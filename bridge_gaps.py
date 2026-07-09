"""Bridge gaps in a 1-px skeleton to improve line CONTINUITY (run with venv_skimage python).

The borders fragment into hundreds of pieces because faded/broken ink leaves gaps that
isotropic closing (skeletonize.py's binary_closing) can't span without fusing parallel
lines. This links skeleton ENDPOINTS across a gap only when they are close AND roughly
collinear & facing each other -- so a broken border is rejoined but neighbouring lines are
not merged. Two passes:

  1. endpoint <-> endpoint : rejoin a clean break in a single line (most common).
  2. endpoint  -> line     : extend a stub that stops just short of a crossing border
                             (T-junction), casting a short ray along the stub's direction.

Then re-skeletonize (drawn bridges may be 1px-thick already, but junctions need cleanup)
and re-prune spurs. Output is a Byte 0/1 skeleton, same grid as input.

Run (between skeletonize.py and trace_lines.py):
  python bridge_gaps.py in_skel.tif out_skel.tif [max_gap=25] [cos=0.80] [ray=0|1]
"""
import sys
import numpy as np
from osgeo import gdal
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage.draw import line as draw_line

from skeletonize import prune_spurs  # reuse the same de-hairing

gdal.UseExceptions()
OPTS = ["COMPRESS=LZW", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER"]
_8 = np.ones((3, 3), np.uint8)


def endpoints(sk):
    """(r,c) of skeleton pixels with exactly one 8-neighbour."""
    nbr = ndimage.convolve(sk.astype(np.uint8), _8, mode="constant")
    return np.argwhere(sk & (nbr == 2))


def outward_dirs(sk, eps, R=10):
    """Unit vector at each endpoint pointing AWAY from its line (endpoint - local centroid).
    Robust for both straight and gently curved strokes."""
    ys, xs = np.nonzero(sk)
    tree = cKDTree(np.c_[ys, xs])
    dirs = np.zeros((len(eps), 2))
    near = tree.query_ball_point(eps, R)
    for k, (ep, idx) in enumerate(zip(eps, near)):
        pts = np.c_[ys[idx], xs[idx]]
        v = ep - pts.mean(axis=0)
        n = np.hypot(v[0], v[1])
        if n > 1e-6:
            dirs[k] = v / n
    return dirs


def bridge_ep_ep(sk, eps, dirs, max_gap, cos_thr):
    """Pass 1: connect endpoint pairs that face each other and are roughly collinear."""
    pts = eps.astype(float)
    tree = cKDTree(pts)
    cands = []
    for i, j in tree.query_pairs(max_gap):
        g = pts[j] - pts[i]
        d = np.hypot(g[0], g[1])
        if d < 1.5:
            continue
        u = g / d
        ai = dirs[i].dot(u)        # i's stub points toward j
        aj = dirs[j].dot(-u)       # j's stub points toward i
        if ai > cos_thr and aj > cos_thr:
            cands.append(((ai + aj) - d / max_gap, i, j))   # prefer aligned & short
    cands.sort(reverse=True)
    used, nb = set(), 0
    for _, i, j in cands:
        if i in used or j in used:
            continue
        rr, cc = draw_line(int(pts[i, 0]), int(pts[i, 1]), int(pts[j, 0]), int(pts[j, 1]))
        sk[rr, cc] = 1
        used.update((i, j))
        nb += 1
    return nb, used


def bridge_ep_line(sk, eps, dirs, used, max_gap, cos_thr):
    """Pass 2: extend each still-unmatched endpoint along its outward direction; if the ray
    hits another part of the skeleton within max_gap, draw the connector (T-junction)."""
    H, W = sk.shape
    nb = 0
    for k, ep in enumerate(eps):
        if k in used:
            continue
        d = dirs[k]
        if not d.any():
            continue
        r0, c0 = float(ep[0]), float(ep[1])
        for step in range(2, max_gap + 1):
            r, c = int(round(r0 + d[0] * step)), int(round(c0 + d[1] * step))
            if not (0 <= r < H and 0 <= c < W):
                break
            if sk[r, c]:                      # ray reached existing line
                rr, cc = draw_line(int(ep[0]), int(ep[1]), r, c)
                sk[rr, cc] = 1
                nb += 1
                break
    return nb


def run(src, dst, max_gap=25, cos_thr=0.80, ray=True, prune_len=40):
    ds = gdal.Open(src)
    sk = (ds.GetRasterBand(1).ReadAsArray() > 0).astype(np.uint8)
    n0 = int(sk.sum())

    eps = endpoints(sk)
    print("[bridge] %d endpoints" % len(eps))
    if len(eps) >= 2:
        dirs = outward_dirs(sk, eps)
        nb1, used = bridge_ep_ep(sk, eps, dirs, max_gap, cos_thr)
        nb2 = bridge_ep_line(sk, eps, dirs, used, max_gap, cos_thr) if ray else 0
        print("[bridge] linked %d endpoint-endpoint, %d endpoint-line (max_gap=%d, cos=%.2f)"
              % (nb1, nb2, max_gap, cos_thr))

    # NB: do NOT re-skeletonize globally -- skimage.skeletonize is not idempotent on an
    # already-thin skeleton and fragments 8-connected diagonal chains (202 -> 417 endpoints
    # on M15). The connectors from draw_line are already 1-px, so only re-prune any spurs.
    if prune_len:
        sk = prune_spurs(sk, prune_len)
    print("[bridge] skeleton %d -> %d px; endpoints %d -> %d"
          % (n0, int(sk.sum()), len(eps), len(endpoints(sk))))

    out = gdal.GetDriverByName("GTiff").Create(
        dst, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte, options=OPTS)
    out.SetGeoTransform(ds.GetGeoTransform())
    out.SetProjection(ds.GetProjection())
    b = out.GetRasterBand(1)
    b.WriteArray(sk)
    b.SetNoDataValue(0)
    out.FlushCache()
    print("[bridge] wrote %s" % dst)


if __name__ == "__main__":
    a = sys.argv
    max_gap = int(a[3]) if len(a) > 3 else 25
    cos_thr = float(a[4]) if len(a) > 4 else 0.80
    ray = (a[5] != "0") if len(a) > 5 else True
    run(a[1], a[2], max_gap, cos_thr, ray)
