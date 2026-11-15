"""Bridge gaps in a 1-px skeleton to improve line CONTINUITY (run with venv_skimage python).

The borders fragment into hundreds of pieces because faded/broken ink leaves gaps that
isotropic closing (skeletonize.py's binary_closing) can't span without fusing parallel
lines. This links skeleton ENDPOINTS across a gap only when they are close AND roughly
collinear & facing each other -- so a broken border is rejoined but neighbouring lines are
not merged. Two passes:

  1. endpoint <-> endpoint : rejoin a clean break in a single line (most common).
  2. endpoint  -> line     : extend a stub that stops just short of a crossing border
                             (T-junction), casting a short ray along the stub's direction.

The connectors are drawn 1-px and the skeleton is NOT re-skeletonized (skimage.skeletonize is
not idempotent on an already-thin skeleton); only spurs are re-pruned. Output is a Byte 0/1
skeleton, same grid as input.

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
NB8 = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))


def endpoints(sk):
    """(r,c) of skeleton pixels with exactly one 8-neighbour."""
    nbr = ndimage.convolve(sk.astype(np.uint8), _8, mode="constant")
    return np.argwhere(sk & (nbr == 2))


def walk_back(sk, ep, n):
    """Up to n pixels of the endpoint's OWN stroke, walked from the endpoint inwards along the
    skeleton and stopping at the first junction. The endpoint itself is not included."""
    H, W = sk.shape
    prev, cur = None, (int(ep[0]), int(ep[1]))
    seen, path = {cur}, []
    for _ in range(n):
        r, c = cur
        nxt = [(r + dr, c + dc) for dr, dc in NB8
               if 0 <= r + dr < H and 0 <= c + dc < W and sk[r + dr, c + dc]
               and (r + dr, c + dc) not in seen]
        if len(nxt) > 1 and prev is not None:
            # a staircase corner offers the same next pixel twice (orthogonal + diagonal);
            # drop the candidates that merely neighbour the pixel we came from
            far = [q for q in nxt if max(abs(q[0] - prev[0]), abs(q[1] - prev[1])) > 1]
            nxt = far or nxt
        if len(nxt) == 2:
            # exactly two, a 4-neighbour and the diagonal right beside it: one direction, not
            # a fork (three candidates in a row is a T-junction and must stop the walk)
            orth = [q for q in nxt if q[0] == r or q[1] == c]
            if len(orth) == 1 and all(abs(q[0] - orth[0][0]) + abs(q[1] - orth[0][1]) <= 1 for q in nxt):
                nxt = orth
        if len(nxt) != 1:                     # dead end, or a real junction: stop here
            break
        prev, cur = cur, nxt[0]
        seen.add(cur); path.append(cur)
    return path


def outward_dirs(sk, eps, R=10):
    """Unit vector at each endpoint pointing AWAY from its line: endpoint minus the centroid of
    the last R pixels of its OWN stroke.

    The centroid used to be taken over every skeleton pixel within R of the endpoint. That
    includes the line the stub is about to meet: for a stub ending 2-8 px short of a crossing
    border, the border's pixels outweighed the stub's own, the centroid landed BEYOND the
    endpoint and the vector flipped -- so pass 2 cast its ray backwards into the stub, "hit" it
    at step 2, and reported a link without drawing anything. Short T-gaps were never bridged."""
    dirs = np.zeros((len(eps), 2))
    for k, ep in enumerate(eps):
        own = walk_back(sk, ep, R)
        if not own:
            continue
        v = np.asarray(ep, float) - np.asarray(own, float).mean(axis=0)
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
        pr, pc = int(round(r0 + d[0])), int(round(c0 + d[1]))
        for step in range(2, max_gap + 1):
            r, c = int(round(r0 + d[0] * step)), int(round(c0 + d[1] * step))
            if not (0 <= r < H and 0 <= c < W):
                break
            hit = None
            if sk[r, c]:                      # ray reached existing line
                hit = (r, c)
            elif step > 2 and r != pr and c != pc:
                # a diagonal step can slip BETWEEN two pixels of an 8-connected diagonal line
                # without landing on either; the two cells it squeezes past are the line
                if sk[pr, c]:
                    hit = (pr, c)
                elif sk[r, pc]:
                    hit = (r, pc)
            if hit is not None:
                rr, cc = draw_line(int(ep[0]), int(ep[1]), hit[0], hit[1])
                sk[rr, cc] = 1
                nb += 1
                break
            pr, pc = r, c
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
