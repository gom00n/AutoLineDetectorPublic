"""Skeletonize a binary line mask to 1-px centerlines (run with venv_skimage python).

close_r bridges small ink breaks by morphological closing, but it also welds DOUBLED parallel
casings (two lines 2-4 px apart, e.g. red-ink double-line tracks) into one blob whose skeleton
is a single weaving centerline. max_fill keeps the closing local: any connected patch of
pixels the closing ADDED that is bigger than max_fill px is reverted -- a repaired line break
adds a few dozen px, the filled corridor between two casings adds hundreds, so the breaks stay
repaired and the casings stay two separate lines. 0 disables the guard (old behaviour).

Usage: python skeletonize.py in_mask.tif out_skel.tif [close_r] [prune_len] [max_fill]"""
import sys
import numpy as np
from osgeo import gdal
from scipy import ndimage
from skimage.morphology import skeletonize
gdal.UseExceptions()

OPTS = ["COMPRESS=LZW", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER"]
_8 = np.ones((3, 3), np.uint8)


def _disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r


def prune_spurs(sk, max_len, passes=4):
    """Remove 'hair': short side-branches that hang off the line network (a free endpoint
    AND adjacent to a junction). Long lines and standalone short segments are kept."""
    sk = sk.astype(bool)
    for _ in range(passes):
        nbr = ndimage.convolve(sk.astype(np.uint8), _8, mode="constant")
        endpoints = sk & (nbr == 2)            # exactly 1 neighbour
        junctions = sk & (nbr >= 4)            # 3+ neighbours
        arcs = sk & ~junctions                 # split network at junctions
        lab, n = ndimage.label(arcs, structure=np.ones((3, 3), int))
        if n == 0:
            break
        sizes = np.bincount(lab.ravel())
        has_ep = set(np.unique(lab[endpoints])) - {0}
        jadj = set(np.unique(lab[ndimage.binary_dilation(junctions, _8) & arcs])) - {0}
        spur = np.zeros(n + 1, bool)
        for cid in has_ep & jadj:
            if sizes[cid] < max_len:
                spur[cid] = True
        if not spur.any():
            break
        sk &= ~spur[lab]
    return sk.astype(np.uint8)


def run(src, dst, close_r=3, prune_len=40, max_fill=0):
    ds = gdal.Open(src)
    a = ds.GetRasterBand(1).ReadAsArray() > 0
    if close_r:                       # bridge small ink breaks before thinning
        closed = ndimage.binary_closing(a, structure=_disk(close_r))
        if max_fill:                  # ...but un-fill big patches (inter-casing corridors)
            lab, n = ndimage.label(closed & ~a, structure=np.ones((3, 3), int))
            if n:
                sizes = np.bincount(lab.ravel()); sizes[0] = 0
                big = sizes > max_fill
                if big.any():
                    closed &= ~big[lab]
                    print("[skel] reverted %d closing fills > %d px (kept doubled casings apart)"
                          % (int(big.sum()), max_fill))
        a = closed
    sk = skeletonize(a).astype(np.uint8)
    if prune_len:
        before = int(sk.sum())
        sk = prune_spurs(sk, prune_len)
        print("[skel] pruned spurs <%dpx: %d -> %d px" % (prune_len, before, int(sk.sum())))
    out = gdal.GetDriverByName("GTiff").Create(
        dst, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte, options=OPTS)
    out.SetGeoTransform(ds.GetGeoTransform())
    out.SetProjection(ds.GetProjection())
    b = out.GetRasterBand(1)
    b.WriteArray(sk)
    b.SetNoDataValue(0)
    out.FlushCache()
    print("[skel] %d ink -> %d skeleton px  %s" % (int((a > 0).sum()), int(sk.sum()), dst))


if __name__ == "__main__":
    close_r = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    prune_len = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    max_fill = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    run(sys.argv[1], sys.argv[2], close_r, prune_len, max_fill)
