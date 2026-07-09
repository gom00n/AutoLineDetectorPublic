"""QA overlay (any project): detected lines over the faded scan, cropped to the line extent.
Prefers <stem>_lines_snapped.shp when it exists and colors by its `snapped` field --
GREEN = snapped to a reference boundary, RED = kept as ink -- so one glance at the PNG shows
how the snap stage did. Falls back to <stem>_lines.shp (all red). Reads the scan from --srcdir
and the result/preview in --workdir. Rotation-aware (world -> full-res pixels via the .tfwx affine).

Usage: python overlay_9a.py <stem> [--suffix _lines] [--srcdir DIR] [--workdir DIR] [--win cx cy w]
"""
import os, sys, argparse
import numpy as np
from osgeo import gdal, ogr
from PIL import Image, ImageDraw
gdal.UseExceptions()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from extract_9a import source_paths
DEF_SRCDIR, DEF_WORK = config.NINE_SRCDIR, config.NINE_OUT


def gt_from_tfwx(p):
    A, D, B, E, C, F = [float(x) for x in open(p).read().split()][:6]
    return (C - 0.5 * A - 0.5 * B, A, B, F - 0.5 * D - 0.5 * E, D, E)


ap = argparse.ArgumentParser()
ap.add_argument("stem")
ap.add_argument("--suffix", default=None, help="result suffix; default: _lines_snapped if present, else _lines")
ap.add_argument("--srcdir", default=DEF_SRCDIR, help="folder holding <stem>.tif (+ .tfwx)")
ap.add_argument("--workdir", default=None, help="folder holding <stem><suffix>.shp; preview written here")
ap.add_argument("--win", nargs=3, type=int, default=None, metavar=("CX", "CY", "W"))
a = ap.parse_args()

wd = a.workdir or os.path.join(DEF_WORK, a.stem)
suffix = a.suffix
if suffix is None:                                   # prefer the snapped result when it exists
    snapped_shp = os.path.join(wd, a.stem + "_lines_snapped.shp")
    suffix = "_lines_snapped" if os.path.exists(snapped_shp) else "_lines"
shp = os.path.join(wd, a.stem + suffix + ".shp")
src, wf = source_paths(a.srcdir, a.stem)
gt = gt_from_tfwx(wf)
ds = gdal.Open(src); nx, ny = ds.RasterXSize, ds.RasterYSize

# world -> pixel (rotation-aware) via affine inverse
a_, b_, d_, e_ = gt[1], gt[2], gt[4], gt[5]; det = a_ * e_ - b_ * d_


def world_to_px(x, y):
    dx, dy = x - gt[0], y - gt[3]
    return (e_ * dx - b_ * dy) / det, (-d_ * dx + a_ * dy) / det


dl = ogr.Open(shp); ly = dl.GetLayer()
has_snap = ly.GetLayerDefn().GetFieldIndex("snapped") >= 0
pix = []                                             # (points, snapped 0/1)
for f in ly:
    g = f.GetGeometryRef()
    fl = f.GetField("snapped") if has_snap else 0
    pix.append(([world_to_px(x, y) for x, y in g.GetPoints()], 1 if fl else 0))
if not pix:
    print("[overlay] no lines in %s" % shp); sys.exit(0)
allc = np.array([p for ln, _ in pix for p in ln])

if a.win:
    cx, cy, win = a.win
    x0, y0 = max(0, cx - win // 2), max(0, cy - win // 2)
    x1, y1 = min(nx, x0 + win), min(ny, y0 + win)
else:
    x0, y0 = max(0, int(allc[:, 0].min())), max(0, int(allc[:, 1].min()))
    x1, y1 = min(nx, int(allc[:, 0].max())), min(ny, int(allc[:, 1].max()))
W, H = max(1, x1 - x0), max(1, y1 - y0)
F = max(1, max(W, H) // 1600)
ow, oh = max(1, W // F), max(1, H // F)
rd = lambda bnd: ds.GetRasterBand(bnd).ReadAsArray(x0, y0, W, H, buf_xsize=ow, buf_ysize=oh).astype(np.float32)
g = 0.299 * rd(1) + 0.587 * rd(2) + 0.114 * rd(3)
g = np.clip(175 + (g - g.mean()) * 0.40, 120, 255).astype(np.uint8)
img = Image.fromarray(np.dstack([g, g, g]))
dr = ImageDraw.Draw(img)
COL = {0: (225, 20, 20), 1: (20, 160, 40)}           # red = kept as ink, green = snapped
for fl_draw in (0, 1):                               # snapped drawn on top
    for ln, fl in pix:
        if fl != fl_draw:
            continue
        pts = [((c - x0) / F, (r - y0) / F) for c, r in ln]
        if len(pts) >= 2:
            dr.line(pts, fill=COL[fl], width=2)
out = os.path.join(wd, a.stem + "_preview.png")
img.save(out)
n_snap = sum(fl for _, fl in pix)
tag = ("  (%d snapped green / %d ink red)" % (n_snap, len(pix) - n_snap)) if has_snap else ""
print("[overlay] %d lines%s, window (%d,%d)+%dx%d /%d -> %s" % (len(pix), tag, x0, y0, W, H, F, out))
