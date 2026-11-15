"""Motzag-GDB-corridor line extraction for ALL projects (base arcgispro-py3 python).

Two ink methods, same corridor + windowed-read machinery:
  --ink black (black-ink): black borders are dark in ALL channels, so `max(R,G,B) < cap` isolates them
       and rejects colored confusers (red numbers leak 0%; purple ~18%, killed by the B-G exclusion).
  --ink red (red-ink): red ink over a yellow wash reads only ~15-20 DN darker than the wash, so
       a global threshold misses it. An ADAPTIVE local-mean on the Red band (R < local_mean - OFFSET,
       capped) recovers it. `--srcdir` points at the scan folder; output folder is set by --out / GUI.

CORRIDOR comes from a GDB: the user picks which feature layers to scan; in them we take the
polygons whose attribute Motzag == <scan number from the file name> (M8_9076 -> 8), buffer them,
and that defines the zone to detect in. We DON'T read the whole scan: the corridor's bounding box
gives a pixel window and we read only that (~2-3% of the sheet). If no Motzag==scan polygons exist
in the chosen layers (or the source .tif is missing) we NOTIFY and skip the sheet (exit 3) rather
than abort the batch. A genuine config error (e.g. the GDB can't be opened) exits 2 so the GUI
surfaces it as an error instead of a silent skip. Success exits 0.

CRS (Israeli Grid 05/12, EPSG 6991) is stamped on the windowed mask, so the rest of the pipeline
(skeletonize -> bridge_gaps -> trace_lines -> merge_lines -> filter_components) just works.

Run:  python extract_9a.py <stem> [--ink black|red] [--cap N] [--offset 14] [--block 101]
      [--min-px N] [--buf 35] [--gdb PATH] [--layers "L1,L2"] [--srcdir DIR] [--out MASK]
"""
import os, re, sys, argparse
import numpy as np
from osgeo import gdal, ogr, osr
from scipy import ndimage
gdal.UseExceptions()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config                                        # single source of truth for paths/defaults
SRCDIR, WORK = config.NINE_SRCDIR, config.NINE_OUT
GDB, LAYERS = config.NINE_GDB, config.NINE_LAYERS

# --- ink detection methods (both share the Motzag-GDB corridor + per-polygon windowed read) ---
#   "black" (black-ink): a pixel is ink if max(R,G,B) < cap, minus bluish admin lines (B-G > PURPLE_DB).
#   "red"   (red-ink): ADAPTIVE local-mean on the Red band -- ink = R darker than its local
#           background by >= OFFSET and below cap. Recovers red ink drawn over a colored wash, which
#           a global threshold misses.
CAP_BLACK = 145      # default ceiling for black ink (max RGB < cap)
CAP_RED = 212        # default Red ceiling for red ink (excludes bright paper/yellow)
PURPLE_DB = 20       # black only: exclude bluish admin dash-dot (B - G > PURPLE_DB)
BLOCK = 101          # red only: local-background window (px)
OFFSET = 14          # red only: how much darker than local background a line must be
MIN_PX = 30          # drop components smaller than this (coarse 0.85 m/px lines)
BUF_M = 35.0         # buffer the Motzag corridor (catch borders sitting just outside the polygon)
EPSG = 6991
SKIP = 3             # exit code: no target area for this sheet (notify, don't abort the batch)
CREATE_OPTS = ["COMPRESS=LZW", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER"]


def scan_number(stem):
    """M8_9076 -> 8 (the Motzag value selecting this scan's areas of interest)."""
    return int(re.match(r"[Mm](\d+)", stem).group(1))


def geotransform_from_worldfile(tfwx_path):
    with open(tfwx_path) as fh:
        A, D, B, E, C, F = [float(x) for x in fh.read().split()][:6]
    return (C - 0.5 * A - 0.5 * B, A, B, F - 0.5 * D - 0.5 * E, D, E)


def source_paths(sdir, stem):
    """Resolve <stem>'s scan + world file: some sheets are delivered as .tif/.tfwx, others as
    .jpg/.jgwx (same project, mixed formats). Try both, in that order."""
    for ext, wext in ((".tif", ".tfwx"), (".jpg", ".jgwx")):
        src = os.path.join(sdir, stem + ext)
        if os.path.exists(src):
            return src, os.path.join(sdir, stem + wext)
    return os.path.join(sdir, stem + ".tif"), os.path.join(sdir, stem + ".tfwx")


def world_to_px(gt, X, Y):
    """Invert a (possibly rotated) geotransform: world (X,Y) -> (col, row)."""
    a, b, d, e = gt[1], gt[2], gt[4], gt[5]
    det = a * e - b * d
    dx, dy = X - gt[0], Y - gt[3]
    return (e * dx - b * dy) / det, (-d * dx + a * dy) / det


def pixel_size(gt):
    """True ground size of a pixel in metres: sqrt(|det|) of the geotransform's 2x2 part.

    abs(gt[1]) is the pixel size only on a NORTH-UP grid. On a rotated world file it is
    m/px * cos(rotation) -- 0.083 against a true 1.015 m/px on one sheet scanned sideways."""
    return abs(gt[1] * gt[5] - gt[2] * gt[4]) ** 0.5


def size_filter(mask, min_px):
    lab, n = ndimage.label(mask, structure=np.ones((3, 3), int))
    sizes = np.bincount(lab.ravel()); keep = sizes >= min_px; keep[0] = False
    return keep[lab].astype(np.uint8), n, int(keep.sum())


def motzag_geoms(gdb_ds, layers, scan):
    """Cloned Motzag==scan polygon geometries from the chosen GDB layers + combined world extent.
    Warns (doesn't fail) on a missing layer or a layer without a Motzag field."""
    geoms, ext = [], None
    for name in layers:
        lyr = gdb_ds.GetLayerByName(name)
        if lyr is None:
            print("[9A] WARNING: layer '%s' not found in GDB -- skipped" % name); continue
        if lyr.GetLayerDefn().GetFieldIndex("Motzag") < 0:
            print("[9A] WARNING: layer '%s' has no Motzag field -- skipped" % name); continue
        lyr.SetAttributeFilter("Motzag = %d" % scan)
        for f in lyr:
            g = f.GetGeometryRef()
            if g is None:
                continue
            g = g.Clone()
            g.FlattenTo2D()                           # drop Z/M so envelope/rasterize behave
            if g.HasCurveGeometry():                  # true-curve polygon (e.g. MULTISURFACE):
                g = g.GetLinearGeometry()              # RasterizeLayer burns 0px on curves -- linearize
            geoms.append(g)
            e = g.GetEnvelope()                       # (minX, maxX, minY, maxY)
            ext = e if ext is None else (min(ext[0], e[0]), max(ext[1], e[1]),
                                         min(ext[2], e[2]), max(ext[3], e[3]))
        lyr.SetAttributeFilter(None)
    return geoms, ext


def _px_window(gt, ex, buf_m, mrg, bounds):
    """World envelope ex=(minX,maxX,minY,maxY) (+buf) -> clamped pixel window (c0,r0,c1,r1)."""
    cr = [world_to_px(gt, x, y) for x in (ex[0] - buf_m, ex[1] + buf_m)
          for y in (ex[2] - buf_m, ex[3] + buf_m)]
    cs = [c for c, _ in cr]; rs = [r for _, r in cr]
    c0 = max(bounds[0], int(np.floor(min(cs))) - mrg); r0 = max(bounds[1], int(np.floor(min(rs))) - mrg)
    c1 = min(bounds[2], int(np.ceil(max(cs))) + mrg); r1 = min(bounds[3], int(np.ceil(max(rs))) + mrg)
    return c0, r0, c1, r1


def _rasterize_geom(geom, sw, sh, gtsub, proj):
    """Burn a single polygon geometry into a sw x sh bool array on gtsub."""
    vds = ogr.GetDriverByName("Memory").CreateDataSource("m")
    srs = osr.SpatialReference(); srs.ImportFromWkt(proj)
    lyr = vds.CreateLayer("l", srs, ogr.wkbPolygon)
    f = ogr.Feature(lyr.GetLayerDefn()); f.SetGeometry(geom); lyr.CreateFeature(f)
    mem = gdal.GetDriverByName("MEM").Create("", sw, sh, 1, gdal.GDT_Byte)
    mem.SetGeoTransform(gtsub); mem.SetProjection(proj)
    gdal.RasterizeLayer(mem, [1], lyr, burn_values=[1])
    return mem.GetRasterBand(1).ReadAsArray().astype(bool)


def process(stem, ink="black", cap=None, offset=OFFSET, block=BLOCK, purple_db=PURPLE_DB,
            min_px=MIN_PX, buf_m=BUF_M, gdb=GDB, layers=LAYERS, srcdir=None, out_path=None):
    ink = "red" if str(ink).lower().startswith("r") else "black"
    if cap is None:
        cap = CAP_RED if ink == "red" else CAP_BLACK
    sdir = srcdir or SRCDIR
    src, wf = source_paths(sdir, stem)
    if not os.path.exists(src):
        print("[9A] NOTICE: source scan %s not found -- skipping %s (not an error)." % (src, stem))
        return None
    gt = geotransform_from_worldfile(wf)
    sr = osr.SpatialReference(); sr.ImportFromEPSG(EPSG); proj = sr.ExportToWkt()
    dst = out_path or os.path.join(WORK, stem, stem + "_mask.tif")    # intermediate: GUI sends it to temp
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    scan = scan_number(stem)

    ds = gdal.Open(src); nx, ny = ds.RasterXSize, ds.RasterYSize

    # 1) corridor polygons from the chosen GDB layers (Motzag == scan)
    gdb_ds = ogr.Open(gdb)
    if gdb_ds is None:
        raise RuntimeError("cannot open GDB %s -- check the path/config (this is a real error, "
                           "not a 'no target here' skip)" % gdb)
    geoms, ext = motzag_geoms(gdb_ds, layers, scan)
    if not geoms:
        print("[9A] NOTICE: no Motzag=%d polygons in layers %s -- skipping %s (not an error)."
              % (scan, layers, stem))
        return None

    # 2) overall output grid = corridor bbox window; gt_win keeps georef correct
    mrg = 8
    col0, row0, col1, row1 = _px_window(gt, ext, buf_m, mrg, (0, 0, nx, ny))
    ncols, nrows = col1 - col0, row1 - row0
    if ncols <= 0 or nrows <= 0:
        print("[9A] NOTICE: corridor window falls outside %s -- skipping." % stem); return None
    gt_win = (gt[0] + col0 * gt[1] + row0 * gt[2], gt[1], gt[2],
              gt[3] + col0 * gt[4] + row0 * gt[5], gt[4], gt[5])
    print("[9A] %s: ink=%s  %d Motzag=%d polygons; output grid %dx%d"
          % (stem, ink, len(geoms), scan, ncols, nrows))

    # 3) read ONE small window PER polygon (the targets are scattered, so a single bbox is useless);
    #    detect ink (black global OR red adaptive) + clip to the buffered polygon, OR into the mask.
    mask = np.zeros((nrows, ncols), np.uint8)
    # pixel_size, not abs(gt[1]): on a rotated grid the old form turned a 20 m buffer into
    # 241 px instead of 20, growing the corridor to fill each polygon's whole window
    buf_px = max(1, int(round(buf_m / max(pixel_size(gt), 1e-9))))
    struct = ndimage.generate_binary_structure(2, 1)
    read_px = 0
    for g in geoms:
        sc0, sr0, sc1, sr1 = _px_window(gt, g.GetEnvelope(), buf_m, mrg, (col0, row0, col1, row1))
        sw, sh = sc1 - sc0, sr1 - sr0
        if sw <= 0 or sh <= 0:
            continue
        R = ds.GetRasterBand(1).ReadAsArray(sc0, sr0, sw, sh).astype(np.float32)
        if ink == "red":                                    # adaptive local-mean on the Red band
            localm = ndimage.uniform_filter(R, size=block)
            sub = (R < localm - offset) & (R < cap)
        else:                                               # black: dark in all channels, drop purple
            G = ds.GetRasterBand(2).ReadAsArray(sc0, sr0, sw, sh).astype(np.float32)
            B = ds.GetRasterBand(3).ReadAsArray(sc0, sr0, sw, sh).astype(np.float32)
            sub = (np.maximum(np.maximum(R, G), B) < cap) & ~((B - G) > purple_db)
        gtsub = (gt[0] + sc0 * gt[1] + sr0 * gt[2], gt[1], gt[2],
                 gt[3] + sc0 * gt[4] + sr0 * gt[5], gt[4], gt[5])
        cor = _rasterize_geom(g, sw, sh, gtsub, proj)
        if buf_m > 0:
            cor = ndimage.binary_dilation(cor, struct, iterations=buf_px)
        sub = (sub & cor).astype(np.uint8)
        mask[sr0 - row0:sr1 - row0, sc0 - col0:sc1 - col0] |= sub
        read_px += sw * sh
    ds = None
    print("[9A] read %.1f%% of the sheet across %d polygon windows (vs 100%% whole-scan)"
          % (100.0 * read_px / (nx * ny), len(geoms)))

    # 4) size filter + write the windowed mask (georeferenced via gt_win + EPSG)
    mask, ncc, nkeep = size_filter(mask, min_px)
    print("[9A] cap=%d  kept %d/%d components (>= %dpx) -> %.3f%% ink" % (cap, nkeep, ncc, min_px, 100 * mask.mean()))
    out = gdal.GetDriverByName("GTiff").Create(dst, ncols, nrows, 1, gdal.GDT_Byte, options=CREATE_OPTS)
    out.SetGeoTransform(gt_win); out.SetProjection(proj)
    bnd = out.GetRasterBand(1); bnd.WriteArray(mask); bnd.SetNoDataValue(0); out.FlushCache(); out = None
    print("[9A] wrote %s" % dst)
    return dst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stem")
    ap.add_argument("--ink", default="black", choices=["black", "red"],
                    help="black = max(R,G,B)<cap (black-ink); red = adaptive Red-band (red-ink)")
    ap.add_argument("--cap", type=int, default=None, help="ceiling (default 145 black / 212 red)")
    ap.add_argument("--offset", type=int, default=OFFSET, help="red only: darker-than-local margin")
    ap.add_argument("--block", type=int, default=BLOCK, help="red only: local-background window px")
    ap.add_argument("--min-px", type=int, default=MIN_PX, dest="min_px")
    ap.add_argument("--buf", type=float, default=BUF_M)
    ap.add_argument("--gdb", default=GDB)
    ap.add_argument("--layers", default=",".join(LAYERS), help="comma-separated GDB layer names")
    ap.add_argument("--srcdir", default=None, help="folder holding <stem>.tif (default: project srcdir)")
    ap.add_argument("--out", default=None, help="mask output path (default: work/<stem>/<stem>_mask.tif)")
    a = ap.parse_args()
    sel = [x.strip() for x in a.layers.split(",") if x.strip()]
    try:
        res = process(a.stem, ink=a.ink, cap=a.cap, offset=a.offset, block=a.block,
                      min_px=a.min_px, buf_m=a.buf, gdb=a.gdb, layers=sel,
                      srcdir=a.srcdir, out_path=a.out)
    except RuntimeError as e:                       # hard config error (e.g. bad GDB): exit 2, NOT skip
        print("[9A] ERROR: %s" % e)
        sys.exit(2)
    sys.exit(0 if res else SKIP)                    # 0 = wrote mask; 3 = legit skip (no target / no source)
