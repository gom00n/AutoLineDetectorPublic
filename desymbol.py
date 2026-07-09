"""Remove traced map-SYMBOL strokes from the detected lines (base arcgispro-py3: ogr+shapely).

The historical sheets carry symbology ink that survives detection because it is connected to
real borders: tick strokes and X-marks between double-line casings, glyph outlines (survey
triangles, arrows, digits touching a line = the 'dot on a line' + wiggle), and small isolated
ornament clusters. None of these can be removed by a length threshold -- a real border piece
crossing a road is the same length as a tick -- so this stage uses the topology the user's
hand-cleaning follows (validated against the M14 edits of 2026-07-07):

  row     -- symbol strokes come in PERIODIC ROWS: >= row_min short strokes (each <= stroke_len,
             with NO collinear continuation at either end) whose midpoints line up (collinear
             within 3 m, spaced <= row_gap): tick/X rungs along a double track, connector
             strokes of ornament chains. A real border piece crossing the double line is
             excluded: it continues collinearly beyond the casings.
  ring    -- a CLOSED feature with a small fat interior (perimeter <= ring_perim, thickness
             area/perimeter >= 0.5 m): a traced glyph circle/blob. Real parcels of that size
             (a ~10 x 12 m ring) do not exist on these rural sheets.
  cluster -- a connected component that is compact (bbox diagonal <= clus_diag) and sits far
             from the rest of the network: an ornament that min_comp couldn't see because the
             cluster's own length exceeds it (e.g. a traced triangle glyph).

Validated against the user's hand-cleaned M14 (2026-07-07): generic short-stroke topology rules
(dead-end angle, graph redundancy, blob pairing) were tried first and FAILED -- their signatures
overlap real border links completely. Only the row + cluster rules survived with ~zero false
positives; everything else stays manual.

Usage: python desymbol.py in.shp out.shp [--stroke-len 12] [--near-d 2.5] [--blob-sep 2.0]
       [--clus-diag 25]
"""
import os, argparse
import numpy as np
from osgeo import ogr
from shapely import wkt, STRtree
from shapely.geometry import LineString, Point, Polygon, box
ogr.UseExceptions()

RND = 1                    # junction key rounding (0.1 m)
COS_CONT = 0.80            # how collinear a neighbour must be to count as a continuation
END_WIN = 3.0              # tangent window at a feature end (m)


def _key(xy):
    return (round(xy[0], RND), round(xy[1], RND))


def _end_dir(cs, e, win=END_WIN):
    """Unit vector pointing INTO the line from end e (0 or -1)."""
    p = np.array(cs[e], float)
    q = p
    idx = range(1, len(cs)) if e == 0 else range(len(cs) - 2, -1, -1)
    for k in idx:                                        # walk until win metres away
        q = np.array(cs[k], float)
        if np.hypot(*(q - p)) >= win:
            break
    v = q - p
    n = np.hypot(*v)
    return v / n if n > 1e-9 else None


def _pass(lines, dropped, stroke_len, row_min, row_gap, ring_perim, clus_diag, protect):
    """One rule pass over the live features; adds to `dropped`, returns True if it grew."""
    n = len(lines)
    before = len(dropped)
    coords = [list(L.coords) for L in lines]
    live = [i for i in range(n) if i not in dropped]
    node = {}
    for i in live:
        for e in (0, -1):
            node.setdefault(_key(coords[i][e]), []).append((i, e))
    tree = STRtree([lines[i] for i in live])
    lmap = {k: i for k, i in enumerate(live)}

    # --- rings: closed glyph circles/blobs ------------------------------------------------
    for i in live:
        if i in protect or len(coords[i]) < 4:
            continue
        if _key(coords[i][0]) == _key(coords[i][-1]) and lines[i].length <= ring_perim:
            try:
                ring = Polygon(coords[i])
                area = ring.area if ring.is_valid else ring.buffer(0).area
            except Exception:
                continue
            if area / max(lines[i].length, 1e-9) >= 0.5:
                dropped[i] = "ring"

    # --- symbol rows: periodic short no-continuation strokes with collinear midpoints ----
    def has_continuation(i):
        for e in (0, -1):
            din = _end_dir(coords[i], e)
            if din is None:
                continue
            for j, f in node.get(_key(coords[i][e]), []):
                if j == i or j in dropped:
                    continue
                dj = _end_dir(coords[j], f)
                if dj is not None and float(np.dot(din, dj)) < -COS_CONT:
                    return True                          # a neighbour continues my direction
        return False

    cand = [i for i in live
            if i not in protect and i not in dropped
            and lines[i].length <= stroke_len and not has_continuation(i)]
    mids = {i: np.array([lines[i].centroid.x, lines[i].centroid.y]) for i in cand}
    # chain candidates whose midpoints sit close together, then test each chain for
    # collinearity (symbol rows march in a line along the track; corner junk scatters)
    parent = {i: i for i in cand}

    def cfind(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    cl = sorted(cand)
    for ai in range(len(cl)):
        for bi in range(ai + 1, len(cl)):
            a, b = cl[ai], cl[bi]
            if np.hypot(*(mids[a] - mids[b])) <= row_gap:
                parent[cfind(a)] = cfind(b)
    chains = {}
    for i in cand:
        chains.setdefault(cfind(i), []).append(i)
    for members in chains.values():
        if len(members) < row_min:
            continue
        pts = np.array([mids[i] for i in members])
        c = pts.mean(axis=0)
        u, s, vt = np.linalg.svd(pts - c)
        resid = np.abs((pts - c) @ vt[1])                # distance to the fitted row axis
        if float(resid.max()) <= 3.0:
            for i in members:
                dropped[i] = "row"

    # --- clusters: compact connected components far from the rest ------------------------
    gparent = {i: i for i in live}

    def find(x):
        while gparent[x] != x:
            gparent[x] = gparent[gparent[x]]
            x = gparent[x]
        return x

    for members in node.values():
        ms = [j for j, _ in members if j not in dropped]
        if not ms:
            continue
        r0 = find(ms[0])
        for j in ms[1:]:
            gparent[find(j)] = r0
    comp = {}
    for i in live:
        if i not in dropped:
            comp.setdefault(find(i), []).append(i)
    for members in comp.values():
        ms = [i for i in members if i not in protect]
        if not ms or len(ms) != len(members):
            continue
        xs, ys = [], []
        for i in ms:
            b = lines[i].bounds
            xs += [b[0], b[2]]; ys += [b[1], b[3]]
        diag = float(np.hypot(max(xs) - min(xs), max(ys) - min(ys)))
        if diag > clus_diag:
            continue
        away = True                                      # every member far from OTHER components
        for i in ms:
            for c in tree.query(lines[i].buffer(4.0)):
                j = lmap[int(c)]
                if j not in dropped and find(j) != find(ms[0]) and lines[i].distance(lines[j]) <= 4.0:
                    away = False
                    break
            if not away:
                break
        if away:
            for i in ms:
                dropped[i] = "cluster"
    return len(dropped) > before


def desymbol(lines, stroke_len=12.0, row_min=3, row_gap=30.0, ring_perim=45.0,
             clus_diag=25.0, protect=None):
    """Pure core. lines = [LineString]; protect = optional set of indices never dropped.
    Rules iterate until stable (dropping a row of connectors isolates the glyph blobs it
    chained, which the next pass removes). Returns (kept_idx, [(idx, reason), ...])."""
    protect = protect or set()
    dropped = {}
    for _ in range(3):
        if not _pass(lines, dropped, stroke_len, row_min, row_gap, ring_perim,
                     clus_diag, protect):
            break
    kept = [i for i in range(len(lines)) if i not in dropped]
    return kept, sorted(dropped.items())


def run(src, out, stroke_len=12.0, row_min=3, row_gap=30.0, ring_perim=45.0, clus_diag=25.0):
    dl = ogr.Open(src); ll = dl.GetLayer(); sr = ll.GetSpatialRef()
    feats = [(f.GetField("arc") if f.GetFieldIndex("arc") >= 0 else f.GetFID(),
              wkt.loads(f.GetGeometryRef().ExportToWkt())) for f in ll if f.GetGeometryRef()]
    lines = [L for _, L in feats]
    kept, dropped = desymbol(lines, stroke_len, row_min, row_gap, ring_perim, clus_diag)

    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out):
        drv.DeleteDataSource(out)
    dso = drv.CreateDataSource(out); lyo = dso.CreateLayer("lines", sr, ogr.wkbLineString)
    lyo.CreateField(ogr.FieldDefn("arc", ogr.OFTInteger)); defn = lyo.GetLayerDefn()
    for i in kept:
        g = ogr.Geometry(ogr.wkbLineString)
        for x, y in lines[i].coords:
            g.AddPoint_2D(float(x), float(y))
        f = ogr.Feature(defn); f.SetField("arc", int(feats[i][0])); f.SetGeometry(g)
        lyo.CreateFeature(f); f = None
    dso = None
    import collections
    tags = collections.Counter(t for _, t in dropped)
    print("[desymbol] %d lines -> dropped %d (%s) -> %d -> %s"
          % (len(lines), len(dropped), ", ".join("%s=%d" % kv for kv in sorted(tags.items())) or "none",
             len(kept), out))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("out")
    ap.add_argument("--stroke-len", type=float, default=12.0, dest="stroke_len",
                    help="max length (m) of a symbol stroke (ticks/X-marks)")
    ap.add_argument("--row-min", type=int, default=3, dest="row_min",
                    help="min strokes on one host pair to count as a symbol row")
    ap.add_argument("--row-gap", type=float, default=30.0, dest="row_gap",
                    help="max spacing (m) between consecutive strokes of a row")
    ap.add_argument("--ring-perim", type=float, default=45.0, dest="ring_perim",
                    help="max perimeter (m) of a closed glyph ring")
    ap.add_argument("--clus-diag", type=float, default=25.0, dest="clus_diag",
                    help="max bbox diagonal (m) of an isolated compact ornament cluster")
    a = ap.parse_args()
    run(a.src, a.out, a.stroke_len, a.row_min, a.row_gap, a.ring_perim, a.clus_diag)
