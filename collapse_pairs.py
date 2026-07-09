"""Collapse DOUBLED parallel lines into one centerline (base arcgispro-py3 python: ogr+shapely).

Thick double-cased roads on the historical sheets trace as TWO parallel lines a couple of metres
apart, joined by short "rungs" -- a ladder. This step finds pairs of lines that run mutually
parallel within `pair_d` metres over >= `frac` of BOTH their lengths, replaces each pair by its
centerline, and drops the leftover rungs (short lines whose both ends sit on a collapsed pair).

Real narrow strip parcels are much wider than a road casing (>= ~10 m vs ~2-3 m), so the default
pair_d=3.0 does not touch them; the MUTUAL coverage requirement also protects a long border that
merely passes near a short segment. Optional stage, off by default -- check the QA overlay.

Usage: python collapse_pairs.py in.shp out.shp [--pair-d 3.0] [--frac 0.7]
"""
import os, argparse
import numpy as np
from osgeo import ogr
from shapely import wkt, STRtree
from shapely.geometry import LineString, Point, box
from shapely.ops import substring
ogr.UseExceptions()

COS_MIN = 0.8          # how parallel the partner must run (|cos|, ~36 deg)
EPS = 1.0              # tangent window on the partner line (m)
STEP = 1.0             # densify spacing (m)
SIMP = 0.3             # centerline simplify (m)


def _densify(L, step=STEP):
    n = max(1, int(L.length // step))
    return [L.interpolate(k / n, normalized=True) for k in range(n + 1)]


def _tangent(L, s, eps=EPS):
    q0 = L.interpolate(max(0.0, s - eps)); q1 = L.interpolate(min(L.length, s + eps))
    v = np.array([q1.x - q0.x, q1.y - q0.y]); n = np.hypot(*v)
    return v / n if n > 1e-9 else None


def collapse(lines, pair_d=3.0, frac=0.7, step=STEP):
    """Pure core: lines = [LineString]. Returns (out_lines, n_pairs, n_rungs_dropped),
    where out_lines keeps the untouched lines and appends one centerline per collapsed pair."""
    if not lines:
        return [], 0, 0
    tree = STRtree(lines)
    pts = [_densify(L, step) for L in lines]
    cov = {}                                          # (i, j) -> points of i covered by j
    for i, L in enumerate(lines):
        for k, p in enumerate(pts[i]):
            a = pts[i][max(0, k - 1)]; b = pts[i][min(len(pts[i]) - 1, k + 1)]
            v = np.array([b.x - a.x, b.y - a.y]); nv = np.hypot(*v)
            if nv < 1e-9:
                continue
            v /= nv
            best, bd = None, pair_d + 1e-9
            for c in tree.query(box(p.x - pair_d, p.y - pair_d, p.x + pair_d, p.y + pair_d)):
                j = int(c)
                if j == i:
                    continue
                d = p.distance(lines[j])
                if d > pair_d or d >= bd:
                    continue
                t = _tangent(lines[j], lines[j].project(p))
                if t is None or abs(float(np.dot(v, t))) < COS_MIN:
                    continue                          # partner crosses rather than runs alongside
                bd, best = d, j
            if best is not None:
                cov[(i, best)] = cov.get((i, best), 0) + 1

    cand = []                                         # mutually-covered pairs
    for (i, j), c in cov.items():
        if i < j and cov.get((j, i)):
            fi = c / len(pts[i]); fj = cov[(j, i)] / len(pts[j])
            if fi >= frac and fj >= frac:
                cand.append((min(fi, fj), i, j))
    cand.sort(reverse=True)

    consumed, centers = set(), []
    for _, i, j in cand:
        if i in consumed or j in consumed:
            continue
        A, B = (lines[i], lines[j]) if lines[i].length >= lines[j].length else (lines[j], lines[i])
        mids = []
        for p in _densify(A, step):                   # centerline = midway between the two casings
            if p.distance(B) <= pair_d:
                q = B.interpolate(B.project(p))
                mids.append(((p.x + q.x) / 2.0, (p.y + q.y) / 2.0))
        if len(mids) >= 2:
            c = LineString(mids).simplify(SIMP)
            if c.length > 0.5:
                centers.append(c); consumed.add(i); consumed.add(j)

    out, n_rungs = [], 0
    ctree = STRtree(centers) if centers else None
    for k, L in enumerate(lines):
        if k in consumed:
            continue
        if ctree is not None and L.length <= 2.0 * pair_d:     # a ladder rung: both ends on the pair
            p0, p1 = Point(L.coords[0]), Point(L.coords[-1])
            def near(p):
                for c in ctree.query(box(p.x - pair_d, p.y - pair_d, p.x + pair_d, p.y + pair_d)):
                    if p.distance(centers[int(c)]) <= pair_d:
                        return True
                return False
            if near(p0) and near(p1):
                n_rungs += 1; continue
        out.append(L)
    return out + centers, len(centers), n_rungs


def run(src, out, pair_d=3.0, frac=0.7):
    dl = ogr.Open(src); ll = dl.GetLayer(); sr = ll.GetSpatialRef()
    lines = [wkt.loads(f.GetGeometryRef().ExportToWkt()) for f in ll if f.GetGeometryRef()]
    res, n_pairs, n_rungs = collapse(lines, pair_d, frac)

    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out):
        drv.DeleteDataSource(out)
    dso = drv.CreateDataSource(out); lyo = dso.CreateLayer("lines", sr, ogr.wkbLineString)
    lyo.CreateField(ogr.FieldDefn("arc", ogr.OFTInteger)); defn = lyo.GetLayerDefn()
    for k, L in enumerate(res):
        g = ogr.Geometry(ogr.wkbLineString)
        for x, y in L.coords:
            g.AddPoint_2D(float(x), float(y))
        f = ogr.Feature(defn); f.SetField("arc", k); f.SetGeometry(g); lyo.CreateFeature(f); f = None
    dso = None
    print("[collapse] %d lines: %d parallel pairs -> centerlines, %d rungs dropped -> %d lines -> %s"
          % (len(lines), n_pairs, n_rungs, len(res), out))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("out")
    ap.add_argument("--pair-d", type=float, default=3.0, dest="pair_d",
                    help="max separation (m) of a doubled pair; keep BELOW real parcel widths")
    ap.add_argument("--frac", type=float, default=0.7,
                    help="fraction of BOTH lines that must run alongside each other")
    a = ap.parse_args()
    run(a.src, a.out, a.pair_d, a.frac)
