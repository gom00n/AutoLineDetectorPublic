"""Close near-meeting borders and build parcel POLYGONS from the detected lines
(base arcgispro-py3 python: ogr + shapely -- FeatureToPolygon is Basic-locked, polygonize isn't).

Detected borders often stop a metre or two short of the border they should meet, so nothing
closes. Two repairs, then polygonize:
  1. JOIN: two loose ends (dangles) facing each other within `join_d` metres are welded at
     their midpoint (also closes a near-complete ring of the same line).
  2. EXTEND: a remaining loose end is projected straight ahead up to `ext_d` metres; if that ray
     hits another line, the end is extended exactly onto it.
Then unary_union nodes everything, shapely.polygonize builds faces, and polygons smaller than
`min_area` m2 (slivers) are dropped. The lines file is NOT modified -- polygons are a new output
with `area_m2` / `perim_m` fields, ready for QA in Pro.

Usage: python close_parcels.py in.shp out_parcels.shp [--ext 5] [--join 3] [--min-area 30]
"""
import os, argparse
import numpy as np
from osgeo import ogr
from shapely import wkt, STRtree
from shapely.geometry import LineString
from shapely.ops import unary_union, polygonize
ogr.UseExceptions()

RND = 2                     # endpoint welding key rounding (1 cm)


def _key(p):
    return (round(p[0], RND), round(p[1], RND))


def _degrees(lines):
    deg = {}
    for cs in lines:
        for e in (0, -1):
            deg[_key(cs[e])] = deg.get(_key(cs[e]), 0) + 1
    return deg


def _append(cs, e, pt):
    if e == -1:
        cs.append(tuple(pt))
    else:
        cs.insert(0, tuple(pt))


def _end_dir(cs, e):
    """Unit vector pointing OUT of the line at end e (0 = start, -1 = end)."""
    p = np.array(cs[e], float); q = None
    idx = range(len(cs) - 2, -1, -1) if e == -1 else range(1, len(cs))
    for k in idx:
        q = np.array(cs[k], float)
        if np.hypot(*(p - q)) > 0.05:
            break
    if q is None:
        return None
    v = p - q; n = np.hypot(*v)
    return v / n if n > 1e-9 else None


def join_dangles(lines, join_d):
    """Weld pairs of facing loose ends within join_d at their midpoint (in place)."""
    deg = _degrees(lines)
    ends = [[i, e, np.array(cs[e], float)] for i, cs in enumerate(lines) for e in (0, -1)
            if deg[_key(cs[e])] == 1]
    used, n = set(), 0
    for a in range(len(ends)):
        if a in used:
            continue
        i, e, p = ends[a]
        best, bd = None, join_d
        for b in range(a + 1, len(ends)):
            if b in used:
                continue
            j, f, q = ends[b]
            if j == i and LineString(lines[i]).length <= 4 * join_d:
                continue                              # don't fold a short line onto itself
            d = float(np.hypot(*(p - q)))
            if d < bd:
                bd, best = d, b
        if best is None:
            continue
        j, f, q = ends[best]
        m = ((p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0)
        _append(lines[i], e, m); _append(lines[j], f, m)
        used.add(a); used.add(best); n += 1
    return n


def extend_dangles(lines, ext_d):
    """Extend each remaining loose end straight ahead onto the first line it hits (in place)."""
    geoms = [LineString(cs) for cs in lines]
    tree = STRtree(geoms)
    deg = _degrees(lines)
    n = 0
    for i, cs in enumerate(lines):
        for e in (0, -1):
            if deg.get(_key(cs[e])) != 1:
                continue
            v = _end_dir(cs, e)
            if v is None:
                continue
            p = np.array(cs[e], float)
            ray = LineString([tuple(p + v * 0.03), tuple(p + v * ext_d)])
            hit, hd = None, ext_d + 1.0
            for c in tree.query(ray):
                j = int(c)
                if j == i:
                    continue
                x = ray.intersection(geoms[j])
                if x.is_empty:
                    continue
                for g in getattr(x, "geoms", [x]):
                    q = np.array(g.coords[0], float)
                    d = float(np.hypot(*(q - p)))
                    if 0.02 < d < hd:
                        hd, hit = d, q
            if hit is not None:
                _append(cs, e, hit); n += 1
    return n


def close(lines, ext_d=5.0, join_d=3.0):
    """Pure core: lines = [[(x,y),...], ...] coord lists, repaired IN PLACE.
    Returns (n_dangles_before, n_joined_pairs, n_extended)."""
    nd = sum(1 for v in _degrees(lines).values() if v == 1)
    nj = join_dangles(lines, join_d)
    ne = extend_dangles(lines, ext_d)
    return nd, nj, ne


def run(src, out, ext_d=5.0, join_d=3.0, min_area=30.0):
    dl = ogr.Open(src); ll = dl.GetLayer(); sr = ll.GetSpatialRef()
    lines = [list(wkt.loads(f.GetGeometryRef().ExportToWkt()).coords)
             for f in ll if f.GetGeometryRef()]
    nd, nj, ne = close(lines, ext_d, join_d)

    net = unary_union([LineString(cs) for cs in lines if len(cs) >= 2])
    polys = list(polygonize(net))
    kept = [p for p in polys if p.area >= min_area]

    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out):
        drv.DeleteDataSource(out)
    dso = drv.CreateDataSource(out); lyo = dso.CreateLayer("parcels", sr, ogr.wkbPolygon)
    lyo.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    fa = ogr.FieldDefn("area_m2", ogr.OFTReal); fa.SetWidth(18); fa.SetPrecision(1); lyo.CreateField(fa)
    fp = ogr.FieldDefn("perim_m", ogr.OFTReal); fp.SetWidth(18); fp.SetPrecision(1); lyo.CreateField(fp)
    defn = lyo.GetLayerDefn()
    for k, p in enumerate(sorted(kept, key=lambda g: -g.area)):
        g = ogr.CreateGeometryFromWkt(p.wkt)
        f = ogr.Feature(defn); f.SetField("id", k)
        f.SetField("area_m2", float(p.area)); f.SetField("perim_m", float(p.length))
        f.SetGeometry(g); lyo.CreateFeature(f); f = None
    dso = None
    med = float(np.median([p.area for p in kept])) if kept else 0.0
    print("[parcels] %d lines, %d loose ends: joined %d pairs, extended %d; polygonize -> "
          "%d faces, kept %d >= %.0f m2 (median %.0f m2) -> %s"
          % (len(lines), nd, nj, ne, len(polys), len(kept), min_area, med, out))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("out")
    ap.add_argument("--ext", type=float, default=5.0, dest="ext_d",
                    help="max straight-ahead extension (m) of a loose end onto a nearby line")
    ap.add_argument("--join", type=float, default=3.0, dest="join_d",
                    help="max gap (m) between two facing loose ends to weld them")
    ap.add_argument("--min-area", type=float, default=30.0, dest="min_area",
                    help="drop polygons smaller than this (m2) as slivers")
    a = ap.parse_args()
    run(a.src, a.out, a.ext_d, a.join_d, a.min_area)
