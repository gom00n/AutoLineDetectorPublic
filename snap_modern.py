"""Conflate detected historical lines onto reference boundaries -- per-point network snapping.

The reference is ONE noded network: polygon-boundary rings pooled from every --snap-src source
(e.g. the project's block/Gush polygons Poly_ID/Poly_Map/Poly_map2 PLUS the cadastral parcels
PARCEL_ALL from a separate parcels GDB), unioned so shared parcel edges are deduplicated and
junctions are noded. Do NOT include the Motzag sheet-index polygons -- those just define the
scan footprint.

Per line: densify at 1 m and project every point onto the nearest network edge within `tol`
that runs roughly PARALLEL to the line there (perpendicular parcel-dividers are ignored, so
crossing the grid never makes teeth). Contiguous stretches that follow the network for at
least `min_run` metres are replaced by their projected feet; the rest keeps the ink geometry.
A line is SPLIT into parts at the transitions, so one line can be half-snapped, half-historical
-- no more all-or-nothing per line. Lines shorter than min_run stay all-or-nothing at
match_frac. The `snapped` field is 0 = kept as ink, 1 = snapped to a reference boundary.

Sources are read 3D/measured + true-curve safe, spatially filtered to the lines' extent, and
reprojected to the detected CRS. A source can be a File GDB or a shapefile (ogr.Open takes both;
for a shapefile the layer name is its basename without .shp).

Usage (base arcgispro-py3 python: gdal/ogr + shapely + numpy):
  python snap_modern.py in.shp out.shp [--tol 3] [--match-frac 0.6] [--min-run 10]
         [--snap-src GDB1 "Poly_ID,Poly_Map"] [--snap-src GDB2 "PARCEL_ALL"]   (repeatable)
"""
import os, sys, argparse
import numpy as np
from osgeo import ogr, osr
from shapely import wkt, STRtree
from shapely.geometry import LineString, Point, box
from shapely.ops import substring, unary_union, linemerge
ogr.UseExceptions(); osr.UseExceptions()

STEP = 1.0        # densify spacing along a detected line (m)
COS_MIN = 0.5     # reject reference edges running > ~60 deg off the line (no teeth)
GAP_FILL = 4.0    # bridge off-network gaps shorter than this between two snapped runs (m)
SIMP = 0.3        # simplify tolerance for the projected-feet polyline (m)
# run-quality guards (added after the red-ink M14 review, 2026-07-07):
KINK_MAX = 1.0       # a run may not start/end mid-line with a sideways jump bigger than this (m)
STRAIGHT_COS = 0.94  # ...when the ink runs STRAIGHT through the transition (~20 deg); a real
#                      corner in the ink legitimises the jump, a straight line does not.
WIGGLE = 1.2         # max feet-polyline length vs the ink stretch it replaces: oscillating
#                      between near-parallel reference edges (saw teeth) adds length; reject.


def _ct(src_sr, dst_sr):
    if src_sr is None:
        return None
    s = src_sr.Clone(); s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    d = dst_sr.Clone(); d.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return osr.CoordinateTransformation(s, d) if not s.IsSame(d) else None


def _boundary_rings(geom_iter, ct):
    rings = []
    for g in geom_iter:
        if g is None:
            continue
        g = g.Clone()
        g.FlattenTo2D()                          # drop Z/M (GDB 3D-measured parcels) so shapely parses
        if g.HasCurveGeometry():                 # GDB true-curve arcs -> stroke to line segments
            g = g.GetLinearGeometry()            # (shapely can't parse CIRCULARSTRING/COMPOUNDCURVE)
        if ct is not None:
            g.Transform(ct)
        b = g.Boundary()
        if b is None:
            continue
        if b.HasCurveGeometry():                 # boundary of a curve polygon can still be curved
            b = b.GetLinearGeometry()
        s = wkt.loads(b.ExportToWkt())
        for ls in getattr(s, "geoms", [s]):
            if ls.length > 1.0:
                rings.append(ls)
    return rings


def _src_rings(src, layers, dst_sr, bbox, pad=300.0):
    """Boundary rings of ALL polygons in the chosen layers of one OGR source (File GDB or
    shapefile), spatially clipped to the lines' bbox, reprojected to dst_sr."""
    ds = ogr.Open(src)
    if ds is None:
        print("[snap] WARNING: cannot open snap source %s" % src); return []
    rings = []
    for name in layers:
        lyr = ds.GetLayerByName(name)
        if lyr is None:
            print("[snap] WARNING: snap layer '%s' not found in %s" % (name, src)); continue
        lsr = lyr.GetSpatialRef()
        ctf = _ct(dst_sr, lsr)                                  # bbox (line CRS) -> layer CRS to filter
        cs = [(bbox[0], bbox[2]), (bbox[1], bbox[2]), (bbox[1], bbox[3]), (bbox[0], bbox[3])]
        tc = [ctf.TransformPoint(x, y)[:2] for x, y in cs] if ctf else cs
        xs = [c[0] for c in tc]; ys = [c[1] for c in tc]
        lyr.SetSpatialFilterRect(min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad)
        rings += _boundary_rings((f.GetGeometryRef() for f in lyr), _ct(lsr, dst_sr))
        lyr.SetSpatialFilter(None)
    return rings


def build_network(rings):
    """Node + dedupe the pooled boundary rings into junction-to-junction edges + an STRtree.
    unary_union removes the doubled shared parcel edges and nodes all crossings; linemerge
    chains the pieces back up between junctions."""
    lines = []
    if rings:
        net = unary_union(rings)
        lines = [g for g in getattr(net, "geoms", [net]) if g.geom_type == "LineString"]
    if not lines:
        return [], None
    merged = linemerge(lines) if len(lines) > 1 else lines[0]
    edges = [g for g in getattr(merged, "geoms", [merged]) if g.length > 0.1]
    return edges, (STRtree(edges) if edges else None)


def _densify(line, step):
    n = max(1, int(line.length // step))
    return [line.interpolate(k / n, normalized=True) for k in range(n + 1)]


_seg_cache = {}       # id(edge) -> (coords array, cumulative segment lengths) for _local_tangent


def _local_tangent(g, s):
    """Unit direction of the EXACT segment of edge g containing arc-length s. An interpolated
    +-EPS window would smear across corners (linemerge chains L-shaped edges through degree-2
    corner nodes), yielding a diagonal pseudo-tangent that lets points near a parcel corner
    snap onto the perpendicular arm -- the exact segment cannot."""
    cached = _seg_cache.get(id(g))
    if cached is None:
        xy = np.asarray(g.coords, float)
        cum = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(xy, axis=0).T))])
        _seg_cache[id(g)] = cached = (xy, cum)
    xy, cum = cached
    i = min(max(int(np.searchsorted(cum, s, side="right")) - 1, 0), len(xy) - 2)
    v = xy[i + 1] - xy[i]; nv = np.hypot(*v)
    return v / nv if nv > 1e-9 else None


def _nearest_foot(p, edges, tree, tol, tangent=None):
    """Foot of p on the nearest network edge within tol; if `tangent` (unit vector) is given,
    edges running > ~60 deg off it AT THE FOOT are skipped (that's what keeps perpendicular
    dividers from making teeth). Returns None when nothing qualifies."""
    best, bestd = None, tol + 1.0
    for c in tree.query(box(p.x - tol, p.y - tol, p.x + tol, p.y + tol)):
        g = edges[int(c)]
        d = p.distance(g)
        if d > tol or d >= bestd:
            continue
        s = g.project(p)
        if tangent is not None:
            dr = _local_tangent(g, s)
            if dr is not None and abs(float(np.dot(tangent, dr))) < COS_MIN:
                continue
        bestd, best = d, g.interpolate(s)
    return best


def _feet(L, edges, tree, tol, step):
    """Densified points of L + per point its projected foot on the network (or None)."""
    pts = _densify(L, step)
    n = len(pts)
    feet = [None] * n
    for i, p in enumerate(pts):
        a, b = pts[max(0, i - 1)], pts[min(n - 1, i + 1)]
        dl = np.array([b.x - a.x, b.y - a.y]); nl = np.hypot(*dl)
        feet[i] = _nearest_foot(p, edges, tree, tol, dl / nl if nl > 1e-6 else None)
    return pts, feet


def _runs(flags):
    """[(i0, i1_exclusive, flag), ...] maximal same-flag runs, as mutable lists."""
    runs, i = [], 0
    while i < len(flags):
        j = i
        while j < len(flags) and flags[j] == flags[i]:
            j += 1
        runs.append([i, j, flags[i]]); i = j
    return runs


def _merge_runs(runs):
    out = []
    for r in runs:
        if out and out[-1][2] == r[2]:
            out[-1][1] = r[1]
        else:
            out.append(r)
    return out


def _ink_tangent(pts, i, j):
    """Unit vector pts[i] -> pts[j] (indices clamped); None when degenerate."""
    n = len(pts)
    a, b = pts[max(0, min(i, n - 1))], pts[max(0, min(j, n - 1))]
    v = np.array([b.x - a.x, b.y - a.y]); nv = np.hypot(*v)
    return v / nv if nv > 1e-9 else None


def _run_feet(i0, i1, pts, feet, edges, tree, tol):
    """Projected-feet coordinates of one snapped run (bridged gap points fall back to the
    nearest edge with no tangent test, and to the ink point when even that fails)."""
    coords = []
    for k in range(i0, i1):
        f = feet[k] or _nearest_foot(pts[k], edges, tree, tol)
        coords.append((f.x, f.y) if f else (pts[k].x, pts[k].y))
    return coords


def _straight_through(pts, ib, w):
    """True when the ink runs straight through point index ib (~20 deg over a +-w window)."""
    t0 = _ink_tangent(pts, ib - w, ib)
    t1 = _ink_tangent(pts, ib, ib + w)
    return t0 is not None and t1 is not None and float(np.dot(t0, t1)) >= STRAIGHT_COS


def snap_line(L, edges, tree, tol, match_frac=0.6, min_run=10.0, step=STEP):
    """Conflate one detected line onto the reference network.
    Returns [(LineString, snapped), ...] parts covering L in order (snapped: 1/0).
    Stretches following the network within tol for >= min_run metres become the projected feet;
    the rest keeps the ink geometry (original vertices, via substring). Lines shorter than
    min_run stay all-or-nothing at match_frac. Runs that would saw or kink (see _bad_run) are
    kept as ink."""
    if not edges or tree is None or L.length < 1.0:
        return [(L, 0)]
    pts, feet = _feet(L, edges, tree, tol, step)
    n = len(pts)
    if not any(f is not None for f in feet):
        return [(L, 0)]
    seg = L.length / (n - 1) if n > 1 else L.length          # metres between densified points
    flags = [f is not None for f in feet]

    if L.length < min_run:                                   # short line: all-or-nothing
        if sum(flags) / n < match_frac:
            return [(L, 0)]
        coords = [((f.x, f.y) if f else (p.x, p.y)) for p, f in zip(pts, feet)]
        line = LineString(coords).simplify(SIMP)
        if line.length <= 0.5 or line.length > WIGGLE * L.length:      # saw guard
            return [(L, 0)]
        return [(line, 1)]

    runs = _runs(flags)
    for k, r in enumerate(runs):                             # bridge junction jitter between snapped runs
        if (not r[2]) and 0 < k < len(runs) - 1 and (r[1] - r[0]) * seg <= GAP_FILL:
            r[2] = True
    runs = _merge_runs(runs)
    for r in runs:                                           # too short to trust: keep as ink
        if r[2] and (r[1] - r[0]) * seg < min_run:
            r[2] = False
    runs = _merge_runs(runs)

    # --- run-quality guards (saw + kink; see the constants up top) ---------------------------
    # A kinking run is TRIMMED back to where its offset drops below KINK_MAX: a genuine run
    # that merely overhangs a corner loses its tip, while a false parallel-grab (the reference
    # is offset along the WHOLE run -- red-ink pic 1) trims away completely. Saw runs (feet much
    # longer than the ink -- oscillation between near-parallel edges, pic 2) are dropped whole.
    for i0, i1, fl in runs:                                  # sync flags with gap-filled runs
        for k in range(i0, i1):
            flags[k] = fl

    def _off(k):
        return pts[k].distance(feet[k]) if feet[k] is not None else 0.0

    w = max(2, int(round(3.0 / max(seg, 1e-6))))             # ~3 m tangent window in points
    for _ in range(4):                                       # trims/flips converge fast
        changed = False
        for i0, i1, fl in runs:
            if not fl:
                continue
            coords = _run_feet(i0, i1, pts, feet, edges, tree, tol)
            ink_len = max(seg * (i1 - 1 - i0), 1e-6)
            if len(coords) >= 2 and LineString(coords).length > WIGGLE * ink_len:
                for k in range(i0, i1):                      # saw: drop the run entirely
                    flags[k] = False
                changed = True
                continue
            j0, j1 = i0, i1                                  # kink: trim (line ends are welded,
            if i0 > 0 and _off(i0) > KINK_MAX and _straight_through(pts, i0, w):        # not kinked)
                while j0 < i1 and _off(j0) > KINK_MAX:
                    j0 += 1
            if i1 < n and _off(i1 - 1) > KINK_MAX and _straight_through(pts, i1 - 1, w):
                while j1 > j0 and _off(j1 - 1) > KINK_MAX:
                    j1 -= 1
            if (j1 - j0) * seg < min_run:                    # too little left after trimming
                j0 = j1 = i0
            if (j0, j1) != (i0, i1):
                for k in list(range(i0, j0)) + list(range(j1, i1)):
                    flags[k] = False
                changed = True
        if not changed:
            break
        runs = _merge_runs(_runs(flags))

    snapped_coords = []                                      # snapped runs first, so ink runs can
    for i0, i1, fl in runs:                                  # share their exact transition coords
        if not fl:
            snapped_coords.append(None); continue
        coords = _run_feet(i0, i1, pts, feet, edges, tree, tol)
        line = LineString(coords).simplify(SIMP) if len(coords) >= 2 else None
        snapped_coords.append(list(line.coords) if line is not None else None)

    parts = []
    for idx, (i0, i1, fl) in enumerate(runs):
        if fl:
            coords = snapped_coords[idx]
        else:
            d0 = 0.0 if i0 == 0 else i0 * seg                # ink: substring keeps original vertices
            d1 = L.length if i1 >= n else (i1 - 1) * seg
            piece = substring(L, d0, d1)
            coords = [] if piece.is_empty else list(piece.coords)
            if idx > 0 and snapped_coords[idx - 1]:          # share the transition coordinate
                coords = [snapped_coords[idx - 1][-1]] + coords
            if idx < len(runs) - 1 and snapped_coords[idx + 1]:
                coords = coords + [snapped_coords[idx + 1][0]]
        if coords and len(coords) >= 2:
            line = LineString(coords)
            if line.length > 0.05:
                parts.append((line, int(fl)))
    return parts if parts else [(L, 0)]


def run(src, out, tol=3.0, match_frac=0.6, min_run=10.0, sources=None):
    """sources: list of (path, [layer, ...]) -- polygon layers to snap onto, pooled across
    several GDBs/shapefiles (e.g. the project's Poly_ID/Poly_Map + a separate parcels GDB)."""
    dl = ogr.Open(src); ll = dl.GetLayer(); sr = ll.GetSpatialRef(); ext = ll.GetExtent()
    feats = [(f.GetField("arc") if f.GetFieldIndex("arc") >= 0 else f.GetFID(),
              wkt.loads(f.GetGeometryRef().ExportToWkt())) for f in ll]

    rings = []                                            # one combined reference-boundary pool
    for path, layers in (sources or []):
        if path and layers:
            rings += _src_rings(path, layers, sr, ext)
    edges, tree = build_network(rings)
    print("[snap] %d lines; %d boundary rings -> %d network edges; tol=%.2f match_frac=%.2f min_run=%.1f"
          % (len(feats), len(rings), len(edges), tol, match_frac, min_run))

    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out):
        drv.DeleteDataSource(out)
    dso = drv.CreateDataSource(out); lyo = dso.CreateLayer("lines", sr, ogr.wkbLineString)
    lyo.CreateField(ogr.FieldDefn("arc", ogr.OFTInteger))
    lyo.CreateField(ogr.FieldDefn("snapped", ogr.OFTInteger)); defn = lyo.GetLayerDefn()

    results, ends = [], []                                # per line: [[coords, flag], ...] parts
    key = lambda xy: (round(xy[0], 1), round(xy[1], 1))
    for arc, L in feats:
        parts = snap_line(L, edges, tree, tol, match_frac, min_run)
        results.append((arc, [[list(P.coords), fl] for P, fl in parts]))
        ends.append((key(L.coords[0]), key(L.coords[-1])))

    node = {}                                             # weld shared junctions (original line ends)
    for i, (a, b) in enumerate(ends):
        node.setdefault(a, []).append((i, 0)); node.setdefault(b, []).append((i, -1))
    n_weld = 0
    for members in node.values():
        if len(members) < 2:
            continue
        # snapped ends sit ON a reference boundary -- they win; ink ends are pulled onto them.
        # Only an all-ink junction is welded to the plain mean.
        onb = [results[i][1][e][0][e] for i, e in members if results[i][1][e][1]]
        tgt = np.mean(onb if onb else [results[i][1][e][0][e] for i, e in members], axis=0)
        tgt = (float(tgt[0]), float(tgt[1]))
        for i, e in members:
            coords = results[i][1][e][0]                  # part e (0=first / -1=last), coord e
            if tuple(coords[e]) != tgt:
                coords[e] = tgt; n_weld += 1

    n_full = n_part = n_kept = 0
    len_snap = len_tot = 0.0
    for (arc, parts) in results:
        fls = [fl for _, fl in parts]
        if all(fls):
            n_full += 1
        elif any(fls):
            n_part += 1
        else:
            n_kept += 1
        for coords, fl in parts:
            if len(coords) < 2:
                continue
            g = ogr.Geometry(ogr.wkbLineString)
            for x, y in coords:
                g.AddPoint_2D(float(x), float(y))
            ln = g.Length(); len_tot += ln; len_snap += ln if fl else 0.0
            f = ogr.Feature(defn); f.SetField("arc", int(arc))
            f.SetField("snapped", fl); f.SetGeometry(g); lyo.CreateFeature(f); f = None
    dso = None
    print("[snap] %d lines: %d fully snapped, %d partially, %d kept as ink; %.0f%% of length "
          "on a reference boundary; welded %d ends -> %s"
          % (len(results), n_full, n_part, n_kept,
             (100.0 * len_snap / len_tot) if len_tot else 0.0, n_weld, out))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("out")
    ap.add_argument("--tol", type=float, default=3.0, help="max distance (m) line->boundary to snap")
    ap.add_argument("--match-frac", type=float, default=0.6, dest="match_frac",
                    help="short lines (< min-run): fraction that must lie within tol to snap whole")
    ap.add_argument("--min-run", type=float, default=10.0, dest="min_run",
                    help="min stretch (m) that must follow a boundary before that stretch snaps")
    ap.add_argument("--snap-src", nargs=2, action="append", default=[], metavar=("SRC", "LAYERS"),
                    dest="snap_src", help="a GDB or shapefile + comma-separated polygon layer(s) "
                    "to snap onto; repeatable (all sources are pooled into one network)")
    a = ap.parse_args()
    sources = [(g, [x.strip() for x in lyr.split(",") if x.strip()]) for g, lyr in a.snap_src]
    run(a.src, a.out, a.tol, a.match_frac, a.min_run, sources)
