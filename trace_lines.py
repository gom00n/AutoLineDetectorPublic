"""Vectorize a 1-px skeleton to TOPOLOGICALLY CLEAN polylines by tracing the skeleton graph
(base arcgispro-py3 python: numpy + scipy + gdal/ogr).

Why not arcpy RasterToPolyline: on this data it shatters a connected skeleton into thousands
of ~4-px features whose endpoints are offset ~0.5-1 m and never coincide (0 degree-2 nodes),
so the lines can't be merged or closed into parcels. Here each arc runs junction-to-junction
and every junction is collapsed to ONE shared coordinate -> exact topology, naturally long
lines, and pixel-precise hair / loop pruning.

Method (vectorised; only small per-arc walks happen in Python):
  deg     = 8-neighbour count of each skeleton pixel
  junc    = deg >= 3            (collapsed to one centroid node per junction blob)
  arcs    = skeleton minus junctions, labelled into connected components (one per arc)
  each arc component is ordered into a path; its ends are snapped to the adjacent junction
  centroid (shared node) or kept as a free end.

Pruning:
  - loop : arc whose two ends meet the SAME junction (or start==end) and is < MIN_LOOP px.
  - spur : arc with one free end and one junction end, shorter than MIN_SPUR px ("hair").
  - standalone arcs (both ends free) are always kept.

Usage: python trace_lines.py in_skel.tif out_lines.shp [min_spur=12] [min_loop=30]
"""
import os, re, sys, numpy as np
from osgeo import gdal, ogr, osr
from scipy import ndimage
gdal.UseExceptions()


def source_crs(skel_path):
    """The scan is georeferenced by a polynomial/GCP transform (.tfwx), so GDAL surfaces no
    projection -- the CRS lives in the source <stem>.tif.aux.xml ESRI sidecar. Pull its WKID
    (e.g. 6991 = Israeli Grid 05/12) so the output shapefile carries the same CRS as the scan.
    Returns an osr.SpatialReference or None. Layout: skel is ROOT/<stem>/<stem>_skel_*.tif."""
    folder = os.path.dirname(skel_path)
    root, stem = os.path.dirname(folder), os.path.basename(folder)
    aux = os.path.join(root, stem + ".tif.aux.xml")
    if not os.path.exists(aux):
        return None
    with open(aux, encoding="utf-8", errors="ignore") as fh:
        text = fh.read()
    m = re.search(r"<LatestWKID>(\d+)</LatestWKID>", text) or re.search(r"<WKID>(\d+)</WKID>", text)
    sr = osr.SpatialReference()
    if m and sr.ImportFromEPSG(int(m.group(1))) == 0:
        return sr
    w = re.search(r"<WKT>(.*?)</WKT>", text, re.S)               # fall back to the ESRI WKT
    if w and sr.SetFromUserInput(w.group(1).strip()) == 0:
        return sr
    return None

NB8 = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
_3 = np.ones((3, 3), np.uint8)


def _rdp(pts, eps):
    """Ramer-Douglas-Peucker line simplification (removes the 1-px staircase). eps in pixels."""
    pts = np.asarray(pts, float)
    if len(pts) < 3 or eps <= 0:
        return pts
    keep = np.zeros(len(pts), bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        seg = pts[j] - pts[i]
        L = np.hypot(seg[0], seg[1])
        rel = pts[i + 1:j] - pts[i]
        if L == 0:
            dist = np.hypot(rel[:, 0], rel[:, 1])
        else:
            dist = np.abs(seg[0] * rel[:, 1] - seg[1] * rel[:, 0]) / L   # perpendicular distance
        k = int(np.argmax(dist))
        if dist[k] > eps:
            idx = i + 1 + k
            keep[idx] = True
            stack.append((i, idx))
            stack.append((idx, j))
    return pts[keep]


def _order(coords):
    """Order an arc component's pixels (a simple path or small loop) into a sequence."""
    cset = set(map(tuple, coords))

    def nb(p):
        r, c = p
        return [(r + dr, c + dc) for dr, dc in NB8 if (r + dr, c + dc) in cset]

    ends = [p for p in cset if len(nb(p)) == 1]
    start = ends[0] if ends else next(iter(cset))
    path = [start]
    used = {start}
    prev, cur = None, start
    while True:
        nxts = [q for q in nb(cur) if q != prev and q not in used]
        if not nxts:
            break
        nb_ = nxts[0]
        path.append(nb_)
        used.add(nb_)
        prev, cur = cur, nb_
    return path


def _trace(sk, min_spur, min_loop, simplify_px):
    skb = sk.astype(np.uint8)
    deg = (ndimage.convolve(skb, _3, mode="constant") - skb) * sk     # 8-neighbour count
    junc = sk & (deg >= 3)
    jlab, jn = ndimage.label(junc, structure=_3)
    jcent = ndimage.center_of_mass(junc, jlab, range(1, jn + 1)) if jn else []
    jcent = {i + 1: jcent[i] for i in range(jn)}                      # label -> (r,c) centroid

    arcs_mask = sk & ~junc
    alab, an = ndimage.label(arcs_mask, structure=_3)
    objs = ndimage.find_objects(alab)

    def junc_labels_near(p):
        """Distinct junction labels touching pixel p, in neighbour order."""
        r, c = p
        out = []
        for dr, dc in NB8:
            jl = jlab[r + dr, c + dc]
            if jl and jl not in out:
                out.append(jl)
        return out

    out_arcs, n_spur, n_loop = [], 0, 0
    for cid in range(1, an + 1):
        sl = objs[cid - 1]
        if sl is None:
            continue
        sub = alab[sl] == cid
        coords = np.argwhere(sub) + np.array([sl[0].start, sl[1].start])
        path = _order([tuple(p) for p in coords])
        plen = len(path)
        a, b = path[0], path[-1]
        la, lb = junc_labels_near(a), junc_labels_near(b)
        ja, jb = (la[0] if la else 0), (lb[0] if lb else 0)     # 0 = free end
        if ja and ja == jb:
            # A very short arc can touch TWO junction blobs with the same pixel (always, when
            # it is 1 px long: a == b). Taking the first label for both ends made it look
            # like a bubble off one junction, so it was dropped as a loop -- cutting a 3-4 px
            # hole in a border that runs through two T-junctions a few pixels apart. Prefer
            # the reading in which the arc LINKS two different junctions.
            other = [l for l in lb if l != ja]
            if other:
                jb = other[0]
            else:
                other = [l for l in la if l != jb]
                if other:
                    ja = other[0]
        free_a, free_b = (ja == 0 and deg[a] <= 1), (jb == 0 and deg[b] <= 1)

        if ja and jb and ja == jb and plen < min_loop:          # bubble off one junction
            n_loop += 1
            continue
        if (free_a ^ free_b) and (ja or jb) and plen < min_spur:  # hair off a junction
            n_spur += 1
            continue

        spath = _rdp(path, simplify_px)                         # straighten (kill pixel staircase)
        verts = [jcent[ja]] if ja else []                       # snap end to shared node
        verts += [(float(r), float(c)) for r, c in spath]
        if jb:
            verts.append(jcent[jb])
        out_arcs.append(verts)
    return out_arcs, n_spur, n_loop, an


def run(src, out, min_spur=12, min_loop=30, simplify_px=1.5):
    ds = gdal.Open(src)
    sk0 = ds.GetRasterBand(1).ReadAsArray() > 0
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    sk = np.pad(sk0, 1)                                          # border so neighbours stay in-bounds

    arcs, n_spur, n_loop, an = _trace(sk, min_spur, min_loop, simplify_px)

    def to_xy(rc):
        r, c = rc[0] - 1, rc[1] - 1                             # undo pad
        x = gt[0] + (c + 0.5) * gt[1] + (r + 0.5) * gt[2]
        y = gt[3] + (c + 0.5) * gt[4] + (r + 0.5) * gt[5]
        return x, y

    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out):
        drv.DeleteDataSource(out)
    sr = None                       # GDAL surfaces no projection (georef is world-file/GCP only)
    if proj:                        # -> recover the scan's CRS from its .tif.aux.xml ESRI sidecar
        sr = osr.SpatialReference()
        sr.ImportFromWkt(proj)
    else:
        sr = source_crs(src)
    print("[trace] CRS: %s" % (sr.GetName() if sr else "unknown (no .prj)"))
    dso = drv.CreateDataSource(out)
    lyr = dso.CreateLayer("lines", sr, ogr.wkbLineString)
    lyr.CreateField(ogr.FieldDefn("arc", ogr.OFTInteger))
    defn = lyr.GetLayerDefn()
    for i, verts in enumerate(arcs):
        if len(verts) < 2:
            continue
        g = ogr.Geometry(ogr.wkbLineString)
        for rc in verts:
            g.AddPoint_2D(*to_xy(rc))
        f = ogr.Feature(defn)
        f.SetField("arc", i)
        f.SetGeometry(g)
        lyr.CreateFeature(f)
        f = None
    dso = None
    print("[trace] %d arc components; dropped %d spurs(<%dpx) + %d loops(<%dpx) -> %d features"
          % (an, n_spur, min_spur, n_loop, min_loop, len(arcs)))
    print("[trace] wrote %s" % out)


if __name__ == "__main__":
    a = sys.argv
    run(a[1], a[2],
        int(a[3]) if len(a) > 3 else 12,
        int(a[4]) if len(a) > 4 else 30,
        float(a[5]) if len(a) > 5 else 1.5)
