"""Drop disconnected junk by connected-component size (base arcgispro-py3 python: ogr+numpy).

Map legend, the marginal number/area-table, and stray symbols survive detection as SHORT,
ISOLATED line clusters, whereas real parcel borders form large CONNECTED networks (features
sharing endpoint nodes). This keeps only features whose connected component's total length is
>= MIN_COMP_LEN, removing most junk automatically. On M15 (vs the user's hand-cleaned truth)
a 20 m threshold dropped 96% of the margin table and 85% of other junk while keeping 59/64
real borders -- i.e. it turns ~333 manual deletions into a couple dozen.

Tune MIN_COMP_LEN per sheet if real borders are sparse/fragmented (lower it) or junk survives
(raise it). Run AFTER merge_lines (shared nodes must exist).

Also straightens the kept lines: a Douglas-Peucker pass (ogr.SimplifyPreserveTopology, free on
Basic -- arcpy SimplifyLine is license-blocked) removes the 1-px staircase without breaking
topology, so borders read as clean straight segments between corners.

Usage: python filter_components.py in_lines.shp out_lines.shp [min_comp_len_m=20] [simplify_m=0.5]
"""
import os, sys
from collections import defaultdict
from osgeo import ogr, osr
ogr.UseExceptions()


def _components(geoms, prec=1):
    """Union-find features that share an endpoint (rounded to 10^-prec m). Returns parent list."""
    n = len(geoms)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    node2 = defaultdict(list)
    for i, g in enumerate(geoms):
        m = g.GetPointCount()
        for nd in ((round(g.GetX(0), prec), round(g.GetY(0), prec)),
                   (round(g.GetX(m - 1), prec), round(g.GetY(m - 1), prec))):
            node2[nd].append(i)
    for fs in node2.values():
        r0 = find(fs[0])
        for j in fs[1:]:
            parent[find(j)] = r0
    return [find(i) for i in range(n)]


def run(src, out, min_comp_len=20.0, simplify_m=0.5):
    ds = ogr.Open(src)
    lyr = ds.GetLayer()
    sr = lyr.GetSpatialRef()
    geoms = [f.GetGeometryRef().Clone() for f in lyr]
    n = len(geoms)
    roots = _components(geoms)
    comp_len = defaultdict(float)
    for i, g in enumerate(geoms):
        comp_len[roots[i]] += g.Length()
    keep = [i for i in range(n) if comp_len[roots[i]] >= min_comp_len]

    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out):
        drv.DeleteDataSource(out)
    dso = drv.CreateDataSource(out)
    lyo = dso.CreateLayer("lines", sr, ogr.wkbLineString)
    lyo.CreateField(ogr.FieldDefn("arc", ogr.OFTInteger))
    defn = lyo.GetLayerDefn()
    for k, i in enumerate(keep):
        g = geoms[i]
        if simplify_m and simplify_m > 0:                 # straighten: drop the pixel staircase
            g = g.SimplifyPreserveTopology(simplify_m)
        f = ogr.Feature(defn)
        f.SetField("arc", k)
        f.SetGeometry(g)
        lyo.CreateFeature(f)
        f = None
    dso = None
    ncomp = len(set(roots))
    nkeep_c = len({roots[i] for i in keep})
    print("[filter] %d feats / %d components -> kept %d feats / %d components (>= %.0fm)"
          % (n, ncomp, len(keep), nkeep_c, min_comp_len))
    print("[filter] wrote %s" % out)


if __name__ == "__main__":
    a = sys.argv
    run(a[1], a[2],
        float(a[3]) if len(a) > 3 else 20.0,
        float(a[4]) if len(a) > 4 else 0.5)
