"""Remove 'hair' barbs from the final polylines (base arcgispro-py3 python: ogr + numpy).

After merge/filter a few short branches still hang off the border network -- a feature with one
FREE end (a node touched by no other line) and one JUNCTION end (a node shared by >=2 others),
shorter than MAX_HAIR. Those are skeleton hair, not real borders. Long dangles (real dead-end
borders) and through-segments (both ends at junctions) are kept. Iterated, since removing one
barb can expose the next.

Usage: python dehair_lines.py in.shp out.shp [max_hair_m=4] [passes=3]
"""
import os, sys
from collections import Counter
from osgeo import ogr
ogr.UseExceptions()


def _key(xy):
    return (round(xy[0], 1), round(xy[1], 1))


def run(src, out, max_hair=4.0, passes=3):
    ds = ogr.Open(src); lyr = ds.GetLayer(); sr = lyr.GetSpatialRef()
    feats = [(f.GetField("arc") if f.GetFieldIndex("arc") >= 0 else f.GetFID(),
              f.GetGeometryRef().Clone()) for f in lyr]
    alive = [True] * len(feats)
    removed_total = 0
    for _ in range(passes):
        deg = Counter()
        for i, (_, g) in enumerate(feats):
            if not alive[i]:
                continue
            c = g.GetPoints()
            deg[_key(c[0])] += 1; deg[_key(c[-1])] += 1
        drop = []
        for i, (_, g) in enumerate(feats):
            if not alive[i] or g.Length() >= max_hair:
                continue
            c = g.GetPoints()
            da, db = deg[_key(c[0])], deg[_key(c[-1])]
            if min(da, db) == 1 and max(da, db) >= 3:        # free end + junction end = hair
                drop.append(i)
        if not drop:
            break
        for i in drop:
            alive[i] = False
        removed_total += len(drop)

    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(out):
        drv.DeleteDataSource(out)
    dso = drv.CreateDataSource(out); lyo = dso.CreateLayer("lines", sr, ogr.wkbLineString)
    lyo.CreateField(ogr.FieldDefn("arc", ogr.OFTInteger)); defn = lyo.GetLayerDefn()
    kept = 0
    for i, (arc, g) in enumerate(feats):
        if not alive[i]:
            continue
        f = ogr.Feature(defn); f.SetField("arc", int(arc)); f.SetGeometry(g)
        lyo.CreateFeature(f); f = None; kept += 1
    dso = None
    print("[dehair] removed %d hair barbs (<%.0fm); %d -> %d features -> %s"
          % (removed_total, max_hair, len(feats), kept, out))


if __name__ == "__main__":
    a = sys.argv
    run(a[1], a[2],
        float(a[3]) if len(a) > 3 else 4.0,
        int(a[4]) if len(a) > 4 else 3)
