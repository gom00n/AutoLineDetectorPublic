"""Validation metrics: compare a detected-lines result against a hand-cleaned REFERENCE shapefile
so tuning is data-driven instead of eyeballed (base arcgispro-py3 python: ogr + shapely + numpy).

  recall    = % of REFERENCE length whose sampled points lie within TOL of a result line
              (how much of the truth did we detect)
  precision = % of RESULT length whose sampled points lie within TOL of a reference line
              (how much of what we drew is real)
  F1        = harmonic mean.

Both layers must share a CRS (our outputs + hand-cleaned refs are EPSG 6991).
Usage: python metrics.py result.shp reference.shp [tol_m=2.0] [step_m=1.0]
"""
import sys
import numpy as np
from osgeo import ogr
from shapely import wkt
from shapely.ops import unary_union
ogr.UseExceptions()


def load(path):
    ds = ogr.Open(path)
    if ds is None:
        raise SystemExit("cannot open " + path)
    lyr = ds.GetLayer()
    geoms = [wkt.loads(f.GetGeometryRef().ExportToWkt()) for f in lyr if f.GetGeometryRef()]
    return geoms


def sample(lines, step):
    pts = []
    for L in lines:
        n = max(1, int(L.length // step))
        pts += [L.interpolate(k / n, normalized=True) for k in range(n + 1)]
    return pts


def run(res, ref, tol=2.0, step=1.0):
    R, F = load(res), load(ref)
    if not R or not F:
        raise SystemExit("empty result or reference")
    Ru, Fu = unary_union(R), unary_union(F)
    rp, fp = sample(R, step), sample(F, step)
    recall = float(np.mean([p.distance(Ru) <= tol for p in fp]))
    precision = float(np.mean([p.distance(Fu) <= tol for p in rp]))
    f1 = 2 * recall * precision / max(recall + precision, 1e-9)
    rl, fl = sum(x.length for x in R), sum(x.length for x in F)
    print("result    : %d lines, %.1f km" % (len(R), rl / 1000))
    print("reference : %d lines, %.1f km" % (len(F), fl / 1000))
    print("tol=%.1fm  RECALL=%.1f%%  PRECISION=%.1f%%  F1=%.1f%%"
          % (tol, 100 * recall, 100 * precision, 100 * f1))


if __name__ == "__main__":
    a = sys.argv
    if len(a) < 3:
        raise SystemExit(__doc__)
    run(a[1], a[2], float(a[3]) if len(a) > 3 else 2.0, float(a[4]) if len(a) > 4 else 1.0)
