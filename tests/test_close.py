"""Regression tests for the two vector-repair stages on synthetic geometry (no data needed):
  collapse_pairs.collapse -- doubled road casings -> one centerline, rungs dropped, strips kept
  close_parcels.close     -- join/extend loose ends so polygonize can close parcels

Run with the base interpreter:
  "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" tests\\test_close.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from shapely.geometry import LineString
from shapely.ops import unary_union, polygonize
from collapse_pairs import collapse
from close_parcels import close

fails = []
n_checks = 0


def check(name, cond, msg=""):
    global n_checks
    n_checks += 1
    if cond:
        print("  PASS  %s" % name)
    else:
        print("  FAIL  %s  %s" % (name, msg))
        fails.append(name)


# ---------- collapse_pairs ----------
# a doubled road (2.4 m apart) with one rung, plus two REAL strip borders 8 m apart
A = LineString([(0, 0), (100, 0)])
B = LineString([(0, 2.4), (100, 2.4)])
rung = LineString([(50, 0), (50, 2.4)])
C = LineString([(0, 50), (100, 50)])
D = LineString([(0, 58), (100, 58)])
out, n_pairs, n_rungs = collapse([A, B, rung, C, D], pair_d=3.0, frac=0.7)
center = [L for L in out if abs(np.mean([y for _, y in L.coords]) - 1.2) < 0.05]
check("collapse_pair_to_centerline", n_pairs == 1 and len(center) == 1,
      "pairs=%d out=%d" % (n_pairs, len(out)))
check("collapse_drops_rung", n_rungs == 1 and len(out) == 3,
      "rungs=%d out=%d" % (n_rungs, len(out)))
check("collapse_keeps_strips", any(L.equals(C) for L in out) and any(L.equals(D) for L in out))

# nothing doubled -> untouched
out2, p2, r2 = collapse([C, D], pair_d=3.0, frac=0.7)
check("collapse_noop", p2 == 0 and r2 == 0 and len(out2) == 2)

# ---------- close_parcels ----------
# square with a 1.5 m corner gap + a divider stopping 4 m short of the top edge
lines = [
    list(LineString([(0, 0), (100, 0)]).coords),          # bottom
    list(LineString([(100, 0), (100, 100)]).coords),      # right
    list(LineString([(100, 100), (0, 100)]).coords),      # top
    list(LineString([(0, 100), (0, 1.5)]).coords),        # left, 1.5 m short of the corner
    list(LineString([(50, 0), (50, 96)]).coords),         # divider, 4 m short of the top
]
nd, nj, ne = close(lines, ext_d=5.0, join_d=3.0)
polys = [p for p in polygonize(unary_union([LineString(cs) for cs in lines])) if p.area >= 30]
check("close_joined_corner", nj >= 1, "joined=%d" % nj)
check("close_extended_divider", ne >= 1, "extended=%d" % ne)
check("close_two_parcels", len(polys) == 2 and all(4500 < p.area < 5500 for p in polys),
      "polys=%s" % [round(p.area) for p in polys])

# near-closed ring (one line, 2 m gap) -> self-join closes it
ring = [list(LineString([(0, 0), (100, 0), (100, 100), (0, 100), (0, 2)]).coords)]
nd, nj, ne = close(ring, ext_d=5.0, join_d=3.0)
polys = [p for p in polygonize(unary_union([LineString(ring[0])])) if p.area >= 30]
check("close_ring_selfjoin", nj == 1 and len(polys) == 1 and abs(polys[0].area - 10000) < 200,
      "joined=%d polys=%s" % (nj, [round(p.area) for p in polys]))

# far dangles must NOT be joined
farlines = [list(LineString([(0, 0), (10, 0)]).coords), list(LineString([(50, 50), (60, 50)]).coords)]
nd, nj, ne = close(farlines, ext_d=5.0, join_d=3.0)
check("close_far_untouched", nj == 0 and ne == 0)

print()
if fails:
    print("FAILED: %d test(s): %s" % (len(fails), ", ".join(fails)))
    sys.exit(1)
print("all %d close/collapse tests passed" % n_checks)
