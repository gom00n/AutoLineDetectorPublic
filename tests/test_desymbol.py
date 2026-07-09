"""Regression tests for desymbol.desymbol on synthetic geometry (no data needed).
Encodes the M14 hand-cleaning classes of 2026-07-07: tick rows between casings go, a real
border crossing the track stays (collinear continuation), glyph rings and isolated compact
clusters go, sub-row-min strokes stay.

Run with the base interpreter:
  "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" tests\\test_desymbol.py
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from shapely.geometry import LineString
from desymbol import desymbol

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


# world: a double track (two casings 2.5 m apart), a border crossing it at x=50,
# a row of 4 symbol ticks, a glyph ring, an isolated triangle, and 2 lone strokes.
lines = []
lines += [LineString([(0, 0), (48, 0)]), LineString([(52, 0), (100, 0)])]          # casing A
lines += [LineString([(0, 2.5), (48, 2.5)]), LineString([(52, 2.5), (100, 2.5)])]  # casing B
# a real border crossing the track at x=50 (split at both casings)
cross_mid = LineString([(50, 0), (50, 2.5)])
lines += [LineString([(50, -30), (50, 0)]), cross_mid, LineString([(50, 2.5), (50, 30)])]
# symbol ticks between the casings (a periodic row; midpoints collinear, 15 m apart)
ticks = [LineString([(x, 0), (x, 2.5)]) for x in (15, 30, 45, 60, 75)]
lines += ticks
# glyph ring (octagon, perimeter ~ 25 m) sitting off to the side
t = np.linspace(0, 2 * np.pi, 9)
ring = LineString([(200 + 4 * np.cos(a), 50 + 4 * np.sin(a)) for a in t])
lines.append(ring)
# isolated compact triangle glyph
tri = [LineString([(300, 80), (308, 80)]), LineString([(308, 80), (304, 87)]),
       LineString([(304, 87), (300, 80)])]
lines += tri
# only TWO strokes on another double track (below row_min) -> stay; the track keeps them
# non-isolated so the cluster rule does not apply (as post-filter real data guarantees)
lines += [LineString([(0, 50), (40, 50)]), LineString([(0, 53), (40, 53)])]
lone = [LineString([(10, 50), (10, 53)]), LineString([(25, 50), (25, 53)])]
lines += lone

kept, dropped = desymbol(lines)
di = {i for i, _ in dropped}
tags = dict(dropped)


def idx(g):
    return next(i for i, L in enumerate(lines) if L.equals(g))


check("tick_row_dropped", all(idx(t_) in di for t_ in ticks),
      "tags=%s" % [(i, t) for i, t in dropped])
check("crossing_kept", idx(cross_mid) not in di, "crossing tagged %s" % tags.get(idx(cross_mid)))
check("casings_kept", all(idx(lines[k]) not in di for k in range(4)))
check("ring_dropped", idx(ring) in di, "ring not dropped")
check("triangle_cluster_dropped", all(idx(g) in di for g in tri),
      "tags=%s" % [(i, t) for i, t in dropped])
check("lone_pair_kept", all(idx(g) not in di for g in lone))

# protect wins over every rule
kept2, dropped2 = desymbol(lines, protect={idx(t_) for t_ in ticks})
check("protect_respected", all(idx(t_) not in {i for i, _ in dropped2} for t_ in ticks))

print()
if fails:
    print("FAILED: %d test(s): %s" % (len(fails), ", ".join(fails)))
    sys.exit(1)
print("all %d desymbol tests passed" % n_checks)
