"""Regression tests for the snap stage (snap_modern.build_network / snap_line) on synthetic
geometry -- no GDB, no shapefiles. Encodes every failure diagnosed on the M13/M14 sheets so
future snap tuning can't silently regress.

Run with the base interpreter:
  "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" tests\\test_snap.py

Reference world: one row of three 100 x 100 m parcels --
  (0,0)-(100,100) | (100,0)-(200,100) | (200,0)-(300,100)
so y=0 / y=100 are long boundary chains crossing T-junctions at x=100 and x=200, and
x=100 / x=200 are shared (doubled) parcel dividers that the network union must dedupe.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from shapely.geometry import Polygon, LineString
from snap_modern import build_network, snap_line

TOL, FRAC, MIN_RUN = 3.0, 0.6, 10.0

SQUARES = [Polygon([(x, 0), (x + 100, 0), (x + 100, 100), (x, 100)]) for x in (0, 100, 200)]
RINGS = [p.exterior for p in SQUARES]
EDGES, TREE = build_network(RINGS)

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


def snap(L):
    return snap_line(L, EDGES, TREE, TOL, FRAC, MIN_RUN)


def max_dy(coords, y):
    return max(abs(cy - y) for _, cy in coords)


# --- network construction: shared dividers deduped, total length = the true union ---
total = sum(e.length for e in EDGES)
check("network_dedupe", abs(total - 1000.0) < 1e-6,
      "edge length %.1f != 1000 (doubled shared edges not deduped?)" % total)

# --- 1. line 0.5 m parallel to one edge -> snapped onto it (M14: 0.5 m gap unsnapped) ---
parts = snap(LineString([(10, 0.5), (90, 0.5)]))
check("parallel_close", len(parts) == 1 and parts[0][1] == 1
      and max_dy(parts[0][0].coords, 0.0) < 0.05,
      "parts=%s" % [(p.length, f) for p, f in parts])

# --- 2. line following a CHAIN of three parcels' edges -> snapped (the multi-ring failure) ---
parts = snap(LineString([(10, 100.7), (290, 100.7)]))
check("multi_edge_chain", len(parts) == 1 and parts[0][1] == 1
      and max_dy(parts[0][0].coords, 100.0) < 0.05,
      "parts=%s" % [(round(p.length, 1), f) for p, f in parts])

# --- 3. diagonal line far from any boundary -> kept as ink, geometry untouched ---
L = LineString([(30, 30), (70, 70)])
parts = snap(L)
check("diagonal_kept", len(parts) == 1 and parts[0][1] == 0
      and parts[0][0].equals(L), "parts=%s" % [(round(p.length, 1), f) for p, f in parts])

# --- 4. half-along / half-divergent -> SPLIT into a snapped part + an ink part ---
L = LineString([(10, 0.5), (90, 0.5), (140, 40)])
parts = snap(L)
ok = (len(parts) == 2 and parts[0][1] == 1 and parts[1][1] == 0
      and max_dy(parts[0][0].coords, 0.0) < 0.05                       # snapped part on y=0
      and abs(parts[1][0].coords[-1][0] - 140) < 1e-6                  # ink part ends at the ink end
      and abs(parts[1][0].coords[-1][1] - 40) < 1e-6
      and parts[0][0].coords[-1] == parts[1][0].coords[0])             # no gap at the transition
check("partial_split", ok, "parts=%s" % [(round(p.length, 1), f) for p, f in parts])

# --- 5. line crossing a perpendicular divider -> NO teeth, kept exactly as ink ---
L = LineString([(60, 50), (140, 50)])
parts = snap(L)
check("no_teeth_crossing", len(parts) == 1 and parts[0][1] == 0 and parts[0][0].equals(L),
      "parts=%s" % [(round(p.length, 1), f) for p, f in parts])

# --- 6. line hugging the boundary THROUGH a junction -> stays snapped across it
#        (M14: 'snapped fine until it unsnapped near intersection') ---
parts = snap(LineString([(10, 0.5), (190, 0.5)]))
check("junction_continuation", len(parts) == 1 and parts[0][1] == 1
      and max_dy(parts[0][0].coords, 0.0) < 0.05
      and abs(parts[0][0].length - 180) < 5,
      "parts=%s" % [(round(p.length, 1), f) for p, f in parts])

# --- 7. short line (< min_run): all-or-nothing at match_frac ---
parts = snap(LineString([(10, 0.4), (16, 0.4)]))
check("short_all_snap", len(parts) == 1 and parts[0][1] == 1
      and max_dy(parts[0][0].coords, 0.0) < 0.05,
      "parts=%s" % [(round(p.length, 1), f) for p, f in parts])
L = LineString([(40, 40), (46, 44)])
parts = snap(L)
check("short_kept", len(parts) == 1 and parts[0][1] == 0 and parts[0][0].equals(L))

# --- 8. KINK guard (red-ink pic 1): a straight line passing 2.2 m from a boundary that only
#        exists along part of it must NOT grab it -- that would put a Z-jog in a straight
#        historical border. The whole line stays ink. ---
L = LineString([(-80, 2.2), (90, 2.2)])
parts = snap(L)
check("kink_reject_parallel_grab", len(parts) == 1 and parts[0][1] == 0 and parts[0][0].equals(L),
      "parts=%s" % [(round(p.length, 1), f) for p, f in parts])

# --- 9. ...but the same shape at 0.5 m offset IS the boundary (sub-noise jump): split into
#        an ink tail + a snapped stretch, sharing the transition coordinate. ---
L = LineString([(-80, 0.5), (90, 0.5)])
parts = snap(L)
ok = (len(parts) == 2 and parts[0][1] == 0 and parts[1][1] == 1
      and max_dy(parts[1][0].coords, 0.0) < 0.05
      and parts[0][0].coords[-1] == parts[1][0].coords[0])
check("low_offset_still_splits", ok, "parts=%s" % [(round(p.length, 1), f) for p, f in parts])

# --- 10. SAW guard (red-ink pic 2): a line down the middle of a 2.4 m sliver (two parallel
#         edges) must not oscillate between them -- kept as ink, geometry untouched. ---
SLIVER = [Polygon([(0, 0), (100, 0), (100, 2.4), (0, 2.4)]).exterior]
SEDGES, STREE = build_network(SLIVER)
L = LineString([(x, 1.2 + 0.18 * (-1) ** x) for x in range(5, 96)])
parts = snap_line(L, SEDGES, STREE, TOL, FRAC, MIN_RUN)
check("saw_reject_sliver", all(f == 0 for _, f in parts)
      and abs(sum(p.length for p, _ in parts) - L.length) < 0.5,
      "parts=%s" % [(round(p.length, 1), f) for p, f in parts])

# --- 11. wobbly ink along ONE real boundary still snaps (feet are SHORTER than the ink,
#         so the saw guard must not fire on legitimate matches). ---
L = LineString([(x, 0.4 + 0.3 * (-1) ** x) for x in range(5, 96)])
parts = snap_line(L, SEDGES, STREE, TOL, FRAC, MIN_RUN)
check("wobble_still_snaps", len(parts) == 1 and parts[0][1] == 1
      and max_dy(parts[0][0].coords, 0.0) < 0.05,
      "parts=%s" % [(round(p.length, 1), f) for p, f in parts])

print()
if fails:
    print("FAILED: %d test(s): %s" % (len(fails), ", ".join(fails)))
    sys.exit(1)
print("all %d snap tests passed" % n_checks)
