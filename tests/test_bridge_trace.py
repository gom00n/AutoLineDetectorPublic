"""Regression tests for the two raster-topology stages, on synthetic skeletons (no scans needed):
  bridge_gaps -- a stub must be extended onto the border it stops short of, at EVERY gap size
  trace_lines -- an arc linking two junctions a few pixels apart must survive

Both encode bugs found by review on 2026-11-15, not by a failing sheet, so they had been
shaping the output silently:
  * outward_dirs took its centroid over every skeleton pixel near the endpoint, including the
    border the stub was about to meet. For gaps of 2-8 px the vector flipped, pass 2 cast its
    ray backwards into the stub and logged a link it never drew.
  * _trace read a 1-px arc between two junction blobs as a loop off ONE junction and dropped
    it, cutting a hole in the through-border.

Run with the scikit-image venv (bridge_gaps imports skimage):
  venv_skimage\\Scripts\\python.exe tests\\test_bridge_trace.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import bridge_gaps as B
import trace_lines as T

fails = []
n_checks = 0


def check(name, cond, msg=""):
    global n_checks
    n_checks += 1
    print(("  PASS  " if cond else "  FAIL  ") + name + ("" if cond else "  " + msg))
    if not cond:
        fails.append(name)


def bridge(sk, max_gap=25, cos=0.80):
    """Both passes, as bridge_gaps.run does them. Returns (n_ep_ep, n_ep_line)."""
    eps = B.endpoints(sk.astype(bool))
    dirs = B.outward_dirs(sk, eps)
    n1, used = B.bridge_ep_ep(sk, eps, dirs, max_gap, cos)
    return n1, B.bridge_ep_line(sk, eps, dirs, used, max_gap, cos)


# ---------- bridge: T-junction, every gap from 2 px up to max_gap ----------
bad = []
for gap in range(2, 25):
    sk = np.zeros((90, 90), np.uint8)
    sk[60, 5:85] = 1                                  # horizontal border
    tip = 60 - gap - 1                                # `gap` empty pixels below the stub's tip
    sk[tip - 30:tip + 1, 45] = 1                      # vertical stub
    bridge(sk)
    if not sk[tip + 1:60, 45].all():
        bad.append(gap)
check("bridge_T_every_gap_2_to_24px", not bad, "not bridged at gaps %s" % bad)

# the direction itself: a stub 5 px short of a border must point AT the border
sk = np.zeros((90, 90), np.uint8)
sk[60, 5:85] = 1; sk[24:55, 45] = 1
eps = B.endpoints(sk.astype(bool)); dirs = B.outward_dirs(sk, eps)
k = [i for i, e in enumerate(eps) if tuple(e) == (54, 45)][0]
check("bridge_dir_ignores_neighbour_line", dirs[k][0] > 0.99, "dir=%s" % dirs[k])

# beyond max_gap nothing is drawn
sk = np.zeros((90, 90), np.uint8)
sk[80, 5:85] = 1; sk[10:41, 45] = 1                   # 39 px short
before = int(sk.sum()); bridge(sk)
check("bridge_T_beyond_max_gap_untouched", int(sk.sum()) == before)

# a stub aimed at an 8-connected DIAGONAL border must stop at it, whichever pixel parity it has
bad = []
for shift in range(8):
    sk = np.zeros((80, 80), np.uint8)
    for i in range(5, 75):
        sk[i, i] = 1
    tip = (34, 46 + shift)
    for i in range(18):
        sk[tip[0] - i, tip[1] + i] = 1
    before = sk.copy(); bridge(sk)
    added = np.argwhere((sk > 0) & (before == 0))
    if len(added) == 0 or any(r > c for r, c in added):   # nothing drawn, or leaked through
        bad.append(shift)
check("bridge_T_onto_diagonal_no_leak", not bad, "failed at shifts %s" % bad)

# a clean break in one straight line is still rejoined end to end (pass 1 unchanged)
sk = np.zeros((40, 90), np.uint8)
sk[20, 5:40] = 1; sk[20, 52:85] = 1
n1, n2 = bridge(sk)
check("bridge_collinear_break_rejoined", n1 == 1 and sk[20, 5:85].all(), "n1=%d" % n1)

# two parallel lines 4 px apart must NOT be linked to each other (doubled casings stay two)
sk = np.zeros((40, 90), np.uint8)
sk[18, 10:80] = 1; sk[22, 10:80] = 1
before = int(sk.sum()); bridge(sk)
check("bridge_parallel_casings_untouched", int(sk.sum()) == before)

# walk_back stops at a junction instead of wandering down another branch
sk = np.zeros((40, 40), np.uint8)
sk[20, 5:35] = 1; sk[10:20, 20] = 1                   # T: stub from (10,20) down to the bar
own = B.walk_back(sk, (10, 20), 30)
check("walk_back_stops_at_junction", own and all(c == 20 for _, c in own) and len(own) <= 9,
      "walked %s" % own)


# ---------- trace: an arc linking two close junctions survives ----------
def trace(sk, min_spur=4, min_loop=6):
    arcs, n_spur, n_loop, _ = T._trace(np.pad(sk, 1), min_spur, min_loop, 1.0)
    return [[(r - 1, c - 1) for r, c in a] for a in arcs], n_spur, n_loop


def cols_covered(arcs, row):
    out = set()
    for a in arcs:
        for (r0, c0), (r1, c1) in zip(a[:-1], a[1:]):
            if abs(r0 - row) < 0.6 and abs(r1 - row) < 0.6:
                out.update(range(int(round(min(c0, c1))), int(round(max(c0, c1))) + 1))
    return out


# one straight border with a branch UP at col 8 and a branch DOWN at col `c2`
for c2 in (12, 13, 14):
    sk = np.zeros((14, 30), bool)
    sk[6, 1:28] = True; sk[1:6, 8] = True; sk[7:13, c2] = True
    arcs, _, n_loop = trace(sk)
    missing = [c for c in range(1, 28) if c not in cols_covered(arcs, 6)]
    check("trace_border_through_two_junctions_%dpx_apart" % (c2 - 8), not missing and n_loop == 0,
          "uncovered cols %s, loops dropped %d" % (missing, n_loop))

# the linking arc must end on BOTH junction nodes (shared coordinates = real topology)
sk = np.zeros((14, 30), bool)
sk[6, 1:28] = True; sk[1:6, 8] = True; sk[7:13, 12] = True
arcs, _, _ = trace(sk)
nodes = {}
for a in arcs:
    for p in (a[0], a[-1]):
        nodes[(round(p[0], 3), round(p[1], 3))] = nodes.get((round(p[0], 3), round(p[1], 3)), 0) + 1
check("trace_two_junctions_are_both_degree_3", sorted(nodes.values()).count(3) == 2,
      "node degrees %s" % sorted(nodes.values()))

# a genuine small bubble hanging off ONE junction is still dropped as a loop
# (a lollipop: a stem rising from the border, with a 7-px ring closing on the stem's top pixel)
sk = np.zeros((20, 30), bool)
sk[10, 1:28] = True; sk[6:10, 14] = True
for r, c in ((5, 13), (4, 12), (3, 13), (2, 14), (3, 15), (4, 16), (5, 15)):
    sk[r, c] = True
arcs, _, n_loop = trace(sk, min_loop=30)
check("trace_real_bubble_still_dropped", n_loop == 1 and
      not any(r < 4.5 for a in arcs for r, _ in a), "loops=%d" % n_loop)

print()
if fails:
    print("FAILED: %d test(s): %s" % (len(fails), ", ".join(fails)))
    sys.exit(1)
print("all %d bridge/trace tests passed" % n_checks)
