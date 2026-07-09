"""Merge touching polylines into longer continuous lines, then straighten (base arcpy python).

RasterToPolyline splits the border mesh at every junction, so one border becomes many short
features. UnsplitLine dissolves all features that meet at a 2-valent node (a simple pass-through)
into a single polyline -- it stops only at true 3+ junctions, where merging would be ambiguous.
SimplifyLine (POINT_REMOVE) then removes the 1-px staircase so lines read straight.

Both tools are core (Data Management / Cartography) and run on a Basic license.

Usage: python merge_lines.py in_lines.shp out_lines.shp [simplify_m=0.30]
"""
import sys
import arcpy

arcpy.env.overwriteOutput = True


def _stats(fc):
    n = int(arcpy.management.GetCount(fc)[0])
    tot = sum(r[0] for r in arcpy.da.SearchCursor(fc, ["SHAPE@LENGTH"]))
    return n, tot


def run(src, out, simplify_m=0.30):
    n0, len0 = _stats(src)
    merged = "in_memory/_merged"
    # UnsplitLine is NOT licensed on Basic (ERROR 000824); Dissolve(UNSPLIT_LINES) is the
    # core equivalent -- merges runs through 2-valent nodes, splits at true 3+ junctions.
    arcpy.management.Dissolve(src, merged, "", "", "SINGLE_PART", "UNSPLIT_LINES")
    n1, len1 = _stats(merged)
    if simplify_m and simplify_m > 0:
        try:
            arcpy.cartography.SimplifyLine(
                merged, out, "POINT_REMOVE", simplify_m,
                error_resolving_option="RESOLVE_ERRORS",
                collapsed_point_option="NO_KEEP", error_checking_option="NO_CHECK")
        except arcpy.ExecuteError as e:
            print("[merge] SimplifyLine unavailable (%s); keeping unsimplified"
                  % str(e).splitlines()[-1].strip())
            arcpy.management.CopyFeatures(merged, out)
    else:
        arcpy.management.CopyFeatures(merged, out)
    n2, len2 = _stats(out)
    print("[merge] features %d -> %d (unsplit) -> %d (simplify)" % (n0, n1, n2))
    print("[merge] avg length %.1f -> %.1f -> %.1f m;  total length %.0f -> %.0f m"
          % (len0 / max(n0, 1), len1 / max(n1, 1), len2 / max(n2, 1), len0, len2))
    print("[merge] wrote %s" % out)


if __name__ == "__main__":
    a = sys.argv
    run(a[1], a[2], float(a[3]) if len(a) > 3 else 0.30)
