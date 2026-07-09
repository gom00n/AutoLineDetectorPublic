"""Central configuration for the LineDetector pipeline -- the single source of truth for
interpreter paths, project roots, the modern cadastre, and per-project default parameters.
The GUI, the batch runner, and extract_9a import from here so paths aren't duplicated.

ALL projects use the same engine: extract_9a (Motzag-GDB corridor + per-polygon windowed read),
differing only by ink ("black" global threshold / "red" adaptive local-mean) and by their
INPUT (srcdir) and OUTPUT (outdir) folders -- both overridable in the GUI.
Edit the paths below to move the pipeline to another machine or project.
"""
import os, re, glob

HERE = os.path.dirname(os.path.abspath(__file__))

# --- interpreters (two are needed: ArcGIS base + a scikit-image venv) ---
BASE_PY = r"C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3\python.exe"
VENV_PY = os.path.join(HERE, "venv_skimage", "Scripts", "python.exe")

# --- project roots ---
# Point DATA_ROOT at wherever your scans / GDBs live (override via env var).
DATA_ROOT = os.environ.get("LINEDETECTOR_DATA", r"C:\Data\LineDetector")
BLACK_ROOT = os.path.join(DATA_ROOT, "black_ink")     # black-ink 1:20k topo-cadastral sheets
RED1_ROOT = os.path.join(DATA_ROOT, "red_ink_1")      # red ink over a yellow wash
RED2_ROOT = os.path.join(DATA_ROOT, "red_ink_2")      # red ink, double-line casings
EPSG = 6991                                            # Israeli Grid 05/12
# parcels to snap onto: a local GDB of the clean (unmodified) cadastre
PARCELS_GDB = os.path.join(DATA_ROOT, "reference_cadastre", "Layers.gdb")
PARCELS_LAYER = "PARCEL_ALL"

# --- black ink (global threshold) ---
NINE_SRCDIR = os.path.join(BLACK_ROOT, "scans")
NINE_OUT = os.path.join(BLACK_ROOT, "work")
NINE_GDB = os.path.join(BLACK_ROOT, "corridor_project", "corridor.gdb")
NINE_LAYERS = ["PolyMap_black_1", "PolyMap_black_2", "PolyMap_black_1_and_2"]

# --- red ink #1 (adaptive, over wash) ---
FIVE_SRCDIR = RED1_ROOT                                # scans sit at the project root
FIVE_OUT = RED1_ROOT                                   # outputs go to RED1_ROOT\<stem>\ (existing layout)
FIVE_GDB = os.path.join(RED1_ROOT, "ArcGis", "corridor.gdb")   # set its Motzag layers via "Choose..."

# --- red ink #2 (adaptive, double-line) -- scans and outputs under the data root ---
EIGHT_SRCDIR = os.path.join(RED2_ROOT, "Scans")
EIGHT_OUT = RED2_ROOT
EIGHT_GDB = os.path.join(RED2_ROOT, "corridor_project", "corridor.gdb")
EIGHT_LAYERS = ["Poly_map2", "Poly_Map"]               # Motzag-tagged corridor polygon layers (true-curve)

STEM_RE = re.compile(r"^M\d+_\d+$")                    # a sheet stem, e.g. M8_9076 / M11_8457

# shared vector-backend defaults (same for every project; tune per-mode where they differ)
_BACKEND = dict(gap="25", cos="0.80", simp_m="0.5", hair="4", passes="3")


def _defaults(srcdir, outdir, ink, gdb, layers, **over):
    """Build a mode's default-parameter dict; `over` overrides the red-ink baseline below."""
    d = dict(srcdir=srcdir, outdir=outdir,
             cap=("212" if ink == "red" else "145"), offset="14", block="101",
             min_px=("400" if ink == "red" else "30"), buf=("20" if ink == "red" else "35"),
             close_r=("3" if ink == "red" else "1"), prune=("40" if ink == "red" else "15"),
             # max_fill: revert closing fills > N px so doubled casings stay TWO lines
             # (double-line fix, 2026-07-07); 0 = off for black ink (close_r=1 can't weld pairs)
             max_fill=("100" if ink == "red" else "0"),
             spur=("12" if ink == "red" else "4"), loop=("30" if ink == "red" else "6"),
             simp_px=("1.5" if ink == "red" else "1.0"),
             # hair/min_comp validated 2026-07-06 against a hand-cleaned reference sheet: do NOT
             # raise them. hair=8/min_comp=30 ate 19 real lines while catching 2/45 of the junk;
             # the hand-deleted junk is junction-attached symbol strokes (ticks/X-marks between
             # double casings, traced glyphs) that no length threshold can separate from borders.
             min_comp=("20" if ink == "red" else "40"),
             snap_tol=("3" if ink == "red" else "8"), snap_frac=("0.6" if ink == "red" else "0.5"),
             snap_run="10", pair_d="3.0", ext_m="5", join_m="3", min_area="30",
             stroke_len="12", row_gap="30", ring_perim="45",
             gdb=gdb, layers=layers, snap_layers="",
             parcels_gdb=PARCELS_GDB, parcels_layers=PARCELS_LAYER)
    d.update(_BACKEND)
    d.update(over)
    return d


def _mode(ink, srcdir, outdir, gdb, layers, **over):
    return dict(extract="extract_9a.py", overlay="overlay_9a.py", ink=ink, src_ext=".tif",
                has_extract_params=True, defaults=_defaults(srcdir, outdir, ink, gdb, layers, **over))


# --- projects for the GUI / batch runner (all Motzag-GDB based) ---
PROJECTS = {
    "Black ink - global threshold (1:20k)": _mode("black", NINE_SRCDIR, NINE_OUT, NINE_GDB, ",".join(NINE_LAYERS)),
    "Red ink - adaptive (over wash)":       _mode("red",   FIVE_SRCDIR, FIVE_OUT, FIVE_GDB, ""),
    # snap_layers stays EMPTY for the double-line red mode: snap ONLY to the reference parcels
    # (parcels_gdb/PARCEL_ALL) -- mixing in the corridor-GDB block outlines gave near-parallel
    # duplicate reference edges and caused the saw-tooth snapping artifact.
    "Red ink - adaptive, double-line":      _mode("red",   EIGHT_SRCDIR, EIGHT_OUT, EIGHT_GDB, ",".join(EIGHT_LAYERS),
                                            snap_tol="3"),
}


def list_sheets(srcdir, src_ext=".tif"):
    """All sheet stems found in a source folder (filtered to M<n>_<n>). Looks for both .tif and
    .jpg scans -- some projects deliver a few sheets as .jpg (+ .jgwx world file) instead."""
    out = []
    for ext in {src_ext, ".jpg"}:
        for p in glob.glob(os.path.join(srcdir or "", "*" + ext)):
            stem = os.path.splitext(os.path.basename(p))[0]
            if STEM_RE.match(stem):
                out.append(stem)
    return sorted(set(out))


def resolve_stems(text, srcdir, src_ext=".tif"):
    """Turn the GUI 'sheets' field into a list of stems, looked up in `srcdir`. Accepts:
       - a single stem or a comma/space/newline-separated list (glob '*' allowed),
       - a FOLDER path typed straight in (enumerate its sheets, overriding srcdir),
       - 'all' / '*' (every sheet in srcdir)."""
    text = (text or "").strip()
    if not text:
        return []
    if os.path.isdir(text):                           # a folder typed straight into Sheets
        return list_sheets(text, src_ext)
    if text.lower() in ("all", "*"):
        return list_sheets(srcdir, src_ext)
    avail = list_sheets(srcdir, src_ext)
    out = []
    for tok in re.split(r"[,\s]+", text):
        if not tok:
            continue
        if "*" in tok or "?" in tok:                  # glob against available stems
            import fnmatch
            out += [s for s in avail if fnmatch.fnmatch(s, tok)]
        else:
            out.append(tok)
    seen = set()
    return [s for s in out if not (s in seen or seen.add(s))]
