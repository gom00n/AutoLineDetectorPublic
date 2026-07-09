# LineDetector

Extract parcel-border **lines** from scanned historical cadastral map sheets (raster) and
vectorize them to georeferenced polylines, on an ArcGIS Pro **Basic** license (no Spatial
Analyst). Built for Israeli cadastral projects (Israel TM Grid / Gush-Parcel), but the
vector backend is general.

## One engine, two inks

**Every project runs the same engine** (`extract_9a.py`): the target area comes from a GDB — the
polygons whose `Motzag ==` the scan number parsed from the file name (`M8_9076` → 8) — read as one
small window per polygon. Only the ink detection differs, because the map series hide their ink
differently:

- **RED ink** — dark-red lines, often drawn **over a yellow project-area wash**.
  The wash lightens the ink to only ~15–20 DN darker than its background, so a global
  `R < threshold` catches lines on white paper but **misses ~90% of the lines inside the colored
  area** (exactly the area of interest). Fix: an **adaptive (local-mean) threshold** on the Red
  band (`--ink red`).
- **BLACK ink** — 1:20,000 topo-cadastral sheets (~0.85 m/px), with red numbers and
  purple dash-dot admin lines as colored confusers. Black ink is dark in **all** channels, so a
  global `max(R,G,B) < CAP` isolates it and rejects every colored confuser (red leaks 0%; purple
  ~18%, killed by a `B-G` exclusion). No adaptive threshold needed (`--ink black`).

## Pipeline

```
source RGB TIFF
  └─ extract_9a.py        binarize (--ink red: adaptive local-mean | black: global cap)  [base: gdal+scipy]
                          → clip to the buffered Motzag corridor (per-polygon windowed read)
  └─ skeletonize.py       gap-close → skimage.skeletonize → spur-prune ("de-hair")     [venv: scikit-image]
  └─ bridge_gaps.py       link facing+collinear endpoints across gaps (continuity)     [venv: scikit-image]
  └─ trace_lines.py       trace skeleton graph → topo-clean polylines (NOT RasterToPolyline) [base: gdal/ogr]
                          → junctions share one node; spur/loop prune; RDP straighten; recover CRS
  └─ merge_lines.py       Dissolve(UNSPLIT_LINES) → merge runs into longer lines        [base: arcpy]
  └─ filter_components.py drop small isolated components (legend/table/symbols) + simplify [base: ogr]
  → <stem>_lines.shp      (georeferenced polylines, EPSG 6991)
  └─ snap_modern.py       (optional) conflate onto the reference-boundary network        [base: ogr+shapely]
  → <stem>_lines_snapped.shp
  └─ overlay_9a.py        QA PNG over the faded scan (green = snapped, red = ink) → <stem>_preview.png
```

## Running it — the GUI

`run_gui.bat` (or `pipeline_gui.py` under the ArcGIS base python) is a Tkinter front end over all
the scripts. Pick a **project mode** (black-ink / red-ink presets), enter one or many **sheets** (a
single stem like `M8_9076`, a comma/space list, a `*` glob, a folder, or `all`; "Scan" fills in
every sheet), tweak parameters (**hover any field for help**), tick the **stages**, and Run. A
**preflight check** validates interpreters, source sheets, and every named GDB layer before
anything runs. It shells out to the correct interpreter per stage, streams the log, and runs
sheets as a batch.

Paths and per-project defaults live in **`config.py`** — the single source of truth. Edit it to
relocate the pipeline to another machine or project.

**New to the parameters?** Open **`GUI_GUIDE.html`** in a browser — an illustrated, plain-language
guide to every knob in the GUI (what it does, which way to turn it, sensible defaults) plus how to
label reference/artifact data to measure and improve results.

### Outputs (per sheet, in the work folder)
- `<stem>_lines.shp` — the final detected polylines (+ `.shx/.dbf/.prj`).
- `<stem>_lines_snapped.shp` — present only if the Snap stage ran.
- `<stem>_preview.png` — QA overlay (lines over the faded scan).
- `<stem>.log` — the run log for that sheet.

Intermediate rasters (`_mask` / `_skel` / `_skel_bridged`) and the trace/merge shapefiles go to a
**temp folder** by default, so the work folder stays clean. The GUI's **"keep intermediates"**
checkbox writes them into the work folder instead (handy when re-running single stages to tune).

### Core scripts (CLI, if you bypass the GUI)
- **`extract_9a.py <stem> [--ink black|red] [--cap N] [--offset 14] [--block 101] [--min-px N]
  [--buf N] [--gdb PATH] [--layers "L1,L2"] [--srcdir DIR] [--out MASK]`** — ink mask clipped to
  the buffered **Motzag** corridor. Exit **0** = wrote mask, **3** = legit skip (no Motzag
  polygons / missing source — batch continues), **2** = real config error (e.g. GDB can't open).
- **`skeletonize.py in.tif out.tif [close_r=3] [prune_len=40]`** *(venv)* — close tiny gaps, thin
  to 1-px, prune spurs.
- **`bridge_gaps.py in_skel.tif out_skel.tif [max_gap=25] [cos=0.80] [ray=1]`** *(venv)* — rejoin
  broken borders by linking endpoints that are close AND collinear/facing; does not re-skeletonize.
- **`trace_lines.py in_skel.tif out.shp [min_spur=12] [min_loop=30] [simplify_px=1.5]`** —
  skeleton-graph vectoriser giving exact topology; stamps the scan CRS.
- **`merge_lines.py in.shp out.shp [simplify_m=0]`** — `Dissolve(UNSPLIT_LINES)`.
- **`filter_components.py in.shp out.shp [min_comp_len=20] [simplify_m=0.5]`** — junk removal by
  connected-component size + `SimplifyPreserveTopology`.
- **`dehair_lines.py in.shp out.shp [max_hair=4] [passes=3]`** — remove short hair barbs.
- **`snap_modern.py in.shp out.shp [--tol 3] [--match-frac 0.6] [--min-run 10]
  [--snap-src GDB "L1,L2"]...`** — per-point network conflation with per-run splitting (see
  findings). `--snap-src` is repeatable and takes a GDB or shapefile; all sources are pooled into
  one reference network. Writes a `snapped` field: 0 = kept as ink, 1 = on a reference boundary.
  Regression tests: `tests\test_snap.py` (run with the base python; pure-shapely, no data needed).
- **`overlay_9a.py <stem> [--srcdir DIR] [--workdir DIR]`** — QA overlay PNG; prefers
  `_lines_snapped.shp` and colors green (snapped) / red (kept as ink).
- **`metrics.py result.shp reference.shp [tol] [step]`** — recall / precision / F1 vs a hand-
  cleaned reference, for data-driven tuning.

## Key findings (hard-won)
- **Match the detection to the ink.** The red-over-yellow ink needs an **adaptive** threshold (a
  global `R<threshold` is blind over the wash). Black ink is dark in all channels, so a
  **global `max(R,G,B) < CAP`** isolates it and rejects red/purple confusers — simpler, and the
  whole vector backend transfers unchanged.
- **Target only the area of interest.** The corridor comes from a **GDB**: the polygons with
  **`Motzag == <scan #>`** (parsed from the file name, e.g. `M8_9076` → 8) across the chosen
  layers, buffered. Those polygons are scattered, so a single bounding box = the whole sheet;
  instead it reads **one small window per polygon** (a few % of a 206-MP scan). No Motzag
  polygons → notify and skip the sheet, don't abort.
- **Don't vectorise with RasterToPolyline.** It shatters the connected skeleton into thousands of
  ~4-px features whose endpoints never coincide (0 degree-2 nodes) → unmergeable. The skeleton-graph
  tracer (`trace_lines.py`) collapses each junction to one shared node instead.
- **Continuity = directional endpoint-linking** (`bridge_gaps.py`): bridge a gap only when the two
  ends are close AND aligned, so broken borders rejoin without fusing parallel neighbours. Re-running
  `skimage.skeletonize` afterwards is forbidden — it is non-idempotent and shatters a thin skeleton.
- **Junk is separable by connectivity.** Legend / margin number-table / symbols are short ISOLATED
  clusters; real borders form big CONNECTED networks. A `min_comp_len≈20 m` filter drops ~90% of junk
  automatically. (Corridor-tightening fails — it clips boundary-coincident real borders first.)
- **Basic license blocks a lot:** `UnsplitLine`, `SimplifyLine`, `Integrate`, `FeatureToPolygon`,
  all Spatial Analyst (ERROR 000824). Use `Dissolve(UNSPLIT_LINES)`, `ogr.SimplifyPreserveTopology`,
  numpy/scipy/GDAL, and `shapely` (polygonize/conflation) instead.
- **Conflation is per-point on the UNIONED reference network, guarded by a tangent test.** The
  modern cadastre was redrawn, so it can guide but not replace the extraction. Two simpler designs
  failed: naive per-point snapping staircases a line that crosses perpendicular dividers ("teeth"),
  and per-WHOLE-LINE replacement (one boundary ring must own ≥ match_frac of the line) fails on
  tiled parcels — no single ring owns a line that follows a boundary *chain*, and one divergence
  near a junction rejects the whole line. `snap_modern` therefore `unary_union`s all reference
  boundaries into ONE noded network (ring ownership disappears), projects each densified point onto
  the nearest edge running roughly **parallel** there (|cos| ≥ 0.5 — kills the teeth), and snaps
  per contiguous **run** (≥ `min_run` m), splitting a line into snapped + ink parts at the
  transitions. (Red-ink projects tol ≈ 3 m; coarse 1:20k black-ink ≈ 8 m.)
- **scikit-image is not in the ArcGIS python env** — run skeletonize/bridge from a throwaway venv made
  with `--system-site-packages`, and **keep numpy < 2** (base GDAL's `gdal_array` is numpy-1.x ABI).
- **CRS recovery:** GDAL surfaces no projection (georef is polynomial/GCP via `.tfwx`); the real CRS
  (EPSG 6991, Israeli Grid 05/12) lives in `<stem>.tif.aux.xml` and is read by `trace_lines`;
  `extract_9a` also stamps it directly on its windowed mask.
- GDAL Python gotcha: hold the datasource reference — `gdal.Open(p).GetRasterBand(1)` /
  `ogr.Open(p).GetLayer()` get garbage-collected mid-call.

## Environment
Two interpreters (see `requirements.txt`):
1. **ArcGIS Pro python** (`arcgispro-py3`) — numpy, scipy, GDAL/OGR, PIL, arcpy, shapely. Runs
   everything except skeletonize and bridge.
2. **venv** with `scikit-image` (numpy<2), created `--system-site-packages` off the ArcGIS python —
   runs `skeletonize.py` and `bridge_gaps.py`.

## Configuration
**`config.py`** is the single source of truth: interpreter paths, project roots, the parcels GDB,
corridor GDBs / layers, and per-project GUI defaults. Detection tuning constants
(`CAP_BLACK`/`CAP_RED`/`PURPLE_DB`/`BLOCK`/`OFFSET`/`MIN_PX`/`BUF_M`) live at the top of
`extract_9a.py`.
