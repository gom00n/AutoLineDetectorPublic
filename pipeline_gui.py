"""LineDetector pipeline GUI (Tkinter -- ships with the ArcGIS Python, no install).

A form-and-buttons front end over the CLI scripts. It runs each stage as a subprocess with the
correct interpreter (ArcGIS base vs the scikit-image venv), chained, streaming all output into the
log pane. Pick a project mode, a sheet <stem>, tweak parameters, choose which stages to run.

Run it with the ArcGIS base python:
  "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" pipeline_gui.py
"""
import os, re, sys, subprocess, threading, queue, tempfile, glob, shutil
import tkinter as tk
from tkinter import ttk, scrolledtext, filedialog, messagebox

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config                                            # central paths + per-project defaults

HERE = config.HERE
BASE_PY, VENV_PY = config.BASE_PY, config.VENV_PY
MODES = config.PROJECTS
STAGES = ["extract", "skeletonize", "bridge", "trace", "merge", "filter",
          "dehair", "desymbol", "collapse", "snap", "parcels", "overlay"]
OPTIONAL = ("dehair", "desymbol", "collapse", "snap", "parcels")   # off by default
NOISE = ("code page", "Code page", "codepage", "FutureWarning", "warnings.warn",
         "memory leak", "swig/python")

HELP = {
    "mode": "Project preset. ALL projects use the same Motzag-GDB engine and differ only by ink and "
            "default folders. black-ink = BLACK ink (global brightness threshold). red-ink = RED ink "
            "over a wash (adaptive local threshold). The target area is always the GDB polygons whose "
            "Motzag == the scan number. Switching mode reloads that preset's defaults.",
    "stem": "The sheet basename(s), e.g. M8_9076. The number after 'M' is the Motzag value that selects "
            "this scan's target polygons from the GDB. Use 'Pick files…' to tick sheets from the Input folder.",
    "srcdir": "INPUT folder: where the <stem>.tif scans (+ .tfwx / .tif.aux.xml) live. Any drive/path, "
              "including a P: network folder. 'Scan project' and 'Pick files…' list the sheets found here.",
    "outdir": "OUTPUT folder: each sheet's results go to <output>\\<stem>\\ — <stem>_lines.shp, "
              "<stem>_lines_snapped.shp, <stem>_preview.png and <stem>.log.",
    "base": "ArcGIS Pro Python (arcgispro-py3). Runs every stage EXCEPT skeletonize and bridge.",
    "venv": "The scikit-image virtual-env Python. Runs skeletonize and bridge_gaps (which need skimage).",
    # parameters
    "cap": "Ink brightness ceiling. BLACK ink: a pixel is ink if max(R,G,B) < cap (~145). RED ink: the "
           "absolute Red ceiling (~212) on the adaptive threshold. Lower = stricter; higher = catches "
           "fainter strokes but risks grabbing the coloured wash.",
    "offset": "RED ink only: how much darker (in DN) than its local background a pixel must be to count as "
              "ink. Higher = stricter (only strong lines); lower = catches faint strokes but more wash noise.",
    "block": "RED ink only: window size (px) used to estimate the local background for the adaptive "
             "threshold — roughly the spacing between lines. Larger = smoother background. Ignored for black ink.",
    "min_px": "Smallest connected blob (in pixels) to keep in the extracted mask. Drops speckle, dust and "
              "digit fragments. Black ~30; red ~400 (finer scans carry more speckle). Larger = cleaner.",
    "buf": "Buffer (metres) grown around the target (Motzag) polygons before clipping. Catches border lines "
           "sitting just outside the polygon edge. Bigger = more edge borders, but more surrounding junk.",
    "close_r": "Gap-closing radius (px) before thinning: bridges hairline ink breaks so a line isn't split. "
               "0 = off.",
    "max_fill": "Un-fill guard for gap-closing: a connected patch the closing ADDED that is bigger than this "
                "many px is reverted, so DOUBLED parallel casings (double-line tracks) stay two separate lines "
                "while small ink breaks stay repaired. 0 = off (close everything).",
    "stroke_len": "De-symbol: max length (m) of a symbol stroke (tick / X-mark between double casings).",
    "row_gap": "De-symbol: max spacing (m) between consecutive strokes for them to count as one symbol row "
               "(3+ collinear strokes are removed together).",
    "ring_perim": "De-symbol: closed rings up to this perimeter (m) with a fat interior are traced glyph "
                  "circles/blobs and get removed.",
    "prune": "Spur-prune length (px): removes short 'hair' branches left on the skeleton after thinning rough "
             "or thick strokes.",
    "gap": "Bridge max gap (px): the largest break the bridge step will span to rejoin two border ends — only "
           "if they are close AND aligned. Bigger spans dashed lines but risks false joins.",
    "cos": "How aligned two ends must be to bridge them (cosine of the angle; 0.80 ≈ 37°). Higher = stricter, "
           "i.e. more parallel, fewer wrong joins.",
    "spur": "Trace-stage hair removal: drop branches shorter than this (px) that have one free end and meet a "
            "junction.",
    "loop": "Drop small skeleton loops shorter than this (px) — bubbles from thinning doubled lines. Caution: "
            "very small real parcels are also small loops.",
    "simp_px": "Straightening tolerance (px) inside the tracer: removes the pixel staircase so lines run "
               "straight between corners.",
    "min_comp": "Junk filter: keep a line only if its connected network totals at least this length (metres). "
                "Legend / number-table / symbols are small isolated clusters; real borders form big networks.",
    "simp_m": "Final straighten tolerance (metres) — Douglas-Peucker simplification to drop the staircase (~0.5).",
    "hair": "De-hair: remove short barbs up to this length (metres) that have one free end and one junction "
            "end. Real borders and long dangles are kept.",
    "passes": "How many de-hair passes to run (removing one barb can expose the next).",
    "snap_tol": "Snap tolerance (metres): a stretch of a line moves onto a reference boundary only while it "
                "stays within this distance of it. red-ink ≈ 3 m; 1:20k black-ink ≈ 8 m (coarser offset).",
    "snap_frac": "SHORT lines only (shorter than min_run): fraction of the line that must lie within tol "
                 "before the whole line snaps (0.6 = 60%). Longer lines are split per stretch instead.",
    "snap_run": "Minimum stretch (metres) that must follow a boundary before that stretch is snapped. A line "
                "can be split: the following stretch snaps (snapped=1), the diverging rest stays ink "
                "(snapped=0). Larger = more conservative, fewer but surer snaps.",
    "pair_d": "Collapse: max separation (metres) of a DOUBLED line pair (road drawn as two casings). Pairs "
              "running mutually parallel within this distance become ONE centerline and their short 'rungs' "
              "are dropped. Keep well below real parcel widths (default 3).",
    "ext_m": "Parcels: extend a loose line end straight ahead up to this many metres to meet the border it "
             "stops short of. Bigger closes more parcels but risks wrong connections.",
    "join_m": "Parcels: weld two facing loose ends when they are within this many metres of each other "
              "(closes near-meeting corners).",
    "min_area": "Parcels: drop closed polygons smaller than this (m²) — slivers from doubled lines and "
                "skeleton bubbles, not real parcels.",
    "parcels_gdb": "A GDB holding the clean cadastral PARCEL polygons to snap onto (e.g. Cad-GOV\\Layers.gdb). "
                   "Browse to it, then pick the layer(s) in 'Parcels lyr'. Snap stage only.",
    "parcels_layers": "Which layer(s) in the Parcels GDB to snap onto (e.g. PARCEL_ALL). Use 'Choose…'. These "
                      "are pooled with the 'Snap layers' boundaries; a line snaps to whichever it best follows.",
    "gdb": "The geodatabase holding the corridor polygons (required — every project is Motzag-based). "
           "The code reads only the layers you choose and, in them, the polygons with Motzag == the scan "
           "number, to define the detect zone. Use 'Choose…' to tick the layers.",
    "layers": "Which GDB feature layers to scan for the Motzag polygons (comma-separated). Use 'Choose…' "
              "to tick them from the GDB. Only the buffered Motzag==scan polygons are read from the scan.",
    "snap_layers": "Which layer(s) IN THE CORRIDOR GDB to snap detected lines onto — the block/Gush polygons "
                   "(e.g. Poly_ID, Poly_Map, Poly_map2). Pooled with the Parcels GDB layers. The Snap stage "
                   "replaces a detected line by a boundary it closely follows. Use 'Choose…' to pick from the GDB.",
    # stages
    "stage_extract": "Stage 1: turn the RGB scan into a binary mask of just the target ink (black global "
                     "threshold or red adaptive), clipped to the buffered Motzag corridor read from the GDB.",
    "stage_skeletonize": "Thin the ink mask to 1-pixel centrelines (runs in the scikit-image venv).",
    "stage_bridge": "Rejoin broken border ends across small gaps (close + aligned) to improve continuity.",
    "stage_trace": "Vectorise the skeleton into clean polylines with shared junction nodes "
                   "(replaces RasterToPolyline).",
    "stage_merge": "Merge touching segments into longer continuous lines (Dissolve UNSPLIT_LINES).",
    "stage_filter": "Drop junk by connected-component size and straighten. Produces <stem>_lines.shp.",
    "stage_dehair": "Optional: remove leftover short hair barbs from the final lines.",
    "stage_desymbol": "Optional: remove traced map-SYMBOL ink -- periodic tick/X rows along double "
                      "tracks, glyph rings, isolated compact ornament clusters. Review the result: "
                      "validated on M14 but symbol rows can rarely overlap real dense corners.",
    "stage_collapse": "Optional: collapse DOUBLED parallel lines (thick roads traced as two casings + rungs "
                      "= a 'ladder') into one centerline. Runs before Snap so the centerline gets snapped. "
                      "Check the QA overlay after enabling.",
    "stage_parcels": "Optional: close near-meeting borders (join/extend loose ends) and build parcel "
                     "POLYGONS via shapely polygonize (FeatureToPolygon needs Standard — this doesn't). "
                     "Writes <stem>_parcels.shp with area_m2/perim_m; the lines outputs are not modified.",
    "stage_snap": "Optional: conflate lines onto the reference-boundary NETWORK pooled from the 'Snap layers' "
                  "(corridor GDB) + 'Parcels lyr' (parcels GDB). Stretches following a boundary within tol for "
                  ">= min_run metres are replaced by it; a line can be SPLIT into a snapped part and an ink "
                  "part. Writes <stem>_lines_snapped.shp with a 'snapped' field (1 = snapped, 0 = kept as "
                  "ink). Divergent historical lines are left as drawn.",
    "stage_overlay": "Render a QA PNG over the faded scan. After a snap run it colors by result: GREEN = "
                     "snapped to a reference boundary, RED = kept as ink.",
}


class Tooltip:
    """Lightweight hover tooltip for any Tk widget."""
    def __init__(self, widget, text):
        self.widget, self.text, self.tip = widget, text, None
        if text:
            widget.bind("<Enter>", self.show, add="+")
            widget.bind("<Leave>", self.hide, add="+")

    def show(self, _):
        if self.tip or not self.text:
            return
        x = self.widget.winfo_rootx() + 18
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 3
        self.tip = tk.Toplevel(self.widget)
        self.tip.wm_overrideredirect(True)
        self.tip.wm_geometry("+%d+%d" % (x, y))
        tk.Label(self.tip, text=self.text, justify="left", background="#ffffe0",
                 relief="solid", borderwidth=1, wraplength=380, font=("", 9)).pack()

    def hide(self, _):
        if self.tip:
            self.tip.destroy(); self.tip = None


class App:
    def __init__(self, root):
        self.root = root
        root.title("LineDetector pipeline")
        self.q = queue.Queue()
        self.proc = None
        self.worker = None
        self.stop_flag = False
        self.cur_stem = ""
        self.tmp = tempfile.mkdtemp(prefix="ld_gui_")
        self.v = {}

        top = ttk.Frame(root, padding=8); top.pack(fill="x")
        ttk.Label(top, text="Mode").grid(row=0, column=0, sticky="w")
        self.mode = tk.StringVar(value=list(MODES)[0])
        mb = ttk.Combobox(top, textvariable=self.mode, values=list(MODES), state="readonly", width=28)
        mb.grid(row=0, column=1, sticky="w"); mb.bind("<<ComboboxSelected>>", lambda e: self.load_defaults())
        Tooltip(mb, HELP["mode"])
        sl = ttk.Label(top, text="Sheets"); sl.grid(row=0, column=2, sticky="e", padx=(12, 2))
        self.stem = tk.StringVar(value="")
        se = ttk.Entry(top, textvariable=self.stem, width=40); se.grid(row=0, column=3, sticky="we")
        ttk.Button(top, text="Scan", command=self.scan_project).grid(row=0, column=4, padx=2)
        ttk.Button(top, text="Pick files…", command=self.pick_files).grid(row=0, column=5, padx=2)
        sheets_help = ("One or more sheets to process, looked up in the Input folder. Accepts: a single "
                       "stem (M8_9076); a comma/space list; wildcards (M2_*); a FOLDER path; or 'all'. "
                       "'Scan' fills this with every sheet in the Input folder; 'Pick files…' lets you tick them.")
        Tooltip(sl, sheets_help); Tooltip(se, sheets_help)
        ttk.Label(top, text="(hover any field for help)", foreground="#777").grid(row=1, column=3, sticky="w")
        top.columnconfigure(3, weight=1)

        # interpreters (fixed paths from config; no UI needed)
        self.base = tk.StringVar(value=BASE_PY); self.venv = tk.StringVar(value=VENV_PY)

        # input / output folders (overridable per run; defaults come from the mode preset)
        ff = ttk.LabelFrame(root, text="Folders", padding=6); ff.pack(fill="x", padx=8, pady=(4, 0))
        self.v["srcdir"] = tk.StringVar(); self.v["outdir"] = tk.StringVar()
        il = ttk.Label(ff, text="Input (scans)"); il.grid(row=0, column=0, sticky="w")
        ie = ttk.Entry(ff, textvariable=self.v["srcdir"], width=66); ie.grid(row=0, column=1, sticky="we")
        ttk.Button(ff, text="…", width=3, command=self.pick_srcdir).grid(row=0, column=2, padx=2)
        ol = ttk.Label(ff, text="Output (results)"); ol.grid(row=1, column=0, sticky="w")
        oe = ttk.Entry(ff, textvariable=self.v["outdir"], width=66); oe.grid(row=1, column=1, sticky="we")
        ttk.Button(ff, text="…", width=3, command=self.pick_outdir).grid(row=1, column=2, padx=2)
        Tooltip(il, HELP["srcdir"]); Tooltip(ie, HELP["srcdir"])
        Tooltip(ol, HELP["outdir"]); Tooltip(oe, HELP["outdir"])
        ff.columnconfigure(1, weight=1)

        # parameters
        pf = ttk.LabelFrame(root, text="Parameters", padding=6); pf.pack(fill="x", padx=8, pady=4)
        specs = [("Detection", [("cap", "CAP"), ("offset", "offset"), ("block", "block"),
                                ("min_px", "MIN_PX"), ("buf", "buffer m")]),
                 ("Skeleton", [("close_r", "close_r"), ("prune", "prune_len"), ("max_fill", "max_fill")]),
                 ("Bridge", [("gap", "max_gap"), ("cos", "cos")]),
                 ("Trace", [("spur", "min_spur"), ("loop", "min_loop"), ("simp_px", "simplify_px")]),
                 ("Filter", [("min_comp", "min_comp_len"), ("simp_m", "simplify_m")]),
                 ("De-hair", [("hair", "max_hair"), ("passes", "passes")]),
                 ("De-symbol", [("stroke_len", "stroke_len"), ("row_gap", "row_gap"),
                                ("ring_perim", "ring_perim")]),
                 ("Snap", [("snap_tol", "tol"), ("snap_frac", "match_frac"), ("snap_run", "min_run")]),
                 ("Collapse", [("pair_d", "pair_d")]),
                 ("Parcels", [("ext_m", "extend_m"), ("join_m", "join_m"), ("min_area", "min_area")])]
        r = 0
        for group, fields in specs:
            ttk.Label(pf, text=group, font=("", 9, "bold")).grid(row=r, column=0, sticky="w", pady=(4, 0))
            c = 1
            for key, lab in fields:
                lw = ttk.Label(pf, text=lab); lw.grid(row=r, column=c, sticky="e"); c += 1
                self.v[key] = tk.StringVar()
                ew = ttk.Entry(pf, textvariable=self.v[key], width=7); ew.grid(row=r, column=c, sticky="w"); c += 1
                Tooltip(lw, HELP.get(key, "")); Tooltip(ew, HELP.get(key, ""))
            r += 1
        gl = ttk.Label(pf, text="Corridor GDB"); gl.grid(row=r, column=0, sticky="w")
        self.v["gdb"] = tk.StringVar()
        ge = ttk.Entry(pf, textvariable=self.v["gdb"], width=60); ge.grid(row=r, column=1, columnspan=5, sticky="we")
        ttk.Button(pf, text="…", width=3, command=self.pick_gdb).grid(row=r, column=6, sticky="w")
        Tooltip(gl, HELP["gdb"]); Tooltip(ge, HELP["gdb"]); r += 1
        ll = ttk.Label(pf, text="GDB layers"); ll.grid(row=r, column=0, sticky="w")
        self.v["layers"] = tk.StringVar()
        le = ttk.Entry(pf, textvariable=self.v["layers"], width=52); le.grid(row=r, column=1, columnspan=4, sticky="we")
        ttk.Button(pf, text="Choose…", command=self.choose_layers).grid(row=r, column=5, columnspan=2, sticky="w")
        Tooltip(ll, HELP["layers"]); Tooltip(le, HELP["layers"]); r += 1
        sll = ttk.Label(pf, text="Snap layers"); sll.grid(row=r, column=0, sticky="w")
        self.v["snap_layers"] = tk.StringVar()
        sle = ttk.Entry(pf, textvariable=self.v["snap_layers"], width=52); sle.grid(row=r, column=1, columnspan=4, sticky="we")
        ttk.Button(pf, text="Choose…", command=lambda: self.choose_layers("snap_layers")).grid(row=r, column=5, columnspan=2, sticky="w")
        Tooltip(sll, HELP["snap_layers"]); Tooltip(sle, HELP["snap_layers"]); r += 1
        pgl = ttk.Label(pf, text="Parcels GDB"); pgl.grid(row=r, column=0, sticky="w")
        self.v["parcels_gdb"] = tk.StringVar()
        pge = ttk.Entry(pf, textvariable=self.v["parcels_gdb"], width=60); pge.grid(row=r, column=1, columnspan=5, sticky="we")
        ttk.Button(pf, text="…", width=3, command=lambda: self.pick_gdb("parcels_gdb")).grid(row=r, column=6, sticky="w")
        Tooltip(pgl, HELP["parcels_gdb"]); Tooltip(pge, HELP["parcels_gdb"]); r += 1
        pll = ttk.Label(pf, text="Parcels lyr"); pll.grid(row=r, column=0, sticky="w")
        self.v["parcels_layers"] = tk.StringVar()
        ple = ttk.Entry(pf, textvariable=self.v["parcels_layers"], width=52); ple.grid(row=r, column=1, columnspan=4, sticky="we")
        ttk.Button(pf, text="Choose…", command=lambda: self.choose_layers("parcels_layers", "parcels_gdb")).grid(row=r, column=5, columnspan=2, sticky="w")
        Tooltip(pll, HELP["parcels_layers"]); Tooltip(ple, HELP["parcels_layers"])

        # stage checkboxes
        sf = ttk.LabelFrame(root, text="Stages to run", padding=6); sf.pack(fill="x", padx=8)
        self.stage = {}
        for i, s in enumerate(STAGES):
            self.stage[s] = tk.BooleanVar(value=(s not in OPTIONAL))
            cb = ttk.Checkbutton(sf, text=s, variable=self.stage[s]); cb.grid(row=0, column=i, padx=3)
            Tooltip(cb, HELP.get("stage_" + s, ""))

        # actions
        af = ttk.Frame(root, padding=6); af.pack(fill="x")
        self.run_btn = ttk.Button(af, text="▶ Run", command=self.start); self.run_btn.pack(side="left")
        self.stop_btn = ttk.Button(af, text="■ Stop", command=self.stop, state="disabled"); self.stop_btn.pack(side="left", padx=4)
        ob = ttk.Button(af, text="Open output folder", command=self.open_out); ob.pack(side="left", padx=4)
        vb = ttk.Button(af, text="View overlay", command=self.view_overlay); vb.pack(side="left")
        mb2 = ttk.Button(af, text="Metrics…", command=self.run_metrics); mb2.pack(side="left", padx=4)
        self.keep = tk.BooleanVar(value=False)
        kc = ttk.Checkbutton(af, text="keep intermediates", variable=self.keep); kc.pack(side="left", padx=10)
        Tooltip(kc, "Off (default): intermediate mask/skeleton TIFs go to a temp folder, so the work "
                    "folder keeps only <stem>_lines.shp, <stem>_lines_snapped.shp, <stem>_preview.png "
                    "and <stem>.log. On: write them into the work folder (handy when re-running single "
                    "stages while tuning).")
        Tooltip(self.run_btn, "Run the ticked stages in order, with the right interpreter per stage.")
        Tooltip(self.stop_btn, "Terminate the currently running stage.")
        Tooltip(ob, "Open the sheet's working folder (the outputs land there).")
        Tooltip(vb, "Open the QA overlay PNG for this sheet (run the overlay stage first).")
        Tooltip(mb2, "Score this sheet's result against a hand-cleaned REFERENCE shapefile (you pick it): "
                     "prints recall / precision / F1 to the log. Match tolerance = the Snap 'tol' value. "
                     "Reference must be in the same CRS (EPSG 6991). Uses _lines_snapped.shp if present, else _lines.shp.")
        ttk.Button(af, text="Clear log", command=lambda: self.log_w.delete("1.0", "end")).pack(side="right")
        self.prog = tk.StringVar(value="idle")
        ttk.Label(af, textvariable=self.prog, foreground="#080").pack(side="right", padx=10)

        self.log_w = scrolledtext.ScrolledText(root, height=18, font=("Consolas", 9))
        self.log_w.pack(fill="both", expand=True, padx=8, pady=(0, 8))

        self.load_defaults()
        self.root.after(80, self.drain)
        root.protocol("WM_DELETE_WINDOW", self.on_close)

    # ---- helpers ----
    def load_defaults(self):
        d = MODES[self.mode.get()]["defaults"]
        for k, val in d.items():
            self.v[k].set(val)

    def pick_gdb(self, target="gdb"):
        p = filedialog.askdirectory(title="Select a .gdb folder")
        if p:
            self.v[target].set(p)

    def choose_layers(self, target="layers", gdb_key="gdb"):
        """List the polygon layers in the chosen GDB with checkboxes; OK fills the <target> field."""
        from osgeo import ogr
        gdb = self.v[gdb_key].get().strip()
        ds = ogr.Open(gdb)
        if ds is None:
            messagebox.showerror("GDB", "Cannot open:\n" + gdb); return
        polys = [ds.GetLayerByIndex(i).GetName() for i in range(ds.GetLayerCount())
                 if "Polygon" in ogr.GeometryTypeToName(ds.GetLayerByIndex(i).GetGeomType())]
        if not polys:
            messagebox.showinfo("GDB", "No polygon layers found in this GDB."); return
        cur = set(x.strip() for x in self.v[target].get().split(",") if x.strip())
        prompt = ("Tick the cadastral PARCEL layer(s) to snap onto:" if target == "snap_layers"
                  else "Tick the layers to scan for Motzag polygons:")
        win = tk.Toplevel(self.root); win.title("Select GDB layers"); win.transient(self.root); win.grab_set()
        ttk.Label(win, text=prompt, padding=8).pack(anchor="w")
        cbvars = {}
        for name in polys:
            cbvars[name] = tk.BooleanVar(value=(name in cur))
            ttk.Checkbutton(win, text=name, variable=cbvars[name]).pack(anchor="w", padx=14, pady=1)

        def ok():
            self.v[target].set(",".join(n for n, v in cbvars.items() if v.get()))
            win.destroy()
        ttk.Button(win, text="OK", command=ok).pack(pady=8)

    def pick_srcdir(self):
        p = filedialog.askdirectory(title="Input folder (where the scans are)")
        if p:
            self.v["srcdir"].set(p)

    def pick_outdir(self):
        p = filedialog.askdirectory(title="Output folder (where results go)")
        if p:
            self.v["outdir"].set(p)

    def pick_files(self):
        """Tick specific sheets found in the Input folder; fills the Sheets field."""
        src = self.v["srcdir"].get().strip(); ext = MODES[self.mode.get()]["src_ext"]
        sheets = config.list_sheets(src, ext)
        if not sheets:
            messagebox.showinfo("Input folder", "No %s sheets (M<n>_<n>) found in:\n%s" % (ext, src or "(empty)")); return
        cur = set(config.resolve_stems(self.stem.get(), src, ext))
        win = tk.Toplevel(self.root); win.title("Select sheets"); win.transient(self.root); win.grab_set()
        ttk.Label(win, text="Tick the sheets to process:", padding=8).pack(anchor="w")
        grid = ttk.Frame(win); grid.pack(padx=10)
        cbvars = {}
        for i, name in enumerate(sheets):
            cbvars[name] = tk.BooleanVar(value=(name in cur))
            ttk.Checkbutton(grid, text=name, variable=cbvars[name]).grid(row=i % 14, column=i // 14, sticky="w", padx=8)

        def ok():
            self.stem.set(", ".join(n for n, v in cbvars.items() if v.get()))
            win.destroy()
        ttk.Button(win, text="OK", command=ok).pack(pady=8)

    def scan_project(self):
        src = self.v["srcdir"].get().strip(); ext = MODES[self.mode.get()]["src_ext"]
        self.stem.set(", ".join(config.list_sheets(src, ext)) or "(none found)")

    def workdir(self, stem=None):
        s = stem or self.cur_stem or self.stem.get().strip()
        return os.path.join(self.v["outdir"].get().strip(), s)

    def _ref_stem(self):
        if self.cur_stem:
            return self.cur_stem
        src = self.v["srcdir"].get().strip(); ext = MODES[self.mode.get()]["src_ext"]
        stems = config.resolve_stems(self.stem.get(), src, ext)
        return stems[0] if stems else self.stem.get().strip()

    def log(self, line):
        self.log_w.insert("end", line + "\n"); self.log_w.see("end")

    def drain(self):
        try:
            while True:
                self.log(self.q.get_nowait())
        except queue.Empty:
            pass
        self.root.after(80, self.drain)

    # ---- pipeline ----
    def build(self, stem):
        m = MODES[self.mode.get()]; s = stem; g = self.v
        wd = self.workdir(stem); ex = lambda f: os.path.join(HERE, f)
        idir = wd if self.keep.get() else self.tmp           # intermediate TIFs: temp by default
        mask = os.path.join(idir, s + "_mask.tif")
        skel = os.path.join(idir, s + "_skel.tif")
        brd = os.path.join(idir, s + "_skel_bridged.tif")
        final = os.path.join(wd, s + "_lines.shp")           # short, clear names
        snapped = os.path.join(wd, s + "_lines_snapped.shp")
        parcels = os.path.join(wd, s + "_parcels.shp")
        tr, mg, dh, dy, cl = (os.path.join(self.tmp, s + x)
                              for x in ("_tr.shp", "_mg.shp", "_dh.shp", "_dy.shp", "_cl.shp"))
        BASE, VENV = self.base.get(), self.venv.get()
        get = lambda k: g[k].get().strip()
        srcdir = get("srcdir")

        cmds = []  # (label, interpreter, argv)
        if self.stage["extract"].get():                      # unified Motzag extractor (extract_9a)
            argv = [ex(m["extract"]), s, "--ink", m["ink"], "--srcdir", srcdir,
                    "--cap", get("cap"), "--offset", get("offset"), "--block", get("block"),
                    "--min-px", get("min_px"), "--buf", get("buf"), "--out", mask]
            if get("gdb"):
                argv += ["--gdb", get("gdb")]
            if get("layers"):
                argv += ["--layers", get("layers")]
            cmds.append(("extract", BASE, argv))
        if self.stage["skeletonize"].get():
            cmds.append(("skeletonize", VENV, [ex("skeletonize.py"), mask, skel, get("close_r"), get("prune"),
                                               get("max_fill") or "0"]))
        if self.stage["bridge"].get():
            cmds.append(("bridge", VENV, [ex("bridge_gaps.py"), skel, brd, get("gap"), get("cos"), "1"]))
        if self.stage["trace"].get():
            cmds.append(("trace", BASE, [ex("trace_lines.py"), brd, tr, get("spur"), get("loop"), get("simp_px")]))
        if self.stage["merge"].get():
            cmds.append(("merge", BASE, [ex("merge_lines.py"), tr, mg, "0"]))
        if self.stage["filter"].get():
            cmds.append(("filter", BASE, [ex("filter_components.py"), mg, final, get("min_comp"), get("simp_m")]))
        snap_in = final
        if self.stage["dehair"].get():
            cmds.append(("dehair", BASE, [ex("dehair_lines.py"), final, dh, get("hair"), get("passes")]))
            snap_in = dh
        if self.stage["desymbol"].get():                 # traced map-symbol ink (tick rows, glyphs)
            cmds.append(("desymbol", BASE, [ex("desymbol.py"), snap_in, dy,
                                            "--stroke-len", get("stroke_len") or "12",
                                            "--row-gap", get("row_gap") or "30",
                                            "--ring-perim", get("ring_perim") or "45"]))
            snap_in = dy
        if self.stage["collapse"].get():                 # doubled road casings -> one centerline
            cmds.append(("collapse", BASE, [ex("collapse_pairs.py"), snap_in, cl,
                                            "--pair-d", get("pair_d") or "3.0"]))
            snap_in = cl
        if self.stage["snap"].get():
            sargv = [ex("snap_modern.py"), snap_in, snapped, "--tol", get("snap_tol"),
                     "--match-frac", get("snap_frac"), "--min-run", get("snap_run") or "10"]
            if get("gdb") and get("snap_layers"):            # block/Poly_* layers from the corridor GDB
                sargv += ["--snap-src", get("gdb"), get("snap_layers")]
            if get("parcels_gdb") and get("parcels_layers"):  # parcels from a separate GDB
                sargv += ["--snap-src", get("parcels_gdb"), get("parcels_layers")]
            cmds.append(("snap", BASE, sargv))
        if self.stage["parcels"].get():                  # close loose ends + polygonize
            parcels_in = snapped if self.stage["snap"].get() else snap_in
            cmds.append(("parcels", BASE, [ex("close_parcels.py"), parcels_in, parcels,
                                           "--ext", get("ext_m") or "5", "--join", get("join_m") or "3",
                                           "--min-area", get("min_area") or "30"]))
        if self.stage["overlay"].get():
            cmds.append(("overlay", BASE, [ex(m["overlay"]), s, "--srcdir", srcdir, "--workdir", wd]))
        return cmds

    def _preflight(self, stems):
        """Cheap checks BEFORE a long run: interpreters, source sheets, every named GDB layer.
        Catches the silent failures (typo'd layer name -> 0 reference rings after a 10-minute
        run). Returns a list of problem strings; empty = good to go."""
        from osgeo import ogr
        probs = []
        for what, p in (("ArcGIS base python", self.base.get()), ("skimage venv python", self.venv.get())):
            if not os.path.exists(p):
                probs.append("%s not found: %s" % (what, p))
        src = self.v["srcdir"].get().strip()
        if self.stage["extract"].get():
            for s in stems:
                if not os.path.exists(os.path.join(src, s + ".tif")):
                    probs.append("scan not found: %s.tif in %s" % (s, src or "(no Input folder)"))
                elif not os.path.exists(os.path.join(src, s + ".tfwx")):
                    probs.append("world file missing: %s.tfwx in %s" % (s, src))
            if not self.v["layers"].get().strip():
                probs.append("Corridor GDB layers not set (use 'Choose…' next to GDB layers)")
        checks = []                                   # (label, gdb path, comma-separated layers)
        if self.stage["extract"].get():
            checks.append(("Corridor GDB", self.v["gdb"].get().strip(), self.v["layers"].get()))
        if self.stage["snap"].get():
            if self.v["snap_layers"].get().strip():
                checks.append(("Corridor GDB (Snap layers)", self.v["gdb"].get().strip(),
                               self.v["snap_layers"].get()))
            if self.v["parcels_gdb"].get().strip() and self.v["parcels_layers"].get().strip():
                checks.append(("Parcels GDB", self.v["parcels_gdb"].get().strip(),
                               self.v["parcels_layers"].get()))
        for what, gdb, layers in checks:
            if not gdb:
                probs.append("%s: no GDB path set" % what); continue
            ds = ogr.Open(gdb)
            if ds is None:
                probs.append("%s: cannot open %s" % (what, gdb)); continue
            for name in (x.strip() for x in layers.split(",") if x.strip()):
                if ds.GetLayerByName(name) is None:
                    probs.append("%s: layer '%s' not found in %s" % (what, name, os.path.basename(gdb)))
            ds = None
        return probs

    def start(self):
        if self.worker and self.worker.is_alive():
            return
        src = self.v["srcdir"].get().strip(); ext = MODES[self.mode.get()]["src_ext"]
        stems = config.resolve_stems(self.stem.get(), src, ext)
        if not stems:
            messagebox.showinfo("No sheets", "Enter a stem, a list, a folder, or 'all'.\n"
                                "'Scan' / 'Pick files…' list the sheets in the Input folder."); return
        if not any(self.stage[s].get() for s in STAGES):
            messagebox.showinfo("Nothing to run", "Tick at least one stage."); return
        probs = self._preflight(stems)
        if probs:
            messagebox.showerror("Pre-run check failed", "Fix these before running:\n\n- "
                                 + "\n- ".join(probs)); return
        if len(stems) > 1 and not messagebox.askyesno(
                "Batch", "Process %d sheets?\n\n%s" % (len(stems),
                ", ".join(stems[:12]) + (" ..." if len(stems) > 12 else ""))):
            return
        self.stop_flag = False
        self.run_btn["state"] = "disabled"; self.stop_btn["state"] = "normal"
        self.worker = threading.Thread(target=self.run_all, args=(stems,), daemon=True); self.worker.start()

    def emit(self, line, lf):
        self.q.put(line)
        if lf:
            try:
                lf.write(line + "\n"); lf.flush()
            except Exception:
                pass

    def run_all(self, stems):
        import datetime
        n = len(stems)
        for i, stem in enumerate(stems, 1):
            if self.stop_flag:
                break
            self.cur_stem = stem
            wd = self.workdir(stem); os.makedirs(wd, exist_ok=True)
            try:
                lf = open(os.path.join(wd, stem + ".log"), "w", encoding="utf-8")
                lf.write("# %s  %s\n# mode=%s  params=%s\n" % (
                    stem, datetime.datetime.now().isoformat(timespec="seconds"),
                    self.mode.get(), {k: self.v[k].get() for k in self.v}))
            except Exception:
                lf = None
            self.q.put("\n########## [%d/%d] %s ##########" % (i, n, stem))
            self.root.after(0, lambda i=i, s=stem: self.prog.set("sheet %d/%d: %s" % (i, n, s)))
            try:
                cmds = self.build(stem)
            except Exception as e:
                self.emit("!! build error: %s" % e, lf)
                if lf:
                    lf.close()
                continue
            self._run_cmds(cmds, lf)
            if not self.keep.get() and self.stage["filter"].get():   # full run: keep only shp/png/log
                removed = 0
                # any intermediate that landed in the work folder ...
                junk = [os.path.join(wd, stem + suf) for suf in
                        ("_mask.tif", "_skel.tif", "_skel_bridged.tif")]
                # ... plus this sheet's temp intermediates, so they don't pile up across a batch
                for patt in ("_mask.tif", "_skel.tif", "_skel_bridged.tif", "_tr.*", "_mg.*", "_dh.*", "_dy.*", "_cl.*"):
                    junk += glob.glob(os.path.join(self.tmp, stem + patt))
                for fp in junk:
                    try:
                        if os.path.exists(fp):
                            os.remove(fp); removed += 1
                    except Exception:
                        pass
                if removed:
                    self.emit("[clean] removed %d intermediate file(s) (work folder + temp)" % removed, lf)
            if lf:
                lf.close()
                self.q.put("[log] %s.log" % stem)
        self.proc = None
        self.root.after(0, lambda: self.prog.set("done" if not self.stop_flag else "stopped"))
        self.q.put("\n[all done]")
        self.root.after(0, self.finish)

    def _run_cmds(self, cmds, lf):
        for label, interp, argv in cmds:
            if self.stop_flag:
                return
            if not os.path.exists(interp):
                self.emit("!! interpreter not found: %s" % interp, lf); return
            self.emit("\n===== %s =====" % label, lf)
            self.emit("$ %s %s" % (os.path.basename(interp),
                      " ".join('"%s"' % a if " " in a else a for a in argv[1:])), lf)
            try:
                p = subprocess.Popen([interp] + argv, cwd=HERE, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, text=True, bufsize=1)
            except Exception as e:
                self.emit("!! failed to launch: %s" % e, lf); return
            self.proc = p
            for line in p.stdout:
                line = line.rstrip()
                if line and not any(nz in line for nz in NOISE):
                    self.emit(line, lf)
            p.wait()
            if label == "extract" and p.returncode == 3:
                self.emit("-- skipped (no target area / missing sheet); continuing to next sheet.", lf)
                return
            if p.returncode != 0:
                self.emit("!! %s exited with code %d -- skipping rest of this sheet." % (label, p.returncode), lf)
                return

    def finish(self):
        self.run_btn["state"] = "normal"; self.stop_btn["state"] = "disabled"

    def stop(self):
        self.stop_flag = True
        if self.proc and self.proc.poll() is None:
            self.proc.terminate()
        self.q.put("!! stopping...")

    def on_close(self):
        """Stop any running stage and remove the session temp dir before closing."""
        self.stop_flag = True
        if self.proc and self.proc.poll() is None:
            try:
                self.proc.terminate()
            except Exception:
                pass
        shutil.rmtree(self.tmp, ignore_errors=True)
        self.root.destroy()

    def open_out(self):
        wd = self.workdir(self._ref_stem())
        if os.path.isdir(wd):
            os.startfile(wd)
        else:
            messagebox.showinfo("No folder", wd + "\n(run the pipeline first)")

    def view_overlay(self):
        s = self._ref_stem(); wd = self.workdir(s)
        p = os.path.join(wd, s + "_preview.png")
        if os.path.exists(p):
            os.startfile(p); return
        messagebox.showinfo("No overlay", "Run the 'overlay' stage first for %s." % s)

    def run_metrics(self):
        """Score this sheet's result vs a hand-cleaned reference: recall / precision / F1 into the log."""
        if self.worker and self.worker.is_alive():
            messagebox.showinfo("Busy", "Wait for the current run to finish, then try Metrics."); return
        s = self._ref_stem(); wd = self.workdir(s)
        snapped = os.path.join(wd, s + "_lines_snapped.shp")
        res = snapped if os.path.exists(snapped) else os.path.join(wd, s + "_lines.shp")
        if not os.path.exists(res):
            messagebox.showinfo("No result", "No result lines for %s yet -- run the pipeline first." % s); return
        ref = filedialog.askopenfilename(title="Hand-cleaned REFERENCE shapefile for %s (EPSG 6991)" % s,
                                         filetypes=[("Shapefile", "*.shp")])
        if not ref:
            return
        tol = self.v["snap_tol"].get().strip() or "2.0"
        argv = [os.path.join(HERE, "metrics.py"), res, ref, tol]
        self.q.put("\n===== metrics: %s  vs  %s   (tol=%s m) =====" %
                   (os.path.basename(res), os.path.basename(ref), tol))

        def work():
            try:
                p = subprocess.Popen([self.base.get()] + argv, cwd=HERE, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, text=True, bufsize=1)
            except Exception as e:
                self.q.put("!! metrics failed to launch: %s" % e); return
            for line in p.stdout:
                line = line.rstrip()
                if line and not any(nz in line for nz in NOISE):
                    self.q.put(line)
            p.wait()
        threading.Thread(target=work, daemon=True).start()


if __name__ == "__main__":
    root = tk.Tk()
    try:
        root.call("tk", "scaling", 1.2)
    except tk.TclError:
        pass
    App(root)
    root.mainloop()
