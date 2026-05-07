"""
Side-by-side viewer for original vs. standard-setting unit cells.
Shows equivalent Miller planes, d-spacings, 2-theta, and predicted intensity.

Usage:
    python cell_compare.py [path/to/file.cif]
    (or double-click the .exe)
"""

import sys
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import gemmi
import spglib
import itertools

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

WAVELENGTHS = {
    "Cu Ka1": 1.54056,
    "Cu Ka2": 1.54439,
    "Mo Ka1": 0.70930,
    "Mo Ka2": 0.71359,
    "Co Ka1": 1.78897,
    "Ag Ka1": 0.55941,
    "Cr Ka1": 2.28970,
}
DEFAULT_WAVELENGTH = "Cu Ka1"

# Colour palette
BG          = "#1e1e2e"
BG_PANEL    = "#252536"
BG_HEADER   = "#2d2d44"
BG_ENTRY    = "#33334d"
FG          = "#e0e0e0"
FG_DIM      = "#888899"
FG_ACCENT   = "#7aa2f7"
FG_STRONG   = "#ffffff"
FG_ABSENT   = "#f7768e"
FG_WEAK     = "#636379"
HIGHLIGHT   = "#3d59a1"
ROW_EVEN    = "#252536"
ROW_ODD     = "#2a2a3e"
BORDER      = "#3b3b55"
BTN_BG      = "#33334d"
BTN_FG      = "#c0caf5"
BTN_ACTIVE  = "#414168"

# ---------------------------------------------------------------------------
# CIF loading
# ---------------------------------------------------------------------------

def load_cif(path):
    doc = gemmi.cif.read(path)
    block_names = [block.name for block in doc]
    if not block_names:
        raise ValueError("CIF file contains no data blocks.")
    if len(block_names) > 1:
        chosen = choose_block(block_names)
        if chosen is None:
            return None
    else:
        chosen = block_names[0]
    st = gemmi.read_small_structure(path)
    if st.cell.a == 0:
        raise ValueError("Could not read cell parameters from CIF.")
    return st, chosen


def choose_block(names):
    result = [None]
    win = tk.Toplevel()
    win.title("Choose data block")
    win.configure(bg=BG)
    win.resizable(False, False)
    tk.Label(win, text="Multiple data blocks found.\nSelect one:",
             font=("Segoe UI", 12), bg=BG, fg=FG, padx=16, pady=12).pack()
    lb = tk.Listbox(win, font=("Consolas", 12), height=min(len(names), 10),
                    bg=BG_ENTRY, fg=FG, selectbackground=HIGHLIGHT,
                    selectforeground=FG_STRONG, relief=tk.FLAT, bd=0,
                    highlightthickness=1, highlightcolor=BORDER)
    for n in names:
        lb.insert(tk.END, n)
    lb.selection_set(0)
    lb.pack(padx=16, pady=6)

    def on_ok():
        sel = lb.curselection()
        if sel:
            result[0] = names[sel[0]]
        win.destroy()

    tk.Button(win, text="OK", command=on_ok, font=("Segoe UI", 11),
              width=10, bg=BTN_BG, fg=BTN_FG, activebackground=BTN_ACTIVE,
              relief=tk.FLAT, cursor="hand2").pack(pady=12)
    win.grab_set()
    win.wait_window()
    return result[0]


# ---------------------------------------------------------------------------
# Lattice math
# ---------------------------------------------------------------------------

def lattice_matrix(a, b, c, alpha, beta, gamma):
    al, be, ga = np.radians(alpha), np.radians(beta), np.radians(gamma)
    bx = b * np.cos(ga)
    by = b * np.sin(ga)
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = np.sqrt(max(c**2 - cx**2 - cy**2, 0.0))
    return np.array([[a, 0, 0], [bx, by, 0], [cx, cy, cz]])


def cell_params_from_lattice(L):
    av, bv, cv = L[0], L[1], L[2]
    a = np.linalg.norm(av)
    b = np.linalg.norm(bv)
    c = np.linalg.norm(cv)
    alpha = np.degrees(np.arccos(np.clip(np.dot(bv, cv) / (b * c), -1, 1)))
    beta  = np.degrees(np.arccos(np.clip(np.dot(av, cv) / (a * c), -1, 1)))
    gamma = np.degrees(np.arccos(np.clip(np.dot(av, bv) / (a * b), -1, 1)))
    return a, b, c, alpha, beta, gamma


def reciprocal_metric_tensor(a, b, c, alpha, beta, gamma):
    al, be, ga = np.radians(alpha), np.radians(beta), np.radians(gamma)
    V = a * b * c * np.sqrt(
        1 - np.cos(al)**2 - np.cos(be)**2 - np.cos(ga)**2
        + 2 * np.cos(al) * np.cos(be) * np.cos(ga))
    sa, sb, sg = np.sin(al), np.sin(be), np.sin(ga)
    a_s = b * c * sa / V
    b_s = a * c * sb / V
    c_s = a * b * sg / V
    cas = (np.cos(be) * np.cos(ga) - np.cos(al)) / (sb * sg)
    cbs = (np.cos(al) * np.cos(ga) - np.cos(be)) / (sa * sg)
    cgs = (np.cos(al) * np.cos(be) - np.cos(ga)) / (sa * sb)
    return np.array([
        [a_s**2,          a_s * b_s * cgs, a_s * c_s * cbs],
        [a_s * b_s * cgs, b_s**2,          b_s * c_s * cas],
        [a_s * c_s * cbs, b_s * c_s * cas, c_s**2]])


def d_spacing(hkl, G_star):
    h = np.array(hkl, dtype=float)
    inv_d2 = h @ G_star @ h
    if inv_d2 <= 1e-12:
        return float('inf')
    return 1.0 / np.sqrt(inv_d2)


def two_theta(d, wavelength):
    if d <= 0 or d == float('inf'):
        return None
    ratio = wavelength / (2 * d)
    if abs(ratio) >= 1.0:
        return None
    return 2 * np.degrees(np.arcsin(ratio))


# ---------------------------------------------------------------------------
# Cell standardisation via spglib
# ---------------------------------------------------------------------------

def expand_asu(st):
    positions = []
    numbers = []
    ops = st.spacegroup.operations()
    for site in st.sites:
        z = gemmi.Element(site.type_symbol).atomic_number
        for op in ops:
            pos = op.apply_to_xyz([site.fract.x, site.fract.y, site.fract.z])
            pos = [p % 1.0 for p in pos]
            dup = False
            for ex in positions:
                diff = [min(abs(pos[i] - ex[i]), 1 - abs(pos[i] - ex[i]))
                        for i in range(3)]
                if all(d < 0.02 for d in diff):
                    dup = True
                    break
            if not dup:
                positions.append(pos)
                numbers.append(z)
    return positions, numbers


def standardize_cell(st):
    cell = st.cell
    L = lattice_matrix(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    positions, numbers = expand_asu(st)
    cell_tuple = (L, positions, numbers)
    ds = None
    for sp in [1e-3, 0.01, 0.05, 0.1, 0.5]:
        ds = spglib.get_symmetry_dataset(cell_tuple, symprec=sp)
        if ds is not None and ds.number > 1:
            break
    if ds is None or ds.number <= 1:
        return None
    std_params = cell_params_from_lattice(ds.std_lattice)
    P = ds.transformation_matrix
    std_sg = ds.international if hasattr(ds, 'international') else str(ds.number)
    return std_params, P, std_sg


# ---------------------------------------------------------------------------
# HKL helpers
# ---------------------------------------------------------------------------

def transform_hkl(hkl, P):
    result = P.T @ np.array(hkl, dtype=float)
    return tuple(int(round(x)) for x in result)


def laue_equivalents(hkl, sg):
    equivs = set()
    h_arr = list(hkl)
    for op in sg.operations():
        h_new = op.apply_to_hkl(h_arr)
        equivs.add(tuple(h_new))
        equivs.add(tuple(-x for x in h_new))
    return equivs


def canonical_hkl(equivs):
    return max(equivs, key=lambda hkl: (hkl[0], hkl[1], hkl[2]))


def generate_merged_reflections(max_idx, sg, G_star, sf_calc, st, wavelength):
    seen = set()
    reflections = []
    for h, k, l in itertools.product(range(-max_idx, max_idx + 1), repeat=3):
        if h == 0 and k == 0 and l == 0:
            continue
        hkl = (h, k, l)
        if hkl in seen:
            continue
        equivs = laue_equivalents(hkl, sg)
        for eq in equivs:
            seen.add(eq)
        rep = canonical_hkl(equivs)
        d = d_spacing(rep, G_star)
        tt = two_theta(d, wavelength)
        sf = sf_calc.calculate_sf_from_small_structure(st, list(rep))
        intensity = abs(sf) ** 2
        reflections.append({
            'hkl': rep, 'equivs': equivs, 'd': d,
            'two_theta': tt, 'intensity': intensity,
            'multiplicity': len(equivs)})
    reflections.sort(key=lambda r: (r['two_theta'] is None, r['two_theta'] or 999))
    return reflections


# ---------------------------------------------------------------------------
# Build paired table data
# ---------------------------------------------------------------------------

def build_table(st, max_idx, wavelength):
    cell = st.cell
    orig_params = (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    orig_sg = st.spacegroup
    G_star_orig = reciprocal_metric_tensor(*orig_params)
    sf_calc = gemmi.StructureFactorCalculatorX(cell)

    std_result = standardize_cell(st)
    if std_result is None:
        std_params, P, std_sg_name = orig_params, np.eye(3), orig_sg.hm
    else:
        std_params, P, std_sg_name = std_result

    G_star_std = reciprocal_metric_tensor(*std_params)
    orig_refls = generate_merged_reflections(
        max_idx, orig_sg, G_star_orig, sf_calc, st, wavelength)

    rows = []
    for refl in orig_refls:
        hkl_orig = refl['hkl']
        hkl_std_raw = transform_hkl(hkl_orig, P)
        neg = tuple(-x for x in hkl_std_raw)
        hkl_std = neg if neg > hkl_std_raw else hkl_std_raw
        rows.append({
            'hkl_orig': hkl_orig, 'hkl_std': hkl_std,
            'd_orig': refl['d'], 'd_std': d_spacing(hkl_std, G_star_std),
            'two_theta_orig': refl['two_theta'],
            'two_theta_std': two_theta(d_spacing(hkl_std, G_star_std), wavelength),
            'intensity': refl['intensity'],
            'multiplicity': refl['multiplicity']})

    max_int = max((r['intensity'] for r in rows), default=1.0)
    for r in rows:
        r['intensity_norm'] = round(100 * r['intensity'] / max_int) if max_int > 0 else 0

    info = {'orig_params': orig_params, 'std_params': std_params,
            'orig_sg': orig_sg.hm, 'std_sg': std_sg_name, 'P': P}
    return rows, info


# ---------------------------------------------------------------------------
# Tkinter UI
# ---------------------------------------------------------------------------

def fmt_hkl(hkl):
    return f"({hkl[0]:2d} {hkl[1]:2d} {hkl[2]:2d})"

def fmt_intensity(raw, norm, threshold=0.5):
    if raw < threshold:
        return "✗  0"
    return str(norm)

def parse_hkl(text):
    text = text.strip().replace("(", "").replace(")", "").replace(",", " ")
    parts = text.split()
    if len(parts) != 3:
        return None
    try:
        return tuple(int(x) for x in parts)
    except ValueError:
        return None


class CellCompareApp:
    ABSENT_THRESHOLD = 0.5

    def __init__(self, root, cif_path=None):
        self.root = root
        self.root.title("Cell Compare")
        self.root.geometry("1440x860")
        self.root.minsize(1000, 520)
        self.root.configure(bg=BG)

        self.st = None
        self.rows = []
        self.info = {}
        self.font_size = 14
        self.wavelength_var = tk.StringVar(value=DEFAULT_WAVELENGTH)
        self.max_idx_var = tk.IntVar(value=1)

        self._apply_theme()
        self._build_ui()

        if cif_path and os.path.isfile(cif_path):
            self._load_file(cif_path)

    # ---- Theme ----

    def _apply_theme(self):
        style = ttk.Style()
        style.theme_use("clam")

        style.configure("Treeview",
                        background=ROW_EVEN, foreground=FG,
                        fieldbackground=ROW_EVEN,
                        font=("Consolas", self.font_size),
                        rowheight=self.font_size + 14,
                        borderwidth=0, relief=tk.FLAT)
        style.configure("Treeview.Heading",
                        background=BG_HEADER, foreground=FG_ACCENT,
                        font=("Segoe UI Semibold", self.font_size - 1),
                        borderwidth=0, relief=tk.FLAT, padding=(0, 6))
        style.map("Treeview.Heading",
                  background=[("active", BG_HEADER)])
        style.map("Treeview",
                  background=[("selected", HIGHLIGHT)],
                  foreground=[("selected", FG_STRONG)])

        style.configure("TScrollbar",
                        background=BG_PANEL, troughcolor=BG,
                        borderwidth=0, arrowsize=14)
        style.map("TScrollbar",
                  background=[("active", BTN_ACTIVE), ("!disabled", BG_PANEL)])

        style.configure("TCombobox",
                        fieldbackground=BG_ENTRY, background=BG_ENTRY,
                        foreground=FG, arrowcolor=FG_DIM,
                        borderwidth=1, relief=tk.FLAT)
        style.map("TCombobox",
                  fieldbackground=[("readonly", BG_ENTRY)],
                  selectbackground=[("readonly", BG_ENTRY)],
                  selectforeground=[("readonly", FG)])

    # ---- UI construction ----

    def _build_ui(self):
        # ---- toolbar ----
        toolbar = tk.Frame(self.root, bg=BG_HEADER, pady=6)
        toolbar.pack(fill=tk.X)

        self._make_btn(toolbar, "Open CIF", self._open_file).pack(
            side=tk.LEFT, padx=(12, 6))

        sep1 = tk.Frame(toolbar, width=1, bg=BORDER)
        sep1.pack(side=tk.LEFT, fill=tk.Y, padx=8, pady=2)

        tk.Label(toolbar, text="max |h|,|k|,|l|", font=("Segoe UI", 11),
                 bg=BG_HEADER, fg=FG_DIM).pack(side=tk.LEFT, padx=(4, 4))
        spin = tk.Spinbox(toolbar, from_=1, to=10, width=3,
                          textvariable=self.max_idx_var,
                          font=("Consolas", 12), command=self._recompute,
                          bg=BG_ENTRY, fg=FG, buttonbackground=BTN_BG,
                          insertbackground=FG, relief=tk.FLAT, bd=2,
                          highlightthickness=0)
        spin.pack(side=tk.LEFT)
        spin.bind("<Return>", lambda e: self._recompute())

        sep2 = tk.Frame(toolbar, width=1, bg=BORDER)
        sep2.pack(side=tk.LEFT, fill=tk.Y, padx=8, pady=2)

        tk.Label(toolbar, text="Source", font=("Segoe UI", 11),
                 bg=BG_HEADER, fg=FG_DIM).pack(side=tk.LEFT, padx=(4, 4))
        wl_menu = ttk.Combobox(toolbar, textvariable=self.wavelength_var,
                               values=list(WAVELENGTHS.keys()),
                               state="readonly", width=9,
                               font=("Segoe UI", 11))
        wl_menu.pack(side=tk.LEFT)
        wl_menu.bind("<<ComboboxSelected>>", lambda e: self._recompute())

        self._make_btn(toolbar, " + ", self._zoom_in).pack(side=tk.RIGHT, padx=(2, 12))
        self._make_btn(toolbar, " – ", self._zoom_out).pack(side=tk.RIGHT, padx=2)

        # ---- dual header ----
        hdr_frame = tk.Frame(self.root, bg=BG)
        hdr_frame.pack(fill=tk.X, padx=0)
        hdr_frame.columnconfigure(0, weight=1)
        hdr_frame.columnconfigure(1, weight=1)

        self.orig_hdr = self._make_cell_header(hdr_frame, 0, "original")
        self.std_hdr  = self._make_cell_header(hdr_frame, 1, "standard")

        # ---- search row ----
        search_row = tk.Frame(self.root, bg=BG_PANEL, pady=4)
        search_row.pack(fill=tk.X)
        search_row.columnconfigure(0, weight=1)
        search_row.columnconfigure(1, weight=1)

        self.search_left  = self._make_search(search_row, 0, "left")
        self.search_right = self._make_search(search_row, 1, "right")

        # ---- tables ----
        table_frame = tk.Frame(self.root, bg=BG)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=0, pady=0)
        table_frame.columnconfigure(0, weight=1)
        table_frame.columnconfigure(1, weight=1)
        table_frame.rowconfigure(0, weight=1)

        self.left_tree, self.left_scroll = self._make_tree(table_frame, 0)
        self.right_tree, self.right_scroll = self._make_tree(table_frame, 1)

        self.left_tree.configure(yscrollcommand=self._sync_left)
        self.right_tree.configure(yscrollcommand=self._sync_right)
        self.left_scroll.configure(command=self._scroll_both)
        self.right_scroll.configure(command=self._scroll_both)

        # ---- status bar ----
        self.status = tk.Label(self.root, text="No file loaded. Click Open CIF.",
                               font=("Segoe UI", 10), bg=BG_HEADER, fg=FG_DIM,
                               anchor="w", padx=12, pady=3)
        self.status.pack(fill=tk.X, side=tk.BOTTOM)

        self.root.bind("<Escape>", lambda e: self._clear_search())
        # Drag-and-drop support (basic: accepts command-line paths)
        self.root.drop_target_register = lambda *a: None  # placeholder

    def _make_btn(self, parent, text, cmd):
        return tk.Button(parent, text=text, command=cmd,
                         font=("Segoe UI", 11), bg=BTN_BG, fg=BTN_FG,
                         activebackground=BTN_ACTIVE, activeforeground=FG_STRONG,
                         relief=tk.FLAT, cursor="hand2", padx=10, pady=2,
                         highlightthickness=0, bd=0)

    def _make_cell_header(self, parent, col, label):
        frame = tk.Frame(parent, bg=BG_PANEL, padx=16, pady=8,
                         highlightbackground=BORDER, highlightthickness=1)
        frame.grid(row=0, column=col, sticky="nsew", padx=(0 if col else 0, 0))
        lbl_sg = tk.Label(frame, text="", font=("Segoe UI Semibold", 13),
                          bg=BG_PANEL, fg=FG_ACCENT, anchor="w")
        lbl_sg.pack(fill=tk.X)
        lbl_cell = tk.Label(frame, text="", font=("Consolas", 12),
                            bg=BG_PANEL, fg=FG, anchor="w")
        lbl_cell.pack(fill=tk.X)
        return lbl_sg, lbl_cell

    def _make_search(self, parent, col, side):
        frame = tk.Frame(parent, bg=BG_PANEL)
        frame.grid(row=0, column=col, sticky="ew", padx=12, pady=2)
        tk.Label(frame, text="Search (h k l):", font=("Segoe UI", 11),
                 bg=BG_PANEL, fg=FG_DIM).pack(side=tk.LEFT)
        entry = tk.Entry(frame, font=("Consolas", 13), width=12,
                         bg=BG_ENTRY, fg=FG, insertbackground=FG_ACCENT,
                         relief=tk.FLAT, bd=4, highlightthickness=1,
                         highlightcolor=FG_ACCENT, highlightbackground=BORDER)
        entry.pack(side=tk.LEFT, padx=6)
        entry.bind("<KeyRelease>", lambda e, s=side: self._on_search(s))
        self._make_btn(frame, "Clear", self._clear_search).pack(side=tk.LEFT, padx=2)
        return entry

    def _make_tree(self, parent, col):
        frame = tk.Frame(parent, bg=BG)
        frame.grid(row=0, column=col, sticky="nsew", padx=(0, 1 if col == 0 else 0))

        cols = ("hkl", "d", "2theta", "I")
        tree = ttk.Treeview(frame, columns=cols, show="headings", selectmode="browse")
        tree.heading("hkl", text="(h k l)")
        tree.heading("d", text="d (Å)")
        tree.heading("2theta", text="2θ (°)")
        tree.heading("I", text="I")
        tree.column("hkl", width=110, anchor="center", minwidth=80)
        tree.column("d", width=100, anchor="center", minwidth=70)
        tree.column("2theta", width=100, anchor="center", minwidth=70)
        tree.column("I", width=65, anchor="center", minwidth=50)

        scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self._configure_tags(tree)
        return tree, scroll

    def _configure_tags(self, tree):
        tree.tag_configure("strong",
                           foreground=FG_STRONG, background=ROW_EVEN,
                           font=("Consolas", self.font_size, "bold"))
        tree.tag_configure("medium",
                           foreground=FG, background=ROW_EVEN,
                           font=("Consolas", self.font_size))
        tree.tag_configure("weak",
                           foreground=FG_WEAK, background=ROW_EVEN,
                           font=("Consolas", self.font_size))
        tree.tag_configure("absent",
                           foreground=FG_ABSENT, background=ROW_EVEN,
                           font=("Consolas", self.font_size))
        tree.tag_configure("strong_alt",
                           foreground=FG_STRONG, background=ROW_ODD,
                           font=("Consolas", self.font_size, "bold"))
        tree.tag_configure("medium_alt",
                           foreground=FG, background=ROW_ODD,
                           font=("Consolas", self.font_size))
        tree.tag_configure("weak_alt",
                           foreground=FG_WEAK, background=ROW_ODD,
                           font=("Consolas", self.font_size))
        tree.tag_configure("absent_alt",
                           foreground=FG_ABSENT, background=ROW_ODD,
                           font=("Consolas", self.font_size))
        tree.tag_configure("highlight",
                           background=HIGHLIGHT, foreground=FG_STRONG)

    # ---- Scroll sync ----

    def _sync_left(self, *args):
        self.left_scroll.set(*args)
        self.right_tree.yview_moveto(args[0])

    def _sync_right(self, *args):
        self.right_scroll.set(*args)
        self.left_tree.yview_moveto(args[0])

    def _scroll_both(self, *args):
        self.left_tree.yview(*args)
        self.right_tree.yview(*args)

    # ---- File handling ----

    def _open_file(self):
        path = filedialog.askopenfilename(
            title="Select a CIF file",
            filetypes=[("CIF files", "*.cif"), ("All files", "*.*")])
        if path:
            self._load_file(path)

    def _load_file(self, path):
        try:
            result = load_cif(path)
            if result is None:
                return
            self.st, block = result
            self.root.title(f"Cell Compare — {os.path.basename(path)}")
            self.status.config(text=f"Loaded: {os.path.basename(path)}  [{block}]")
            self._recompute()
        except Exception as e:
            messagebox.showerror("Error loading CIF", str(e))

    # ---- Computation ----

    def _recompute(self):
        if self.st is None:
            return
        max_idx = self.max_idx_var.get()
        wl = WAVELENGTHS[self.wavelength_var.get()]
        self.root.config(cursor="watch")
        self.root.update_idletasks()
        try:
            self.rows, self.info = build_table(self.st, max_idx, wl)
        except Exception as e:
            self.root.config(cursor="")
            messagebox.showerror("Computation error", str(e))
            return
        self._update_headers()
        self._populate_trees()
        self._clear_search()
        self.status.config(
            text=f"{len(self.rows)} reflections  |  "
                 f"max |hkl|={max_idx}  |  "
                 f"λ={wl:.5f} Å ({self.wavelength_var.get()})")
        self.root.config(cursor="")

    def _update_headers(self):
        def fmt_sg(sg, label):
            return f"{sg}  ({label})"
        def fmt_cell(params):
            a, b, c, al, be, ga = params
            return (f"a = {a:.4f}   b = {b:.4f}   c = {c:.4f}   "
                    f"α = {al:.2f}°   β = {be:.2f}°   "
                    f"γ = {ga:.2f}°")

        self.orig_hdr[0].config(text=fmt_sg(self.info['orig_sg'], "original"))
        self.orig_hdr[1].config(text=fmt_cell(self.info['orig_params']))
        self.std_hdr[0].config(text=fmt_sg(self.info['std_sg'], "standard"))
        self.std_hdr[1].config(text=fmt_cell(self.info['std_params']))

    def _populate_trees(self):
        for tree in (self.left_tree, self.right_tree):
            tree.delete(*tree.get_children())

        for i, row in enumerate(self.rows):
            base_tag = self._intensity_tag(row['intensity'], row['intensity_norm'])
            tag = base_tag + ("_alt" if i % 2 else "")

            d_o = f"{row['d_orig']:.4f}" if row['d_orig'] < 1e6 else "—"
            tt_o = f"{row['two_theta_orig']:.3f}" if row['two_theta_orig'] else "—"
            d_s = f"{row['d_std']:.4f}" if row['d_std'] < 1e6 else "—"
            tt_s = f"{row['two_theta_std']:.3f}" if row['two_theta_std'] else "—"
            i_txt = fmt_intensity(row['intensity'], row['intensity_norm'])

            iid = str(i)
            self.left_tree.insert("", tk.END, iid=iid,
                values=(fmt_hkl(row['hkl_orig']), d_o, tt_o, i_txt), tags=(tag,))
            self.right_tree.insert("", tk.END, iid=iid,
                values=(fmt_hkl(row['hkl_std']), d_s, tt_s, i_txt), tags=(tag,))

    def _intensity_tag(self, raw, norm):
        if raw < self.ABSENT_THRESHOLD:
            return "absent"
        if norm >= 90:
            return "strong"
        if norm >= 20:
            return "medium"
        return "weak"

    # ---- Search ----

    def _on_search(self, side):
        entry = self.search_left if side == "left" else self.search_right
        other = self.search_right if side == "left" else self.search_left
        text = entry.get().strip()
        if not text:
            self._clear_highlight()
            other.config(state=tk.NORMAL)
            return
        other.delete(0, tk.END)
        other.config(state=tk.DISABLED)
        hkl = parse_hkl(text)
        if hkl is None:
            self._clear_highlight()
            return
        self._highlight_hkl(hkl, side)

    def _highlight_hkl(self, hkl, side):
        self._clear_highlight()
        target = hkl
        neg = tuple(-x for x in hkl)

        # Direct / Friedel match
        for i, row in enumerate(self.rows):
            key = 'hkl_orig' if side == "left" else 'hkl_std'
            if row[key] == target or row[key] == neg:
                self._set_highlight(i)
                return

        # Laue-equivalent search
        if self.st:
            sg = self.st.spacegroup
            if side == "left":
                equivs = laue_equivalents(hkl, sg)
            else:
                P = self.info.get('P', np.eye(3))
                hkl_orig = tuple(int(round(x))
                                 for x in np.linalg.inv(P).T @ np.array(hkl))
                equivs = laue_equivalents(hkl_orig, sg)

            for i, row in enumerate(self.rows):
                h = row['hkl_orig']
                if h in equivs or tuple(-x for x in h) in equivs:
                    self._set_highlight(i)
                    return

    def _set_highlight(self, idx):
        iid = str(idx)
        for tree in (self.left_tree, self.right_tree):
            tree.item(iid, tags=("highlight",))
            tree.see(iid)

    def _clear_highlight(self):
        for i, row in enumerate(self.rows):
            base_tag = self._intensity_tag(row['intensity'], row['intensity_norm'])
            tag = base_tag + ("_alt" if i % 2 else "")
            iid = str(i)
            for tree in (self.left_tree, self.right_tree):
                tree.item(iid, tags=(tag,))

    def _clear_search(self, event=None):
        self.search_left.config(state=tk.NORMAL)
        self.search_right.config(state=tk.NORMAL)
        self.search_left.delete(0, tk.END)
        self.search_right.delete(0, tk.END)
        self._clear_highlight()

    # ---- Zoom ----

    def _zoom_in(self):
        self.font_size = min(self.font_size + 2, 28)
        self._apply_font_size()

    def _zoom_out(self):
        self.font_size = max(self.font_size - 2, 8)
        self._apply_font_size()

    def _apply_font_size(self):
        style = ttk.Style()
        style.configure("Treeview",
                        font=("Consolas", self.font_size),
                        rowheight=self.font_size + 14)
        style.configure("Treeview.Heading",
                        font=("Segoe UI Semibold", self.font_size - 1))
        for tree in (self.left_tree, self.right_tree):
            self._configure_tags(tree)
        # Re-tag to refresh
        self._clear_highlight()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    cif_path = sys.argv[1] if len(sys.argv) > 1 else None
    root = tk.Tk()
    CellCompareApp(root, cif_path)
    root.mainloop()


if __name__ == "__main__":
    main()
