# -*- coding: utf-8 -*-
"""
PBI (Propagation-Based Imaging) Phase Retrieval GUI tab.

Supports:
 - Multi-distance TIE phase retrieval
 - Multi-distance CTF phase retrieval
 - Single-distance Paganin retrieval (homogeneous-material)

Includes live preview of loaded TIFF images and result maps with W/L controls.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter.filedialog import asksaveasfilename
import os
import numpy as np
import tifffile
import threading
import json
import io
import zipfile
import datetime
import traceback
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib import rcParams

from GUI.ui.styles import Styles as stl
from GUI.ui.widgets import Widget as wg, DropZone, VerticalScrolledFrame as vsf, ToggleButton
from GUI.ui.tooltips import ToolTip
from GUI.utils import resource_path
import GUI.i18n as i18n
from src.PBIRec.phase_retrieval import (
    normalize_image,
    paganin_single_distance,
    tie_multi_distance,
    ctf_multi_distance,
)


ALGO_OPTIONS = ["TIE 多距离", "CTF 多距离", "Paganin 单距离"]

# Shared dark-theme matplotlib params
_MPL_DARK = {
    "text.color": "white",
    "xtick.color": "white",
    "ytick.color": "white",
    "axes.grid": False,
    "axes.labelcolor": "white",
    "axes.edgecolor": "#666666",
}


class PBIRec_GUI:
    """GUI panel for PBI (inline) phase retrieval from experimental data."""

    def __init__(self, master, status_var=None):
        self.master = master
        self.status_var = status_var

        stl.configure_style()

        # --- Internal state ---
        self.intensity_images = []       # list of 2-D arrays (raw)
        self.intensity_paths = []        # display names
        self.flat_image = None
        self.dark_image = None
        self.phase_result = None
        self.absorption_result = None

        # Preview tracking
        self._preview_img_data = None
        self._preview_im_handle = None
        self._result_handles = {}        # {"phase": im, "abs": im}

        # --- Tk variables ---
        self.energy_var = tk.DoubleVar(value=30.0)
        self.pixel_var = tk.DoubleVar(value=50.0)
        self.dso_var = tk.DoubleVar(value=100.0)
        self.algo_var = tk.StringVar(value=ALGO_OPTIONS[0])
        self.delta_beta_var = tk.DoubleVar(value=1000.0)
        self.alpha_var = tk.DoubleVar(value=1e-3)
        self.pad_var = tk.IntVar(value=0)
        self.beam_var = tk.StringVar(value="Conical")
        self.zip_var = tk.BooleanVar(value=False)

        self._build_ui()

    # ==================================================================
    #  UI construction
    # ==================================================================
    def _build_ui(self):
        self.master.grid_columnconfigure(0, weight=0, minsize=340)
        self.master.grid_columnconfigure(1, weight=1)
        self.master.grid_rowconfigure(0, weight=1)

        # ---- Left panel: parameters & file loading ----
        left = ttk.Frame(self.master, style="TFrame")
        left.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)
        left.grid_columnconfigure(0, weight=1)
        left.grid_columnconfigure(1, weight=1)

        row = 0

        # Drop-zone: intensity images (multi-file)
        self.img_dropzone = DropZone(
            left,
            title=i18n.PBI_DZ_IMG_TITLE,
            description=i18n.PBI_DZ_IMG_DESC,
            button_text=i18n.PBI_DZ_IMG_BTN,
            extensions=[".tif", ".tiff"],
            filetypes=[("TIFF images", "*.tif *.tiff")],
            multiple=True,
            on_files_dropped=self._on_images_loaded,
        )
        self.img_dropzone.grid(row=row, column=0, columnspan=2, sticky="ew",
                               padx=5, pady=(6, 4), ipady=2)
        row += 1

        # Drop-zone: flat
        self.flat_dropzone = DropZone(
            left,
            title=i18n.PBI_DZ_FLAT_TITLE,
            description=i18n.PBI_DZ_FLAT_DESC,
            button_text=i18n.PBI_DZ_FLAT_BTN,
            extensions=[".tif", ".tiff"],
            filetypes=[("TIFF images", "*.tif *.tiff")],
            multiple=False,
            on_files_dropped=self._on_flat_loaded,
        )
        self.flat_dropzone.grid(row=row, column=0, columnspan=2, sticky="ew",
                                padx=5, pady=(4, 4), ipady=2)
        row += 1

        # Drop-zone: dark
        self.dark_dropzone = DropZone(
            left,
            title=i18n.PBI_DZ_DARK_TITLE,
            description=i18n.PBI_DZ_DARK_DESC,
            button_text=i18n.PBI_DZ_DARK_BTN,
            extensions=[".tif", ".tiff"],
            filetypes=[("TIFF images", "*.tif *.tiff")],
            multiple=False,
            on_files_dropped=self._on_dark_loaded,
        )
        self.dark_dropzone.grid(row=row, column=0, columnspan=2, sticky="ew",
                                padx=5, pady=(4, 4), ipady=2)
        row += 1

        # --- Distance table ---
        dist_frame = ttk.LabelFrame(left, text=i18n.PBI_LBL_DIST_TABLE)
        dist_frame.grid(row=row, column=0, columnspan=2, sticky="ew", padx=5, pady=4)
        row += 1

        self.dist_text = tk.Text(dist_frame, height=4, width=38, bg="gray12",
                                 fg="white", font=("Consolas", 9), wrap=tk.WORD)
        self.dist_text.pack(fill="x", padx=4, pady=4)
        self.dist_text.insert("1.0", i18n.PBI_DIST_PLACEHOLDER)

        # --- Parameters ---
        lbl_e, _ = wg.create_label_entry(left, i18n.PBI_LBL_ENERGY, row, 0,
                                          textvariable=self.energy_var, padx=10)
        row += 1
        lbl_px, _ = wg.create_label_entry(left, i18n.PBI_LBL_PIXEL, row, 0,
                                           textvariable=self.pixel_var, padx=10)
        row += 1
        lbl_dso, _ = wg.create_label_entry(left, i18n.PBI_LBL_DSO, row, 0,
                                            textvariable=self.dso_var, padx=10)
        row += 1
        lbl_beam, _ = wg.create_label_combobox(left, i18n.PBI_LBL_BEAM,
                                                ["Conical", "Plane"], row, 0,
                                                textvariable=self.beam_var)
        row += 1
        lbl_algo, _ = wg.create_label_combobox(left, i18n.PBI_LBL_ALGO,
                                                ALGO_OPTIONS, row, 0,
                                                textvariable=self.algo_var)
        row += 1
        lbl_db, _ = wg.create_label_entry(left, i18n.PBI_LBL_DELTA_BETA, row, 0,
                                           textvariable=self.delta_beta_var, padx=10)
        row += 1
        lbl_alpha, _ = wg.create_label_entry(left, i18n.PBI_LBL_ALPHA, row, 0,
                                              textvariable=self.alpha_var, padx=10)
        row += 1
        lbl_pad, _ = wg.create_label_entry(left, i18n.PBI_LBL_PAD, row, 0,
                                            textvariable=self.pad_var, padx=10)
        row += 1

        ToggleButton(left, text=i18n.TOGGLE_CREATE_ZIP,
                     variable=self.zip_var).grid(row=row, column=0, pady=5)
        row += 1

        self.run_btn = wg.create_button(left, i18n.BTN_RUN, row, 0,
                                         command=self._run_retrieval)
        row += 1

        wg.create_button(left, i18n.PBI_BTN_SAVE_PHASE, row, 0,
                          command=lambda: self._save_image(self.phase_result, "phase"))
        wg.create_button(left, i18n.PBI_BTN_SAVE_ABS, row, 1,
                          command=lambda: self._save_image(self.absorption_result, "absorption"))
        row += 1

        # Tooltips
        ToolTip(lbl_e, i18n.PBI_TT_ENERGY)
        ToolTip(lbl_px, i18n.PBI_TT_PIXEL)
        ToolTip(lbl_dso, i18n.PBI_TT_DSO)
        ToolTip(lbl_beam, i18n.PBI_TT_BEAM)
        ToolTip(lbl_algo, i18n.PBI_TT_ALGO)
        ToolTip(lbl_db, i18n.PBI_TT_DELTA_BETA)
        ToolTip(lbl_alpha, i18n.PBI_TT_ALPHA)
        ToolTip(lbl_pad, i18n.PBI_TT_PAD)

        # ---- Right panel: preview + results ----
        self._build_right_panel()

    # ------------------------------------------------------------------
    def _build_right_panel(self):
        """Build the right panel with notebook tabs for input preview and results."""
        right = ttk.Frame(self.master, style="TFrame")
        right.grid(row=0, column=1, sticky="nsew", padx=4, pady=4)
        right.grid_rowconfigure(0, weight=1)
        right.grid_columnconfigure(0, weight=1)

        self.right_notebook = ttk.Notebook(right)
        self.right_notebook.grid(row=0, column=0, sticky="nsew")

        # --- Tab 1: Input preview ---
        self.preview_tab = ttk.Frame(self.right_notebook, style="TFrame")
        self.right_notebook.add(self.preview_tab, text=i18n.PBI_PREVIEW_TAB_INPUT)
        self._build_preview_tab()

        # --- Tab 2: Results ---
        self.results_tab = ttk.Frame(self.right_notebook, style="TFrame")
        self.right_notebook.add(self.results_tab, text=i18n.PBI_PREVIEW_TAB_RESULT)
        self._build_results_tab()

    # ------------------------------------------------------------------
    def _build_preview_tab(self):
        """Input-image preview tab with thumbnail list + big preview + toolbar."""
        self.preview_tab.grid_rowconfigure(1, weight=1)
        self.preview_tab.grid_columnconfigure(0, weight=0, minsize=140)
        self.preview_tab.grid_columnconfigure(1, weight=1)

        # --- Image list (left column) ---
        list_frame = ttk.LabelFrame(self.preview_tab, text=i18n.PBI_PREVIEW_TITLE)
        list_frame.grid(row=0, column=0, rowspan=3, sticky="nsew", padx=(4, 2), pady=4)
        list_frame.grid_rowconfigure(0, weight=1)
        list_frame.grid_columnconfigure(0, weight=1)

        self.img_listbox = tk.Listbox(
            list_frame, bg="#1e1e1e", fg="white", selectbackground="#0078d4",
            font=("Consolas", 9), width=22, activestyle="none",
        )
        self.img_listbox.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        self.img_listbox.bind("<<ListboxSelect>>", self._on_listbox_select)

        sb = ttk.Scrollbar(list_frame, orient="vertical",
                           command=self.img_listbox.yview)
        sb.grid(row=0, column=1, sticky="ns")
        self.img_listbox.configure(yscrollcommand=sb.set)

        # --- Image info label ---
        self.preview_info_var = tk.StringVar(value=i18n.PBI_PREVIEW_NONE)
        info_lbl = tk.Label(
            self.preview_tab, textvariable=self.preview_info_var,
            bg="#252525", fg="#cccccc", font=("Consolas", 9),
            anchor="w", justify="left",
        )
        info_lbl.grid(row=0, column=1, sticky="ew", padx=(2, 4), pady=(4, 0), ipady=2)

        # --- Big preview canvas ---
        self.preview_canvas_frame = ttk.Frame(self.preview_tab, style="TFrame")
        self.preview_canvas_frame.grid(row=1, column=1, sticky="nsew",
                                       padx=(2, 4), pady=(0, 4))
        self.preview_canvas_frame.grid_rowconfigure(0, weight=1)
        self.preview_canvas_frame.grid_columnconfigure(0, weight=1)

        self._init_preview_canvas()

    # ------------------------------------------------------------------
    def _build_results_tab(self):
        """Results tab with phase + absorption display and W/L controls."""
        self.results_tab.grid_rowconfigure(0, weight=1)
        self.results_tab.grid_rowconfigure(1, weight=0)
        self.results_tab.grid_columnconfigure(0, weight=1)

        # Canvas area
        self.results_canvas_frame = ttk.Frame(self.results_tab, style="TFrame")
        self.results_canvas_frame.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)
        self.results_canvas_frame.grid_rowconfigure(0, weight=1)
        self.results_canvas_frame.grid_columnconfigure(0, weight=1)

        self._init_results_canvas()

        # W/L controls at bottom
        wl_frame = ttk.LabelFrame(self.results_tab, text=i18n.PBI_PREVIEW_WL_TITLE)
        wl_frame.grid(row=1, column=0, sticky="ew", padx=4, pady=(0, 4))

        # Phase W/L
        tk.Label(wl_frame, text=f"{i18n.PBI_PLOT_PHASE}:", bg="#333333",
                 fg="white", font=("Segoe UI", 9)).grid(
            row=0, column=0, padx=(6, 2), pady=2, sticky="w")

        self.phase_wl_min = tk.DoubleVar(value=0)
        self.phase_wl_max = tk.DoubleVar(value=1)

        tk.Label(wl_frame, text="min:", bg="#333333", fg="#aaa",
                 font=("Segoe UI", 8)).grid(row=0, column=1, padx=2)
        tk.Entry(wl_frame, textvariable=self.phase_wl_min, width=10,
                 bg="#1e1e1e", fg="white", font=("Consolas", 9)).grid(
            row=0, column=2, padx=2)
        tk.Label(wl_frame, text="max:", bg="#333333", fg="#aaa",
                 font=("Segoe UI", 8)).grid(row=0, column=3, padx=2)
        tk.Entry(wl_frame, textvariable=self.phase_wl_max, width=10,
                 bg="#1e1e1e", fg="white", font=("Consolas", 9)).grid(
            row=0, column=4, padx=2)

        # Absorption W/L
        tk.Label(wl_frame, text=f"{i18n.PBI_PLOT_ABS}:", bg="#333333",
                 fg="white", font=("Segoe UI", 9)).grid(
            row=1, column=0, padx=(6, 2), pady=2, sticky="w")

        self.abs_wl_min = tk.DoubleVar(value=0)
        self.abs_wl_max = tk.DoubleVar(value=1)

        tk.Label(wl_frame, text="min:", bg="#333333", fg="#aaa",
                 font=("Segoe UI", 8)).grid(row=1, column=1, padx=2)
        tk.Entry(wl_frame, textvariable=self.abs_wl_min, width=10,
                 bg="#1e1e1e", fg="white", font=("Consolas", 9)).grid(
            row=1, column=2, padx=2)
        tk.Label(wl_frame, text="max:", bg="#333333", fg="#aaa",
                 font=("Segoe UI", 8)).grid(row=1, column=3, padx=2)
        tk.Entry(wl_frame, textvariable=self.abs_wl_max, width=10,
                 bg="#1e1e1e", fg="white", font=("Consolas", 9)).grid(
            row=1, column=4, padx=2)

        ttk.Button(wl_frame, text=i18n.PBI_PREVIEW_AUTO_WL,
                   command=self._auto_wl).grid(
            row=0, column=5, padx=(8, 4), pady=2)
        ttk.Button(wl_frame, text="Apply",
                   command=self._apply_wl).grid(
            row=1, column=5, padx=(8, 4), pady=2)

        # Pixel info readout
        self.result_info_var = tk.StringVar(value="")
        tk.Label(wl_frame, textvariable=self.result_info_var, bg="#333333",
                 fg="#cccccc", font=("Consolas", 8), anchor="w").grid(
            row=2, column=0, columnspan=6, sticky="ew", padx=6, pady=(0, 2))

    # ==================================================================
    #  Canvas initialisation
    # ==================================================================
    def _init_preview_canvas(self):
        """Create the initial empty preview matplotlib canvas."""
        rcParams.update(_MPL_DARK)

        self._preview_fig = Figure(figsize=(5, 5))
        self._preview_fig.set_facecolor("#333333")
        self._preview_ax = self._preview_fig.add_subplot(111)
        self._preview_ax.set_facecolor("#2a2a2a")
        self._preview_ax.set_title(i18n.PBI_PREVIEW_NONE, color="white",
                                   fontsize=10)
        self._preview_ax.tick_params(colors="white")

        self._preview_canvas = FigureCanvasTkAgg(self._preview_fig,
                                                 master=self.preview_canvas_frame)
        self._preview_canvas.draw()
        self._preview_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # Toolbar for zoom / pan / save
        tb_frame = tk.Frame(self.preview_canvas_frame, bg="#333333")
        tb_frame.grid(row=1, column=0, sticky="ew")
        self._preview_toolbar = NavigationToolbar2Tk(self._preview_canvas,
                                                     tb_frame)
        self._preview_toolbar.update()

        # Mouse-motion for pixel readout
        self._preview_canvas.mpl_connect("motion_notify_event",
                                         self._on_preview_mouse_move)

    # ------------------------------------------------------------------
    def _init_results_canvas(self):
        """Create the initial empty results matplotlib canvas (side-by-side)."""
        rcParams.update(_MPL_DARK)

        self._results_fig = Figure(figsize=(9, 4))
        self._results_fig.set_facecolor("#333333")

        self._phase_ax = self._results_fig.add_subplot(1, 2, 1)
        self._abs_ax = self._results_fig.add_subplot(1, 2, 2)
        for ax in (self._phase_ax, self._abs_ax):
            ax.set_facecolor("#2a2a2a")
            ax.tick_params(colors="white")

        self._phase_ax.set_title(i18n.PBI_PLOT_PHASE, color="white", fontsize=10)
        self._abs_ax.set_title(i18n.PBI_PLOT_ABS, color="white", fontsize=10)
        self._results_fig.tight_layout()

        self._results_canvas = FigureCanvasTkAgg(self._results_fig,
                                                 master=self.results_canvas_frame)
        self._results_canvas.draw()
        self._results_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # Toolbar
        tb_frame = tk.Frame(self.results_canvas_frame, bg="#333333")
        tb_frame.grid(row=1, column=0, sticky="ew")
        self._results_toolbar = NavigationToolbar2Tk(self._results_canvas,
                                                     tb_frame)
        self._results_toolbar.update()

        # Mouse-motion for pixel readout
        self._results_canvas.mpl_connect("motion_notify_event",
                                         self._on_result_mouse_move)

    # ==================================================================
    #  File loading callbacks
    # ==================================================================
    def _on_images_loaded(self, paths):
        """Load one or more TIFF files as intensity images."""
        self.intensity_images = []
        self.intensity_paths = []
        for p in sorted(paths):
            try:
                data = tifffile.imread(p).astype(np.float64)
                if data.ndim == 3:
                    for i in range(data.shape[0]):
                        self.intensity_images.append(data[i])
                        self.intensity_paths.append(
                            f"{os.path.basename(p)} [{i}]")
                else:
                    self.intensity_images.append(data)
                    self.intensity_paths.append(os.path.basename(p))
            except Exception as e:
                messagebox.showerror(i18n.ERR_TITLE,
                                     f"{i18n.PBI_ERR_LOAD_IMG}\n{e}")
                return

        n = len(self.intensity_images)
        self.img_dropzone.set_status_text(f"{n} " + i18n.PBI_IMAGES_LOADED)
        self.set_status(f"{n} {i18n.PBI_IMAGES_LOADED}")

        # Populate the list & auto-preview first image
        self._populate_image_list()
        if n > 0:
            self.img_listbox.selection_clear(0, tk.END)
            self.img_listbox.selection_set(0)
            self._show_preview(0)

    def _on_flat_loaded(self, paths):
        try:
            self.flat_image = tifffile.imread(paths[0]).astype(np.float64)
            if self.flat_image.ndim == 3:
                self.flat_image = self.flat_image.mean(axis=0)
            self.flat_dropzone.set_status_text(os.path.basename(paths[0]))
            self._populate_image_list()
        except Exception as e:
            messagebox.showerror(i18n.ERR_TITLE, str(e))

    def _on_dark_loaded(self, paths):
        try:
            self.dark_image = tifffile.imread(paths[0]).astype(np.float64)
            if self.dark_image.ndim == 3:
                self.dark_image = self.dark_image.mean(axis=0)
            self.dark_dropzone.set_status_text(os.path.basename(paths[0]))
            self._populate_image_list()
        except Exception as e:
            messagebox.showerror(i18n.ERR_TITLE, str(e))

    # ==================================================================
    #  Image list & preview
    # ==================================================================
    def _populate_image_list(self):
        """Refresh the listbox with loaded image names."""
        self.img_listbox.delete(0, tk.END)
        for idx, name in enumerate(self.intensity_paths):
            self.img_listbox.insert(tk.END, f"#{idx + 1}  {name}")
        if self.flat_image is not None:
            self.img_listbox.insert(tk.END,
                                    f"── {i18n.PBI_PREVIEW_FLAT} ──")
        if self.dark_image is not None:
            self.img_listbox.insert(tk.END,
                                    f"── {i18n.PBI_PREVIEW_DARK} ──")

    def _on_listbox_select(self, event):
        """Handle click in the image listbox."""
        sel = self.img_listbox.curselection()
        if not sel:
            return
        idx = sel[0]
        n_int = len(self.intensity_images)

        if idx < n_int:
            self._show_preview(idx)
        elif self.flat_image is not None and idx == n_int:
            self._show_preview_data(self.flat_image, i18n.PBI_PREVIEW_FLAT)
        else:
            if self.dark_image is not None:
                self._show_preview_data(self.dark_image, i18n.PBI_PREVIEW_DARK)

    def _show_preview(self, index: int):
        """Preview the index-th intensity image."""
        if index < 0 or index >= len(self.intensity_images):
            return
        data = self.intensity_images[index]
        title = i18n.PBI_PREVIEW_IMG_N.format(index + 1)
        self._show_preview_data(data, title)

    def _show_preview_data(self, data: np.ndarray, title: str):
        """Display a 2-D array on the preview canvas with auto-contrast."""
        rcParams.update(_MPL_DARK)
        self._preview_img_data = data

        ax = self._preview_ax

        # Clear everything and rebuild the axes to avoid colorbar accumulation
        self._preview_fig.clear()
        self._preview_ax = self._preview_fig.add_subplot(111)
        ax = self._preview_ax
        ax.set_facecolor("#2a2a2a")
        ax.tick_params(colors="white")

        vmin = float(np.nanpercentile(data, 1))
        vmax = float(np.nanpercentile(data, 99))
        if vmin == vmax:
            vmax = vmin + 1.0

        self._preview_im_handle = ax.imshow(data, cmap="gray",
                                            vmin=vmin, vmax=vmax,
                                            aspect="equal")
        ax.set_title(title, color="white", fontsize=10)
        self._preview_fig.colorbar(self._preview_im_handle, ax=ax,
                                   fraction=0.046, pad=0.04)
        self._preview_fig.tight_layout()
        self._preview_canvas.draw_idle()

        # Update info text
        ny, nx = data.shape
        info = (i18n.PBI_PREVIEW_SHAPE.format(ny, nx) + "   " +
                i18n.PBI_PREVIEW_RANGE.format(float(np.nanmin(data)),
                                               float(np.nanmax(data))))
        self.preview_info_var.set(info)

        # Switch to preview tab
        self.right_notebook.select(self.preview_tab)

    def _on_preview_mouse_move(self, event):
        """Show pixel (x, y) and value on mouse hover over preview."""
        if event.inaxes is None or self._preview_img_data is None:
            return
        x = int(round(event.xdata))
        y = int(round(event.ydata))
        data = self._preview_img_data
        ny, nx = data.shape
        if 0 <= x < nx and 0 <= y < ny:
            val = data[y, x]
            self.preview_info_var.set(
                i18n.PBI_PREVIEW_PIXEL_INFO.format(x, y, val) + "   " +
                i18n.PBI_PREVIEW_SHAPE.format(ny, nx))

    def _on_result_mouse_move(self, event):
        """Show pixel value on mouse move over result images."""
        if event.inaxes is None:
            return
        x = int(round(event.xdata))
        y = int(round(event.ydata))

        parts = []
        if event.inaxes == self._phase_ax and self.phase_result is not None:
            ny, nx = self.phase_result.shape
            if 0 <= x < nx and 0 <= y < ny:
                parts.append(
                    f"{i18n.PBI_PLOT_PHASE} "
                    f"{i18n.PBI_PREVIEW_PIXEL_INFO.format(x, y, self.phase_result[y, x])}")
        elif event.inaxes == self._abs_ax and self.absorption_result is not None:
            ny, nx = self.absorption_result.shape
            if 0 <= x < nx and 0 <= y < ny:
                parts.append(
                    f"{i18n.PBI_PLOT_ABS} "
                    f"{i18n.PBI_PREVIEW_PIXEL_INFO.format(x, y, self.absorption_result[y, x])}")

        if parts:
            self.result_info_var.set("  |  ".join(parts))

    # ==================================================================
    #  Distance parsing
    # ==================================================================
    def _parse_distances(self):
        """Parse DOD values (cm) from the text widget."""
        raw = self.dist_text.get("1.0", "end").strip()
        if not raw:
            return []
        tokens = raw.replace(",", " ").split()
        distances = []
        for t in tokens:
            try:
                distances.append(float(t))
            except ValueError:
                pass
        return distances

    # ==================================================================
    #  Run retrieval
    # ==================================================================
    def _run_retrieval(self):
        algo = self.algo_var.get()

        if not self.intensity_images:
            messagebox.showwarning(i18n.VAL_INVALID_PARAMS,
                                   i18n.PBI_ERR_NO_IMAGES)
            return

        distances = self._parse_distances()
        energy = self.energy_var.get()
        pixel = self.pixel_var.get()
        dso = self.dso_var.get() if self.beam_var.get() == "Conical" else 0.0
        pad = self.pad_var.get()

        if energy <= 0 or pixel <= 0:
            messagebox.showwarning(i18n.VAL_INVALID_PARAMS, i18n.PBI_ERR_PARAMS)
            return

        # Normalise
        imgs_norm = [normalize_image(im, self.flat_image, self.dark_image)
                     for im in self.intensity_images]

        def _work():
            nonlocal imgs_norm, distances
            self.set_status(i18n.PBI_ST_RUNNING.format(algo))
            try:
                if algo == "Paganin 单距离":
                    if len(distances) < 1:
                        raise ValueError(i18n.PBI_ERR_NEED_1_DIST)
                    thickness, phase = paganin_single_distance(
                        imgs_norm[0],
                        energy_keV=energy,
                        pixel_size_um=pixel,
                        dod_cm=distances[0],
                        delta_beta_ratio=self.delta_beta_var.get(),
                        dso_cm=dso,
                        pad=pad,
                    )
                    self.phase_result = phase
                    self.absorption_result = thickness

                elif algo == "TIE 多距离":
                    if (len(distances) < 2 or
                            len(distances) != len(imgs_norm)):
                        raise ValueError(
                            i18n.PBI_ERR_DIST_MISMATCH.format(
                                len(imgs_norm), len(distances)))
                    phase, absorption, _ = tie_multi_distance(
                        imgs_norm, distances,
                        energy_keV=energy,
                        pixel_size_um=pixel,
                        dso_cm=dso,
                        alpha=self.alpha_var.get(),
                        pad=pad,
                    )
                    self.phase_result = phase
                    self.absorption_result = absorption

                elif algo == "CTF 多距离":
                    if (len(distances) < 2 or
                            len(distances) != len(imgs_norm)):
                        raise ValueError(
                            i18n.PBI_ERR_DIST_MISMATCH.format(
                                len(imgs_norm), len(distances)))
                    phase, absorption = ctf_multi_distance(
                        imgs_norm, distances,
                        energy_keV=energy,
                        pixel_size_um=pixel,
                        dso_cm=dso,
                        alpha_phase=self.alpha_var.get(),
                        alpha_abs=self.alpha_var.get(),
                        pad=pad,
                    )
                    self.phase_result = phase
                    self.absorption_result = absorption
                else:
                    raise ValueError(f"Unknown algorithm: {algo}")

                self.master.after(0, self._show_results)
                self.master.after(0,
                                 lambda: self.set_status(i18n.PBI_ST_FINISHED))
                if self.zip_var.get():
                    self.master.after(0, self._export_zip)

            except Exception as e:
                tb = traceback.format_exc()
                self.master.after(
                    0, lambda: messagebox.showerror(
                        i18n.ERR_TITLE, f"{e}\n\n{tb}"))
                self.master.after(0,
                                 lambda: self.set_status(i18n.ST_ERROR))

        threading.Thread(target=_work, daemon=True).start()

    # ==================================================================
    #  Results display
    # ==================================================================
    def _show_results(self):
        """Refresh the results canvas with phase + absorption maps."""
        rcParams.update(_MPL_DARK)

        # Fully rebuild the figure to avoid colorbar accumulation
        self._results_fig.clear()
        self._phase_ax = self._results_fig.add_subplot(1, 2, 1)
        self._abs_ax = self._results_fig.add_subplot(1, 2, 2)

        for ax in (self._phase_ax, self._abs_ax):
            ax.set_facecolor("#2a2a2a")
            ax.tick_params(colors="white")

        self._result_handles = {}

        if self.phase_result is not None:
            p = self.phase_result
            vmin_p = float(np.nanpercentile(p, 1))
            vmax_p = float(np.nanpercentile(p, 99))
            if vmin_p == vmax_p:
                vmax_p = vmin_p + 1.0
            im1 = self._phase_ax.imshow(p, cmap="gray",
                                        vmin=vmin_p, vmax=vmax_p,
                                        aspect="equal")
            self._results_fig.colorbar(im1, ax=self._phase_ax,
                                       fraction=0.046, pad=0.04)
            self._result_handles["phase"] = im1
            self.phase_wl_min.set(round(vmin_p, 6))
            self.phase_wl_max.set(round(vmax_p, 6))
        self._phase_ax.set_title(i18n.PBI_PLOT_PHASE, color="white",
                                 fontsize=10)

        if self.absorption_result is not None:
            a = self.absorption_result
            vmin_a = float(np.nanpercentile(a, 1))
            vmax_a = float(np.nanpercentile(a, 99))
            if vmin_a == vmax_a:
                vmax_a = vmin_a + 1.0
            im2 = self._abs_ax.imshow(a, cmap="gray",
                                      vmin=vmin_a, vmax=vmax_a,
                                      aspect="equal")
            self._results_fig.colorbar(im2, ax=self._abs_ax,
                                       fraction=0.046, pad=0.04)
            self._result_handles["abs"] = im2
            self.abs_wl_min.set(round(vmin_a, 6))
            self.abs_wl_max.set(round(vmax_a, 6))
        self._abs_ax.set_title(i18n.PBI_PLOT_ABS, color="white", fontsize=10)

        self._results_fig.tight_layout()
        self._results_canvas.draw_idle()

        # Switch to results tab
        self.right_notebook.select(self.results_tab)

    # ------------------------------------------------------------------
    def _auto_wl(self):
        """Reset W/L entries to auto 1–99 percentile range."""
        if self.phase_result is not None:
            self.phase_wl_min.set(
                round(float(np.nanpercentile(self.phase_result, 1)), 6))
            self.phase_wl_max.set(
                round(float(np.nanpercentile(self.phase_result, 99)), 6))
        if self.absorption_result is not None:
            self.abs_wl_min.set(
                round(float(np.nanpercentile(self.absorption_result, 1)), 6))
            self.abs_wl_max.set(
                round(float(np.nanpercentile(self.absorption_result, 99)), 6))
        self._apply_wl()

    def _apply_wl(self):
        """Apply current W/L entries to result images (no re-computation)."""
        try:
            if "phase" in self._result_handles:
                self._result_handles["phase"].set_clim(
                    self.phase_wl_min.get(), self.phase_wl_max.get())
            if "abs" in self._result_handles:
                self._result_handles["abs"].set_clim(
                    self.abs_wl_min.get(), self.abs_wl_max.get())
            self._results_canvas.draw_idle()
        except Exception:
            pass

    # ==================================================================
    #  Save / Export
    # ==================================================================
    def _save_image(self, data, label="image"):
        if data is None:
            messagebox.showinfo(i18n.PBI_NO_RESULT_TITLE,
                                i18n.PBI_NO_RESULT_MSG)
            return
        path = asksaveasfilename(
            defaultextension=".tif",
            filetypes=[("TIFF", "*.tif"), ("All", "*.*")])
        if path:
            tifffile.imwrite(path, data.astype(np.float32))
            self.set_status(f"{label} → {path}")

    def _export_zip(self):
        path = asksaveasfilename(
            defaultextension=".zip",
            filetypes=[("ZIP", "*.zip"), ("All", "*.*")])
        if not path:
            return

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        config = {
            "algorithm": self.algo_var.get(),
            "energy_keV": self.energy_var.get(),
            "pixel_um": self.pixel_var.get(),
            "DSO_cm": self.dso_var.get(),
            "beam": self.beam_var.get(),
            "delta_beta": self.delta_beta_var.get(),
            "alpha": self.alpha_var.get(),
            "pad": self.pad_var.get(),
            "distances_cm": self._parse_distances(),
            "n_images": len(self.intensity_images),
            "timestamp": timestamp,
        }

        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("config.json", json.dumps(config, indent=2))
            if self.phase_result is not None:
                buf = io.BytesIO()
                tifffile.imwrite(buf, self.phase_result.astype(np.float32))
                zf.writestr("phase.tif", buf.getvalue())
            if self.absorption_result is not None:
                buf = io.BytesIO()
                tifffile.imwrite(buf,
                                 self.absorption_result.astype(np.float32))
                zf.writestr("absorption.tif", buf.getvalue())
            zf.writestr("README.txt",
                         f"XPCIpy PBI Phase Retrieval\n"
                         f"Timestamp: {timestamp}\n")

        self.set_status(i18n.PBI_ST_ZIP_SAVED)

    # ==================================================================
    #  Helpers
    # ==================================================================
    def set_status(self, text):
        if self.status_var is not None:
            self.status_var.set(text)
            try:
                self.master.update_idletasks()
            except Exception:
                pass
