# napari_dea/_widget.py
"""
Yasas Wijesekara, 2025
University Medicine Greifswald

Stratified differential expression widget for napari-spatialdata.
Supports per-group DE (e.g. celltype, condition) via Wilcoxon or DESeq2,
with violin-plot popups and background threading.
"""

from typing import Optional, List
import numpy as np
import pandas as pd
from napari.utils.notifications import show_info, show_error
from qtpy.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QComboBox,
    QPushButton,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QLineEdit,
    QProgressBar,
    QSpinBox,
    QFrame,
)
from qtpy.QtCore import QThread

from spatialdata._types import ArrayLike
from geopandas import GeoDataFrame
from spatialdata.transformations import Identity
from spatialdata.models import force_2d
from napari.viewer import Viewer
from napari.layers import Shapes, Points
import scipy.sparse as sp
from spatialdata import polygon_query
from shapely.ops import unary_union
from shapely.geometry import Polygon
from shapely import make_valid

from napari_dea._worker import DEAWorker
from napari_dea._plot_dialog import ViolinPlotDialog


class TwoRegionDEWidget(QWidget):
    """
    Dock widget for differential expression analysis between two spatial regions.

    Assumptions:
    - A SpatialData object is available as ``viewer.sdata``
      or as ``layer.metadata["sdata"]`` on some layer (see ``_get_sdata()``).
    - For the selected table, there is an AnnData in ``sdata.tables[table_name]``.
    """

    def __init__(self, viewer: Viewer):
        super().__init__()
        self.viewer = viewer

        self._last_de: Optional[pd.DataFrame] = None
        self._last_adata = None
        self._thread: Optional[QThread] = None
        self._worker: Optional[DEAWorker] = None

        self._build_ui()
        self._populate_tables()
        self._populate_shape_layers()
        self._update_groupby_options()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _build_ui(self):
        from qtpy.QtWidgets import QGroupBox

        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        # ═════════════════════════════════════════════════════════════════
        # Data Selection
        # ═════════════════════════════════════════════════════════════════
        group_data = QGroupBox("Data Selection")
        layout_data = QVBoxLayout()
        group_data.setLayout(layout_data)

        row_data = QHBoxLayout()

        col_tbl = QVBoxLayout()
        col_tbl.addWidget(QLabel("Table (AnnData in SpatialData):"))
        self.table_combo = QComboBox()
        self.table_combo.currentTextChanged.connect(self._update_groupby_options)
        col_tbl.addWidget(self.table_combo)
        row_data.addLayout(col_tbl, stretch=2)

        col_a = QVBoxLayout()
        col_a.addWidget(QLabel("Region A (Shapes):"))
        self.regionA_combo = QComboBox()
        col_a.addWidget(self.regionA_combo)
        row_data.addLayout(col_a, stretch=1)

        col_b = QVBoxLayout()
        col_b.addWidget(QLabel("Region B (Shapes):"))
        self.regionB_combo = QComboBox()
        col_b.addWidget(self.regionB_combo)
        row_data.addLayout(col_b, stretch=1)

        layout_data.addLayout(row_data)
        main_layout.addWidget(group_data)

        # ═════════════════════════════════════════════════════════════════
        # DE Parameters
        # ═════════════════════════════════════════════════════════════════
        group_params = QGroupBox("DE Parameters")
        layout_params = QVBoxLayout()
        group_params.setLayout(layout_params)

        row_params = QHBoxLayout()

        col_g = QVBoxLayout()
        col_g.addWidget(QLabel("Group DE by (obs column):"))
        self.groupby_combo = QComboBox()
        col_g.addWidget(self.groupby_combo)
        row_params.addLayout(col_g, stretch=2)

        col_m = QVBoxLayout()
        col_m.addWidget(QLabel("DE Method:"))
        self.method_combo = QComboBox()
        self.method_combo.addItems(["Wilcoxon (scanpy)", "DESeq2 (pydeseq2)"])
        col_m.addWidget(self.method_combo)
        row_params.addLayout(col_m, stretch=1)

        col_min = QVBoxLayout()
        col_min.addWidget(QLabel("Min cells per region:"))
        self.min_cells_spin = QSpinBox()
        self.min_cells_spin.setRange(1, 10000)
        self.min_cells_spin.setValue(5)
        col_min.addWidget(self.min_cells_spin)
        row_params.addLayout(col_min, stretch=1)

        layout_params.addLayout(row_params)
        main_layout.addWidget(group_params)

        # ═════════════════════════════════════════════════════════════════
        # Execution
        # ═════════════════════════════════════════════════════════════════
        group_exec = QGroupBox("Execution")
        layout_exec = QHBoxLayout()
        group_exec.setLayout(layout_exec)

        run_frame = QFrame()
        run_frame.setFrameShape(QFrame.StyledPanel)
        run_frame_layout = QHBoxLayout()
        run_frame_layout.setContentsMargins(4, 4, 4, 4)
        run_frame.setLayout(run_frame_layout)

        self.run_button = QPushButton("Run DE: Region A vs Region B")
        self.run_button.setMinimumHeight(32)
        self.run_button.clicked.connect(self._on_run_clicked)
        run_frame_layout.addWidget(self.run_button)
        layout_exec.addWidget(run_frame)

        prog_frame = QFrame()
        prog_frame.setFrameShape(QFrame.StyledPanel)
        prog_frame_layout = QHBoxLayout()
        prog_frame_layout.setContentsMargins(4, 4, 4, 4)
        prog_frame.setLayout(prog_frame_layout)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        prog_frame_layout.addWidget(self.progress_bar)
        layout_exec.addWidget(prog_frame, stretch=1)

        main_layout.addWidget(group_exec)

        # ═════════════════════════════════════════════════════════════════
        # Visualization & Results
        # ═════════════════════════════════════════════════════════════════
        group_vis = QGroupBox("Visualization & Results")
        layout_vis = QVBoxLayout()
        group_vis.setLayout(layout_vis)

        # Top row: Color region + Filter + Gene
        row_vis_top = QHBoxLayout()

        col_cr = QVBoxLayout()
        col_cr.addWidget(QLabel("Color region:"))
        self.color_region_combo = QComboBox()
        self.color_region_combo.addItems(
            ["Whole slide", "Only A", "Only B", "Union (A ∪ B)"]
        )
        col_cr.addWidget(self.color_region_combo)
        row_vis_top.addLayout(col_cr)

        col_filt = QVBoxLayout()
        col_filt.addWidget(QLabel("Filter results by group:"))
        self.group_filter_combo = QComboBox()
        self.group_filter_combo.addItem("All")
        self.group_filter_combo.currentTextChanged.connect(self._on_filter_changed)
        col_filt.addWidget(self.group_filter_combo)
        row_vis_top.addLayout(col_filt)

        col_gene = QVBoxLayout()
        col_gene.addWidget(QLabel("Gene to visualize:"))
        gene_frame = QFrame()
        gene_frame.setFrameShape(QFrame.StyledPanel)
        gene_frame_layout = QHBoxLayout()
        gene_frame_layout.setContentsMargins(4, 4, 4, 4)
        gene_frame.setLayout(gene_frame_layout)
        self.gene_line = QLineEdit()
        self.gene_line.setPlaceholderText("e.g. ACTB")
        gene_frame_layout.addWidget(self.gene_line)
        col_gene.addWidget(gene_frame)
        row_vis_top.addLayout(col_gene, stretch=2)

        layout_vis.addLayout(row_vis_top)

        # Color button row
        btn_row = QHBoxLayout()
        btn_row.addStretch(1)
        self.show_gene_button = QPushButton("Color by Gene")
        self.show_gene_button.clicked.connect(self._on_show_gene_clicked)
        btn_row.addWidget(self.show_gene_button)
        layout_vis.addLayout(btn_row)

        layout_vis.addWidget(QLabel("DE results (Region A vs Region B):"))
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(5)
        self.results_table.setHorizontalHeaderLabels(
            ["group", "gene", "logFC", "p_adj", "score"]
        )
        self.results_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch
        )
        self.results_table.setAlternatingRowColors(False)
        self.results_table.setStyleSheet("""
            QTableWidget {
                background-color: #2b2b2b;
                color: #ffffff;
                gridline-color: #3d3d3d;
            }
            QTableWidget::item {
                background-color: #2b2b2b;
                color: #ffffff;
                padding: 4px;
            }
            QTableWidget::item:selected {
                background-color: #4a4a4a;
                color: #ffffff;
            }
            QHeaderView::section {
                background-color: #3d3d3d;
                color: #ffffff;
                padding: 4px;
                border: 1px solid #555555;
            }
        """)
        layout_vis.addWidget(self.results_table)
        self.results_table.cellClicked.connect(self._on_result_click)
        self.results_table.cellDoubleClicked.connect(self._on_result_double_click)

        main_layout.addWidget(group_vis, stretch=1)
        main_layout.addStretch(0)

    # ------------------------------------------------------------------
    # SpatialData + layers discovery
    # ------------------------------------------------------------------
    def _get_sdata(self):
        """
        Retrieve the SpatialData object.
        """
        if hasattr(self.viewer, "sdata"):
            return self.viewer.sdata

        for layer in self.viewer.layers:
            if isinstance(layer.metadata, dict) and "sdata" in layer.metadata:
                return layer.metadata["sdata"]

        return None

    def _populate_tables(self):
        sdata = self._get_sdata()
        if sdata is None:
            return
        self.table_combo.clear()
        self._table_names = list(getattr(sdata, "tables", {}).keys())
        self.table_combo.addItems(self._table_names)

    def _populate_shape_layers(self):
        """Fill combos with all Shapes layers and Points layers."""
        self.regionA_combo.clear()
        self.regionB_combo.clear()

        shapes_names: List[str] = []
        points_names: List[str] = []
        for layer in self.viewer.layers:
            if isinstance(layer, Shapes):
                shapes_names.append(layer.name)
            if isinstance(layer, Points):
                points_names.append(layer.name)

        self.regionA_combo.addItems(shapes_names)
        self.regionB_combo.addItems(shapes_names)

    def _update_groupby_options(self):
        """Populate the group-by dropdown from the selected table's obs columns."""
        self.groupby_combo.clear()
        sdata = self._get_sdata()
        if sdata is None:
            return
        table_name = self.table_combo.currentText()
        if not table_name or table_name not in getattr(sdata, "tables", {}):
            return
        adata = sdata.tables[table_name]
        candidates = []
        for col in adata.obs.columns:
            dtype = adata.obs[col].dtype
            is_cat = hasattr(dtype, "categories")
            is_obj = pd.api.types.is_object_dtype(dtype)
            is_bool = pd.api.types.is_bool_dtype(dtype)
            nuniq = adata.obs[col].nunique()
            if is_cat or is_obj or is_bool or nuniq <= 200:
                candidates.append(col)
        self.groupby_combo.addItems(candidates)

    # ------------------------------------------------------------------
    # Helpers for polygon_query
    # ------------------------------------------------------------------
    def _indices_from_polygon_query(
        self,
        layer: Shapes,
        table_name: str,
    ) -> np.ndarray:
        """
        Use spatialdata.polygon_query to get indices of cells in the table that
        fall inside the union of napari-drawn polygons.
        """
        if "sdata" not in layer.metadata:
            raise ValueError("Shapes layer has no 'sdata' in metadata.")
        if "_current_cs" not in layer.metadata:
            raise ValueError("Shapes layer has no '_current_cs' in metadata.")

        sdata = layer.metadata["sdata"]
        coordinate_system = layer.metadata["_current_cs"]

        if len(layer.data) == 0:
            raise ValueError("Cannot query ROI: shapes layer has no shapes.")

        # Convert napari shapes to world coords, then to shapely Polygons
        coords = [
            np.array([layer.data_to_world(xy) for xy in shape._data])
            for shape in layer._data_view.shapes
        ]

        def _fix_coords(coords: ArrayLike) -> ArrayLike:
            # drop Z if present and flip (row, col) -> (x, y) if needed
            remove_z = coords.shape[1] == 3
            first_index = 1 if remove_z else 0
            coords = coords[:, first_index::]
            return np.fliplr(coords)

        polygons: list[Polygon] = [Polygon(_fix_coords(p)) for p in coords]

        if not polygons:
            return np.array([], dtype=int)

        # Fix invalid geometries (self-intersections, etc.) before unioning
        polygons = [make_valid(p) if not p.is_valid else p for p in polygons]
        polygons = [p for p in polygons if not p.is_empty]

        if not polygons:
            return np.array([], dtype=int)

        gdf = GeoDataFrame({"geometry": polygons})
        force_2d(gdf)

        # Merge all shapes into a single Polygon/MultiPolygon
        try:
            poly_union = unary_union(polygons)
        except Exception as e:
            show_error(f"Failed to merge ROI polygons: {e}. Try redrawing the shapes.")
            return np.array([], dtype=int)

        try:
            roi_sdata = polygon_query(
                sdata,
                polygon=poly_union,
                target_coordinate_system=coordinate_system,
                filter_table=True,
                clip=False,
            )
        except Exception as e:
            show_error(f"polygon_query failed: {e}")
            return np.array([], dtype=int)

        # No data found
        if roi_sdata is None or table_name not in getattr(roi_sdata, "tables", {}):
            return np.array([], dtype=int)

        # Map ROI obs back to the full table via obs_names
        full_table = sdata.tables[table_name]
        roi_table = roi_sdata.tables[table_name]

        mask = full_table.obs_names.isin(roi_table.obs_names)
        idx = np.where(mask)[0]
        return idx

    # ------------------------------------------------------------------
    # reset state
    # ------------------------------------------------------------------
    def reset_state(self):
        """Reset widget UI and internal state."""
        # Clear stored data
        self._last_de = None
        self._last_adata = None
        if hasattr(self, "_full_adata"):
            self._full_adata = None
        if hasattr(self, "_keep_mask"):
            self._keep_mask = None
        if hasattr(self, "_keep_idx"):
            self._keep_idx = None
        if hasattr(self, "_last_table_name"):
            self._last_table_name = None

        # Clear UI elements
        self.gene_line.clear()
        self.results_table.setRowCount(0)
        self.group_filter_combo.clear()
        self.group_filter_combo.addItem("All")
        self.progress_bar.setValue(0)

        # Repopulate dropdowns from current viewer state
        self._populate_tables()
        self._populate_shape_layers()
        self._update_groupby_options()

    def showEvent(self, event):
        """Called when the dock widget is shown (toggled on)."""
        super().showEvent(event)

        # Napari toggles visibility instead of destroying the widget.
        # Each time it's shown again, reset its state.
        show_info("resetting DEWidget state")
        self.reset_state()

    # ------------------------------------------------------------------
    # DE pipeline
    # ------------------------------------------------------------------
    def _get_layer(self, layer_name: str) -> Shapes:
        layer = self.viewer.layers[layer_name]
        if not isinstance(layer, Shapes):
            raise ValueError(f"Layer {layer_name} is not a Shapes layer.")
        return layer

    def _on_run_clicked(self):
        try:
            if self._thread is not None and self._thread.isRunning():
                show_info("DE analysis is already running.")
                return
        except RuntimeError:
            # C++ object has been deleted (previous thread finished)
            self._thread = None
            self._worker = None

        show_info("Please wait ...")
        sdata = self._get_sdata()
        if sdata is None:
            show_error("No SpatialData object found.")
            return

        table_name = self.table_combo.currentText()
        regionA_layer = self.regionA_combo.currentText()
        regionB_layer = self.regionB_combo.currentText()
        groupby_col = self.groupby_combo.currentText()

        if not table_name or not regionA_layer or not regionB_layer:
            show_error("Please select table and both region layers.")
            return
        if not groupby_col:
            show_error("Please select a group-by column.")
            return

        method_text = self.method_combo.currentText()
        method = "wilcoxon" if "Wilcoxon" in method_text else "deseq2"

        if method == "deseq2":
            try:
                import pydeseq2  # noqa: F401
            except ImportError:
                show_error(
                    "pydeseq2 is not installed. Install it with: pip install pydeseq2"
                )
                return

        # --- get the selected layers ---
        roiA = self._get_layer(regionA_layer)
        roiB = self._get_layer(regionB_layer)

        # --- Use spatialdata.polygon_query to get indices for A and B ---
        idxA = self._indices_from_polygon_query(layer=roiA, table_name=table_name)
        idxB = self._indices_from_polygon_query(layer=roiB, table_name=table_name)

        if len(idxA) == 0 or len(idxB) == 0:
            show_error("One of the regions contains no cells/spots.")
            return

        # ------------------------------------------------------------------
        # Build region labels on the FULL AnnData
        # ------------------------------------------------------------------
        full_adata = sdata.tables[table_name].copy()
        n = full_adata.n_obs
        region = np.full(n, "other", dtype=object)
        region[idxA] = "regionA"
        region[idxB] = "regionB"
        full_adata.obs["__region__"] = region

        # Save full context for later visualization
        self._full_adata = full_adata
        self._last_table_name = table_name
        self._idxA = idxA
        self._idxB = idxB
        self._keep_mask = np.isin(region, ["regionA", "regionB"])
        self._keep_idx = np.where(self._keep_mask)[0]

        min_cells = self.min_cells_spin.value()

        self.run_button.setEnabled(False)
        self.progress_bar.setValue(0)

        # ------------------------------------------------------------------
        # Run DE in background thread
        # ------------------------------------------------------------------
        self._thread = QThread()
        self._worker = DEAWorker(
            full_adata, idxA, idxB, groupby_col, min_cells, method
        )
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.progress.connect(self._on_worker_progress)
        self._worker.finished.connect(self._on_worker_finished)
        self._worker.error.connect(self._on_worker_error)
        self._worker.finished.connect(self._thread.quit)
        self._worker.finished.connect(self._worker.deleteLater)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.start()

    def _on_worker_progress(self, current: int, total: int):
        if total > 0:
            self.progress_bar.setValue(int(100 * current / total))

    def _on_worker_finished(self, de_df: pd.DataFrame):
        self.run_button.setEnabled(True)
        self._last_de = de_df
        self._populate_results_table(de_df)
        show_info("DE analysis complete.")

    def _on_worker_error(self, msg: str):
        self.run_button.setEnabled(True)
        show_error(msg)

    def _populate_results_table(self, de_df: pd.DataFrame, n_top: int = 200):
        self.results_table.setRowCount(0)
        if de_df is None or len(de_df) == 0:
            return

        de_df = de_df.sort_values("p_adj").head(n_top)

        self.results_table.setRowCount(len(de_df))
        for i, (_, row) in enumerate(de_df.iterrows()):
            self.results_table.setItem(i, 0, QTableWidgetItem(str(row["group"])))
            self.results_table.setItem(i, 1, QTableWidgetItem(str(row["gene"])))
            self.results_table.setItem(
                i, 2, QTableWidgetItem(f"{row['logFC']:.3f}")
            )
            self.results_table.setItem(
                i, 3, QTableWidgetItem(f"{row['p_adj']:.2e}")
            )
            self.results_table.setItem(
                i, 4, QTableWidgetItem(f"{row['score']:.3f}")
            )

        # Update filter combo (block signals to avoid recursion)
        self.group_filter_combo.blockSignals(True)
        current_filter = self.group_filter_combo.currentText()
        self.group_filter_combo.clear()
        self.group_filter_combo.addItem("All")
        groups = sorted(de_df["group"].unique())
        self.group_filter_combo.addItems([str(g) for g in groups])
        if current_filter in [str(g) for g in groups] or current_filter == "All":
            self.group_filter_combo.setCurrentText(current_filter)
        else:
            self.group_filter_combo.setCurrentText("All")
        self.group_filter_combo.blockSignals(False)

    def _on_filter_changed(self):
        filter_text = self.group_filter_combo.currentText()
        if not filter_text or self._last_de is None:
            return
        if filter_text == "All":
            self._populate_results_table(self._last_de)
        else:
            filtered = self._last_de[self._last_de["group"] == filter_text]
            self._populate_results_table(filtered)

    def _on_result_click(self, row: int, column: int):
        """Single-click on a DE result row → populate gene line edit."""
        item_gene = self.results_table.item(row, 1)
        if item_gene is not None:
            self.gene_line.setText(item_gene.text())

    def _on_result_double_click(self, row: int, column: int):
        """Double-click on a DE result row → open violin plot popup."""
        item_group = self.results_table.item(row, 0)
        item_gene = self.results_table.item(row, 1)
        if item_group is None or item_gene is None:
            return
        group = item_group.text()
        gene = item_gene.text()

        if not hasattr(self, "_full_adata") or self._full_adata is None:
            show_error("Run DE first (or provide an AnnData context).")
            return

        groupby_col = self.groupby_combo.currentText()
        dialog = ViolinPlotDialog(
            self._full_adata, gene, group, groupby_col, parent=self
        )
        dialog.show()

    # ------------------------------------------------------------------
    # Visualization: color a points layer by gene expression
    # ------------------------------------------------------------------
    def _get_target_layer(self):
        """Auto-select the best Points/Shapes layer to color."""
        if hasattr(self, "_last_table_name") and self._last_table_name:
            for layer in self.viewer.layers:
                if isinstance(layer, (Points, Shapes)):
                    table_names = layer.metadata.get("table_names", [])
                    if self._last_table_name in table_names:
                        return layer
        for layer in self.viewer.layers:
            if isinstance(layer, (Points, Shapes)):
                return layer
        return None

    def _map_instance_ids_to_positions(self, pl, full_adata):
        """
        Map layer points to full_adata positions.
        Returns (positions, valid_mask) where:
          - positions[i] = table row index for layer point i (only meaningful where valid)
          - valid_mask[i] = True if point i has a corresponding table row
        Returns (None, None) if no mapping is possible.
        """
        n_points = len(pl.data)

        if "indices" not in pl.metadata:
            if n_points == full_adata.n_obs:
                return np.arange(n_points), np.ones(n_points, dtype=bool)
            return None, None

        instance_ids = pl.metadata["indices"]
        if len(instance_ids) != n_points:
            return None, None

        positions = np.full(n_points, -1, dtype=int)
        valid_mask = np.zeros(n_points, dtype=bool)

        # Strategy 1: SpatialData table keys (instance_key + region_key)
        try:
            from spatialdata.models import get_table_keys

            _, region_key, instance_key = get_table_keys(full_adata)
            element_name = pl.metadata.get("name")

            if (
                region_key in full_adata.obs.columns
                and instance_key in full_adata.obs.columns
            ):
                if element_name is not None:
                    mask = full_adata.obs[region_key] == element_name
                    element_positions = np.where(mask)[0]
                    element_ids = full_adata.obs.loc[mask, instance_key].values
                else:
                    element_positions = np.arange(full_adata.n_obs)
                    element_ids = full_adata.obs[instance_key].values

                id_to_pos = pd.Series(element_positions, index=element_ids)
                mapped = id_to_pos.reindex(instance_ids).values
                valid = ~np.isnan(mapped)
                if valid.any():
                    positions[valid] = mapped[valid].astype(int)
                    valid_mask[:] = valid
                    return positions, valid_mask
        except Exception:
            pass

        # Strategy 2: Match against obs_names
        obs_pos_map = pd.Series(range(full_adata.n_obs), index=full_adata.obs_names)
        mapped = obs_pos_map.reindex(instance_ids).values
        valid = ~np.isnan(mapped)
        if valid.any():
            positions[valid] = mapped[valid].astype(int)
            valid_mask[:] = valid
            return positions, valid_mask

        # Strategy 3: Treat as positional indices (last resort)
        mapped = np.array(instance_ids, dtype=int)
        valid = (mapped >= 0) & (mapped < full_adata.n_obs)
        if valid.any():
            positions[valid] = mapped[valid]
            valid_mask[:] = valid
            return positions, valid_mask

        return None, None

    def _on_show_gene_clicked(self):
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors

        gene = self.gene_line.text().strip()
        if not gene:
            return

        if not hasattr(self, "_full_adata") or self._full_adata is None:
            show_error("Run DE first (or provide an AnnData context).")
            return

        full_adata = self._full_adata

        if gene not in full_adata.var_names:
            show_error(f"Gene {gene} not in var_names.")
            return

        # Get full expression vector (one value per obs in the full table)
        x = full_adata[:, gene].X
        if sp.issparse(x):
            expr_full = x.toarray().ravel()
        else:
            expr_full = np.asarray(x).ravel()

        # Auto-detect target layer
        pl = self._get_target_layer()
        if pl is None:
            show_error("No Points or Shapes layer found to color.")
            return

        n_points = len(pl.data)

        # Map layer points to full_adata positions (may be partial)
        positions, valid_mask = self._map_instance_ids_to_positions(pl, full_adata)
        if positions is None:
            show_error(
                f"Cannot map layer '{pl.name}' to the table. "
                f"({n_points} points vs {full_adata.n_obs} table rows). "
                f"Coloring aborted."
            )
            return

        # Build per-point expression array (NaN = no mapping or no data)
        expr = np.full(n_points, np.nan, dtype=float)
        expr[valid_mask] = expr_full[positions[valid_mask]]

        # Build region mask on the full table
        region_mode = self.color_region_combo.currentText()
        n_full = len(expr_full)
        region_mask_full = np.zeros(n_full, dtype=bool)

        if region_mode == "Only A":
            if hasattr(self, "_idxA"):
                region_mask_full[self._idxA] = True
        elif region_mode == "Only B":
            if hasattr(self, "_idxB"):
                region_mask_full[self._idxB] = True
        elif region_mode == "Union (A ∪ B)":
            if hasattr(self, "_keep_mask"):
                region_mask_full = self._keep_mask.copy()
        else:  # Whole slide
            region_mask_full[:] = True

        # Map full-table mask to layer points
        layer_mask = np.zeros(n_points, dtype=bool)
        layer_mask[valid_mask] = region_mask_full[positions[valid_mask]]

        valid_expr = expr[valid_mask]
        if len(valid_expr) == 0 or not np.any(np.isfinite(valid_expr)):
            show_error("No valid expression values to color.")
            return

        # Start with gray for ALL points
        rgba = np.tile([0.5, 0.5, 0.5, 1.0], (n_points, 1))

        # Compute colormap only for valid points
        vmin = float(np.nanmin(valid_expr))
        vmax = float(np.nanmax(valid_expr))
        if vmin == vmax:
            vmax = vmin + 1e-6

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cmap = cm.get_cmap("viridis")
        valid_rgba = cmap(norm(valid_expr))

        # Within valid points, gray out those outside the selected region
        valid_selected = layer_mask[valid_mask]
        valid_rgba[~valid_selected] = [0.5, 0.5, 0.5, 1.0]

        rgba[valid_mask] = valid_rgba

        # Assign direct colors
        pl.face_color = rgba
        if hasattr(pl, "properties"):
            pl.properties["expr"] = expr

        n_mapped = int(valid_mask.sum())
        n_selected = int(layer_mask.sum())
        show_info(
            f"Colored '{pl.name}' by {gene} ({region_mode}). "
            f"Mapped {n_mapped}/{n_points} points; {n_selected} inside region."
        )


def two_region_de_widget(viewer: Viewer) -> TwoRegionDEWidget:
    """Factory function used by napari to create the dock widget."""
    return TwoRegionDEWidget(viewer)
