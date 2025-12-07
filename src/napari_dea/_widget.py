# napari_dea/_widget.py
"""
Yasas Wijesekara, 205
University Medicine Greifswald
"""

from typing import Optional, List
import numpy as np
import pandas as pd
from napari.utils.notifications import show_info, show_error
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QComboBox, QPushButton,
    QLabel, QTableWidget, QTableWidgetItem, QHeaderView, QLineEdit
)

from spatialdata._types import ArrayLike
from geopandas import GeoDataFrame
from spatialdata.transformations import Identity
from spatialdata.models import force_2d
from napari.viewer import Viewer
from napari.layers import Shapes, Points
import scanpy as sc
import scipy.sparse as sp
from spatialdata import polygon_query
from shapely.ops import unary_union
from shapely.geometry import Polygon


class TwoRegionDEWidget(QWidget):
    """
    Dock widget for differential expression between two spatial regions.

    Assumptions:
    - A SpatialData object is available as `viewer.sdata`
      or as `layer.metadata["sdata"]` on some layer (see _get_sdata()).
    - For the selected table, there is an AnnData in `sdata.tables[table_name]`.
    """

    def __init__(self, viewer: Viewer):
        super().__init__()
        self.viewer = viewer

        self._last_de: Optional[pd.DataFrame] = None
        self._last_adata = None

        self._build_ui()
        self._populate_tables()
        self._populate_shape_layers()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _build_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        # --- Table selection (AnnData) ---
        layout.addWidget(QLabel("Table (AnnData in SpatialData):"))
        self.table_combo = QComboBox()
        layout.addWidget(self.table_combo)

        # --- Region A shapes layer ---
        layout.addWidget(QLabel("Region A (Shapes layer):"))
        self.regionA_combo = QComboBox()
        layout.addWidget(self.regionA_combo)

        # --- Region B shapes layer ---
        layout.addWidget(QLabel("Region B (Shapes layer):"))
        self.regionB_combo = QComboBox()
        layout.addWidget(self.regionB_combo)

        # --- Optional: name of points layer or shapes layer to color ---
        layout.addWidget(QLabel("Points layer to color (optional):"))
        self.points_layer_combo = QComboBox()
        layout.addWidget(self.points_layer_combo)

        # --- Run button ---
        self.run_button = QPushButton("Run DE: Region A vs Region B")
        layout.addWidget(self.run_button)
        self.run_button.clicked.connect(self._on_run_clicked)

        # --- Gene selection and visualization ---
        hl = QHBoxLayout()
        hl.addWidget(QLabel("Gene to visualize:"))
        self.gene_line = QLineEdit()
        hl.addWidget(self.gene_line)
        self.show_gene_button = QPushButton("Color points/shapes by gene")
        hl.addWidget(self.show_gene_button)
        layout.addLayout(hl)
        self.show_gene_button.clicked.connect(self._on_show_gene_clicked)

        # --- Results table ---
        layout.addWidget(QLabel("DE results (Region A vs Region B):"))
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(4)
        self.results_table.setHorizontalHeaderLabels(
            ["gene", "logFC", "p_adj", "score"]
        )
        self.results_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch
        )
        layout.addWidget(self.results_table)

        # Click a row → populate gene line
        self.results_table.cellDoubleClicked.connect(self._on_result_double_click)

    # ------------------------------------------------------------------
    # SpatialData + layers discovery
    # ------------------------------------------------------------------
    def _get_sdata(self):
        """
        Retrieve the SpatialData object.

        You MUST adapt this to how you integrate SpatialData with napari.
        Examples:
        - viewer has attribute: self.viewer.sdata
        - or it's stored in a special layer's metadata
        """
        # EXAMPLE 1: direct attribute (your own app)
        if hasattr(self.viewer, "sdata"):
            return self.viewer.sdata

        # EXAMPLE 2: stored in metadata of a specific layer
        for layer in self.viewer.layers:
            if isinstance(layer.metadata, dict) and "sdata" in layer.metadata:
                return layer.metadata["sdata"]

        # If nothing found, return None
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
        self.points_layer_combo.clear()

        shapes_names: List[str] = []
        points_names: List[str] = []
        for layer in self.viewer.layers:
            if isinstance(layer, Shapes):
                shapes_names.append(layer.name)
            if isinstance(layer, Points):
                points_names.append(layer.name)

        self.regionA_combo.addItems(shapes_names)
        self.regionB_combo.addItems(shapes_names)
        self.points_layer_combo.addItem("")  # optional blank
        self.points_layer_combo.addItems(points_names + shapes_names)

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
        gdf = GeoDataFrame({"geometry": polygons})

        # enforce 2D geometry
        force_2d(gdf)

        if not polygons:
            return np.array([], dtype=int)

        # Merge all shapes into a single Polygon/MultiPolygon
        poly_union = unary_union(polygons)

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

        # Repopulate dropdowns from current viewer state
        self._populate_tables()
        self._populate_shape_layers()
        

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
        show_info("Please wait ...")
        sdata = self._get_sdata()
        if sdata is None:
            show_error("No SpatialData object found.")
            return

        table_name = self.table_combo.currentText()
        regionA_layer = self.regionA_combo.currentText()
        regionB_layer = self.regionB_combo.currentText()

        if not table_name or not regionA_layer or not regionB_layer:
            show_error("Please select table and both region layers.")
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
        self._keep_mask = np.isin(region, ["regionA", "regionB"])
        self._keep_idx = np.where(self._keep_mask)[0]

        # ------------------------------------------------------------------
        # Subset to regionA ∪ regionB for DE
        # ------------------------------------------------------------------
        adata = full_adata[self._keep_mask].copy()

        # Make __region__ categorical with fixed categories
        adata.obs["__region__"] = pd.Categorical(
            adata.obs["__region__"],
            categories=["regionA", "regionB"],
        )

        # ------------------------------------------------------------------
        # Run DE: regionA vs regionB using scanpy
        # ------------------------------------------------------------------
        sc.tl.rank_genes_groups(
            adata,
            groupby="__region__",
            groups=["regionA"],
            reference="regionB",
            method="wilcoxon",
        )
        de_df = sc.get.rank_genes_groups_df(adata, group="regionA")
        self._last_de = de_df
        self._last_adata = adata

        self._populate_results_table(de_df)

    def _populate_results_table(self, de_df: pd.DataFrame, n_top: int = 100):
        if de_df is None or len(de_df) == 0:
            return

        de_df = de_df.sort_values("pvals_adj").head(n_top)

        self.results_table.setRowCount(len(de_df))
        # cols = ["names", "logfoldchanges", "pvals_adj", "scores"]

        for i, (_, row) in enumerate(de_df.iterrows()):
            gene = str(row["names"])
            logfc = f"{row['logfoldchanges']:.3f}"
            padj = f"{row['pvals_adj']:.2e}"
            score = f"{row['scores']:.3f}"

            self.results_table.setItem(i, 0, QTableWidgetItem(gene))
            self.results_table.setItem(i, 1, QTableWidgetItem(logfc))
            self.results_table.setItem(i, 2, QTableWidgetItem(padj))
            self.results_table.setItem(i, 3, QTableWidgetItem(score))

    def _on_result_double_click(self, row: int, column: int):
        """Double-click on a DE result row → put gene name into line edit."""
        item = self.results_table.item(row, 0)
        if item is not None:
            self.gene_line.setText(item.text())

    # ------------------------------------------------------------------
    # Visualization: color a points layer by gene expression
    # ------------------------------------------------------------------
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

        # Get full expression vector (one value per obs)
        x = full_adata[:, gene].X
        if sp.issparse(x):
            expr_full = x.toarray().ravel()
        else:
            expr_full = np.asarray(x).ravel()

        # Optional restriction to regionA/B
        if hasattr(self, "_keep_mask") and self._keep_mask is not None:
            mask = self._keep_mask
            if len(mask) == len(expr_full):
                expr = np.zeros_like(expr_full, dtype=float)
                expr[mask] = expr_full[mask]
            else:
                expr = expr_full.astype(float)
        else:
            expr = expr_full.astype(float)

        expr = np.asarray(expr, dtype=float)

        if not np.any(np.isfinite(expr)):
            show_error("Expression contains no finite values.")
            return

        points_layer_name = self.points_layer_combo.currentText()
        if not points_layer_name:
            show_error("No points layer selected to color.")
            return

        pl = self.viewer.layers[points_layer_name]
        if not isinstance(pl, Points) and not isinstance(pl, Shapes):
            show_error("Selected layer is not a Points/Shapes layer.")
            return

        # Length check
        n_points = len(pl.data)
        if len(expr) != n_points:
            show_error(
                f"Length mismatch: {len(expr)} expression values vs "
                f"{n_points} points in layer '{points_layer_name}'. "
                "Coloring aborted."
            )
            return

        # Create RGBA colors from colormap
        vmin = float(np.nanmin(expr))
        vmax = float(np.nanmax(expr))
        if vmin == vmax:
            vmax = vmin + 1e-6

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cmap = cm.get_cmap("viridis")
        rgba = cmap(norm(expr))  # shape (N, 4)

        # Assign direct colors
        pl.face_color = rgba
        pl.properties["expr"] = expr

        show_info("Colored points using explicit RGBA colormap.")


def two_region_de_widget(viewer: Viewer) -> TwoRegionDEWidget:
    """Factory function used by napari to create the dock widget."""
    return TwoRegionDEWidget(viewer)
