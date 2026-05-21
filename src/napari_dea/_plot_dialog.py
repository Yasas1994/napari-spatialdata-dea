"""Popup dialog showing a violin plot for a single gene in a single group."""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from qtpy.QtWidgets import QDialog, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class ViolinPlotDialog(QDialog):
    """Non-modal dialog with a matplotlib violin plot."""

    def __init__(
        self,
        adata,
        gene: str,
        group: str,
        groupby_col: str,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle(f"DE Violin Plot: {gene} in {group}")
        self.resize(550, 420)

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.fig = Figure(figsize=(5.5, 4.2))
        self.canvas = FigureCanvas(self.fig)
        layout.addWidget(self.canvas)

        self._plot(adata, gene, group, groupby_col)

    def _plot(self, adata, gene: str, group: str, groupby_col: str) -> None:
        if gene not in adata.var_names:
            ax = self.fig.add_subplot(111)
            ax.text(
                0.5,
                0.5,
                f"Gene '{gene}' not found in var_names.",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_axis_off()
            self.fig.tight_layout()
            self.canvas.draw()
            return

        mask = (adata.obs[groupby_col] == group) & adata.obs[
            "__region__"
        ].isin(["regionA", "regionB"])
        sub = adata[mask]

        if sub.n_obs == 0:
            ax = self.fig.add_subplot(111)
            ax.text(
                0.5,
                0.5,
                "No data available\nfor this group/region combination.",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_axis_off()
            self.fig.tight_layout()
            self.canvas.draw()
            return

        x = sub[:, gene].X
        if sp.issparse(x):
            expr = x.toarray().ravel()
        else:
            expr = np.asarray(x).ravel()

        regions = sub.obs["__region__"].values
        df = pd.DataFrame({"expression": expr, "region": regions})

        data_a = df[df["region"] == "regionA"]["expression"].values
        data_b = df[df["region"] == "regionB"]["expression"].values

        ax = self.fig.add_subplot(111)

        # Violin plot
        parts = ax.violinplot(
            [data_a, data_b],
            positions=[1, 2],
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        colors = ["#1f77b4", "#ff7f0e"]
        for pc, color in zip(parts["bodies"], colors):
            pc.set_facecolor(color)
            pc.set_alpha(0.5)

        # Jittered strip overlay
        rng = np.random.default_rng(42)
        jitter_a = rng.normal(1.0, 0.04, size=len(data_a))
        jitter_b = rng.normal(2.0, 0.04, size=len(data_b))
        ax.scatter(
            jitter_a,
            data_a,
            alpha=0.6,
            s=12,
            color=colors[0],
            edgecolors="white",
            linewidths=0.3,
        )
        ax.scatter(
            jitter_b,
            data_b,
            alpha=0.6,
            s=12,
            color=colors[1],
            edgecolors="white",
            linewidths=0.3,
        )

        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Region A", "Region B"])
        ax.set_ylabel("Expression")
        ax.set_title(
            f"{gene} in {group}\n(n_A={len(data_a)}, n_B={len(data_b)})"
        )

        self.fig.tight_layout()
        self.canvas.draw()
