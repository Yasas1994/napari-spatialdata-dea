# napari-spatialdata-dea

A napari plugin for **stratified differential expression analysis** between spatial regions of interest.

![GitHub](https://img.shields.io/github/license/Yasas1994/napari-spatialdata-dea)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/Yasas1994/napari-spatialdata-dea/main?color=8a35da)

> [!NOTE]
> napari-spatialdata-dea [documentation](https://napari-spatialdata-dea.readthedocs.io/en/latest/) will be online soon!

## Features

- **Stratified DEA** — run differential expression per group (e.g. cell type, treatment, condition) by selecting any `obs` column in your AnnData table.
- **Two statistical methods** — Wilcoxon rank-sum (via `scanpy`) or DESeq2 (via optional `pydeseq2`).
- **Multiple ROIs per condition** — draw several polygons in a single Shapes layer; they are automatically merged before querying.
- **Smart visualization** — color only the bins inside Region A, Region B, their union, or the whole slide.
- **Violin plots** — double-click any DE result to open a non-modal violin plot comparing expression in Region A vs Region B for that group and gene.
- **Background threading** — DE computation runs in a `QThread` so napari stays responsive.
- **Auto-detect target layer** — the plugin automatically finds the correct Points/Shapes layer to color based on the selected table.

## Installation

Make sure you have both [spatialdata and napari-spatialdata](https://spatialdata.scverse.org/en/latest/installation.html) installed in the same environment, then run:

```bash
pip install git+https://github.com/Yasas1994/napari-spatialdata-dea.git@main
```

### Optional: DESeq2 support

To enable the DESeq2 method, install the optional dependency:

```bash
pip install pydeseq2
```

## Quick start

1. **Load your data** in napari via `napari-spatialdata`.
2. **Draw regions of interest** as Shapes layers (see the [napari ROI tutorial](https://spatialdata.scverse.org/en/latest/tutorials/notebooks/notebooks/examples/napari_rois.html)). Each Shapes layer can contain multiple polygons — they will be merged automatically.
3. **Open the plugin** from the napari Plugins menu: **Differential Expression Analysis**.
4. **Configure the analysis** in the four sections of the widget:
   - **Data Selection** — pick the table and the two ROI layers.
   - **DE Parameters** — choose the grouping column, statistical method, and minimum cells per region.
   - **Execution** — click **Run DE** and wait for the progress bar to fill.
   - **Visualization & Results** — browse results, filter by group, and visualize genes.

### Interacting with results

| Action | Effect |
|--------|--------|
| **Single-click** a row | Copies the gene name to the *Gene to visualize* field |
| **Double-click** a row | Opens a violin plot for that (group, gene) pair |
| **Color region** dropdown | Choose *Whole slide*, *Only A*, *Only B*, or *Union (A ∪ B)* before clicking **Color by Gene** |
| **Filter results by group** | Narrow the table to one group without re-running DE |

## Architecture

The plugin consists of three modules under `src/napari_dea/`:

| File | Purpose |
|------|---------|
| `_widget.py` | Main `TwoRegionDEWidget` — UI layout, thread orchestration, and layer coloring |
| `_worker.py` | `DEAWorker` — background `QThread` that loops over groups and runs `rank_genes_groups` or `pydeseq2` |
| `_plot_dialog.py` | `ViolinPlotDialog` — non-modal matplotlib popup with jittered violin plots |

## Requirements

- Python ≥ 3.10
- napari ≥ 0.4
- spatialdata ≥ 0.2
- napari-spatialdata
- scanpy
- pandas, numpy, shapely, geopandas
- matplotlib (bundled with napari)
- **Optional:** pydeseq2 (for DESeq2 method)

## License

BSD 3-Clause
