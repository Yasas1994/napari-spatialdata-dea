"""Background worker for stratified differential expression analysis."""

from qtpy.QtCore import QObject, Signal
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


class DEAWorker(QObject):
    """Run DEA per group in a background QThread."""

    progress = Signal(int, int)
    finished = Signal(object)  # pd.DataFrame
    error = Signal(str)

    def __init__(
        self,
        full_adata,
        idxA: np.ndarray,
        idxB: np.ndarray,
        groupby_col: str,
        min_cells: int,
        method: str = "wilcoxon",
    ):
        super().__init__()
        self.full_adata = full_adata
        self.idxA = idxA
        self.idxB = idxB
        self.groupby_col = groupby_col
        self.min_cells = min_cells
        self.method = method.lower().strip()

    def run(self) -> None:
        try:
            adata = self.full_adata.copy()
            n = adata.n_obs
            region = np.full(n, "other", dtype=object)
            region[self.idxA] = "regionA"
            region[self.idxB] = "regionB"
            adata.obs["__region__"] = region

            if self.groupby_col not in adata.obs.columns:
                self.error.emit(
                    f"Group-by column '{self.groupby_col}' not found in adata.obs."
                )
                return

            groups = adata.obs[self.groupby_col].dropna().unique()
            total = len(groups)
            results = []

            for i, group in enumerate(groups):
                self.progress.emit(i, total)
                mask = (adata.obs[self.groupby_col] == group) & adata.obs[
                    "__region__"
                ].isin(["regionA", "regionB"])
                if mask.sum() == 0:
                    continue

                sub = adata[mask].copy()
                counts = sub.obs["__region__"].value_counts()
                if (
                    counts.get("regionA", 0) < self.min_cells
                    or counts.get("regionB", 0) < self.min_cells
                ):
                    continue

                sub.obs["__region__"] = pd.Categorical(
                    sub.obs["__region__"], categories=["regionA", "regionB"]
                )

                if self.method == "wilcoxon":
                    sc.tl.rank_genes_groups(
                        sub,
                        groupby="__region__",
                        groups=["regionA"],
                        reference="regionB",
                        method="wilcoxon",
                    )
                    df = sc.get.rank_genes_groups_df(sub, group="regionA")
                    df = df.rename(
                        columns={
                            "names": "gene",
                            "logfoldchanges": "logFC",
                            "pvals_adj": "p_adj",
                            "scores": "score",
                        }
                    )
                elif self.method == "deseq2":
                    df = self._run_deseq2(sub)
                    if df is None:
                        return
                else:
                    self.error.emit(f"Unknown DE method: {self.method}")
                    return

                df["group"] = str(group)
                results.append(df)

            self.progress.emit(total, total)

            if not results:
                self.error.emit(
                    "No group had enough cells in both regions."
                )
                return

            result_df = pd.concat(results, ignore_index=True)
            self.finished.emit(result_df)
        except Exception as e:
            self.error.emit(str(e))

    def _run_deseq2(self, sub) -> pd.DataFrame | None:
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
        except ImportError:
            self.error.emit(
                "pydeseq2 is not installed. Install it with: pip install pydeseq2"
            )
            return None

        X = sub.X
        if sp.issparse(X):
            X = X.toarray()
        X = np.asarray(X)

        if not np.all(X >= 0):
            self.error.emit(
                "DESeq2 requires non-negative counts. Switch to Wilcoxon or provide raw counts."
            )
            return None

        if not np.issubdtype(X.dtype, np.integer):
            if not np.allclose(X, X.astype(int)):
                self.error.emit(
                    "DESeq2 requires integer counts. Switch to Wilcoxon or provide raw counts."
                )
                return None

        # pydeseq2 expects counts as genes (rows) x samples (columns)
        counts_df = pd.DataFrame(
            X.T, index=sub.var_names, columns=sub.obs_names
        )
        obs_df = sub.obs[["__region__"]].copy()

        dds = DeseqDataSet(
            counts=counts_df,
            metadata=obs_df,
            design="~ __region__",
        )
        dds.deseq()

        stat_res = DeseqStats(
            dds, contrast=["__region__", "regionA", "regionB"]
        )
        stat_res.summary()
        df = (
            stat_res.results_df.reset_index()
            .rename(
                columns={
                    "index": "gene",
                    "log2FoldChange": "logFC",
                    "padj": "p_adj",
                    "stat": "score",
                }
            )
        )
        return df
