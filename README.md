# scRNA-seq Analysis Pipeline

Scripts for statistics and plotting figures in the single-cell RNA-seq study of brain tissues and gut bulk RNA-seq, covering differential expression analysis, cell-type marker visualization, GO enrichment, Ro/e tissue distribution analysis, and volcano plots.

## Repository Contents

- /fig1-6 # R scripts, Rmarkdown scripts and output figures
- figX/figX.Rmd # X is number 1-6, including the reproducible R scripts for each panel in figure
- figX/figX.html # Readability report by R markdown, include annotations, scripts and figures
- /script # General R scripts used in this study (shared themes, group orders and colors)
- environment.yml # Conda environment definition (name: `DEA`)
- install_packages.R # Installs the GitHub-hosted R packages

| File | Description |
| --- | --- |
| `fig1/plot_cell_correlation.R` | Cell-type average-expression correlation (Pearson/Spearman, corrplot + pheatmap) |
| `fig2/GVis.R` | `FindAllMarkers` + `ClusterGVis` heatmap / line-plot visualization with GO annotation |
| `fig3/DEA.R` | Differential expression analysis (`MMT_D30` vs `GF_D30`) per tissue and cell type using Seurat `FindMarkers` (Wilcoxon test) |
| `fig3/plot_Volcano.R` | Volcano plots of significant DEGs per tissue using `scMarkerViz` |
| `fig4/go_enrichment.R` | GO Biological Process enrichment of DEGs from selected neuronal cell types |
| `fig5/DEG_count_analysis.R` | edgeR differential expression of gut bulk RNA-seq per intestinal segment; bidirectional DEG-count bar plot |
| `fig6/plot_RoE.R` | Ro/e tissue-distribution analysis (`STARTRAC`) with heatmaps and a bubble plot |
| `script/stat_plot_functions.R` | Shared ggplot theme, group orders and color palettes for fig1-6 |

## Environment Setup

All analyses were performed in the conda environment `DEA`:

```bash
conda env create -f environment.yml
conda activate DEA
Rscript install_packages.R   # installs GitHub-hosted R packages
```

## Input Data

Place the following Seurat objects in the working directory
(adjust paths in the scripts if needed):

- `Brain.rds` — metadata columns required: `subtype`, `group`, `Tissue`
- `Brain_new.rds` — metadata columns required: `manual_celltype`, `group`
- `gut_bulkRNA_seq_TPM_20260801.csv` — TPM matrix of gut bulk RNA-seq (fig5)

## Workflow

1. **fig1/fig1.Rmd** — average expression per cell type, Pearson/Spearman correlation, mixed plots and clustered heatmaps.
2. **fig2/fig2.Rmd** — per-cell-type `FindAllMarkers`, `ClusterGVis` visualizations, and GO enrichment; outputs one folder per cell type.
3. **fig3/fig3.Rmd** — runs `FindMarkers` comparing `MMT_D30` vs `GF_D30` for every cell type within each tissue, then volcano plots of significant DEGs; outputs `MMT_D30_vs_GF_D30_<Tissue>_all_DEGs.csv`.
4. **fig4/fig4.Rmd** — GO BP enrichment of significant DEGs for GABAergic / Glutamatergic / Immature neurons; outputs `GO_results/`.
5. **fig5/fig5.Rmd** — edgeR exactTest of gut bulk RNA-seq per intestinal segment (`D30C vs D0`, `D30T vs D30C`); outputs `DEG_summary_by_segment.csv`, `DEG_table_perSegment.csv` and the bidirectional bar plot.
6. **fig6/fig6.Rmd** — Ro/e distribution analysis of glial subtypes (`Oligodendrocytes`, `Microglia`, `Astrocytes`); outputs `RoE_table.csv`.

To render the HTML reports (R and rmarkdown required):

```bash
cd fig1 && Rscript -e "rmarkdown::render('fig1.Rmd')"
```

## Notes

- `scMarkerViz` (used by `fig3/plot_Volcano.R` and `fig3/fig3.Rmd`) is not on CRAN or Bioconductor; it is installed from <https://github.com/DaoqiuWang/scMarkerViz> by `install_packages.R`.
- `jjAnno` (used by `fig2`) is installed from <https://github.com/junjunlab/jjAnno>.
- The scripts expect output files to be written to the working directory.
