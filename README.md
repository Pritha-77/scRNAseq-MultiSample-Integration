# MultiSample-scRNA-seq-Workflow

This repository contains an R script for comprehensive analysis of multi-sample single-cell RNA-seq (scRNA-seq) data using Seurat and complementary tools. The workflow supports data loading, QC, normalization, ambient RNA correction (DecontX), doublet detection, data integration (Harmony, RPCA, SCTransform), clustering, marker gene identification, cell type annotation, and differential expression analysis.

## Features

- Load and merge multiple samples from GEO dataset GSE180665
- Quality control, including mitochondrial content filtering and feature scatterplots
- Ambient RNA contamination correction using DecontX per patient
- Doublet detection with scDblFinder and filtering
- Batch correction and integration using Harmony, RPCA, or SCTransform methods
- Dimensionality reduction (PCA, UMAP) and clustering with clustree-based resolution optimization
- Identification of cluster markers and cell type annotation using canonical marker genes
- Differential gene expression analysis between cell types/clusters
- Visualization with violin plots, heatmaps, feature plots, and barplots

## Requirements

- R (version 4.x or higher)
- Packages: Seurat, SeuratDisk, DoubletFinder, scDblFinder, SingleR, celldex, patchwork, tidyverse, ggplot2, gridExtra, SoupX, scran, ggridges, decontX, celda, SeuratData, presto, DropletUtils, clustree, pheatmap, and dependencies

## Usage

1. Clone or download the repository:  
