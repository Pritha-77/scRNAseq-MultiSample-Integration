setwd("E:/scRNA-seq/scRNA")

library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(DoubletFinder)
library(scDblFinder)
library(SingleR)
library(celldex)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(SoupX)
library(scran)
library(ggridges)
library(decontX)
library(celda)
library(SeuratData)
library(presto)

#--------------------------#
# 1. Load Data+ Merge
#--------------------------#
#tar -xvf GSE180665_RAW.tar
#for i in *.gz; do tar -xvzf $i ; done

# get data location
dirs <- list.dirs(path = 'E:/scRNA-seq/trail.1', recursive = F, full.names = F)

dirs

for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  cts <- ReadMtx(mtx = paste0('E:/scRNA-seq/trail.1/',x,'/matrix.mtx.gz'),
                 features = paste0('E:/scRNA-seq/trail.1/',x,'/features.tsv.gz'),
                 cells = paste0('E:/scRNA-seq/trail.1/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}

ls()

# merge datasets

seu <- merge(HB17_background, y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor, HB53_background,
                                              HB53_tumor),
                       add.cell.ids = ls()[3:9],
                       project = 'HB')


seu

View(seu@meta.data)
# create a sample column
seu$sample <- rownames(seu@meta.data)

# split sample column
seu@meta.data <- separate(seu@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcode'), 
                                    sep = '_')

View(seu@meta.data)

#Join layers (Seurat v5)

seu <- JoinLayers(seu)


#--------------------------#
# 2. QC: Mitochondrial content
#--------------------------#

mt_gene <- c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
             "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB")
mt_gene_det <- mt_gene[mt_gene %in% rownames(seu)]
seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mt_gene_det)
summary(seu$percent.mt)


seu[["percent.mt2"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
summary(seu$percent.mt2)

library(ggplot2)
dat <- data.frame(
  byPattern = seu$percent.mt,
  byGene    = seu$percent.mt2,
  stringsAsFactors = FALSE
)
cor_val <- cor.test(dat$byPattern, dat$byGene, method = "spearman")
ggplot(dat, aes(x = byPattern, y = byGene)) + geom_point() + geom_smooth() + labs(x = "% of MT, est b
 y genes",
                                                                                  y = "% of MT,est by pattern", subtitle = paste0("rho=", round(cor_val$estimate,
                                                                                                                                             3), "; p-value=", cor_val$p.value[1])) + theme_classic()

#_________________OR_______________________#  

# calculate mitochondrial percentage
seu$percent.mt <- PercentageFeatureSet(seu, pattern='^MT-')

# explore QC
# Calculate % mitochondrial gene content
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.mt"]]

#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts


# Visualize QC metrics as a violin plot
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



# 1. Data Extraction and Transformation
# Fetch the features and the current cell identity (ident) from the Seurat object.
# 'ident' is the current clustering/grouping Seurat is using (e.g., sample ID, cluster ID).
plot_data <- FetchData(
  object = seu,
  vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "ident")
) %>%
  # Rename 'ident' for clarity in the plot
  rename(Cell_Identity = ident) %>%
  # Convert data from wide format (one column per feature) to long format.
  # This is essential for using facet_wrap in ggplot2.
  pivot_longer(
    cols = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    names_to = "Feature",
    values_to = "Value"
  )

# 2. ggplot2 Plot Generation
# Create the ggplot object
qc_plot <- ggplot(
  plot_data,
  aes(x = Cell_Identity, y = Value, fill = Cell_Identity)
) +
  # Create the violin plot. 'scale = "width"' ensures the maximum width is the same.
  # 'draw_quantiles = c(0.5)' adds the median line.
  geom_violin(scale = "width", draw_quantiles = c(0.5), alpha = 0.7) +
  
  # Optionally add a jittered point layer to show individual cells
  geom_jitter(
    height = 0, 
    width = 0.1, 
    size = 0.1, 
    alpha = 0.5,
    color = "black" # Use a neutral color for jittered points
  ) +
  
  # Arrange the plots in a grid with 3 columns.
  # 'scales = "free_y"' allows each feature to have its own appropriate y-axis scale.
  facet_wrap(~ Feature, scales = "free_y", ncol = 3) +
  
  # Apply a clean theme and set labels
  theme_minimal(base_size = 14) +
  labs(
    x = "Cell Identity / Cluster",
    y = "Observed Value",
    title = "Quality Control Metrics Distribution",
    subtitle = "Facet-wrapped Violin Plots replicating Seurat's VlnPlot"
  ) +
  
  # Adjust aesthetics
  theme(
    legend.position = "none", # Hide legend since fill is redundant with x-axis
    # Rotate x-axis labels for readability if Cell_Identity names are long
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold")
  )

# Print the plot
print(qc_plot)


## FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter low-quality cells (QC thresholds) #thrshold Seurat Tutorial
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#--------------------------#
# 3. Normalization & Feature Selection
#--------------------------#
seu <- NormalizeData(seu, normalization.method = "LogNormalize")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu, resolution = 0.5)


P1 <- DimPlot(seu, group.by = "Patient", pt.size = 0.2)
P2 <- DimPlot(seu, group.by = "Type", pt.size = 0.2)

P1+P2

#--------------------------#
# 4. DecontX (Ambient RNA correction)
#--------------------------#


var_genes <- VariableFeatures(seu)

sce <- as.SingleCellExperiment(seu)
sce <- sce[var_genes, ]

keep <- rowSums(counts(sce) > 0) >= 5
sce <- sce[keep, ]


# --- RUN DECONTX PER PATIENT ---

patients <- unique(seu$Patient)
patient_cells <- split(colnames(seu), seu$Patient)

sce_list <- list()

for (p in patients) {
  message("Running DecontX for patient: ", p)
  
  cells <- patient_cells[[p]]
  sce_sub <- sce[, cells]
  
  z_sub <- seu$seurat_clusters[cells]
  
  sce_decont_sub <- decontX(sce_sub, z = z_sub)
  
  # ensure correct naming
  colnames(sce_decont_sub) <- cells
  
  # save into list
  sce_list[[p]] <- sce_decont_sub
}

sce_list
lapply(sce_list, ncol)

contam_vec <- unlist(
  lapply(sce_list, function(sce_obj) {
    vec <- colData(sce_obj)$decontX_contamination
    names(vec) <- colnames(sce_obj)   # <-- keep real cell names
    vec
  })
)

head(names(contam_vec))
head(colnames(seu))

clean_names <- sub("^[^.]+\\.", "", names(contam_vec))
names(contam_vec) <- clean_names

head(names(contam_vec))

seu$decontX_contamination <- contam_vec[colnames(seu)]



# Combine decontX counts
dec_mat_list <- lapply(sce_list, function(sce_obj) {
  mat <- decontXcounts(sce_obj)
  colnames(mat) <- colnames(sce_obj)
  mat
})

decont_counts <- do.call(cbind, dec_mat_list)
decont_counts <- decont_counts[, colnames(seu)]  # reorder to match Seurat


#Clean RNA assay and rebuild assays container

rna_assay <- seu@assays$RNA
seu@assays <- list(RNA = rna_assay)   # reset assays to clean state


# 3. Add decontX assay

seu[["decontX"]] <- CreateAssayObject(counts = decont_counts)

# Verify
names(seu@assays)
Assays(seu)


# 4. Normalize, find HVGs, scale

DefaultAssay(seu) <- "decontX"
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)



############################################################
### 1) Define cancer markers (Tumor + Immune + Stromal)  ###
############################################################

markers <- c(
  # Malignant / epithelial tumor cells
  "EPCAM","KRT8","KRT18","KRT19","KRT5","MKI67","TOP2A",
  
  # T cells
  "CD3D","CD3E","CD3G","CD4","CD8A","CD8B",
  
  # B cells / Plasma cells
  "MS4A1","CD79A","CD79B","MZB1",
  
  # NK cells
  "GNLY","NKG7","KLRD1",
  
  # Myeloid / Macrophages / Monocytes
  "LYZ","S100A8","S100A9","FCGR3A","MS4A7","CD14","APOE",
  
  # Endothelial cells
  "PECAM1","VWF","KDR",
  
  # Fibroblasts / CAFs
  "COL1A1","COL1A2","COL3A1","ACTA2","FAP"
)


############################################################
### 2) Helper function: % cells expressing each marker    ###
############################################################

percent_expr_by_cluster <- function(seu, assay_name = "RNA", markers) {
  DefaultAssay(seu) <- assay_name
  expr <- GetAssayData(seu, slot = "counts")
  clust <- seu$seurat_clusters
  
  res <- lapply(split(seq_len(ncol(expr)), clust), function(idx) {
    sapply(markers, function(g) mean(expr[g, idx] > 0, na.rm = TRUE))
  })
  
  mat <- do.call(cbind, res)
  colnames(mat) <- paste0("C", sort(unique(clust)))
  rownames(mat) <- markers
  return(mat)
}


############################################################
### 3) Percent expressing: Before vs After DecontX        ###
############################################################

pct_before <- percent_expr_by_cluster(seu, assay_name = "RNA", markers = markers)

markers_exist <- markers[markers %in% rownames(seu[["decontX"]])]
markers_missing <- setdiff(markers, rownames(seu[["decontX"]]))

markers_exist
markers_missing

pct_after <- percent_expr_by_cluster(seu, assay_name = "decontX", markers = markers_exist)

library(dplyr)
library(tidyr)
library(ggplot2)

df_before <- as.data.frame(pct_before) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cluster", values_to = "pct") %>%
  mutate(stage = "before")

df_after <- as.data.frame(pct_after) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cluster", values_to = "pct") %>%
  mutate(stage = "after")

df_comb <- bind_rows(df_before, df_after)

ggplot(df_comb, aes(x = cluster, y = pct, fill = stage)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  theme_minimal() +
  labs(y = "Fraction of cells expressing (>0 UMI)", x = "Cluster")


############################################################
### 4) Violin plots for selected cancer markers           ###
############################################################

plot_markers <- c("CD8B","KRT5","TOP2A","NKG7","PECAM1","FAP")

df_before <- FetchData(seu, vars = c("seurat_clusters", plot_markers)) %>%
  rownames_to_column("cell") %>%
  pivot_longer(cols = all_of(plot_markers), names_to = "gene", values_to = "expr") %>%
  mutate(stage = "before")

df_after <- FetchData(seu, vars = c("seurat_clusters", plot_markers), assay = "decontX") %>%
  rownames_to_column("cell") %>%
  pivot_longer(cols = all_of(plot_markers), names_to = "gene", values_to = "expr") %>%
  mutate(stage = "after")

df_both <- bind_rows(df_before, df_after)

ggplot(df_both, aes(x = seurat_clusters, y = expr, fill = stage)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_classic() +
  labs(y = "Expression", x = "Cluster")


#--------------------------#
# 5. Cell Cycle Scoring
#--------------------------#

feat_s <- cc.genes$s.genes
feat_g2m <- cc.genes$g2m.genes
feat_s
feat_g2m

# The CellCycleScoring() function the scores each cell based on speci c features for S-phase/G2M-phase. For a given cell with signifcant high S.Score or G2M.Score they are assigned as S/G2M. 
#Cells with low both S.Score and G2M.Score were assigned as G1.
#This will overestimate cells in S-/G2M-phases in the tissues with low cell cycle i.e. neurons.

DefaultAssay(seu) <- "RNA"
seu <- CellCycleScoring(seu, s.features = feat_s, g2m.features = feat_g2m)

dat_s <- data.frame(cell_id = Cells(seu), cat = "S_Score", Phase = seu$Phase,
                    score = seu$S.Score)
dat_g2m <- data.frame(cell_id = Cells(seu), cat = "G2M_Score", Phase = seu$Phase,
                      score = seu$G2M.Score)
dat <- rbind(dat_s, dat_g2m)
dat$Phase <- factor(dat$Phase, levels = c("G1", "S", "G2M"))

ggplot(dat, aes(x = Phase, y = score, fill = Phase)) + geom_violin() + labs(x = "",
                                                                            y = "Score", fill = "Phase") + facet_wrap(~cat) + theme_classic()

#--------------------------#
# 6. Doublet Detection
#--------------------------#

library(scDblFinder)
library(SingleCellExperiment)

patients <- unique(seu$Patient)
seu$doublet <- NA
seu$doublet_score <- NA

for (p in patients) {
  cells <- colnames(seu)[seu$Patient == p]
  sce_sub <- as.SingleCellExperiment(seu[, cells])
  sce_sub <- scDblFinder(sce_sub, clusters = sce_sub$seurat_clusters, samples = sce_sub$Patient)
  
  seu$doublet[cells] <- colData(sce_sub)$scDblFinder.class
  seu$doublet_score[cells] <- colData(sce_sub)$scDblFinder.score
}

table(seu$doublet == "TRUE" | seu$percent.mt >= 10)

# Create a logical column for singlets
seu$singlet <- seu$doublet == "1"   # TRUE for singlets

# Filter cells: singlets and percent.mt < 10%
seu <- subset(seu, subset = singlet & percent.mt < 10) #or 5 (standard seruat protocol)

# Check
seu
table(seu$singlet)
summary(seu$percent.mt)

##____________USE ANY ONE OF THE INTEGRATION METHOD Depending on your dataset- 7.A/7.B/7.C___________________________##


#--------------------------#
# 7.A Harmony Integration
#--------------------------#

# run Harmony -----------
seu_integrated <- seu %>%
  RunHarmony(group.by.vars = 'Patient', plot_convergence = FALSE)

seu_integrated@reductions

seu_integrated.embed <- Embeddings(seu_integrated, "harmony")
seu_integrated.embed[1:10,1:10]


#--------------------------#
# 7.B RPCA Integration
#--------------------------#
# Split your merged Seurat object back to list of samples if needed
seu_list <- SplitObject(seu, split.by = "Patient")

# Normalize, find HVGs within each object
for (i in 1:length(seu_list)) {
  seu_list[[i]] <- NormalizeData(seu_list[[i]])
  seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method = "vst", nfeatures = 3000)
}

# Select features for integration
features <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = 3000)

# Prepare objects for RPCA integration
seu_list <- lapply(seu_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find anchors using RPCA
anchors <- FindIntegrationAnchors(object.list = seu_list,
                                  anchor.features = features,
                                  reduction = "rpca")  # key step: use RPCA

# Integrate data
seu_integrated <- IntegrateData(anchorset = anchors)

# Switch default assay to integrated
DefaultAssay(seu_integrated) <- "integrated"

# Scale data (optional)
seu_integrated <- ScaleData(seu_integrated, verbose = FALSE)

# PCA, UMAP, clustering
seu_integrated <- RunPCA(seu_integrated, npcs = 50)
seu_integrated <- FindNeighbors(seu_integrated, dims = 1:30)
seu_integrated <- FindClusters(seu_integrated, resolution = 0.5)
seu_integrated <- RunUMAP(seu_integrated, dims = 1:30)

# Visualize
DimPlot(seu_integrated, group.by = "Patient")
DimPlot(seu_integrated, group.by = "Type")


#--------------------------#
# 7.c SCTransform
#--------------------------#

# ---- Recommended: SCT-based integration (Seurat v5 compatible) ----
# Split by patient and run SCTransform per sample

obj.list <- SplitObject(seu, split.by = 'Patient')

for (i in seq_along(obj.list)) {
  # SCTransform is robust for batch correction and compatible with v5 assays
  obj.list[[i]] <- SCTransform(obj.list[[i]], verbose = FALSE)
}

# Select features for integration (increase if needed; 3000 is a common choice)
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)

# Prepare objects for SCT integration (required in v5)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

# Find anchors using RPCA for memory-efficiency; reduce dims to 20 to save memory
anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  dims = 1:20
)

# Integrate (SCT) using the anchors (also specify dims)
seu_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)


#--------------------------#
# 8. Regress out confounders
#--------------------------#

pot_conf <- c(
  "percent.mt",
  "decontX_contamination",
  "doublet_score",
  "S.Score",
  "G2M.Score"
)

seu_filt <- ScaleData(seu_integrated, vars.to.regress = pot_conf)


#--------------------------#
# 9. Dimensionality Reduction
#--------------------------#
# Perform PCA with RunPCA() function in Seurat. Set Principle components to 50.
set.seed(1001)
DefaultAssay(seu_filt) <- "RNA"
seu_filt <- RunPCA(seu_filt, assay = "RNA", npcs = 50)

#Evaluate dimensionality with elbow plot
ElbowPlot(seu_filt, ndims = 50, reduction = "pca")

#We can use cutoffs to derive this computationally based on: * How low variation is. * How much variation has already been described. * The rate of change in variation between PCs.

library(dplyr)
# Determine percent of variation associated with each PC
pct <- seu_filt[["pca"]]@stdev/sum(seu_filt[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % 
#variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC, last
# point where change of % of variation is more than 0.1%.
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] +
  1
pc <- min(co1, co2)
pc

plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pc)) + geom_text() +
  geom_vline(xintercept = 90, color = "grey") + geom_hline(yintercept = min(pct[pct >
                                                                                  5]), color = "grey") + theme_bw()
#Make non-linear dimensional reduction
seu_filt <- FindNeighbors(seu_filt, dims = 1:pc, reduction = "pca")
seu_filt <- RunUMAP(seu_filt, dims = 1:pc, reduction = "pca")


# UMAP
DimPlot(seu_filt)

#Clustering

seu_filt <- FindClusters(seu_filt, resolution = 0.5)

seu_filt[["cluster_byDefault"]] <- seu_filt$seurat_clusters


DimPlot(seu_filt, group.by = "seurat_clusters", label = TRUE, pt.size = 0.2) + NoLegend()




P1 <- DimPlot(seu_filt, group.by = "Patient", pt.size = 0.2)
P2 <- DimPlot(seu_filt, group.by = "Type", pt.size = 0.2)

P1+P2


#Testing resolutions
library(clustree)
reso <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
reso_res <- lapply(1:length(reso), function(x, seu_filt, reso) {
  seu_filt <- FindClusters(seu_filt, resolution = reso[x])
  clust <- setNames(seu_filt$seurat_clusters, Cells(seu_filt))
  return(clust)
}, seu_filt, reso)


names(reso_res) <- paste0("k", 1:length(reso))

k_tab <- do.call(cbind, reso_res)
k_dat <- as.data.frame(k_tab)
head(k_dat, 2)

clustree(k_dat, prefix = "k", node_colour = "sc3_stability")

seu_filt <- FindClusters(seu_filt, resolution = 0.4)

DimPlot(seu_filt, group.by = "seurat_clusters", label = TRUE, pt.size = 0.2) + NoLegend() +
  ggtitle("Optimized Clusters")

DimPlot(seu_filt, group.by = "cluster_byDefault", label = TRUE, pt.size = 0.2) +
  NoLegend() + ggtitle("Default Clusters")


P3 <- DimPlot(seu_filt, group.by = "Patient", pt.size = 0.2)
P4 <- DimPlot(seu_filt, group.by = "Type", pt.size = 0.2)

P3+P4

library(ggplot2)
library(dplyr)

# Gene of interest
gene_of_interest <- "EPCAM"  # replace with your gene

# Extract metadata + gene expression
df <- FetchData(seu_filt, vars = c("Patient", "Type", gene_of_interest))

# Summarize mean expression per patient and Type
df_summary <- df %>%
  group_by(Patient, Type) %>%
  summarise(mean_expr = mean(get(gene_of_interest)), .groups = "drop")

# Bar plot
ggplot(df_summary, aes(x = Patient, y = mean_expr, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +  # side-by-side bars
  theme_classic() +
  labs(x = "Patient", y = paste("Mean expression of", gene_of_interest), fill = "Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#_______Bar plot of cells per cluster per Type


# Make a table of cell counts per cluster and per Type
tbl <- table(seu_filt$Type, seu_filt$seurat_clusters)

# Barplot
barplot(tbl, beside = TRUE,
        main = "Cell numbers in each cluster of each Type",
        xlab = "Clusters",
        ylab = "Number of cells",
        col = rainbow(nrow(tbl)),
        legend = rownames(tbl))


#==========================#
# 10. Identify cluster-specfic marker genes
#==========================#
#we are selecting only positive markers (only.pos=TRUE), genes that are found in 0.25 of cells and have a logFC 0.25 di erence.
markers <- FindAllMarkers(seu_filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)

# Plot marker expression
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(seu_filt, features = top10$gene) + NoLegend()


#select the top 2 genes for each cluster by log2FC
top_genes <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
head(top_genes)

DoHeatmap(seu_filt, features = top_genes$gene) + NoLegend()


#Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(seu_filt, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(seu_filt, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seu_filt.markers <- FindAllMarkers(seu_filt, only.pos = TRUE)
seu_filt.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


# Identify markers for all clusters
seu_filt.markers <- FindAllMarkers(seu_filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Save marker results
write.table(seu_filt.markers,file="Markers_info.csv",
            sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

cluster0.markers <- FindMarkers(seu_filt, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(seu_filt, features = c("SCGB2A2","H4C3","OR1A1"))
VlnPlot(seu_filt, features = c("OR1A1"))

FeaturePlot(seu_filt, features = "OR1A1")
FeaturePlot(seu_filt, features = "SCGB2A2")
FeaturePlot(seu_filt, features = "ENPP1")


FeaturePlot(seu_filt, features = c("OR1A1", "SCGB2A2", "ENPP1"))


seu_filt.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seu_filt, features = top10$gene) + NoLegend()



# Example: GAPDH in Seurat- In case of ENSMBL ID
ens_id <- rowData(sce2)$ID[rowData(sce2)$Symbol == "GAPDH"]
FeaturePlot(seu_filt, features = ens_id) + ggtitle("After QC")

## Example: GAPDH in Seurat
FeaturePlot(seu_filt, features = "GAPDH-a1") + ggtitle("Apter QC")


library(ggplot2)
library(dplyr)

# Gene of interest
gene_of_interest <- "EPCAM"  # replace with your gene

# Extract metadata + gene expression
df <- FetchData(seu_filt, vars = c("Patient", "cellType_byClust", gene_of_interest))

# Summarize mean expression per patient and cell type
df_summary <- df %>%
  group_by(Patient, cellType_byClust) %>%
  summarise(mean_expr = mean(get(gene_of_interest)), .groups = "drop")

# Bar plot
ggplot(df_summary, aes(x = Patient, y = mean_expr, fill = cellType_byClust)) +
  geom_bar(stat = "identity", position = "dodge") +  # side-by-side bars
  theme_classic() +
  labs(x = "Patient", y = paste("Mean expression of", gene_of_interest), fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#==========================#
# 11. Annotate clusters to known Markers
#==========================#

#________________Evaluate known marker genes expression________________#

# Known marker genes in cancer datasets

## Tumor cells: EPCAM, KRT8, KRT18
## Fibroblasts: COL1A1, PDGFRA
## Endothelial cells: PECAM1, VWF
## T-cells: CD3D, CD3E
## B-cells: MS4A1
## NK cells: GNLY, NKG7
## Monocytes / Macrophages: CD14, FCGR3A, CD68
## Dendritic cells: FCER1A, CST3
## Platelets: PPBP

# Evaluate marker gene expression: Tumor cells
known_marker <- c("EPCAM", "KRT8", "KRT18")
FeaturePlot(seu_filt, features = known_marker)

# Evaluate marker gene expression: Fibroblasts
known_marker <- c("COL1A1", "PDGFRA")
FeaturePlot(seu_filt, features = known_marker)

# Match markers to clusters
known_marker <- c("EPCAM", "KRT8", "KRT18", "COL1A1", "PDGFRA",
                  "PECAM1", "VWF", "CD3D", "CD3E", "MS4A1", "GNLY",
                  "NKG7", "CD14", "FCGR3A", "CD68", "FCER1A", "CST3", "PPBP")
mat <- GetAssayData(seu_filt, assay = "RNA", slot = "data")
mat <- mat[known_marker, ]
mat <- as.matrix(mat)

clust <- unique(seu_filt$seurat_clusters)
clust <- as.character(clust)
avgExp_byClust <- lapply(clust, function(clust, seu, known) {
  sub <- subset(seu, subset = seurat_clusters == clust)
  mat <- GetAssayData(sub, assay = "RNA", slot = "data")
  mat <- mat[known, ]
  mat <- as.matrix(mat)
  avg <- rowMeans(mat)
  return(avg)
}, seu_filt, known_marker)
names(avgExp_byClust) <- paste0("C", clust)

library(pheatmap)
avgExp_mat <- do.call(cbind, avgExp_byClust)
pheatmap(avgExp_mat, scale = "row")

# Assign cell types by clusters (update based on your dataset)
seu_filt[["cellType_byClust"]] <- NA
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(0,1)] <- "Tumor cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(2)] <- "Fibroblasts"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(3)] <- "Endothelial cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(4,5)] <- "T-cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(6)] <- "B-cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(7)] <- "NK cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(8)] <- "Monocytes/Macrophages"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(9)] <- "Dendritic cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(10)] <- "Platelets"
seu_filt$cellType_byClust[is.na(seu_filt$cellType_byClust)] <- "Unspecified"

table(seu_filt$cellType_byClust)

DimPlot(seu_filt, group.by = c("seurat_clusters", "cellType_byClust"), label = TRUE,
        pt.size = 0.2) + NoLegend()


DimPlot(seu_filt, group.by = c("seurat_clusters", "cellType_byClust"), 
        label = TRUE, pt.size = 0.2)


DimPlot(seu_filt, group.by = c("seurat_clusters", "cellType_byClust"), 
        label = FALSE, pt.size = 0.2) +
  theme(legend.position = "right")  # legend on the side


# Violin plot for specific markers (example: Tumor vs T-cell)
mark_gene <- c("EPCAM", "CD3E")
DefaultAssay(seu_filt) <- "RNA"

df <- FetchData(seu_filt, vars = c("seurat_clusters", mark_gene))
library(tidyr)

df_long <- pivot_longer(df, cols = all_of(mark_gene),
                        names_to = "gene", values_to = "expression")

library(ggplot2)
ggplot(df_long, aes(x = seurat_clusters, y = expression, fill = seurat_clusters)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~gene, scales = "free_y") +
  theme_classic()



#==========================#
#12. Differential gene expression
#==========================#

#___________________Differential gene expression between cell types_________________________#

#comparing Tumor cells to Fibroblasts
deg <- FindMarkers(seu_filt, group.by = "cellType_byClust", ident.1 = "Tumor cells",
                   ident.2 = "Fibroblasts", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1)

deg$Symbol <- rownames(deg)
deg <- deg[order(-deg$avg_log2FC), ]
head(deg)

# avg_log2FC: average log2_expression of ident.1 - average log2 expression of ident.2
#pct.1: percent of cells with a given gene expression (expression level > 0) in ident.1
#pct.2: percent of cells with a given gene expression (expression level > 0) in ident.2
#p_val: p value estimated by using wilcox method
#p_val_adj: adjust p value calculated with FDR


#Volcano Plot
deg$sig <- NA
deg$sig[deg$avg_log2FC >= 0.585 & deg$p_val_adj < 0.05] <- "Tumor cells"
deg$sig[deg$avg_log2FC <= -0.585 & deg$p_val_adj < 0.05] <- "Fibroblasts"
deg$sig[is.na(deg$sig)] <- "nc"
deg$sig <- factor(deg$sig, levels = c("Tumor cells", "nc", "Fibroblasts"))
table(deg$sig)
ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) + geom_point() +
  scale_color_manual(values = c("red", "grey", "blue")) + labs(x = "log2_Fold-Changes",
                                                               y = "-log10_adjPV", color = "") + theme_classic()

