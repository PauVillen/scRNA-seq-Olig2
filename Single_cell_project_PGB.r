#Olig 2 PGB Project - Single-cell part
setwd("/Users/pauvillen14/Desktop/BIOINFO/PGB/PROJECT/single-cell")

# ---------------------------------
# ---- Import needed libraries ----
# ---------------------------------
library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(readr)
library(data.table)


# --------------------------------------------
# ---- Read raw data and prepare metadata ----
# --------------------------------------------
#We have to reed only a random number of cells from the count matrix for memory issues
raw.data <- read.csv("matrix_subset.csv", row.names = 1)

raw.data <- Matrix(as.matrix(raw.data), sparse = TRUE)
raw.data <- t(raw.data) #now: 50281 genes, 15000 cells

meta.data <- read.csv("metadata.csv")
cell.names <- colnames(raw.data)
cell.meta.data <- meta.data[meta.data$sample_name %in% cell.names, ] #keep the same cells as in raw
cell.meta.data <- cell.meta.data[match(cell.names, cell.meta.data$sample_name), ] #put in the same order
rownames(cell.meta.data) <- cell.meta.data$sample_name

cell.meta.data$subclass_color[cell.meta.data$subclass_label == "Oligo"] <- "#000000"
cell.meta.data$subclass_color[cell.meta.data$subclass_label == "OPC"] <- "#32CD32"


# --------------------------------------------------------------
# ---- Check that reduced metadata has the same proportion -----
# ---- of cells across subclasses as the original metadata -----
# --------------------------------------------------------------

tab_full <- table(meta.data$subclass_label)
tab_sub  <- table(cell.meta.data$subclass_label)

# Absolute counts
data.frame(subclass = names(tab_full),
           full = as.numeric(tab_full[names(tab_full)]),
           sub  = as.numeric(tab_sub[names(tab_full)]))

# Relative proportions
prop_full <- prop.table(tab_full)
prop_sub  <- prop.table(tab_sub)
data.frame(subclass = names(prop_full),
           prop_full = round(100*as.numeric(prop_full),2),
           prop_sub  = round(100*as.numeric(prop_sub),2),
           diff_pct  = round(100*(as.numeric(prop_sub) - as.numeric(prop_full))/as.numeric(prop_full), 1)) 

# Barplot comparison
library(ggplot2)
df <- data.frame(
  subclass = rep(names(tab_full), 2),
  type = rep(c("Original Metadata","Reduced Metadata"), each = length(tab_full)),
  count = c(as.numeric(tab_full), as.numeric(tab_sub[names(tab_full)]))
)
ggplot(df, aes(x = reorder(subclass, -count), y = count, fill = type)) +
  geom_col(position="dodge") + theme(axis.text.x = element_text(angle=60, hjust=1)) +
  ylab("Cell count") + xlab("Subclass")



# --------------------------------------------
# ---- Create and prepare Seurat Object ----
# --------------------------------------------
tiss <- CreateSeuratObject(counts = raw.data)
tiss <- AddMetaData(object = tiss, cell.meta.data)

VlnPlot(tiss, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#filter:
tiss <- subset(tiss, subset = nFeature_RNA > 2000 & nFeature_RNA < 9000 &
                 nCount_RNA > 10000 & nCount_RNA < 80000)

VlnPlot(tiss, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

tiss <- NormalizeData(object = tiss, scale.factor = 1e4)


# --------------------------------------------
# ---- PCA ----
# --------------------------------------------
tiss <- FindVariableFeatures(object = tiss)
tiss <- ScaleData(object = tiss, features = VariableFeatures(tiss)) #scale only variable features to reduce memory usage
top10 <- head(VariableFeatures(tiss), 10)
top100 <- head(VariableFeatures(tiss), 100)
plot1 <- VariableFeaturePlot(tiss)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

tiss <- RunPCA(tiss, features = VariableFeatures(object = tiss))
ElbowPlot(tiss, ndims = 50)

# --------------------------------------------
# ---- Find Clusters ----
# --------------------------------------------
tiss <- FindNeighbors(tiss, dims = 1:40) 
tiss <- FindClusters(tiss, resolution = 0.5) 


# --------------------------------------------
# ---- UMAP to visualize clusters ----
# --------------------------------------------
colors <- tiss@meta.data %>%
  distinct(subclass_label, subclass_color) %>%
  filter(!is.na(subclass_color))

color_vec <- setNames(colors$subclass_color, colors$subclass_label)

set.seed(345)
tiss <- RunUMAP(tiss, dims = 1:40)
p1 <- DimPlot(tiss, reduction = "umap", group.by = 'class_label')
p2 <- DimPlot(tiss, reduction = "umap", group.by = 'seurat_clusters', label = TRUE)
p3 <- FeaturePlot(tiss, features = "OLIG2")
p2 + p3
VlnPlot(tiss, features="OLIG2", group.by = "seurat_clusters")
DotPlot(tiss, features= "OLIG2", group.by = "seurat_clusters")


# --------------------------------------------
# ---- Find marker genes ----
# --------------------------------------------
cluster2.markers <- FindMarkers(tiss, ident.1 = 29, min.pct = 0.25, test.use= )

# --------------------------------------------
# ---- Find out real annotation ----
# --------------------------------------------
tiss$cell_type_annotation <- tiss@meta.data$cell_type_alias_label
tiss$cell_type_annotation[is.na(tiss$cell_type_annotation)] <- "unknown"
tiss$cell_type_annotation <- as.factor(tiss$cell_type_annotation)
DimPlot(tiss, reduction = "umap", group.by = "cell_type_annotation", label = FALSE) +
  theme(legend.text = element_text(size = 5))



# --------------------------------------------
# ---- AUCell ----
# --------------------------------------------
library(GSEABase)
exprMatrix <- GetAssayData(tiss, slot = "counts")
exprMatrix <- as(exprMatrix, "dgCMatrix")

genes <- read.csv("adjusted_genes_with_ensembl.csv")
gene_ids <- unique(genes$geneID)

gene_ids_upper <- toupper(gene_ids)
rownames_upper <- toupper(rownames(exprMatrix))

# Get intersection
intersect_genes <- intersect(gene_ids_upper, rownames_upper)
# Match the uppercase intersection back to the original rownames
matched_genes <- rownames(exprMatrix)[toupper(rownames(exprMatrix)) %in% intersect_genes]

geneSet <- GeneSet(matched_genes, setName = "OLIG2_related_genes")
library(AUCell)
cells_AUC <- AUCell_run(exprMatrix, geneSet)
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)
set.seed(333)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
auc_values <- getAUC(cells_AUC)
tiss$OLIG2_AUC <- as.numeric(auc_values["OLIG2_related_genes", ])
FeaturePlot(tiss, features = "OLIG2_AUC") +
  ggtitle("AUCell enrichment scores for OLIG2-related gene set") +
  theme(plot.title = element_text(hjust = 0.5))
VlnPlot(tiss, features="OLIG2_AUC", group.by = "seurat_clusters")
DotPlot(tiss, features= "OLIG2_AUC", group.by = "seurat_clusters")

