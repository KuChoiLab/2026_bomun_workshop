#1. Loading necessary packages and practice data
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(patchwork)
library(EnhancedVolcano)
library(msigdbr)
library(clusterProfiler)

hs09 <- readRDS("hs09.rds")
hs12 <- readRDS("hs12.rds")
hs13 <- readRDS("hs13.rds")
hs254 <- readRDS("hs254.rds")
hs255 <- readRDS("hs255.rds")
hs266 <- readRDS("hs266.rds")

#2. Remove Doublets

#2-1. Convert your Seurat objects to SingleCellExperiment
sce09 <- as.SingleCellExperiment(hs09, assay = "RNA")
sce12 <- as.SingleCellExperiment(hs12, assay = "RNA")
sce13 <- as.SingleCellExperiment(hs13, assay = "RNA")
sce254 <- as.SingleCellExperiment(hs254, assay = "RNA")
sce255 <- as.SingleCellExperiment(hs255, assay = "RNA")
sce266 <- as.SingleCellExperiment(hs266, assay = "RNA")

#2-2. Run scDblFinder
sce09 <- scDblFinder(sce09)
sce12 <- scDblFinder(sce12)
sce13 <- scDblFinder(sce13)
sce254 <- scDblFinder(sce254)
sce255 <- scDblFinder(sce255)
sce266 <- scDblFinder(sce266)

#2-3. Add the results back to your Seurat object
hs09$doublet_score <- sce09$scDblFinder.score
hs09$is_doublet <- sce09$scDblFinder.class

hs12$doublet_score <- sce12$scDblFinder.score
hs12$is_doublet <- sce12$scDblFinder.class 

hs13$doublet_score <- sce13$scDblFinder.score
hs13$is_doublet <- sce13$scDblFinder.class 

hs254$doublet_score <- sce254$scDblFinder.score
hs254$is_doublet <- sce254$scDblFinder.class 

hs255$doublet_score <- sce255$scDblFinder.score
hs255$is_doublet <- sce255$scDblFinder.class 

hs266$doublet_score <- sce266$scDblFinder.score
hs266$is_doublet <- sce266$scDblFinder.class 

#2-4 Check how many cells were detected as doublets 
table(hs09$is_doublet)
table(hs12$is_doublet)
table(hs13$is_doublet)
table(hs254$is_doublet)
table(hs255$is_doublet)
table(hs266$is_doublet)

#2-5 Filtering doublet cells
hs09 <- subset(hs09, subset = is_doublet == "singlet")
hs12 <- subset(hs12, subset = is_doublet == "singlet")
hs13 <- subset(hs13, subset = is_doublet == "singlet")
hs254 <- subset(hs254, subset = is_doublet == "singlet")
hs255 <- subset(hs255, subset = is_doublet == "singlet")
hs266 <- subset(hs266, subset = is_doublet == "singlet")

#2-6 Merge separate datas into one for convenience in further analysis
data <- merge(
  x = hs09, 
  y = c(hs12, hs13, hs254, hs255, hs266),
  add.cell.ids = c("Hs09", "Hs12", "Hs13", "Hs254", "Hs255", "Hs266"),
)

#3. Quality control

#3-1. Take a look at such metrics within our data
VlnPlot(data, group.by = 'sample',feature = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
range(data$nCount_RNA)
range(data$nFeature_RNA)
range(data$percent.mt)

#3-2. Filtering cells with high mitochondrial percentage
data <- subset(data, subset = percent.mt < 20 & nCount_RNA >= 200)

#4. Normalization
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

#5. Feature selection & Scaling
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

#5-1. Plotting variable features
top10 <- head(VariableFeatures(data), 10)
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#6. Scaling the data
data <- ScaleData(data)

#7. Run linear dimensionality reduction(PCA)
data <- RunPCA(data, verbose = FALSE)

#7-1. Examine and visualize PCA results in a few different ways
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

#7-2. Elbow plot (Determining # of PCs to use for further analysis)
#7-2-1. Extract the data from your Seurat object
stdev <- data[["pca"]]@stdev
pct <- stdev / sum(stdev) * 100
cum_pct <- cumsum(pct)
pcs <- 1:length(stdev)
plot_df <- data.frame(PC = pcs, stdev = stdev, pct = pct, cum_pct = cum_pct)

# 7-2-2. Identify thresholds 
pc_90 <- which(cum_pct >= 90)[1]
pc_01 <- which(diff(pct) > -0.1)[1]

# 7-2-3. Create the Main Elbow Plot
p_main <- ggplot(plot_df, aes(x = PC, y = stdev)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = pc_90, color = "royalblue", linetype = "solid") +
  geom_vline(xintercept = pc_01, color = "red", linetype = "dashed") +
  theme_classic() +
  labs(title = "Method for selecting the appropriate number of PCs:\n1) 90% variance (Blue line)\n2) <0.1% drop (Red line)",
       y = "Standard Deviation")

# 7-2-4. Create the Cumulative Contribution Subplots
# Subplot 1: 90% Threshold
p1 <- ggplot(plot_df, aes(x = cum_pct, y = stdev)) +
  geom_point(aes(color = PC <= pc_90), size = 1) +
  geom_text(aes(label = PC, color = PC <= pc_90), vjust = -1, size = 2, alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "cyan", "FALSE" = "tomato")) +
  theme_bw() +
  labs(title = paste0("PC = ", pc_90, " <- 90% reached"),
       x = "Cumulative % Contribution", y = "SD (%)") +
  theme(legend.position = "none")

# 7-2-5. Subplot 2: 0.1% Drop Threshold
p2 <- ggplot(plot_df, aes(x = cum_pct, y = stdev)) +
  geom_point(aes(color = PC <= pc_01), size = 1) +
  geom_text(aes(label = PC, color = PC <= pc_01), vjust = -1, size = 2, alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "cyan", "FALSE" = "tomato")) +
  theme_bw() +
  labs(title = paste0("PC = ", pc_01, " <- <0.1% drop"),
       x = "Cumulative % Contribution", y = "SD (%)") +
  theme(legend.position = "none")

p1 + p2

#8. Clustering the cell
data <- FindNeighbors(data, dims = 1:10,reduction='pca')
data <- FindClusters(data, resolution = 0.1)

#9. Non-linear dimensionality reduction(UMAP)
data <- RunUMAP(data, dims = 1:10)

#9-1. Examine our clusters
DimPlot(data,reduction="umap",label = TRUE)
FeaturePlot(data, features = c("ADIPOQ"),label=TRUE) #Adipocyte cell marker
FeaturePlot(data, features = c("PTPRC"),label=TRUE) #Immune cell marker
FeaturePlot(data, features = c("PROX1"),label=TRUE) #Lymphatic endothelial cell marker
FeaturePlot(data, features = c("STEAP4"),label=TRUE) #Pericyte cell marker
FeaturePlot(data,features = c("JAM2"),label=TRUE) #Endothelial cell marker
FeaturePlot(data,features = c("DCLK1"),label=TRUE) #Adipose stem and progenitor cell marker

VlnPlot(data, features=c('ADIPOQ','PTPRC','PROX1','STEAP4','JAM2','DCLK1'), group.by='seurat_clusters', pt.size=0) & theme(aspect.ratio=1)
                                                                                  
#10. Annotation
#10-1. Finding differentially expressed genes
data <- JoinLayers(data)

cluster2_markers <- FindMarkers(data, ident.1 = 2)
cluster2_markers %>%
  dplyr::filter(avg_log2FC > 1)

#10-2. Annotation
new.cluster.ids <- c("unannotated1", "ASPC", "Adipose cells", "Immunce cells", "Endothelial cells", "Immune cells",
                     "LEC", "Pericyte", "unannotated2")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
data$cell_type <- Idents(data)

#10-3. Check newly annotated clusters
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#11. Differential abundance analysis between 'normal' vs 'obese'
#11-1. Label cells with BMI <= 25 'normal' and cells with BMI > 25 'obese'
data$condition <- ifelse(data$bmi <=25,'normal','obese')

#11-2. Create a frequency table of cell type vs Condition
counts <- as.data.frame(table(data$cell_type, data$condition))
colnames(counts) <- c("Cell_type", "Condition", "Freq")

#11-3. Calculate percentages for each condition
counts <- counts %>%
  group_by(Condition) %>%
  mutate(Percentage = (Freq / sum(Freq)) * 100)

#11-4. Draw the Stacked Barplot
ggplot(counts, aes(x = Condition, y = Percentage, fill = Cell_type)) +
  geom_bar(stat = "identity", color = "white") +
  theme_minimal() +
  labs(title = "Cell Composition by BMI Condition",
       x = "Condition (BMI Status)",
       y = "Percentage of Total Cells (%)",
       fill = "Cell Type") +
  scale_fill_brewer(palette = "Set3") + # Change palette as needed
  theme(axis.text = element_text(size = 12),
        legend.title = element_text(face = "bold"))

#12. Differential expression analysis(gene expression difference in adipose cells of normal vs obese)
#12-1. Run differential expression analysis
adipocyte <- subset(data,subset=cell_type=='Adipose cells')
Idents(adipocyte) <- "condition"
adipocyte_de <- FindMarkers(adipocyte, 
                                    ident.1 = "obese", 
                                    ident.2 = "normal")
#12-3. Visualization via volcano plot

EnhancedVolcano(adipocyte_de,
                lab = rownames(adipocyte_de),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                drawConnectors = TRUE)

#13. Pathway analysis
#13-1. Create a ranking metric
adipocyte_deg$ranking_metric <- sign(adipocyte_deg$avg_log2FC) * (-log10(adipocyte_deg$p_val_adj + 1e-300))

#13-2. Create the named vector for clusterProfiler
gene_list <- adipocyte_deg$ranking_metric
names(gene_list) <- rownames(adipocyte_deg)

#13-3. Sort the list in descending order (CRITICAL for GSEA)
gene_list <- sort(gene_list, decreasing = TRUE)

#13-4. Download Hallmark gene sets for Humans
h_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

#13-5. Run GSEA
gsea_res <- GSEA(gene_list, TERM2GENE = h_df, pvalueCutoff = 0.05)

#13-6. Visualize result
dotplot(gsea_res, showCategory = 10, split = ".sign") + 
  facet_grid(.~.sign) +
  theme_minimal() +
  labs(title = "GSEA: Obese vs. Normal Adipocytes")


