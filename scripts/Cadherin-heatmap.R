
library(Seurat)

counts <- readRDS(file = 'deprecated/superior colliculus/Data/SC_20200103.rds')
# counts <- readRDS(file = 'deprecated/superior colliculus/Data/SC.rds')

DimPlot(counts, reduction = 'pca', dims = c(1,2), group.by = 'orig.ident', shuffle = TRUE)
DimPlot(counts, reduction = 'pca', dims = c(1,3), group.by = 'orig.ident', shuffle = TRUE)
DimPlot(counts, reduction = 'pca', dims = c(3,4), group.by = 'orig.ident', shuffle = TRUE)
DimPlot(counts, reduction = 'pca', dims = c(3,5), group.by = 'orig.ident', shuffle = TRUE)
DimPlot(counts, reduction = 'pca', dims = c(3,6), group.by = 'orig.ident', shuffle = TRUE)
ElbowPlot(counts)

cadherins <- grep(pattern = 'Cdh', ignore.case = TRUE, x = rownames(counts@assays$RNA), value = TRUE)

cadmat <- DietSeurat(counts, features = cadherins, assays = 'RNA')
cadmat <- NormalizeData(cadmat)
cadmat <- ScaleData(cadmat, features = rownames(cadmat))
cadmat <- RunPCA(cadmat, features = rownames(cadmat))


cadmat <- FindNeighbors()
cadmat <- RunTSNE(cadmat, )

d <- dist(x = cadmat@assays$RNA@scale.data)
c <- hclust(d = d, method = 'ward.D2')
ComplexHeatmap::Heatmap(
  matrix = cadmat@assays$RNA@scale.data,
  show_column_names = FALSE,
  show_row_names = FALSE,
  # cluster_rows = c,
  use_raster = TRUE
)

DefaultAssay(counts) <- 'RNA'
VlnPlot(counts, cadherins[1:3])

neurons <- counts[,counts$CellType %in% c('Excitatory','Inhibitory')]
DefaultAssay(neurons) <- 'RNA'
