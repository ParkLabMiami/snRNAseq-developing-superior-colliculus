

# DecontX test ------------------------------------------------------------

# BiocManager::install('celda')
library(patchwork)
library(topGO)
library(org.Mm.eg.db)
library(scran)
library(scuttle)
library(DESeq2)
library(DropletUtils)
library(Seurat)
library(ggplot2)
library(dplyr)
library(enrichR)
library(simplifyEnrichment)
library(celda)

samplesheet <- read.csv(file = 'data/samplesheet.csv')
sample_names <- samplesheet$SampleName
sc <- readRDS(file = 'data/sc.rds')
raw.counts <- vector(mode = 'list', length = length(samplesheet$raw_feature_bc_matrix_path))
names(raw.counts) <- samplesheet$SampleName
for (i in seq_along(samplesheet$raw_feature_bc_matrix_path)) {
  raw.counts[[i]] <- Read10X_h5(filename = samplesheet$raw_feature_bc_matrix_path[i])
}

# P4 test case ------------------------------------------------------------

# object setup
p4 = sc[, sc$time == 'P4']
p4.counts = SingleCellExperiment(assays = list(counts = p4@assays$RNA@counts))
colData(p4.counts) = cbind(colData(p4.counts), p4@meta.data)
p4.raw = SingleCellExperiment(assays = list(counts = raw.counts$P4))
colnames(p4.raw) = paste(colnames(p4.raw), 'P4', sep = '_')
all(colnames(p4.counts) %in% colnames(p4.raw))

# Run
p4.counts = decontX(x = p4.counts, background = p4.raw, varGenes = 2500, z = p4.counts$celltype)
decontXcounts(p4.counts) = round(decontXcounts(p4.counts))
p4.counts = logNormCounts(x = p4.counts, assay.type = 'counts', name = 'logcounts')
p4.counts = logNormCounts(x = p4.counts, assay.type = 'decontXcounts', name = 'logdecontXcounts')
assays(p4.counts)
geneFit.logcounts = modelGeneVar(p4.counts, assay.type = 'logcounts', density.weights = FALSE)
geneFit.logdecontXcounts = modelGeneVar(p4.counts, assay.type = 'logdecontXcounts', density.weights = FALSE)

# Verify mean-variance relationship in rounded, DecontX output numeric matrix
{
par(mfrow = c(1,2))
plot(x = geneFit.logcounts@metadata$mean, 
     y = geneFit.logcounts@metadata$var,
     xlab = 'Mean of log-expression',
     ylab = 'Variance of log-expression',
     main = 'log-counts')
curve(geneFit.logcounts@metadata$trend(x), col = 'dodgerblue', add = TRUE, lwd = 2)
plot(x = geneFit.logdecontXcounts@metadata$mean, 
     y = geneFit.logdecontXcounts@metadata$var,
     xlab = 'Mean of log-expression',
     ylab = 'Variance of log-expression',
     main = 'log-decontXcounts')
curve(geneFit.logdecontXcounts@metadata$trend(x), col = 'dodgerblue', add = TRUE, lwd = 2)
par(mfrow = c(1,1))
}

## Conclusion: mean-variance relationship is preserved after DecontX

# Aggregate results
p4@meta.data$decontX_contamination = p4.counts$decontX_contamination
p4@meta.data$decontX_clusters = p4.counts$decontX_clusters
p4[['decontXcounts']] = CreateAssayObject(counts = round(decontXcounts(p4.counts)))

# standard analysis
DefaultAssay(p4) = 'decontXcounts'
p4 = p4 %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:15) %>% 
  RunUMAP(dims = 1:15)
p.count = DimPlot(p4, group.by = 'celltype', pt.size = 2, shuffle = TRUE, label.size = 5)

DefaultAssay(p4) = 'RNA'
p4 = p4 %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:15) %>% 
  RunUMAP(dims = 1:15)
p.decontX = DimPlot(p4, group.by = 'celltype', pt.size = 2, shuffle = TRUE, label.size = 5)

# compare UMAPs
p.count | p.decontX

rm(p4, p4.raw, p4.counts, geneFit.logcounts, geneFit.logdecontXcounts)


# Run DecontX all samples -------------------------------------------------

# object setup
sc.counts = SingleCellExperiment(assays = list(counts = sc@assays$RNA@counts))
colData(sc.counts) = cbind(colData(sc.counts), sc@meta.data)
for (i in 1:length(raw.counts)) {
  colnames(raw.counts[[i]]) = paste(colnames(raw.counts[[i]]), names(raw.counts)[i], sep = '_')
}
sc.raw = Reduce(f = cbind, x = raw.counts)
if (!all(colnames(sc.counts) %in% colnames(sc.raw))) {
  stop('unmatched cell barcodes')
}
sc.raw = SingleCellExperiment(assays = list(counts = sc.raw))
sc.raw$time = sapply(X = strsplit(x = colnames(sc.raw), split = '_'), `[[`, 2)
table(sc.raw$time)

# Run
sc.counts = decontX(
  x = sc.counts, 
  batch = sc.counts$time,
  background = sc.raw,
  bgBatch = sc.raw$time,
  z = sc.counts$celltype
)

# Aggregate results
sc@meta.data$decontX_contamination = sc.counts$decontX_contamination
sc@meta.data$decontX_clusters = sc.counts$decontX_clusters
sc[['decontXcounts']] = CreateAssayObject(counts = round(decontXcounts(sc.counts)))

# standard analysis
DefaultAssay(sc) = 'decontXcounts'
sc = sc %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:15) %>% 
  RunUMAP(dims = 1:15)
p.celltype = DimPlot(sc, group.by = 'celltype', pt.size = 1, shuffle = TRUE, label.size = 5, label = TRUE)
p.time = DimPlot(sc, group.by = 'time', pt.size = 1, shuffle = TRUE)
p.celltype | p.time


Idents(sc) = 'time'
de.rna = FindAllMarkers(
  object = sc[, sc$celltype == 'Excitatory Neuron'],
  assay = 'RNA', 
  slot = 'data',
  only.pos = TRUE
)
de.decontx = FindAllMarkers(
  object = sc[, sc$celltype == 'Excitatory Neuron'],
  assay = 'decontXcounts', 
  slot = 'data',
  only.pos = TRUE
)
colnames(de.rna)[1:5] = paste(colnames(de.rna)[1:5], 'rna', sep = '.')
colnames(de.decontx)[1:5] = paste(colnames(de.decontx)[1:5], 'decontx', sep = '.')

de = merge(x = de.decontx, y = de.rna, by = c('gene','cluster'), all.x = TRUE)
de %>% 
  mutate(
    logAdjPval.rna = ifelse(
      test = is.infinite(-log10(p_val_adj.rna)),
      yes = 304,
      no = -log10(p_val_adj.rna)
    ),
    logAdjPval.decontx = ifelse(
      test = is.infinite(-log10(p_val_adj.decontx)),
      yes = 304,
      no = -log10(p_val_adj.decontx)
    )
  ) %>% 
  mutate(sig.diff = logAdjPval.decontx - logAdjPval.rna,
         fc.diff = avg_log2FC.decontx - avg_log2FC.rna) %>% 
  ggplot(mapping = aes(x = fc.diff, y = sig.diff)) + 
  geom_point() +
  geom_text(mapping = aes(label = ifelse(sig.diff > 10 | fc.diff > 0.25,
                                         gene,
                                         ''))) +
  facet_wrap(. ~ cluster) + 
  theme_bw()

# Results: For samples where decontX did not successfully remove contaminant genes, DE test results were affected the most. Case: Stmn1 is globally detected in all cell-types and global expression values decrease with time. After DecontX, Stmn1 is "removed" from neuronal cells in samples P4, P8, and P21 but not in E19. DE test comparing E19 and P4 in excitatory neurons shows Cst3 as top gene whose results differ between the two count matrices.

de[de$gene == 'Cst3',]

DefaultAssay(sc) = 'RNA'
p.cst3.rna = FeaturePlot(sc, 'Cst3', order = TRUE, split.by = 'time')
DefaultAssay(sc) = 'decontXcounts'
p.cst3.decontx = FeaturePlot(sc, 'Cst3', order = TRUE, split.by = 'time')
p.cst3.rna / p.cst3.decontx

DefaultAssay(sc) = 'RNA'
p.stmn1.rna = FeaturePlot(sc, 'Stmn1', order = TRUE, split.by = 'time')
DefaultAssay(sc) = 'decontXcounts'
p.stmn1.decontx = FeaturePlot(sc, 'Stmn1', order = TRUE, split.by = 'time')
p.stmn1.rna / p.stmn1.decontx

DefaultAssay(sc) = 'RNA'
p.plp1.rna = FeaturePlot(sc, 'Plp1', split.by = 'time', order = TRUE,)
DefaultAssay(sc) = 'decontXcounts'
p.plp1.decontx = FeaturePlot(sc, 'Plp1', split.by = 'time', order = TRUE,)
p.plp1.rna / p.plp1.decontx

DefaultAssay(sc) = 'RNA'
p.0610040F04Rik.rna = FeaturePlot(sc, '0610040F04Rik', split.by = 'time', order = TRUE,)
DefaultAssay(sc) = 'decontXcounts'
p.0610040F04Rik.decontx = FeaturePlot(sc, '0610040F04Rik', split.by = 'time', order = TRUE,)
p.0610040F04Rik.rna / p.0610040F04Rik.decontx



# Summary of findings:
# DecontX can successfully remove ambient RNAs from individual samples. However, with the experimental current design of the superior colliculus snRNAseq data, there are other non-ambient RNA sources of batch effects that prevents a time-course analysis without the use of batch correction/data integration methods.

# FindMarkers methods test ------------------------------------------------

DefaultAssay(sc) = 'RNA'
Idents(sc) = 'time'
de = FindAllMarkers(
  object = sc[, sc$celltype == 'Excitatory Neuron'],
  assay = 'RNA',
  slot = 'counts',
  logfc.threshold = 0.25,
  only.pos = TRUE,
  test.use = 'wilcox'
)
View(de)
tmp = DotPlot(sc, assay = 'RNA', features = c('Nnat',''))


