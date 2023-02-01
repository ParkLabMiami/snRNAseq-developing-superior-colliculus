

### GEO data prepatation of counts matrix + metadata

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

outs = 'data/GEO-submission/'
dir.create(path = outs)

e19 = Read10X_h5(filename = 'data/raw_feature_bc_matrix/E19_raw_feature_bc_matrix.h5', use.names = FALSE)
p4 = Read10X_h5(filename = 'data/raw_feature_bc_matrix/P4_raw_feature_bc_matrix.h5', use.names = FALSE)
p8 = Read10X_h5(filename = 'data/raw_feature_bc_matrix/P8_raw_feature_bc_matrix.h5', use.names = FALSE)
p21 = Read10X_h5(filename = 'data/raw_feature_bc_matrix/P21_raw_feature_bc_matrix.h5', use.names = FALSE)
sc = readRDS(file = 'data/sc.rds')
neuron = readRDS(file = 'data/neuron.rds')


# Most of this is similar to the preprocessing for the web portal

whichNeuronMeta = c('subtype','celltype_adjusted','SingleR_Transseq', 'SingleR_Vectorseq', 'UMAP_1', 'UMAP_2')

# Neuron metadata
neuron.metadata = FetchData(
  object = neuron,
  vars = whichNeuronMeta
)
colnames(neuron.metadata) = gsub(pattern = 'UMAP_', replacement = 'neuron_UMAP_', x = colnames(neuron.metadata))

whichSCMeta = c('celltype','celltype_subtype','time','integrated_snn_res.0.8','nCount_RNA','nFeature_RNA','percent_mt','percent_rp','UMAP_1','UMAP_2')

# Full superior colliculus metadata
metadata = FetchData(
  object = sc,
  vars = whichSCMeta
)
colnames(metadata) = gsub(pattern = 'UMAP_', replacement = 'sc_UMAP_', x = colnames(metadata))
rownames(metadata) = gsub(pattern = '\\.', replacement = '_', x = rownames(metadata))
sc.column.names = gsub(pattern = '\\.', replacement = '_', colnames(sc))

# Extract subtypes identified from neuron analysis and transfer to full metadata
metadata$subtype = as.character(metadata$celltype_subtype)
metadata$subtype[match(rownames(neuron.metadata), rownames(metadata))] = as.character(neuron.metadata$subtype)

# Extract adjusted celltypes identified from neuron analysis and transfer to full metadata
metadata$celltype = as.character(metadata$celltype)
metadata$celltype[match(rownames(neuron.metadata), rownames(metadata))] = as.character(neuron.metadata$celltype_adjusted)

# Extract SingleR predictions from neuron analysis
metadata$SingleR_Transseq = NA
metadata$SingleR_Transseq[match(rownames(neuron.metadata), rownames(metadata))] = as.character(neuron.metadata$SingleR_Transseq)
metadata$SingleR_Vectorseq = NA
metadata$SingleR_Vectorseq[match(rownames(neuron.metadata), rownames(metadata))] = as.character(neuron.metadata$SingleR_Vectorseq)

# Combine subtype and adjusted celltype information
metadata = base::merge(
  x = metadata, 
  y = neuron.metadata[grepl('neuron', colnames(neuron.metadata))],
  by = 'row.names', 
  all.x = TRUE
) %>% 
  tibble::column_to_rownames(var = 'Row.names')
metadata$celltype_subtype = NULL
sc.column.names = gsub(pattern = '\\.', replacement = '_', colnames(sc))
metadata = metadata[sc.column.names,]

# Reintroduce factor levels
metadata$celltype = factor(
  x = metadata$celltype,
  levels = levels(sc$celltype)
)

metadata$subtype = factor(
  x = metadata$subtype,
  levels = c(
    levels(neuron$subtype), 
    levels(sc$celltype_subtype)[!grepl('Neuron', levels(sc$celltype_subtype), ignore.case = TRUE)]
  )
)

metadata$time = factor(
  x = metadata$time,
  levels = levels(sc$time)
)

metadata$SingleR_Transseq = factor(
  x = metadata$SingleR_Transseq,
  levels = levels(neuron$SingleR_Transseq)
)

metadata$SingleR_Vectorseq = factor(
  x = metadata$SingleR_Vectorseq,
  levels = levels(neuron$SingleR_Vectorseq)
)

# count matrix stuff
counts = sc[['RNA']]@counts
genes = data.frame(
  'Ensembl_ID' = rownames(e19),
  'GeneName' = rownames(counts)
)
barcodes = colnames(counts)

# write out
writeMM(obj = counts, file = paste0(outs, 'sc_mat.mtx'))
write.table(x = genes, file = paste0(outs, 'genes.tsv'), sep = '\t', row.names = FALSE, col.names = FALSE)
write(x = barcodes, file = paste0(outs, 'barcodes.tsv'))
write.table(x = metadata, file = paste0(outs, 'metadata.tsv'), sep = '\t', row.names = TRUE, col.names = NA)
