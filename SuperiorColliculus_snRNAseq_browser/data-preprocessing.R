
# Data processing script to generate Superior Colliculus snRNAseq web browser

# This script prepares gene expression data and associated metadata. 


# Libraries ---------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(plotly)
library(dplyr)
library(tibble)
library(tidyr)
library(HDF5Array)


# Directory setup ---------------------------------------------------------

# app.dir <- 'SuperiorColliculus_snRNAseq_browser/'
# dir.create(path = app.dir)

setwd('SuperiorColliculus_snRNAseq_browser/')

sc <- readRDS(file = '../data/sc.rds')
neuron <- readRDS(file = '../data/neuron.rds')


# Gather metadata ---------------------------------------------------------

whichNeuronMeta = c('subtype','celltype_adjusted','SingleR_Transseq', 'SingleR_Vectorseq', 'UMAP_1', 'UMAP_2')

# Neuron metadata
neuron.metadata <- FetchData(
  object = neuron,
  vars = whichNeuronMeta
)
colnames(neuron.metadata) <- gsub(pattern = 'UMAP_', replacement = 'neuron_UMAP_', x = colnames(neuron.metadata))

whichSCMeta = c('celltype','celltype_subtype','time','integrated_snn_res.0.8','nCount_RNA','nFeature_RNA','percent_mt','percent_rp','UMAP_1','UMAP_2')

# Full superior colliculus metadata
metadata <- FetchData(
  object = sc,
  vars = whichSCMeta
)
colnames(metadata) <- gsub(pattern = 'UMAP_', replacement = 'sc_UMAP_', x = colnames(metadata))
rownames(metadata) <- gsub(pattern = '\\.', replacement = '_', x = rownames(metadata))
sc.column.names <- gsub(pattern = '\\.', replacement = '_', colnames(sc))

# Extract subtypes identified from neuron analysis and transfer to full metadata
metadata$subtype <- as.character(metadata$celltype_subtype)
metadata$subtype[match(rownames(neuron.metadata), rownames(metadata))] <- as.character(neuron.metadata$subtype)

# Extract adjusted celltypes identified from neuron analysis and transfer to full metadata
metadata$celltype <- as.character(metadata$celltype)
metadata$celltype[match(rownames(neuron.metadata), rownames(metadata))] <- as.character(neuron.metadata$celltype_adjusted)

# Extract SingleR predictions from neuron analysis
metadata$SingleR_Transseq = NA
metadata$SingleR_Transseq[match(rownames(neuron.metadata), rownames(metadata))] = as.character(neuron.metadata$SingleR_Transseq)
metadata$SingleR_Vectorseq = NA
metadata$SingleR_Vectorseq[match(rownames(neuron.metadata), rownames(metadata))] = as.character(neuron.metadata$SingleR_Vectorseq)

# Combine subtype and adjusted celltype information
metadata <- base::merge(
  x = metadata, 
  y = neuron.metadata[grepl('neuron', colnames(neuron.metadata))],
  by = 'row.names', 
  all.x = TRUE
) %>% 
  column_to_rownames(var = 'Row.names')
metadata$celltype_subtype <- NULL
sc.column.names <- gsub(pattern = '\\.', replacement = '_', colnames(sc))
metadata <- metadata[sc.column.names,]

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


# Data setup for browser UI -----------------------------------------------

# dataset key-value pairs for subsetting data. Values should match prefixes from above e.g. 'neuron' for 'neuron_UMAP_1'
dataset_dict <- c(
  SuperiorColliculus = 'sc',
  Neurons = 'neuron'
)

library_size = metadata$nCount_RNA
feature_size = metadata$nFeature_RNA

# Categorical variables by which cells can be grouped.
categoricalVars <- colnames(metadata)[sapply(metadata, function(x) !is.numeric(x))]

# All possible genes to query
all_features <- rownames(sc[['RNA']]@data)

# Dimplot label coordinates
dimReducVars <- grep(pattern = 'PC|tSNE|UMAP', colnames(metadata), value = TRUE)
label_coords <- metadata[c(categoricalVars, dimReducVars)] %>%
  pivot_longer(cols = !contains(match = c('PC','UMAP','tSNE')), 
               names_to = 'variables', 
               values_to = 'values') %>%
  group_by(variables, values) %>% 
  summarise(across(.cols = all_of(dimReducVars), 
                   .fns = function(x) median(x, na.rm = TRUE))) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'variables', values_from = 'values') %>% 
  as.data.frame()


# Website details ---------------------------------------------------------

window_title = "snRNAseq of the developing mouse superior colliculus"
title_link_text = "by Park lab"
title_link_url = "https://www.parklabmiami.com"


# Save data ---------------------------------------------------------------

# Write raw counts matrix as h5 file (for 66k cells and 33k genes, ~ 300Mb)
#
# NOTE: By taking the transpose, "barcodes" and "genes" are swapped in the .h5 file. This is so that columns correspond to features/genes and rows correspond to observations/cells. 
counts <- Matrix::t(sc[['RNA']]@counts)

# Note: If `sc.h5` already exists, this throws an error. For purpose of recompiling the web portal, delete `sc.h5` before running this.
HDF5Array::writeTENxMatrix(
  x = counts,
  filepath = 'sc.h5', 
  verbose = TRUE
)

save(metadata, categoricalVars, dataset_dict, dimReducVars, all_features, label_coords, window_title, title_link_text, title_link_url, library_size, feature_size, file = 'appdata.RData')
