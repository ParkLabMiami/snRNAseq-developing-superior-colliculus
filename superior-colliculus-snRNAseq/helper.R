
#### Helper scripts and functions (e.g. plotting) for the website ###


# Env variables -----------------------------------------------------------

# Pre-computed values and expression matrix h5
load(file = 'appdata.RData')
f <- h5file(filename = 'sc.h5', mode = 'r')

# Variables for storing expression values
# 2022-08-22: under-the-hood data structures were change to avoid retrieving epression data from non-matrices. 
t_init <- Sys.time()
log_x_sc <- matrix(NA, nrow = 1, ncol = length(rownames(metadata)))
colnames(log_x_sc) <- rownames(metadata)
last_query <- 1
log_times <- .POSIXct(xx = rep(t_init, 1000))
names(log_times) <- character(length(log_times))

# Colors
time.cols <- c("#dd4539","#d9b123","#68c545","#514fd1")
celltype.cols <- RColorBrewer::brewer.pal(n = 12, name = 'Paired')
celltype.cols[11] <- 'gold'
neuron.cols <- c("#c2444c","#e14327","#be6231","#de8d26","#debc22","#b99f3e","#a6c336","#6f9a3e","#63d135","#59c251","#4db873","#36dbbc","#7b85dc","#6074eb","#554ea8","#6934c4","#9253ea","#412070","#542298","#ac69d3","#ce4ce2","#d980cc","#cc4bb5","#df35bd","#92337c","#d54681")
combined.subtype.cols = c("#bea3de","#66db40","#7048dc","#a7c738","#bc4bdf","#66cb5b","#bc3eaf","#54d592","#e143a4","#488b33","#6a75e2","#d9b83a","#7a4eb0","#b4c072","#d980d2","#83be84","#e14275","#6bcab4","#e44630","#71c8e7","#e38836","#5c86cc","#b04c26","#4699aa","#b63c4b","#3c805e","#a34877","#767427","#775e90","#a7782d","#366e8e","#e17e74","#637745","#da95ad","#daae7b","#955d47")
SingleR_Transseq.cols = c("#c86f3c","#6588cd","#aeae43","#a863c3","#56b980","#d5447f","#c06171","#5e7d36")
SingleR_Vectorseq.cols = c("#bea3de","#66db40","#7048dc","#a7c738","#bc4bdf","#66cb5b","#bc3eaf","#e143a4","#488b33","#6a75e2","#d9b83a","#7a4eb0","#b4c072","#d980d2","#83be84","#e14275","#6bcab4","#e44630","#71c8e7","#e38836","#5c86cc","#b04c26","#4699aa","#b63c4b","#3c805e","#a34877","#767427","#775e90","#a7782d","#366e8e","#e17e74","#637745","#da95ad","#daae7b","#955d47")

my_colors <- list(
  sc = list(
    celltype = celltype.cols,
    time = time.cols,
    subtype = combined.subtype.cols,
    SingleR_Transseq = SingleR_Transseq.cols,
    SingleR_Vectorseq = SingleR_Vectorseq.cols
  ),
  neuron = list(
    celltype = celltype.cols,
    time = time.cols,
    subtype = neuron.cols,
    SingleR_Transseq = SingleR_Transseq.cols,
    SingleR_Vectorseq = SingleR_Vectorseq.cols
  )
)

# Helper functions --------------------------------------------------------

readFeatureValue <- function(feature) {
  # for loop makes it so function can take as input integer(length = 2)
  for ( feat in feature ) {
    # if feature is not already in expression vector
    if (!feat %in% rownames(log_x_sc)) {
      # get expression vector slot with earliest retrieval time
      # returns sequentially since all have t_init at first
      last_query <<- which.min(log_times)
      
      # extract non-zero values `g.i`, column `c`, and row `i` values from h5 connection. indexing in h5 file is 0-based. 
      gene_index <- which(f[[names(f)[1]]][['barcodes']][1:length(all_features)] == feat)
      i <- c(f[[names(f)[1]]][['indptr']][gene_index:(gene_index + 1)]) + 1
      nonzero_vals <- f[[names(f)[1]]][['data']][i[1]:(i[2]-1)]
      c <- f[[names(f)[1]]][['indices']][i[1]:(i[2]-1)] + 1
      tmp <- rep.int(0, times = nrow(metadata))
      tmp[c] <- nonzero_vals
      
      # Normalize by library size
      tmp <- log1p(tmp/library_size * 10000)
      log_x_sc <<- rbind(log_x_sc, tmp)
      rownames(log_x_sc)[nrow(log_x_sc)] <<- feat
      log_times[last_query] <<- Sys.time()
      names(log_times)[last_query] <<- feat
    }
  }
}

getGroupColorPalette <- function(dataset_value, groupby) {
  return(my_colors[[dataset_value]][[groupby]])
}

drawDimPlot <- function(
  dataset_value, 
  groupby,
  reduction = 'UMAP',
  draw_labels
) {
  dims <- paste(reduction, 1:2, sep = '_')
  plot_these <- c(paste(dataset_value, dims, sep = '_'), groupby)
  umap_df <- metadata[plot_these]
  colnames(umap_df) <- c('dim1', 'dim2', 'groupby')
  umap_df <- umap_df[!is.na(umap_df$dim1),]
  if (is.factor(umap_df$groupby)) {
    umap_df$groupby <- droplevels(umap_df$groupby)
  }
  umap_df <- umap_df[sample(1:nrow(umap_df), nrow(umap_df)),]
  label_df <- label_coords[plot_these]
  colnames(label_df) <- c('dim1', 'dim2', 'label')
  label_df <- label_df[!is.na(label_df$dim1) & !is.na(label_df$label),]
  label_df$label <- factor(label_df$label, levels = levels(umap_df$groupby))
  
  cols <- getGroupColorPalette(dataset_value = dataset_value, groupby = groupby)
  
  p <- ggplot(data = umap_df) +
    geom_point(mapping = aes(x = dim1, y = dim2, color = groupby),
               size = 1.5) +
    scale_color_manual(values = cols) + 
    theme_bw() +
    xlab(paste(reduction, 1, sep = '_')) +
    ylab(paste(reduction, 2, sep = '_')) +
    theme(axis.text = element_blank(),
          axis.title = element_text(size = 14),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    guides(color = guide_legend(title = groupby,
                                override.aes = list(size = 6)))
  if (draw_labels) {
    p <- p +
      geom_text(
        data = label_df,
        mapping = aes(x = dim1, y = dim2, label = label),
        size = 5
      )
  }
  return(p)
}


drawFeaturePlot <- function(
  dataset_value, 
  feature,
  reduction = 'UMAP'
) {
  dims <- paste(reduction, 1:2, sep = '_')
  umap_df <- metadata[paste(dataset_value, dims, sep = '_')]
  readFeatureValue(feature = feature)
  umap_df$feature <- log_x_sc[feature,]
  colnames(umap_df) <- c('dim1','dim2','feature')
  umap_df <- umap_df[order(umap_df$feature),]
  
  p <- ggplot(data = umap_df) +
    geom_point(mapping = aes(x = dim1, y = dim2, color = feature),
               size = 1.5) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(n = 9, name = 'YlOrRd')) +
    theme_bw() +
    xlab(paste(reduction, 1, sep = '_')) +
    ylab(paste(reduction, 2, sep = '_')) +
    theme(axis.text = element_blank(),
          axis.title = element_text(size = 14),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    guides(color = guide_colorbar(frame.colour = 'black',
                                  ticks.colour = 'black',
                                  title = feature))
  return(p)
}


drawSplitFeaturePlot <- function(
  dataset_value, 
  feature,
  reduction = 'UMAP'
) {
  dims <- paste(reduction, 1:2, sep = '_')
  umap_df <- metadata[c(paste(dataset_value, dims, sep = '_'), 'time')]
  readFeatureValue(feature = feature)
  umap_df$feature <- log_x_sc[feature,]
  umap_df <- umap_df[order(umap_df$feature),]
  colnames(umap_df) <- c('dim1','dim2','split', 'feature')
  
  p <- ggplot(data = umap_df) +
    geom_point(mapping = aes(x = dim1, y = dim2, color = feature),
               size = 1.5) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(n = 9, name = 'YlOrRd')) +
    facet_wrap(. ~ split, ncol = 4) + 
    theme_bw() +
    xlab(paste(reduction, 1, sep = '_')) +
    ylab(paste(reduction, 2, sep = '_')) +
    theme(axis.text = element_blank(),
          axis.title = element_text(size = 14),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    guides(color = guide_colorbar(frame.colour = 'black',
                                  ticks.colour = 'black',
                                  title = feature))
  return(p)
}


drawDotPlot <- function(
  dataset_value,
  feature,
  groupby
) {
  plot_these <- c(paste0(dataset_value, '_UMAP_1'), groupby)
  data_df <- metadata[plot_these]
  readFeatureValue(feature = feature)
  if (length(feature) >= 2) {
    data_df <- cbind(data_df, t(log_x_sc[feature,]))
  } else {
    data_df <- cbind(data_df, log_x_sc[feature,])
    colnames(data_df)[3] <- feature
  }
  colnames(data_df) <- c('dataset', 'groupby', feature)
  data_df <- data_df[!is.na(data_df$groupby) & 
                       !is.na(data_df$dataset),]
  pct <- data_df %>% 
    group_by(groupby) %>% 
    summarise(across(all_of(feature), .fns = function(a) mean(a > 0))) %>%
    ungroup() %>% 
    reshape2::melt(id.vars = c('groupby'))
  colnames(pct) <- c('groupby','feature','pct')
    
  avg <- data_df %>% 
    mutate(across(all_of(feature), .fns = scale)) %>% 
    group_by(groupby) %>% 
    summarise(across(all_of(feature), .fns = function(a) mean(a))) %>%
    ungroup() %>% 
    reshape2::melt(id.vars = c('groupby'))
  colnames(avg) <- c('groupby','feature','avg')
  
  lims <- quantile(avg$avg, probs = c(0.05, 0.95))
  avg$avg[avg$avg < lims[1]] <- lims[1]
  avg$avg[avg$avg > lims[2]] <- lims[2]
  
  p <- merge(pct, avg) %>% 
    ggplot(mapping = aes(x = groupby, y = feature)) +
    geom_point(mapping = aes(size = pct, fill = avg), pch = 21) +
    # scale_fill_viridis_c(option = 'A') +
    scale_fill_gradient(low = 'grey90', high = 'blue') +
    scale_size(range = c(0, 6), limits = c(0, 1)) +
    theme_bw() +
    xlab(label = groupby) +
    theme(axis.text.x = element_text(size = 16, angle = 65, hjust = 1),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                                      size = 16),
          legend.title.align = 1,
          legend.text = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title.x = element_blank(),
          legend.spacing.y = unit(1, "cm")) +
    guides(fill = guide_colorbar(title = 'Average scaled\nexpression',
                                 title.position = 'left',
                                 frame.colour = 'black',
                                 title.hjust = 0.5,
                                 title.vjust = 0.5,
                                 ticks = TRUE,
                                 ticks.colour = 'black'),
           size = guide_legend(title = 'Percent detected',
                               title.position = 'left',
                               override.aes = list(color = 'black',
                                                   fill = 'black')))
  return(p)
}

drawSplitDotPlot <- function(
  dataset_value,
  feature,
  groupby,
  splitby = 'time',
  ncol = 2
) {
  plot_these <- c(paste0(dataset_value, '_UMAP_1'), groupby, splitby)
  data_df <- metadata[plot_these]
  readFeatureValue(feature = feature)
  if (length(feature) >= 2) {
    data_df <- cbind(data_df, t(log_x_sc[feature,]))
  } else {
    data_df <- cbind(data_df, log_x_sc[feature,])
    colnames(data_df)[4] <- feature
  }
  colnames(data_df)[1:3] <- c('dataset', 'groupby', 'splitby')
  data_df <- data_df[!is.na(data_df$groupby) & 
                       !is.na(data_df$dataset),]
  pct <- data_df %>% 
    group_by(groupby, splitby) %>% 
    summarise(across(all_of(feature), .fns = function(a) mean(a > 0))) %>%
    reshape(varying = feature,
            timevar = 'feature',
            idvar = c('groupby','splitby'),
            v.names = 'pct',
            times = feature,
            direction = 'long')
    
  avg <- data_df %>% 
    mutate(across(all_of(feature), .fns = scale)) %>% 
    group_by(groupby, splitby) %>% 
    summarise(across(all_of(feature), .fns = function(a) mean(a))) %>%
    reshape(varying = feature,
            timevar = 'feature',
            idvar = c('groupby','splitby'),
            v.names = 'avg',
            times = feature,
            direction = 'long')
  
  lims <- quantile(avg$avg, probs = c(0.05, 0.95))
  avg$avg[avg$avg < lims[1]] <- lims[1]
  avg$avg[avg$avg > lims[2]] <- lims[2]
  
  p <- merge(pct, avg) %>% 
    ggplot(mapping = aes(x = groupby, y = feature)) +
    geom_point(mapping = aes(size = pct, fill = avg), pch = 21) +
    facet_wrap(. ~ splitby, ncol = ncol) +
    # scale_fill_viridis_c(option = 'A') +
    scale_fill_gradient(low = 'grey90', high = 'blue') +
    scale_size(range = c(0, 6), limits = c(0, 1)) +
    theme_bw() +
    xlab(label = groupby) +
    theme(axis.text.x = element_text(size = 16, angle = 65, hjust = 1),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                                      size = 16),
          legend.title.align = 1,
          legend.text = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_blank(),
          legend.spacing.y = unit(1, "cm"),
          strip.text = element_text(size = 16)) +
    guides(fill = guide_colorbar(title = 'Normalized expression',
                                 title.position = 'left',
                                 frame.colour = 'black',
                                 title.hjust = 0.5,
                                 title.vjust = 0.5,
                                 ticks = TRUE,
                                 ticks.colour = 'black'),
           size = guide_legend(title = 'Percent detected',
                               title.position = 'left',
                               override.aes = list(color = 'black',
                                                   fill = 'black')))
  return(p)
}


drawFeatureViolinPlot <- function(
  dataset_value,
  feature,
  groupby
) {
  plot_these <- c(paste0(dataset_value, '_UMAP_1'), groupby)
  data_df <- metadata[plot_these]
  readFeatureValue(feature = feature)
  data_df$feature <- log_x_sc[feature,]
  colnames(data_df) <- c('dataset', 'groupby', 'feature')
  data_df <- data_df[!is.na(data_df$groupby) & 
                       !is.na(data_df$feature) &
                       !is.na(data_df$dataset),]
  data_df$groupby <- factor(
    x = data_df$groupby,
    levels = rev(
      x = levels(droplevels(metadata[[groupby]]))
    )
  )
  cols <- getGroupColorPalette(dataset_value = dataset_value, groupby = groupby)
  cols <- rev(cols)
  p <- ggplot(data = data_df) +
    geom_violin(mapping = aes(x = groupby, y = feature, fill = groupby),
                scale = 'width') +
    # geom_beeswarm(mapping = aes(x = groupby, y = feature)) +
    geom_quasirandom(mapping = aes(x = groupby, y = feature), bandwidth = 0.5, size = 0.5, alpha = 0.3, varwidth = TRUE, width = 0.4) +
    scale_fill_manual(values = cols) +
    xlab(label = groupby) +
    ylab(label = feature) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(size = 18, angle = 65, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 14))
  return(p)
}
