

results_out = 'notes/ambientContribMaximum-patterns/'
dir.create(results_out)

# Aggregate count setup ---------------------------------------------------

sc_sce <- SingleCellExperiment(
  assays = list(counts = sc[['RNA']]@counts)
)
colData(sc_sce) <- cbind(colData(sc_sce), sc@meta.data)
sc_summed <- aggregateAcrossCells(
  x = sc_sce,
  ids = paste(sc_sce$celltype, sc_sce$time, sep = '_')
)
sc_summed$ids <- factor(x = sc_summed$ids)
sc_summed$ids <- factor(
  x = sc_summed$ids,
  levels = apply(
    X = expand.grid(levels(sc$time), levels(sc$celltype)), 
    MARGIN = 1, 
    FUN = function(x) paste(x[2],x[1],sep = '_')
  )
)
# 2023-01-19: discovered potential bug whereby timepoints in `ambient` did not correspond to matching celltype_time ids in `sc_summed`. According to ?ambientContribMaximum, these needs to be corresponding. Prior to date, it was not. 
sc_summed = sc_summed[, order(sc_summed$ids)]



# Comparing threshold= values in ambientContribMaximum --------------------

# Note: NA/NaN values indicate values for which expected counts from ambient profile exceeds observed counts in aggregate counts

# p-value threshold = 0.1 (discard fewer genes)
celltypes <- levels(sc$celltype)
max.contrib.1 <- list()
for (c in celltypes) {
  c.pull <- grepl(pattern = c, x = sc_summed$ids)
  summed <- sc_summed[, c.pull]
  max.contrib.1[[c]] <- ambientContribMaximum(
    y = counts(summed),
    ambient = ambient,
    threshold = 0.1,
    mode = 'proportion'
  )
  message('Done with: ', c)
}

# p-value threshold = 0.01 (discard more genes)
celltypes <- levels(sc$celltype)
max.contrib.01 <- list()
for (c in celltypes) {
  c.pull <- grepl(pattern = c, x = sc_summed$ids)
  summed <- sc_summed[, c.pull]
  max.contrib.01[[c]] <- ambientContribMaximum(
    y = counts(summed),
    ambient = ambient,
    threshold = 0.01, 
    mode = 'proportion'
  )
  message('Done with: ', c)
}

# Plot per-gene contrib values for the two thresholds
e.1 = max.contrib.1$`Excitatory Neuron`
colnames(e.1) = paste(colnames(e.1), 'p0.1', sep = '_')
e.01 = max.contrib.01$`Excitatory Neuron`
colnames(e.01) = paste(colnames(e.01), 'p0.01', sep = '_')
e = cbind(e.1, e.01)
head(e)
diffData = as.data.frame(e) %>% 
  tibble::rownames_to_column(var = 'gene') %>%
  tidyr::pivot_longer(cols = -gene) %>% 
  tidyr::separate(col = 'name', into = c('celltype', 'time', 'pval'), sep = '_') %>% 
  tidyr::pivot_wider(names_from = pval, values_from = value) %>% 
  mutate(time = factor(time, levels = levels(sc$time)),
         contribDiff = p0.01 - p0.1)
p.contribution = diffData %>% 
  ggplot(mapping = aes(x = p0.1, y = p0.01)) +
  geom_abline(slope = 1, intercept = 0, lty = 'dashed', color = 'indianred') +
  geom_point(size = 0.5, alpha = 0.7) + 
  facet_grid(. ~ time) +
  labs(y = '`threshold=0.01`', x = '`threshold=0.1`') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 65, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))
p.contribution
ggsave(filename = paste0(results_out, 'ambient-scaling-affects-contribution.tiff'), plot = p.contribution, device = 'tiff', height = 2.5, width = 7, dpi = 320)

p.contributionDiff = diffData %>% 
  ggplot(mapping = aes(x = p0.1, y = contribDiff)) +
  geom_point(size = 0.5, alpha = 0.7) + 
  facet_grid(. ~ time)  +
  labs(y = 'Difference between\n`threshold=0.01`-`threshold=0.1`', x = '`threshold=0.1`') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 65, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))
p.contributionDiff
ggsave(filename = paste0(results_out, 'ambient-scaling-affects-already-high-contributors.tiff'), plot = p.contributionDiff, device = 'tiff', height = 2.5, width = 7, dpi = 320)

# Background: The `threshold=` argument effectively controls how much to scale the ambient RNA profile, which is then compared against the observed gene count. More stringent (numerically smaller)  `threshold=` values indicates greater upwards scaling of the ambient profile before rejecting the joint-null and effectively increases the maximum possible ambient contribution proportion. For a given contribution proportion threshold, decreasing `threshold=` yields greater number of genes above the said proportion threshold.

# Results: Using different `threshold=` values in `ambientContribMaximum` yields a linear curve with some genes plateauing at the upper ends of maximum contribution values. The resulting monotonic curve suggests that per-gene proportions are independent of each other (as expected from documentation `?ambientContribMaximum`). Plotting the difference in contribution proportions between the two thresholds shows that proportions can shift by up to 0.1. However, genes for which the shift was greater than 0.1 are also genes for which ambient contribution proportions are above 0.6, which is already a high proportion.


# Identify suitable contribution threshold --------------------------------

celltypes <- levels(sc$celltype)
max.contrib <- list()
for (c in celltypes) {
  c.pull <- grepl(pattern = c, x = sc_summed$ids)
  summed <- sc_summed[, c.pull]
  max.contrib[[c]] <- ambientContribMaximum(
    y = counts(summed),
    ambient = ambient,
    threshold = 0.1,
    mode = 'proportion'
  )
  message('Done with: ', c)
}

# max.contrib$`Excitatory Neuron`[order(max.contrib$`Excitatory Neuron`[,1], decreasing = TRUE, na.last = TRUE),][500:530,1:2]

# Plot of ranked contributions
p.geneRank = data.frame(
  rank = rank(max.contrib$`Excitatory Neuron`[,1], na.last = TRUE),
  contrib = max.contrib$`Excitatory Neuron`[,1]
) %>% 
  ggplot() + 
  geom_point(mapping = aes(x = rank, y = contrib)) +
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  labs(x = 'Gene rank',
       y = 'Ambient contribution proportion') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, hjust = 1, angle = 65),
        axis.title = element_text(size = 12))
p.geneRank
ggsave(filename = paste0(results_out, 'ambient-contribution-gene-rank.tiff'), plot = p.geneRank, device = 'tiff', height = 3.5, width = 3.5, dpi = 320)

# Plot of gene percent detection rate versus whether gene would have been removed due to higher than threshold contribution value
c = celltypes[1]
pct.exp = list()
for (t in levels(sc$time)) {
  pct.exp[[t]] = sparseMatrixStats::rowMeans2(
    x = (sc@assays$RNA@counts > 0), 
    cols = which(sc$celltype == c & sc$time == t)
  )
}
pct.exp = data.frame(
  gene = rownames(sc),
  as.data.frame(pct.exp)
) %>% 
  tidyr::pivot_longer(cols = -gene) %>% 
  dplyr::rename(pct.exp = value,
                time = name) %>% 
  mutate(time = factor(time, levels = levels(sc$time)))
contrib = max.contrib[[c]] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'gene') %>% 
  tidyr::pivot_longer(cols = -gene) %>% 
  tidyr::separate(col = name, into = c('celltype', 'time'), sep = '_') %>% 
  dplyr::rename(contrib = value) %>% 
  mutate(time = factor(time, levels = levels(sc$time)))
contribData = dplyr::left_join(pct.exp, contrib, by = c('gene','time'))

thresholds = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5)
tmp = list()
for (i in 1:length(thresholds)) {
  tmp[[i]] = contribData %>% 
    mutate(remove = ifelse(contrib > thresholds[i], 'remove', 'keep'),
           threshold = thresholds[i])
}
p.thresholdRange = Reduce(f = rbind, x = tmp) %>%
  filter(!is.na(contrib)) %>% 
  ggplot() + 
  # geom_density(mapping = aes(x = pct.exp, fill = remove), alpha = 0.4) +
  geom_histogram(mapping = aes(x = pct.exp, fill = remove), 
                 bins = 50, 
                 position = 'identity',
                 alpha = 0.4) +
  facet_grid(threshold ~ time) +
  scale_y_continuous(trans = 'log1p',
                     labels = scales::label_number(),
                     breaks = 10^seq(0, 6, 1)) +
  scale_fill_discrete(name = 'Retention') +
  labs(x = 'Percent detection', y = 'Number of genes', 
       title = 'Retention of genes across ambientContribMaximum thresholds') +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 65, hjust = 1),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
p.thresholdRange
ggsave(filename = paste0(results_out, 'contribution-threshold-gene-retention-histogram.tiff'), plot = p.thresholdRange, device = 'tiff', height = 6, width = 6.5, dpi = 320)
  

# Results: As contribution threshold increases, i.e. as genes are allowed to have higher ambient contribution to the observed count, genes at the upper end of percent detection are retained at increasing proportions. Conversely, genes at the lower end of percent detection are also retained at increasing proportions but at a lower rate i.e. highly detected genes are first to be retained as contribution threshold increases. This implies that genes for which ambient contribution is significant usually also have low-mid detection rates.

FeaturePlot(sc, order = TRUE, split.by = 'time', 'Src')


# Conclusions -------------------------------------------------------------

# Here I compared how changed `threshold=` in `ambientContribMaximum` controls the estimated contribution proportion. `threshold=` I also compared how various thresholds on the estimated contribution proportion affect which genes are kept and removed. I recommend the following:
# 1) use `threshold=0.1` in `ambientContribMaximum` for identifying a suitable scale factor for ambient RNA profile. Genes with high contribution proportions are affected the most by various this parameter, which would not impact downstream DEG filtering because of filtering by explicit proportion thresholds.
# 2) Instead, vary the contribution proportion threshold above which genes are removed from DEG results. 

# rom empirically checking various genes with estimated maximum contribution proportions near a threshold of 0.1, I see that, for Excitatory Neurons, some genes are specific to non-Excitatory Neurons. The density distribution of 