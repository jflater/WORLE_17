library(phyloseq)
library(phylosmith)
library(tidyverse)
library(viridis)

phy <- readRDS("worle_with_meta.RDS")
soil_depth_1 <- subset_samples(phy, matrix == "soil" & depth == "d1") %>%
  rarefy_even_depth(sample.size = 10000, rngseed = 3242343, verbose = F) %>%
  filter_taxa(function(x) sum(x) >= 3, T)
rm(phy)

tsne_soil_depth_1 <- tsne_phyloseq_ggplot(soil_depth_1, treatment = c('treatment'), perplexity = 10, circle = TRUE, colors = 'default') +
  scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw() +
  guides(fill=guide_legend(title="Treatment"))

tsne_soil_depth_1

only_prairie_samples <- subset_samples(soil_depth_1, in_plot_location %in% c("s6", "s7", "s8", "s9"))
only_prairie_samples_tsne <- tsne_phyloseq_ggplot(only_prairie_samples, treatment = c('treatment'), perplexity = 10, circle = TRUE, colors = 'default') +
  scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw() +
  guides(fill=guide_legend(title="Treatment"))

only_prairie_samples_tsne

only_crop_samples <- subset_samples(soil_depth_1, in_plot_location %in% c("s1", "s2", "s3", "s4", "s5"))
only_crop_samples_tsne <- tsne_phyloseq_ggplot(only_crop_samples, treatment = c('treatment'), perplexity = 10, circle = TRUE, colors = 'default') +
  scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw() +
  guides(fill=guide_legend(title="Treatment"))

only_crop_samples_tsne

