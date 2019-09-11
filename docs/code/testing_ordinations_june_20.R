setwd("Documents/STRIPs/WORLE_17")

library(phyloseq)
phy <- readRDS("worle_with_meta.RDS")

# Testing Ordinations using multiple cores
# Loading required library and displaying core configuration
pkgTest <- function(x) {
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

pkgTest('doParallel')
pkgTest('phylosmith')
library(doParallel)
library(phylosmith)
library(tidyverse)
library(viridisLite)
library(viridis)

# Use multiple cores
detectCores(all.tests=TRUE)
# 48

# Setting up and registering the cluster
cl <- makeCluster(detectCores(all.tests=TRUE))
registerDoParallel(cl)

water <- subset_samples(phy, matrix == "water") %>%
  filter_taxa(function(x) sum(x) >= 15, T)
              
water_tsne <- tsne_phyloseq_ggplot(water, treatment = c('treatment'), perplexity = 10, circle = T, colors = 'default') +
  scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw() +
  guides(fill = guide_legend(title = "Treatment"))

water_tsne

soil_depth_1 <- subset_samples(phy, matrix == "soil" & depth == "d1" & day %in% c("t0", "t2")) %>%
  filter_taxa(function(x) sum(x) > 20, T)

tsne_soil_depth_1 <- tsne_phyloseq_ggplot(soil_depth_1, treatment = c('treatment'), perplexity = 10, circle = TRUE, colors = 'default') +
  scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw() +
  guides(fill=guide_legend(title="Treatment"))

tsne_soil_depth_1

nmds_soil_depth_1 <- nmds_phyloseq_ggplot(soil_depth_1, treatment = c('treatment'), circle = T, colors = 'default') +
  scale_fill_viridis(discrete = T, option = "viridis") + ggplot2::theme_bw() +
  guides(fill = guide_legend(title = "Treatment"))

nmds_soil_depth_1

# Shut down cluster
stopCluster(cl)
rm(cl)
registerDoSEQ()
