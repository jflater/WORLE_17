library(tidyverse)
library(phyloseq)

phy <- readRDS("data/usda.RDS")
taxa_names(phy) <- paste0('ASV_', seq_along(taxa_names(phy))) #this renames the "OTU" to ASV1....ASV1000...
sample_data(phy) <- read.delim(sep="\t", file="data/meta.txt", row.names=1, header=TRUE)
colnames(sample_data(phy)) <- c("PI", "source", "date", "concentration")

sample_names(phy)

worle <- subset_samples(phy, PI == "Jared")

worle_meta <- data.frame(phyloseq::sample_data(worle)) %>%
  rownames_to_column() %>%
  select(sample_id = rowname)

write_excel_csv2(worle_meta, path = "data/worle_meta.txt", col_names = T)
