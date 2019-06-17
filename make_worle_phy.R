list.of.packages <- c("phyloseq", "tidyverse", "phylosmith")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

rm(list.of.packages)
rm(new.packages)

library(phyloseq)
library(tidyverse)
library(phylosmith)

phy <- readRDS("data/usda.RDS")
taxa_names(phy) <- paste0('ASV_', seq_along(taxa_names(phy))) #this renames the "OTU" to ASV1....ASV1000...
sample_data(phy) <- read.delim(sep="\t", file="data/meta.txt", row.names=1, header=TRUE)
colnames(sample_data(phy)) <- c("PI", "source", "date", "concentration")

worle <- subset_samples(phy, PI == "Jared")
worle
rm(phy)

sample_names(worle)
# Four sample names are wrong and need to be changed, these should be d2, not d1. That is the reason that a 2 was tagged on the end because they are duplicate names
#"P1-s1-d1-t2-2"
#"P1-s2-d1-t2-2"
#"P1-s3-d1-t2-2"

sample_names(worle)[9] <- "P1-s1-d2-t2"
sample_names(worle)[9]
sample_names(worle)[21] <- "P1-s2-d2-t2"
sample_names(worle)[32] <- "P1-s3-d2-t2"
sample_names(worle)[44] <- "P1-s4-d2-t2"

meta <- as.data.frame(read_csv(file = "data/working_on_meta.csv")) 

sample_names(worle) == meta$id
test <- data.frame(sample_data(worle)) 

row.names(meta) <- meta$id
row.names(test) == row.names(meta)

sample_data(worle) <- meta
head(sample_data(worle))

saveRDS(worle, file = "worle_with_meta.RDS")
