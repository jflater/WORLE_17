library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(data.table)

# Wafergen text file, results from the facility
dat <- read.table("data/wafergen_results/STRIPs.txt", sep="\t", header=T)[,c("Assay","Sample","Ct")]

# User generated metadata file, be sure the your Sample column here matches the Sample column in your results
meta <- read.csv("data/meta_data.txt", sep = "\t", header = T)

# Merge based on Sample
merged_data <- merge(dat, meta, by.x = "Sample", by.y = "Sample", all.x=T)

# Filter Ct values that are more than 28, these are poor reads and we will not consider them
hits <- merged_data %>%
  filter(Ct != "NA", Ct <= 28) %>%
  group_by(ID, Assay) %>%
  summarize(mean_ct = mean(Ct, na.rm = T))

hits <- merged_data %>%
  filter(Ct != "NA", Ct <= 28) %>%
  group_by(ID, Assay)

mean_ct <- hits %>%
  summarize(mean_ct = mean(Ct, na.rm = T))

df %>%
  group_by(Label, Code) %>%
  mutate_if(is.numeric, sum) %>%
  mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>%
  distinct()
# We don't really care about the marker genes and the primers for specific bacteria, remove those with the ! added to the column of interest
arg_hits <- hits %>%
  filter(!Assay %in% c("16S new 2", "16S old 1", "Bacteroidetes", "Firmicutes"))

levels(arg_hits$ID)
IDS

soil_hits <- arg_hits %>%
  filter(ID %in% c("P3_S5", "P4_S9", "P8_S2")) %>%
  select(Assay) 

Crop_Prairie_Interface_soil <- hits %>%
  filter(Treatment == "Crop_Prairie_Interface")

Swine_manure <- hits %>%
  filter(Treatment == "Swine") 
  

soil <- soil %>%
  filter(Ct <= 28) %>%
  arrange(Ct)

soil$Assay <- as.vector(soil$Assay)
soil$Assay <- factor(soil$Assay)

ggplot(Swine_manure, 
       aes(x = Assay, y = Ct)) +
  geom_boxplot() +
  facet_grid(~Block) +
  rotate_x_text(angle = 45) 
