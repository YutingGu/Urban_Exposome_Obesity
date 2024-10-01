library(glmnet)
library(fake)
library(sharp)
library(dplyr)
library(tidyverse)
library(devtools)
library(ggplot2)
library(rstudioapi)
library(igraph)

path_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir = path_dir)
rm(path_dir)

# loading data
clusters <- read_csv("/rds/general/project/hda_23-24/live/CompEpi/group1/Dataset/cluster_id.csv")
clusters <- data.frame(clusters)
cleaned_ukb_with_protein <- readRDS("/rds/general/project/hda_23-24/live/CompEpi/group1/Dataset/cleaned_ukb_with_protein.rds")
protein <- readRDS("/rds/general/project/hda_23-24/live/CompEpi/group1/Dataset/Proteins_nd_imputed.rds")
protein <- data.frame(protein)

# merging cluster membership and ukbk data
ukb_cluster <- merge(cleaned_ukb_with_protein, clusters, by = "eid")
protein_names <- colnames(protein)
variables <- append(protein_names, c("age","sex","Cluster")) #age to be added
ukb_cluster <- ukb_cluster[,variables]

# creating binary variable for cluster membership
ukb_cluster$cluster0 <- ifelse(ukb_cluster$Cluster == 0, 1,0)
ukb_cluster$cluster1 <- ifelse(ukb_cluster$Cluster == 1, 1,0)
ukb_cluster$cluster2 <- ifelse(ukb_cluster$Cluster == 2, 1,0)
ukb_cluster$cluster3 <- ifelse(ukb_cluster$Cluster == 3, 1,0)

saveRDS(ukb_cluster, "/rds/general/project/hda_23-24/live/CompEpi/group1/Dataset/ukb_cluster.rds")

# stability selection

penalty_factor <- c(rep(1, length(protein_names)),0,0)
variables_to_include <- append(protein_names, c("age","sex"))
X <- ukb_cluster[,variables_to_include]

# CLUSTER 0 - stability selection 

Y_0 <- ukb_cluster$cluster0

# running stability selection
t0 <- Sys.time()
out_cluster0 <- VariableSelection(
  xdata = X,
  ydata = Y_0,
  verbose = FALSE,
  penalty.factor = penalty_factor,
  family = "binomial"
)
t1 <- Sys.time()
print(t1 - t0)
saveRDS(out_cluster0, "/rds/general/project/hda_23-24/live/CompEpi/group1/Results/stability_selection_output/stab_select_cluster0.rds")

# CLUSTER 1 - stability selection 

Y_1 <- ukb_cluster$cluster1

# running stability selection
t0 <- Sys.time()
out_cluster1 <- VariableSelection(
  xdata = X,
  ydata = Y_1,
  verbose = FALSE,
  penalty.factor = penalty_factor,
  family = "binomial"
)
t1 <- Sys.time()
print(t1 - t0)
saveRDS(out_cluster1, "/rds/general/project/hda_23-24/live/CompEpi/group1/Results/stability_selection_output/stab_select_cluster1.rds")

# CLUSTER 2 - stability selection 

Y_2 <- ukb_cluster$cluster2

# running stability selection
t0 <- Sys.time()
out_cluster2 <- VariableSelection(
  xdata = X,
  ydata = Y_2,
  verbose = FALSE,
  penalty.factor = penalty_factor,
  family = "binomial"
)
t1 <- Sys.time()
print(t1 - t0)
saveRDS(out_cluster2, "/rds/general/project/hda_23-24/live/CompEpi/group1/Results/stability_selection_output/stab_select_cluster2.rds")

# CLUSTER 3 - stability selection 

Y_3 <- ukb_cluster$cluster3

# running stability selection
t0 <- Sys.time()
out_cluster3 <- VariableSelection(
  xdata = X,
  ydata = Y_3,
  verbose = FALSE,
  penalty.factor = penalty_factor,
  family = "binomial"
)
t1 <- Sys.time()
print(t1 - t0)
saveRDS(out_cluster3, "/rds/general/project/hda_23-24/live/CompEpi/group1/Results/stability_selection_output/stab_select_cluster3.rds")







