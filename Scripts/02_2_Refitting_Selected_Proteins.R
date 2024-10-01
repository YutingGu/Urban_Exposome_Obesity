rm(list=ls())

library(glmnet)
library(igraph)
library(pheatmap)
library(fake)
library(sharp)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(rstudioapi)
library(gplots)

stab_cluster0 <- readRDS("/rds/general/project/hda_23-24/live/CompEpi/group1/Results/stability_selection_output/stab_select_cluster0.rds")
stab_cluster1 <- readRDS("/rds/general/project/hda_23-24/live/CompEpi/group1/Results/stability_selection_output/stab_select_cluster1.rds")
stab_cluster2 <- readRDS("/rds/general/project/hda_23-24/live/CompEpi/group1/Results/stability_selection_output/stab_select_cluster2.rds")
stab_cluster3 <- readRDS("/rds/general/project/hda_23-24/live/CompEpi/group1/Results/stability_selection_output/stab_select_cluster3.rds")
# loading ukb + protein + cluster data
ukb_cluster <- readRDS("/rds/general/project/hda_23-24/live/CompEpi/group1/Dataset/ukb_cluster.rds")
# proteins
protein <- readRDS("/rds/general/project/hda_23-24/live/CompEpi/group1/Dataset/Proteins_nd_imputed.rds")
protein <- data.frame(protein)


# Cluster 0

# calibration plot
CalibrationPlot(stab_cluster0)

# calibrated parameters 
hat_params0 <- Argmax(stab_cluster0)
selection_threshold0 <- hat_params0[2]
print(hat_params0)

# selection proportion with calibrated model
selprop0 <- SelectionProportions(stab_cluster0)
print(selprop0)
selected_proteins0 <- selprop0[selprop0 >= selection_threshold0]
selected_proteins0 <- sort(selected_proteins0, decreasing = TRUE) #ordered by descending values of proportion selection
nb_selected_proteins0 <- length(selected_proteins0)

# visualisation of selection proportions
selprop0 <- sort(selprop0, decreasing = TRUE)
par(mar = c(10, 5, 1, 1))
plot(selprop0,
     type = "h", lwd = 3, las = 1, xlab = "Proteins", ylab = "Selection Proportion", xaxt = "n",
     col = ifelse(selprop0 >= hat_params0[2], yes = "red", no = "grey"), cex.lab = 1.5
)
abline(h = hat_params0[2], lty = 2, col = "darkred")
#for (i in 1:length(selprop0)) {
#  axis(
#    side = 1, at = i, labels = names(selprop0)[i], las = 2,
#    col = ifelse(selprop0[i] >= hat_params0[2], yes = "red", no = "transparent"),
#    col.axis = ifelse(selprop0[i] >= hat_params0[2], yes = "red", no = "transparent")
#  )
#}

# running logistic regression with selected proteins cluster 0
selected_proteins_conf0 <- append(names(selected_proteins0), c("age","sex","cluster0"))
ukb_cluster0 <- ukb_cluster[,selected_proteins_conf0 ]
log_reg0 <- glm(cluster0 ~ . , family = binomial(link = "logit"), ukb_cluster0)

coeff0 <- exp(log_reg0$coefficients) # coefficients
coeff0 <- sort(coeff0, decreasing = TRUE)

# Cluster 1

# calibration plot
CalibrationPlot(stab_cluster1)

# calibrated parameters 
hat_params1 <- Argmax(stab_cluster1)
selection_threshold1 <- hat_params1[2]
print(hat_params1)

# selection proportion with calibrated model
selprop1 <- SelectionProportions(stab_cluster1)
#print(selprop1)
selected_proteins1 <- selprop1[selprop1 >= selection_threshold1]
selected_proteins1 <- sort(selected_proteins1, decreasing = TRUE) #ordered by descending values of proportion selection
nb_selected_proteins1 <- length(selected_proteins1)

# visualisation of selection proportions
selprop1 <- sort(selprop1, decreasing = TRUE)
par(mar = c(10, 5, 1, 1))
plot(selprop1,
     type = "h", lwd = 3, las = 1, xlab = "", ylab = "Selection Proportion", xaxt = "n",
     col = ifelse(selprop1 >= hat_params1[2], yes = "red", no = "grey"), cex.lab = 1.5
)
abline(h = hat_params1[2], lty = 2, col = "darkred")
for (i in 1:length(selprop1)) {
  axis(
    side = 1, at = i, labels = names(selprop1)[i], las = 2,
    col = ifelse(selprop1[i] >= hat_params1[2], yes = "red", no = "grey"),
    col.axis = ifelse(selprop1[i] >= hat_params1[2], yes = "red", no = "grey")
  )
}

# running logistic regression with selected proteins cluster 1
selected_proteins_conf1 <- append(names(selected_proteins1), c("age","sex","cluster1"))
ukb_cluster1 <- ukb_cluster[,selected_proteins_conf1 ]
log_reg1 <- glm(cluster1 ~ . , family = binomial(link = "logit"), ukb_cluster1)

coeff1 <- exp(log_reg1$coefficients) # coefficients
coeff1 <- sort(coeff1, decreasing = TRUE)

# Cluster 2

# calibration plot
CalibrationPlot(stab_cluster2)

# calibrated parameters 
hat_params2 <- Argmax(stab_cluster2)
selection_threshold2 <- hat_params2[2]
print(hat_params2)

# selection proportion with calibrated model
selprop2 <- SelectionProportions(stab_cluster2)
#print(selprop1)
selected_proteins2 <- selprop2[selprop2 >= selection_threshold2]
selected_proteins2 <- sort(selected_proteins2, decreasing = TRUE) #ordered by descending values of proportion selection
nb_selected_proteins2 <- length(selected_proteins2)

# visualisation of selection proportions
par(mar = c(10, 5, 1, 1))
plot(selprop2,
     type = "h", lwd = 3, las = 1, xlab = "", ylab = "Selection Proportion", xaxt = "n",
     col = ifelse(selprop2 >= hat_params2[2], yes = "red", no = "grey"), cex.lab = 1.5
)
abline(h = hat_params2[2], lty = 2, col = "darkred")
for (i in 1:length(selprop2)) {
  axis(
    side = 1, at = i, labels = names(selprop2)[i], las = 2,
    col = ifelse(selprop2[i] >= hat_params2[2], yes = "red", no = "grey"),
    col.axis = ifelse(selprop2[i] >= hat_params2[2], yes = "red", no = "grey")
  )
}

# running logistic regression with selected proteins cluster 2
selected_proteins_conf2 <- append(names(selected_proteins2), c("age","sex","cluster2"))
ukb_cluster2 <- ukb_cluster[,selected_proteins_conf2 ]
log_reg2 <- glm(cluster2 ~ . , family = binomial(link = "logit"), ukb_cluster2)

coeff2 <- exp(log_reg2$coefficients)
coeff2 <- sort(coeff2, decreasing = TRUE)

# Cluster 3

# calibration plot
CalibrationPlot(stab_cluster3)

# calibrated parameters 
hat_params3 <- Argmax(stab_cluster3)
selection_threshold3 <- hat_params3[2]
print(hat_params3)

# selection proportion with calibrated model
selprop3 <- SelectionProportions(stab_cluster3)
#print(selprop1)
selected_proteins3 <- selprop3[selprop3 >= selection_threshold3]
selected_proteins3 <- sort(selected_proteins3, decreasing = TRUE) #ordered by descending values of proportion selection
nb_selected_proteins3 <- length(selected_proteins3)


# visualisation of selection proportions
par(mar = c(10, 5, 1, 1))
plot(selprop3,
     type = "h", lwd = 3, las = 1, xlab = "", ylab = "Selection Proportion", xaxt = "n",
     col = ifelse(selprop3 >= hat_params3[2], yes = "red", no = "grey"), cex.lab = 1.5
)
abline(h = hat_params3[2], lty = 2, col = "darkred")
for (i in 1:length(selprop3)) {
  axis(
    side = 1, at = i, labels = names(selprop3)[i], las = 2,
    col = ifelse(selprop3[i] >= hat_params3[2], yes = "red", no = "grey"),
    col.axis = ifelse(selprop3[i] >= hat_params3[2], yes = "red", no = "grey")
  )
}

# running logistic regression with selected proteins cluster 3
selected_proteins_conf3 <- append(names(selected_proteins3), c("age","sex","cluster3"))
ukb_cluster3 <- ukb_cluster[,selected_proteins_conf3]
log_reg3 <- glm(cluster3 ~ . , family = binomial(link = "logit"), ukb_cluster3)

coeff3 <- exp(log_reg3$coefficients)
coeff3 <- sort(coeff3, decreasing =TRUE)


# Heatmap 

# initializing heatmap
heatmap <- data.frame(matrix(nrow = length(protein), ncol = 4))
colnames(heatmap) <- c("cluster0","cluster1","cluster2","cluster3")
heatmap$protein <- colnames(protein)

# coefficients from model 0
beta_0 <- as.data.frame(coeff0)
beta_0$protein <- names(coef(log_reg0))
heatmap <- merge(heatmap, beta_0, by = "protein", all.x = TRUE)
heatmap$cluster0 <- heatmap$coeff0
heatmap <- select(heatmap, select = -c("coeff0"))

# coefficients from model 1
beta_1 <- as.data.frame(coeff1)
beta_1$protein <- names(coef(log_reg1))
heatmap <- merge(heatmap, beta_1, by = "protein", all.x = TRUE)
heatmap$cluster1 <- heatmap$coeff1
heatmap <- select(heatmap, select = -c("coeff1"))

# coefficients from model 2
beta_2 <- as.data.frame(coeff2)
beta_2$protein <- names(coef(log_reg2))
heatmap <- merge(heatmap, beta_2, by = "protein", all.x = TRUE)
heatmap$cluster2 <- heatmap$coeff2
heatmap <- select(heatmap, select = -c("coeff2"))

# coefficients from model 3
beta_3 <- as.data.frame(coeff3)
beta_3$protein <- names(coef(log_reg3))
heatmap <- merge(heatmap, beta_3, by = "protein", all.x = TRUE)
heatmap$cluster3 <- heatmap$coeff3
heatmap <- select(heatmap, select = -c("coeff3"))

# plot coeff 
beta_0 <- beta_0[!(beta_0$protein %in% c("age","sex","(Intercept)")),]
beta_1 <- beta_1[!(beta_1$protein %in% c("age","sex","(Intercept)")),]
beta_2 <- beta_2[!(beta_2$protein %in% c("age","sex","(Intercept)")),]
beta_3 <- beta_3[!(beta_3$protein %in% c("age","sex","(Intercept)")),]


par(mfrow = c(2, 2))
plot(beta_0$coeff0, 
     xlab = "Protein Cluster 0", ylab = "Beta values")
abline(h = 1, col = "red")
plot(beta_1$coeff1, 
     xlab = "Protein Cluster 1 ", ylab = "Beta values")
abline(h = 1, col = "red")
plot(beta_2$coeff2, 
     xlab = "Protein Cluster 2", ylab = "Beta values")
abline(h = 1, col = "red")
plot(beta_3$coeff3, 
     xlab = "Protein Cluster 3", ylab = "Beta values")
abline(h = 1, col = "red")


# replacing NA with 0
heatmap$cluster0 <- ifelse(is.na(heatmap$cluster0) == TRUE, 0, heatmap$cluster0)
heatmap$cluster1 <- ifelse(is.na(heatmap$cluster1) == TRUE, 0, heatmap$cluster1)
heatmap$cluster2 <- ifelse(is.na(heatmap$cluster2) == TRUE, 0, heatmap$cluster2)
heatmap$cluster3 <- ifelse(is.na(heatmap$cluster3) == TRUE, 0, heatmap$cluster3)

# number of times a protein was selected
heatmap$nb_times_selected <- rowSums(heatmap[, c("cluster0", "cluster1", "cluster2", "cluster3")] != 0)
heatmap <- heatmap[order(heatmap$nb_times_selected, decreasing = TRUE), ]
table_nb_selected_proteins <- table(heatmap$nb_times_selected)

# plotting heatmap coefficient
map_selected_protein <- heatmap[heatmap$nb_times_selected != 0,] #removing proteins never selected across the 4 clusters
protein_names <- map_selected_protein$protein
map_selected_protein <- map_selected_protein[, c(2:5)] #removing 'protein' and 'nb_selected_protein' columns for heatmap
map_selected_protein <- as.matrix(map_selected_protein)


custom_palette <- function(value) {
  if (value == 0) {
    return("white")  # White for value = 0
  } else if (value > 0 & value <= 2.5) {
    return(colorRampPalette(c("white", "blue"))(256))  # Custom palette for values between 0 and 2.5
  } else {
    return("blue")  # Default color for values above 2.5
  }
}
# Use custom color palette in heatmap
heatmap(map_selected_protein, 
        Rowv = NA, 
       # Colv = as.dendrogram(as.dist(1 - cor(t(map_selected_protein)))), 
       col = sapply(map_selected_protein, custom_palette), 
        scale = "column", 
        ylab = "Proteins", 
        main = "Betas of selected proteins per cluster", 
        labRow = protein_names)

heatmap(map_selected_protein, 
        Rowv = NA, 
        # Colv = as.dendrogram(as.dist(1 - cor(t(map_selected_protein)))), 
        col = ifelse(map_selected_protein < 0.000001, "white", 
                     ifelse(map_selected_protein > 0.001 & map_selected_protein <= 3, 
                            colorRampPalette(c("red","white"))(256), "pink")), 
        scale = "column", 
        ylab = "Proteins", 
        main = "Betas of selected proteins per cluster", 
        labRow = protein_names)


# % of commonly selected proteins across the clusters
num_clusters <- 4
common_proteins <- matrix(0, ncol = num_clusters, nrow = num_clusters, dimnames = list(paste0("cluster", 0:(num_clusters-1)), paste0("cluster", 0:(num_clusters-1))))
# Calculate the number of common proteins for each pair of clusters
for (i in 1:num_clusters) {
  for (j in 1:num_clusters) {
    common_proteins[i, j] <- sum(heatmap[, paste0("cluster", i-1)] != 0 & heatmap[, paste0("cluster", j-1)] != 0)
  }
}
# Normalizing by dividing each column by its corresponding diagonal value
diagonal_values <- diag(common_proteins)
for (i in 1:num_clusters) {
  common_proteins[,i] <- 100 * common_proteins[,i] / diagonal_values[i]
}
common_proteins <- round(common_proteins, 0)
# Define the order of rows and columns based on the order of clusters in common_proteins matrix
row_order <- colnames(common_proteins)
col_order <- rownames(common_proteins)
# Create heatmap with specified row and column orders
heatmap.2(common_proteins, 
          trace = "none",
          dendrogram = "none",
          Colv = FALSE,
          Rowv = list(row_order),
          scale = "none",
          key = TRUE, keysize = 1.5,
          key.title = NA,
          key.xlab = "Percentage",
          key.ylab = NA,
          density.info = "none",
          cexRow = 0.8, cexCol = 0.8, margins = c(10,10),
          main = "Commonly selected proteins percentage ",
          labRow = paste0("Cluster", 0:(num_clusters-1)),
          labCol = paste0("Cluster", 0:(num_clusters-1)),
          cellnote = common_proteins, # Add cell annotations
          notecol = "black", notecex = 0.5, # Customize cell annotation properties
          cexCell = 3,
          col = colorRampPalette(c("lightyellow","orange", "red"))(100)
)

# Investigation 

# Cluster 2 -> link with obesity 
par(mfrow = c(1, 1))
boxplot(ukb_cluster[ukb_cluster$Cluster == 2, "CEACAM1"],ukb_cluster[ukb_cluster$Cluster == 0, "CEACAM1"])
t.test(ukb_cluster[ukb_cluster$Cluster == 2, "CEACAM1"], ukb_cluster[ukb_cluster$Cluster == 3, "CEACAM1"])
ukb_cluster$cluster2 <- factor(ukb_cluster$cluster2)
obesity <- lm(CEACAM1 ~ cluster2, ukb_cluster)
summary(obesity)



