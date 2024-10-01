## Part 3 - Univariate analysis Obesity + Covars ~ clusters 

library(tidyverse)
library(ggplot2)
library(rstudioapi)
library(readr)
library(dplyr)
library(tidyr)
library(lme4)
library(reshape2)

# set working directory
#setwd("/rds/general/project/hda_23-24/live/CompEpi/group1/")

path_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir = path_dir)
rm(path_dir)


# load dataset
cluster_id <- read.csv("../Dataset/cluster_id.csv")
cluster_id$eid <- as.character(cluster_id$eid)
cleaned_ukb_with_protein <- readRDS("../Dataset/cleaned_ukb_with_protein.rds")

# join the label with the 
data <- cleaned_ukb_with_protein %>% left_join(cluster_id, join_by(eid==eid))

# convert the cluster_id to one-hot encode
data$Cluster_0 <- ifelse(data$Cluster==0, 1, 0)
data$Cluster_1 <- ifelse(data$Cluster==1, 1, 0)
data$Cluster_2 <- ifelse(data$Cluster==2, 1, 0)
data$Cluster_3 <- ifelse(data$Cluster==3, 1, 0)

Covars <- names(data)[17:23]
data_Ob = data[c("eid", Covars,"Cluster_0","Cluster_1","Cluster_2","Cluster_3","age")]

data_Ob$sex <- as.factor(data_Ob$sex)


############ missing values ###############################################
sum(is.na(data_Ob))
colSums(is.na(data_Ob))
aggregate(. ~ cluster, data = df, FUN = function(x) sum(is.na(x)))

# create covar columns
covariate_columns <- c("sex", "waist_circ", "hh_income", "bmi", "hba1c", "bp_sys", "ldl")

# Create a list to hold the results
missing_values_by_cluster <- list()

# Calculate missing values for each cluster
missing_values_by_cluster$Cluster_0 <- colSums(is.na(data_Ob[data_Ob$Cluster_0 == 1, covariate_columns]))
missing_values_by_cluster$Cluster_1 <- colSums(is.na(data_Ob[data_Ob$Cluster_1 == 1, covariate_columns]))
missing_values_by_cluster$Cluster_2 <- colSums(is.na(data_Ob[data_Ob$Cluster_2 == 1, covariate_columns]))
missing_values_by_cluster$Cluster_3 <- colSums(is.na(data_Ob[data_Ob$Cluster_3 == 1, covariate_columns]))

# Print the results
print(missing_values_by_cluster)

# Calculate the percentage of missing values for each cluster and covariate
calculate_percent_missing <- function(cluster_column) {
  colSums(is.na(data_Ob[data_Ob[[cluster_column]] == 1, covariate_columns])) / 
    sum(data_Ob[[cluster_column]] == 1) * 100
}


# Calculate percentages
percent_missing_cluster_0 <- calculate_percent_missing("Cluster_0")
percent_missing_cluster_1 <- calculate_percent_missing("Cluster_1")
percent_missing_cluster_2 <- calculate_percent_missing("Cluster_2")
percent_missing_cluster_3 <- calculate_percent_missing("Cluster_3")

# Create a data frame for plotting
df_missing_percent <- data.frame(
  Cluster_0 = percent_missing_cluster_0,
  Cluster_1 = percent_missing_cluster_1,
  Cluster_2 = percent_missing_cluster_2,
  Cluster_3 = percent_missing_cluster_3,
  Covariate = covariate_columns
)

# Melt the data frame for ggplot
df_melted <- melt(df_missing_percent, id.vars = "Covariate", variable.name = "Cluster", value.name = "Percent_Missing")

ggplot(df_melted, aes(x = Covariate, y = Percent_Missing, fill = Cluster)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Covariate", y = "Percentage Missing", title = "Percentage of Missing Data by Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Larger font for x-axis labels
        axis.text.y = element_text(size = 14),                         # Larger font for y-axis labels
        axis.title.x = element_text(size = 16),                        # Larger font for x-axis title
        axis.title.y = element_text(size = 16),                        # Larger font for y-axis title
        plot.title = element_text(size = 18, hjust = 0.5)) +           # Larger and centered plot title
  scale_y_continuous(limits = c(0, 20))  # Set y-axis limits from 0 to 20%

## remove missing data
data_Ob <- na.omit(data_Ob)
dim(data_Ob)

################################################################################
## Logistic model 1:  
## (m0) Cluster(i) ~ age + gender 
## (m1) Cluster(i) ~ age + gender + Obesity_Covar

model0_0 <- glm(Cluster_0 ~ age + sex, data = data_Ob, family = binomial())
model0_1 <- glm(Cluster_0 ~ age + sex + bmi, data = data_Ob, family = binomial())

summary(model0_0)
summary(model0_1)

model1_0 <- glm(Cluster_1 ~ age + sex, data = data_Ob, family = binomial())
model1_1 <- glm(Cluster_1 ~ age + sex + bmi, data = data_Ob, family = binomial())

summary(model1_0)
summary(model1_1)


model2_0 <- glm(Cluster_2 ~ age + sex, data = data_Ob, family = binomial())
model2_1 <- glm(Cluster_2 ~ age + sex + bmi, data = data_Ob, family = binomial())

summary(model2_0)
summary(model2_1)

model3_0 <- glm(Cluster_3 ~ age + sex, data = data_Ob, family = binomial())
model3_1 <- glm(Cluster_3 ~ age + sex + bmi, data = data_Ob, family = binomial())

summary(model3_0)
summary(model3_1)

################################################################################
## Regression model 2:  
## (m1) BMI ~ Cluster + age + gender + covars

## stepwise cluster addition
lm_bmi_covars <- lm(bmi ~ age + sex + hh_income, data = data_Ob)
lm_bmi_ClusCovar <- lm(bmi ~ Cluster_0 + age + sex + hh_income, data = data_Ob)
summary(lm_bmi_covars)
summary(lm_bmi_ClusCovar)

## stepwise covar addition to test model improvement 
lm_bmi_c0_1 <- lm(bmi ~ Cluster_0, data = data_Ob)
summary(lm_bmi_c0_1)

lm_bmi_c0_2 <- lm(bmi ~ Cluster_0 + age + sex, data = data_Ob)
summary(lm_bmi_c0_2)

lm_bmi_c0_3 <- lm(bmi ~ Cluster_0 + age + sex + hh_income, data = data_Ob) 
summary(lm_bmi_c0_3)

## bmi ~ clusters + covars 
lm_bmi_c0 <- lm(bmi ~ Cluster_0 + age + sex , data = data_Ob) 
lm_bmi_c1 <- lm(bmi ~ Cluster_1 + age + sex , data = data_Ob) 
lm_bmi_c2 <- lm(bmi ~ Cluster_2 + age + sex , data = data_Ob) 
lm_bmi_c3 <- lm(bmi ~ Cluster_3 + age + sex , data = data_Ob) 

summary(lm_bmi_c0)
summary(lm_bmi_c1)
summary(lm_bmi_c2)
summary(lm_bmi_c3)


## bmi ~ clusters + covars + income 
lm_bmi_income_c0 <- lm(bmi ~ Cluster_0 + age + sex + hh_income, data = data_Ob) 
lm_bmi_income_c1 <- lm(bmi ~ Cluster_1 + age + sex + hh_income, data = data_Ob) 
lm_bmi_income_c2 <- lm(bmi ~ Cluster_2 + age + sex + hh_income, data = data_Ob) 
lm_bmi_income_c3 <- lm(bmi ~ Cluster_3 + age + sex + hh_income, data = data_Ob) 

summary(lm_bmi_income_c0)
summary(lm_bmi_income_c1)
summary(lm_bmi_income_c2)
summary(lm_bmi_income_c3)

 
## Hba1c ~ clusters + covars 
lm_a1c_c0 <- lm(hba1c ~ Cluster_0 + age + sex , data = data_Ob) 
lm_a1c_c1 <- lm(hba1c ~ Cluster_1 + age + sex , data = data_Ob) 
lm_a1c_c2 <- lm(hba1c ~ Cluster_2 + age + sex , data = data_Ob) 
lm_a1c_c3 <- lm(hba1c ~ Cluster_3 + age + sex , data = data_Ob) 

summary(lm_a1c_c0)
summary(lm_a1c_c1)
summary(lm_a1c_c2)
summary(lm_a1c_c3)

## Hba1c ~ clusters + covars + income
lm_a1c_income_c0 <- lm(hba1c ~ Cluster_0 + age + sex + hh_income, data = data_Ob) 
lm_a1c_income_c1 <- lm(hba1c ~ Cluster_1 + age + sex + hh_income, data = data_Ob) 
lm_a1c_income_c2 <- lm(hba1c ~ Cluster_2 + age + sex + hh_income, data = data_Ob) 
lm_a1c_income_c3 <- lm(hba1c ~ Cluster_3 + age + sex + hh_income, data = data_Ob) 

summary(lm_a1c_income_c0)
summary(lm_a1c_income_c1)
summary(lm_a1c_income_c2)
summary(lm_a1c_income_c3)

## bp_sys ~ clusters + covars + income
lm_bp_income_c0 <- lm(bp_sys ~ Cluster_0 + age + sex + hh_income, data = data_Ob) 
lm_bp_income_c1 <- lm(bp_sys ~ Cluster_1 + age + sex + hh_income, data = data_Ob) 
lm_bp_income_c2 <- lm(bp_sys ~ Cluster_2 + age + sex + hh_income, data = data_Ob) 
lm_bp_income_c3 <- lm(bp_sys ~ Cluster_3 + age + sex + hh_income, data = data_Ob) 

summary(lm_bp_income_c0)
summary(lm_bp_income_c1)
summary(lm_bp_income_c2)
summary(lm_bp_income_c3)

## ldl ~ clusters + covars + income

lm_ldl_income_c0 <- lm(ldl ~ Cluster_0 + age + sex + hh_income, data = data_Ob) 
lm_ldl_income_c1 <- lm(ldl ~ Cluster_1 + age + sex + hh_income, data = data_Ob) 
lm_ldl_income_c2 <- lm(ldl ~ Cluster_2 + age + sex + hh_income, data = data_Ob) 
lm_ldl_income_c3 <- lm(ldl ~ Cluster_3 + age + sex + hh_income, data = data_Ob) 

summary(lm_ldl_income_c0)
summary(lm_ldl_income_c1)
summary(lm_ldl_income_c2)
summary(lm_ldl_income_c3)


# heatmap with ORs 

# Define outcomes and clusters
outcomes <- c("bmi", "hba1c", "bp_sys", "ldl")
cluster_suffixes <- c("c0", "c1", "c2", "c3")  # Correct suffixes for cluster identification

# Initialize a results dataframe
results <- data.frame(Outcome = character(), Cluster = character(), Coefficient = numeric(), pValue = numeric(), stringsAsFactors = FALSE)

# Function to extract coefficients and p-values
extract_coef <- function(outcome, cluster_suffix) {
  # Construct the model name
  model_name <- paste("lm", outcome, "income", cluster_suffix, sep = "_")
  
  # Check if the model exists and if so, proceed with extraction
  if (exists(model_name, envir = .GlobalEnv)) {
    model <- get(model_name)
    coef_info <- summary(model)$coefficients
    
    # Extract coefficient and p-value for the cluster, assuming coefficient row matches exactly as shown in summary
    coef_label <- paste("Cluster", gsub("c", "", cluster_suffix), sep = "_")  # This should now correctly form labels like "Cluster_0"
    
    if (coef_label %in% rownames(coef_info)) {
      coef <- coef_info[coef_label, "Estimate"]
      pval <- coef_info[coef_label, "Pr(>|t|)"]
      return(data.frame(Outcome = outcome, Cluster = coef_label, Coefficient = coef, pValue = pval))
    } else {
      print(paste("Coefficient name", coef_label, "not found in the model summary for", model_name))
    }
  } else {
    print(paste("Model does not exist:", model_name))
  }
  return(data.frame(Outcome = outcome, Cluster = cluster_suffix, Coefficient = NA, pValue = 1))
}

# Populate the results dataframe using loops
for (outcome in outcomes) {
  for (cluster_suffix in cluster_suffixes) {
    results <- rbind(results, extract_coef(outcome, cluster_suffix))
  }
}

# Print results to check the correct extraction
print(results)

# Filter results to include only significant coefficients (pValue < 0.05)
significant_results <- results[results$pValue < 0.05, ]
significant_results$Outcome <- factor(significant_results$Outcome,
                                      levels = c("bp" ,"ldl" ,"a1c", "bmi"))  # Set the levels in the correct order

# Create the heatmap
heatmap_plot <- ggplot(significant_results, aes(x = Cluster, y = Outcome, fill = Coefficient)) +
  geom_tile(color = "white", linewidth = 0.1) + # Defines the tile borders
  geom_text(aes(label = sprintf("%.2f", Coefficient)), color = "black", size = 3) +  # Add coefficient numbers
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       na.value = "white") +  # No limits set, color scales naturally with data
  labs(title = "Heatmap of Significant Coefficients", x = "", y = "", fill = "Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),  # Center x labels under each column
        axis.text.y = element_text(hjust = 1),  # Ensure y labels are aligned to the left of the plot
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        axis.ticks.x = element_blank(),  # Remove x-axis ticks to clean up the plot
        axis.ticks.y = element_blank()) +  # Remove y-axis ticks for a clean look
  scale_x_discrete(position = "top") +  # Ensure x-axis labels are at the top of the plot
  coord_fixed(ratio = 1)  # Ensure that tiles are square

## bigger fonts 
heatmap_plot <- ggplot(significant_results, aes(x = Cluster, y = Outcome, fill = Coefficient)) +
  geom_tile(color = "white", linewidth = 0.1) + # Defines the tile borders
  geom_text(aes(label = sprintf("%.2f", Coefficient)), color = "black", size = 6) +  # Add coefficient numbers
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       na.value = "white") +  # No limits set, color scales naturally with data
  labs(title = "Heatmap of Significant Coefficients", x = "", y = "", fill = "Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"),  # Center x labels under each column
        axis.text.y = element_text(hjust = 1, size = 16, face = "bold"),  # Ensure y labels are aligned to the left of the plot
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        axis.ticks.x = element_blank(),  # Remove x-axis ticks to clean up the plot
        axis.ticks.y = element_blank()) +  # Remove y-axis ticks for a clean look
  scale_x_discrete(position = "top") +  # Ensure x-axis labels are at the top of the plot
  coord_fixed(ratio = 1)  # Ensure that tiles are square

# Print plot directly to view interactively
print(heatmap_plot)

########################################################################
## Obesity + traits 

# Assuming your data frame is named data_Ob
data_Ob$bmi_hi <- ifelse(data_Ob$bmi > 25, 1, 0)
data_Ob$hba1c_hi <- ifelse(data_Ob$hba1c > 42, 1, 0)
data_Ob$bp_sys_hi <- ifelse(data_Ob$bp_sys > 140, 1, 0)
data_Ob$ldl_hi <- ifelse(data_Ob$ldl > 3, 1, 0)

### logistic models
# bmi_hi ~ clusters + covars + income 
glm_c0_bmi_hi <- glm(bmi_hi ~ Cluster_0 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c1_bmi_hi <- glm(bmi_hi ~ Cluster_1 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c2_bmi_hi <- glm(bmi_hi ~ Cluster_2 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c3_bmi_hi <- glm(bmi_hi ~ Cluster_3 + age + sex + hh_income, data = data_Ob, family = binomial())

summary(glm_c0_bmi_hi)
summary(glm_c1_bmi_hi)
summary(glm_c2_bmi_hi)
summary(glm_c3_bmi_hi)

c0_bmi_OR <- exp(coef(glm_c0_bmi_hi))
c1_bmi_OR <- exp(coef(glm_c1_bmi_hi))
c2_bmi_OR <- exp(coef(glm_c2_bmi_hi))
c3_bmi_OR <- exp(coef(glm_c3_bmi_hi))

c0_bmi_OR 
c1_bmi_OR 
c2_bmi_OR 
c3_bmi_OR 

# conf_int <- confint(glm_c0_bmihi)  # Default 95% CI
# exp_conf_int <- exp(conf_int)  # Exponentiating the confidence intervals

# hba1c_hi ~ clusters + covars + income 
glm_c0_hba1c_hi <- glm(hba1c_hi ~ Cluster_0 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c1_hba1c_hi <- glm(hba1c_hi ~ Cluster_1 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c2_hba1c_hi <- glm(hba1c_hi ~ Cluster_2 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c3_hba1c_hi <- glm(hba1c_hi ~ Cluster_3 + age + sex + hh_income, data = data_Ob, family = binomial())

summary(glm_c0_hba1c_hi)
summary(glm_c1_hba1c_hi)
summary(glm_c2_hba1c_hi)
summary(glm_c3_hba1c_hi)

c0_a1c_OR <- exp(coef(glm_c0_hba1c_hi))
c1_a1c_OR <- exp(coef(glm_c1_hba1c_hi))
c2_a1c_OR <- exp(coef(glm_c2_hba1c_hi))
c3_a1c_OR <- exp(coef(glm_c3_hba1c_hi))

c0_a1c_OR 
c1_a1c_OR 
c2_a1c_OR 
c3_a1c_OR

# ldl_hi ~ clusters + covars + income 
glm_c0_ldl_hi <- glm(ldl_hi ~ Cluster_0 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c1_ldl_hi <- glm(ldl_hi ~ Cluster_1 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c2_ldl_hi <- glm(ldl_hi ~ Cluster_2 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c3_ldl_hi <- glm(ldl_hi ~ Cluster_3 + age + sex + hh_income, data = data_Ob, family = binomial())

summary(glm_c0_ldl_hi)
summary(glm_c1_ldl_hi)
summary(glm_c2_ldl_hi)
summary(glm_c3_ldl_hi)

c0_ldl_OR <- exp(coef(glm_c0_ldl_hi))
c1_ldl_OR <- exp(coef(glm_c1_ldl_hi))
c2_ldl_OR <- exp(coef(glm_c2_ldl_hi))
c3_ldl_OR <- exp(coef(glm_c3_ldl_hi))

c0_ldl_OR 
c1_ldl_OR 
c2_ldl_OR 
c3_ldl_OR


# bp_sys_hi ~ clusters + covars + income 
glm_c0_bp_sys_hi <- glm(bp_sys_hi ~ Cluster_0 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c1_bp_sys_hi <- glm(bp_sys_hi ~ Cluster_1 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c2_bp_sys_hi <- glm(bp_sys_hi ~ Cluster_2 + age + sex + hh_income, data = data_Ob, family = binomial())
glm_c3_bp_sys_hi <- glm(bp_sys_hi ~ Cluster_3 + age + sex + hh_income, data = data_Ob, family = binomial())

summary(glm_c0_bp_sys_hi)
summary(glm_c1_bp_sys_hi)
summary(glm_c2_bp_sys_hi)
summary(glm_c3_bp_sys_hi)

c0_bp_OR <- exp(coef(glm_c0_bp_sys_hi))
c1_bp_OR <- exp(coef(glm_c1_bp_sys_hi))
c2_bp_OR <- exp(coef(glm_c2_bp_sys_hi))
c3_bp_OR <- exp(coef(glm_c3_bp_sys_hi))

c0_bp_OR 
c1_bp_OR 
c2_bp_OR 
c3_bp_OR

### heatmap with ORs 

# Define binary outcomes and clusters
outcomes <- c("bmi_hi", "hba1c_hi", "bp_sys_hi", "ldl_hi")
cluster_suffixes <- c("c0", "c1", "c2", "c3")  # Cluster suffixes for model identification

# Initialize a results dataframe
results <- data.frame(Outcome = character(), Cluster = character(), OddsRatio = numeric(), pValue = numeric(), stringsAsFactors = FALSE)

# Function to extract odds ratios and p-values
extract_odds_ratios <- function(outcome, cluster_suffix) {
  # Construct the model name
  model_name <- paste("glm", cluster_suffix, outcome, sep = "_")
  
  # Check if the model exists and if so, proceed with extraction
  if (exists(model_name, envir = .GlobalEnv)) {
    model <- get(model_name)
    coef_info <- summary(model)$coefficients
    
    # Extract odds ratio and p-value for the cluster
    cluster_label <- paste("Cluster", gsub("c", "", cluster_suffix), sep = "_")
    if (cluster_label %in% rownames(coef_info)) {
      odds_ratio <- exp(coef_info[cluster_label, "Estimate"])
      pval <- coef_info[cluster_label, "Pr(>|z|)"]
      return(data.frame(Outcome = outcome, Cluster = cluster_label, OddsRatio = odds_ratio, pValue = pval))
    } else {
      print(paste("Coefficient name", cluster_label, "not found in the model summary for", model_name))
    }
  } else {
    print(paste("Model does not exist:", model_name))
  }
  return(data.frame(Outcome = outcome, Cluster = cluster_suffix, OddsRatio = NA, pValue = 1))
}

# Populate the results dataframe using loops
for (outcome in outcomes) {
  for (cluster_suffix in cluster_suffixes) {
    results <- rbind(results, extract_odds_ratios(outcome, cluster_suffix))
  }
}

# Filter results to include only significant odds ratios (pValue < 0.05)
significant_results <- results[results$pValue < 0.05, ]
significant_results$Outcome <- factor(significant_results$Outcome,
                                      levels = outcomes)  # Set the levels in the correct order


#Desired order of outcomes from top to bottom
desired_order <- c("High LDL", "Hypertension", "High HbA1c", "High BMI")

# Map original outcome variable names to descriptive labels
outcome_labels <- c(
  "bmi_hi" = "High BMI",
  "hba1c_hi" = "High HbA1c",
  "bp_sys_hi" = "Hypertension",
  "ldl_hi" = "High LDL"
)

# Update the 'Outcome' column to have factor levels in the specific order
significant_results$Outcome <- factor(significant_results$Outcome,
                                      levels = names(outcome_labels)[match(desired_order, outcome_labels)],
                                      labels = desired_order)


# Create the heatmap
heatmap_plot <- ggplot(significant_results, aes(x = Cluster, y = Outcome, fill = OddsRatio)) +
  geom_tile(color = "white", linewidth = 0.1) + # Defines the tile borders
  geom_text(aes(label = sprintf("%.2f", OddsRatio)), color = "black", size = 4) +  # Add odds ratio numbers
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1,
                       na.value = "white", limits = c(0.5, 2)) +  # Adjust color scale based on expected OR range
  labs(title = "Heatmap of Significant Odds Ratios", x = "", y = "", fill = "Odds Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, face = "bold"),  # Center x labels under each column
        axis.text.y = element_text(hjust = 1, size = 12, face = "bold"),  # Ensure y labels are aligned to the left of the plot
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        axis.ticks.x = element_blank(),  # Remove x-axis ticks to clean up the plot
        axis.ticks.y = element_blank()) +  # Remove y-axis ticks for a clean look
  scale_x_discrete(position = "top") +  # Ensure x-axis labels are at the top of the plot
  coord_fixed(ratio = 1)  # Ensure that tiles are square

# Print plot directly to view interactively
print(heatmap_plot)

##############################################################################################

## bmi ~ clusters + covars + income 
lm_bmi_Cluster_0 <- lm(bmi ~ Cluster_0 + age + sex + hh_income, data = data_Ob) 
lm_bmi_Cluster_1 <- lm(bmi ~ Cluster_1 + age + sex + hh_income, data = data_Ob) 
lm_bmi_Cluster_2 <- lm(bmi ~ Cluster_2 + age + sex + hh_income, data = data_Ob) 
lm_bmi_Cluster_3 <- lm(bmi ~ Cluster_3 + age + sex + hh_income, data = data_Ob) 

# Directly specify cluster variables based on known model naming
cluster_variables <- c("Cluster_0", "Cluster_1", "Cluster_2", "Cluster_3")

# Apply the extraction function to each model
model_data <- do.call(rbind, lapply(seq_along(models), function(i) {
  model <- models[[i]]
  cluster_variable <- cluster_variables[i]
  extract_coef_data(model, cluster_variable)
}))

# Print the collected model data
print(model_data)

library(grid)  # Ensure this is loaded for unit functionality

# Create the forest plot with corrected lineheight
forestplot(labeltext = label_matrix,
           mean = model_data$Coefficient,
           lower = model_data$CI_Lower,
           upper = model_data$CI_Upper,
           title = "Forest Plot of Cluster Coefficients",
           xlab = "Effect Size and 95% Confidence Intervals",
           zero = 0,
           lineheight = unit(1.5, "lines"),  # Correctly set using unit
           clip = range(c(model_data$CI_Lower, model_data$CI_Upper), na.rm = TRUE),
           graphwidth = unit(2, "inches"))


# Adjust the levels of 'Cluster' to reverse the order
plot_data$Cluster <- factor(plot_data$Cluster, levels = rev(unique(plot_data$Cluster)))

# Create the forest plot with ggplot2 in reversed order
ggplot(plot_data, aes(x = Coefficient, y = Cluster)) +
  geom_point(size = 4, color = "blue") +  # Coefficient as a point
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, color = "blue") +  # Confidence intervals
  labs(title = "Forest Plot of Cluster Coefficients",
       x = "Coefficient and 95% Confidence Intervals",
       y = "") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  # Zero line for reference

# Display the plot
ggsave("forest_plot_reversed.png", width = 8, height = 6, dpi = 300)





