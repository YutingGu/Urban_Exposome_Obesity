library(readr)
library(dplyr)
library(tidyr)

setwd("/rds/general/project/hda_23-24/live/CompEpi/group1/Scripts")


# load dataset
cluster_id <- read_csv("../Dataset/cluster_id.csv")
cluster_id$eid <- as.character(cluster_id$eid)
cleaned_ukb_with_protein <- readRDS("../Dataset/cleaned_ukb_with_protein.rds")

# join the label with the 
data <- cleaned_ukb_with_protein %>% left_join(cluster_id, join_by(eid==eid))

# convert the cluster_id to one-hot encode
data$Cluster_0 <- ifelse(data$Cluster==0, 1, 0)
data$Cluster_1 <- ifelse(data$Cluster==1, 1, 0)
data$Cluster_2 <- ifelse(data$Cluster==2, 1, 0)
data$Cluster_3 <- ifelse(data$Cluster==3, 1, 0)

table(data$Cluster)

################################################################################
exposome_list <- names(data)[2:16]
data_exposome = data[c(exposome_list,"Cluster_0","Cluster_1","Cluster_2","Cluster_3","sex")]


# "no2_airpol"       "no_airpol"        "pm10"             "pm2_5"            "pm2_5_ab"        
# "trf_int_rd"       "inv_dist_rd"      "trf_int_maj_rd"   "inv_dist_maj_rd"  "grnspc_buf1000"  
# "dom_grdn_buf1000" "water_buf1000"    "nat_env_buf1000"  "dis_coast"        "avg_sndpol"      

data$Cluster = as.factor(data$Cluster)

### significant ones
ggplot(data, aes(x = log(inv_dist_maj_rd), fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  ggtitle("Inverse distance to the nearest major road")

  
ggplot(data, aes(x = pm2_5_ab, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  ggtitle("Particulate matter air pollution PM2.5 Absorbance")


ggplot(data, aes(x = pm2_5, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  ggtitle("Particulate matter air pollution PM2.5")


ggplot(data, aes(x = pm10, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  ggtitle("Particulate matter air pollution PM10")


ggplot(data, aes(x = avg_sndpol, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()  +
  ggtitle("Average 24-hour sound level of noise pollution")

ggplot(data, aes(x = log(inv_dist_rd), fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()  +
  ggtitle("Inverse distance to the nearest road")


### others
ggplot(data, aes(x = no2_airpol, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(data, aes(x = no_airpol, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()
  
ggplot(data, aes(x = log(trf_int_rd), fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(data, aes(x = trf_int_maj_rd, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(data, aes(x = grnspc_buf1000, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(data, aes(x = dom_grdn_buf1000, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(data, aes(x = log(water_buf1000), fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(data, aes(x = nat_env_buf1000, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(data, aes(x = dis_coast, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()


################################################################################
ggplot(data, aes(x = log(inv_dist_maj_rd), y = pm2_5_ab, fill = Cluster, colour = Cluster)) +
  geom_point() +
  labs(x = "log(Inverse distance to the nearest major road)", 
       y = "Particulate matter air pollution PM2.5 Absorbance") +
  theme_bw()

ggplot(data, aes(x = log(inv_dist_maj_rd), y = avg_sndpol, fill = Cluster, colour = Cluster)) +
  geom_point() +
  labs(x = "log(Inverse distance to the nearest major road)", 
       y = "Average 24-hour sound level of noise pollution") +
  theme_bw()

ggplot(data, aes(x = pm2_5_ab, y = pm10, fill = Cluster, colour = Cluster)) +
  geom_point() +
  labs(x = "Particulate matter air pollution PM2.5 Absorbance", 
       y = "Particulate matter air pollution PM10") +
  theme_bw()



################################################################################
model_result_unadj <- data.frame()
for (cid in c("Cluster_0","Cluster_1","Cluster_2","Cluster_3")){
  for (exposome in exposome_list){
    formula <- as.formula(paste(cid, '~', exposome))
    model <- glm(formula, data = data_exposome, family = binomial(link = "logit"))
    model_summary <- summary(model)
    exposome_coef <- model_summary$coefficients[exposome,]
    OR_CI = confint(model)
    
    OR <- exp(exposome_coef['Estimate'])
    OR_lb <- exp(OR_CI[exposome,1])
    OR_ub <- exp(OR_CI[exposome,2])
    p_value <- exposome_coef['Pr(>|z|)']
    
    new_row <- c(cid, exposome, OR, OR_lb, OR_ub, p_value)  
    model_result_unadj <- rbind(model_result_unadj, new_row)
  }
}

names(model_result_unadj) <- c('cluster', 'exposome', 'OR', 'OR_0.025', 'OR_0.975', 'p_value')


# save the results
saveRDS(model_result_unadj, '../Results/01_univariate_results/01_Univariate_Analysis_Exposome_Results_unadj.rds')
model_result_unadj <- readRDS('../Results/01_univariate_results/01_Univariate_Analysis_Exposome_Results_unadj.rds')

################################################################################
# Plot the valcano plot
Volcano = function(results, annot = NULL, thr = 0.05) {
  par(mar = c(4.5, 4.5, 1, 1))
  plot(results$OR, -log10(results$p_value), pch = 19,
       las = 1, cex = 0.5, xlab = expression(beta),
       ylab = expression(-log[10](p[value])), 
       col = ifelse(p.adjust(results$p_value,method = "BH") < 0.05, yes = "tomato",no = "darkgrey"))
  if (!is.null(annot)) {
    annot = annot[rownames(results), ]
    text(results$OR, -log10(results$p_value), pos = 3,
         offset = 0.2, cex = 0.5, labels = ifelse(results$p_value < thr, yes = annot$Symbol, no = ""))
  }
  abline(v = 1, lty = 3)
  abline(h = -log10(0.05/nrow(results)), lty = 2,col = "darkred")
  legend("topleft", col = c("darkred", "tomato","darkgrey"), lty = c(2, NA, NA), pch = c(NA,19, 19), 
         cex = 0.7, legend = c("Bonferroni threshold at 0.05","FDR significant hits", "Not significant"))
}

model_result_unadj$OR = as.numeric(model_result_unadj$OR)
model_result_unadj$OR_0.025 = as.numeric(model_result_unadj$OR_0.025)
model_result_unadj$OR_0.975 = as.numeric(model_result_unadj$OR_0.975)
model_result_unadj$p_value = as.numeric(model_result_unadj$p_value)

# plot volcano plot for each cluster
Volcano(model_result_unadj[model_result_unadj$cluster == 'Cluster_0',])
Volcano(model_result_unadj[model_result_unadj$cluster == 'Cluster_1',])
Volcano(model_result_unadj[model_result_unadj$cluster == 'Cluster_2',])
Volcano(model_result_unadj[model_result_unadj$cluster == 'Cluster_3',])

# Bonferroni (FWER) Correction
# 0.0033
0.05/15
model_result_unadj$p_value_significance <- model_result_unadj$p_value < (0.05/15)

# save the results to CSV
write_csv(model_result_unadj, '../Results/01_univariate_results/01_Univariate_Analysis_Exposome_Results_unadj.csv')

# subset data of cluster
result_0 <- model_result_unadj[model_result_unadj$cluster == 'Cluster_0',]
result_1 <- model_result_unadj[model_result_unadj$cluster == 'Cluster_1',]
result_2 <- model_result_unadj[model_result_unadj$cluster == 'Cluster_2',]
result_3 <- model_result_unadj[model_result_unadj$cluster == 'Cluster_3',]


# 








