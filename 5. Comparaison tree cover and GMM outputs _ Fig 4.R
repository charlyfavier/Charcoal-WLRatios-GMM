################################################################################
###### Comparison of Treecover with modelled proportions of woody charcoals ####
################################################################################
# Ref : Cornet et al. in prep
# Date: 2025-08-27
# Purpose: Compute Vegetation Class areas and Fire Occurence per vegetation 
#          class in all sites 
#          Computes Tree Cover from two remote sensing data inside radii from
#          1000 m to 5000 m around lake in Savannas, Savannas + Agricultural
#          Land, Area burned at least once
################################################################################

# Load required packages
library(ggplot2) 
library(ggpubr)
library(corrplot) 
library(ggrepel)
library(ggpubr) 
library(dendextend)
library(dplyr)
library(tidyverse)

################################################################################
# Load data
################################################################################

# Load GMM outputs 
GMM_outputs_30yrs <- read.csv("Output/Bayesian_GMM_outputs_30_yrs.csv",header=T, dec=".",check.names = F)
# Load TC estimates
TC_data_brut <- read.csv("Output/TC vegclass disks.csv",header=T,  dec=".",check.names = F)

GMM_outputs_30yrs = GMM_outputs_30yrs %>%
  rename(
    P_woody      = `theta[2]50%`,
    P_woody_low  = `theta[2]25%`,
    P_woody_high = `theta[2]75%`
  )

TC_data_brut =  TC_data_brut %>%
  column_to_rownames(var = "Site")
################################################################################
# Parameters for visulalization
################################################################################

site_colors <- c("DKL" = "#70e000", "FBA" = "#3f37c9", "NGO" = "#d62828", "ELI" = "#f8961e", "GRO" = "#00f5d4")

x_max = max(GMM_outputs_30yrs$P_woody_high)
x_min = min(GMM_outputs_30yrs$P_woody_low)


################################################################################
# Check of colinearity
################################################################################
# Calculation of the correlation matrix
# cor_matrix <- cor(TC_data[,10:length(TC_data)], use = "pairwise.complete.obs") # voir avec Charly
cor_matrix <- cor(TC_data_brut, use = "pairwise.complete.obs")

# Displaying the correlation matrix with a colour map
corrplot(cor_matrix, method = "circle", type = "upper",  order = "hclust",
         col = colorRampPalette(c("#00b4d8", "white", "#a7c957"))(200), 
         tl.cex = 0.6)


# Compute dendrogram
hc <- hclust(as.dist(1 - abs(cor_matrix)))  # distance = 1 - |corr|

png("Output/Figure S16 dendrogram_TC.png", width = 4000, height = 3500, res = 600)

par(mar = c(10, 4, 4, 2))  
as.dendrogram(hc) %>%
  set("branches_k_color", k = 6) %>%   # colorier en 4 groupes
  set("labels_cex", 0.7) %>%           # réduire la taille des labels
  plot(main = "Hierarchical Clustering of Variables")

dev.off()

# FILTER DATA
kept_variables = c("TC1m_5000_burned_weighted_mean", "TC1m_5000_sav","TC1m_1000_savagr","TC10m_5000_sav","TC10m_5000_savagr","TC10m_5000_burned_weighted_mean")

data_filtered <- GMM_outputs_30yrs %>% 
  dplyr::select(Site, P_woody,P_woody_low,P_woody_high) %>% 
  left_join(TC_data_brut  %>%
              dplyr::select(all_of(kept_variables)) %>%  
              rownames_to_column(var = "Site"),
            by = join_by(Site)) %>%
  mutate(color = ifelse(
    Site %in% names(site_colors),                # Si le site est dans la liste
    site_colors[Site],                           # Prend la couleur correspondante
    "#a2999e"                                    # Sinon couleur par défaut
  ))
data_filtered <- data_filtered[order(data_filtered$color == "#a2999e", decreasing = TRUE), ]   
cor_matrix <- cor(data_filtered[,kept_variables],use = "pairwise.complete.obs")

# Display the correlation matrix with a colour map
corrplot(cor_matrix, method = "circle", type = "upper", 
         col = colorRampPalette(c("#00b4d8", "white", "#a7c957"))(200), 
         tl.cex = 0.8)




# Loop on kept variables
for (colname in kept_variables) {
  cat(colname)
  model <- lm(get(colname) ~ P_woody+0, data = data_filtered)
  model_summary <- summary(model)
  p_value <- coef(model_summary)[nrow(coef(model_summary)), 4]
  r_squared = model_summary$adj.r.squared

  p <- ggplot(data_filtered, aes(x = P_woody, y = get(colname), label = Site)) +
    
    geom_segment(aes(x = P_woody_low, xend = P_woody_high,
                     y = get(colname), yend = get(colname), color = color), size = 0.5) +
    geom_point(aes(color = color), size = 3) +
    # Linear regression
    geom_smooth(method = "lm",formula = y ~ 0 + x, se = FALSE, color = "red") +
    
    # Labels
    geom_text_repel(size = 3.5, max.overlaps = 20) +
    xlim(0,1)+
  # Statistics
    annotate("text",
           x = max(data_filtered$P_woody, na.rm = TRUE) * 0.8,
           y = max(data_filtered[[colname]], na.rm = TRUE) * 0.2,
           label = paste("R² =", round(r_squared, 4), "\np-value =", signif(p_value, 3)),
           color = "black", size = 4) +
    
    # Formatting
    labs(title = paste("Linear regression : P_woody vs", colname),
         x = "Propotion of woody charcoal",
         y = "TC") +
    theme_pubr() +
    theme(legend.position = "none") +
    scale_color_identity()
  
  ggsave(paste0("Output/Figure Comparison P_ligneous with ",colname,".png"), plot = p, width = 8, height = 8, dpi = 600)
}
 

### Figure 4 : regression of TC10m_5000_sav~P_ligneous

colname = "TC10m_5000_sav"
model <- lm(get(colname) ~ P_woody+0, data = data_filtered)
model_summary <- summary(model)
p_value <- coef(model_summary)[nrow(coef(model_summary)), 4]
r_squared = model_summary$adj.r.squared
# Filtering based on significant p-values
plot_Pwoody_TC10m_5000_sav <- ggplot(data_filtered) +
    
    geom_segment(aes(x = P_woody_low, xend = P_woody_high,
                     y = get(colname), yend = get(colname), color = color), size = 0.5) +
    geom_point(aes(x = P_woody, y = get(colname), color = color), size = 3) +
    # Linear regression
    geom_smooth(aes(x = P_woody, y = get(colname)),method = "lm",formula = y ~ 0 + x, se = F, color = "red") +
    
    scale_color_identity()+
  
    # Labels
    geom_text_repel(aes(x = P_woody, y = get(colname),label = Site), size = 5, max.overlaps = 20,box.padding = 0.4) +
    xlim(0,1)+
    # Statistics
    annotate("text",
             x = max(data_filtered$P_woody, na.rm = TRUE) * 0.8,
             y = max(data_filtered[[colname]], na.rm = TRUE) * 0.2,
             label = paste("R² =", round(r_squared, 4), "\np-value =", signif(p_value, 3)),
             color = "black",size=6) +
    
    # Formatting
    labs(title = "",
         x = "Propotion of woody charcoal",
         y = "JRC TC estimates in savannas") +
    theme_pubr() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 20),  
      axis.title.y = element_text(size = 20),  
      axis.text = element_text(size = 15)     
    )

plot_Pwoody_TC10m_5000_sav

ggsave("Output/Figure 4 Comparison P_ligneous-TC10m_5000_sav.png", plot = plot_Pwoody_TC10m_5000_sav, width = 8, height = 8, dpi = 600)




