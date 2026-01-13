################################################################################
################## Homogeneous Fuel Analysis: Charcoal Data ####################
################################################################################
# Ref : Cornet et al. in prep
# Date: 2025-08-27
# Purpose: Find optimal lambda in boxcox logit tranform to normalize WL ratios
#          using pure herbaceous and pure woody fuel sedimentary charcoal data
#          Generate figures
################################################################################

# Load required packages
library(EnvStats)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(patchwork)

if (!dir.exists("Output")) {
  dir.create("Output")
}
################################################################################
# Transformation functions
################################################################################

# Box-Cox logit transform and its inverse
BoxCoxLogit = function(x,lambda=-.07) {
  if(lambda==0) log(x/(1-x)) else ((x/(1-x))^lambda-1)/lambda
} 
InvBoxCoxLogit = function(y,lambda = -.07) {
  if(lambda==0) exp(y)/(1+exp(y)) else (lambda*y+1)^(1/lambda)/(1+(lambda*y+1)^(1/lambda)) 
}

################################################################################
# Load and preprocess data
################################################################################

# Load charcoal data
charcoal_data = read.csv("Input/Charcoal_all_sites.csv", header=T, sep=";", dec=".",check.names = F)

# Sites to consider
site_list = c("DINA","ELI","ESB","NDA","DKL")
type_fuel = c(DINA = "Grass",ELI = "Grass", ESB = "Grass", NDA = "Grass", DKL = "Tree")

# Herbaceous and woody fuel selection
# Factorize site
# Filter out extreme WLRatio values (keep between 0.04 and 0.96)
# Create transformed WL variable : WL/-1-WL)
charcoal_data = charcoal_data %>% 
  filter(Site %in% site_list) %>% 
  filter(!(Site == "DKL" & Depth < 12)) %>%
  mutate(Site = factor(recode(Site, "DKL" = "DKL(b)"), levels = c("ELI","ESB","DINA","NDA","DKL(b)")))%>% 
  filter(WLRatio > .04 &  WLRatio < .96) %>%
  mutate(WL_1mWL = WLRatio/(1-WLRatio)) %>%
  mutate(Fuel = type_fuel[Site]) 

################################################################################
# Box-Cox transformation per site
################################################################################

# Define range of lambda to test
lambda_test =  seq(-2, 2, by = 0.01) # Lambda values to test

table_boxcox <- data.frame(lambda = lambda_test)

# Compute Box-Cox objective function for each site
for (site_name in levels(charcoal_data$Site)) {
  wl_vals <- charcoal_data %>%
    filter(Site == site_name) %>%
    pull(WL_1mWL)
  table_boxcox[[site_name]] <- boxcox(wl_vals, lambda = lambda_test)$objective
}

# Convert to long format for plotting
table_boxcox_long = table_boxcox %>%
  pivot_longer(-lambda, names_to = "Site", values_to = "objective") %>%
  mutate(Site = factor(Site,levels=c("ELI","ESB","DINA","NDA","DKL(b)")))

################################################################################
# Weighted lambda optimization
################################################################################

weights = c(DINA = 1,ELI = 1, ESB = 1, NDA = 1, "DKL(b)" = 4)

table_boxcox_long_weighted <- table_boxcox_long %>%
  mutate(weight = weights[table_boxcox_long$Site]) %>%
  group_by(lambda) %>%
  summarise(weighted_mean = weighted.mean(objective, w = weight), .groups = "drop")

lambda_optim_weighted <- table_boxcox_long_weighted %>%
  filter(weighted_mean == max(weighted_mean)) %>%
  pull(lambda)# Plot

################################################################################
# Define site-specific color palette
################################################################################

palette_all_sites <- c(
  "ELI"  = "#FFFB03",  
  "ESB"  = "#FFD903",  
  "DINA" = "#FFB703",  
  "NDA"  = "#FF9603",  
  "DKL(b)"  = "#006400"   
)

linewidth_all_sites <- c(
  "ELI"    = 0.7,
  "ESB"    = 0.7,
  "DINA"   = 0.7,
  "NDA"    = 0.7,
  "DKL(b)" = 1.4
)

################################################################################
# Plot: Box-Cox results and optimal lambda
################################################################################

plot_lambda_boxcox = ggplot(data = table_boxcox_long) +
  geom_line(aes(x = lambda, y = objective, color = Site, linewidth = Site)) +
  geom_vline(xintercept = lambda_optim_weighted, linetype = "dashed", linewidth = 1, color = "black") +
  annotate("text", x = lambda_optim_weighted + 0.05, y = 1.05,
           label =paste0("lambda == ", round(lambda_optim_weighted, 2)),parse=T,
           hjust = 0, vjust = 1.2, size = 4.5, fontface = "italic") +
  scale_color_manual(values = palette_all_sites) +
    scale_linewidth_manual(values = linewidth_all_sites)+
  guides(linewidth = "none") +  
  labs(
    x = expression("Box-Cox parameter " (lambda)),
    y = "Log-likelihood (objective function)",
    color = "Study sites") +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )



plot_lambda_boxcox
ggsave(filename = "Output/Figure S1 Optimal Lambda for BCLogit.png", plot = plot_lambda_boxcox, width = 6, height = 6, dpi = 600)

################################################################################
# Density plots: Original vs Transformed WL Ratio
################################################################################

# Plot A: WLRatio (original scale)
plot_char_treeherb_1 <- ggplot(data = charcoal_data) +
  geom_density(aes(x = WLRatio, color = Site, fill = Site, linewidth = Site),
               alpha = 0.4, position = "identity") +
  labs(x = "WL ratio", y = "Density", color = "Study sites", fill = "Study sites") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = palette_all_sites) +
  scale_fill_manual(values = palette_all_sites) +
  scale_linewidth_manual(values = linewidth_all_sites)+
  guides(linewidth = "none") +  
  theme_pubr() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),  # ← Titre de légende en gras
    legend.text  = element_text(size = 11),
    plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  )

# Plot B: Box_Cox logit transformed WLRatio
plot_char_treeherb_2 <- ggplot(data = charcoal_data) +
  geom_density(aes(x = BoxCoxLogit(WLRatio, lambda_optim_weighted), color = Site, fill = Site, linewidth = Site),
               alpha = 0.4, position = "identity") +
  labs(x = "Box-Cox logit-transformed WL ratio", 
       y = "Density", color = "Study sites", fill = "Study sites") +
  scale_x_continuous(limits = c(-4.5, 4.5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = palette_all_sites) +
  scale_fill_manual(values = palette_all_sites) +
  scale_linewidth_manual(values = linewidth_all_sites)+
  guides(linewidth = "none") +  
  theme_pubr() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),  # ← Titre de légende en gras
    legend.text  = element_text(size = 11),
    legend.position = "none",
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  )


# Combine both plots vertically
plot_treeherb <- plot_treeherb <- plot_char_treeherb_1 / plot_char_treeherb_2 +   
  plot_layout(ncol = 1, nrow = 2, guides = "collect") +      
  plot_annotation(tag_levels = "A") &                        
  theme(legend.position = "bottom") 


print(plot_treeherb)
ggsave(filename = "Output/Figure 2 Charcoal_distribution_homogeneous_fuel.png", plot = plot_treeherb, width = 6, height = 6, dpi = 600)


################################################################################
# Descriptive statistics by site for original and transformed WL Ratio and transformd
# (mode,mean, std, median, Q1 et Q3)
################################################################################

# Mode function
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Function to summarize statistics
get_stats <- function(data, variable, transformation_label) {
    data %>%
      group_by(Fuel, Site) %>%
      summarise(
        mu = mean({{ variable }}, na.rm = TRUE),
        mediane = median({{ variable }}, na.rm = TRUE),
        Q1 = quantile({{ variable }}, 0.25, na.rm = TRUE),
        Q3 = quantile({{ variable }}, 0.75, na.rm = TRUE),
        mode = get_mode(round({{ variable }}, 3)),  
        sigma = sd({{ variable }}, na.rm = TRUE)
      ) %>%
      mutate(transformation = transformation_label)
  }

# Compute stats for original and transformed data
stats_original <- get_stats(charcoal_data, WLRatio, "Original")
stats_transformed <- get_stats(charcoal_data, BoxCoxLogit(WLRatio), "Transformed")

# Combined final table
combined_stats <- bind_rows(stats_original, stats_transformed) %>%
  arrange(Fuel, transformation, Site)%>%
  mutate(across(where(is.numeric), ~ round(., 3)))

print(combined_stats)

write.csv(combined_stats, "Output/Synth_homogeneous_fuel_data.csv", row.names = FALSE)






