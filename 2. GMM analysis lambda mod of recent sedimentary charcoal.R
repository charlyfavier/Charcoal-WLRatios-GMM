################################################################################
###### GMM model for sedimentary charcoal < 30 years in sediment cores #########
################################################################################
# Ref : Cornet et al. in prep
# Date: 2025-08-27
# Purpose: Estimate proportions from herbaceous vs woody charcoal from  
# Stan Bayesian Gaussian Mixture Model (GMM) with additional constraints
################################################################################


# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggpattern)
source("Functions/GMM analysis function.R")

################################################################################
# Load charcoal and age model data, keep records from 1990, filter out extreme values
################################################################################
# Charcoal data
charcoal_data = read.csv("Input/Charcoal_all_sites.csv", header=T, sep=";", dec=".",check.names = F)
# Age models
age_models = read.csv("Input/Agemodels_all_sites.csv", header=T, sep=";", dec=".",check.names = F)
# Associate both tables and filter
charcoal_data = charcoal_data %>% 
  left_join(age_models,by=join_by(Site,Depth)) %>%
  mutate(Site = ifelse(Site=="DKL",ifelse(Depth<=10,"DKL","DKL(b)"),Site))%>%
  filter(WLRatio > .01 &  WLRatio < .99) %>% 
  filter(Age > 1990)

#Ensure equitable charcoal reparition in cas Volume changes in sequence
charcoal_data <- charcoal_data %>% group_by(Site) %>%
  mutate(keep = runif(n())<min(Volume)/Volume)

write.csv(charcoal_data, "Output/Charcoal_data_kept.csv", row.names = FALSE)


charcoal_data <- charcoal_data %>% 
  filter(keep) %>%
  dplyr::select(-keep, Volume)



################################################################################
# Loop over sites to apply GMM
################################################################################
# Initialize result table
table_GMM_results <- data.frame()

# List of sites
sites <- unique(charcoal_data$Site)

# Main Loop over Sites
for (i in seq_along(sites)) {
  study_site <- sites[i]
  message(sprintf("[%d/%d] Analyse en cours pour le site : %s", i, length(sites), study_site))
  
  ################################################################################
  # GMM model fit
  ################################################################################
  # Prepare data for Stan GMM model
  WLseries <- charcoal_data %>% filter(Site == study_site) %>% pull(WLRatio)
  
  # GMM model fit
  line_param_studysite = analyze_WLseries_GMM(WLseries, study_site = study_site)
  
  # Add site results to global table
  table_GMM_results <- rbind(table_GMM_results, c(Site = study_site,line_param_studysite))
}

table_GMM_results = table_GMM_results %>%
  rename_with(~ c("Site",names(line_param_studysite)))%>%  
  mutate(across(-Site, as.numeric)) %>% 
  left_join(charcoal_data %>% group_by(Site) %>% summarise(Prop.WL.gt.0.5 = sum(WLRatio > 0.5) / n(), meanWL = mean(WLRatio)) , by = "Site")

# Save results
write.csv(table_GMM_results, "Output/Bayesian_GMM_outputs_30_yrs.csv", row.names = FALSE)

################################################################################
# Graphs and Visualization of GMM Results
################################################################################
table_GMM_results = read.csv("Output/Bayesian_GMM_outputs_30_yrs.csv", h=T, check.names = F)
# Plot the lambda mode median and and interquartile ranges (25%-75%) 
plot_lambda = ggplot(table_GMM_results) +
  geom_segment(
    aes(x = Site,xend = Site,y = `lambda25%`,yend = `lambda75%`),
    color = "gray",   
    linewidth = 1.2
  ) +
  geom_hline(yintercept = -.07,linetype="dashed")+
  geom_point(aes(x = Site, y = `lambda50%`),size = 3)+
  labs(x = "Site",y = "Box-Cox parameter",color=NULL) +
  theme_pubr(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1),
        axis.text.x = element_text(angle = 60, hjust = 1))

# Plot the WLratio mode median and and interquartile ranges (25%-75%) 
# for woody vs herbaceous charcoal distributions for each site 
plot_modes_WL = ggplot(table_GMM_results) +
  geom_segment(
    aes(x = `WL[1]50%`,xend = `WL[1]50%`,y = `WL[2]25%`,yend = `WL[2]75%`),
    color = "#006400",   
    linewidth = 1.2
  ) +
  geom_segment(
    aes(x = `WL[1]25%`,xend = `WL[1]75%`,y = `WL[2]50%`,yend = `WL[2]50%`),
    color = "#ffb703",   
    linewidth = 1.2
  ) +
  geom_point(aes(x = `WL[1]50%`, y = `WL[2]50%`),size = 3)+
  geom_text_repel(
    aes(x = `WL[1]50%`, y = `WL[2]50%`, label = Site),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.8
  ) +
  xlim(0.1,0.37)+
  ylim(0.53,0.80)+
  labs(x = "Herbaceous charcoal WL Ratio mode",y = "Woody charcoal WL Ratio mode",color=NULL) +
  theme_pubr(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1))

# Plot the BC logit WLratio sd median and and interquartile ranges (25%-75%) 
# for woody vs herbaceous charcoal distributions for each site 
plot_sd_WL = ggplot(table_GMM_results) +
  geom_segment(
    aes(x = `sigma[1]50%`,xend = `sigma[1]50%`,y = `sigma[2]25%`,yend = `sigma[2]75%`),
    color = "#006400",   
    linewidth = 1.2
  ) +
  geom_segment(
    aes(x = `sigma[1]25%`,xend = `sigma[1]75%`,y = `sigma[2]50%`,yend = `sigma[2]50%`),
    color = "#ffb703",   
    linewidth = 1.2
  ) +
  geom_point(aes(x = `sigma[1]50%`, y = `sigma[2]50%`),size = 3)+
  geom_text_repel(
    aes(x = `sigma[1]50%`, y = `sigma[2]50%`, label = Site),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.8
  ) +
  xlim(0.5,1.15)+
  ylim(0.45,1.1)+
  labs(x = "Herbaceous charcoal SD",y = "Woody charcoal SD",color=NULL) +
  theme_pubr(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1))

# Plot the median and and interquartile ranges (25%-75%) 
# of sd BoxCoxlogit WLratio vs  mode WLratio vs 
# for herbaceous charcoal distributions for each site 
plot_sdmode_WL_herbs = ggplot(table_GMM_results) +
  geom_segment(
    aes(x = `WL[1]50%`,xend = `WL[1]50%`,y = `sigma[1]25%`,yend = `sigma[1]75%`),
    color = "#ffb703",   
    linewidth = 1.2
  ) +
  geom_segment(
    aes(x = `WL[1]25%`,xend = `WL[1]75%`,y = `sigma[1]50%`,yend = `sigma[1]50%`),
    color = "#ffb703",   
    linewidth = 1.2
  ) +
  geom_point(aes(x = `WL[1]50%`, y = `sigma[1]50%`),size = 3)+
  geom_text_repel(
    aes(x = `WL[1]50%`, y = `sigma[1]50%`, label = Site),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.8
  ) +
  xlim(0.1,0.37)+
  ylim(0.5,1.15)+
  labs(x = "Herbaceous charcoal WL mode",y = "Herbaceous charcoal SD",color=NULL) +
  theme_pubr(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1))

# Plot the median and and interquartile ranges (25%-75%) 
# of sd BoxCoxlogit WLratio vs  mode WLratio vs 
# for woody charcoal distributions for each site 
plot_sdmode_WL_trees = ggplot(table_GMM_results) +
  geom_segment(
    aes(x = `WL[2]50%`,xend = `WL[2]50%`,y = `sigma[2]25%`,yend = `sigma[2]75%`),
    color = "#006400",   
    linewidth = 1.2
  ) +
  geom_segment(
    aes(x = `WL[2]25%`,xend = `WL[2]75%`,y = `sigma[2]50%`,yend = `sigma[2]50%`),
    color = "#006400",   
    linewidth = 1.2
  ) +
  geom_point(aes(x = `WL[2]50%`, y = `sigma[2]50%`),size = 3)+
  geom_text_repel(
    aes(x = `WL[2]50%`, y = `sigma[2]50%`, label = Site),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.8
  ) +
  ylim(0.45,1.1)+
  xlim(0.53,0.80)+
  labs(x = "Woody charcoal WL mode",y = "Woody charcoal SD",color=NULL) +
  theme_pubr(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1))



# Prepare data for fuel proportion plots
data_to_plot =table_GMM_results %>% 
  arrange(`theta[2]50%`) %>%
  mutate(p_H_low = `theta[1]25%`,
         p_H_high = `theta[1]75%` - `theta[1]25%`,
         p_W = 1-`theta[1]75%`) %>% 
  dplyr::select(Site,p_H_low,p_H_high,p_W) %>%
  pivot_longer( cols = c(p_H_low,p_H_high,p_W),
                names_to = "Type",
                values_to = "Value")%>%
  left_join( data.frame(Type = c("p_H_low",  "p_H_high", "p_W"), Fuel = c("herbaceous", "herbaceous", "Tree")), by = "Type") %>%
  mutate(Type = factor(Type, levels = c("p_H_low",  "p_H_high", "p_W")),
         Site = factor(Site, levels = unique(Site)))



# Plot fuel proportions as cumulated bars per site with pattern fill for
# uncertainty regions
plot_proportions =  ggplot() +
  #geom_col(data = data_to_plot, aes(x = Site, y = Value, fill = Type), width = 0.6) +
  geom_col_pattern(
    data = data_to_plot,
    aes(
      x = Site,
      y = Value,
      pattern_fill = Type,
      pattern_fill2 = Type,
    ),
    pattern = "gradient" ,
    width = 0.6,
    pattern_angle = 0,
  )+ 
  geom_segment(data = table_GMM_results %>% 
                 arrange(`theta[2]50%`) %>%
                 mutate(Site = factor(Site, levels = unique(Site))), 
               aes(x = as.numeric(Site) - 0.6/2,
                   xend = as.numeric(Site) + 0.6/2,
                   y = `theta[2]50%`, yend = `theta[2]50%`),
               color = "black", linewidth = 1) +
  scale_pattern_fill_manual(
    values = c(
      "p_H_low"  = "#ffb703",
      "p_H_high" = "#006400",
      "p_W"      = "#006400"
    ),
    guide="none")+
  scale_pattern_fill2_manual(
    values = c(
      "p_H_low"  = "#ffb703",
      "p_H_high" = "#ffb703",
      "p_W"      = "#006400"
    ),
    guide="none")+
  labs(
    x = "Site",
    y = "Porportion of fuel",
    fill = "Fuel"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_rect(data = data.frame(x=c(0,0),y=c(0,0),Fuel=c("Herbaceous","Woody")),aes(xmin=x,xmax=x, ymin=y,ymax=y, fill = Fuel))+
  scale_fill_manual(
    values = c("Herbaceous" = "#ffb703","Woody" = "#006400"),name="Fuel type")+
  guides(fill = guide_legend(override.aes = list(fill = c("#ffb703", "#006400")))) +
  theme_pubr(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1),
  )


################################################################################
# Combine plots and save figures
################################################################################

# Combine plot_proportions and plot_modes_WL
figure3 <- plot_proportions/plot_modes_WL+
  plot_annotation(tag_levels = 'A')+ 
  plot_layout(heights = c(1, 1)) +
  theme_pubr(base_size = 14)

#figure3

# Save Figure 3
ggsave("Output/Figure 3 PWoody WLmodesd.pdf", figure3, width = 7, height = 14, units = "in")


# Combine plot_sd_WL, plot_sdmode_WL_herbs and plot_sdmode_WL_trees
figureS12 <- (plot_sd_WL | plot_lambda) / (plot_sdmode_WL_herbs|plot_sdmode_WL_trees)+
  plot_annotation(tag_levels = 'A')+ 
  plot_layout(heights = c(1, 1)) +
  theme_pubr(base_size = 14)

#figureS2

# Save Figure S2
ggsave("Output/Figure S12 other comparisons.pdf", figureS12, width = 14, height = 14, units = "in")

