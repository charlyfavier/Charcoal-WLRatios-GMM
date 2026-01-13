################################################################################
###### GMM model for sedimentary charcoal data along Sele sediment core ########
################################################################################
# Reference : Cornet et al. (in prep)
# Date      : 2025-08-27
# Purpose   : Estimate proportions of grass vs. tree charcoal using a 
#             Bayesian Gaussian Mixture Model (GMM) with constraints.
################################################################################

#===============================================================================
# Load Required Libraries
#===============================================================================

library(progress)     
library(dplyr)        
library(ggplot2)      
library(ggpubr)       
library(locfit)      
library(MASS)
source("Functions/GMM analysis function.R")

#===============================================================================
# Load Data
#===============================================================================

# Load charcoal particle data
charcoal_data = read.csv(
  "Input/Sele_charcoal.csv",
  header = T,
  sep = ";",
  dec = ".",
  check.names = F
)

# Load age-depth model
age_model <- read.csv(
  "Input/Sele_2015_age_depth.csv",
  header = T,
  sep = ";",
  dec = ".",
  check.names = F
)

#===============================================================================
# Compute Age for each sample
#===============================================================================

# Merge age-depth information into charcoal data
charcoal_data <- charcoal_data %>%
  left_join(age_model, by="Depth")

#===============================================================================
# Compute Concentrations and WL Ratio statistics per sample
#===============================================================================

# Compute Summary Statistics per Depth Level
#Sample volume
samples_volume = 0.5

charcoal_per_level_Sele = charcoal_data %>% 
  group_by(Depth) %>% 
  summarise(Age = mean(Age),
            Counts = n(), 
            Concentration = n()/samples_volume, 
            AreaConcentration = sum(ProjArea)/samples_volume, 
            meanWL = mean(WLRatio),
            Prop.WL.gt.0.5 = sum(WLRatio > 0.5) / n())

# Build Complete Depth Series
p_fuels = data.frame(Depth = seq(min(charcoal_data$Depth),max(charcoal_data$Depth),1)) %>%
  left_join(age_model, by="Depth") %>%
  left_join(charcoal_per_level_Sele %>% dplyr::select(-Age), by="Depth") %>%
  mutate(Counts = ifelse(is.na(Counts),0,Counts),
         Concentration = ifelse(is.na(Concentration),0,Concentration),
         AreaConcentration = ifelse(is.na(AreaConcentration),0,AreaConcentration))


# RUN model on slices 
gmm_output = analyze_WLcore_GMM(WLRatio_series  = charcoal_data$WLRatio,
                                depth_series = charcoal_data$Depth,
                                age_serie = charcoal_data$Age,
                                stan_file       = "Stan/Herb_Woody_Charcoal_GMM.stan",
                                n_bins        = 30,
                                sampling_type = "regular", # "adaptive"
                                graph.plot      = TRUE,
                                N_samples     = 1000,
                                progress      = TRUE,
                                chains          = 5,
                                parallel_chains = 5,
                                iter_warmup     = 1000,
                                iter_sampling   = 2000,
                                refresh         = 0,
                                show_messages   = FALSE,
                                mc.cores        = 5,  
                                check_toolchain = TRUE)

all_quantiles = gmm_output$parameter_quantiles

#===============================================================================
# Visualization of Posterior Parameters
#===============================================================================

#-----------------------------------------
# Lambda parameter
#-----------------------------------------
plot(ggplot(all_quantiles %>% dplyr::filter(param == "lambda"),aes(x=depth,group=param))+
       geom_point(aes(y=q0.5))+
       geom_smooth(aes(y=q0.5),span=0.1)+
       geom_linerange(aes(ymin = q0.05, ymax = q0.95), linewidth = 1) +
       ylab("lambda")+
       theme_pubclean())


#-----------------------------------------
# Mean parameters (muH, muW)
#-----------------------------------------
plot(ggplot(all_quantiles %>% dplyr::filter(param %in% c("muH","muW")),aes(x=depth,group=param,color=param))+
       geom_point(aes(y=q0.5))+
       geom_smooth(aes(y=q0.5),span=0.1)+
       geom_linerange(aes(ymin = q0.05, ymax = q0.95), linewidth = 1) +
       scale_color_manual(
         values = c("muH" = "#ffb703", "muW" = "#006400")
       ) +
       scale_fill_manual(
         values = c("muH" = "#ffb703", "muW" = "#006400")
       ) +
       ylab("mu")+
       theme_pubclean())

#-----------------------------------------
# Standard deviation parameters (sigmaH, sigmaW)
#-----------------------------------------
plot(ggplot(all_quantiles %>% dplyr::filter(param %in% c("sigmaH","sigmaW")),aes(x=depth,group=param,color=param))+
       geom_point(aes(y=q0.5))+
       geom_linerange(aes(ymin = q0.05, ymax = q0.95), linewidth = 1) +
       geom_smooth(aes(y=q0.5, ymax = q0.95),span=0.1, linewidth = 1) +
       scale_color_manual(
         values = c("sigmaH" = "#ffb703", "sigmaW" = "#006400")
       ) +
       scale_fill_manual(
         values = c("sigmaH" = "#ffb703", "sigmaW" = "#006400")
       ) +
       ylab("sigma")+
       theme_pubclean())

# RUN INTERPOLATION OF THE PROPORTION OF CHARCOAL WITH WLRatio > 0.5 and 
# APPLY LOCAL TRANSFORM TO PORPOTION OF WOODY CHARCOAL
pp= analyze_WLcore_proportions(
  gmm_output,
  p_fuels$Depth,
  p_fuels$Age,
  nn_locfit = 0.03,
  h_locfit  = 20,
  progress  = TRUE
)

# Combine sample characteristics with model outputs
p_fuels = p_fuels %>% left_join(pp , by = "Age" )

# Quick diagnostic plot: observed per-level proportion vs reconstructed P_H
ggplot() +
  geom_point(data = charcoal_per_level_Sele, aes(x = Age, y = Prop.WL.gt.0.5)) +
  geom_line(data = p_fuels, aes(x = Age, y = P_W), color = "red") +
  labs(
    x = "Age",
    y = "Observed Prop WL > 0.5"
  )+
  scale_x_reverse()+
  theme_pubclean()

# Save p_fuels table with computed probabilities and summary statistics
write.csv(p_fuels, "Output/Charcoal_statistics_Sele.csv", row.names = FALSE)


#===============================================================================
# Comparison Figures: WL-based proxies vs GMM-derived proportions
#===============================================================================

# Load previously saved p_fuel ans GMM summary table 
table_GMM_results  =read.csv("Output/Bayesian_GMM_outputs_30_yrs.csv",h=T,check.names = FALSE)
p_fuels  =read.csv("Output/Charcoal_statistics_Sele.csv",h=T,check.names = FALSE)

# Plot 1: Proportion WL > 0.5 (empirical) vs GMM ligneous proportion (theta[2]50%)
plot_PropWLsup.PLigneous = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",linewidth=1,color = "gray") +
  geom_point(data = p_fuels,aes(x = Prop.WL.gt.0.5, y = 1-P_H,color= Concentration), size = 3) +
  geom_point(data = table_GMM_results, aes(x = Prop.WL.gt.0.5, y = `theta[2]50%`), color = "#0B1F3A", size = 3) +
  geom_smooth(data = p_fuels,aes(x = Prop.WL.gt.0.5, y = 1-P_H),se=F,method=rlm, linetype = "dashed",linewidth=1) +
  labs(x = "Proportion of WL Ratio > 0.5", y = "Proportion of ligneous-derived charcoal", title = "") +
  xlim(0,1)+
  ylim(0,1)+
  scale_color_continuous(low="#fee4e2",high="#F8766D", guide = "none") +
  theme_pubr(base_size = 14)

# Plot 2: Mean WL Ratio vs GMM ligneous proportion
plot_meanWL.PLigneous = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",linewidth=1,color = "gray") +
  geom_point(data = p_fuels,aes(x = meanWL, y = 1-P_H,color=Concentration),  size = 3) +
  geom_point(data = table_GMM_results, aes(x = meanWL, y = `theta[2]50%` ), color = "#0B1F3A", size = 3) +
  geom_smooth(data = p_fuels,aes(x = meanWL, y = 1-P_H),se=F,method=rlm, linetype = "dashed",linewidth=1) +
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "Mean WL Ratio", y = "Proportion of ligneous-derived charcoal", title = "") +
  scale_color_continuous(low="#fee4e2",high="#F8766D",guide="none") +
  theme_pubr(base_size = 14)

# Plot 3: Proportion WL > 0.5 vs Mean WL Ratio (proxy agreement)
plot_PropWLsup.meanWL = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",linewidth=1,color = "gray") +
  geom_point(data = p_fuels,aes(x = Prop.WL.gt.0.5, y = meanWL,color=Concentration),  size = 3) +
  geom_point(data = table_GMM_results, aes(x = Prop.WL.gt.0.5, y = meanWL ), color = "#0B1F3A", size = 3) +
  geom_smooth(data = p_fuels,aes(x = Prop.WL.gt.0.5, y = meanWL),se=F,method=rlm, linetype = "dashed",linewidth=1) +
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "Proportion of WL Ratio > 0.5", y = "Mean WL Ratio", title = "") +
  scale_color_continuous(low="#fee4e2",high="#F8766D",guide="none") +
  theme_pubr(base_size = 14)

# Combine the three plots vertically
plot_comparison <- plot_PropWLsup.PLigneous / 
  plot_meanWL.PLigneous / 
  plot_PropWLsup.meanWL

# Save final comparison figure to disk
ggsave("Output/Figure S13 Comparison proxies lambdamod.pdf", plot_comparison, width = 5, height = 15, units = "in")
