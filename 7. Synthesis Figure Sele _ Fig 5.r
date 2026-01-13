################################################################################
################## Synthesis figure: pollen & charcoal SELE ####################
################################################################################
# Ref : Cornet et al. in prep
# Date: 2025-08-27
# Purpose: Comparison of the proportion of compact charcoal
#          produced by the model with pollen records over time
################################################################################

# Load required libraries
library(dplyr)
library(ggplot2)
library(MASS)
library(patchwork)
library(ggpubr)
library(tidyr)

################################################################################
# Define themes for figures
################################################################################

theme_without_xaxis = theme_pubr()+
   theme(
     axis.title.x = element_blank(),
     axis.text.x = element_blank(),
     axis.ticks.x = element_blank(),
     axis.line.x = element_blank(),
     axis.title.y = element_text(size = 9),
     axis.text.y  = element_text(size = 9),
     axis.title.y.left  = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
     axis.title.y.right = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
     plot.margin = unit(c(0.5, 1.5, 0, 0.5), "lines"),
     panel.background = element_rect(fill = "transparent", colour = NA), 
     plot.background  = element_rect(fill = "transparent", colour = NA)  
   )
 theme_with_xaxis = theme_pubr()+
   theme(
     axis.title.x = element_text(size = 10),
     axis.text.x  = element_text(size = 9),
     axis.ticks.x = element_line(),
     axis.line.x  = element_line(),
     axis.title.y = element_text(size = 10),
     axis.text.y  = element_text(size = 9),
     axis.title.y.left  = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
     axis.title.y.right = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
     plot.margin = unit(c(0.5, 1.5, 1, 0.5), "lines"),  
     panel.background = element_rect(fill = "transparent", colour = NA), 
     plot.background  = element_rect(fill = "transparent", colour = NA)  
   )
 
 # Function to adjust tag position in plots
 tag_side_theme <- function(side = c("left", "right"), hjust_offset = 0.3) {
   side <- match.arg(side)
   
   if (side == "left") {
     pos <- c(0, 1)               
     hjust_val <- 1 + hjust_offset        
   } else {
     pos <- c(1, 1)               
     hjust_val <- -hjust_offset           
    }
   
   theme(
     plot.tag.position = pos,
     plot.tag = element_text(hjust = hjust_val)
   )
 }
################################################################################
# Load pollen, charcoal data and age model data from SELE
################################################################################

# Pollen data
pollen_sele <- read.csv("Input/Pollen_Sele.csv", header = TRUE, sep = ";", dec = ".", check.names = FALSE)
# Age model for pollen core
 age_model_pollen <- read.csv("Input/Sele_pollen_age_depth.csv", header = TRUE, sep = ";", dec = ".", check.names = FALSE)
# Charcoal data
 charcoal_data <- read.csv("Input/Sele_charcoal.csv", header = TRUE, sep = ";", dec = ".", check.names = FALSE)
# Age model for charcoal core
 age_model_charcoal <- read.csv("Input/Sele_2015_age_depth.csv", header = TRUE, sep = ";", dec = ".", check.names = FALSE)
# output GMM model
 p_fuels <- read.csv("Output/Charcoal_statistics_Sele.csv", header = TRUE, dec = ".", check.names = FALSE)


# Set age range for plotting
 MaxAge = 3500
 MinAge = -70

# Apply age models and filter data within age range
charcoal_data <- charcoal_data %>% left_join(age_model_charcoal, by="Depth") %>% filter(Age <= MaxAge & !is.na(Age))
pollen_sele <- pollen_sele  %>% left_join(age_model_pollen, by="Depth") %>% filter(Age <= MaxAge & !is.na(Age))
p_fuels <- p_fuels %>% filter((Age <= MaxAge & !is.na(Age))) 

################################################################################
# Compute pollen group percentages
################################################################################

# Savanna pollen types
savanna_tree <- c(
  "Anacardiaceae undiff.", "Apocynaceae undiff.", "Caesalpiniaceae undiff.",
  "Aphania senegalensis", "Polygala", "Periplocaceae cf. Zacateza pedicellata", "Verbenaceae", "Cordia",
  "Erythrococca", "Lannea", "Croton", "Cissus", "Rubiaceae undiff.",
  "Malvaceae undiff.", "Flabellaria", "Fabaceae undiff.", "Capparis",
  "Afzelia", "Bombax", "Albizia-type", "Borassus",
  "Hymenocardia", "Vitex", "Glyphaea brevis", "Loranthaceae", "Cassia",
  "Labiatae", "Mimosaceae undiff.", "Desmodium", "Mitragyna inermis",
  "Convolvulaceae undiff."
)
# Forest pollen types
forest_tree <- c(
  "Adiantaceae", "Anthocleista", "Adenia", "Annonaceae undiff.",
  "Antidesma", "Baissea", "Burseraceae undiff.", "Landolphia",
  "Berlinia","Ficus", "Irvingia", "Triplochiton scleroxylon", "Cynometra",
  "Mansonia", "Pycnanthus", "Rothmannia", "Sapindaceae undiff.", "Saba",
  "Drypetes", "Diospyros", "Piptadeniastrum africanum", "Sapotaceae undiff.",
  "Cola", "Podocarpus", "Passifloraceae undiff.", "Ceiba pentandra",
  "Morelia senegalensis", "Canarium", "Dialium", "Mussaenda", "Parinari",
  "Crudia", "Araliaceae undiff.", "Cnestis", "Dalbergia", "Nauclea-type",
  "Pterocarpus", "Tetrapleura", "Pentadesma", "Syzygium", "Celtis", "Eugenia",
  "Gaertnera", "Meliaceae undiff.", "Morinda", "Menispermaceae undiff.",
  "Flacourtiaceae", "Lasiodiscus", "Ericaceae","Zanthoxylum","Musanga", "Macaranga","Elaeis guineensis","Daniellia","Holoptelea grandis","Mallotus","Alchornea")
pioneer_tree <- c(
  "Musanga", "Anthocleista", "Elaeis guineensis","Mallotus","Alchornea")

# Compute percentages for each group
pollen_sele = pollen_sele %>%
       mutate(pourc_savanna_tree = rowSums(across(all_of(savanna_tree))/ Total_terrestrial_pollen*100, na.rm = TRUE),
              pourc_forest_tree = rowSums(across(all_of(forest_tree))/ Total_terrestrial_pollen *100, na.rm = TRUE),
              pourc_pioneer_tree = rowSums(across(all_of(pioneer_tree))/ Total_terrestrial_pollen *100, na.rm = TRUE),
              pourc_poaceae = `Poaceae undiff.`/Total_terrestrial_pollen *100)


################################################################################
# Plot pollen percentages
################################################################################

### % savanna tree
plot_savannatree_pollen <- ggplot(pollen_sele, aes(x = Age, y = pourc_savanna_tree)) + 
  geom_line(color = "#f85e00", linewidth = 1) +  
  geom_smooth(color = "#ae2012", fill = "#ae2012", se = FALSE, method = "loess", linetype = "dashed") +  
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), position = "left") +
  ylab("Savanna Tree pollen (%)")+
  theme_without_xaxis +
  tag_side_theme("left")

### % savanna tree / poaceae
plot_savannatree_poaceae_pollen <- ggplot(pollen_sele, aes(x = Age, y = pourc_savanna_tree/pourc_poaceae)) + 
  geom_line(color = "#f85e00", linewidth = 1) +  
  geom_smooth(color = "#ae2012", fill = "#ae2012", se = FALSE, method = "loess", linetype = "dashed") +  
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), position = "left") +
  ylab("Savanna Tree pollen (%)")+
  theme_without_xaxis +
  tag_side_theme("left")

### % Poaceae
plot_poaceae_pollen <- ggplot(pollen_sele, aes(x = Age, y = pourc_poaceae)) +
  geom_line(color = "#ffb703", linewidth = 1) +
  geom_smooth(color = "#ae2012", fill = "#ae2012", se = FALSE, method = "loess", linetype = "dashed") +
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(45, 100), position = "right") +
  ylab("Poaceae pollen (%)")+
  theme_without_xaxis +
  tag_side_theme("right")

### % Forest trees
plot_foresttree_pollen <- ggplot(pollen_sele, aes(x = Age, y = pourc_forest_tree)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_smooth(color = "#ae2012", fill = "#ae2012", se = FALSE, method = "loess", linetype = "dashed") +
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,55), position = "right") +
  ylab("Forest trees pollen (%)")+
  theme_without_xaxis+
  tag_side_theme("right")

### % Pioneer trees
plot_pioneertree_pollen <- ggplot(pollen_sele, aes(x = Age, y = pourc_pioneer_tree )) + 
  geom_line(color = "#38b000", linewidth = 1) +  
  geom_smooth(color = "#006400", fill = "#006400", se = FALSE, method = "loess", linetype = "dashed") +  
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), position = "left") +
  ylab("Pioneer trees pollen (%)")+
  theme_without_xaxis +
  tag_side_theme("left")


# % all trees
plot_arboreal_pollen <- ggplot(pollen_sele, aes(x = Age, y = Pourc_Arbres)) +
  geom_line(color = "#2b9348", linewidth = 1) +
  geom_smooth(color = "#006400", fill = "#006400", se = FALSE, method = "loess", linetype = "dashed") +
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 55), position = "right") +
  ylab("Arboreal pollen (%)")+
  theme_without_xaxis+
  tag_side_theme("right")

# % evergreen
plot_evergreen_pollen <- ggplot(pollen_sele, aes(x = Age, y = Pourc_Evergreen)) +
  geom_line(color = "#2b9348", linewidth = 1) +
  geom_smooth(color = "#006400", fill = "#006400", se = FALSE, method = "loess", linetype = "dashed") +
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(45, NA), position = "left") +
  ylab("Evergreen trees pollen (%)")+
  theme_without_xaxis+
  tag_side_theme("left")

################################################################################
# Charcoal WLRatio: density and mean over time
################################################################################

## Compute density with kde2d
dens_WL = kde2d(charcoal_data$Age,charcoal_data$WLRatio ,lims =c(min(charcoal_data$Age),max(charcoal_data$Age),0,1),n=200)

# Transform to dataframe and normalize
dens_df <- expand.grid(
  Age = dens_WL$x,
  WLRatio = dens_WL$y
) %>%
  mutate(z = as.vector(dens_WL$z))%>%
  group_by(Age) %>%
  mutate(z = z / max(z)) %>%
  ungroup()


# Plot density + smoothed mean
plot_kde_charcoal = ggplot() +
  geom_raster(data = dens_df, aes(x = Age, y = WLRatio, fill = z), interpolate = TRUE) +
  scale_fill_gradient(low = "white", high = "black") +
  geom_point(data = charcoal_data, aes(x = Age, y = WLRatio),
             color = "red", alpha = 0.05, size = 0.3) +
  geom_smooth(data = charcoal_data, aes(x = Age, y = WLRatio), 
              method = "loess", span=0.02,
              color = "#4A90E2", se = FALSE, linewidth = 1.3) +
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),position = "right") +
  labs(y = "Charcoal W/L ratio") +
  theme_without_xaxis+
  theme(
    legend.position = "none"
  )+
  tag_side_theme("right")

################################################################################
# Plot charcoal proportions 
################################################################################

# Tree-derived charcoal proportion
plot_treecharcoal_proportion <- ggplot(p_fuels, aes(x = Age, y = 1-P_G)) +
  geom_line(color = "#006400", linewidth = 1) +
  geom_smooth(color = "#006400", fill = "#006400", se = FALSE, method = "loess", linetype = "dashed") +
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), position = "left") + 
  # ylab("# charcoal/cm3")+
  ylab("Woody charcoal\nproportion")+
  theme_without_xaxis+
  tag_side_theme("left")

# Prepare long dataframe for stacked plot (grass/tree)
p_fuels_long <- p_fuels %>% 
  dplyr::select(Depth,Age,P_H)%>%
  mutate(P_W = 1- P_H) %>%
  pivot_longer(cols = c(P_H, P_W),
               names_to = "Vegetation",
               values_to = "Proportion") %>%
  mutate(Vegetation = recode(Vegetation,
                             "P_H" = "Herbaceous",
                             "P_W" = "Woody"))%>%
  mutate(Vegetation = factor(Vegetation, levels = c("Herbaceous", "Woody")))

# Stacked area plot for charcoal proportion
plot_charcoal_proportion = ggplot(p_fuels_long, aes(x = Age, y = Proportion, fill = Vegetation)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = c("Herbaceous" = "#ffb703", "Woody" = "#006400"),guide="none") +
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), position = "left", expand = c(0, 0)) + 
  # ylab("Proportion of Woody / Herbaceous charcoal")+
  ylab("Proportion of Woody/\nHerbaceous charcoal")+
  theme_without_xaxis+
  tag_side_theme("left")+
  theme(legend.position = "right")

################################################################################
# Plot Charcoal concentration / cm3
################################################################################


# Total charcoal concentration
plot_charcoal_concentration <- ggplot(p_fuels, aes(x = Age, y = Concentration)) +
  geom_area(fill = "#f85e00", alpha = 0.4) +  # colouring under the curve
  geom_line(color = "#f85e00", linewidth = 1) +
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 200), position = "left") + 
  # ylab("# charcoal/cm3")+
  ylab("Charcoal\nconcentration (#/cm3)")+
  theme_without_xaxis+
  tag_side_theme("left")

# Tree-derived charcoal concentration 
plot_tree_charcoal_concentration <- ggplot(p_fuels,aes(x = Age, y = Concentration*P_W)) +
  geom_area(fill = "#006400", alpha = 0.4) +  # colouring under the curve
  geom_line(color = "#006400", linewidth = 1) + 
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), position = "right", expand = c(0, 0)) + 
  # ylab("Tree charcoal") +
  labs(x = "Age cal. BP", y = "Woody charcoal\nconcentration (#/cm3)") +
  theme_with_xaxis+
  tag_side_theme("right")

# Grass-derived charcoal concentration / cm
plot_grass_charcoal_concentration <- ggplot(p_fuels,aes(x = Age, y = Concentration*P_H)) +
  geom_area(fill = "#ffb703", alpha = 0.4) +
  geom_line(color = "#ffb703", linewidth = 1) +
  scale_x_reverse(limits = c(MaxAge, MinAge), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), position = "left", expand = c(0, 0)) + 
  # labs(x = "Age cal. BP", y = "Grass charcoal") +
  ylab("Herbaceous charcoal\nconcentration (#/cm3)") +
  theme_without_xaxis+
  tag_side_theme("left")

################################################################################
# Combine all plots into final figure
################################################################################

final_plot <-   plot_charcoal_concentration /
  plot_spacer() /
  plot_kde_charcoal /
  plot_spacer()/
  plot_charcoal_proportion / 
  plot_spacer()/ 
  plot_poaceae_pollen /
  plot_spacer()/
  plot_grass_charcoal_concentration /
  plot_spacer()/
  plot_arboreal_pollen  / 
  plot_spacer()/
  plot_pioneertree_pollen / 
  plot_spacer()/
  plot_tree_charcoal_concentration + 
  plot_layout(axes = "collect")+
  plot_annotation(tag_levels = 'A')+
  plot_layout(heights = c(1, 
                          -.1, 1,
                          -.1,1,
                          -.1,1,
                          -.8,1,
                          -.1,1,
                          -.8,1,
                          -.4,1,
                          -.5,1))


# Display plot
final_plot

# Save plot
ggsave("Output/Figure 5 Synthesis Sele pollen charcoal.pdf", final_plot, width = 5.16, height = 9.17, units = "in")



