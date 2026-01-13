################################################################################
################## Computation of Tree Covers in lanscape subsets ##############
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
library(raster)       
library(ggplot2)      
library(ggpubr)      
library(progress)
library(dplyr)
library(tidyr)
library(purrr)
library(ggrepel)

if (!dir.exists("Output")) {
  dir.create("Output")
}

################################################################################
# Utilities
################################################################################

# Define Land Cover Colors and Labels
landcover_colors <- c(
  "1" = "#3B9AB2",   # WATER
  "2" = "#C19A6B",   # BARE SOIL
  "3" = "#C0C0C0",   # URBAN AREA
  "4" = "#FFFF99",   # FIELDS
  "5" = "#F4A460",   # SAVANNA
  "6" = "#ABDDA4",   # DENSER SAVANNA
  "7" = "#66A61E",   # TREE CROPS
  "8" = "#7FCDBB",   # WETLAND
  "9" = "#006766",   # INUNDATED FOREST
  "10" = "#006400"   # FOREST
)
landcover_labels <- c(
  "Water", "Bare Soil", "Urban Area", "Fields", "Savanna",
  "Dense Savanna", "Tree Crops", "Wetland", "Inundated forest", "Forest"
)
colors_legend = data.frame(category = landcover_labels,color = landcover_colors) %>%
  mutate(category = factor(category, levels = landcover_labels))

landcover_colors_cat = landcover_colors
names(landcover_colors_cat) = landcover_labels

# Create Combination of Parameters for Tree Cover computations : 
# TC product resolution, radii, 
# masks (savanna, savanna + agricultural land, burned areas)
resolutions <- c("TC1m", "TC10m")
distances <- c(1000, 2000, 3000, 4000, 5000)
masks <- c("classif_buffer_savanna", "classif_buffer_savanna_agr", "baa_classif_t_1_buf","baa_classif_t_resclassif")
codes = c("sav","savagr","burned")

raster_combinations = expand.grid(
  distance = distances,
  resolution = resolutions,
  mask = masks,
  stringsAsFactors = FALSE
) %>%
  mutate(code = recode(mask,
                       "classif_buffer_savanna" = "sav",
                       "classif_buffer_savanna_agr" = "savagr",
                       "baa_classif_t_1_buf" = "burned_once",
                       "baa_classif_t_resclassif" = "burned_weighted_mean"))

# Function to compute mean tree cover for a given combination
calc_mean <- function(resolution, distance, mask) {
  expr <- parse(text = paste0(resolution, "_classif*disk_masks$radius_", distance, "*", mask))
  expr2 <- parse(text = paste0("disk_masks$radius_", distance, "*", mask))
  return(sum(values(eval(expr)), na.rm = TRUE)/sum(values(eval(expr2)), na.rm = TRUE))
}

# Site list
sites = c("NON","BOC","NIAK","NAP","NIN","DTO","SLE","AZL","FBA","GRO","DINA",
          "ELI","ESB","NDA","R1","NGO","DKL","NGG","GBL")

################################################################################
# Computation of landscape, vegetation and fire statistics
################################################################################

# Initialize Result Tables
results_fires = data.frame()
results_treecover = data.frame()
results_tc_class = data.frame()
results_tc_fires_per_class = data.frame()

# Progress Bar
pb <- progress_bar$new(
  total = length(sites),
  format = "Processing sites [:bar] :percent eta: :eta"
)

# Main loop over sites
for (site in sites) {
  
  ################################################################################
  # Load rasters 
  ################################################################################
  
  # Load vegetation classification
  classif = raster(paste0("Sites/",site,"/",site, "_Classif.tif"))
  
  # Compute site center coordinates
  x <- xmin(classif) + (xmax(classif) - xmin(classif)) / 2
  y <- ymin(classif) + (ymax(classif) - ymin(classif)) / 2
  centre = SpatialPoints(cbind(x, y), proj4string = crs(classif))
  
  # Define UTM projection for the site
  utm_crs = CRS(paste0("+proj=utm +zone=",
                       floor((x + 180)/6) + 1 , 
                       ifelse(y >= 0, " +north", " +south"),
                       " +datum=WGS84 +units=m +no_defs"))
  
  # Compute a raster with distance from the centre
  center_utm <- spTransform(centre, utm_crs)
  dist_raster = raster(
    extent(x-.05, x+.05, y-.05, y+.05),
    crs = crs(classif),
    res = res(classif)
  ) %>% 
    projectRaster(crs = utm_crs) %>%
    distanceFromPoints(center_utm) %>%
    projectRaster( crs = crs(classif)) %>%
    resample(classif, method = "bilinear") %>%
    crop(extent(classif))
  
  # Create binary masks for different radii
  disk_masks <- setNames(lapply(distances, function(d) {
    r <- dist_raster
    values(r) <- ifelse(values(r) <= d, 1, NA)
    r
  }), paste0("radius_", distances))
  
  ##############################################################################
  # Load buffered vegetation rasters and apply 5 km mask
  ##############################################################################
  # Vegetation classification without holes with buffer 20m around each class
  classif_buffer = raster(paste0("Sites/",site,"/",site, "_Classif_smooth2_buffer.tif"))*disk_masks$radius_5000
  # Mask for class 5&6 with buffer 20m
  classif_buffer_savanna = raster(paste0("Sites/",site,"/",site, "_Classif_smooth_buffer_savanna.tif"))*disk_masks$radius_5000
  # Mask for class 4 to 7 with buffer 20m
  classif_buffer_savanna_agr = raster(paste0("Sites/",site,"/",site, "_Classif_smooth_buffer_savanna_agr.tif"))*disk_masks$radius_5000
  
  
  ##############################################################################
  # Fire rasters: at least 1 fire, 2 fires, and number of fires per pixel
  ##############################################################################
  # Enveloppes burnt at least once with buffer
  baa_classif_t_1_buf = raster(paste0("Sites/",site,"/",site, "_BA_GABAM_1985_2021_1fire_resclassif_buffer.tif"))*disk_masks$radius_5000
  # Enveloppes burnt at least twice with buffer
  baa_classif_t_2_buf = raster(paste0("Sites/",site,"/",site, "_BA_GABAM_1985_2021_2fires_resclassif_buffer.tif"))*disk_masks$radius_5000
  # Create mask for regions burnt at least once
  baa_classif_t_1_buf_0 = raster(baa_classif_t_1_buf)
  values(baa_classif_t_1_buf_0) <- ifelse(is.na(values(baa_classif_t_1_buf)), 0, 1)
  # Number of fires in each pixel inside buffer
  baa_classif_t_resclassif = raster(paste0("Sites/",site,"/",site, "_BA_GABAM_1985_2021_sum_resclassif.tif"))*baa_classif_t_1_buf_0
  
  ##############################################################################
  # Tree cover rasters
  ##############################################################################
  # TC JRC
  TC10m_classif = raster(paste0("Sites/",site,"/",site, "_TC_JRC_GFC2020_V1.tif"))*disk_masks$radius_5000
  
  # TREE COVER 1m
  TC1m_classif = raster(paste0("Sites/",site,"/",site, "_TC_1m_2019_resclassif.tif"))*disk_masks$radius_5000
  
  ##############################################################################
  # Compute tree cover per vegetation class over 5km
  ##############################################################################
  tc_per_class_site = data.frame(zonal(TC1m_classif,classif_buffer, fun = 'mean', na.rm = TRUE)) %>% 
    rename(TC1m = mean) %>%
    left_join(data.frame(zonal(TC10m_classif,classif_buffer, fun = 'mean', na.rm = TRUE)),by=join_by(zone))%>%
    rename(TC10m = mean, class = zone)
  
  ##############################################################################
  # Compute fire characteristics per vegetation class over 5km
  ##############################################################################
  ba_per_class_site = data.frame(freq(classif_buffer,useNA = "no")) %>% 
    rename(totalarea = count) %>%
    left_join(data.frame(freq(baa_classif_t_1_buf*classif_buffer,useNA = "no")), by = "value") %>%
    rename(area_1fire = count) %>%
    mutate(prop_1fire = area_1fire/totalarea) %>%
    left_join(data.frame(freq(baa_classif_t_2_buf*classif_buffer,useNA = "no")), by = "value") %>%
    rename(area_2fires = count) %>%
    mutate(prop_2fires = area_2fires/totalarea) %>%
    left_join(data.frame(zonal(baa_classif_t_resclassif,classif_buffer, fun = 'mean', na.rm = TRUE)),by=join_by(value==zone)) %>%
    rename(Weighted_mean_fire = mean, class=value)
  # Merge Results
  results_tc_fires_per_class= results_tc_fires_per_class %>% bind_rows(
    left_join(ba_per_class_site,tc_per_class_site,by=join_by(class)) %>%
      mutate(Site=site)%>%
      dplyr::select(Site, everything()))
  
  
  ##############################################################################
  # Compute tree cover for different resolutions, radii, and masks
  ##############################################################################
  results_treecover_site <- raster_combinations %>%
    mutate(value = pmap_dbl(list(resolution, distance, mask), calc_mean),
           col_name = paste(resolution, distance, code, sep = "_")) %>%
    dplyr::select(col_name, value) %>%
    pivot_wider(names_from = col_name, values_from = value) 
  # Merge Results
  results_treecover = results_treecover %>% 
    bind_rows(results_treecover_site%>%
                mutate(Site = site)%>%
                dplyr::select(Site, everything()))
  
  # Update progress bar
  pb$tick()
}

# Save table
write.csv(results_treecover,"Output/TC vegclass disks.csv",row.names = F)

################################################################################
# Post-processing results_tc_fires_per_class 
################################################################################

# Complete missing classes, add labels, compute relative areas
legend = data.frame(class = 1:10, category = landcover_labels)
results_tc_fires_per_class <- results_tc_fires_per_class %>%
  left_join(legend , by=join_by(class)) %>%
  mutate(class = factor(class),
         category = factor(category, levels = landcover_labels), 
         Site = as.factor(Site), 
         prop_1fire = ifelse(is.na(prop_1fire)& !is.na(totalarea),0,prop_1fire),
         prop_2fires = ifelse(is.na(prop_2fires)& !is.na(totalarea),0,prop_2fires)
  )%>%
  complete(Site, category, fill = list(mean = NA))


# Save table
 write.csv(results_tc_fires_per_class,"Output/TC and fires per site and vegclass.csv",row.names = F)

 

# Compute relative areas of vegetation classes with/without fire
 results_tc_fires_per_class <- read.csv("Output/TC and fires per site and vegclass.csv",header=T, dec=".",check.names = F)
 results_fire = results_tc_fires_per_class %>% 
  group_by(Site) %>% 
  mutate(area_rep = totalarea/sum(totalarea, na.rm = T)) %>% 
  mutate(area_1fire = prop_1fire*area_rep, area_nofire = (1-prop_1fire)*area_rep) %>%
  dplyr::select(Site,category,area_1fire,area_nofire) %>%
  pivot_longer(cols = c("area_1fire","area_nofire"), names_to = "fire", values_to = "reparea" ) %>%
  mutate(fire = factor(fire, levels = c( "area_nofire","area_1fire"))) %>%
  left_join(colors_legend, by=join_by(category))%>%
  mutate(color = ifelse(fire == "area_nofire", color, "red"))



################################################################################
# Visualization
################################################################################

# Compare TC1m vs TC10m for vegetation classes 4 to 7
plot_compTC = ggplot(results_tc_fires_per_class %>% filter(class %in% c(4,5,6,7)) )+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey80",linewidth=1) +
  geom_point(aes(x=TC10m,y=TC1m,color=category,shape = category),size=3)+
  geom_smooth(aes(x = TC10m, y = TC1m, group = category,color=category),method = "lm", se = FALSE, show.legend = F)+
  scale_color_manual(values = landcover_colors_cat) +  # assign colors per category
  scale_shape_manual(values = 15:18) +
  xlim(0,1)+
  ylim(0,1)+
  labs(
    title = "",
    x = "Tree Cover JRC (10m resolution)",
    y = "Tree cover Reiner et al. (1m resolution)"
  ) +
  coord_fixed(ratio = 1)+
  theme_pubr()+
  theme(legend.position = "right", legend.title = element_blank())
plot(plot_compTC)
ggsave("Output/Figure S15 comparision TC 1m 10m.png", plot = plot_compTC, width = 6, height = 6, dpi = 600)


# Compare TC10m in savannas vs savannas + agricultural land
plot_comp_sav_savagr = 
  ggplot(results_treecover)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey80",linewidth=1) +
  geom_point(aes(x=TC10m_5000_sav,y=TC10m_5000_savagr,color=Site))+
  geom_text_repel(aes(x = TC10m_5000_sav, y = TC10m_5000_savagr, label = Site), 
                  max.overlaps = Inf, # allow all labels, adjust dynamically,
                  nudge_x = 0.03,nudge_y = 0.0,
                  min.segment.length = 0, seed = 42, box.padding = 0.1,
                  point.padding = 0.0,     # spacing around points
                  size = 3) +
  coord_fixed(ratio = 1)+
  xlim(0,1)+
  ylim(0,1)+
  labs(
    title = "Tree Cover JRC (10m resolution)",
    x = "TC in savannas",
    y = "TC in savannas and agriculatural land"
  ) +
  theme_pubr()+
  theme(legend.position = "none", legend.title = element_blank()) 

plot(plot_comp_sav_savagr)
ggsave("Output/Comparision TC in savanna and sav_agr.png", plot = plot_comp_sav_savagr, width = 6, height = 6, dpi = 600)


