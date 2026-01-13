################################################################################
##################### Final summary multiplot for each site ####################
################################################################################
# Ref : Cornet et al. in prep
# Date: 2025-08-27
# Purpose: Plot data for each site in one panel
################################################################################

# Load libraries
library(grid)
library(gridExtra)
library(png)
library(progress)
library(patchwork)
library(raster)       
library(ggplot2)
library(ggspatial)
library(ggnewscale)
library(dplyr)
library(tidyr)
library(ggpubr)

# Charcoal data
charcoal_data = read.csv("Output/Charcoal_data_kept.csv", header=T,  dec=".",check.names = F)


# Age models
age_models = read.csv("Input/Agemodels_all_sites.csv", header=T, sep=";", dec=".",check.names = F)

age_models = age_models %>% 
  group_by(Site) %>%
  arrange(Depth) %>%
  mutate(
    # Différence centrée pour points internes
    Accrate = -(lead(Depth) - lag(Depth))/(lead(Age) - lag(Age)) ,
    # Bord gauche : forward difference
    Accrate = ifelse(row_number() == 1, -(lead(Depth) - Depth)/(lead(Age) - Age) , Accrate),
    # Bord droit : backward difference
    Accrate = ifelse(row_number() == n(), -(Depth - lag(Depth))/(Age - lag(Age)) , Accrate)
  ) %>%
  ungroup() %>%
  arrange(Site,Depth)


# Site metadata
site_metadata = read.csv("Input/Sites metadata.csv", header=T, sep=";", dec=".",check.names = F)

# Associate both tables and filter 
#Ensure equitable charcoal reparition in cas Volume changes in sequence

charcoal_data_KB = charcoal_data %>%  dplyr::select(-Age) %>% 
  mutate(Site = ifelse(Site=="DKL(b)","DKL",Site))%>%
  group_by(Site,Depth) %>%
  summarise(N = n(), Areatot = sum(ProjArea), WLmoy = mean(WLRatio)) %>%
  group_by(Site) %>%
  complete(Depth = 0:max(Depth)) %>%
  arrange(Site, Depth) %>%
  mutate(N = replace_na(N, 0),
         Areatot = replace_na(Areatot, 0)) %>%
  left_join(age_models,by=join_by(Site,Depth)) %>%
  mutate(CharConc = N/Volume, AreaConc = N/Volume, CharInflux = CharConc * Accrate, AreaInflux = AreaConc * Accrate) %>%
  filter(Age >= 1990) 

charcoal_stats = charcoal_data %>% 
  filter(keep) %>% 
  group_by(Site,Depth) %>% 
  summarise(n=n()) %>% 
  group_by(Site) %>% 
  summarise(Charcoal_total = sum(n), Mean=mean(n),Median = median(n),StandardDeviataion = sd(n)) %>%
  left_join(site_metadata %>% dplyr::select(Site,Country), by="Site") %>%
  arrange(Country, Site)

write.csv(charcoal_stats, "Output/Charcoal_site_statistics.csv", row.names = FALSE)

table_GMM_results = read.csv("Output/Bayesian_GMM_outputs_30_yrs.csv",h=T,check.names=F)

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
# Site names
################################################################################
sites = unique(charcoal_data_KB$Site)

################################################################################
# Plots figures and maps
################################################################################

# Progress Bar
pb <- progress_bar$new(
  total = length(sites),
  format = "Processing sites [:bar] :percent eta: :eta"
)


# Loop over sites
for (site in sites) {
  
  list_files <- list.files(path = paste0("Sites/",site), pattern = "\\.png$", ignore.case = TRUE, full.names = TRUE)
  AgeModel_file = list_files[grepl(paste0("^", site, "_.*_AM\\.png$"), basename(list_files), ignore.case = TRUE)]
  
  
  ######
  #  Age model image 
  #####
  
  if(length(AgeModel_file)==0){
    plot_agemodel = plot_spacer()
  }else{
    agemodel <- readPNG(AgeModel_file)
    grob_agemodel <- rasterGrob(agemodel, interpolate = TRUE)
    plot_agemodel <- wrap_elements(full = grob_agemodel)+
      theme(plot.margin = ggplot2::margin(),        
            plot.tag = element_text(face = "bold", size = 16))
  }
  # Site charcoal data
  site_data <- subset(charcoal_data_KB, Site == site)
  
  #Change y concentration scale according to concentraiton range
  if (max(site_data$CharConc, na.rm = TRUE) < 50){
    ylim_values <- c(0, 50)              # par défaut
    line_color <- "#E35D5D"  
  }else if ( (max(site_data$CharConc, na.rm = TRUE) < 200)){
    ylim_values <- c(0, 200)              # par défaut
    line_color <- "#F3B94D"  
  }else{
    ylim_values <- c(0, 500)              # par défaut
    line_color <- "#A04E98"  
  }
  
  #Scale factor concentration/influx
  scale_factor <- max(site_data$CharConc, na.rm = TRUE) / max(site_data$CharInflux, na.rm = TRUE)
  
  ######
  #  Plot concentration and influx / Time
  #####
  if(length(AgeModel_file)==0){
    plot_charconc <- ggplot(site_data, aes(x = Age)) +
      geom_point(aes(y = CharConc,color = "Concentration"), size = 2) +      # couleur définie ici
      scale_color_manual(name=NULL, values = c("Concentration" = line_color))+
      labs(
        x = "Age CE",
        y = "Charcoal concentration (#/cm3)") +
      scale_x_continuous(limits = c(1990, 2025))+
      coord_cartesian(ylim = ylim_values) +
      theme_pubr(base_size = 14)+
      theme(axis.title.y.right = element_text(angle = 90, vjust = 0.5) )+
      theme(legend.position = c(0.7, 0.95),
            legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
            legend.key = element_blank(),
            legend.text = element_text(size = 12),
            legend.spacing.y = unit(0, "cm"),
            plot.tag = element_text(face = "bold", size = 16))
  }else{
    plot_charconc <- ggplot(site_data, aes(x = Age)) +
      geom_line(aes(y = CharConc,color = "Concentration"), linewidth = 1) +      # couleur définie ici
      geom_line(aes(y = CharInflux * scale_factor, color = "Influx"), linewidth = 1, linetype = "dashed") +    
      scale_color_manual(name=NULL, values = c("Concentration" = line_color, "Influx" = "#3A8486"))+
      labs(
        x = "Age CE",
        y = "Charcoal concentration (#/cm3)") +
      scale_x_continuous(limits = c(1990, 2025))+
      guides(color = guide_legend(nrow = 2))+
      scale_y_continuous(sec.axis = sec_axis(~./scale_factor, name = "Charcoal influx (#/cm2/yr)"))+
      coord_cartesian(ylim = ylim_values) +
      theme_pubr(base_size = 14)+
      theme(axis.title.y.right = element_text(angle = 90, vjust = 0.5) )+
      theme(legend.position = c(0.7, 0.95),
            legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
            legend.key = element_blank(),
            legend.text = element_text(size = 12),
            legend.spacing.y = unit(0, "cm"),
            plot.tag = element_text(face = "bold", size = 16))
    
  }
#  ggsave(filename = paste0("Sites/",site,"/CharConc_Influx_", site, ".png"),
#         plot = plot_charconc, width = 6, height = 6/8*6, dpi = 600, bg = "white")
  
  
  ######
  #  Plot WL / Time 
  #####
  if(length(AgeModel_file)==0){
    plot_WL <- ggplot(site_data%>% filter(!is.na(WLmoy)), aes(x = Age)) +
      geom_point(aes(y = WLmoy),color = "red", size = 2) +      
      geom_hline(yintercept = 0.5,linetype="dashed")+
      labs(
        x = "Age CE",
        y = "Mean WL ratio per sample") +
      theme_pubr() +
      # scale_x_reverse(limits = c(2020, 1990)) +
      scale_x_continuous(limits = c(1990, 2025))+
      scale_y_continuous(limits = c(0,1),expand=c(0,0)) +
      theme_pubr(base_size = 14)+
      theme(plot.tag = element_text(face = "bold", size = 16))
  }else{
    plot_WL <- ggplot(site_data%>% filter(!is.na(WLmoy)), aes(x = Age)) +
      geom_line(aes(y = WLmoy),color = "red", linewidth = 1) +      
      geom_hline(yintercept = 0.5,linetype="dashed")+
      labs(
        x = "Age CE",
        y = "Mean WL ratio per sample") +
      theme_pubr() +
      # scale_x_reverse(limits = c(2020, 1990)) +
      scale_x_continuous(limits = c(1990, 2025))+
      scale_y_continuous(limits = c(0,1),expand=c(0,0)) +
      theme_pubr(base_size = 14)+
      theme(plot.tag = element_text(face = "bold", size = 16))
  }
#  ggsave(filename = paste0("Sites/",site,"/WLRatio_mean_", site, ".png"),
#         plot = plot_WL, width = 6, height = 6/8*6, dpi = 600, bg = "white")
  
  ######
  #  Plot distributions WLratio and Box Cox logit transformed WLratio
  ##### 
  
  # Compute Distributions
  table_site = table_GMM_results %>% filter(Site==site)
  
  xx <- seq(1e-4, to = 1 - 1e-4, length.out = 512)
  dx <- diff(xx)[1]
  
  distributions <- data.frame(
    x = xx[-1] - dx/2,
    herbaceous = table_site$`theta[1]50%` * diff(pnorm(BoxCoxLogit(xx, table_site$`lambda50%`), table_site$`mu[1]50%`, table_site$`sigma[1]50%`)) / dx,
    woody = table_site$`theta[2]50%` * diff(pnorm(BoxCoxLogit(xx, table_site$`lambda50%`), table_site$`mu[2]50%`, table_site$`sigma[2]50%`)) / dx
  )
  distributions$total <- distributions$herbaceous + distributions$woody
  
  ######
  # Plot Raw WL Ratio distribution
  #####
  
  charcoal_data_site = charcoal_data %>% filter(Site==site & keep)
  
  plot_GMM_WLratio <- ggplot(data = charcoal_data_site) +
    geom_density(aes(x = WLRatio), color = "#1e88e5", fill = "#1e88e5", alpha = 0.3, adjust=1) +
    geom_line(data = distributions, aes(x = x, y = herbaceous), color = "#e01e37", size = 1) +
    geom_line(data = distributions, aes(x = x, y = woody), color = "#a7c957", size = 1) +
    geom_line(data = distributions, aes(x = x, y = total), color = "black", linetype = "dashed", size = 1) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    annotate("text",
             x = Inf, y = Inf,               
             label = paste0("N = ",nrow(charcoal_data_site)),
             hjust = 1.1, vjust = 2,          
             size = 6, color = "black")+
    theme_pubr(base_size = 14) +
    labs(x = "WL Ratio", y = "Density")+
    theme(plot.caption = element_text(hjust = 1),
          plot.tag = element_text(face = "bold", size = 16))
  
  
  ######
  # Plot Transformed WL Ratio 
  #####
  plot_GMM_BCL_WLratio <- ggplot(data = charcoal_data_site) +
    geom_density(aes(x = BoxCoxLogit(WLRatio, lambda=table_site$`lambda50%`)), color = "#1e88e5", fill = "#1e88e5", alpha = 0.3, adjust=1) +
    stat_function(
      fun = function(x) table_site$`theta[1]50%` * dnorm(x, table_site$`mu[1]50%`, table_site$`sigma[1]50%`),
      color = "#e01e37", size = 1) +
    stat_function(
      fun = function(x)  table_site$`theta[2]50%` * dnorm(x, table_site$`mu[2]50%`, table_site$`sigma[2]50%`),
      color = "#a7c957", size = 1) +
    stat_function(
      fun = function(x) table_site$`theta[1]50%` * dnorm(x, table_site$`mu[1]50%`, table_site$`sigma[1]50%`) +
        table_site$`theta[2]50%` * dnorm(x, table_site$`mu[2]50%`, table_site$`sigma[2]50%`),
      color = "black", linetype = "dashed", size = 1) +
    annotate("text",
             x = Inf, y = Inf,               
             label = paste0("N = ",nrow(charcoal_data_site)),
             hjust = 1.1, vjust = 2,          
             size = 6, color = "black")+
    scale_x_continuous(limits = c(-4.5, 4.5), expand = c(0, 0)) +
    theme_pubr(base_size = 14) +
    labs(x = "Box-Cox logit-transformed WL Ratio", y = "Density")+
  theme(plot.tag = element_text(face = "bold", size = 16))

  ######
  # Plot relative distributions of vegetation classes with and without fire per site
  ######
  
    plot_vegclass_distr_fire = ggplot(results_fire%>%filter(Site==site), aes(x = category, y = reparea, fill = color)) +
      geom_bar(stat = "identity") +
      scale_fill_identity(
        guide = "legend",
        name = NULL,
        labels = c("red" = "At least 1 fire detected"),
        breaks = c("red")
      )+
      labs(
        x = "",
        y = "Relative proportion"
      ) +
      scale_y_continuous(expand = c(0, 0))+
      theme_pubr(base_size = 14)+
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5,margin = ggplot2::margin(b = 50)),
            legend.title = element_blank(),
        plot.tag = element_text(face = "bold", size = 16))
    
  
  ################################################################################
  # Load and mask vegetation and fire rasters
  ################################################################################
  
    # Load land cover classification raster
    classif = raster(paste0("Sites/",site,"/",site, "_Classif.tif"))
    
    # Get center coordinates of the raster
    x <- xmin(classif) + (xmax(classif) - xmin(classif)) / 2
    y <- ymin(classif) + (ymax(classif) - ymin(classif)) / 2
    centre = SpatialPoints(cbind(x, y), proj4string = crs(classif))
    
    # Define UTM projection for site
    utm_crs = CRS(paste0("+proj=utm +zone=",
                         floor((x + 180)/6) + 1 , 
                         ifelse(y >= 0, " +north", " +south"),
                         " +datum=WGS84 +units=m +no_defs"))
    
    # Compute a raster with distance from the centre
    center_utm <- spTransform(centre, utm_crs)
    disk_mask_5000 = raster(
      extent(x-.05, x+.05, y-.05, y+.05),
      crs = crs(classif),
      res = res(classif)
    ) %>% projectRaster(crs = utm_crs) %>%
      distanceFromPoints(center_utm) %>%
      projectRaster( crs = crs(classif)) %>%
      resample(classif, method = "bilinear") %>%
      crop(extent(classif))
    
    # Create 5 km mask (1 inside disk, NA outside)
    values(disk_mask_5000) <- ifelse(values(disk_mask_5000) <= 5000, 1, NA)
    
    
    # Apply 5km-mask to vegetation classification
    classif = classif*disk_mask_5000
    
    # Load number of fires in each pixel and apply 5km-mask
    baa_classif_t_resclassif = raster(paste0("Sites/",site,"/",site, "_BA_GABAM_1985_2021_sum_resclassif.tif"))*disk_mask_5000
    
    # Load TC JRC
    tc_classif = raster(paste0("Sites/",site,"/",site, "_TC_JRC_GFC2020_V1.tif"))*disk_mask_5000
    
    # Load TREE COVER 1m
    tc1m_classif_resclassif = raster(paste0("Sites/",site,"/",site, "_TC_1m_2019_resclassif.tif"))*disk_mask_5000
    
    ################################################################################
    # Define color ramp for fire intensity visualization
    ################################################################################
    max_val <- ceiling(max(values(baa_classif_t_resclassif), na.rm = TRUE))
    colors_fire <- c(rgb(1,1,1,0),colorRampPalette(c("#FFB6C160", "#FF000060", "#8B000060"),alpha = T)(max_val))  
  
    
    ######
    # Plot Landcover
    ######
    
    plot_landcover = ggplot() +
           layer_spatial(data = classif,aes(fill = as.factor(after_stat(band1))))+
      scale_fill_manual(breaks=1:10,values=landcover_colors,  na.value = NA,labels = landcover_labels, name="")+
      labs(title = paste0(site, " / Land Cover Map"))+
      theme_void()+
      theme(
        plot.title = element_text(hjust = 0.5,size = 16),
        legend.position = "right",
        plot.tag = element_text(face = "bold", size = 16)
      )
    ######
    # Plot Land Cover + Fire Map
    ######
    
    plot_fire_map = ggplot() +
      layer_spatial(data = classif,aes(fill = as.factor(after_stat(band1))))+
      scale_fill_manual(breaks=1:10,values=landcover_colors,  na.value = NA,labels = landcover_labels, name="",guide="none")+
      layer_spatial(data = baa_classif_t_resclassif,aes(fill = as.factor(after_stat(band1))))+
      new_scale_fill() +
      layer_spatial(
        data = baa_classif_t_resclassif,
        aes(fill = after_stat(band1))
      ) +
      scale_fill_gradientn(
        colors = colors_fire,
        limits=c(0,max_val),
        breaks = 0:(max_val),
        labels = c(0,rep("",max_val-1),max_val),
        name = "Number of fires",
        na.value = NA
      )+
      labs(title = paste0(site, " / Burned Areas"))+
      theme_void(base_size = 14)+
      theme(
        plot.title = element_text(hjust = 0.5,size = 14),
        legend.position = "right",
        plot.tag = element_text(face = "bold", size = 16)
      )
    ######
    # Tree Cover 10m
    ######
    plot_tc10m = ggplot() +
      layer_spatial(data = tc_classif,aes(fill = as.factor(after_stat(band1))))+
      scale_fill_manual(breaks=c(1,0),values=(terrain.colors(2)),  na.value = NA, name="")+
      labs(title = paste0(site, " / TC 10m"))+
      theme_void()+
      theme(
        plot.title = element_text(hjust = 0.5,size = 14),
        legend.position = "right",
        plot.tag = element_text(face = "bold", size = 16)  # gras + taille
      )
    
    ######
    # Tree Cover 1m
    ######
    
    plot_tc1m = ggplot() +
      layer_spatial(data = tc1m_classif_resclassif,aes(fill = after_stat(band1)))+
      scale_fill_gradientn(
        colors = rev(terrain.colors(100)),  # palette continue
        limits = c(0, 1), 
        na.value = NA,
        name = "TC"
      )+
      labs(title = paste0(site, " / TC 1m"))+
      theme_void()+
      theme(
        plot.title = element_text(hjust = 0.5,size = 16),
        legend.position = "right",
        plot.tag = element_text(face = "bold", size = 16)  # gras + taille
      )
    
    ######
    # Patchwork of plots
    ######
    
    plot_fig_AF =  (plot_agemodel | plot_charconc | plot_WL) /  
   ( plot_GMM_WLratio | plot_GMM_BCL_WLratio | plot_vegclass_distr_fire   ) +
      plot_layout( heights = c(1, 1))+
      plot_annotation(tag_levels = "A") 
    
    plot_fig_GJ =  (plot_landcover | plot_fire_map) /
   (plot_tc10m | plot_tc1m)+
      plot_layout( heights = c(1,1),tag_level = 'new')+
      plot_annotation(tag_levels = list(c("G","H","I","J"))) & 
      theme(plot.margin = margin(0, 0, 0, 0)) 
  
 
    # Sauvegarde du graphique
    ggsave(paste0("Sites/",site,"/",site, "_Fig_Appendix_S7_",site,"_AF.png"),plot_fig_AF,width=18,height=8.5,dpi = 600)
    ggsave(paste0("Sites/",site,"/",site, "_Fig_Appendix_S7_",site,"_GJ.png"),plot_fig_GJ,width=8.5,height=8.5,dpi = 600)
    

  # Update progress bar
  pb$tick()
  
}
