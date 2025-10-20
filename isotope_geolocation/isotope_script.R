# Isotope geolocation script for 'An epigenetic aging clock in butterfies'
# author: Megan S Reich
# contact: mreich@uottawa.ca; meganreich13@gmail.com
# description: R script to apply hydrogen and strontium isotope-based geographic
  # assignment, combined with a monthly species distribution model, to estimate
  # the natal origin of 5 painted lady butterflies captured in Benin and Senegal.
  # From each, estimates of the  minimum distance-traveled and maximum distance
  # traveled are calculated.
# comments: run in R version 4.5.1
# Sr Isoscape manuscript citation:
#   Reich, M. S., Ghouri, S., Zabudsky, S., Hu, L., Le Corre, M., Ng’iru, I., Benyamini, D., Shipilina, D., Collins, S. C., Martins, D. J., Vila, R., Talavera, G., & Bataille, C. P. (2024). Trans-Saharan migratory patterns in Vanessa cardui and evidence for a southward leapfrog migration. iScience, 27(12), 111342. https://doi.org/10.1016/j.isci.2024.111342
#
# Sr Isoscape data citation and location:
#   Reich, M. S., Ghouri, S., Zabudsky, S., Hu, L., Le Corre, M., Ng’iru, I., … bataille, clement. (2024, October 30). Isotope geolocation. Retrieved from osf.io/rnfmp
#
# d2H Calibration data citation:
#   Ghouri, S., Reich, M. S., Lopez‐Mañas, R., Talavera, G., Bowen, G. J., Vila, R., Talla, V. N. K., Collins, S. C., Martins, D. J., & Bataille, C. P. (2024). A hydrogen isoscape for tracing the migration of herbivorous lepidopterans across the Afro‐Palearctic range. Rapid Communications in Mass Spectrometry, 38(3), e9675. https://doi.org/10.1002/rcm.9675

# species distribution model:
#   Talavera, G., García-Berro, A., Talla, V. N. K., Ng’iru, I., Bahleman, F., Kébé, K., Nzala, K. M., Plasencia, D., Marafi, M. A. J., Kassie, A., Goudégnon, E. O. A., Kiki, M., Benyamini, D., Reich, M. S., López-Mañas, R., Benetello, F., Collins, S. C., Bataille, C. P., Pierce, N. E., … Vila, R. (2023). The Afrotropical breeding grounds of the Palearctic-African migratory painted lady butterflies (Vanessa cardui). Proceedings of the National Academy of Sciences, 120(16), e2218280120. https://doi.org/10.1073/pnas.2218280120

# download package if missing
list.of.packages <- c("tidyverse","osfr", "assignR", "terra",  "sf", "geosphere", "tidyterra","spatialEco""ggspatial","raster")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load libraries
library(tidyverse)
library(osfr)
library(terra)
library(assignR)
library(sf)
library(raster)
library(geosphere)
library(spatialEco)
library(tidyterra)
library(ggspatial)

### projection
WGS84 <- "+init=epsg:4326"
laea<-"+proj=laea +lat_0=13.88 +lon_0=22.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" #equal area

# countries outline
extent3<-st_read("isotope_geolocation/extent3.shp") # from rnaturalearth

# preliminary version of Table S4
unknown <- read.csv("isotope_geolocation/isotope_data.csv")

# select long-distance individuals
unknown2<- unknown %>%
  filter(status=="long-distance")

# download isoscapes from OSF (Reich et al 2024)
proj <- osf_retrieve_node("dzh2g")
dir.create("isotope_geolocation/data_osf_dzh2g")
proj %>%
  osf_ls_files() %>%
  filter(name=="isoscapes") %>%
  osf_download(path="isotope_geolocation/data_osf_dzh2g", conflicts="overwrite")

# import isoscapes
rf_sr<-rast("isotope_geolocation/data_osf_dzh2g/isoscapes/Sr_mlr_mean_c.tif")
rf_sr_err<-rast("isotope_geolocation/data_osf_dzh2g/isoscapes/Sr_mlr_sd001_c.tif")
d2hw_mean<-rast("isotope_geolocation/data_osf_dzh2g/isoscapes/d2h_wing_mean_c.tif")
d2hw_sd<-rast("isotope_geolocation/data_osf_dzh2g/isoscapes/d2h_wing_sd_c.tif")

#august species distribution model from Talavera et al., 2023
Aug<-terra::rast("isotope_geolocation/occRep_ensemble_TSSbin_08.asc")
tot <- terra::global(Aug, sum, na.rm=T)
Aug_P<- Aug/tot[[1]]
Aug_P_cm<-terra::project(Aug_P,rf_sr,method="bilinear")
Aug_P_cm <- mask(Aug_P_cm, rf_sr)

#match the isoscapes
rf_sr <- mask(rf_sr, Aug_P_cm)
rf_sr_err <- mask(rf_sr_err, Aug_P_cm)
d2hw_mean <- mask(d2hw_mean, Aug_P_cm)
d2hw_sd <- mask(d2hw_sd, Aug_P_cm)

#stack isoscapes for geographic assignment
d2Hw_iso<-c(d2hw_mean, d2hw_sd)
sr_iso<-c(rf_sr, rf_sr_err)
iso<-assignR::isoStack(d2Hw_iso,sr_iso)

#prepare isotope data for assignR
iso.data_iso<-unknown2 %>%
  dplyr::select("Sample_ID","Delta2H_permil","X87_86_Kr_Rb_corr_xc") %>%
  as.data.frame()

#geographic assignment
dir.create("isotope_geolocation/output")
d_assign<-assignR::pdRaster(iso, unknown = iso.data_iso,  outDir = "isotope_geolocation/output", prior = Aug_P_cm, genplot=TRUE)

#average map
avg <- mean(d_assign, na.rm=T)

#capture location
indxy_laea <- unknown2 %>%
  sf::st_as_sf(coords = c("Lon", "Lat"), crs=WGS84) %>% #WGS84
  st_transform(indxy, crs = laea)

# plot Figure 3A ####
avg <- terra::setMinMax(avg, force =T)
Uniques<-terra::minmax(avg)
Uniques.max <-Uniques[[2]]
Uniques.min <-Uniques[[1]]
avg_plot <-ggplot()+
  geom_spatraster(data = avg) +
  scale_fill_gradient2(limits = c(Uniques.min, Uniques.max), low = "#F2F2F2",mid= "#FAE08C", high = "#3B4CC0", midpoint=((Uniques.max-Uniques.min)/2),
                       na.value = NA, breaks = seq(Uniques.min, Uniques.max, length.out = 5), labels = scales::label_scientific(digits =2)) +
  labs(fill= "Average probability")+
  geom_sf(data=extent3,color="grey30",fill="transparent", lwd=0.3)+ #countries
  geom_sf(data=indxy_laea, col="#E63D3A",fill="transparent", pch=21, cex=1.3)+
  annotation_scale(
    location = "bl",         # top-left
    width_hint = 0.2,        # controls width of scale bar
    unit_category = "metric",
    bar_cols = c("grey60", "white"),
    text_cex = 0.6,
    line_width = 0.3 ,
    text_col = "grey30"
  )+
  coord_sf(expand = FALSE)+
  theme_void()+
  theme(plot.margin = margin(0, 0, .2, .2, "cm"),
        legend.margin=margin(0,.2,.2,0),
        plot.background = element_rect(fill=NA, color=NA),
        panel.border=element_rect(fill=NA, color="grey30"),
        legend.background = element_rect(fill = NA, color=NA),
        legend.key.height = unit(1, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(plot=avg_plot, filename=paste0("isotope_geolocation/output/epi_spdist_avg_",nrow(unknown2),".png"), device="png", height=4, width=4, units="in",bg="white")
ggsave(plot=avg_plot, filename=paste0("isotope_geolocation/output/epi_spdist_avg_",nrow(unknown2),".pdf"), device="pdf", height=4, width=4, units="in",bg="transparent")

# combined rose plot modified from assignR
plot.wDist_MR_combined = function(x, ..., bin = 20, pty = "bear", indices = NULL, wedge_color = "blue"){
  # Check if x is of class "wDist"
  if(!inherits(x, "wDist")){
    stop("x must be a wDist object")
  }

  # Check if x is empty
  n = length(x)
  if(n == 0){
    stop("x is empty")
  }

  # Determine the number of indices to plot automatically if indices are not provided
  if (is.null(indices)) {
    indices <- 1:n
  }

  # Ensure indices are within bounds
  valid_indices <- indices[indices >= 1 & indices <= n]
  if(length(valid_indices) == 0){
    stop("No valid indices selected")
  }

  # Determine the number of indices to plot
  np = length(valid_indices)

  # Check if bin is numeric, within valid range, and a factor of 360
  if(!is.numeric(bin)){
    stop("bin must be numeric")
  }
  if(length(bin) > 1){
    stop("bin must be length 1")
  }
  if(bin <=0 | bin > 90){
    stop("bin must be a value between 0 and 90")
  }
  if(360 %% bin != 0){
    stop("bin should be a factor of 360")
  }

  # Check if pty is a valid value
  if(!(pty %in% c("both", "dist", "bear"))){
    stop("pty not valid for plot.xist")
  }

  # Combine bearing data and density values for selected indices
  combined_bearing_data <- numeric(0)
  combined_density_values <- numeric(0)

  for(i in valid_indices){
    combined_bearing_data <- c(combined_bearing_data, x[[i]]$b.dens$x)
    combined_density_values <- c(combined_density_values, x[[i]]$b.dens$y)
  }

  # Create a single combined rose plot
  opar = par(no.readonly = TRUE)
  on.exit(par(opar))

  if(pty %in% c("both", "bear")){
    # Define functions for creating arcs and wedges
    arc = function(a1, a2, b){
      a = seq(a1, a2, by = 0.5)
      r = 2 * pi * a / 360
      x = sin(r) * b
      y = cos(r) * b
      return(cbind(x, y))
    }

    wedge = function(a1, a2, b){
      xy = arc(a1, a2, b)
      xy = rbind(c(0,0), xy, c(0,0))
      return(xy)
    }

    # Create bins for bearing values
    bins = seq(-180, 179.9, by = bin)
    vals = numeric(length(bins))

    # Set up the layout of the combined plot
    mfr = 1
    mfc = 1
    par(mfrow = c(mfr, mfc), mar = c(1,1,2,1))

    # Ensure bearing values are within -180 to 180 degrees
    for(j in seq_along(combined_bearing_data)){
      if(combined_bearing_data[j] < -180){
        combined_bearing_data[j] = combined_bearing_data[j] + 360
      } else if(combined_bearing_data[j] >= 180){
        combined_bearing_data[j] = combined_bearing_data[j] - 360
      }
    }

    # Calculate values for each bin
    for(j in seq_along(bins)){
      vals[j] = sum(combined_density_values[combined_bearing_data >= bins[j] & combined_bearing_data < bins[j] + bin])
    }

    # Determine the maximum value for scaling
    b.max = max(vals)

    # Create and customize the plot
    xy = arc(-180, 180, b.max)
    plot(xy, type = "l", col = "dark grey", axes = FALSE,
         ylim = c(-b.max, 1.05 * b.max),
         xlab = "", ylab = "", asp = 1, main = "")
    lines(arc(-180, 180, b.max/2), col = "dark grey")
    for(j in c(-180, -90, 0, 90)){
      lines(wedge(j, j, b.max * 1.05), col = "dark grey")
    }

    # Fill the wedges with user-specified color
    for(j in seq_along(bins)){
      xy = wedge(bins[j], bins[j] + bin, vals[j])
      polygon(xy, col = wedge_color)
    }
  }

  # Restore original plotting parameters
  par(opar)
  return() }

# plot Figure 3A inset ####
indxy_sv <- vect(cbind(unknown2$Lon,unknown2$Lat),atts=unknown2,crs="+proj=longlat +datum=WGS84 +no_defs") # Spatvect
indxy_sv_aea <- project(indxy_sv, y = laea)

listy <- paste0("isotope_geolocation/output/",unknown2$Sample_ID[1:nrow(unknown2)],"_like.tif")
stackyp <- rast(listy)
names(stackyp)<-unknown2$Sample_ID
bear<-assignR::wDist(stackyp, indxy_sv_aea)

#plot
png(file=paste0("isotope_geolocation/output/bearings.png"))
plot(bear, wedge_color="#3B4CC0", bin=10)
dev.off()

#plot
png(file=paste0("isotope_geolocation/output/bearing_combined.png"))
plot.wDist_MR_combined(bear, wedge_color="#3B4CC0", bin=10)
dev.off()

# plot individual maps (Figure S4) ####
i720<-rast("isotope_geolocation/output/18B720_like.tif")
#capture location
indxy_laea <- unknown2 %>%
  filter(Sample_ID=="18B720") %>%
  sf::st_as_sf(coords = c("Lon", "Lat"), crs=WGS84) %>% #WGS84
  st_transform(indxy, crs = laea)
#map
i720 <- terra::setMinMax(i720, force =T)
Uniques<-terra::minmax(i720)
Uniques.max <-Uniques[[2]]
Uniques.min <-Uniques[[1]]
i720_plot <-ggplot()+
  geom_spatraster(data = i720) +
  scale_fill_gradient2(limits = c(Uniques.min, Uniques.max), low = "#F2F2F2",mid= "#FAE08C", high = "#3B4CC0", midpoint=((Uniques.max-Uniques.min)/2),
                       na.value = NA, breaks = seq(Uniques.min, Uniques.max, length.out = 5), labels = scales::label_scientific(digits =2)) +
  labs(fill= "Probability")+
  geom_sf(data=extent3,color="grey30",fill="transparent", lwd=0.3)+ #countries
  geom_sf(data=indxy_laea, col="#E63D3A",fill="transparent", pch=21, cex=1.3)+
  coord_sf(expand = FALSE)+
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.margin=margin(0,0,.2,0), plot.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color="grey30"),legend.background = element_rect(fill = NA, color=NA),legend.key.height = unit(1, "cm"), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot=i720_plot, filename=paste0("isotope_geolocation/output/18B720_plot.png"), device="png", height=4, width=4, units="in",bg="white")

i721<-rast("isotope_geolocation/output/18B721_like.tif")
#capture location
indxy_laea <- unknown2 %>%
  filter(Sample_ID=="18B721") %>%
  sf::st_as_sf(coords = c("Lon", "Lat"), crs=WGS84) %>% #WGS84
  st_transform(indxy, crs = laea)
#map
i721 <- terra::setMinMax(i721, force =T)
Uniques<-terra::minmax(i721)
Uniques.max <-Uniques[[2]]
Uniques.min <-Uniques[[1]]
i721_plot <-ggplot()+
  geom_spatraster(data = i721) +
  scale_fill_gradient2(limits = c(Uniques.min, Uniques.max), low = "#F2F2F2",mid= "#FAE08C", high = "#3B4CC0", midpoint=((Uniques.max-Uniques.min)/2),
                       na.value = NA, breaks = seq(Uniques.min, Uniques.max, length.out = 5), labels = scales::label_scientific(digits =2)) +
  labs(fill= "Probability")+
  geom_sf(data=extent3,color="grey30",fill="transparent", lwd=0.3)+ #countries
  geom_sf(data=indxy_laea, col="#E63D3A",fill="transparent", pch=21, cex=1.3)+
  coord_sf(expand = FALSE)+
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.margin=margin(0,0,.2,0), plot.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color="grey30"),legend.background = element_rect(fill = NA, color=NA),legend.key.height = unit(1, "cm"), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot=i721_plot, filename=paste0("isotope_geolocation/output/18B721_plot.png"), device="png", height=4, width=4, units="in",bg="white")


i773<-rast("isotope_geolocation/output/18B773_like.tif")
#capture location
indxy_laea <- unknown2 %>%
  filter(Sample_ID=="18B773") %>%
  sf::st_as_sf(coords = c("Lon", "Lat"), crs=WGS84) %>% #WGS84
  st_transform(indxy, crs = laea)
#map
i773 <- terra::setMinMax(i773, force =T)
Uniques<-terra::minmax(i773)
Uniques.max <-Uniques[[2]]
Uniques.min <-Uniques[[1]]
i773_plot <-ggplot()+
  geom_spatraster(data = i773) +
  scale_fill_gradient2(limits = c(Uniques.min, Uniques.max), low = "#F2F2F2",mid= "#FAE08C", high = "#3B4CC0", midpoint=((Uniques.max-Uniques.min)/2),
                       na.value = NA, breaks = seq(Uniques.min, Uniques.max, length.out = 5), labels = scales::label_scientific(digits =2)) +
  labs(fill= "Probability")+
  geom_sf(data=extent3,color="grey30",fill="transparent", lwd=0.3)+ #countries
  geom_sf(data=indxy_laea, col="#E63D3A",fill="transparent", pch=21, cex=1.3)+
  coord_sf(expand = FALSE)+
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.margin=margin(0,0,.2,0), plot.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color="grey30"),legend.background = element_rect(fill = NA, color=NA),legend.key.height = unit(1, "cm"), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot=i773_plot, filename=paste0("isotope_geolocation/output/18B773_plot.png"), device="png", height=4, width=4, units="in",bg="white")

i116<-rast("isotope_geolocation/output/19H116_like.tif")
#capture location
indxy_laea <- unknown2 %>%
  filter(Sample_ID=="19H116") %>%
  sf::st_as_sf(coords = c("Lon", "Lat"), crs=WGS84) %>% #WGS84
  st_transform(indxy, crs = laea)
#map
i116 <- terra::setMinMax(i116, force =T)
Uniques<-terra::minmax(i116)
Uniques.max <-Uniques[[2]]
Uniques.min <-Uniques[[1]]
i116_plot <-ggplot()+
  geom_spatraster(data = i116) +
  scale_fill_gradient2(limits = c(Uniques.min, Uniques.max), low = "#F2F2F2",mid= "#FAE08C", high = "#3B4CC0", midpoint=((Uniques.max-Uniques.min)/2),
                       na.value = NA, breaks = seq(Uniques.min, Uniques.max, length.out = 5), labels = scales::label_scientific(digits =2)) +
  labs(fill= "Probability")+
  geom_sf(data=extent3,color="grey30",fill="transparent", lwd=0.3)+ #countries
  geom_sf(data=indxy_laea, col="#E63D3A",fill="transparent", pch=21, cex=1.3)+
  coord_sf(expand = FALSE)+
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.margin=margin(0,0,.2,0), plot.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color="grey30"),legend.background = element_rect(fill = NA, color=NA),legend.key.height = unit(1, "cm"), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot=i116_plot, filename=paste0("isotope_geolocation/output/19H116_plot.png"), device="png", height=4, width=4, units="in",bg="white")

i128<-rast("isotope_geolocation/output/19H128_like.tif")
#capture location
indxy_laea <- unknown2 %>%
  filter(Sample_ID=="19H128") %>%
  sf::st_as_sf(coords = c("Lon", "Lat"), crs=WGS84) %>% #WGS84
  st_transform(indxy, crs = laea)
#map
i128 <- terra::setMinMax(i128, force =T)
Uniques<-terra::minmax(i128)
Uniques.max <-Uniques[[2]]
Uniques.min <-Uniques[[1]]
i128_plot <-ggplot()+
  geom_spatraster(data = i128) +
  scale_fill_gradient2(limits = c(Uniques.min, Uniques.max), low = "#F2F2F2",mid= "#FAE08C", high = "#3B4CC0", midpoint=((Uniques.max-Uniques.min)/2),
                       na.value = NA, breaks = seq(Uniques.min, Uniques.max, length.out = 5), labels = scales::label_scientific(digits =2)) +
  labs(fill= "Probability")+
  geom_sf(data=extent3,color="grey30",fill="transparent", lwd=0.3)+ #countries
  geom_sf(data=indxy_laea, col="#E63D3A",fill="transparent", pch=21, cex=1.3)+
  annotation_scale(location = "bl",         # top-left
    width_hint = 0.2,        # controls width of scale bar
    unit_category = "metric", bar_cols = c("grey60", "white"), text_cex = 0.6,
    line_width = 0.3 , text_col = "grey30")+
  coord_sf(expand = FALSE)+ theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.margin=margin(0,0,.2,0), plot.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color="grey30"),legend.background = element_rect(fill = NA, color=NA),legend.key.height = unit(1, "cm"), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot=i128_plot, filename=paste0("isotope_geolocation/output/19H128_plot.png"), device="png", height=4, width=4, units="in",bg="white")

# distance metrics ####
metrics <- data.frame(ID=NA,capx=NA,capy=NA, HSr_local=NA,HSr_distmin=NA, HSr_distminy=NA, HSr_distminy=NA, HSr_distmax=NA, HSr_distcentx=NA, HSr_distcenty=NA, HSr_area=NA,HSr_per=NA)[numeric(0), ]

for(k in 1:nrow(unknown2)){
  ind<-unknown2[k,]
  #location of capture
  capxy<-ind
  capxy <- sf::st_as_sf(capxy, coords = c("Lon", "Lat"), crs=WGS84) #WGS84
  capxy_eck <- sf::st_transform(capxy, crs = laea)
  capxy_sp <- as(capxy, "Spatial")
  #bring in raster
  x<-rast(paste("isotope_geolocation/output/",ind$Sample_ID[1],"_like.tif",sep=""))
  x333<-qtlRaster(x, .333, thresholdType="prob",genplot=F)
  Tarea<-ncell(x333)-freq(x333, value=NA) #land cells
  HSr_area <- zonal(x333, x333, fun='sum') #square kilometers - equal area projection so this method should be okay
  local<-terra::extract(x333, capxy_eck, xy=TRUE,method="simple") #is the capture location within the 2:1 odds ratio threshold? layer: TRUE=local, FALSE=migrant
  #binary
  x333_bi <- x333*1
  x333_bi[x333_bi <= 0] <- NA
  #sieve out single cells
  patch <- terra::patches(x333_bi, directions = 8)
  patchFreq <- freq(patch) #frequency table
  excludeID <- patchFreq$value[which(patchFreq$count == 1)]
  if (length(excludeID) > 0) {x333_bi[patch %in% excludeID] <- NA}
  #polygon
  lin <- as.polygons(x333_bi)
  #coordinates of centroid
  cent <- centroids(lin, inside=FALSE) #center does not have to be in the polygon
  cent_df<-geom(cent) #get coordinates in laea
  #minimum distance from capture point to polygon
  pol <- terra::project(lin, capxy)
  pol_sf <- sf::st_as_sf(pol)
  pol_sp <- as(pol, "Spatial")
  dist_min <- geosphere::dist2Line(p = c(ind$Lon[1], ind$Lat[1]), line = pol_sp, distfun=distHaversine)
  #get vertices of polygon
  vertices <- spatialEco::extract.vertices(x = pol_sf, join=F)
  vertices_sf <- sf::st_as_sf(vertices)
  st_crs(vertices_sf) <- 4326
  #maximum distance
  max_dist <- st_distance(vertices_sf, capxy, by_element = F)
  maximumd <- as.numeric(max(max_dist)/1000)
  #tally
  indi <- data.frame(ID=ind$Sample_ID,
                     capx=ind$Lon,
                     capy=ind$Lat, #coord of capture site in WGS84
                     HSr_local=local$layer, #is capture location in 2:1 odds ratio area?
                     HSr_distmin=dist_min[1]/1000,
                     HSr_distminx=dist_min[2],
                     HSr_distminy=dist_min[3],
                     HSr_distmax=maximumd,
                     HSr_distcentx=cent_df[3],
                     HSr_distcenty=cent_df[4],
                     HSr_area=HSr_area$layer.1[2],#area of 2:1 odds ratio polygon
                     HSr_per= (HSr_area$layer.1[2]/Tarea$count*100)  )
  metrics<-rbind(metrics,indi)
  print(paste("finished", k, "of", nrow(unknown2)))
}

# Table S4 ####
tableS4 <-unknown %>%
left_join(metrics %>%
    mutate(HSr_distmin = round(HSr_distmin, -1), HSr_distmax = round(HSr_distmax, -1)) %>%
      dplyr::select(ID, HSr_distmin, HSr_distmax), by = c("Sample_ID" = "ID"))
print(tableS4)
write.csv(metrics, file = "isotope_geolocation/output/TableS4.csv")