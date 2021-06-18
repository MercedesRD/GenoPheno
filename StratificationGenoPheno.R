#############################################################################################################################################
###  Preparation of a layer describing the degree of human modification of pre-1750 vegetation
###  Similar to the VAST framework proposed by Keith and Simpson (2008)
###
###  In this script: Align and resample raster layers, reclassify into categories with different degree of expected human pressure on soils

######## 1. Align, resample, crop

### Desired extent: NSW
### Resolution: 30 m
### CRS: EPSG=4326
### Align rasters with the DEM 30m

###  Author: Mercedes Roman Dobarco
###  Date: 31/10/2019

####### Load packages
library(rgdal)
library(sp)
library(sf)
library(gstat)
library(raster)
library(lattice)
library(ggplot2)
library(dplyr)
library(gdalUtils)
library(foreach)
library(parallel)
library(doParallel)
library(snow)
library(doSNOW)
library(rasterVis)
library(viridis)
library(scales)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggplot2) # tidyverse data visualization package
library(viridis) # color palettes
library(scales)
library(rasterVis)
library(lattice)
library(gridExtra)
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(mapview) # for interactive maps
library(shiny)   # for web applications
library(RColorBrewer)


### Set working directory
### Set the home directory where all files and subfolders are stored
HomeDir <- "C:/Covariates/"  ### Change this with your Home directory 

# ### 1. Get the boundary for NSW analysis --------------------------------
setwd(paste0(HomeDir,"Administrative"))
### Look into the shapefile for NSW
nsw <- st_read("nsw.shp")
nsw_buffer <- st_read(paste0(HomeDir,"Administrative/", "nsw_buffer.shp"))
nsw_buff_geom <- st_geometry(nsw_buffer)
st_bbox(nsw_buffer) ## This is what we will be using, yep.
st_bbox(nsw) 

#### Function to load rasters
load_raster <- function (x) {
  maps <- list()
  for (rast in 1:length(x)) {  
    maps[[rast]] <- raster(x[rast])
  }
  return(maps)
}

# ### 2. Align, resample, project the vegetation layers -----------------------

###  NSW Native Vegetation Extent 5m Raster v1.2
# Classes for native vegetation
# .	Tree cover- 1
# .	Woodland matrix - 5
# .	Candidate native grasslands - 2
# .	Forestry plantations - 3
# .	Non-native areas - 0
# .	Water - 4
setwd(paste0(HomeDir,"Vegetation/nswnativevegetationextentv1p25m2017"))
NatVeg5m <- raster("NSW_Native_Vegetation_Extent_v1p2_5m_2017.tif")
gdalwarp(srcfile = paste0(HomeDir,"Vegetation/nswnativevegetationextentv1p25m2017/NSW_Native_Vegetation_Extent_v1p2_5m_2017.tif"),
         dstfile = gsub(".tif","_r30m.tif", 
                        paste0(HomeDir,"Vegetation/nswnativevegetationextentv1p25m2017/NSW_Native_Vegetation_Extent_v1p2_5m_2017.tif")),
         s_srs = '+proj=lcc +lat_1=-30.75 +lat_2=-35.75 +lat_0=-33.25 +lon_0=147 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +units=m +no_defs', 
         t_srs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0',
         te = c(140.916, -37.58847, 153.7221, -28.07375),
         tr=c(0.0002777778,0.0002777778), ## 30m in decimal degrees
         r="mode", ## Because it is a categorical variable and I will aggregate different pixels, I choose the most frequest
         overwrite=TRUE,
         verbose=TRUE)
rm(NatVeg5m)

### CLUM Update 12/2018 at 50m
setwd(paste0(HomeDir,"LandUse/geotiff_clum_50m1218m"))
clum50 <- raster("clum_50m1218m.tif")
gdal_setInstallation()
gdalwarp(srcfile = paste0(HomeDir,"LandUse/geotiff_clum_50m1218m/clum_50m1218m.tif"),
         dstfile = gsub("clum_50m1218m.tif","clum_30m1218.tif", 
                        paste0(HomeDir,"LandUse/geotiff_clum_50m1218m/clum_50m1218m.tif")),
         s_srs = '+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 
         t_srs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0',
         te = c(140.916, -37.58847, 153.7221, -28.07375),
         tr=c(0.0002777778,0.0002777778), ## 30m in decimal degrees
         r="near", ## Because it is a categorical variable and I will reduce from 50 to 30 m resolution.
         overwrite=TRUE,
         verbose=TRUE)
rm(clum50)


### Extant native vegetation (Keith and Simpson, 2006).
### Binary map which discriminates between remnant native vegetation (intact grasslands and woody vegetation)
### and cleared vegetation (includes non-native and secondary grasslands of native vegetation).
### the layer that interests me most is the secondary grassland layer, which I will add to the mask of non-intact (or non-native) vegetation at 5m
setwd(paste0(HomeDir,"Vegetation/NSWExtantNativeVegetationV2/"))
sgrass002 <- raster("sgrass002.tif")
gdalwarp(srcfile = paste0(HomeDir,"Vegetation/NSWExtantNativeVegetationV2/sgrass002.tif"),
         dstfile = gsub("sgrass002.tif","sgrass30m.tif", 
                        paste0(HomeDir,"Vegetation/NSWExtantNativeVegetationV2/sgrass002.tif")),
         s_srs = '+proj=lcc +lat_1=-30.75 +lat_2=-35.75 +lat_0=-33.25 +lon_0=147 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 
         t_srs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0',
         te = c(140.916, -37.58847, 153.7221, -28.07375),
         tr=c(0.0002777778,0.0002777778), ## 30m in decimal degrees
         r="near", ## Because it is a categorical variable and I will reduce from 250 to 30 m resolution.
         overwrite=TRUE,
         verbose=TRUE)
rm(sgrass002)


### but I am also thinking of the extant native vegetation, layer extveg002.tif
setwd(paste0(HomeDir,"Vegetation/NSWExtantNativeVegetationV2/"))
extveg002 <- raster("extveg002.tif")
gdalwarp(srcfile = paste0(HomeDir,"Vegetation/NSWExtantNativeVegetationV2/extveg002.tif"),
         dstfile = gsub("extveg002.tif","extveg30m.tif", 
                        paste0(HomeDir,"Vegetation/NSWExtantNativeVegetationV2/extveg002.tif")),
         s_srs = '+proj=lcc +lat_1=-30.75 +lat_2=-35.75 +lat_0=-33.25 +lon_0=147 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 
         t_srs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0',
         te = c(140.916, -37.58847, 153.7221, -28.07375),
         tr=c(0.0002777778,0.0002777778), ## 30m in decimal degrees
         r="near", ## Because it is a categorical variable and I will reduce from 250 to 30 m resolution.
         overwrite=TRUE,
         verbose=TRUE)
rm(extveg002)



#### Vegetation structure @ 30 m  - Forest structure categories (plant cover fraction and canopy height)
setwd(paste0(HomeDir,"VegetationStructure/"))
VegSt2009 <- raster("alpsbk_aust_y2009_sf1a2.tif")
gdalwarp(srcfile = paste0(HomeDir,"VegetationStructure/alpsbk_aust_y2009_sf1a2.tif"),
         dstfile = paste0(HomeDir,"VegetationStructure/VegStr2009.tif"),
         s_srs = '+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs', 
         t_srs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0',
         te = c(140.916, -37.58847, 153.7221, -28.07375),
         tr=c(0.0002777778,0.0002777778), ## 30m in decimal degrees
         r="near", ## Because it is a categorical variable.
         overwrite=TRUE,
         verbose=TRUE)
rm(VegSt2009)


### Load resampled layers
nsw_buffer <- st_read( paste0(HomeDir,"Administrative/", "nsw_buffer.shp")) # read polygon

setwd(paste0(HomeDir,"Vegetation/nswnativevegetationextentv1p25m2017"))
NatVeg30m <- raster("NSW_Native_Vegetation_Extent_v1p2_5m_2017_r30m.tif")
setwd(paste0(HomeDir,"LandUse/geotiff_clum_50m1218m"))
clum30m <- raster("clum_30m1218.tif")
setwd(paste0(HomeDir,"Vegetation/NSWExtantNativeVegetationV2/"))
extveg30m <- raster("extveg30m.tif")
sgrass30m <- raster("sgrass30m.tif")
setwd(paste0(HomeDir,"VegetationStructure/"))
VegStr2009 <- raster("VegStr2009.tif")
vegxlu <- stack(NatVeg30m,clum30m,extveg30m,VegStr2009)
names(vegxlu) <- c("Veg_Ext","clum_30m", "Native", "VegStr")

plot(vegxlu)
#vegxlu <-mask(vegxlu, mask = nsw_buffer)

par(mfrow=c(2,2), las=1,xaxt='n',yaxt="n",c(5, 4, 4, 4)+0.1)
### Plot extent native vegetation
plot(NatVeg30m,breaks=c(-1,0,1,2,3,4,5),
     col=c("grey","#1F9E89FF","#FDE725FF","grey","cadetblue1","#6DCD59FF"),legend=NULL)

plot(extveg30m, breaks=c(-1,0,1),
     col=c("grey","#31688EFF"),legend=NULL)

plot(mask(clum30m, nsw),
     col=hcl.colors(178, palette = "TealRose", alpha = NULL, rev = TRUE, fixup = TRUE), legend=NULL)

plot(mask(VegStr2009, nsw))

# setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/VASTlayer")
# GenoPheno <- raster("GenoPheno.tif")
# plot(GenoPheno, breaks=c(0,1,2,5,7,10,31,32,33,41,42,43,54,99),
#      col=c("#2A788EFF","#1F988BFF","#E34E65FF","#A6317DFF","#FED99BFF",
#            "#FDE725FF","#D2E21BFF","#A5DB36FF","#7AD151FF","#54C568FF","#35B779FF",
#            "#C63C74FF","#440154FF"),legend=NULL)


#################################################################################################################################

# ### 5. Reclassify vegetation x Land Use into VAST-LU classes ------------
setwd(paste0(HomeDir,"Vegetation/NSWExtantNativeVegetationV2/"))

### Create a layer that separates potentially native vegetation from non-native (and water)
NonNat_reclass <- function(s,...){
  ### First, let's create an empty vector and fill the values
  NonNat <- rep(NA, length(s[1]))
  ###  NSW Native Vegetation Extent 5m Raster v1.2
  # Classes for native vegetation
  # .	Tree cover- 1
  # .	Woodland matrix - 5
  # .	Candidate native grasslands - 2
  # .	Forestry plantations - 3
  # .	Non-native areas - 0
  # .	Water - 4
  ### If there is a missing value for Native vegetation extent, assign NA
  ### I also want to exclude from the analysis those areas that are water bodies
  NonNat <- ifelse(test=(s[1] == 0 | s[1]==3),yes = 0,
                   ifelse(test=(s[1] == 1 | s[1]==2| s[1]==5),yes = 1,
                          ifelse(test = s[1]==4,4,NA)))
  ### NonNat Classes
  ### Non-native and forestry plantations - 0
  ### Potentially native - 1
  ### Water bodies - 4
  return(NonNat)
}

ff <- function(x) calc(x, NonNat_reclass)
beginCluster(7)
NonNat2017 <- clusterR(vegxlu, ff, export = list("NonNat_reclass"),
                       filename = "NonNat2017.tif", format = "GTiff",
                       na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(NonNat2017)

### replace NatVeg30m by this new layer better
#vegxlu <- stack(NonNat2017, extveg30m, clum30m, VegStr2009)

### Differentiate between areas that have been cleared, and areas that have not been cleared
### Cleared vegetation may be native, however, it may be the result of intensive grazing, burning, clearing forests, etc.
### Cross potentially native vegetation with the extant vegetation layer by Keith and Simpson, and after, with Land Use

### Cleared native, remnant native (potentially), 
Cleared_reclass <- function(s,...){
  ### First, let's create an empty vector and fill the values
  cleared <- rep(NA, length(s[1]))
  ### NonNat2017 Classes s[1]
  #   Non-native and forestry plantations - 0
  #   Potentially native - 1
  #   Water bodies - 4
  ### extveg30m classes s[2]
  # .	Intact- 1
  # .	Cleared - 0
  
  ### If there is a missing value for Native vegetation extent, assign NA
  ### I also want to exclude from the analysis those areas that are water bodies
  cleared <- ifelse(test=(s[1]==4), yes=4,                               # Keep a class for water bodies
                    ifelse(test=(s[1]==1 & s[2]==1), yes = 1,            # Potential genosoils, intact and native vegetation.
                           ifelse(test=(s[1]==1 & s[2]==0), yes = 2,     
                                  # Potential native vegetation, cleared. At least modified or transformed vegetation, phenosoils.
                                  ifelse(test=(s[1]==1 & is.na(s[2])), yes=1,
                                  # Potential native vegetation, where there is no data for cleared or non-cleared, will be kept as potential genosoil
                                  ifelse(s[1]==0, yes=3, NA))))) # Non-native areas will be later classified into phenosoils.
  return(cleared)
  ### VegxCleared Classes
  ### Potentially native remnant - 1
  ### Native cleared - 2
  ### Non-native and forestry plantations - 3
  ### water bodies - 4
}

ff <- function(x) calc(x, Cleared_reclass)
beginCluster(7)
VegxCleared <- clusterR(vegxlu, ff, export = list("Cleared_reclass"),
                       filename = "VegxCleared2017.tif", format = "GTiff",
                       na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

#setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/VASTlayer")
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/VASTlayer")
VegxCleared <- raster("VegxCleared2017.tif")
plot(VegxCleared, breaks =c(0,1,2,3,4),
     col=c("mediumseagreen","olivedrab2","mediumvioletred","cadetblue1"))

library(raster)
library(rasterVis)
library(mapview)
library(RColorBrewer)
my_rst <- ratify(VegxCleared)

# Add some custom labels for each class. Add two kinds to show that
# rasterVis::levelplot() can select one of them.
levels(my_rst)[[1]]$PotGeno <- c("Native remnant", "Native cleared", "Non-native", "Water bodies")
levels(my_rst)
# [[1]]
# ID    PotGeno
# 1  1     NatRem
# 2  2 NatCleared
# 3  3  NonNative
# 4  4      Water

# Custom palette
my_palette <- c("mediumseagreen","olivedrab2","gold3","cadetblue1")

# Plot with rasterVis::levelplot(). Can choose a variable/column from the
# 'Raster Attribute Table' (RAT) with the `att` argument; att = 1L by default.
# see https://oscarperpinan.github.io/rastervis/#factor
levelplot(my_rst, col.regions = my_palette)

# Plot with mapview::mapView(). Only the integers from ID column of RAT are
# displayed.
mapView(my_rst, alpha.regions=0.6)

## Add to cluster
#vegxlu <- stack(NonNat2017, extveg30m, clum30m, VegStr2009)
vegxlu <- stack(VegxCleared, clum30m)
plot(vegxlu)
#vegxlu <- stack(vegxlu,HMcNSW30m)


#####################################################################################################################################

### what are the unique combinations of the raster stack?
u <- unique(vegxlu,na.last=TRUE, progress="text")
u.2 <- u;rm(u)
lu <- unique(clum30m,na.last=TRUE, progress="text")

### Intensive uses, 
### Intensive horticulture, intensive animal production, industrial or manufacturing
### imply aprofound change in the original vegetation, 
d <- as.data.frame(u)
d.2 <- as.data.frame(u.2)
colnames(d) <- c("NatVeg30m2017","clum30m2018","extveg30m2008","VegStr30m2009")

### I ignore vegetation structure for now
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/VASTlayer")

### Let's focus on cases with potentially native non-cleared vegetation.
### Potential GENOSOILS
potGenosoils <- filter(.data=d,
       !is.na(NatVeg30m2017),                                                ## Exclude areas where there are NA for extant native vegetation
       !is.na(clum30m2018),                                                  ## Exclude pixels with NA for land use
       extveg30m2008 == 1,                                                   ## Intact vegetation
       NatVeg30m2017 == 1 | NatVeg30m2017 == 2 | NatVeg30m2017 == 5) %>%     ## Potentially native vegetation 
distinct(clum30m2018,.keep_all = TRUE) %>% arrange(clum30m2018)
write.csv(potGenosoils, file = "NativeRemnantLU.csv")

PhenosoilsLow <- filter(.data=d,
                       !is.na(NatVeg30m2017),                             ## Exclude areas where there are NA for extant native vegetation
                       !is.na(clum30m2018),                               ## Exclude pixels with NA for land use
                       extveg30m2008 == 0,                                ## Cleared vegetation
                       NatVeg30m2017 == 1 | NatVeg30m2017 == 2 | NatVeg30m2017 == 5) %>%     ## Potentially native vegetation 
  distinct(clum30m2018,.keep_all = TRUE) %>% arrange(clum30m2018)
write.csv(PhenosoilsLow, file = "NativeClearedLU.csv")

NonNativeCleared <- filter(.data=d,
                        !is.na(NatVeg30m2017),                             ## Exclude areas where there are NA for extant native vegetation
                        !is.na(clum30m2018),                               ## Exclude pixels with NA for land use
                        extveg30m2008 == 0,                                ## Cleared vegetation
                        NatVeg30m2017 == 0 | NatVeg30m2017 == 3) %>%       ## Non-native and forestry plantations
  distinct(clum30m2018,.keep_all = TRUE) %>% arrange(clum30m2018)
write.csv(NonNativeCleared, file = "NonNativeClearedLU.csv")


##################################################################################################################

### Potential REMNANT PEDOGENONS
u <- unique(vegxlu, na.last=TRUE, progress="text")

potGenosoils.2 <- filter(.data=d.2,
                       !is.na(VegxCleared),                                             ## Exclude areas where there are NA for extant native vegetation
                       !is.na(clum_30m1218),                                                 ## Exclude pixels with NA for land use
                       VegxCleared2017 == 1) %>%                                            ## Potentially native intact vegetation 
  distinct(clum_30m1218,.keep_all = TRUE) %>% arrange(clum_30m1218)

PhenosoilsLow.2 <- filter(.data=d.2,
                         !is.na(VegxCleared2017),                                             ## Exclude areas where there are NA for extant native vegetation
                         !is.na(clum_30m1218),                                                 ## Exclude pixels with NA for land use
                         VegxCleared2017 == 2) %>%                                            ## Potentially native intact vegetation 
  distinct(clum_30m1218,.keep_all = TRUE) %>% arrange(clum_30m1218)


### LU classes present

## I don't want to consider intensive uses because they will be automatically clasiffied as areas where the vegetation has been removed

### Create a layer that separates potentially native vegetation from non-native (and water)
VegxLU_reclass <- function(s,...){
  
  ### First, let's create an empty vector and fill the values with 99 to signal anomalies
  GenPhen <- rep(99, length(s[1]))
  ### VegxCleared Classes
  ### Potentially native remnant - 1
  ### Native cleared - 2
  ### Non-native and forestry plantations - 3
  ### water bodies - 4
  
  ### If there is a missing value for Native vegetation extent or Land Use assign NA
  ### I also want to exclude from the analysis those areas that are water bodies
  GenPhen[is.na(s[1])|s[1]==4|is.na(s[2])] <- NA 
  # GenPhen[is.na(s[1])|is.na(s[2])] <- NA 
  
  ### If the vegetation is potentially native and remnant, use the first classification I made
  Genosoil.Remnant <- c(110:117,120:125,130,131,133,
                        610,611,614, # Lake, lake conservation, lake saline
                        630,631, # River, river conservation
                        650,651,654, # Marsh conservation and saline
                        660,661)
  Genosoil.II <- c(132, # Stock route
                   200:222, # Near natural production
                   314,414, # Environmental forest plantation
                   612,632,652,662) # Lake production, river production, marsh production, estuary production
  Phenosoil.Forestry <- c(310:313,410:413)
  Phenosoil.Dry.Grazing <- c(320:325)
  Phenosoil.Dry.Cropping <- c(330:353,360,362:365,542) # including rural residential with agriculture
  Phenosoil.Irrigated.Grazing <- c(420:424)
  Phenosoil.Irrigated.Cropping <- c(430:460,462:465)
  Phenosoil.Degraded <- c(361,461)
  Phenosoil.Recovery <- 134
  Exclude <- c(500:541,543:595,613,620:623,633,640:643,653,663) # Intensive uses + lake intensive, Dam/reservoir, river intensive, channel/aqueducts
  
  #Phenosoil.Residential <- c(540,541,543:545,553)
  # Phenosoil.Intensive <- c(510,511,512,513,515,520,521,522,523,524,525,526,527,528,530,531,532,533,534,535,537,
  #                          538,550,551,552,554,555,560,561,562,563,564,565,566,567,570,571,572,573,574,575,580,
  #                          581,582,583,584,590,591,592,593,595,
  #                          613,620:623,633,640:643,653) # lake intensive, Dam/reservoir, river intensive, channel/aqueducts
  
  ### Assign classes for Genosoils on native remnant
  GenPhen[s[1]==1 & s[2] %in% Genosoil.Remnant] <- 1 # Remnant Genosoils
  GenPhen[s[1]==1 & s[2] %in% Genosoil.II] <- 2 # Genosoil II, low modification
  
  ### Assign classes for Phenosoil Low in native cleared
  GenPhen[s[1]==2 & s[2] %in% Genosoil.Remnant] <- 3 # Phenosoil low modification
  GenPhen[s[1]==2 & s[2] %in% Genosoil.II] <- 3 # Phenosoil low modification
  
  ### Assign classes for Phenosoil Low in non-native
  GenPhen[s[1]==3 & s[2] %in% Genosoil.Remnant] <- 3 # Phenosoil low modification
  GenPhen[s[1]==3 & s[2] %in% Genosoil.II] <- 3 # Phenosoil low modification
  
  ### Assign Phenosoil classes
  GenPhen[s[1] ==3  & s[2] == 314] <- 4 # Phenosoil.Forestry
  GenPhen[s[1] ==3  & s[2] == 414] <- 4 # Phenosoil.Forestry
  GenPhen[s[2] %in% Phenosoil.Forestry] <- 4 # Phenosoil.Forestry
  GenPhen[s[2] %in% Phenosoil.Dry.Grazing]  <- 5 # Phenosoil.Dry.Grazing
  GenPhen[s[2] %in% Phenosoil.Dry.Cropping] <- 6 # Phenosoil.Dry.Cropping
  GenPhen[s[2] %in% Phenosoil.Irrigated.Grazing]  <- 7 # Phenosoil.Irrigated.Grazing
  GenPhen[s[2] %in% Phenosoil.Irrigated.Cropping] <- 8 # Phenosoil.Irrigated.Cropping
  GenPhen[s[2] %in% Phenosoil.Degraded] <- 9 # Phenosoil.Degraded
  GenPhen[s[2] %in% Phenosoil.Recovery] <- 10 # Phenosoil.Recovery
  GenPhen[s[2] %in% Exclude] <- NA # Phenosoil.Intensive - Excluded
  ## Return this value
  GenPhen
}

#### First version of genosoils/phenosoils classification
# GenoPheno <- calc(vegxlu, fun = VegxLU_reclass, 
#                   filename = "GenoPheno.032020.tif", format = "GTiff",
#                   na.rm=T, inf.rm=T,overwrite = T)

ff <- function(x) calc(x, VegxLU_reclass)
beginCluster(2)
GenoPheno <- clusterR(vegxlu, ff, export = list("VegxLU_reclass"),
                      filename = "GenoPheno.032020.tif", format = "GTiff",
                      na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
gc()
GenoPheno <- raster("GenoPheno.032020.tif")


###############################################################################################################
plot(GenoPheno)
GenoPheno.levels <- ratify(GenoPheno)

# Add some custom labels for each class. Add two kinds to show that
# rasterVis::levelplot() can select one of them.
levels(GenoPheno.levels)[[1]]$GenoPheno <- c("Remnant Genosoils", #1
                                             "Genosoil II", #2
                                             "Phenosoil Low", #3
                                             "Phenosoil Forestry", #4
                                             "Phenosoil Dry Grazing", #5 
                                             "Phenosoil Dry Cropping", #6
                                             "Phenosoil Irrigated Grazing", #7
                                             "Phenosoil Irrigated Cropping", #8
                                             "Phenosoil Degraded", #9
                                             "Phenosoil Rehabilitation", #10
                                             "Unclassified") #99
levels(GenoPheno.levels)

# Custom palette
my_palette <- c("mediumseagreen","olivedrab2","orchid1", 
                "blueviolet","darkkhaki", "gold1",
                "darkseagreen","darkcyan",
                "brown1","antiquewhite2","azure")
cols2 <- c("#2A788EFF","#1F988BFF","#E34E65FF","#A6317DFF","#FED99BFF",
           "#FDE725FF","#D2E21BFF","#A5DB36FF","#7AD151FF","#54C568FF","#35B779FF",
           "#C63C74FF","#440154FF")

HomeDir <- "C:/Covariates/"  ### Home directory
setwd(paste0(HomeDir,"Administrative"))
nsw <- st_read("nsw.shp")
spnsw <- as_Spatial(nsw)

### Mask - only NSW
GenoPheno.levels <- mask(GenoPheno.levels, nsw)
writeRaster(GenoPheno.levels, filename = "GenoPheno.rat.tif",format = "GTiff", overwrite = T )

par(mfrow=c(1,1))
levelplot(GenoPheno.levels, col.regions = my_palette)+
  layer(sp.polygons(spnsw,col = "black"))

plot(GenoPheno.levels, breaks=c(0,1,2,3,4,5,6,7,8,9,10,99), col=my_palette)

### Attach the other raster layers and examine combinations, and specifically, those 99 values
vegxlu <- mask(vegxlu, nsw)
writeRaster(vegxlu[[1]], filename = "VegxCleared2017.nsw.tif",format = "GTiff", overwrite = T )
writeRaster(vegxlu[[2]], filename = "clum30.2018.nsw.tif",format = "GTiff", overwrite = T )

### Reload and Stack together
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/VASTlayer")
GenoPheno <- raster("GenoPheno.rat.tif")
VegxCleared <- raster("VegxCleared2017.nsw.tif")
clum30 <- raster("clum30.2018.nsw.tif")

geno.pheno.examine <- stack(GenoPheno, VegxCleared,clum30)
names(geno.pheno.examine) <- c("GenoPheno", "VegxCleared", "clum30")
plot(geno.pheno.examine)
#geno.pheno.examine <- trim(geno.pheno.examine)

### Use unique function
u <- unique(geno.pheno.examine, na.last=TRUE, progress="text")
u <- as.data.frame(u)

unclassified <- filter(.data=u,
                       !is.na(VegxCleared),                                                ## Exclude areas where there are NA for extant native vegetation
                       !is.na(clum30),                                                  ## Exclude pixels with NA for land use
                       GenoPheno == 99) %>%   
  distinct(clum30,.keep_all = TRUE) %>% arrange(VegxCleared, clum30)
write.csv(unclassified, file = "unclassified.csv")


### Repeat on 04/05/2020 for non-native and natural CLUM
vegxlu <- stack(VegxCleared, clum30)
names(vegxlu) <- c( "VegxCleared", "clum30")

ff <- function(x) calc(x, VegxLU_reclass)
beginCluster(6)
GenoPheno <- clusterR(vegxlu, ff, export = list("VegxLU_reclass"),
                      filename = "GenoPheno.04052020.tif", format = "GTiff",
                      na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
gc()

GenoPheno <- raster("GenoPheno.04052020.tif")
plot(GenoPheno)

GenoPheno.levels <- ratify(GenoPheno)

# Add some custom labels for each class. Add two kinds to show that
# rasterVis::levelplot() can select one of them.
levels(GenoPheno.levels)[[1]]$GenoPheno <- c("Remnant Genosoils", #1
                                             "Genosoil II", #2
                                             "Phenosoil Low", #3
                                             "Phenosoil Forestry", #4
                                             "Phenosoil Dry Grazing", #5 
                                             "Phenosoil Dry Cropping", #6
                                             "Phenosoil Irrigated Grazing", #7
                                             "Phenosoil Irrigated Cropping", #8
                                             "Phenosoil Degraded", #9
                                             "Phenosoil Rehabilitation")#, #10
                                            # "Unclassified") #99
levels(GenoPheno.levels)

# Custom palette
my_palette <- c("mediumseagreen","olivedrab2","orchid1", 
                "blueviolet","darkkhaki", "gold1",
                "darkseagreen","darkcyan",
                "brown1","antiquewhite2")

writeRaster(GenoPheno.levels, filename = "GenoPheno.04052020.rat.tif",format = "GTiff", overwrite = T )

par(mfrow=c(1,1))
levelplot(GenoPheno.levels, col.regions = my_palette)+
  layer(sp.polygons(spnsw,col = "black"))

#plot(GenoPheno.levels, breaks=c(0,1,2,3,4,5,6,7,8,9,10,99), col=my_palette)

########################################################################################################################################################


# ### Simplified PEDOPHENON classes ----------------------------------------

### Reload and Stack together
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/VASTlayer")
VegxCleared <- raster("VegxCleared2017.nsw.tif")
clum30 <- raster("clum30.2018.nsw.tif")
vegxlu <- stack(VegxCleared, clum30m)
vegxlu <- mask(vegxlu, nsw)
plot(vegxlu)

### Create a layer that separates potentially native vegetation from non-native (and water)
VegxLU_reclass.simple <- function(s,...){
  
  ### First, let's create an empty vector and fill the values with 99 to signal anomalies
  GenPhen <- rep(99, length(s[1]))
  ### VegxCleared Classes
  ### Potentially native remnant - 1
  ### Native cleared - 2
  ### Non-native and forestry plantations - 3
  ### water bodies - 4
  
  ### If there is a missing value for Native vegetation extent or Land Use assign NA
  ### I also want to exclude from the analysis those areas that are water bodies
  GenPhen[is.na(s[1])|s[1]==4|is.na(s[2])] <- NA 
  # GenPhen[is.na(s[1])|is.na(s[2])] <- NA 
  
  ### If the vegetation is potentially native and remnant, use the first classification I made
  Genosoil.Remnant <- c(110:117,120:125,130,131,133,
                        610,611,614, # Lake, lake conservation, lake saline
                        630,631, # River, river conservation
                        650,651,654, # Marsh conservation and saline
                        660,661)
  Genosoil.II <- c(132, # Stock route
                   200:222, # Near natural production
                   314,414, # Environmental forest plantation
                   612,632,652,662) # Lake production, river production, marsh production, estuary production
  Phenosoil.Forestry <- c(310:313,410:413)
  Phenosoil.Grazing <- c(320:325, 420:424)
  Phenosoil.Cropping <- c(330:353,360,362:365,542,430:460, 462:465,361,461) # including rural residential with agriculture
 
  #Phenosoil.Irrigated.Cropping <- c()
  #Phenosoil.Degraded <- c()
  #Phenosoil.Recovery <- 134
  Exclude <- c(500:541,543:595,613,620:623,633,640:643,653,663, 134) # Intensive uses + lake intensive, Dam/reservoir, river intensive, channel/aqueducts
  
  #Phenosoil.Residential <- c(540,541,543:545,553)
  # Phenosoil.Intensive <- c(510,511,512,513,515,520,521,522,523,524,525,526,527,528,530,531,532,533,534,535,537,
  #                          538,550,551,552,554,555,560,561,562,563,564,565,566,567,570,571,572,573,574,575,580,
  #                          581,582,583,584,590,591,592,593,595,
  #                          613,620:623,633,640:643,653) # lake intensive, Dam/reservoir, river intensive, channel/aqueducts
  
  ### Assign classes for Genosoils on native remnant
  GenPhen[s[1]==1 & s[2] %in% Genosoil.Remnant] <- 1 # Remnant Genosoils
  GenPhen[s[1]==1 & s[2] %in% Genosoil.II] <- 2 # Genosoil II, low modification
  
  ### Assign classes for Phenosoil Low in native cleared
  GenPhen[s[1]==2 & s[2] %in% Genosoil.Remnant] <- 3 # Phenosoil low modification
  GenPhen[s[1]==2 & s[2] %in% Genosoil.II] <- 3 # Phenosoil low modification
  
  ### Assign classes for Phenosoil Low in non-native
  GenPhen[s[1]==3 & s[2] %in% Genosoil.Remnant] <- 3 # Phenosoil low modification
  GenPhen[s[1]==3 & s[2] %in% Genosoil.II] <- 3 # Phenosoil low modification
  
  ### Assign Phenosoil classes
  GenPhen[s[1] ==3  & s[2] == 314] <- 4 # Phenosoil.Forestry
  GenPhen[s[1] ==3  & s[2] == 414] <- 4 # Phenosoil.Forestry
  GenPhen[s[2] %in% Phenosoil.Forestry] <- 4 # Phenosoil.Forestry
  GenPhen[s[2] %in% Phenosoil.Grazing]  <- 5 # Phenosoil.Dry.Grazing
  GenPhen[s[2] %in% Phenosoil.Cropping] <- 6 # Phenosoil.Dry.Cropping
  #GenPhen[s[2] %in% Phenosoil.Irrigated.Grazing]  <- 7 # Phenosoil.Irrigated.Grazing
  #GenPhen[s[2] %in% Phenosoil.Irrigated.Cropping] <- 8 # Phenosoil.Irrigated.Cropping
  #GenPhen[s[2] %in% Phenosoil.Degraded] <- 9 # Phenosoil.Degraded
  #GenPhen[s[2] %in% Phenosoil.Recovery] <- 10 # Phenosoil.Recovery
  GenPhen[s[2] %in% Exclude] <- NA # Phenosoil.Intensive - Excluded
  ## Return this value
  GenPhen
}

ff <- function(x) calc(x, VegxLU_reclass.simple)
beginCluster(6)
GenoPheno <- clusterR(vegxlu, ff, export = list("VegxLU_reclass.simple"),
                      filename = "GenoPheno.simple.tif", format = "GTiff",
                      na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
gc()
GenoPheno <- raster("GenoPheno.simple.tif")

GenoPheno.levels <- ratify(GenoPheno)
# Add some custom labels for each class. Add two kinds to show that
# rasterVis::levelplot() can select one of them.
levels(GenoPheno.levels)[[1]]$GenoPheno <- c("Remnant Genosoils", #1
                                             "Genosoil II", #2
                                             "Phenosoil Low", #3
                                             "Phenosoil Forestry", #4
                                             "Phenosoil Grazing", #5 
                                             "Phenosoil Cropping") #6
levels(GenoPheno.levels)
# Custom palette
my_palette <- c("mediumseagreen","olivedrab2","orchid1", 
                "blueviolet","darkkhaki", "gold1")
par(mfrow=c(1,1))
levelplot(GenoPheno.levels, col.regions = my_palette)+
  layer(sp.polygons(spnsw,col = "black"))
writeRaster(GenoPheno.levels, filename = "GenoPheno.simple.rat.tif",format = "GTiff", overwrite = T )


### end of script