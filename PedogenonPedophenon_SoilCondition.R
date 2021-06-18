####################################################################################################################################
### Assess changes in soil condition for different pedogenon/pedophenon subclasses
### Author: Mercedes Roman
### Date: 31/08/2020

##Load packages
#library(Rtools)
library(rlang)
library(ClusterR)
library(rgdal)
library(gdalUtils)
library(raster)
library(sp)
library(sf)
library(dplyr)
library(tidyverse)
library(ggmap)
library(ggplot2) 
library(viridis) # color palettes
library(scales)
library(rasterVis)
library(lattice)
library(gridExtra)
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(mapview) # for interactive maps
library(shiny)   # for web applications
library(foreach)
library(doParallel)
library(geosphere)
library(dendsort)
library(gplots)
library(dendextend)
library(colorspace)
library(elsa)
library(vegan)
library(multcomp)
library(multcompView)
library(nlme)

### Load functions for examination of study areas
source("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/R_Code/Kmeans/FunctionsMapsStudyAreas.R")

# ### 1. Load maps and models to examine -------------------------------------

### Load maps and models to examine
### Load NSW
nsw <- read_sf("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Farms/nsw.shp")
st_crs(nsw)

#### Load the Pedogenon map created withwith 18 continuous covariates and vegetation PCs (25 environmental covariates)
OutDir <-"C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/3.CovSubset18contVegPCs"
setwd(OutDir)
K1000.18ContVeg <- raster("K1000.tif")
K1000.18ContVeg.3857 <-raster("K1000.18ContVeg.3857.tif")
#crs(K1000.18ContVeg.3857) <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#### Load the kmeans model
OutDir <-"C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/2.CovSubsetVeg/"
setwd(OutDir)
load(paste0(OutDir,"kmeans_clorpt.Veg.k1000.RData"))
kmeans_clorpt.Veg.k1000 <- kmeans_clorpt
rm(kmeans_clorpt)

###  Load SCaRP data
load("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/C_sequestration/Output/SCaRP_MIR/SCaRP.data.31082020.RData")

### Plot the locations of SCaRP samples on a map
SCaRP.sp <- SpatialPointsDataFrame(data = SCaRP.Soil.Harmonized.0_30[!is.na(SCaRP.Soil.Harmonized.0_30$Lat), ],
                                   coords = SCaRP.Soil.Harmonized.0_30[!is.na(SCaRP.Soil.Harmonized.0_30$Lat),3:4 ], 
                                   proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

### Subset only those inside NSW
nsw <- as(nsw,Class = "Spatial")
nsw <- spTransform(nsw, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

library(rgeos)
nsw_intersects <-  gIntersects(nsw, SCaRP.sp, byid=TRUE) 
# what kind of object is this?
class(nsw_intersects)
# subset
SCaRP.nsw.sp <- SCaRP.sp[as.vector(nsw_intersects),] ### 4422 locations
# plot
plot (nsw, lwd = 2)
points (SCaRP.sp, col="#aaaaaa")
points (SCaRP.nsw.sp,  col="red") 
### Transform into a dataframe
SCaRP.nsw.df <- as.data.frame(SCaRP.nsw.sp)
rm(nsw_intersects,new.packages,list.of.packages, SCaRP.Soil.MIR )

#####################################################################################################################

setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII")

# ### 2. NSW Genosoil maps for 18 continuous variables + Vegetation P --------

# ### Choose optimal number of branches based on Objective measure --------------------------

viz.map.legend.hclust.Obj <- function(kmodel,Min.nc,Max.nc, indexC = c("silhouette","dunn")) {
  
  ### Extract centroids from model
  centroids <- kmodel$centroids
  ### Extract the index of the centroids that are na/nan/Inf
  Kcent <- as.data.frame(centroids)
  Kcent.nan <- which(apply(Kcent, MARGIN = 1, FUN = function(x) {any(is.na(x))}))
  ### Exclude these clusters from everywhere
  if(length(Kcent.nan) >0) {
    Kcent.exist <- Kcent[-Kcent.nan,]
  } else if(length(Kcent.nan) == 0) {
    Kcent.exist <- Kcent
  }
  # Kcent.exist <- Kcent[-Kcent.nan,]
  ### Hierarchical clustering
  #hc <- hclust(dist(Kcent.exist), method="ward.D2")
  ### Determine the optimal number of clsuters with the Dunn or Shilhhouete index
  hc.op <-NbClust(data = Kcent.exist, distance = "euclidean",
                  min.nc = Min.nc, max.nc = Max.nc, method = "ward.D2", index=indexC)
  # plot(hc.op, main="Hierarchical clustering of kmeans centroids", sub="", xlab="")
  return(hc.op)
}


# ### 3. Figure 2 for paper ------------------------------------------------------------
library(NbClust)
set.seed(2011)
hc.k1000.Veg.op.S<- viz.map.legend.hclust.Obj(kmodel = kmeans_clorpt.Veg.k1000, Max.nc = 30, Min.nc = 2, indexC ="silhouette" );hc.k1000.Veg.op.S$Best.nc
plot(x=2:30, y=hc.k1000.Veg.op.S$All.index)

set.seed(2011)
hc.k1000.Veg.op.D <- viz.map.legend.hclust.Obj(kmodel = kmeans_clorpt.Veg.k1000, Max.nc = 30,Min.nc = 2, indexC ="dunn" )
plot(hc.k1000.Veg.op.D)
plot(x=2:30, y=hc.k1000.Veg.op.D$All.index)
str(hc.k1000.Veg.op.D)
hc.k1000.Veg.op.D$Best.nc

branches.choice <- data.frame(cbind(cluster=c(2:30, 2:30), 
                                    IndexValue =(c(as.numeric(hc.k1000.Veg.op.S$All.index), as.numeric(hc.k1000.Veg.op.D$All.index))),
                                    Index = c(rep("Silhouette", 29), rep("Dunn", 29))))
branches.choice$IndexValue <- as.numeric(branches.choice$IndexValue)
branches.choice$cluster <- as.numeric(branches.choice$cluster)
branches.choice$Index <- as.factor(branches.choice$Index)

ggplot()+
  geom_point(data=branches.choice, aes(x=cluster, y=IndexValue, colour=Index),size=2)+
  xlab("Number of clusters")+
  ylab("Index")+ theme(legend.text=element_text(size=16),
                       legend.title = element_text(size=16),axis.text =element_text(size=14),
                       axis.title =element_text(size=14))+
  theme(legend.position=c(0.85,0.85))
rm(branches.choice)
  
#### Let's map the k=1000, including pre-1750s vegetation
### check hierarchical histogram
hc.k1000.Veg <- viz.map.legend.hclust(kmodel = kmeans_clorpt.Veg.k1000)
plot(hc.k1000.Veg)
### Choose number of branches
plot.branches.k1000Veg <- viz.branches(hc.k1000.Veg, 20)

library("colorspace")
hcl_palettes(plot = TRUE)

### choose palettes
my_palette <- c("Turku","PurpOr","TealGrn","OrYel","Burg","RdPu","Greens",
                 "Peach","GnBu","Lajolla",
                "OrRd", "BurgYl",  "Blues","Heat 2","Dark Mint",
                "SunsetDark", "PuBuGn", "Viridis", "Heat","YlOrRd", "Terrain")
# my_palette <- c("PurpOr","OrYel","TealGrn","BurgYl","RdPu",
#                 "GnBu","YlOrRd","Peach","Turku","Lajolla",
#                 "OrRd", "Greens", "Burg", "Heat 2", "Dark Mint",
#                 "Blues", "SunsetDark", "PuBuGn", "Viridis", "Heat","Terrain")
# my_palette17 <- c("Plasma","PurpOr","RdPu","Burg","TealGrn","Turku","PuBu","OrRd",
#                   "Greens","GnBu","Lajolla","Peach",
#                   "OrYel", "Heat 2","BurgYl", "Dark Mint",
#                   "Blues", "SunsetDark", "PuBuGn", "Viridis", "Heat", "Terrain")
# my_palette17 <- c("Plasma","TealGrn","Greens","Turku","OrRd",
#                   "GnBu","Lajolla","Peach",
#                   "OrYel", "Heat 2","BurgYl", "Dark Mint",
#                   "Blues", "SunsetDark", "PuBuGn", "Viridis", "Heat", "Terrain")
k1000Veg.out <- viz.map.legend.pal(kmodel =  kmeans_clorpt.Veg.k1000,
                                   branchN = 20,
                                   pal.names = my_palette,
                                   legend.name =  "k1000Veg",
                                   kmap = K1000.18ContVeg.3857,
                                   need.proj = FALSE) 
### Plot the legend
plot(k1000Veg.out$legend.plot)
k1000Veg.out$map.out %>%
  addScaleBar( position = "bottomright",
               options = scaleBarOptions(maxWidth = 200, metric = TRUE, imperial = FALSE,
                                         updateWhenIdle = TRUE))
k1000Veg.out$branch.centroids.ord


# ### 4. Genosoils and Phenosoils - simplified ----------------------------------------

#### Number of Remnant Genosoils and Phenosoil classes per CLORPTon?
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/VASTlayer")
GenoPheno <- raster("GenoPheno.simple_r90m.tif")
plot(GenoPheno)
### Create masks
# mask.Geno <- calc(GenoPheno, fun= function(x){ifelse((x==1|x==2),1,NA) })
# mask.GenoRemnant <- calc(GenoPheno, fun= function(x){ifelse(x==1,1,NA) })
# writeRaster(mask.Geno, 
#             filename = "C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/mask.Geno.tif")
# writeRaster(mask.GenoRemnant, 
#             filename = "C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/mask.GenoRemnant.tif")
mask.Geno <-raster("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/mask.Geno.tif")
mask.GenoRemnant <-raster("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/mask.GenoRemnant.tif")

#### stack GenoPheno and GENON maps
CLORPTon.GPh.s <- stack(K1000.18ContVeg ,GenoPheno)
names(CLORPTon.GPh.s) <- c("K1000.18ContVeg", "GenoPheno")
areaR <-raster::area(K1000.18ContVeg, na.rm=TRUE)
CLORPTon.GPh.s <- stack(CLORPTon.GPh.s,areaR)
names(CLORPTon.GPh.s) <- c( "K1000.18ContVeg", "GenoPheno", "Area")
plot(CLORPTon.GPh.s)

### summarize area per combination CLORPTon and GenoPheno
CLORPTon.GPH.df <- getValues(CLORPTon.GPh.s)
CLORPTon.GPH.df <- as.data.frame(CLORPTon.GPH.df)
colnames(CLORPTon.GPH.df) <- c("K1000.18ContVeg", "GenoPheno", "Area")

setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII")

### Summary statistics by GenoPheno, across pedogenons
summary.GenoPheno <-  CLORPTon.GPH.df %>% 
  filter(., !is.na(K1000.18ContVeg), !is.na(GenoPheno)) %>%
  group_by(., as.factor(GenoPheno), .drop=TRUE) %>% ### Group by final GenoPheno
  summarise(., AreaClass = sum(Area, na.rm=TRUE))
NSW_AreaT <- sum(summary.GenoPheno$AreaClass, na.rm=TRUE)
summary.GenoPheno <- summary.GenoPheno %>%
  mutate(., Area_GenoPheno_Perc = AreaClass/NSW_AreaT*100)
summary.GenoPheno <- summary.GenoPheno %>%
  arrange(., desc(Area_GenoPheno_Perc))

CLORPTon.GPH.summary.18contVeg <- CLORPTon.GPH.df %>% 
  filter(., !is.na(K1000.18ContVeg), !is.na(GenoPheno)) %>%
  group_by(., as.factor(K1000.18ContVeg), as.factor(GenoPheno), .drop=TRUE) %>% ### Group by final cluster and GenoPheno
  summarise(., AreaClass = sum(Area, na.rm=TRUE), 
            N = n())
CLORPTon.GPH.summary.18contVeg <- as.data.frame(CLORPTon.GPH.summary.18contVeg)
colnames(CLORPTon.GPH.summary.18contVeg) <- c("K1000.18ContVeg", "GenoPheno", "Area", "N pixels")
CLORPTon.GPH.summary.18contVeg <- CLORPTon.GPH.summary.18contVeg %>% 
  distinct(K1000.18ContVeg,GenoPheno,Area, .keep_all = TRUE) %>% arrange(K1000.18ContVeg, GenoPheno)
CLORPTon.GPH.summary.18contVeg$K1000.18ContVeg <- as.numeric(as.character(CLORPTon.GPH.summary.18contVeg$K1000.18ContVeg))
CLORPTon.GPH.summary.18contVeg$GenoPheno <- as.numeric(as.character(CLORPTon.GPH.summary.18contVeg$GenoPheno))
#write.csv(CLORPTon.GPH.summary.18contVeg, file = "CLORPTon.Cont18Veg.GenoPhenoSimple.Area.csv")

CLORPTon.GPH.summary.18contVeg %>% filter(., GenoPheno == 1) %>% dim(.)
length(kmeans_clorpt.Veg.k1000$obs_per_cluster[kmeans_clorpt.Veg.k1000$obs_per_cluster != 0])

summary(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno == 1,]$Area)
round(summary(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno %in% c(1),]$Area), digits=2)
hist(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno == 1,]$Area, 
     breaks=100, xlab= "Area",
     main= "Area of remnant genosoil by CLORPTon class")
hist(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno == 2,]$Area, breaks=100,
     xlab= "Area",
     main= "Area of genosoil II by CLORPTon class")
summary(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno == 2,]$Area)
hist(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno %in% c(1,2),]$Area, breaks=100)
round(summary(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno %in% c(2),]$Area), digits=2)
round(summary(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno %in% c(1,2),]$Area), digits=2)

#K1000.18ContVeg.Geno <- mask(K1000.18ContVeg, mask.Geno, overwrite=TRUE,
#                             filename="C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/K1000.18ContVeg.Geno.tif")
K1000.18ContVeg.Geno <- raster("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/K1000.18ContVeg.Geno.tif")

# K1000.18ContVeg.GenoRemnant <- mask(K1000.18ContVeg, mask.GenoRemnant,overwrite=TRUE,
#                                     filename="C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/K1000.18ContVeg.GenoRemnant.tif")
K1000.18ContVeg.GenoRemnant <-raster("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/K1000.18ContVeg.GenoRemnant.tif")

length(unique(getValues(K1000.18ContVeg.Geno))) ### 998 without NA
length(unique(getValues(K1000.18ContVeg.GenoRemnant))) ### 998 without NA

### choose palettes
my_palette <- c("Turku","PurpOr","TealGrn","OrYel","Burg","RdPu","Greens",
                "Peach","GnBu","Lajolla",
                "OrRd", "BurgYl",  "Blues","Heat 2","Dark Mint",
                "SunsetDark", "PuBuGn", "Viridis", "Heat","YlOrRd", "Terrain")

# K1000.18ContVeg.Geno.3857 <- projectRaster(K1000.18ContVeg.Geno, crs=CRS("+init=epsg:3857"), method = "ngb",overwrite=TRUE,
#                                         filename="C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/K1000.18ContVeg.Geno.3857.tif")
K1000.18ContVeg.Geno.3857 <- raster("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/K1000.18ContVeg.Geno.3857.tif")

Geno.18ContVeg.NSW.out <- viz.map.legend.pal(kmodel =  kmeans_clorpt.Veg.k1000, branchN = 20, need.proj = FALSE,
                                             pal.names = my_palette,
                                             legend.name =  "k1000.18ContVeg", kmap = K1000.18ContVeg.Geno.3857 ) 
Geno.18ContVeg.NSW.out$map.out %>%
  addCircleMarkers(radius = 0.2, color ="black" , opacity = 1,
                   lng = SCaRP.nsw.df$Longitude,
                   lat = SCaRP.nsw.df$Latitude)

# K1000.18ContVeg.GenoRemnant.3857 <- projectRaster(K1000.18ContVeg.GenoRemnant, crs=CRS("+init=epsg:3857"), method = "ngb",overwrite=TRUE,
#                                                filename="C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/K1000.18ContVeg.GenoRemnant.3857.tif")
K1000.18ContVeg.GenoRemnant.3857<- raster("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/K1000.18ContVeg.GenoRemnant.3857.tif")
GenoRemnant.18ContVeg.NSW.out <- viz.map.legend.pal(kmodel =  kmeans_clorpt.Veg.k1000, branchN = 20, need.proj = TRUE,
                                                    pal.names = my_palette,
                                                    legend.name =  "k1000.18ContVeg", kmap = K1000.18ContVeg.GenoRemnant ) 
GenoRemnant.18ContVeg.NSW.out$map.out%>%
  addScaleBar( position = "bottomright",
               options = scaleBarOptions(maxWidth = 300, metric = TRUE, imperial = FALSE,
                                         updateWhenIdle = TRUE))

### Create a couple zoom for a couple of areas
plot(k1000Veg.out$legend.plot)
k1000Veg.out$map.out %>%
  addScaleBar( position = "bottomright",
               options = scaleBarOptions(maxWidth = 300, metric = TRUE, imperial = FALSE,
                                         updateWhenIdle = TRUE))%>%
  #fitBounds(lng1=149, lat1=-29.9,  lng2=152.20 , lat2=-31.058) %>%
  addRectangles(color = "red",
                lng1=149.5, lat1=-29.9,
                lng2=152.20 , lat2=-31.94,
                fillColor = "transparent" ) %>%
  addRectangles(color = "red",
                lng1=141.6, lat1=-32.89,
                lng2=142.94 , lat2=-33.84,
                fillColor = "transparent" )

### Compute summary statistics for the Pedogenon/ Genosoil/Phenosoil classes
CLORPTon.GPH.summary.18contVeg %>% filter(., GenoPheno == 1) %>% summary()

genopheno.summary <- CLORPTon.GPH.summary.18contVeg[,c( "GenoPheno", "Area")] %>% 
  group_by(.,GenoPheno) %>% ### Group by GenoPheno
  summarise_all(., list(~min(.), ~quantile(., probs=0.25),~median(.), ~mean(.), ~quantile(., probs=0.75), ~max(.)))
genopheno.summary <- as.data.frame(genopheno.summary)

### change names
colnames(genopheno.summary) <- c("Genosoil/Phenosoil", "Min Area", "Q25 Area","Median Area", "Mean Area", "Q75 Area", "Max Area" )
genopheno.summary[,3:7] <- apply(genopheno.summary[,3:7], MARGIN = 2, FUN = function(x) round(x,1))
genopheno.summary$`Min Area` <- round(genopheno.summary$`Min Area`,3)
#write.csv(genopheno.summary, file = "CLORPTon.Cont18Veg.genopheno.summary.csv")

### Maybe better as histograms for each phenosoil
ggplot(CLORPTon.GPH.summary.18contVeg, aes(Area, fill = as.factor(GenoPheno)))+
geom_histogram(binwidth = 30)
  
ggplot(CLORPTon.GPH.summary.18contVeg, aes(Area))+
  geom_histogram(binwidth = 20)+
  facet_wrap(~as.factor(GenoPheno)) 
summary(CLORPTon.GPH.summary.18contVeg)

### How many unique combinations do we have?
CLORPTon.GPH.summary.18contVeg$combi <- paste0(CLORPTon.GPH.summary.18contVeg$K1000.18ContVeg,"_",CLORPTon.GPH.summary.18contVeg$GenoPheno)
length(unique(CLORPTon.GPH.summary.18contVeg$combi))
length(unique(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno==2,]$combi))
length(unique(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno==3,]$combi))
length(unique(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno==4,]$combi))
length(unique(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno==5,]$combi))
length(unique(CLORPTon.GPH.summary.18contVeg[CLORPTon.GPH.summary.18contVeg$GenoPheno==6,]$combi))

### what percentage of each pedogenon area represent the genosoils?
pedogenon.area.NSW <- CLORPTon.GPH.df %>% 
  filter(., !is.na(K1000.18ContVeg)) %>%
  group_by(., as.factor(K1000.18ContVeg), .drop=TRUE) %>% ### Group by final cluster and GenoPheno
  summarise(., AreaPedogenon = sum(Area, na.rm=TRUE))
pedogenon.area.NSW <- as.data.frame(pedogenon.area.NSW)
colnames(pedogenon.area.NSW) <- c("K1000.18ContVeg", "AreaPedogenon")
pedogenon.area.NSW$K1000.18ContVeg <- as.numeric(as.character(pedogenon.area.NSW$K1000.18ContVeg))

str(pedogenon.area.NSW)
str(CLORPTon.GPH.summary.18contVeg)
CLORPTon.GPH.summary.18contVeg_AreaTotal <- merge(CLORPTon.GPH.summary.18contVeg, pedogenon.area.NSW, by = "K1000.18ContVeg" )
CLORPTon.GPH.summary.18contVeg_AreaTotal$Area_p <- CLORPTon.GPH.summary.18contVeg_AreaTotal$Area/CLORPTon.GPH.summary.18contVeg_AreaTotal$AreaPedogenon*100

### what do the Genosoils represent from their respective pedogenons?
summary(CLORPTon.GPH.summary.18contVeg_AreaTotal[CLORPTon.GPH.summary.18contVeg_AreaTotal$GenoPheno==1,]$Area_p)
hist(CLORPTon.GPH.summary.18contVeg_AreaTotal[CLORPTon.GPH.summary.18contVeg_AreaTotal$GenoPheno==1,]$Area_p, breaks=50)

### Table with summary statistics
genopheno.perc.summary <- CLORPTon.GPH.summary.18contVeg_AreaTotal[,c( "GenoPheno", "Area_p")] %>% 
  group_by(.,GenoPheno) %>% ### Group by GenoPheno
  summarise_all(., list(~min(.), ~quantile(., probs=0.25),~median(.), ~mean(.), ~quantile(., probs=0.75), ~max(.)))
genopheno.perc.summary <- as.data.frame(genopheno.perc.summary)

### change names
colnames(genopheno.perc.summary) <- c("Genosoil/Phenosoil", "Min %", "Q25 %","Median %", "Mean %", "Q75 %", "Max %" )
genopheno.perc.summary[,3:7] <- apply(genopheno.perc.summary[,3:7], MARGIN = 2, FUN = function(x) round(x,1))
genopheno.perc.summary[,2] <- round(genopheno.perc.summary[,2],3)
#write.csv(genopheno.perc.summary, file = "CLORPTon.Cont18Veg.genopheno.summary.PercAreaPedogenon.csv")

ggplot(CLORPTon.GPH.summary.18contVeg_AreaTotal,
       aes(Area, after_stat(density), colour = as.factor(GenoPheno)))+
  geom_freqpoly(binwidth = 50, lwd=1 )

ggplot(CLORPTon.GPH.summary.18contVeg_AreaTotal,
       aes(Area, colour = as.factor(GenoPheno)))+
  geom_freqpoly(binwidth = 10, lwd=1 )

ggplot(CLORPTon.GPH.summary.18contVeg_AreaTotal,
       aes(Area, after_stat(density), fill = as.factor(GenoPheno)))+
  geom_histogram(binwidth = 30)

ggplot(CLORPTon.GPH.summary.18contVeg_AreaTotal,
       aes(Area, fill = as.factor(GenoPheno)))+
  geom_histogram(binwidth = 20)

ggplot(CLORPTon.GPH.summary.18contVeg_AreaTotal,
       aes(Area_p, fill = as.factor(GenoPheno)))+
  geom_histogram(binwidth =5)

ggplot(CLORPTon.GPH.summary.18contVeg_AreaTotal,
       aes(Area_p, after_stat(density), colour = as.factor(GenoPheno)))+
  geom_freqpoly(binwidth = 5, lwd=1.5 ) +
  labs(colour = "Genosoil / Phenosoil")+
  scale_colour_manual(name="Genosoil / Phenosoil",
                          breaks=c(1,2,3,4,5,6),
                        values= c("mediumseagreen","olivedrab2","orchid1","blueviolet","darkkhaki", "gold1"),
                          labels=c("Remnant Genosoil", "Genosoil II", "Genosoil Cleared",
                                   "Phenosoil Forestry", "Phenosoil Grazing", "Phenosoil Cropping"))+
  xlab("Contribution to the pedogenon of origin (% area)") +
  ylab("Density")
#+  scale_color_viridis_d(labels=c("Remnant Genosoil", "Genosoil II", "Genosoil Cleared",
  #                               "Phenosoil Forestry", "Phenosoil Grazing", "Phenosoil Cropping"))

# ### 5. At what level do I present the pedogenon classes? -------------------
library(wesanderson)
#pal <-wes_palette("Darjeeling1", 5, type = c("discrete"))
binpal <- colorBin(palette = c("mediumseagreen","olivedrab2","orchid1","blueviolet","darkkhaki", "gold1"),
                   bins = c(1,2,3,4,5,6,7),
                   na.color = "transparent")
SCaRP.at.GenoPheno <- leaflet() %>%
  # Base groups
  addTiles(group="OSM (default)") %>%
  addProviderTiles("Esri.WorldImagery", group = "World Imagery") %>% # , group = "World Imagery"
  addRasterImage(GenoPheno, opacity = 1, colors=binpal, project=TRUE, 
                 maxBytes = 300000000, group = "Pedogenons and Pedophenons") %>%
  fitBounds(lng1=140, lat1=-38, lng2=154, lat2=-28) %>%
  leafem::addMouseCoordinates() %>%
  addLayersControl(
    baseGroups = c("OSM (default)","World Imagery"),
    overlayGroups = c("Pedogenons and Pedophenons"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addCircleMarkers(radius = 0.2, color ="black" , opacity = 1,
                   lng = SCaRP.nsw.df$Longitude,
                   lat = SCaRP.nsw.df$Latitude)
SCaRP.at.GenoPheno

rm(pal,areaR, hc.k1000.Veg.op.D, hc.k1000.Veg.op.S)

### Save image
save.image("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/20052021.RData")
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII")

### Merge with SCaRP data
SCaRP.sp_clorpton <- raster::extract(CLORPTon.GPh.s, SCaRP.nsw.sp, sp=TRUE)
SCaRP.clorpton.df <-as.data.frame(SCaRP.sp_clorpton)
SCaRP.clorpton.df <- SCaRP.clorpton.df[,1:57]
str(SCaRP.clorpton.df)
SCaRP.clorpton.df$combi <- paste0(SCaRP.clorpton.df$K1000.18ContVeg,"_",SCaRP.clorpton.df$GenoPheno)

### eliminate observations without GenoPheno or CLORPTon defined
SCaRP.clorpton.df <- SCaRP.clorpton.df[!is.na(SCaRP.clorpton.df$K1000.18ContVeg) & !is.na(SCaRP.clorpton.df$GenoPheno),]
### Also, with estimates of the soil properties
SCaRP.clorpton.df <- SCaRP.clorpton.df[complete.cases(SCaRP.clorpton.df[,c("BDMean","OC","MIR_POC.st","MIR_HOC.st","MIR_ROC.st","MIR_TN",
                                                                           "MIR_pH","MIR_clay","MIR_Si","MIR_Al","MIR_Fe","POC_density",
                                                                           "ROC_density","HOC_density","SOC_density")]),]
### Correct soil data
### MIR_POC.st, MIR_ROC.st, MIR_HOC.st have been already standardized to sum up to the total OC concentration
### Similarly, the densities, "POC_density","ROC_density","HOC_density","SOC_density" are also corrected.
colnames(SCaRP.clorpton.df)
summary(SCaRP.clorpton.df[,c("BDMean","OC","MIR_POC.st","MIR_HOC.st","MIR_ROC.st","MIR_TN",
                             "MIR_pH","MIR_clay","MIR_Si","MIR_Al","MIR_Fe","POC_density",
                             "ROC_density","HOC_density","SOC_density")])
### Set negative values to 0
SCaRP.clorpton.df$MIR_Fe <- ifelse(SCaRP.clorpton.df$MIR_Fe < 0, 0, SCaRP.clorpton.df$MIR_Fe)
### Set upper limit for clay and silica
SCaRP.clorpton.df$MIR_clay <- ifelse(SCaRP.clorpton.df$MIR_clay > 1000, 1000, SCaRP.clorpton.df$MIR_clay)
SCaRP.clorpton.df$MIR_Si <- ifelse(SCaRP.clorpton.df$MIR_Si > 1000, 1000, SCaRP.clorpton.df$MIR_Si)
dim(SCaRP.clorpton.df)

### Count how many observations are per pedogenon for 0-10 cm. Summary statistics.
### only for topsoils 0-10 cm
length(unique(SCaRP.clorpton.df[SCaRP.clorpton.df$UDepth==0,]$K1000.18ContVeg))
k1000.Pedogenon.N.topsoil <- SCaRP.clorpton.df %>% 
  filter(., !is.na(K1000.18ContVeg), !is.na(GenoPheno)) %>%
  filter(., UDepth ==0) %>%
  ### Group by Pedogenon
  group_by(., as.factor(K1000.18ContVeg), .drop=TRUE) %>%
  summarise(., N = n())
summary(k1000.Pedogenon.N.topsoil$N)
sd(k1000.Pedogenon.N.topsoil$N)

### Calculate summary statistics for each soil property, count the number of observations in each GenoPheno category
SCaRP.clorpton.df.summary <- SCaRP.clorpton.df %>% 
  filter(., UDepth ==0) %>%
  filter(., !is.na(K1000.18ContVeg), !is.na(GenoPheno)) %>%
  ### Group by GenoPheno
  group_by(., as.factor(GenoPheno), .drop=TRUE) %>%
  summarise(.,  N = n(),
            mean.BD = mean(BDMean,na.rm=TRUE),
            mean.OC = mean(OC,na.rm=TRUE),
            mean.POC = mean(MIR_POC.st,na.rm=TRUE),
            mean.ROC = mean(MIR_ROC.st,na.rm=TRUE),
            mean.HOC = mean(MIR_HOC.st,na.rm=TRUE),
            mean.TN = mean(MIR_TN,na.rm=TRUE),
            mean.pH = mean(MIR_pH,na.rm=TRUE),
            mean.clay = mean(MIR_clay,na.rm=TRUE),
            mean.Si = mean(MIR_Si,na.rm=TRUE),
            mean.Fe = mean(MIR_Fe,na.rm=TRUE),
            mean.Al = mean(MIR_Al,na.rm=TRUE),
            sd.BD = sd(BDMean,na.rm=TRUE),
            sd.OC = sd(OC,na.rm=TRUE),
            sd.POC = sd(MIR_POC.st,na.rm=TRUE),
            sd.ROC = sd(MIR_ROC.st,na.rm=TRUE),
            sd.HOC = sd(MIR_HOC.st,na.rm=TRUE),
            sd.TN = sd(MIR_TN,na.rm=TRUE),
            sd.pH = sd(MIR_pH,na.rm=TRUE),
            sd.clay = sd(MIR_clay,na.rm=TRUE),
            sd.Si = sd(MIR_Si,na.rm=TRUE),
            sd.Fe = sd(MIR_Fe,na.rm=TRUE),
            sd.Al = sd(MIR_Al,na.rm=TRUE))
colnames(SCaRP.clorpton.df.summary)

SCaRP.clorpton.df.summary2 <- SCaRP.clorpton.df.summary %>%
  mutate(., BD = paste0(round(mean.BD,2)," \u00b1",round(sd.BD,2))) %>%
  mutate(., OC = paste0(round(mean.OC,2)," \u00b1",round(sd.OC,2))) %>%
  mutate(., POC = paste0(round(mean.POC,2)," \u00b1",round(sd.POC,2))) %>%
  mutate(., HOC = paste0(round(mean.HOC,2)," \u00b1",round(sd.HOC,2))) %>%
  mutate(., TN = paste0(round(mean.TN,2)," \u00b1",round(sd.TN,2))) %>%
  mutate(., pH = paste0(round(mean.pH,2)," \u00b1",round(sd.pH,2))) %>%
  mutate(., clay = paste0(round(mean.clay,2)," \u00b1",round(sd.clay,2))) %>%
  mutate(., Si = paste0(round(mean.Si,2)," \u00b1",round(sd.Si,2))) %>%
  mutate(., Fe = paste0(round(mean.Fe,2)," \u00b1",round(sd.Fe,2))) %>%
  mutate(., Al = paste0(round(mean.Al,2)," \u00b1",round(sd.Al,2))) 
  
SCaRP.clorpton.df.summary2 <- as.data.frame(SCaRP.clorpton.df.summary2)
SCaRP.clorpton.df.summary2 <- SCaRP.clorpton.df.summary2[,c(1,2,25:34)]
#write.csv(SCaRP.clorpton.df.summary2, file = "SCaRP.genopheno.df.summary2.csv")

### What is the relationship between precipitation and GenoPheno type?
summary(SCaRP.clorpton.df$MIR_pH)
SCaRP.clorpton.df %>%
  group_by(., as.factor(GenoPheno)) %>%
  summarize(., mean.MAP= round(mean(X30ARain, na.rm=TRUE),1),
            sd.MAP= round(sd(X30ARain, na.rm=TRUE),1),
            mean.MAT= round(mean(X30ATemp, na.rm=TRUE),1),
            sd.MAT= round(sd(X30ATemp, na.rm=TRUE),1)) %>%
  mutate(., MAT = paste0(round(mean.MAT,1)," \u00b1",round(sd.MAT,1))) %>%
  mutate(., MAP = paste0(round(mean.MAP,0)," \u00b1",round(sd.MAP,0)))%>% as.data.frame() %>% write.csv(., "SCaRP.GenoPhenoMAT.MAP.csv")

ggplot()+geom_boxplot(aes(x=as.factor(GenoPheno), y=X30ARain),data=SCaRP.clorpton.df)
ggplot()+geom_boxplot(aes(x=as.factor(GenoPheno), y=MIR_pH),data=SCaRP.clorpton.df)
ggplot()+geom_point(aes(colour=as.factor(GenoPheno), shape = as.factor(GenoPheno),
                        x=X30ARain, y=MIR_pH, alpha=1/10), data=SCaRP.clorpton.df)

### Select only topsoil observations
SCaRP.clorpton.df.topsoil <- SCaRP.clorpton.df %>% 
  filter(., !is.na(K1000.18ContVeg), !is.na(GenoPheno)) %>% 
  filter(., LDepth <= 10)
dim(SCaRP.clorpton.df.topsoil)

### Now, count number of observations per subclass
SCaRP.clorpton.df.summary <- SCaRP.clorpton.df %>% 
  filter(., !is.na(K1000.18ContVeg), !is.na(GenoPheno)) %>%
  filter(., UDepth ==0) %>%
  ### Group by final cluster and GenoPheno
  group_by(., as.factor(K1000.18ContVeg), as.factor(GenoPheno), .drop=TRUE) %>%
  summarise(.,  N = n())
SCaRP.clorpton.df.summary <- as.data.frame(SCaRP.clorpton.df.summary)
colnames(SCaRP.clorpton.df.summary) <- c("K1000.18ContVeg", "GenoPheno","N")
SCaRP.clorpton.df.summary <- SCaRP.clorpton.df.summary %>% 
  distinct(K1000.18ContVeg, GenoPheno,  .keep_all = TRUE) %>% arrange(desc(N), K1000.18ContVeg)
SCaRP.clorpton.df.summary <- as.data.frame(SCaRP.clorpton.df.summary)
head(SCaRP.clorpton.df.summary)
summary(SCaRP.clorpton.df.summary)
summary(SCaRP.clorpton.df.summary$N)
hist(SCaRP.clorpton.df.summary$N, breaks=50)
table(SCaRP.clorpton.df.summary$GenoPheno)


# ### 6.1 Preliminary analysis for RDA --------------------------------------------
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII")
library(vegan)
SCaRP.clorpton.df.topsoil$K1000.18ContVeg <- as.factor(SCaRP.clorpton.df.topsoil$K1000.18ContVeg)
SCaRP.clorpton.df.topsoil$GenoPheno <- as.factor(SCaRP.clorpton.df.topsoil$GenoPheno)

Response <- SCaRP.clorpton.df.topsoil[,c("BDMean","OC","MIR_POC.st" ,"MIR_HOC.st","MIR_ROC.st","MIR_TN","MIR_pH","MIR_clay","MIR_Si","MIR_Al","MIR_Fe")]
Response <- SCaRP.clorpton.df.topsoil[,c("BDMean","MIR_POC.st","MIR_pH","MIR_clay","MIR_Si","MIR_Al")]
Response <- as.data.frame(Response)

### Check correlations and transformations
pairs(Response)
library(corrplot)
Response %>% 
  cor(., use = "pairwise.complete.obs") %>% 
  corrplot.mixed(.,upper = "ellipse", lower = "number", number.cex=0.75, tl.cex=0.70,tl.col = "black")
summary(Response)
library(psych)
pairs.panels(Response, smooth = TRUE, hist.col = "darkcyan", stars = TRUE,
             labels= c("BD", "OC", "POC", "HOC", "ROC", "TN", "pH", "Clay", "Si", "Al", "Fe") )

### Eliminate OC if I use fractions
library(faraway)
options(scipen=999)
vif(Response[,c("BDMean","OC","MIR_pH","MIR_clay","MIR_Si","MIR_Al")])
decorana(Response)

### Do I need transformations?
library(AID)
Response$BD_t <-  boxcoxnc(Response$BDMean, method = "pt", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
         alpha = 0.05, verbose = TRUE)$tf.data
Response$POC_t <-boxcoxnc(Response$MIR_POC.st, method = "pt", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
         alpha = 0.05, verbose = TRUE)$tf.data
boxcoxnc(Response$MIR_pH, method = "pt", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
         alpha = 0.05, verbose = TRUE) ### Bi-modal distribution
boxcoxnc(Response$MIR_clay, method = "pt", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
         alpha = 0.05, verbose = TRUE)
boxcoxnc(Response$MIR_clay, method = "pt", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
         alpha = 0.05, verbose = TRUE)
Response$Clay_t <- boxcoxnc(Response$MIR_clay, method = "pt", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
                            alpha = 0.05, verbose = TRUE)$tf.data
Response$Si_t <- boxcoxnc(Response$MIR_Si, method = "pt", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
                            alpha = 0.05, verbose = TRUE)$tf.data
Response$Al_t <- boxcoxnc(Response$MIR_Al, method = "pt", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
                          alpha = 0.05, verbose = TRUE)$tf.data

Response_t <- Response[, c("BD_t","POC_t","MIR_pH","Clay_t","Si_t","Al_t")]
colnames(Response_t) <- c("BD","POC","pH","Clay","Si","Al")
### Better bind
SCaRP.clorpton.df.topsoil<- cbind(SCaRP.clorpton.df.topsoil,Response_t)


# ### 6.2 RDA on soil properties including the cumulative index of terrestrial pressure, and individual GenoPheno data layers ----------

### Let's try on the whole dataset
SCaRP.clorpton.df.topsoil
### Plot the locations of SCaRP samples on a map
SCaRP.clorpton.topsoil.sp <- SpatialPointsDataFrame(data = SCaRP.clorpton.df.topsoil,
                                   coords = SCaRP.clorpton.df.topsoil[,3:4], 
                                   proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

### Load the layers that create the GenoPheno layer
NatVeg5m <- raster("C:/Covariates/Vegetation/nswnativevegetationextentv1p25m2017/NSW_Native_Vegetation_Extent_v1p2_5m_2017.tif")
plot(NatVeg5m)
clum50 <- raster("C:/Covariates/LandUse/geotiff_clum_50m1218m/clum_50m1218m.tif")
plot(clum50)
extveg002 <- raster("C:/Covariates/Vegetation/NSWExtantNativeVegetationV2/extveg002.tif")
plot(extveg002)
### and the Human Modification index
HMc <- raster("C:/Covariates/HMc/gHM/gHM/gHM.tif")
par(mfrow=c(1,1))
plot(HMc)
SCaRP.clorpton.topsoil.sp.p <- spTransform(x = SCaRP.clorpton.topsoil.sp, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
points(SCaRP.clorpton.topsoil.sp.p)
SCaRP.clorpton.topsoil.sp.p <- raster::extract(HMc, SCaRP.clorpton.topsoil.sp.p)
HMi <- SCaRP.clorpton.topsoil.sp.p
### Attach to the old df
SCaRP.clorpton.df.topsoil <- cbind(SCaRP.clorpton.df.topsoil, HMi)

### Repeat for the other layers
SCaRP.clorpton.topsoil.sp.p <- spTransform(x = SCaRP.clorpton.topsoil.sp, 
                                           CRS("+proj=lcc +lat_0=-33.25 +lon_0=147 +lat_1=-30.75 +lat_2=-35.75 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +units=m +no_defs"))
plot(NatVeg5m)
points(SCaRP.clorpton.topsoil.sp.p)
SCaRP.clorpton.topsoil.sp.p <- raster::extract(NatVeg5m, SCaRP.clorpton.topsoil.sp.p)
NatVeg <- SCaRP.clorpton.topsoil.sp.p; str(NatVeg)
### Attach to the old df
SCaRP.clorpton.df.topsoil <- cbind(SCaRP.clorpton.df.topsoil, NatVeg)

SCaRP.clorpton.topsoil.sp.p <- spTransform(x = SCaRP.clorpton.topsoil.sp, 
                                           CRS("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))
plot(clum50)
points(SCaRP.clorpton.topsoil.sp.p)
SCaRP.clorpton.topsoil.sp.p <- raster::extract(clum50, SCaRP.clorpton.topsoil.sp.p)
clum <- SCaRP.clorpton.topsoil.sp.p; str(clum)
### Attach to the old df
SCaRP.clorpton.df.topsoil <- cbind(SCaRP.clorpton.df.topsoil, clum)

SCaRP.clorpton.topsoil.sp.p <- spTransform(x = SCaRP.clorpton.topsoil.sp, 
                                           CRS("+proj=lcc +lat_0=-33.25 +lon_0=147 +lat_1=-30.75 +lat_2=-35.75 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +units=m +no_defs"))
plot(extveg002)
points(SCaRP.clorpton.topsoil.sp.p)
SCaRP.clorpton.topsoil.sp.p <- raster::extract(extveg002, SCaRP.clorpton.topsoil.sp.p)
ExtVeg <- SCaRP.clorpton.topsoil.sp.p; str(ExtVeg)
### Attach to the old df
SCaRP.clorpton.df.topsoil <- cbind(SCaRP.clorpton.df.topsoil, ExtVeg)

### RDA analysis. All factors
SCaRP.clorpton.df.topsoil$K1000.18ContVeg <- factor(SCaRP.clorpton.df.topsoil$K1000.18ContVeg)
SCaRP.clorpton.df.topsoil$GenoPheno <- factor(SCaRP.clorpton.df.topsoil$GenoPheno)
SCaRP.clorpton.df.topsoil$ManageClass2 <- factor(SCaRP.clorpton.df.topsoil$ManageClass2)
SCaRP.clorpton.df.topsoil$NatVeg <- as.factor(SCaRP.clorpton.df.topsoil$NatVeg)
SCaRP.clorpton.df.topsoil$clum <- as.factor(SCaRP.clorpton.df.topsoil$clum)
SCaRP.clorpton.df.topsoil$ExtVeg <- as.factor(SCaRP.clorpton.df.topsoil$ExtVeg)
Response <- SCaRP.clorpton.df.topsoil[,c("BD","POC","pH","Clay","Si","Al")]

save.image("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/25052021.RData")

# ### 6.3 RDA on those suclasses with at least 5 observations -----------------

combi.5N <- SCaRP.clorpton.df.topsoil %>% 
  filter(., !is.na(K1000.18ContVeg), !is.na(GenoPheno)) %>%
  ### Group by combi
  group_by(., as.factor(combi), .drop=TRUE) %>%
  summarise(.,  N = n())

### selct those combis with at least 5N
combi.5N <- combi.5N[combi.5N$N >= 5,]
combi.select <- as.character(combi.5N$`as.factor(combi)`)
SCaRP.clorpton.df.topsoil.combi5N <- SCaRP.clorpton.df.topsoil[SCaRP.clorpton.df.topsoil$combi %in% combi.select,]

hist(combi.5N$N, breaks = 30)
SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg <- factor(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)
SCaRP.clorpton.df.topsoil.combi5N$GenoPheno <- factor(SCaRP.clorpton.df.topsoil.combi5N$GenoPheno)
SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 <- factor(SCaRP.clorpton.df.topsoil.combi5N$ManageClass2)
SCaRP.clorpton.df.topsoil.combi5N$NatVeg <- factor(SCaRP.clorpton.df.topsoil.combi5N$NatVeg)
SCaRP.clorpton.df.topsoil.combi5N$clum <- factor(SCaRP.clorpton.df.topsoil.combi5N$clum)
SCaRP.clorpton.df.topsoil.combi5N$ExtVeg <- factor(SCaRP.clorpton.df.topsoil.combi5N$ExtVeg)

### simplify
SCaRP.clorpton.df.topsoil.combi5N$ManageSimple <- NA
pastures <- c("Pasture", "Permanent Pasture", "pasture", "improved pasture",
              "unimproved pasture", "Low input perennial Pasture",
              "Pasture-continuous grazing (native-introduced)",
              "Native/naturalised pasture", "Native pasture", "Introduced pasture",
              "Rotational Grazing","Low grazed system- TSR or road side reserve",
              "Low input permanent pasture (native/naturalised)",
              "High input perennial Pasture",
              "Permanent pasture - System 17a comparison",
              "Pasture-rotational grazing (native-introduced)",
              "Perennial tropical pasture", "Pasture (native)",
              "Pasture (introduced)","Pasture (irrigated perennial)",
              "Pasture (irrigated annual)", "Pasture (dryland)",
              "Pasture (lucerne)", "Perennial pasture",
              "Pasture - dairy", "Pasture - sheep beef",
              "Annual pasture + Beef", "Annual pasture",
              "Annual pasture + Dairy",
              "Perennial pasture + Dairy",
              "Remanent vegetation",
              "Grazing Cattle (High)",
              "Mxd Gzng 20% shp & 80%Catle",
              "Mxd Grazing Ctle & Horses",   
              "Mxd Grazg 50 sheep 50 cattle",
              "Grzg Cattle Low","Beef & Dairy High Prodn",
              "Grazing Cattle (Low)")
setdiff(pastures, unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="Pasture",]$ManageClass))
setdiff(unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="Pasture",]$ManageClass),pastures)

cropping <- c("Cropping", "Intensive vegetable production","Horticulture","Irrigated Hort",
              "Continuous cropping", "rotational cropping","mixed cropping",
              "No-till/min-till cereal dominant cropping",
              "Irrigated cotton in rotation","Crop (wheat)",
              "Crop (barley)", "Mixed cropping", "cropping",
              "Fallow (mechanical)", "Fallow (chemical)", "Fallow (other)")
setdiff(cropping,unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="Cropping",]$ManageClass))
setdiff(unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="Cropping",]$ManageClass),cropping)
setdiff(unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="Horticulture",]$ManageClass),cropping)

#fallow <- c("Fallow (mechanical)", "Fallow (chemical)", "Fallow (other)")
pasture_cropping <- c("Crop-Pasture Rotation", 
                      "Crop-pasture rotation",
                      "Mixed",
                      "Mxd Grazing & Cropg",
                      "Organic ammedments on cropping & pasture",
                      "Emerging farming systems: Carbon farming",
                      "Emerging farming systems: Crop-pasture rotation (pasture)",
                      "Emerging farming systems: Crop-pasture rotation (crop)",
                      "Emerging farming systems: Pasture-cropping",
                      "Conventional farming- fenceline comparision w 10a",
                      "intermittent")
setdiff(pasture_cropping,unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="CarbonFarming",]$ManageClass))
setdiff(unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="CarbonFarming",]$ManageClass),pasture_cropping)
setdiff(pasture_cropping,unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="Mixed",]$ManageClass))
setdiff(unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="Mixed",]$ManageClass),pasture_cropping)
setdiff(unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="Unknown",]$ManageClass),pasture_cropping)
setdiff(unique(SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 =="Horticulture",]$ManageClass),pasture_cropping)

SCaRP.clorpton.df.topsoil.combi5N$ManageSimple <- 
  ifelse(is.na(SCaRP.clorpton.df.topsoil.combi5N$ManageClass), NA,
         ifelse(SCaRP.clorpton.df.topsoil.combi5N$ManageClass %in% cropping, "Cropping",
                ifelse(SCaRP.clorpton.df.topsoil.combi5N$ManageClass %in% pastures, "Pasture",
                       #ifelse(SCaRP.clorpton.df.topsoil.combi5N$ManageClass %in% fallow, "Fallow",
                       ifelse(SCaRP.clorpton.df.topsoil.combi5N$ManageClass %in% pasture_cropping, "Pasture-cropping", -9999))))
table(SCaRP.clorpton.df.topsoil.combi5N$ManageSimple)
SCaRP.clorpton.df.topsoil.combi5N$ManageSimple <- as.factor(SCaRP.clorpton.df.topsoil.combi5N$ManageSimple)

### ### ### ### ### ### ### ### ### ### ###
rm(Response, Response_t)
Response_6 <- SCaRP.clorpton.df.topsoil.combi5N[,c("BD","POC","pH","Clay","Si","Al")]
Response_3 <- SCaRP.clorpton.df.topsoil.combi5N[,c("BD","POC","pH")]

### Suitability of the rule based algorithm
SCaRP.5N.rda.3layers <- rda(formula = Response_3 ~ NatVeg+ExtVeg+clum,
                    data= SCaRP.clorpton.df.topsoil.combi5N,
                    scale=TRUE)
RsquareAdj(SCaRP.5N.rda.3layers)
SCaRP.5N.rda.GenoPheno <- rda(formula = Response_3 ~ GenoPheno,
                            data= SCaRP.clorpton.df.topsoil.combi5N,
                            scale=TRUE)
RsquareAdj(SCaRP.5N.rda.GenoPheno)

RsquareAdj(SCaRP.5N.rda.3layers)$adj.r.squared -
RsquareAdj(SCaRP.5N.rda.GenoPheno)$adj.r.squared
SCaRP.5N.rda.3layers <- rda(formula = Response_6 ~ NatVeg+ExtVeg+clum,
                            data= SCaRP.clorpton.df.topsoil.combi5N,
                            scale=TRUE)
SCaRP.5N.rda.GenoPheno <- rda(formula = Response_6 ~ GenoPheno,
                              data= SCaRP.clorpton.df.topsoil.combi5N,
                              scale=TRUE)
RsquareAdj(SCaRP.5N.rda.3layers)
RsquareAdj(SCaRP.5N.rda.GenoPheno)
RsquareAdj(SCaRP.5N.rda.3layers)$adj.r.squared -
  RsquareAdj(SCaRP.5N.rda.GenoPheno)$adj.r.squared
### We lose ~ 4% variance explained in the reclassification step.

### Now with my variables K1000.18ContVeg & GenoPheno
SCaRP.5N.rda <- rda(formula = Response_6 ~ K1000.18ContVeg+GenoPheno,
                    data= SCaRP.clorpton.df.topsoil.combi5N,
                    scale=TRUE)
RsquareAdj(SCaRP.5N.rda)

SCaRP.5N.rda.1 <- rda(formula = Response_6 ~ K1000.18ContVeg+GenoPheno ,
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)
SCaRP.5N.rda.2 <- rda(formula = Response_6 ~ K1000.18ContVeg ,
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)
SCaRP.5N.rda.3 <- rda(formula = Response_6 ~ GenoPheno ,
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)
SCaRP.5N.rda.4 <- rda(formula = Response_6 ~ K1000.18ContVeg + Condition(GenoPheno),
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)
SCaRP.5N.rda.5 <- rda(formula = Response_6 ~ GenoPheno + Condition(K1000.18ContVeg),
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)
rm(SCaRP.5N.rda.1,SCaRP.5N.rda.2,SCaRP.5N.rda.3,SCaRP.5N.rda.4,SCaRP.5N.rda.5)

##3 with varpart
mod.varpart <- varpart(Y = Response_6, ~K1000.18ContVeg,~GenoPheno,
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)
plot(mod.varpart)

summary(SCaRP.5N.rda)
par(mfrow=c(1,1))
plot(SCaRP.5N.rda, scaling=1)

## Test statistic significance
anova.cca(SCaRP.5N.rda)
anova.cca(SCaRP.5N.rda, by="axis")
par(mfrow=c(1,1))
plot(SCaRP.5N.rda, scaling=2, type="text")


# Figure 7 from Revised article -------------------------------------------


### Create nice plots
plot(SCaRP.5N.rda, scaling=1, main="Triplot RDA scaling 1 - wa scores",type="text")
var.sc <- scores(SCaRP.5N.rda, choices=1:2, scaling=1, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
ordilabel(SCaRP.5N.rda, dis="species", cex=1, font=2, fill="white", col="red", scaling=1)

plot(SCaRP.5N.rda, scaling=2, main="Triplot RDA scaling 2")
var.sc <- scores(SCaRP.5N.rda, choices=1:2, scaling=2, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="darkcyan")
ordilabel(SCaRP.5N.rda, dis="species", cex=1, font=2, scaling=2, fill="white", col="darkcyan")

# Create groups
pch.group <- c(22,21,1,2)
colvec  <- viridis_pal(option = "C")(length(levels(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)))
colvec  <- viridis_pal(option = "C")(length(levels(SCaRP.clorpton.df.topsoil.combi5N$GenoPheno)))

#colvec  <- viridis_pal(option = "D")(length(levels(SCaRP.clorpton.df.topsoil.combi5N.sub$GenoPheno)))

### Getting the color from the same color palette as the pedogenon map
SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg <- factor(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)
pedogenons.scarp <- as.numeric(as.character(levels(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)))
pal.rda.5N <- k1000Veg.out$branch.centroids.ord[k1000Veg.out$branch.centroids.ord$Centroid %in% pedogenons.scarp,]$colors
colvec  <- pal.rda.5N

with(SCaRP.clorpton.df.topsoil.combi5N, levels(K1000.18ContVeg))
with(SCaRP.clorpton.df.topsoil.combi5N, levels(GenoPheno))
scl <- 1 ## scaling = 1
plot(SCaRP.5N.rda, type = "n", scaling = scl, 
     xlim=c(-0.6,0.6), 
     ylim=c(-0.6,0.6),
     # xlim=c(-1,1), 
     # ylim=c(-1,1),
     xlab="RDA Axis 1",
     ylab="RDA Axis 2")
# plot(SCaRP.5N.rda, type = "n", scaling = scl,
#      xlim=c(-1,1),
#      xlab="RDA Axis 1",
#      ylab="RDA Axis 2")
with(SCaRP.clorpton.df.topsoil.combi5N, 
     points(SCaRP.5N.rda, display = "sites", 
            col = colvec[K1000.18ContVeg],
            scaling = scl,
            pch = pch.group[GenoPheno],
            bg = colvec[K1000.18ContVeg]))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      points(SCaRP.5N.rda, display = "sites", 
#             col = colvec[GenoPheno],
#             scaling = scl,
#             pch = pch.group[GenoPheno],
#             bg = colvec[GenoPheno]))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("topleft",
            legend = c("Remnant Genosoil", #"Genosoil II", 
                       "Phenosoil Cleared", "Phenosoil Grazing",
                       "Phenosoil Cropping"),
            bty = "n",
            col="black",
            pch = pch.group,
            pt.bg  = "black"))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      legend("topleft",
#             legend = levels(GenoPheno),
#             bty = "n",
#             col=colvec,
#             pch = pch.group,
#             pt.bg  = colvec))


#text(SCaRP.5N.rda, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
var.sc <- scores(SCaRP.5N.rda, choices=1:2, scaling=scl, display="sp")
SCaRP.5N.rda.sc <- scores(SCaRP.5N.rda, choices=1:2, scaling=scl)$species
#arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
arrows(0,0,0.12*SCaRP.5N.rda.sc[,1], 0.12*SCaRP.5N.rda.sc[,2], length=0, lty=1, col="blue")
ordilabel(0.12*SCaRP.5N.rda.sc, dis="species", cex=1, font=2, fill="white", col="blue")
# arrows(0,0,SCaRP.5N.rda.sc[,1], SCaRP.5N.rda.sc[,2], length=0, lty=1, col="blue")
# ordilabel(SCaRP.5N.rda.sc, dis="species", cex=1, font=2, fill="white", col="blue")

### Map only those pedogenons present
K1000.18ContVeg.3857.scarptop5N <- calc(K1000.18ContVeg.3857, fun= function(x) {ifelse(x %in% pedogenons.scarp, x, NA)})

k1000Veg.cont.SCaRP.Topsoil5N <- viz.map.legend.pal(kmodel =  kmeans_clorpt.Veg.k1000,
                                                    branchN = 20,
                                                    pal.names = my_palette,
                                                    legend.name =  "k1000VegSCaRP_Topsoil",
                                                    kmap = K1000.18ContVeg.3857.scarptop5N,
                                                    need.proj = FALSE) 
binpal <- colorBin(palette = k1000Veg.cont.SCaRP.Topsoil5N$branch.centroids.ord$colors,
                   bins = c(k1000Veg.cont.SCaRP.Topsoil5N$branch.centroids.ord$Centroid,
                            tail(k1000Veg.cont.SCaRP.Topsoil5N$branch.centroids.ord$Centroid,1)+1),
                   na.color = "transparent")

leaflet() %>%
  # Base groups
  addTiles(group="OSM (default)") %>%
  #addProviderTiles("Esri.WorldImagery", group = "World Imagery") %>% # , group = "World Imagery"
  addProviderTiles("Jawg.Matrix", group = "Jawg.Matrix") %>%
  #addProviderTiles("Stamen.Toner", group = "Stamen.Toner") %>%
  #addProviderTiles("Thunderforest.MobileAtlas", group = "Thunderforest.MobileAtlas")%>%
  addRasterImage(K1000.18ContVeg.3857.scarptop5N, opacity = 1, colors=binpal, project=FALSE, 
                 maxBytes = 300000000, group = "Pedogenons") %>%
  #fitBounds(lng1=140, lat1=-38, lng2=154, lat2=-28) %>%
   leafem::addMouseCoordinates() %>%
  addLayersControl(
  baseGroups = c("OSM (default)","Jawg.Matrix"), #"Stamen.TerrainBackground", "Thunderforest.MobileAtlas"),
  overlayGroups = c("Pedogenons"),
  options = layersControlOptions(collapsed = FALSE))%>%
  addScaleBar( position = "bottomright",
               options = scaleBarOptions(maxWidth = 200, metric = TRUE, imperial = FALSE,
                                         updateWhenIdle = TRUE)) %>%
  addCircleMarkers(radius = 0.1, color ="white", opacity = 1,
                   lng = SCaRP.clorpton.df.topsoil.combi5N$Longitude,
                   lat = SCaRP.clorpton.df.topsoil.combi5N$Latitude )

k1000Veg.cont.SCaRP.Topsoil5N$map.out %>%
  # addProviderTiles("Jawg.Light", group = "Jawg.Light") %>%
  # addLayersControl(
  #   baseGroups = c("OSM (default)","World Imagery", "Topo Map", "Jawg.Light"),
  #   overlayGroups = c("Pedogenons"),
  #   options = layersControlOptions(collapsed = FALSE)
  # ) %>%
  addCircleMarkers(radius = 0.1, color ="black", opacity = 1,
                   lng = SCaRP.clorpton.df.topsoil.combi5N$Longitude,
                   lat = SCaRP.clorpton.df.topsoil.combi5N$Latitude )%>%
  addScaleBar( position = "bottomright",
               options = scaleBarOptions(maxWidth = 200, metric = TRUE, imperial = FALSE,
                                         updateWhenIdle = TRUE))
k1000Veg.cont.SCaRP.Topsoil5N$map.out %>%
  # addProviderTiles("Thunderforest.MobileAtlas", group = "Thunderforest.MobileAtlas") %>%
  # addLayersControl(
  #   baseGroups = c("OSM (default)","World Imagery", "Topo Map", "Thunderforest.MobileAtlas"),
  #   overlayGroups = c("Pedogenons"),
  #   options = layersControlOptions(collapsed = FALSE)
  # ) %>%
  addCircleMarkers(radius = 0.1, color ="white", opacity = 0.8,
                   lng = SCaRP.clorpton.df.topsoil.combi5N$Longitude,
                   lat = SCaRP.clorpton.df.topsoil.combi5N$Latitude)%>%
  addScaleBar( position = "bottomright",
               options = scaleBarOptions(maxWidth = 200, metric = TRUE, imperial = FALSE,
                                         updateWhenIdle = TRUE))

######################################################################################

### how about soil type?
SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass <- as.factor(SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass)
SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 <- as.factor(SCaRP.clorpton.df.topsoil.combi5N$ManageClass2)

SCaRP.5N.rda.SoilType <- rda(formula = Response_6 ~ Aust_SoilClass+GenoPheno,
                    data= SCaRP.clorpton.df.topsoil.combi5N,
                    scale=TRUE)
RsquareAdj(SCaRP.5N.rda.SoilType)

mod.varpart <- varpart(Y = Response_6, ~Aust_SoilClass,~GenoPheno,
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)
plot(mod.varpart)
# SCaRP.5N.rda.SoilType <- rda(formula = Response_6 ~ Aust_SoilClass+HMi,
#                              data= SCaRP.clorpton.df.topsoil.combi5N,
#                              scale=TRUE)


## Test statistic significance
anova.cca(SCaRP.5N.rda.SoilType)
anova.cca(SCaRP.5N.rda.SoilType, by="axis")
par(mfrow=c(1,1))
plot(SCaRP.5N.rda.SoilType, scaling=2, type="text")


# Figure 8 from revised paper ---------------------------------------------

### Create nice plots
plot(SCaRP.5N.rda.SoilType, scaling=1, main="Triplot RDA scaling 1 - wa scores",type="text")
var.sc <- scores(SCaRP.5N.rda.SoilType, choices=1:2, scaling=1, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
ordilabel(SCaRP.5N.rda.SoilType, dis="species", cex=1, font=2, fill="white", col="red", scaling=1)

plot(SCaRP.5N.rda.SoilType, scaling=2, main="Triplot RDA scaling 2")
var.sc <- scores(SCaRP.5N.rda.SoilType, choices=1:2, scaling=2, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="darkcyan")
ordilabel(SCaRP.5N.rda.SoilType, dis="species", cex=1, font=2, scaling=2, fill="white", col="darkcyan")

# Create groups
pch.group <- c(22,2,3,21)
pch.group <- c(22,21,1,2)
colvec  <- viridis_pal(option = "B", direction = -1)(length(levels(SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass)))
#colvec  <- brewer.pal(length(levels(SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass)), "Paired")

with(SCaRP.clorpton.df.topsoil.combi5N, levels(Aust_SoilClass))
with(SCaRP.clorpton.df.topsoil.combi5N, levels(GenoPheno))
scl <- 1 ## scaling = 1
plot(SCaRP.5N.rda.SoilType, type = "n", scaling = scl, 
     xlim=c(-0.6,0.6), 
     ylim=c(-0.6,0.8),
     # xlim=c(-1,1), 
     # ylim=c(-1,1),
     xlab="RDA Axis 1",
     ylab="RDA Axis 2")
# plot(SCaRP.5N.rda, type = "n", scaling = scl,
#      xlim=c(-1,1),
#      xlab="RDA Axis 1",
#      ylab="RDA Axis 2")
with(SCaRP.clorpton.df.topsoil.combi5N, 
     points(SCaRP.5N.rda.SoilType, display = "sites", 
            col = colvec[Aust_SoilClass],
            scaling = scl,
            pch = pch.group[GenoPheno],
            bg = colvec[Aust_SoilClass]))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      points(SCaRP.5N.rda.SoilType, display = "sites", 
#             col = colvec[Aust_SoilClass],
#             scaling = scl,
#             pch = 19,
#             bg = colvec[Aust_SoilClass]))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      points(SCaRP.5N.rda, display = "sites", 
#             col = colvec[GenoPheno],
#             scaling = scl,
#             pch = pch.group[GenoPheno],
#             bg = colvec[GenoPheno]))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("topleft",
            legend = c("Remnant Genosoil", #"Genosoil II", 
                       "Phenosoil Cleared", "Phenosoil Grazing",
                       "Phenosoil Cropping"),
            bty = "n",
            col="black",
            pch = pch.group,
            pt.bg  = "black"))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("bottomleft",
            legend = levels(Aust_SoilClass),
            bty = "n",
            col=colvec,
            pch = 19,
            pt.bg  = colvec))

#text(SCaRP.5N.rda, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
var.sc <- scores(SCaRP.5N.rda.SoilType, choices=1:2, scaling=scl, display="sp")
SCaRP.5N.rda.SoilType.sc <- scores(SCaRP.5N.rda.SoilType, choices=1:2, scaling=scl)$species
#arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
arrows(0,0,0.12*SCaRP.5N.rda.SoilType.sc[,1], 0.12*SCaRP.5N.rda.SoilType.sc[,2], length=0, lty=1, col="blue")
ordilabel(0.12*SCaRP.5N.rda.SoilType.sc, dis="species", cex=1, font=2, fill="white", col="blue")
# arrows(0,0,SCaRP.5N.rda.sc[,1], SCaRP.5N.rda.sc[,2], length=0, lty=1, col="blue")
# ordilabel(SCaRP.5N.rda.sc, dis="species", cex=1, font=2, fill="white", col="blue")
# 
# library(RColorBrewer)
# par(mar=c(3,4,2,2))
# display.brewer.all()



### Variance partitioning
SCaRP.5N.rda.SoilType1 <- rda(formula = Response_6 ~ Aust_SoilClass+GenoPheno ,
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)
SCaRP.5N.rda.SoilType2 <- rda(formula = Response_6 ~ Aust_SoilClass ,
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)
SCaRP.5N.rda.SoilType3 <- rda(formula = Response_6 ~ GenoPheno ,
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)
SCaRP.5N.rda.SoilType4 <- rda(formula = Response_6 ~ Aust_SoilClass + Condition(GenoPheno),
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)
SCaRP.5N.rda.SoilType5 <- rda(formula = Response_6 ~ GenoPheno + Condition(Aust_SoilClass),
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)


### what about the branches???
str(k1000Veg.out$branch.centroids.ord)
k1000Veg.out$branch.centroids.ord$Centroid <- as.factor(k1000Veg.out$branch.centroids.ord$Centroid)
SCaRP.clorpton.df.topsoil.combi5N <- merge(SCaRP.clorpton.df.topsoil.combi5N, k1000Veg.out$branch.centroids.ord,
                                           by.x ="K1000.18ContVeg", by.y = "Centroid"  )
# SCaRP.clorpton.df.topsoil.combi5N <- SCaRP.clorpton.df.topsoil.combi5N[,1:70]
# colnames(SCaRP.clorpton.df.topsoil.combi5N)[[69]] <- "Branch"
# colnames(SCaRP.clorpton.df.topsoil.combi5N)[[70]] <- "colors"
SCaRP.clorpton.df.topsoil.combi5N$Branch <- as.factor(SCaRP.clorpton.df.topsoil.combi5N$Branch)

SCaRP.5N.rda.Branch <- rda(formula = Response_6 ~ Branch + GenoPheno,
                             data= SCaRP.clorpton.df.topsoil.combi5N,
                             scale=TRUE)

RsquareAdj(SCaRP.5N.rda.Branch)

mod.varpart <- varpart(Y = Response_3, ~Branch,~GenoPheno,
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)
plot(mod.varpart)

## Test statistic significance
anova.cca(SCaRP.5N.rda.Branch)
anova.cca(SCaRP.5N.rda.Branch, by="axis")
par(mfrow=c(1,1))
plot(SCaRP.5N.rda.Branch, scaling=2, type="text")

### Create nice plots
plot(SCaRP.5N.rda.Branch, scaling=1, main="Triplot RDA scaling 1 - wa scores",type="text")
var.sc <- scores(SCaRP.5N.rda.Branch, choices=1:2, scaling=1, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
ordilabel(SCaRP.5N.rda.Branch, dis="species", cex=1, font=2, fill="white", col="red", scaling=1)

plot(SCaRP.5N.rda.Branch, scaling=2, main="Triplot RDA scaling 2")
var.sc <- scores(SCaRP.5N.rda.Branch, choices=1:2, scaling=2, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="darkcyan")
ordilabel(SCaRP.5N.rda.Branch, dis="species", cex=1, font=2, scaling=2, fill="white", col="darkcyan")

# Create groups
pch.group <- c(22,2,3,21)
colvec  <- viridis_pal(option = "B", direction = 1)(length(levels(SCaRP.clorpton.df.topsoil.combi5N$Branch)))
# colvec[3] <- "red"
# colvec[5] <- "deeppink"

with(SCaRP.clorpton.df.topsoil.combi5N, levels(Branch))
with(SCaRP.clorpton.df.topsoil.combi5N, levels(GenoPheno))
scl <- 1 ## scaling = 1
plot(SCaRP.5N.rda.Branch, type = "n", scaling = scl, 
     xlim=c(-0.6,0.6), 
     ylim=c(-0.8,0.6),
     # xlim=c(-1,1), 
     # ylim=c(-1,1),
     xlab="RDA Axis 1",
     ylab="RDA Axis 2")
# plot(SCaRP.5N.rda, type = "n", scaling = scl,
#      xlim=c(-1,1),
#      xlab="RDA Axis 1",
#      ylab="RDA Axis 2")
with(SCaRP.clorpton.df.topsoil.combi5N, 
     points(SCaRP.5N.rda.Branch, display = "sites", 
            col = colvec[Branch],
            scaling = scl,
            pch = pch.group[GenoPheno],
            bg = colvec[Branch]))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      points(SCaRP.5N.rda, display = "sites", 
#             col = colvec[GenoPheno],
#             scaling = scl,
#             pch = pch.group[GenoPheno],
#             bg = colvec[GenoPheno]))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("topleft",
            legend = c("Remnant Genosoil", #"Genosoil II", 
                       "Phenosoil Cleared", "Phenosoil Grazing",
                       "Phenosoil Cropping"),
            bty = "n",
            col="black",
            pch = pch.group,
            pt.bg  = "black"))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("bottomleft",
            legend = levels(Branch),
            bty = "n",
            col=colvec,
            pch = 19,
            pt.bg  = colvec))

#text(SCaRP.5N.rda, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
var.sc <- scores(SCaRP.5N.rda.Branch, choices=1:2, scaling=scl, display="sp")
SCaRP.5N.rda.Branch.sc <- scores(SCaRP.5N.rda.Branch, choices=1:2, scaling=scl)$species
#arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
arrows(0,0,0.12*SCaRP.5N.rda.Branch.sc[,1], 0.12*SCaRP.5N.rda.Branch.sc[,2], length=0, lty=1, col="blue")
ordilabel(0.12*SCaRP.5N.rda.Branch.sc, dis="species", cex=1, font=2, fill="white", col="blue")
# arrows(0,0,SCaRP.5N.rda.sc[,1], SCaRP.5N.rda.sc[,2], length=0, lty=1, col="blue")
# ordilabel(SCaRP.5N.rda.sc, dis="species", cex=1, font=2, fill="white", col="blue")
save.image("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/25052021.RData")
load("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/25052021.RData")
# ### 9. Permanova --------------------------------------------------------

# 9.1 Permanova on the K1000ContVeg + GenoPheno on >= 5N ----------------------
Response_6 <- SCaRP.clorpton.df.topsoil.combi5N[,c("BD","POC","pH","Clay","Si","Al")]

### PERMANOVA
set.seed(1984)
permanova.SCARP.5N.Stable.k1000 <- adonis2(Response_6 ~ K1000.18ContVeg*GenoPheno,
                                           data=SCaRP.clorpton.df.topsoil.combi5N,
                                           method = "mahalanobis", by="terms",
                                           #strata = as.numeric(as.character(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)),
                                           permutations=9999)
permanova.SCARP.5N.Stable.k1000
densityplot(permustats(permanova.SCARP.5N.Stable.k1000))

scarp.dist <- vegdist(x = Response_6, method = "mahalanobis")

### Check for pairwise comparison of centroids.
install.packages('devtools')
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(vegan)
SCaRP.clorpton.df.topsoil.combi5N$ii <- droplevels(interaction(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg, 
                                                               SCaRP.clorpton.df.topsoil.combi5N$GenoPheno))
set.seed(1234)
test.pairwise <- pairwise.adonis(x = scarp.dist, 
                                 factors = SCaRP.clorpton.df.topsoil.combi5N$ii, perm=9999)
scarp.pairwise.df <- as.data.frame(test.pairwise)

scarp.pairwise.df <- separate(data=scarp.pairwise.df, pairs,
                                into=c("Subclass1", "Subclass2"), sep = " vs ",
                                remove = FALSE)
str(scarp.pairwise.df)
#ph.emm.ii.pairs.df2$Subclass1 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass1)
#ph.emm.ii.pairs.df2$Subclass2 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass2)
scarp.pairwise.df <- separate(data=scarp.pairwise.df, col = Subclass1,
                                into=c("Pedogenon1", "GenoPheno1"), 
                                remove = FALSE)
scarp.pairwise.df <- separate(data=scarp.pairwise.df, col = Subclass2,
                                into=c("Pedogenon2", "GenoPheno2"), 
                                remove = FALSE)
scarp.pairwise.df.sub <- scarp.pairwise.df[scarp.pairwise.df$Pedogenon1 ==  scarp.pairwise.df$Pedogenon2,]

list.7.pairwise <- scarp.pairwise.df.sub$Pedogenon1


PdGn <- unique(scarp.pairwise.df.sub$Pedogenon1)

list.adj.p.stable <- list()
for(i in 1:length(unique(scarp.pairwise.df.sub$Pedogenon1))){
  
  SCaRP.clorpton.df.topsoil.combi5N.i <- SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg==PdGn[[i]],]
  scarp.dist.i <- vegdist(x = SCaRP.clorpton.df.topsoil.combi5N.i[,c("BD","POC","pH","Clay","Si","Al")], method = "mahalanobis")
  
  SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno <- factor(SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno)
  set.seed(i+1908)
  test.pairwise.i <- pairwise.adonis(x = scarp.dist.i, 
                                     factors = SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno, perm=999)
  list.adj.p.stable[[i]] <-test.pairwise.i 
  print(test.pairwise.i)
}

names(list.adj.p.stable) <- PdGn


meandist(scarp.dist.i, grouping =SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno )

bd.GenoPheno.i <- betadisper(d =scarp.dist.i, 
                             group=SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno)
set.seed(1984)
plot(bd.GenoPheno.i)



### Get that 274 out
SCaRP.clorpton.df.topsoil.combi5N.274 <- SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg=="274",]
scarp.dist.274 <- vegdist(x = SCaRP.clorpton.df.topsoil.combi5N.274[,c("BD","POC","pH","Clay","Si","Al")], method = "mahalanobis")

SCaRP.clorpton.df.topsoil.combi5N.274$GenoPheno <- factor(SCaRP.clorpton.df.topsoil.combi5N.274$GenoPheno)
test.pairwise.274 <- pairwise.adonis(x = scarp.dist.274, 
                                     factors = SCaRP.clorpton.df.topsoil.combi5N.274$GenoPheno, perm=999)
rm(list.7.pairwise)


set.seed(1984)
bd.Pedogenon <- betadisper(d =scarp.dist,
                           group=SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)
boxplot(bd.Pedogenon)
#TukeyHSD(bd.Pedogenon)
### Null hypothesis of no difference in dispersion between groups
set.seed(1984)
permutest(bd.Pedogenon)
## we reject the null hypothesis, and accept the alternative hypothesis. BUUU
#anova(bd.Pedogenon)
set.seed(1984)
bd.GenoPheno <- betadisper(d =scarp.dist,
                           group=SCaRP.clorpton.df.topsoil.combi5N$GenoPheno)
set.seed(1984)
boxplot(bd.GenoPheno)
#TukeyHSD(bd.GenoPheno)
### Null hypothesis of no difference in dispersion between groups
permutest(bd.GenoPheno)
## we reject the null hypothesis, and accept the alternative hypothesis. BUUU

# set.seed(1984)
# permanova.SCARP.5N <- adonis2(Response_6 ~ K1000.18ContVeg + GenoPheno,
#                               data=SCaRP.clorpton.df.topsoil.combi5N,
#                               method = "mahalanobis", by="terms",
#                               #strata = as.numeric(as.character(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)),
#                               permutations=9999)
# permanova.SCARP.5N
# densityplot(permustats(permanova.SCARP.5N))

### Subset some
SCaRP.clorpton.df.topsoil.combi5N.7 <- SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg %in% list.7.pairwise,]
  
Response6.7 <- SCaRP.clorpton.df.topsoil.combi5N.7[,c("BD","POC","pH")]
Response6.7.dist <- vegdist(x = Response6.7, method = "mahalanobis")
SCaRP.clorpton.df.topsoil.combi5N.7$ii <- factor(SCaRP.clorpton.df.topsoil.combi5N.7$ii)
SCaRP.clorpton.df.topsoil.combi5N.7$K1000.18ContVeg <- factor(SCaRP.clorpton.df.topsoil.combi5N.7$K1000.18ContVeg)

Response6.7.dist.groups <- meandist(Response6.7.dist, grouping =SCaRP.clorpton.df.topsoil.combi5N.7$K1000.18ContVeg )
Response6.7.dist.ii <- meandist(Response6.7.dist, grouping =SCaRP.clorpton.df.topsoil.combi5N.7$ii )
dimnames(Response6.7.dist.ii)

ii <- c("340.1 - 340.5", "719.3 - 719.5", "936.3 - 936.6", "274.5 - 274.6", "462.5 - 462.6", "661.5 - 661.6", "796.5 - 796.6")
comparisonGenoPgeno <- c("1 - 5", "3 - 5", "3 - 6", "5 - 6", "5 - 6", "5 - 6", "5 - 6")
comp.Dist <- data.frame(ii = ii, GenoPheno.c =comparisonGenoPgeno, Dist.pair = rep(NA,7) )
comp.Dist[1,3]<-Response6.7.dist.ii["340.1","340.5"]
comp.Dist[2,3]<-Response6.7.dist.ii["719.3","719.5"]
comp.Dist[3,3]<-Response6.7.dist.ii["936.3","936.6"]
comp.Dist[4,3]<-Response6.7.dist.ii["274.5","274.6"]
comp.Dist[5,3]<-Response6.7.dist.ii["462.5","462.6"]
comp.Dist[6,3]<-Response6.7.dist.ii["661.5","661.6"]
comp.Dist[7,3]<-Response6.7.dist.ii["796.5","796.6"]

boxplot(comp.Dist$Dist.pair ~ comp.Dist$GenoPheno.c)
summary(comp.Dist$Dist.pair)


# 9.2 Permanova on the Soil Type + GenoPheno on >= 5N ----------------------

### PERMANOVA
set.seed(1984)
permanova.SCARP.5N.Stable.SoilType <- adonis2(Response_6 ~ Aust_SoilClass*GenoPheno,
                                       data=SCaRP.clorpton.df.topsoil.combi5N,
                                       method = "mahalanobis", by="terms",
                                       #strata = (as.character(SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass)),
                                       permutations=999)
permanova.SCARP.5N.Stable.SoilType
densityplot(permustats(permanova.SCARP.5N.Stable.SoilType))

scarp.dist <- vegdist(x = Response_6, method = "mahalanobis")
bd.Aust_SoilClass <- betadisper(d =scarp.dist,
                                group=SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass)
boxplot(bd.Aust_SoilClass)
#TukeyHSD(bd.Pedogenon)
### Null hypothesis of no difference in dispersion between groups
set.seed(1984)
permutest(bd.Aust_SoilClass, permutations=999)

#### ### ### 
SCaRP.clorpton.df.topsoil.combi5N$Aus.ii <- droplevels(interaction(SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass, 
                                                               SCaRP.clorpton.df.topsoil.combi5N$GenoPheno))
SCaRP.clorpton.df.topsoil.combi5N$Aus.ii <- factor(SCaRP.clorpton.df.topsoil.combi5N$Aus.ii)

set.seed(1234)
test.pairwise <- pairwise.adonis(x = scarp.dist, 
                                 factors = SCaRP.clorpton.df.topsoil.combi5N$Aus.ii, perm=999)
scarp.pairwise.df <- as.data.frame(test.pairwise)

scarp.pairwise.df <- separate(data=scarp.pairwise.df, pairs,
                              into=c("Subclass1", "Subclass2"), sep = " vs ",
                              remove = FALSE)
str(scarp.pairwise.df)
#ph.emm.ii.pairs.df2$Subclass1 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass1)
#ph.emm.ii.pairs.df2$Subclass2 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass2)
scarp.pairwise.df <- separate(data=scarp.pairwise.df, col = Subclass1,
                              into=c("ASC1", "GenoPheno1"), 
                              remove = FALSE)
scarp.pairwise.df <- separate(data=scarp.pairwise.df, col = Subclass2,
                              into=c("ASC2", "GenoPheno2"), 
                              remove = FALSE)
scarp.pairwise.df.sub <- scarp.pairwise.df[scarp.pairwise.df$ASC1 ==  scarp.pairwise.df$ASC2,]

scarp.pairwise.df.sub$ASC1
nrow(scarp.pairwise.df.sub)

ASC <- unique(scarp.pairwise.df.sub$ASC1)
names(list.adj.p.dynamic) <- ASC
list.adj.p.dynamic <- list()
for(i in 1:length(unique(scarp.pairwise.df.sub$ASC1))){
  
  SCaRP.clorpton.df.topsoil.combi5N.i <- SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass==ASC[[i]],]
  scarp.dist.i <- vegdist(x = SCaRP.clorpton.df.topsoil.combi5N.i[,c("BD","POC","pH","Clay","Si","Al")], method = "mahalanobis")
  
  SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno <- factor(SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno)
  set.seed(i+1908)
  test.pairwise.i <- pairwise.adonis(x = scarp.dist.i, 
                                     factors = SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno, perm=999)
  list.adj.p.dynamic[[i]] <-test.pairwise.i 
  print(test.pairwise.i)
}


# 9.3 Permanova on Branch + GenoPheno on >= 5N ----------------------

### PERMANOVA
set.seed(1984)
permanova.SCARP.5N.Stable.Branch <- adonis2(Response_6 ~ Branch*GenoPheno,
                                       data=SCaRP.clorpton.df.topsoil.combi5N,
                                       method = "mahalanobis", by="terms",
                                       #strata = as.numeric(as.character(SCaRP.clorpton.df.topsoil.combi5N$Branch)),
                                       permutations=9999)
permanova.SCARP.5N.Stable.Branch
densityplot(permustats(permanova.SCARP.5N.Stable.Branch))

scarp.dist <- vegdist(x = Response_6, method = "mahalanobis")
bd.Branch <- betadisper(d =scarp.dist,
                        group=SCaRP.clorpton.df.topsoil.combi5N$Branch)
boxplot(bd.Branch)
#TukeyHSD(bd.Pedogenon)
### Null hypothesis of no difference in dispersion between groups
permutest(bd.Branch)

#### ### ### 
SCaRP.clorpton.df.topsoil.combi5N$Branch.ii <- droplevels(interaction(SCaRP.clorpton.df.topsoil.combi5N$Branch, 
                                                                   SCaRP.clorpton.df.topsoil.combi5N$GenoPheno))
SCaRP.clorpton.df.topsoil.combi5N$Branch.ii <- factor(SCaRP.clorpton.df.topsoil.combi5N$Branch.ii)

set.seed(1234)
test.pairwise <- pairwise.adonis(x = scarp.dist, 
                                 factors = SCaRP.clorpton.df.topsoil.combi5N$Branch.ii, perm=999)
scarp.pairwise.df <- as.data.frame(test.pairwise)

scarp.pairwise.df <- separate(data=scarp.pairwise.df, pairs,
                              into=c("Subclass1", "Subclass2"), sep = " vs ",
                              remove = FALSE)
str(scarp.pairwise.df)
#ph.emm.ii.pairs.df2$Subclass1 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass1)
#ph.emm.ii.pairs.df2$Subclass2 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass2)
scarp.pairwise.df <- separate(data=scarp.pairwise.df, col = Subclass1,
                              into=c("Branch1", "GenoPheno1"), 
                              remove = FALSE)
scarp.pairwise.df <- separate(data=scarp.pairwise.df, col = Subclass2,
                              into=c("Branch2", "GenoPheno2"), 
                              remove = FALSE)
scarp.pairwise.df.sub <- scarp.pairwise.df[scarp.pairwise.df$Branch1 ==  scarp.pairwise.df$Branch2,]

scarp.pairwise.df.sub$Branch1
nrow(scarp.pairwise.df.sub)

Branch <- unique(scarp.pairwise.df.sub$Branch1)
names(list.adj.p.stable) <- Branch
list.adj.p.stable <- list()
for(i in 1:length(Branch)){
  
  SCaRP.clorpton.df.topsoil.combi5N.i <- SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$Branch==Branch[[i]],]
  scarp.dist.i <- vegdist(x = SCaRP.clorpton.df.topsoil.combi5N.i[,c("BD","POC","pH","Clay","Si","Al")], method = "mahalanobis")
  
  SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno <- factor(SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno)
  set.seed(i+1908)
  test.pairwise.i <- pairwise.adonis(x = scarp.dist.i, 
                                     factors = SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno, perm=999)
  list.adj.p.stable[[i]] <-test.pairwise.i 
  print(test.pairwise.i)
}



# ### 10. Dynamic soil properties only. -----------------------------------

# ### 10.1 RDA on dynamic soil properties - k1000 + GenoPheno (>5N) -------
load("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/25052021.RData")
Response_3 <- SCaRP.clorpton.df.topsoil.combi5N[,c("BD","POC","pH")]

SCaRP.5N.rda.Dynamic.k1000 <- rda(formula = Response_3 ~ K1000.18ContVeg+GenoPheno,
                            data= SCaRP.clorpton.df.topsoil.combi5N,
                            scale=TRUE)
RsquareAdj(SCaRP.5N.rda.Dynamic.k1000)

SCaRP.5N.rda.1.dynamic <- rda(formula = Response_3 ~ K1000.18ContVeg+GenoPheno ,
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)
SCaRP.5N.rda.2.dynamic <- rda(formula = Response_3 ~ K1000.18ContVeg ,
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)
SCaRP.5N.rda.3.dynamic <- rda(formula = Response_3 ~ GenoPheno ,
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)
SCaRP.5N.rda.4.dynamic <- rda(formula = Response_3 ~ K1000.18ContVeg + Condition(GenoPheno),
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)
SCaRP.5N.rda.5.dynamic <- rda(formula = Response_3 ~ GenoPheno + Condition(K1000.18ContVeg),
                      data= SCaRP.clorpton.df.topsoil.combi5N,
                      scale=TRUE)
### with varpart
mod.varpart.dynamic <- varpart(Y = Response_3, ~ K1000.18ContVeg,~GenoPheno,
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)

summary(SCaRP.5N.rda.dynamic)
par(mfrow=c(1,1))
plot(SCaRP.5N.rda.dynamic, scaling=1)

## Test statistic significance
anova.cca(SCaRP.5N.rda.dynamic)
anova.cca(SCaRP.5N.rda.dynamic, by="axis")
par(mfrow=c(1,1))
plot(SCaRP.5N.rda.dynamic, scaling=2, type="text")

### Create nice plots
plot(SCaRP.5N.rda.dynamic, scaling=1, main="Triplot RDA scaling 1 - wa scores",type="text")
var.sc <- scores(SCaRP.5N.rda.dynamic, choices=1:2, scaling=1, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
ordilabel(SCaRP.5N.rda.dynamic, dis="species", cex=1, font=2, fill="white", col="red", scaling=1)

plot(SCaRP.5N.rda.dynamic, scaling=2, main="Triplot RDA scaling 2")
var.sc <- scores(SCaRP.5N.rda.dynamic, choices=1:2, scaling=2, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="darkcyan")
ordilabel(SCaRP.5N.rda.dynamic, dis="species", cex=1, font=2, scaling=2, fill="white", col="darkcyan")

# Create groups
pch.group <- c(22,2,3,21)
colvec  <- viridis_pal(option = "C")(length(levels(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)))
colvec  <- viridis_pal(option = "C")(length(levels(SCaRP.clorpton.df.topsoil.combi5N$GenoPheno)))

#colvec  <- viridis_pal(option = "D")(length(levels(SCaRP.clorpton.df.topsoil.combi5N.sub$GenoPheno)))

### Getting the color from the same color palette as the pedogenon map
SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg <- factor(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)
pedogenons.scarp <- as.numeric(as.character(levels(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)))
pal.rda.5N <- k1000Veg.out$branch.centroids.ord[k1000Veg.out$branch.centroids.ord$Centroid %in% pedogenons.scarp,]$colors
colvec  <- pal.rda.5N

with(SCaRP.clorpton.df.topsoil.combi5N, levels(K1000.18ContVeg))
with(SCaRP.clorpton.df.topsoil.combi5N, levels(GenoPheno))
scl <- 1 ## scaling = 1
plot(SCaRP.5N.rda.dynamic, type = "n", scaling = scl, 
     xlim=c(-0.6,1), 
     ylim=c(-0.6,0.6),
     # xlim=c(-1,1), 
     # ylim=c(-1,1),
     xlab="RDA Axis 1",
     ylab="RDA Axis 2")
# plot(SCaRP.5N.rda.dynamic, type = "n", scaling = scl,
#      xlim=c(-1,1),
#      xlab="RDA Axis 1",
#      ylab="RDA Axis 2")
with(SCaRP.clorpton.df.topsoil.combi5N, 
     points(SCaRP.5N.rda.dynamic, display = "sites", 
            col = colvec[K1000.18ContVeg],
            scaling = scl,
            pch = pch.group[GenoPheno],
            bg = colvec[K1000.18ContVeg]))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      points(SCaRP.5N.rda.dynamic, display = "sites", 
#             col = colvec[GenoPheno],
#             scaling = scl,
#             pch = pch.group[GenoPheno],
#             bg = colvec[GenoPheno]))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("topright",
            legend = c("Remnant Genosoil", #"Genosoil II", 
                       "Phenosoil Cleared", "Phenosoil Grazing",
                       "Phenosoil Cropping"),
            bty = "n",
            col="black",
            pch = pch.group,
            pt.bg  = "black"))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      legend("topleft",
#             legend = levels(GenoPheno),
#             bty = "n",
#             col=colvec,
#             pch = pch.group,
#             pt.bg  = colvec))


#text(SCaRP.5N.rda.dynamic, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
var.sc <- scores(SCaRP.5N.rda.dynamic, choices=1:2, scaling=scl, display="sp")
SCaRP.5N.rda.dynamic.sc <- scores(SCaRP.5N.rda.dynamic, choices=1:2, scaling=scl)$species
#arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
arrows(0,0,0.10*SCaRP.5N.rda.dynamic.sc[,1], 0.102*SCaRP.5N.rda.dynamic.sc[,2], length=0, lty=1, col="blue")
ordilabel(0.10*SCaRP.5N.rda.dynamic.sc, dis="species", cex=1, font=2, fill="white", col="blue")
# arrows(0,0,SCaRP.5N.rda.dynamic.sc[,1], SCaRP.5N.rda.dynamic.sc[,2], length=0, lty=1, col="blue")
# ordilabel(SCaRP.5N.rda.dynamic.sc, dis="species", cex=1, font=2, fill="white", col="blue")

### ### ### ### ### ### ### ### ### ### ### ### ### ### #### #### ### ### ### ### ### ### ### ### 


# ### 10.2 RDA on dynamic soil properties - Aus Soil Order+ GenoPheno (>5N) -------
### how about soil type?
SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass <- as.factor(SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass)
SCaRP.clorpton.df.topsoil.combi5N$ManageClass2 <- as.factor(SCaRP.clorpton.df.topsoil.combi5N$ManageClass2)
SCaRP.clorpton.df.topsoil.combi5N$ManageSimple <- as.factor(SCaRP.clorpton.df.topsoil.combi5N$ManageSimple)

SCaRP.5N.rda.SoilType <- rda(formula = Response_3 ~ Aust_SoilClass+GenoPheno,
                             data= SCaRP.clorpton.df.topsoil.combi5N,
                             scale=TRUE)
# SCaRP.5N.rda.SoilType <- rda(formula = Response_3 ~ Aust_SoilClass+ManageSimple,
#                              data= SCaRP.clorpton.df.topsoil.combi5N,
#                              scale=TRUE)

RsquareAdj(SCaRP.5N.rda.SoilType)

mod.varpart <- varpart(Y = Response_3, ~Aust_SoilClass,~GenoPheno,
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)

# SCaRP.5N.rda.SoilType <- rda(formula = Response_3 ~ Aust_SoilClass+HMi,
#                              data= SCaRP.clorpton.df.topsoil.combi5N,
#                              scale=TRUE)


## Test statistic significance
anova.cca(SCaRP.5N.rda.SoilType)
anova.cca(SCaRP.5N.rda.SoilType, by="axis")
par(mfrow=c(1,1))
plot(SCaRP.5N.rda.SoilType, scaling=2, type="text")

### Create nice plots
plot(SCaRP.5N.rda.SoilType, scaling=1, main="Triplot RDA scaling 1 - wa scores",type="text")
var.sc <- scores(SCaRP.5N.rda.SoilType, choices=1:2, scaling=1, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
ordilabel(SCaRP.5N.rda.SoilType, dis="species", cex=1, font=2, fill="white", col="red", scaling=1)

plot(SCaRP.5N.rda.SoilType, scaling=2, main="Triplot RDA scaling 2")
var.sc <- scores(SCaRP.5N.rda.SoilType, choices=1:2, scaling=2, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="darkcyan")
ordilabel(SCaRP.5N.rda.SoilType, dis="species", cex=1, font=2, scaling=2, fill="white", col="darkcyan")

# Create groups
pch.group <- c(22,2,3,21)
colvec  <- viridis_pal(option = "B", direction = -1)(length(levels(SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass)))

with(SCaRP.clorpton.df.topsoil.combi5N, levels(Aust_SoilClass))
with(SCaRP.clorpton.df.topsoil.combi5N, levels(GenoPheno))
scl <- 1 ## scaling = 1
plot(SCaRP.5N.rda.SoilType, type = "n", scaling = scl,
      # xlim=c(-1,1), 
      # ylim=c(-1,1),
     xlim=c(-0.6,0.6),
     ylim=c(-0.6,0.8),
     # xlim=c(-1,1), 
     # ylim=c(-1,1),
     xlab="RDA Axis 1",
     ylab="RDA Axis 2")
# plot(SCaRP.5N.rda, type = "n", scaling = scl,
#      xlim=c(-1,1),
#      xlab="RDA Axis 1",
#      ylab="RDA Axis 2")
with(SCaRP.clorpton.df.topsoil.combi5N, 
     points(SCaRP.5N.rda.SoilType, display = "sites", 
            col = colvec[Aust_SoilClass],
            scaling = scl,
            pch = pch.group[GenoPheno],
            bg = colvec[Aust_SoilClass]))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     points(SCaRP.5N.rda.SoilType, display = "sites", 
            col = colvec[Aust_SoilClass],
            scaling = scl,
            pch = 19,
            bg = colvec[Aust_SoilClass]))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      points(SCaRP.5N.rda, display = "sites", 
#             col = colvec[GenoPheno],
#             scaling = scl,
#             pch = pch.group[GenoPheno],
#             bg = colvec[GenoPheno]))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("topleft",
            legend = c("Remnant Genosoil", #"Genosoil II", 
                       "Phenosoil Cleared", "Phenosoil Grazing",
                       "Phenosoil Cropping"),
            bty = "n",
            col="black",
            pch = pch.group,
            pt.bg  = "black"))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("bottomleft",
            legend = levels(Aust_SoilClass),
            bty = "n",
            col=colvec,
            pch = 19,
            pt.bg  = colvec))

#text(SCaRP.5N.rda, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
var.sc <- scores(SCaRP.5N.rda.SoilType, choices=1:2, scaling=scl, display="sp")
SCaRP.5N.rda.SoilType.sc <- scores(SCaRP.5N.rda.SoilType, choices=1:2, scaling=scl)$species
#arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
arrows(0,0,0.12*SCaRP.5N.rda.SoilType.sc[,1], 0.12*SCaRP.5N.rda.SoilType.sc[,2], length=0, lty=1, col="blue")
ordilabel(0.12*SCaRP.5N.rda.SoilType.sc, dis="species", cex=1, font=2, fill="white", col="blue")
# arrows(0,0,SCaRP.5N.rda.sc[,1], SCaRP.5N.rda.sc[,2], length=0, lty=1, col="blue")
# ordilabel(SCaRP.5N.rda.sc, dis="species", cex=1, font=2, fill="white", col="blue")

### Variance partitioning
SCaRP.5N.rda.SoilType1 <- rda(formula = Response_3 ~ Aust_SoilClass+GenoPheno ,
                              data= SCaRP.clorpton.df.topsoil.combi5N,
                              scale=TRUE)
SCaRP.5N.rda.SoilType2 <- rda(formula = Response_3 ~ Aust_SoilClass ,
                              data= SCaRP.clorpton.df.topsoil.combi5N,
                              scale=TRUE)
SCaRP.5N.rda.SoilType3 <- rda(formula = Response_3 ~ GenoPheno ,
                              data= SCaRP.clorpton.df.topsoil.combi5N,
                              scale=TRUE)
SCaRP.5N.rda.SoilType4 <- rda(formula = Response_3 ~ Aust_SoilClass + Condition(GenoPheno),
                              data= SCaRP.clorpton.df.topsoil.combi5N,
                              scale=TRUE)
SCaRP.5N.rda.SoilType5 <- rda(formula = Response_3 ~ GenoPheno + Condition(Aust_SoilClass),
                              data= SCaRP.clorpton.df.topsoil.combi5N,
                              scale=TRUE)

# ### 10.3 RDA on dynamic soil properties - Pedogenon branch + GenoPheno (>5N) -------
### what about the branches???
str(k1000Veg.out$branch.centroids.ord)
k1000Veg.out$branch.centroids.ord$Centroid <- as.factor(k1000Veg.out$branch.centroids.ord$Centroid)
# SCaRP.clorpton.df.topsoil.combi5N <- merge(SCaRP.clorpton.df.topsoil.combi5N, k1000Veg.out$branch.centroids.ord,
#                                            by.x ="K1000.18ContVeg", by.y = "Centroid"  )
# # SCaRP.clorpton.df.topsoil.combi5N <- SCaRP.clorpton.df.topsoil.combi5N[,1:70]
# colnames(SCaRP.clorpton.df.topsoil.combi5N)[[69]] <- "Branch"
# colnames(SCaRP.clorpton.df.topsoil.combi5N)[[70]] <- "colors"
SCaRP.clorpton.df.topsoil.combi5N$Branch <- factor(SCaRP.clorpton.df.topsoil.combi5N$Branch)

SCaRP.5N.rda.Branch <- rda(formula = Response_3 ~ Branch + GenoPheno,
                           data= SCaRP.clorpton.df.topsoil.combi5N,
                           scale=TRUE)

RsquareAdj(SCaRP.5N.rda.Branch)

mod.varpart <- varpart(Y = Response_6, ~Branch,~GenoPheno,
                       data= SCaRP.clorpton.df.topsoil.combi5N,
                       scale=TRUE)
plot(mod.varpart)

## Test statistic significance
anova.cca(SCaRP.5N.rda.Branch)
anova.cca(SCaRP.5N.rda.Branch, by="axis")
par(mfrow=c(1,1))
plot(SCaRP.5N.rda.Branch, scaling=2, type="text")

### Create nice plots
plot(SCaRP.5N.rda.Branch, scaling=1, main="Triplot RDA scaling 1 - wa scores",type="text", xlim=c(-0.6,0.6), 
     ylim=c(-0.8,0.6))
var.sc <- scores(SCaRP.5N.rda.Branch, choices=1:2, scaling=1, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
ordilabel(SCaRP.5N.rda.Branch, dis="species", cex=1, font=2, fill="white", col="red", scaling=1)

plot(SCaRP.5N.rda.Branch, scaling=2, main="Triplot RDA scaling 2")
var.sc <- scores(SCaRP.5N.rda.Branch, choices=1:2, scaling=2, display="sp")
arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="darkcyan")
ordilabel(SCaRP.5N.rda.Branch, dis="species", cex=1, font=2, scaling=2, fill="white", col="darkcyan")

# Create groups
pch.group <- c(22,2,3,21)
#pch.group <- c(21,22,23,24,25,8,11)
colvec  <- viridis_pal(option = "C", direction = -1)(length(levels(SCaRP.clorpton.df.topsoil.combi5N$Branch)))
#colvec  <- viridis_pal(option = "B", direction = -1)(length(levels(SCaRP.clorpton.df.topsoil.combi5N$GenoPheno)))

with(SCaRP.clorpton.df.topsoil.combi5N, levels(Branch))
with(SCaRP.clorpton.df.topsoil.combi5N, levels(GenoPheno))
scl <- 1 ## scaling = 1
plot(SCaRP.5N.rda.Branch, type = "n", scaling = scl, 
     xlim=c(-0.6,0.6), 
     ylim=c(-0.8,0.6),
     # xlim=c(-1,1), 
     # ylim=c(-1,1),
     xlab="RDA Axis 1",
     ylab="RDA Axis 2")
# plot(SCaRP.5N.rda, type = "n", scaling = scl,
#      xlim=c(-1,1),
#      xlab="RDA Axis 1",
#      ylab="RDA Axis 2")
with(SCaRP.clorpton.df.topsoil.combi5N, 
     points(SCaRP.5N.rda.Branch, display = "sites", 
            col = colvec[Branch],
            scaling = scl,
            pch = pch.group[GenoPheno],
            bg = colvec[Branch]))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      points(SCaRP.5N.rda.Branch, display = "sites", 
#             col = colvec[Branch],
#             scaling = scl,
#             pch = pch.group[Branch],
#             bg = colvec[Branch]))
with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("bottomleft",
            legend = levels(Branch),
            bty = "n",
            col=colvec,
            pch = 19,
            pt.bg  = colvec))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      points(SCaRP.5N.rda, display = "sites", 
#             col = colvec[GenoPheno],
#             scaling = scl,
#             pch = pch.group[GenoPheno],
#             bg = colvec[GenoPheno]))

with(SCaRP.clorpton.df.topsoil.combi5N, 
     legend("topleft",
            legend = c("Remnant Genosoil", #"Genosoil II", 
                       "Phenosoil Cleared", "Phenosoil Grazing",
                       "Phenosoil Cropping"),
            bty = "n",
            col="black",
            pch = pch.group,
            pt.bg  = "black"))

# with(SCaRP.clorpton.df.topsoil.combi5N, 
#      legend("bottomleft",
#             legend = levels(Branch),
#             bty = "n",
#             col=colvec,
#             pch = 19,
#             pt.bg  = colvec))

#text(SCaRP.5N.rda, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
var.sc <- scores(SCaRP.5N.rda.Branch, choices=1:2, scaling=scl, display="sp")
SCaRP.5N.rda.Branch.sc <- scores(SCaRP.5N.rda.Branch, choices=1:2, scaling=scl)$species
#arrows(0,0,var.sc[,1], var.sc[,2], length=0, lty=1, col="red")
arrows(0,0,0.12*SCaRP.5N.rda.Branch.sc[,1], 0.12*SCaRP.5N.rda.Branch.sc[,2], length=0, lty=1, col="blue")
ordilabel(0.12*SCaRP.5N.rda.Branch.sc, dis="species", cex=1, font=2, fill="white", col="blue")
# arrows(0,0,SCaRP.5N.rda.sc[,1], SCaRP.5N.rda.sc[,2], length=0, lty=1, col="blue")
# ordilabel(SCaRP.5N.rda.sc, dis="species", cex=1, font=2, fill="white", col="blue")


# ### 11. PERMANOVA on dynamic soil properties ----------------------------

### PERMANOVA
set.seed(1984)
permanova.SCARP.5N.Dynamic.k1000 <- adonis2(Response_3 ~ K1000.18ContVeg*GenoPheno,
                              data=SCaRP.clorpton.df.topsoil.combi5N,
                              method = "mahalanobis", by="terms",
                              #strata = as.numeric(as.character(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)),
                              permutations=9999)
permanova.SCARP.5N.Dynamic.k1000
densityplot(permustats(permanova.SCARP.5N.Dynamic.k1000))

scarp.dist <- vegdist(x = Response_3, method = "mahalanobis")

SCaRP.clorpton.df.topsoil.combi5N$ii <- droplevels(interaction(SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg, 
                                                               SCaRP.clorpton.df.topsoil.combi5N$GenoPheno))
set.seed(1234)
test.pairwise <- pairwise.adonis(x = scarp.dist, 
                                 factors = SCaRP.clorpton.df.topsoil.combi5N$ii, perm=9999)
scarp.pairwise.df <- as.data.frame(test.pairwise)

scarp.pairwise.df <- separate(data=scarp.pairwise.df, pairs,
                              into=c("Subclass1", "Subclass2"), sep = " vs ",
                              remove = FALSE)
str(scarp.pairwise.df)
#ph.emm.ii.pairs.df2$Subclass1 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass1)
#ph.emm.ii.pairs.df2$Subclass2 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass2)
scarp.pairwise.df <- separate(data=scarp.pairwise.df, col = Subclass1,
                              into=c("Pedogenon1", "GenoPheno1"), 
                              remove = FALSE)
scarp.pairwise.df <- separate(data=scarp.pairwise.df, col = Subclass2,
                              into=c("Pedogenon2", "GenoPheno2"), 
                              remove = FALSE)
scarp.pairwise.df.sub <- scarp.pairwise.df[scarp.pairwise.df$Pedogenon1 ==  scarp.pairwise.df$Pedogenon2,]

list.7.pairwise <- scarp.pairwise.df.sub$Pedogenon1

list.adj.p.dynamic <- list()
for(i in 1:length(list.7.pairwise)){
  
  SCaRP.clorpton.df.topsoil.combi5N.i <- SCaRP.clorpton.df.topsoil.combi5N[SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg==list.7.pairwise[[i]],]
  scarp.dist.i <- vegdist(x = SCaRP.clorpton.df.topsoil.combi5N.i[,c("BD","POC","pH")], method = "mahalanobis")
  
  SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno <- factor(SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno)
  set.seed(i+1908)
  test.pairwise.i <- pairwise.adonis(x = scarp.dist.i, 
                                       factors = SCaRP.clorpton.df.topsoil.combi5N.i$GenoPheno, perm=999)
  list.adj.p.dynamic[[i]] <-test.pairwise.i 
  print(test.pairwise.i)
}



set.seed(1984)
bd.Pedogenon <- betadisper(d =scarp.dist,
                           group=SCaRP.clorpton.df.topsoil.combi5N$K1000.18ContVeg)
boxplot(bd.Pedogenon)
#TukeyHSD(bd.Pedogenon)
### Null hypothesis of no difference in dispersion between groups
set.seed(1984)
permutest(bd.Pedogenon)
## we reject the null hypothesis, and accept the alternative hypothesis. BUUU
#anova(bd.Pedogenon)
set.seed(1984)
bd.GenoPheno <- betadisper(d =scarp.dist,
                           group=SCaRP.clorpton.df.topsoil.combi5N$GenoPheno)
set.seed(1984)
boxplot(bd.GenoPheno)
#TukeyHSD(bd.GenoPheno)
### Null hypothesis of no difference in dispersion between groups
permutest(bd.GenoPheno)


### PERMANOVA
set.seed(1984)
permanova.SCARP.5N.Dynamic.SoilType <- adonis2(Response_3 ~ Aust_SoilClass*GenoPheno,
                              data=SCaRP.clorpton.df.topsoil.combi5N,
                              method = "mahalanobis", by="terms",
                              #strata = as.numeric(as.character(SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass)),
                              permutations=9999)
permanova.SCARP.5N.Dynamic.SoilType
densityplot(permustats(permanova.SCARP.5N.Dynamic.SoilType))

scarp.dist <- vegdist(x = Response_3, method = "mahalanobis")
bd.Aust_SoilClass <- betadisper(d =scarp.dist,
                                group=SCaRP.clorpton.df.topsoil.combi5N$Aust_SoilClass)
boxplot(bd.Aust_SoilClass)
TukeyHSD(bd.Aust_SoilClass)



### PERMANOVA
set.seed(1984)
permanova.SCARP.5N.Dynamic.Branch <- adonis2(Response_3 ~ Branch*GenoPheno,
                              data=SCaRP.clorpton.df.topsoil.combi5N,
                              method = "mahalanobis", by="terms",
                              #strata = as.numeric(as.character(SCaRP.clorpton.df.topsoil.combi5N$Branch)),
                              permutations=9999)
permanova.SCARP.5N.Dynamic.Branch
densityplot(permustats(permanova.SCARP.5N.Dynamic.Branch))

scarp.dist <- vegdist(x = Response_3, method = "mahalanobis")
bd.Branch <- betadisper(d =scarp.dist,
                        group=SCaRP.clorpton.df.topsoil.combi5N$Branch)
boxplot(bd.Branch)
TukeyHSD(bd.Branch)


# ### 12. Assess changes in pH ------------------------------------------------

### Load ph data
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/TERN_Data")

#ph_4A1 <- read_csv("lab_ph_4A1_data_splined_dsm_extract_prepared.csv")
ph_4B1 <- read_csv("lab_ph_4B1_data_splined_dsm_extract_prepared.csv")
#ph_4C1 <- read_csv("lab_ph_4C1_data_splined_dsm_extract_prepared.csv")
#ph_4NR <- read_csv("lab_ph_4NR_data_splined_dsm_extract_prepared.csv")

setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII")

### Extract pH data on pedogenons & GenoPheno
### Plot the locations of SCaRP samples on a map
ph_4B1 <- as.data.frame(ph_4B1)
ph_4B1.sp <- SpatialPointsDataFrame(data = ph_4B1,
                                    coords = ph_4B1[,4:5],
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
#plot(ph_4A1$Longitude,ph_4A1$Latitude)
plot(ph_4B1$Longitude,ph_4B1$Latitude)
#plot(ph_4C1$Longitude,ph_4C1$Latitude)
#plot(ph_4NR$Longitude,ph_4NR$Latitude)

plot(CLORPTon.GPh.s[[2]])
points(ph_4B1.sp)

ph_4B1.sp <- raster::extract(CLORPTon.GPh.s, ph_4B1.sp, sp=TRUE)
ph_4B1.df <-as.data.frame(ph_4B1.sp)

## Subset those with data on pedogenon and Genosoil/Phenosoil
ph_4B1.df <- ph_4B1.df[!is.na(ph_4B1.df$K1000.18ContVeg),] 
## Subset those with data on pedogenon and Genosoil/Phenosoil
ph_4B1.df <- ph_4B1.df[!is.na(ph_4B1.df$GenoPheno),] 
str(ph_4B1.df) ### We have about 5000 observations

### Set the -9999 into NA
summary(ph_4B1.df)
ph_4B1.df$X0.5.cm <- ifelse(ph_4B1.df$X0.5.cm == -9999, NA,ph_4B1.df$X0.5.cm )
ph_4B1.df$X5.15.cm <- ifelse(ph_4B1.df$X5.15.cm == -9999, NA,ph_4B1.df$X5.15.cm )
ph_4B1.df$X15.30.cm <- ifelse(ph_4B1.df$X15.30.cm == -9999, NA,ph_4B1.df$X15.30.cm )
ph_4B1.df$X30.60.cm <- ifelse(ph_4B1.df$X30.60.cm == -9999, NA,ph_4B1.df$X30.60.cm )
ph_4B1.df$X60.100.cm <- ifelse(ph_4B1.df$X60.100.cm == -9999, NA,ph_4B1.df$X60.100.cm )
ph_4B1.df$X60.100.cm <- ifelse(ph_4B1.df$X60.100.cm < 0, NA,ph_4B1.df$X60.100.cm )
ph_4B1.df$X100.200.cm <- ifelse(ph_4B1.df$X100.200.cm == -9999, NA,ph_4B1.df$X100.200.cm )

### Join data from Branch
str(k1000Veg.out$branch.centroids.ord)
str(ph_4B1.df)
k1000Veg.out$branch.centroids.ord$Centroid <- factor(k1000Veg.out$branch.centroids.ord$Centroid)
k1000Veg.Levels <- as.character(k1000Veg.out$branch.centroids.ord$Centroid)
ph_4B1.df$K1000.18ContVeg <- factor(ph_4B1.df$K1000.18ContVeg, levels = k1000Veg.Levels)
#ph_4B1.df$K1000.18ContVeg <- factor(ph_4B1.df$K1000.18ContVeg, levels = )

ph_4B1.df <- merge(ph_4B1.df,k1000Veg.out$branch.centroids.ord, by.x="K1000.18ContVeg", by.y="Centroid")
str(ph_4B1.df)

### Extract precipitation

#### Create spatialpointsdataframe
ph_4B1.sp <- SpatialPointsDataFrame(data = ph_4B1.df,
                                    coords = ph_4B1.df[,5:6],
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
#plot(ph_4A1$Longitude,ph_4A1$Latitude)
plot(ph_4B1$Longitude,ph_4B1$Latitude)
#plot(ph_4C1$Longitude,ph_4C1$Latitude)
#plot(ph_4NR$Longitude,ph_4NR$Latitude)
plot(CLORPTon.GPh.s[[2]])
points(ph_4B1.sp)
clim_PTA <- raster("C:/Covariates/90m/NSW/clim_PTA_nsw.tif")
ph_4B1.sp <- raster::extract(clim_PTA, ph_4B1.sp, sp=TRUE)
ph_4B1.env <- as.data.frame(ph_4B1.sp)
boxplot(ph_4B1.env$clim_PTA ~ ph_4B1.env$GenoPheno)

ggplot(data=ph_4B1.env) +  
  geom_boxplot(aes(x =as.factor(GenoPheno), y= clim_PTA_nsw)) +   
  # scale_colour_manual(values= viridis_pal(option = "B", direction = 1)(6),
  #                     name="Genoform / Phenoform",
  #                     breaks=c("1", "2", "3", "4", "5", "6"),
  #                     labels=c("Remnant Genoform", "Genoform II",
  #                              "Phenoform Cleared", "Phenoform Forestry", 
  #                              "Phenoform Grazing", "Phenoform Cropping") ) +
  ylab(label="Annual precipitation (mm)") +
  xlab(label="Genoform / Phenoform") +
  theme_bw()+
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6"),
                   c("Remnant Pedogenon", "Quasi-remnant pedogenon",
                     "Pedophenon Cleared", "Pedophenon Forestry", 
                     "Pedophenon Grazing", "Pedophenon Cropping")) 


boxplot(ph_4B1.env$X5.15.cm ~ ph_4B1.env$GenoPheno)
ggplot(data=ph_4B1.env) +  
  geom_point(aes(y =X5.15.cm, x= clim_PTA_nsw,
                 colour=as.factor(GenoPheno),
                 shape=as.factor(GenoPheno)),
             alpha = 0.8, size=2) +   
  scale_colour_manual(values= viridis_pal(option = "B", direction = 1)(6),
                      name="Genoform / Phenoform",
                      breaks=c("1", "2", "3", "4", "5", "6"),
                      labels=c("Remnant Genoform", "Genoform II",
                               "Phenoform Cleared", "Phenoform Forestry", 
                               "Phenoform Grazing", "Phenoform Cropping") ) +
  xlab(label="Annual precipitation (mm)") +
  ylab(label="pH (5-15 cm)")


# Figure 6 in revised article ---------------------------------------------
pch.group <- c(21,22,25,3,7)

pedogenons.pH <- unique(ph_4B1.df$K1000.18ContVeg)#(as.numeric(as.character(ph_4B1.df$K1000.18ContVeg)))
pal.pedogenons.all <- k1000Veg.out$branch.centroids.ord[k1000Veg.out$branch.centroids.ord$Centroid %in% pedogenons.pH,]$colors

pedogenons.pH2 <- unique(as.numeric(as.character(ph_4B1.df$K1000.18ContVeg)))
pal.pedogenons.all2 <- k1000Veg.out$branch.centroids.ord[as.numeric(as.character(k1000Veg.out$branch.centroids.ord$Centroid)) %in% pedogenons.pH2,]$colors
all.equal(pal.pedogenons.all, pal.pedogenons.all2)

set.seed(1234)
plot(ph_4B1.env$clim_PTA ~ jitter(as.numeric(as.character(ph_4B1.env$GenoPheno)),1.5),
     col= pal.pedogenons.all[ph_4B1.env$K1000.18ContVeg],cex=1,pch=1, 
     xlab="Pedogenon / pedophenon class",
     ylab="Annual precipitation (mm)")

plot(ph_4B1.env$X5.15.cm ~ ph_4B1.env$clim_PTA,
     col= pal.pedogenons.all[ph_4B1.env$K1000.18ContVeg],cex=1,
     pch= pch.group[ph_4B1.env$GenoPheno],
     ylab="pH",
     xlab="Annual precipitation (mm)")

### What levels are present in each factor?
unique(ph_4B1.df$K1000.18ContVeg)
### Can I add more?
#k1000Veg.Levels <- as.character(k1000Veg.out$branch.centroids.ord$Centroid)
#ph_4B1.df$K1000.18ContVeg <- factor(ph_4B1.df$K1000.18ContVeg, levels = k1000Veg.Levels)
#levels(ph_4B1.df$K1000.18ContVeg)

#ph_4B1.df$Branch <- as.factor(ph_4B1.df$Branch)
#ph_4B1.df$GenoPheno <- as.factor(ph_4B1.df$GenoPheno)
ph_4B1.df$combi <- paste0(ph_4B1.df$K1000.18ContVeg,"_",ph_4B1.df$GenoPheno)

### Work only with 5-15 cm data
ph_4B1.df <- ph_4B1.df[!is.na(ph_4B1.df$X5.15.cm),] 
length(unique(ph_4B1.df$K1000.18ContVeg))
dim(ph_4B1.df)

### Map them
k1000Veg.out$map.out %>%
  addCircleMarkers(radius = 0.2, color ="black" , opacity = 0.5,
                   lng = ph_4B1.df$Longitude,
                   lat = ph_4B1.df$Latitude) %>%
  addScaleBar(position = "bottomright",
              options = scaleBarOptions(maxWidth = 300, metric = TRUE, imperial = FALSE,
                                        updateWhenIdle = TRUE))

### Can I plot?
ph_4B1.k1000.summary <- ph_4B1.df %>%
  group_by(., K1000.18ContVeg, GenoPheno) %>%
  summarise(., ph.5.15 = mean(X5.15.cm, na.rm=TRUE),
            count=n())
hist(ph_4B1.k1000.summary$count, breaks=50)
summary(ph_4B1.k1000.summary$count)

which.ph <- ph_4B1.k1000.summary %>% 
  filter(., count >=5)
colnames(which.ph) <- c("K1000.18ContVeg", "GenoPheno", "ph.5.15", "count")
which.ph <- as.data.frame(which.ph)
which.ph <- merge(which.ph,k1000Veg.out$branch.centroids.ord, by.x="K1000.18ContVeg", by.y="Centroid")
which.ph$combi <- paste0(which.ph$K1000.18ContVeg,"_",which.ph$GenoPheno)
#levels(which.ph$K1000.18ContVeg)

### subset those individual pH observations from subclasses with at least 5 observations per combinations
ph_4B1.df.5N <- ph_4B1.df[ph_4B1.df$combi %in% which.ph$combi, ]
ph_4B1.df.5N$K1000.18ContVeg <- factor(ph_4B1.df.5N$K1000.18ContVeg)
ph_4B1.df.5N$GenoPheno <- factor(ph_4B1.df.5N$GenoPheno)

### How many have  >= 5?
sum(which.ph$count) ### 3224 observations
length(unique(which.ph$K1000.18ContVeg))

### Plot, with basic R plot functions
which.ph$GenoPheno <- as.numeric(as.character(which.ph$GenoPheno))
# plot(which.ph$GenoPheno, which.ph$ph.5.15, "n", xlab="Genosoil / Phenosoil class",
#      ylab="pH (5-15 cm)")
#which.ph$K1000.18ContVeg <- factor(which.ph$K1000.18ContVeg)

### Plot all the data, with their colors from the pedogenon map
### Get colour scale
#pedogenons.pH <- unique(as.numeric(as.character(ph_4B1.df$K1000.18ContVeg)))
# pedogenons.pH <- unique(ph_4B1.df$K1000.18ContVeg)
# pal.pedogenons.all <- k1000Veg.out$branch.centroids.ord[k1000Veg.out$branch.centroids.ord$Centroid %in% pedogenons.pH,]$colors

# set.seed(1234)
# plot(ph_4B1.df$X5.15.cm ~ jitter(as.numeric(as.character(ph_4B1.df$GenoPheno)),1.5),
#      col= pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],cex=1, pch=1,
#      xlab="Genosoil / Phenosoil class",
#      ylab="pH (5-15 cm)")

set.seed(1234)
plot(ph_4B1.df$X5.15.cm ~ jitter(ph_4B1.df$GenoPheno,1.5),
     col= pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],cex=1, pch=1,
     xlab="",
     ylab="pH (5-15 cm)")

## Now plot the means for those classes with >= 5N
#which.ph$GenoPheno <- factor(which.ph$GenoPheno)
#which.ph$GenoPheno <- as.numeric(as.character(which.ph$GenoPheno))
# plot(which.ph$ph.5.15 ~ jitter(which.ph$GenoPheno, 1),
#        col = pal.pH.3N[which.ph$K1000.18ContVeg],
#        pch = 1, cex=2)
myK <- levels(which.ph$K1000.18ContVeg)
myK <- unique(which.ph$K1000.18ContVeg)
### Color palette
#pedogenons.pH.5N <- unique(as.numeric(as.character(which.ph$K1000.18ContVeg)))
pedogenons.pH.5N <- unique(which.ph$K1000.18ContVeg)
pal.pH.5N <- k1000Veg.out$branch.centroids.ord[k1000Veg.out$branch.centroids.ord$Centroid %in% pedogenons.pH.5N,]$colors
#colvec  <- pal.pH
for(i in 1:length(myK)){
  which.ph.i <- which.ph[which.ph$K1000.18ContVeg== myK[[i]],]
  #all.p.i <- ph_4B1.df[ph_4B1.df$K1000.18ContVeg==myK[[i]],]
  # points(as.numeric(as.character(all.p.i$GenoPheno)), all.p.i$X5.15.cm,
  #       col = unique(which.ph.i$colors),
  #       pch=20, cex=1)
  points(which.ph.i$GenoPheno, which.ph.i$ph.5.15,
         col = which.ph.i$colors[1],
         pch=20, cex=3)
  lines(which.ph.i$GenoPheno, which.ph.i$ph.5.15,
          col = which.ph.i$colors[1],
          lty=1, lwd=1)
} 

### Better in jitter and without lines
### Figure for paper

set.seed(1234)
par(las=1, oma=c(0,0,0,0), mar=c(4, 4, 1, 1) + 0.1)
plot(ph_4B1.df$X5.15.cm ~ jitter(as.numeric(as.character(ph_4B1.df$GenoPheno)),1.5),
     col= pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],cex=1, pch=1,
     xlab="",
     ylab="pH (5-15 cm)")

set.seed(1234)
points( which.ph$ph.5.15 ~ jitter(which.ph$GenoPheno,1.5),
       col = "black",
       bg=  which.ph$colors,
       # k1000Veg.out$branch.centroids.ord$colors[which.ph$K1000.18ContVeg],
       pch=21, cex=2)

### Other plots
# plot(which.ph$K1000.18ContVeg, which.ph$ph.5.15, 
#      pch=21, col="black", bg=  which.ph$colors,
#      #"n",
#      xlab="Pedogenon class", #which.ph$colors,
#      ylab="pH (5-15 cm)")

boxplot(ph_4B1.df$X5.15.cm ~ ph_4B1.df$Branch,
    # pch=21, col=ph_4B1.df$colors,
     #bg= ph_4B1.df$colors,cex=1,
     xlab="Pedogenon branch",
     ylab="pH (5-15 cm)")

plot(ph_4B1.df$X5.15.cm ~ jitter(ph_4B1.df$Branch,1.5),
        pch=21, col=ph_4B1.df$colors,
        bg= ph_4B1.df$colors,cex=1,
        xlab="Pedogenon branch",
        ylab="pH (5-15 cm)")
points(which.ph$ph.5.15 ~jitter(which.ph$Branch,1.5), cex=3,
       pch=21, col="black", bg=  which.ph$colors)

plot(ph_4B1.df$X5.15.cm ~ ph_4B1.df$Branch,
     pch=21, col=ph_4B1.df$colors,
     bg= ph_4B1.df$colors,cex=1,
     xlab="Pedogenon branch",
     ylab="pH (5-15 cm)")

plot(which.ph$Branch, which.ph$ph.5.15,
     pch=21, col="black", bg=  which.ph$colors,
     #"n",
     xlab="Pedogenon branch", #which.ph$colors,
     ylab="pH (5-15 cm)")

plot(as.factor(which.ph$Branch), which.ph$ph.5.15,
     xlab="Pedogenon branch",
     ylab="pH (5-15 cm)")


# ### 13. Assess changes in pH with GLS model -----------------------------
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII")
### lm model
hist((ph_4B1.df$X5.15.cm), breaks=50)
library(AID)
boxcoxnc(ph_4B1.df$X5.15.cm[1:5000],method = "sw", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, alpha = 0.05, verbose = TRUE)
par(mfrow=c(1,1))
ph_4B1.df$K1000.18ContVeg <- factor(ph_4B1.df$K1000.18ContVeg)
ph_4B1.df$GenoPheno <- factor(ph_4B1.df$GenoPheno)

## Let's try luck with this
lm.pH <- lm(formula= X5.15.cm ~ K1000.18ContVeg * GenoPheno, data=ph_4B1.df)
summary(lm.pH)
anova(lm.pH)

mod2 <- lm(formula= X5.15.cm ~ K1000.18ContVeg + GenoPheno, data=ph_4B1.df)
summary(mod2)

# Likelihood ratio test
anova(lm.pH, mod2, test="LRT") ### The test is significant. Therefore I leave the interaction
# Analysis of Variance Table
# 
# Model 1: X5.15.cm ~ K1000.18ContVeg * GenoPheno
# Model 2: X5.15.cm ~ K1000.18ContVeg + GenoPheno
# Res.Df    RSS   Df Sum of Sq  Pr(>Chi)    
# 1   3724 1562.4                             
# 2   4465 1985.9 -741   -423.54 1.568e-10 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
anova( lm.pH,mod2)
AIC(lm.pH, mod2) ## Just the AIC would say the model with interaction is best

colnames(ph_4B1.df)
length(unique(ph_4B1.df$nloc_fid))
summary(lm.pH)
plot(lm.pH)
anova(lm.pH)
acf(lm.pH$residuals) ## No autocorrelation

plot(ph_4B1.df[!is.na(ph_4B1.df$X5.15.cm),]$K1000.18ContVeg, resid(lm.pH))
plot(ph_4B1.df[!is.na(ph_4B1.df$X5.15.cm),]$GenoPheno, resid(lm.pH))
### Heterogeneous variance by GenoPheno and by Pedogenon

summary(lm.pH) ### Adjusted R-squared:  0.5867 
rm(lm.pH, mod2)

ph_4B1.df.5N

# ### GLS model in a subset of observations ---------------------------------------

ph_4B1.df.5N <- ph_4B1.df[ph_4B1.df$combi %in% unique(which.ph$combi), ] ### 3224 obs in 311 combinations
pal.pH.5N <- k1000Veg.out$branch.centroids.ord[k1000Veg.out$branch.centroids.ord$Centroid %in%
                                                 unique(as.numeric(as.character(ph_4B1.df.5N$K1000.18ContVeg))),]$colors
### update factor levels
ph_4B1.df.5N$K1000.18ContVeg <- factor(ph_4B1.df.5N$K1000.18ContVeg) ### 176 levels
ph_4B1.df.5N$GenoPheno <- factor(ph_4B1.df.5N$GenoPheno) ### 5 levels
str(ph_4B1.df.5N)
save(ph_4B1.df.5N, pal.pH.5N, file="ph_4B1.df.5N.RData")

### GLS model only in those groups with at least 5 obs per class
library(nlme)
setwd("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII")
load("ph_4B1.df.5N.RData")

lm.pH <- lm(X5.15.cm ~ K1000.18ContVeg * GenoPheno, data=ph_4B1.df.5N)
plot(lm.pH)
#e0 <- residuals(lm.pH, type = "pearson")
plot(ph_4B1.df.5N$K1000.18ContVeg,resid(lm.pH), col=pal.pH.5N )
### Normality of residuals
hist(resid(lm.pH), breaks=50) ### Seem normal

### Homogeneity of variance
bartlett.test(resid(lm.pH)~ K1000.18ContVeg, data = cbind(ph_4B1.df.5N,resid(lm.pH)) )
# Bartlett test of homogeneity of variances
# 
# data:  resid(lm.pH) by K1000.18ContVeg
# Bartlett's K-squared = 521.33, df = 175, p-value < 2.2e-16
bartlett.test(resid(lm.pH)~ GenoPheno, data = cbind(ph_4B1.df.5N,resid(lm.pH)) )

plot(lm.pH, which = c(1), col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     add.smooth = TRUE, caption = "")
plot(ph_4B1.df.5N$X5.15.cm_t, lm.pH$fitted, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     add.smooth = TRUE, caption = "")
abline(0,1)
boxplot(resid(lm.pH) ~ ph_4B1.df.5N$K1000.18ContVeg, col = pal.pH.5N)
boxplot(resid(lm.pH) ~ ph_4B1.df.5N$GenoPheno)
boxplot(resid(lm.pH) ~ ph_4B1.df.5N$combi)

##Try a transformation
library(AID)
ph_4B1.df.5N$X5.15.cm_t <- boxcoxnc(ph_4B1.df.5N$X5.15.cm,method = "mle", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, alpha = 0.05, verbose = TRUE)$tf.data
lm.pH <- lm(X5.15.cm_t ~ K1000.18ContVeg * GenoPheno, data=ph_4B1.df.5N)
plot(lm.pH)
write.csv(ph_4B1.df.5N, file="ph_4B1.df.5N.csv")
### same problem

### Fit the beyond optimal GLS model
gls.full <- gls(X5.15.cm ~ K1000.18ContVeg * GenoPheno, data=ph_4B1.df.5N)
gls.full <- gls(X5.15.cm ~ K1000.18ContVeg + GenoPheno + K1000.18ContVeg:GenoPheno, data=ph_4B1.df.5N)
# Error in glsEstimate(glsSt, control = glsEstControl) : 
#   computed "gls" fit is singular, rank 312

gls1 <- gls(X5.15.cm ~ I(droplevels(interaction(K1000.18ContVeg,GenoPheno))),
            data=ph_4B1.df.5N,
            na.action = na.omit)
gls2 <- gls(X5.15.cm ~ K1000.18ContVeg + GenoPheno, data=ph_4B1.df.5N)
gls3 <- gls(X5.15.cm ~ K1000.18ContVeg , data=ph_4B1.df.5N)
gls4 <- gls(X5.15.cm ~ GenoPheno, data=ph_4B1.df.5N)
AIC(gls1,gls2,gls3,gls4)

# df      AIC
# gls1 312 7104.160
# gls2 181 6956.717
# gls3 177 6978.263
# gls4   6 8496.956
### Would suggest model without interaction

### Can I check if the interaction is significant?
anova(gls1)
anova(gls2)
anova(update(gls1, method = "ML"),
      update(gls2, method = "ML"))
#                             Model  df      AIC      BIC    logLik   Test  L.Ratio p-value
# update(gls1, method = "ML")     1 312 6714.902 8611.356 -3045.451                        
# update(gls2, method = "ML")     2 181 6664.667 7764.853 -3151.333 1 vs 2 211.7643  <.0001

anova(update(gls1, method = "ML"),
      update(gls3, method = "ML"))
#                             Model  df      AIC      BIC    logLik   Test  L.Ratio p-value
# update(gls1, method = "ML")     1 312 6714.902 8611.356 -3045.451                        
# update(gls3, method = "ML")     2 177 6706.483 7782.356 -3176.241 1 vs 2 261.5806  <.0001

anova(update(gls1, method = "ML"),
      update(gls4, method = "ML"))
#                             Model  df      AIC      BIC    logLik   Test  L.Ratio p-value
# update(gls1, method = "ML")     1 312 6714.902 8611.356 -3045.451                        
# update(gls4, method = "ML")     2   6 8472.883 8509.354 -4230.442 1 vs 2 2369.981  <.0001

anova(update(gls2, method = "ML"),
      update(gls4, method = "ML"))
anova(update(gls2, method = "ML"),
      update(gls3, method = "ML"))
### all the terms are significant
### ### ### ###
plot(gls1, which = c(1), col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     add.smooth = TRUE, caption = "")
plot(ph_4B1.df.5N$X5.15.cm, gls1$fitted, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     add.smooth = TRUE, caption = "")
abline(0,1)
plot(gls1$fitted,resid(gls1), col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     add.smooth = TRUE, caption = "")
abline(0,0)
boxplot(resid(gls1) ~ ph_4B1.df.5N$K1000.18ContVeg, col = pal.pH.5N)
boxplot(resid(gls1) ~ ph_4B1.df.5N$GenoPheno)
boxplot(resid(gls1) ~ ph_4B1.df.5N$combi)
plot(gls1)
rm(gls1,gls2,gls3,gls4)


### 2. Fit GLS model with interactions

# ### 13.2 Fit GLS model with interactions ----------------------------------
### or, equivalently first define the interaction variable:
ph_4B1.df.5N$ii <- droplevels(interaction(ph_4B1.df.5N$K1000.18ContVeg, ph_4B1.df.5N$GenoPheno))
table(ph_4B1.df.5N$ii, ph_4B1.df.5N$K1000.18ContVeg)

e.gls <- gls(X5.15.cm ~ ii, data = ph_4B1.df.5N, na.action = na.omit, method = "REML")
logLik(e.gls)
anova(e.gls)
summary(e.gls)
# Generalized least squares fit by REML
# Model: X5.15.cm ~ ii 
# Data: ph_4B1.df.5N 
#     AIC      BIC   logLik
# 7104.16 8968.965 -3240.08

### Check graphically the residuals
e0 <- residuals(e.gls, type = "normalized")
plot(e.gls, which = c(1), col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],caption = "")
plot(ph_4B1.df.5N$X5.15.cm, e.gls$fitted, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Fitted pH")
abline(0,1)
plot(ph_4B1.df.5N$X5.15.cm, e0, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Normalized residuals")
abline(0,0)
plot(e.gls$fitted,e0, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     xlab= "Fitted pH", ylab="Normalized residuals")
abline(0,0)
boxplot(e0 ~ ph_4B1.df.5N$K1000.18ContVeg, col = pal.pH.5N)
boxplot(e0 ~ ph_4B1.df.5N$GenoPheno)
boxplot(e0 ~ ph_4B1.df.5N$combi)
plot(e.gls)
#save.image("C:/Users/mrom8073/Desktop/USydney/Postdoc/Projects/Genosoils/Output/7.Article_PartII/29102020.RData")

### 3. Fit variance structure

# ### 13.3  Fit variance structure in GLS model ---------------------------

### 3.1 Variance structure in the model without interaction
gls.1.var <- gls(X5.15.cm ~ K1000.18ContVeg + GenoPheno, weights = varIdent(form = ~ 1 | K1000.18ContVeg),
                 data = ph_4B1.df.5N, na.action = na.omit, method = "REML")
logLik(gls.1.var)
anova(gls.1.var)
summary(gls.1.var)
# Generalized least squares fit by REML
# Model: X5.15.cm ~ ii 
# Data: ph_4B1.df.5N 
#     AIC      BIC   logLik
# 7104.16 8968.965 -3240.08

### Check graphically the residuals
e0 <- residuals(gls.1.var, type = "normalized")
hist(resid(gls.1.var), breaks=30)
plot(gls.1.var, which = c(1), col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],caption = "")
plot(ph_4B1.df.5N$X5.15.cm, gls.1.var$fitted, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Fitted pH")
abline(0,1)
plot(ph_4B1.df.5N$X5.15.cm, e0, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Normalized residuals")
abline(0,0)
plot(gls.1.var$fitted,e0, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     xlab= "Fitted pH", ylab="Normalized residuals")
abline(0,0)
boxplot(e0 ~ ph_4B1.df.5N$K1000.18ContVeg, col = pal.pH.5N)
boxplot(e0 ~ ph_4B1.df.5N$GenoPheno)
boxplot(e0 ~ ph_4B1.df.5N$combi)

# gls.2.var2 <- gls(X5.15.cm ~ K1000.18ContVeg + GenoPheno, weights = varIdent(form = ~ 1 | GenoPheno),
#                  data = ph_4B1.df.5N, na.action = na.omit, method = "REML")
# 
# AIC(gls.2, gls.2.var2, gls.2.var)
# #             df      AIC
# # gls.2      181 6956.717
# # gls.2.var2 185 6923.190
# # gls.2.var  356 6756.928 ### Preferred the model with variance by month

### Clean
rm(e0,e.gls, e.gls.2.var2, lm.pH)

### Model with interaction no variance structure (like lm)
#gls.1 <- gls(X5.15.cm ~ ii, data = ph_4B1.df.5N, na.action = na.omit, method = "REML")

### Try to fit it for the model with interaction AND variance structure
gls.2.var <- gls(X5.15.cm ~ ii, data = ph_4B1.df.5N,
                 weights = varIdent(form = ~ 1 | K1000.18ContVeg),
                 method = "REML", na.action = na.omit)
gls.1.var
### Check graphically the residuals
e0 <- residuals(gls.2.var, type = "normalized")
hist(resid(gls.2.var), breaks=30)
plot(gls.2.var, which = c(1), col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],caption = "")
plot(ph_4B1.df.5N$X5.15.cm, gls.2.var$fitted, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Fitted pH")
abline(0,1)
plot(ph_4B1.df.5N$X5.15.cm, e0, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Normalized residuals")
abline(0,0)
plot(gls.2.var$fitted,e0, col = pal.pH.5N[ph_4B1.df.5N$K1000.18ContVeg],
     xlab= "Fitted pH", ylab="Normalized residuals")
abline(0,0)
boxplot(e0 ~ ph_4B1.df.5N$K1000.18ContVeg, col = pal.pH.5N)
boxplot(e0 ~ ph_4B1.df.5N$GenoPheno)
boxplot(e0 ~ ph_4B1.df.5N$combi)
plot(gls.2.var)
qqnorm(gls.2.var, abline = c(0,1)) ### Pretty good

### Check the significance of the interaction
AIC(gls.2.var,gls.1.var) ## AIC alone says the better model does not have interaction

anova(update(gls.1.var, method = "ML"),
      update(gls.2.var, method = "ML"))
# Model  df      AIC      BIC    logLik   Test  L.Ratio p-value
# update(e.var.gls, method = "ML")     1 487 6476.749 9436.920 -2751.375                        
# update(gls.2.var, method = "ML")     2 356 6413.234 8577.136 -2850.617 1 vs 2 198.4845   1e-04

save(gls.2.var, gls.1.var, file="gls.var.RData")

#### Extract the estimated means
library(multcomp)
library(emmeans)
summary(ref_grid(gls.2.var))
ph.emm.ii <- emmeans(gls.2.var, "ii")

### Separate into the two factors
ph.emm.ii.2 <- separate(data=as.data.frame(ph.emm.ii), col = ii,
                        into=c("Pedogenon", "GenoPheno"), 
                        remove = FALSE)

### Plot it
myK <- sort(as.numeric(unique(ph.emm.ii.2$Pedogenon)))
myK <- as.character(myK)
### Color palette
pal.pH.emmeans <- k1000Veg.out$branch.centroids.ord[k1000Veg.out$branch.centroids.ord$Centroid %in% as.numeric(myK),]$colors
pal.pH.emmeans <- rainbow(n= length(myK))
pal.pH.5N
ph.emm.ii.2$GenoPheno <- as.numeric(ph.emm.ii.2$GenoPheno)
#colvec  <- pal.pH
par(mfrow=c(1,1))
plot(as.numeric(as.character(ph.emm.ii.2$GenoPheno)), ph.emm.ii.2$emmean,
     "n", xlab="Genoform / Phenoform class",
     ylab=" Estimated marginal means pH (5-15 cm)")
for(i in 1:length(myK)){
  which.ph.i <- ph.emm.ii.2[ph.emm.ii.2$Pedogenon== myK[[i]],]
  #all.p.i <- ph_4B1.df[ph_4B1.df$K1000.18ContVeg==myK[[i]],]
  # points(as.numeric(as.character(all.p.i$GenoPheno)), all.p.i$X5.15.cm,
  #       col = unique(which.ph.i$colors),
  #       pch=20, cex=1)
  points(which.ph.i$GenoPheno, which.ph.i$emmean,
         col = pal.pH.5N[[i]],
         pch=20, cex=3)
  lines(which.ph.i$GenoPheno, which.ph.i$emmean,
        col = pal.pH.5N[[i]],
        lty=1, lwd=2)
} 

ph.emm.ii.2$Pedogenon <- factor(ph.emm.ii.2$Pedogenon, levels = myK)

set.seed(1234)
plot(ph.emm.ii.2$emmean ~ jitter(as.numeric(as.character(ph.emm.ii.2$GenoPheno)),1),
     col= pal.pH.5N[ph.emm.ii.2$Pedogenon],cex=3, pch=20,
     xlab="Genosoil / Phenosoil class",
     ylab="pH (5-15 cm)")

ph.emm.ii.size <-eff_size(ph.emm.ii, sigma = sigma(gls.2.var), edf = 3224)

ph.emm.ii.size.df <- as.data.frame(ph.emm.ii.size)
ph.emm.ii.size.df <- separate(data=ph.emm.ii.size.df, contrast,
                                into=c("Subclass1", "Subclass2"), sep = " - ",
                                remove = FALSE)
ph.emm.ii.size.df <- separate(data=ph.emm.ii.size.df, col="Subclass1",
                              into=c("Pedogenon1", "GenoPheno1"),convert=TRUE,
                              remove = FALSE)
ph.emm.ii.size.df <- separate(data=ph.emm.ii.size.df, col="Subclass2",
                              into=c("Pedogenon2", "GenoPheno2"),convert=TRUE,
                              remove = FALSE)

### pairwise differences
ph.emm.ii.pairs <- pairs(ph.emm.ii)
ph.emm.ii.pairs.df <- as.data.frame(ph.emm.ii.pairs)

library(tidyverse)
ph.emm.ii.pairs.df2 <- separate(data=ph.emm.ii.pairs.df, contrast,
                                into=c("Subclass1", "Subclass2"), sep = " - ",
                                remove = FALSE)
str(ph.emm.ii.pairs.df2)
#ph.emm.ii.pairs.df2$Subclass1 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass1)
#ph.emm.ii.pairs.df2$Subclass2 <- gsub(pattern=" ", replacement = "", x = ph.emm.ii.pairs.df2$Subclass2)
ph.emm.ii.pairs.df3 <- separate(data=ph.emm.ii.pairs.df2, col = Subclass1,
                                into=c("Pedogenon1", "GenoPheno1"), 
                                remove = FALSE)
ph.emm.ii.pairs.df4 <- separate(data=ph.emm.ii.pairs.df3, col = Subclass2,
                                into=c("Pedogenon2", "GenoPheno2"), 
                                remove = FALSE)

rm(ph.emm.ii.pairs.df4.sig,ph.emm.ii.pairs.df3)
ph.emm.ii.pairs.df4.sig <- ph.emm.ii.pairs.df4[ph.emm.ii.pairs.df4$p.value < 0.05,]

### By pedogenon
ph.emm.ii.pairs.df.K1000 <- ph.emm.ii.pairs.df4[ph.emm.ii.pairs.df4$Pedogenon1 == ph.emm.ii.pairs.df4$Pedogenon2,]
summary(abs(ph.emm.ii.pairs.df.K1000$estimate), breaks=30)

ph.emm.ii.pairs.df.K1000$GenoPheno1 <-  as.numeric(ph.emm.ii.pairs.df.K1000$GenoPheno1)
ph.emm.ii.pairs.df.K1000$GenoPheno2 <-  as.numeric(ph.emm.ii.pairs.df.K1000$GenoPheno2)
ph.emm.ii.pairs.df.K1000$Pedogenon1 <-  as.numeric(ph.emm.ii.pairs.df.K1000$Pedogenon1)
ph.emm.ii.pairs.df.K1000$Pedogenon2 <-  as.numeric(ph.emm.ii.pairs.df.K1000$Pedogenon2)

ph.emm.ii.pairs.df.K1000 <- arrange(ph.emm.ii.pairs.df.K1000, Pedogenon1, GenoPheno1, GenoPheno2  )
write.csv(ph.emm.ii.pairs.df.K1000,file="ph.emm.ii.pairs.df.K1000.csv")



ph.emm.ii.pairs.df.K1000.summary <- ph.emm.ii.pairs.df.K1000 %>%
  group_by(., GenoPheno1, GenoPheno2 ) %>%
  summarize(.,
            min.diff = round(min(estimate),3),
            q25.diff = round(quantile(estimate,0.25 ),3),
            mean.diff = round(mean(estimate),3),
            q75.diff = round(quantile(estimate,0.75 ),3),
            max.diff = round(max(estimate),3),
           # mean.se = round(mean(SE),3),
            sd_diff = round(sd(estimate),3),
            count=n())

ph.emm.ii.pairs.df.K1000.summary <- as.data.frame(ph.emm.ii.pairs.df.K1000.summary)
write.csv(ph.emm.ii.pairs.df.K1000.summary, file="ph.emm.ii.pairs.df.K1000.summary.csv")

ph.emm.ii.pairs.df.K1000.sig <- ph.emm.ii.pairs.df.K1000[ph.emm.ii.pairs.df.K1000$p.value<0.05,]
ph.emm.ii.pairs.df.K1000.larger <- ph.emm.ii.pairs.df.K1000[abs(ph.emm.ii.pairs.df.K1000$estimate) >0.4,]

ph.emm.ii.pairs.df.K1000.larger.classes <- ph.emm.ii.pairs.df.K1000.larger %>%
  group_by(., GenoPheno1, GenoPheno2) %>%
  summarise(., mean_diff = mean(estimate),
            sd_diff = sd(estimate),
            count=n())
ph.emm.ii.pairs.df.K1000.larger.classes <- as.data.frame(ph.emm.ii.pairs.df.K1000.larger.classes)


# ### 14. GLS model with branch -----------------------------------------------

ph_4B1.df$GenoPheno <- factor(ph_4B1.df$GenoPheno)
ph_4B1.df$Branch.ii <- droplevels(interaction(ph_4B1.df$Branch, ph_4B1.df$GenoPheno))
table(ph_4B1.df$Branch.ii , ph_4B1.df$Branch)

### 3.1 Variance structure in the model without interaction
gls.1.var <- gls(X5.15.cm ~ Branch + GenoPheno, weights = varIdent(form = ~ 1 | Branch),
                 data = ph_4B1.df, na.action = na.omit, method = "REML")
logLik(gls.1.var)
anova(gls.1.var)
summary(gls.1.var)
ph_4B1.df$K1000.18ContVeg <- factor(ph_4B1.df$K1000.18ContVeg)

### Check graphically the residuals
e0 <- residuals(gls.1.var, type = "normalized")
hist(resid(gls.1.var), breaks=30)
plot(gls.1.var, which = c(1), caption = "")
plot(ph_4B1.df$X5.15.cm, gls.1.var$fitted, col = pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Fitted pH")
abline(0,1)
plot(ph_4B1.df$X5.15.cm, e0, col = pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Normalized residuals")
abline(0,0)
plot(gls.1.var$fitted,e0, col = pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],
     xlab= "Fitted pH", ylab="Normalized residuals")
abline(0,0)
boxplot(e0 ~ ph_4B1.df$Branch)
boxplot(e0 ~ ph_4B1.df$GenoPheno)
boxplot(e0 ~ ph_4B1.df$Branch.Branch.ii)

### Clean
rm(e0,e.gls, e.gls.2.var2, lm.pH)

### Model with interaction no variance structure (like lm)
#gls.1 <- gls(X5.15.cm ~ ii, data = ph_4B1.df, na.action = na.omit, method = "REML")

### Try to fit it for the model with interaction AND variance structure
gls.2.var <- gls(X5.15.cm ~ Branch.ii, data = ph_4B1.df,
                 weights = varIdent(form = ~ 1 | Branch),
                 method = "REML", na.action = na.omit)
gls.1.var
anova(gls.1.var)
anova(gls.2.var)
### Check graphically the residuals
e0 <- residuals(gls.2.var, type = "normalized")
hist(resid(gls.2.var), breaks=30)
plot(gls.2.var, which = c(1), col = pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],caption = "")
plot(ph_4B1.df$X5.15.cm, gls.2.var$fitted, col = pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Fitted pH")
abline(0,1)
plot(ph_4B1.df$X5.15.cm, e0, col = pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],
     xlab= "Observed pH", ylab="Normalized residuals")
abline(0,0)
plot(gls.2.var$fitted,e0, col = pal.pedogenons.all[ph_4B1.df$K1000.18ContVeg],
     xlab= "Fitted pH", ylab="Normalized residuals")
abline(0,0)
boxplot(e0 ~ ph_4B1.df$K1000.18ContVeg, col = pal.pedogenons.all)
boxplot(e0 ~ ph_4B1.df$GenoPheno)
boxplot(e0 ~ ph_4B1.df$Branch.ii)
boxplot(e0 ~ ph_4B1.df$Branch)
plot(gls.2.var)
qqnorm(gls.2.var, abline = c(0,1)) ### Pretty good

### Check the significance of the interaction
AIC(gls.2.var,gls.1.var) ## AIC alone says the better model does not have interaction

anova(update(gls.1.var, method = "ML"),
      update(gls.2.var, method = "ML"))

#### Extract the estimated means
library(multcomp)
library(emmeans)
summary(ref_grid(gls.2.var))
ph.emm.Branch.ii <- emmeans(gls.2.var, "Branch.ii")

### Separate into the two factors
ph.emm.Branch.ii.2 <- tidyr::separate(data=as.data.frame(ph.emm.Branch.ii), col = Branch.ii,
                        into=c("Branch", "GenoPheno"), 
                        remove = FALSE)

### Plot it
myK <- sort(as.numeric(unique(ph.emm.Branch.ii.2$Branch)))
myK <- as.character(myK)
### Color palette
pal.pH.emmeans <- viridis_pal(option = "D", direction = 1)(length(myK))
ph.emm.Branch.ii.2$GenoPheno <- as.numeric(ph.emm.Branch.ii.2$GenoPheno)

par(mfrow=c(1,1))
plot(as.numeric(as.character(ph.emm.Branch.ii.2$GenoPheno)), ph.emm.Branch.ii.2$emmean,
     "n", xlab="Pedogenon / Pedophenon class",
     ylab=" Estimated marginal means pH (5-15 cm)")
for(i in 1:length(myK)){
  which.ph.i <- ph.emm.Branch.ii.2[ph.emm.Branch.ii.2$Branch== myK[[i]],]
  points(which.ph.i$GenoPheno, which.ph.i$emmean,
         col = pal.pH.emmeans[[i]],
         pch=20, cex=3)
  lines(which.ph.i$GenoPheno, which.ph.i$emmean,
         col = pal.pH.emmeans[[i]],
         pch=20, cex=3)
} 

ggplot()+ geom_boxplot(data=ph_4B1.df,
                       aes(x=GenoPheno, y=X5.15.cm, fill=GenoPheno))+ facet_wrap(~Branch)
  


ph.emm.Branch.ii.2$Branch <- factor(ph.emm.Branch.ii.2$Branch, levels = myK)

set.seed(1234)
plot(ph.emm.Branch.ii.2$emmean ~ jitter(as.numeric(as.character(ph.emm.Branch.ii.2$GenoPheno)),1),
     col= pal.pH.emmeans[ph.emm.Branch.ii.2$Branch],cex=3, pch=20,
     xlab="Pedogenon / Pedophenon class",
     ylab="pH (5-15 cm)")

ph.emm.Branch.ii.size <-eff_size(ph.emm.Branch.ii, sigma = sigma(gls.2.var), edf = 3224)

ph.emm.Branch.ii.size.df <- as.data.frame(ph.emm.Branch.ii.size)
ph.emm.Branch.ii.size.df <- separate(data=ph.emm.Branch.ii.size.df, contrast,
                              into=c("Subclass1", "Subclass2"), sep = " - ",
                              remove = FALSE)
ph.emm.Branch.ii.size.df <- separate(data=ph.emm.Branch.ii.size.df, col="Subclass1",
                              into=c("Pedogenon1", "GenoPheno1"),convert=TRUE,
                              remove = FALSE)
ph.emm.Branch.ii.size.df <- separate(data=ph.emm.Branch.ii.size.df, col="Subclass2",
                              into=c("Pedogenon2", "GenoPheno2"),convert=TRUE,
                              remove = FALSE)

### pairwise differences
ph.emm.Branch.ii.pairs <- pairs(ph.emm.Branch.ii)
ph.emm.Branch.ii.pairs.df <- as.data.frame(ph.emm.Branch.ii.pairs)

library(tidyverse)
ph.emm.Branch.ii.pairs.df2 <- separate(data=ph.emm.Branch.ii.pairs.df, contrast,
                                into=c("Subclass1", "Subclass2"), sep = " - ",
                                remove = FALSE)
str(ph.emm.Branch.ii.pairs.df2)
#ph.emm.Branch.ii.pairs.df2$Subclass1 <- gsub(pattern=" ", replacement = "", x = ph.emm.Branch.ii.pairs.df2$Subclass1)
#ph.emm.Branch.ii.pairs.df2$Subclass2 <- gsub(pattern=" ", replacement = "", x = ph.emm.Branch.ii.pairs.df2$Subclass2)
ph.emm.Branch.ii.pairs.df3 <- separate(data=ph.emm.Branch.ii.pairs.df2, col = Subclass1,
                                into=c("Pedogenon1", "GenoPheno1"), 
                                remove = FALSE)
ph.emm.Branch.ii.pairs.df4 <- separate(data=ph.emm.Branch.ii.pairs.df3, col = Subclass2,
                                into=c("Pedogenon2", "GenoPheno2"), 
                                remove = FALSE)

rm(ph.emm.Branch.ii.pairs.df4.sig,ph.emm.Branch.ii.pairs.df3)
ph.emm.Branch.ii.pairs.df4.sig <- ph.emm.Branch.ii.pairs.df4[ph.emm.Branch.ii.pairs.df4$p.value < 0.05,]

### By pedogenon
ph.emm.Branch.ii.pairs.df.K1000 <- ph.emm.Branch.ii.pairs.df4[ph.emm.Branch.ii.pairs.df4$Pedogenon1 == ph.emm.Branch.ii.pairs.df4$Pedogenon2,]
summary(abs(ph.emm.Branch.ii.pairs.df.K1000$estimate), breaks=30)

ph.emm.Branch.ii.pairs.df.K1000$GenoPheno1 <-  as.numeric(ph.emm.Branch.ii.pairs.df.K1000$GenoPheno1)
ph.emm.Branch.ii.pairs.df.K1000$GenoPheno2 <-  as.numeric(ph.emm.Branch.ii.pairs.df.K1000$GenoPheno2)
ph.emm.Branch.ii.pairs.df.K1000$Pedogenon1 <-  as.numeric(ph.emm.Branch.ii.pairs.df.K1000$Pedogenon1)
ph.emm.Branch.ii.pairs.df.K1000$Pedogenon2 <-  as.numeric(ph.emm.Branch.ii.pairs.df.K1000$Pedogenon2)

ph.emm.Branch.ii.pairs.df.K1000 <- arrange(ph.emm.Branch.ii.pairs.df.K1000, Pedogenon1, GenoPheno1, GenoPheno2  )
write.csv(ph.emm.Branch.ii.pairs.df.K1000,file="ph.emm.Branch.ii.pairs.df.K1000.csv")



ph.emm.Branch.ii.pairs.df.K1000.summary <- ph.emm.Branch.ii.pairs.df.K1000 %>%
  group_by(., GenoPheno1, GenoPheno2 ) %>%
  summarize(.,
            min.diff = round(min(estimate),3),
            q25.diff = round(quantile(estimate,0.25 ),3),
            mean.diff = round(mean(estimate),3),
            q75.diff = round(quantile(estimate,0.75 ),3),
            max.diff = round(max(estimate),3),
            # mean.se = round(mean(SE),3),
            sd_diff = round(sd(estimate),3),
            count=n())

ph.emm.Branch.ii.pairs.df.K1000.summary <- as.data.frame(ph.emm.Branch.ii.pairs.df.K1000.summary)
write.csv(ph.emm.Branch.ii.pairs.df.K1000.summary, file="ph.emm.Branch.ii.pairs.df.K1000.summary.csv")

ph.emm.Branch.ii.pairs.df.K1000.sig <- ph.emm.Branch.ii.pairs.df.K1000[ph.emm.Branch.ii.pairs.df.K1000$p.value<0.05,]
ph.emm.Branch.ii.pairs.df.K1000.larger <- ph.emm.Branch.ii.pairs.df.K1000[abs(ph.emm.Branch.ii.pairs.df.K1000$estimate) >0.4,]

ph.emm.Branch.ii.pairs.df.K1000.larger.classes <- ph.emm.Branch.ii.pairs.df.K1000.larger %>%
  group_by(., GenoPheno1, GenoPheno2) %>%
  summarise(., mean_diff = mean(estimate),
            sd_diff = sd(estimate),
            count=n())
ph.emm.Branch.ii.pairs.df.K1000.larger.classes <- as.data.frame(ph.emm.Branch.ii.pairs.df.K1000.larger.classes)



#############################################################################################################################