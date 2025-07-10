# A Novel Approach to Field Data Augmentation with Remote Sensing and Machine Learning in Rangelands
# R Code by Javier Osorio Leyton & Hailey E. Schmidt
# Part I: Creating Feature Set, Model Development, and Accuracy Assessment
# -------------------------------------------------------

#### Data Loading and Set-Up ####
# remove all objects in the workspace
rm(list = ls(all.names = TRUE)) # will clear all objects, including hidden objects
cat("\f"); gc() # free up memory and report memory usage

# load required libraries (Note: if these packages are not installed, then install them first and then load)
library(raster)       # raster: for the manipulation of raster data
library(caret)        # caret: for the machine learning algorithms
library(sp)           # sp: for the manipulation of spatial objects
library(sf)           # sf: Simple Features for R
library(nnet)         # nnet: Artificial Neural Network
library(cluster)      # cluster: "Finding Groups in Data": Cluster Analysis Extended
library(randomForest) # randomForest: Random Forest 
library(kernlab)      # kernlab: Support Vector Machines
library(e1071)        # e1071: provides miscellaneous functions requiered by the caret package
library(tidyr)        # tidyr: Tidy Messy Data
library(doSNOW)       # doSNOW: Foreach provides a parallel backend for the %dopar% function using the snow packag
library(corrplot)     # corrplot: A visualization of a correlation matrix
library(glcm)         # glcm: Calculate Textures from Grey-Level Co-Occurrence Matrices (GLCMs)
library(XML)          # XML: Tools for Parsing and Generating XML Within R and S-Plus
library(xml2)         # xml2: Parse XML
library(dplyr)        # dplyr: data manipulation

# collect files and folders
PathDir <- "C:/Users/hailey.schmidt/OneDrive - Texas A&M AgriLife/Desktop/Gis_Rs/Planet" 
WorkDir <- paste(PathDir, "/WorkDir", sep = "")
RastrDir <- paste(PathDir, "/RasterDir", sep = "")
ShapeDir <- paste(PathDir, "/ShapeDir", sep = "")
ImageDir <- paste(PathDir, "/ImageDir/p_2018_01_31", sep = "")
list.files(path = PathDir, full.names = FALSE, recursive = FALSE) # lists Project
list.files(path = RastrDir, full.names = FALSE, recursive = FALSE) # lists Rasters
list.files(path = ShapeDir, full.names = FALSE, recursive = FALSE) # lists Shapefiles
list.files(path = ImageDir, full.names = FALSE, recursive = FALSE) # lists Images

# setup working directory
setwd(WorkDir); getwd()

#### Get satellite data ####
# Open XML file
# https://rpubs.com/Howetowork/499292
xlmName <- paste(ImageDir, "20180131_164252_102f_3B_AnalyticMS_metadata_clip.xml", sep = '/')
xml_1 <- xmlParse(xlmName)
class(xml_1)  # provides the class of the parsed file
xmlRoot(xml_1)

# looking at the data
xmltop = xmlRoot(xml_1) # gives content of root
xmlName(xmltop) # give name of node, PubmedArticleSet
xmlSize(xmltop) # how many children in node
xmlName(xmltop[[1]]) #name of root's children

# raster name
# <acquisition date>_<acquisition time>_<satellite_id>_<productLevel>_<bandProduct>.<extension> 
# example: 20230207_143613_03_241c_3B_AnalyticMS_SR_8b.tif
rasName <- xmlValue(xmltop[[1]][[1]][[1]])
print(rasName)
dateCOll <- stringr::str_sub(rasName, start = 1, end = 8)
baseName <- paste(Sys.Date(), sep = "_")

# gives content of root
xmltop = xmlRoot(xml_1) 

# collect reflectance coefficients
bandNumber <- NULL
reflectanceCoefficient <- NULL
for (i in 1:8) {
  bandNumber[i] <- xmlValue(xmltop[[5]][[1]][[i+5]][[1]])
  reflectanceCoefficient[i] <- xmlValue(xmltop[[5]][[1]][[i + 5]][[5]])    
}
vecPlanetCoef <- as.numeric(reflectanceCoefficient)
str(vecPlanetCoef); vecPlanetCoef[5]

# open Planet image 
imgName <- paste(ImageDir, "composite.tif", sep = '/') 

# load the Planet stack of the study area
s1data = stack(imgName)

# do this if your image has digital numbers and you need surface reflectance
# bandNumber = 1 - Coastal blue: 431 - 452 nm
s1data[[1]] = s1data[[1]] / 10000        
# bandNumber = 2 - Blue: 465 - 515 nm
s1data[[2]] = s1data[[2]] / 10000
# bandNumber = 3 - Green I:  513 - 549 nm
s1data[[3]] = s1data[[3]] / 10000
# bandNumber = 4 - Green:  547 - 583 nm
s1data[[4]] = s1data[[4]] / 10000

# stack the bands into one composite
names(s1data) = c('blue','green','red','nir')
s1data
nrow(s1data); ncol(s1data)
summary(s1data$blue)

# adjust extent to remove areas of no data (NA)
s2data <- crop(s1data, extent(x = s1data, r1 = 50, r2 = nrow(s1data) - 50, c1 = 50, c2 = ncol(s1data) - 50))
s2data
nrow(s2data); ncol(s2data)
summary(s2data$blue)
plot(s2data, colNA="red")

# save the raster
raster::writeRaster(s2data$blue, 
                    filename = paste(baseName, "s2data.tif", sep ="_"), 
                    datatype = 'FLT4S', # Integer: 'INT2S' or Float: 'FLT4S' (default)
                    format = "GTiff",
                    overwrite = TRUE)  

#### Calculate topographic variables ####
# load DEM LiDAR 1m 
demName <- paste(RastrDir, "martranch_dem1m.tif", sep = '/') 
rstDem1m <- raster(demName)    
rstDem1m
plot(rstDem1m)

# downscale to 3m
tmpDem3m <-  aggregate(x = rstDem1m, fact = 3, fun = mean, expand = FALSE, na.rm = TRUE)
tmpDem3m
plot(tmpDem3m)

# crop raster to reduce working extent
tmpDem3Crop <- crop(x = tmpDem3m,  
                    y = extent(s2data$blue)) 
plot(tmpDem3Crop)

# resample DTM raster to match DSM data
tifDem3m <- resample(x = tmpDem3Crop, 
                     y = s2data$blue, 
                     method = "ngb")
tifDem3m
plot(tifDem3m)
names(tifDem3m) <- c("dem")

# save the raster
raster::writeRaster(tifDem3m, 
                    filename = paste(baseName, "dem3m.tif", sep ="_"), 
                    datatype = 'FLT4S', # Integer: 'INT2S' or Float: 'FLT4S' (default)
                    format = "GTiff",
                    overwrite = TRUE) 

# calculate slope
tifSlope3m = terrain(x = tifDem3m, opt = 'slope', units = 'degrees', neighbors = 8, overwrite = TRUE)        
tifSlope3m
plot(tifSlope3m)
names(tifSlope3m) <- c("slope")

# save the raster
raster::writeRaster(tifSlope3m, 
                    filename = paste(baseName, "slope3m.tif", sep ="_"), 
                    datatype = 'FLT4S', # Integer: 'INT2S' or Float: 'FLT4S' (default)
                    format = "GTiff",
                    overwrite = TRUE) 

# calculate aspect
tifAspect3m = terrain(x = tifDem3m, opt = 'aspect', units = 'degrees', neighbors = 8, overwrite = TRUE)        
tifAspect3m
plot(tifAspect3m) 
names(tifAspect3m) <- c("aspect")

# save the raster
raster::writeRaster(tifAspect3m, 
                    filename = paste(baseName, "aspect3m.tif", sep ="_"), 
                    datatype = 'FLT4S', # Integer: 'INT2S' or Float: 'FLT4S' (default)
                    format = "GTiff",
                    overwrite = TRUE) 

#### POLARIS Soils Data ####
## Variable - Description,Units
# silt - silt percentage, %
# sand - sand percentage, %
# clay - clay percentage, %
# bd - bulk density, g/cm3
# theta_s - saturated soil water content, m3/m3
# theta_r - residual soil water content, m3/m3
# ksat - saturated hydraulic conductivity, log10(cm/hr)
# ph - soil pH in H2O, N/A
# om - organic matter, log10(%)
# 
## Depth from surface
# 1. 0-5 cm
# 2. 5-15 cm
# 3. 15-30 cm
# 4. 30-60 cm
# 5. 60-100 cm
# 6. 100-200 cm
# 
## Statistics provided per layer and variable:
# mean - Arithmetic mean
# mode - Mode
# p50 - Median
# p5 - 5th percentile
# p95 - 95th percentile

## open soils data 
siltName <- paste(RastrDir, "silt_0_5_lat3031_lon-100-99.tif", sep = '/')
sandName <- paste(RastrDir, "sand_0_5_lat3031_lon-100-99.tif", sep = '/')
clayName <- paste(RastrDir, "clay_0_5_lat3031_lon-100-99.tif", sep = '/')
ksatName <- paste(RastrDir, "ksat_0_5_lat3031_lon-100-99.tif", sep = '/')
phName <- paste(RastrDir, "ph_0_5_lat3031_lon-100-99.tif", sep = '/')
omName <- paste(RastrDir, "om_0_5_lat3031_lon-100-99.tif", sep = '/')
bdName <- paste(RastrDir, "bd_0_5_lat3031_lon-100-99.tif", sep = '/')
wpName <- paste(RastrDir, "theta_r_0_5_lat3031_lon-100-99.tif", sep = '/')
fcName <- paste(RastrDir, "theta_s_0_5_lat3031_lon-100-99.tif", sep = '/')

rstSilt <- raster(siltName)
rstSand <- raster(sandName)
rstClay <- raster(clayName)
rstKsat <- raster(ksatName)
rstPh <- raster(phName)
rstOm <- raster(omName)
rstBd <- raster(bdName)
rstWp <- raster(wpName)
rstFc <- raster(fcName)

# load the stack of the study area
sol1data = stack(rstSilt, rstSand, rstClay, rstKsat, rstPh, rstOm, rstBd, rstWp, rstFc)
sol1data

# rename the layers of the stack based on previously saved information
names(sol1data) = c('silt', 'sand', 'clay', 'shc','ph', 'om', 'bd', 'wp', 'fc')
sol1data

# reproject soil raster to match imagery
sol1dataProj <- projectRaster(sol1data,
                              crs = crs(s2data$blue),
                              res = 30, 
                              method = 'ngb')
sol1dataProj

# downscale to 3m to match the imagery resolution
tmpSol3m <-  disaggregate(x = sol1dataProj, fact = 10, method = '')
tmpSol3m
plot(tmpSol3m)
summary(tmpSol3m)

# crop raster to reduce working extent
tmpSol3Crop <- crop(x = tmpSol3m,  
                    y = extent(s2data$blue)) 
tmpSol3Crop
plot(tmpSol3Crop)

# resample DTM raster to match DSM data
tifSol3m <- resample(x = tmpSol3Crop, 
                     y = s2data$blue, 
                     method = "ngb")
tifSol3m[[4]] = 10^(tifSol3m[[4]]) # ksat is in log10 
tifSol3m[[6]] = 10^(tifSol3m[[6]]) # OM is in log10 
tifSol3m
plot(tifSol3m)
names(tifSol3m) <- c("silt","sand","clay","shc","ph","om","bd","wp","fc")

# save the raster
raster::writeRaster(tifSol3m, 
                    filename = paste(baseName, "soil3m.tif", sep ="_"), 
                    datatype = 'FLT4S', # Integer: 'INT2S' or Float: 'FLT4S' (default)
                    format = "GTiff",
                    overwrite = TRUE) 

#### calculate vegetation indices ####
# Calculates the Difference Vegetation Index (DVI) DVI = NIR – Red
tifDvi = s2data$nir - s2data$red
tifDvi                   # View raster structure
plot(tifDvi)             # plot the raster      

# Caculate Enhanced Vegetation Index (EVI) EVI = 2.5 * (NIR - Red)/(NIR + 6 * Red - 7.5 * Blue + 1)
tifEvi <- 2.5 * (s2data$nir - s2data$red) / (s2data$nir + (6 * s2data$red) - (7.5 * s2data$blue) + 1)
tifEvi                   # View raster structure
plot(tifEvi)             # plot the raster  

# Calculates the Green Atmospherically Resistant Index (GARI): GARI = NIR – [Green - gamma(Blue - Red) / NIR – [Green + gamma(Blue - Red)
tifGari = (s2data$nir - (s2data$green - 1.7 * (s2data$blue - s2data$red))) / (s2data$nir + (s2data$green - 1.7 * (s2data$blue - s2data$red)))
tifGari                   # View raster structure
plot(tifGari)             # plot the raster      

# Calculates the Green Leaf Index (GLI): GARI = (Green - Red) + (Green - Blue) / (2 * Green) + Red + Blue    
tifGli = ((s2data$green - s2data$red) + (s2data$green - s2data$blue)) / ((2 * s2data$green) + s2data$red + s2data$blue)
tifGli                   # View raster structure
plot(tifGli)             # plot the raster      

# Calculates the Green Normalized Difference Vegetation Index (GNDVI) GNDVI = (NIR – Green)/(NIR + Green)    
tifGndvi <- (s2data$nir - s2data$green)/(s2data$nir + s2data$green)
tifGndvi                   # View raster structure
plot(tifGndvi)             # plot the raster   

# Calculates the Modified Non-Linear Index (MNLI) MNLI = (NIR^2 – Red)*(1+L)/(NIR^2 + Red + L)
L = 0.5
tifMnli <- ((s2data$nir^2 - s2data$red)*(1 + L))/(s2data$nir^2 + s2data$red + L)
tifMnli                   # View raster structure
plot(tifMnli)             # plot the raster   

# Calculate Modified Soil Adjusted Vegetation Index (MSAVI) MSAVI = (2 * NIR + 1 – sqrt ((2 * NIR + 1)2 – 8 * (NIR - R))) / 2
tifMsavi <- 0.5 * (2 * s2data$nir + 1 - sqrt((2 * s2data$nir + 1)^2 - 8 * (s2data$nir - s2data$red))) 
tifMsavi                   # View raster structure
plot(tifMsavi)             # plot the raster      

## Calculate Nonlinear Vegetation Index (NLI) NDVI = = (NIR^2 – Red)/(NIR^2 + Red)
tifNli <- (s2data$nir^2 - s2data$red)/(s2data$nir^2 + s2data$red)
tifNli                   # View raster structure
plot(tifNli)             # plot the raster   

## Calculate Normalized Difference Vegetation Index (NDVI) NDVI = (NIR – Red)/(NIR + Red)
tifNdvi <- (s2data$nir - s2data$red)/(s2data$nir + s2data$red)
tifNdvi                   # View raster structure
plot(tifNdvi)             # plot the raster   

# Calculate Modified Normalized Differential Greenness Index (NDGI) mNDGI = (Green^2 – Red)/(Green^2 + Red)
tifmNdgi <- (s2data$green * 2 - s2data$red)/(s2data$green * 2 + s2data$red) 
tifmNdgi                   # View raster structure
plot(tifmNdgi)             # plot the raster   

# Calculate Optimized Soil-Adjusted Vegetation Index (OSAVI) SAVI = ((NIR - Red) / (NIR + Red + 0.16)) 
tifOsavi <- (s2data$nir - s2data$red)/(s2data$nir + s2data$red + 0.16)
tifOsavi                   # View raster structure
plot(tifOsavi)             # plot the raster    

# Calculates the Renormalized Difference Vegetation Index (RDVI) RDVI = (NIR – Red)/sqrt(NIR + Red)
tifRdvi <- (s2data$nir - s2data$red)/sqrt(s2data$nir + s2data$red)
tifRdvi                   # View raster structure
plot(tifRdvi)             # plot the raster   

# Calculate the Soil-Adjusted Vegetation Index (SAVI) SAVI = ((NIR - Red) / (NIR + Red + L)) x (1 + L)
L = 0.5 
tifSavi <- (s2data$nir - s2data$red)/(s2data$nir + s2data$red + L)*(1 + L)
tifSavi                   # View raster structure
plot(tifSavi)             # plot the raster    

# Calculates the Transformed Difference Vegetation Index (TDVI) TDVI = 1.5 [(Nir - red) / sqrt(NIR^2 + Red + 0.5)]
tifTdvi <- 1.5 * ((s2data$nir - s2data$red)/sqrt(s2data$nir^2 + s2data$red + 0.5))
tifTdvi                   # View raster structure
plot(tifTdvi)             # plot the raster   

# Calculates the Visible Atmospherically Resistant Index (VARI)  VARI = (Green - Red) / (Green + Red – Blue)
tifVari = (s2data$green - s2data$red) / (s2data$green + s2data$red - s2data$blue)
tifVari                  # View raster structure
plot(tifVari)             # plot the raster      

# Calculates the Visible Atmospherically Resistant Vegetation Index (VARVI)  VARVI = (Green - Red) / (Green + NIR + Red)
tifVarvi = (s2data$green - s2data$red) / (s2data$green + s2data$nir + s2data$red)
tifVarvi                   # View raster structure
plot(tifVarvi)             # plot the raster      

# Calculates the Atmospherically Resistant Vegetation Index (ARVI) ARVI = (NIR – (2 * Red) + Blue) / (NIR + (2 * Red) + Blue)
tifArvi = (s2data$nir - (2 * s2data$red) + s2data$blue) / (s2data$nir + (2 * s2data$red) + s2data$blue)
tifArvi                   # View raster structure
plot(tifArvi)             # plot the raster      

# Calculates the Wide Dynamic Range Vegetation Index (WDRVI) WDRVI = (a * NIR – Red)/(a *NIR + Red)
a = 0.2
tifWdrvi = (a * s2data$nir - s2data$red)/(a * s2data$nir + s2data$red)
tifWdrvi                   # View raster structure
plot(tifWdrvi)             # plot the raster      

# Calculates the Green Ratio Vegetation Index (GRVI) GRVI = NIR / Green
tifGrvi = s2data$nir / s2data$green
tifGrvi                   # View raster structure
plot(tifGrvi)             # plot the raster      

# Calculates the Simple Ratio (SR) SR = NIR / Red
tifSr = s2data$nir / s2data$red
tifSr                   # View raster structure
plot(tifSr)             # plot the raster      

# Calculates the Blue Green Index (BGI)) BGI = Blue / Green
tifBgi = s2data$blue / s2data$green
tifBgi                   # View raster structure
plot(tifBgi)             # plot the raster      

# Calculates the Excess Green Index (ExGI) ExGI = (2 * Green) - (Red + Blue)
tifExgi = (2 * s2data$green) - (s2data$red + s2data$blue)
tifExgi                   # View raster structure
plot(tifExgi)             # plot the raster      

# Calculates the Visible-Band Difference Vegetation Index (VDVI) VDVI = (2 * Green - Red - Blue) / (2 * Green + Red + Blue) 
tifVdvi = ((2 * s2data$green) - s2data$red - s2data$blue) / ((2 * s2data$green) + s2data$red + s2data$blue)
tifVdvi                   # View raster structure
plot(tifVdvi)             # plot the raster      

# Calculates the Red-Green-Blue Vegetation Index (RGBVI) RGBVI = (Green^2 - Red * Blue) / (Green^2 + Red * Blue)
tifRgbvi = (s2data$green^2 - (s2data$red * s2data$blue))/(s2data$green^2 + (s2data$red * s2data$blue))
tifRgbvi                   # View raster structure
plot(tifRgbvi)             # plot the raster      

# Calculates the Green Chlorophyll Index (GCI) GCI = (NIR / Green) - 1
tifGci = (s2data$nir / s2data$green) - 1
tifGci                   # View raster structure
plot(tifGci)             # plot the raster      

# Calculates the Structure Insensitive Pigment Index (SIPI) GCI = (NIR / Green) - 1
tifSipi = (s2data$nir - s2data$blue) / (s2data$nir - s2data$red) 
tifSipi                   # View raster structure
plot(tifSipi)             # plot the raster      

# Calculates the Red-Green-Blue Vegetation Index (NGRDI) NGRDI = (Green - Red) / (Green + Red)
tifNgrdvi = (s2data$green - s2data$red)/(s2data$green + s2data$red)
tifNgrdvi                   # View raster structure
plot(tifNgrdvi)             # plot the raster      

# Calculates the Modified excess green (MExG) Mexg = 1.62 * green - 0.884 * red - 0.311 * blue
tifMexg <- 1.62 * s2data$green - 0.884 * s2data$red - 0.311 * s2data$blue  
tifMexg                   # View raster structure
plot(tifMexg)             # plot the raster      

# Calculate Modified Green Red Vegetation Index (MGVRI) MGVRI = (Green^2 – Red^2)/(Green^2 + Red^2)
tifMgvri <- (s2data$green ^ 2 - s2data$red ^ 2)/(s2data$green ^ 2 + s2data$red ^ 2) 
tifMgvri                   # View raster structure
plot(tifMgvri)             # plot the raster   

# Calculate Color Index of Vegetation (CIVE) CIVE = 0.441 Red − 0.881 Green + 0.385 Blue + 18.78745
tifCive <- 0.441 * s2data$red - 0.881 * s2data$green + 0.385 * s2data$blue + 18.78745
tifCive                   # View raster structure
plot(tifCive)             # plot the raster   
summary(tifCive)

# Calculate Normalized Difference Water Index (NDWI) NDWI = (Green - NIR) / (Green + NIR)
tifNdwi <- (s2data$green - s2data$nir) / (s2data$green + s2data$nir)
tifNdwi                   # View raster structure
plot(tifNdwi)             # plot the raster 

#### calculate textures ####
## Calculate Textures ----
# Grey-Level Co-Occurrence Matrices (GLCMs)
# Calculation rotation-invariant texture features in all 4 directions (0°, 45°, 90° and 135°) and then combined to one rotation-invariant texture. 
glcmNir <- glcm(x = s2data$nir, 
                window = c(5,5), 
                shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)), 
                statistics = c("mean","variance","correlation"), 
                na_opt = "ignore") # ,na_val = 999

names(glcmNir) <- c('txt_mean','txt_variance','txt_correlation')
glcmNir
plot(glcmNir) 

glcmNir[is.infinite(glcmNir)] <- NA 

xbar <- cellStats(glcmNir$txt_mean, stat = mean, na.rm = TRUE) # max, min, mean
glcmNir$txt_mean[is.na(glcmNir$txt_mean[]) | is.nan(glcmNir$txt_mean[]) | is.infinite(glcmNir$txt_mean[])] = xbar
xbar <- cellStats(glcmNir$txt_variance, stat = mean, na.rm = TRUE) # max, min, mean
glcmNir$txt_variance[is.na(glcmNir$txt_variance[]) | is.nan(glcmNir$txt_variance[]) | is.infinite(glcmNir$txt_variance[])] = xbar
xbar <- cellStats(glcmNir$txt_correlation, stat = mean, na.rm = TRUE) # max, min, mean
glcmNir$txt_correlation[is.na(glcmNir$txt_correlation[]) | is.nan(glcmNir$txt_correlation[]) | is.infinite(glcmNir$txt_correlation[])] = xbar
glcmNir

#### create feature stack ####
# Create a raster stack with the original layers and the calculated indexes  
b2data <- stack(s2data,
                tifDvi,tifEvi,tifGari,tifGli,tifGndvi,tifMnli,tifMsavi,tifNli,tifNdvi,tifmNdgi,tifOsavi,tifRdvi,tifSavi,
                tifTdvi,tifVari,tifVarvi,tifArvi,tifWdrvi,tifGrvi,tifSr,tifBgi,tifExgi,tifVdvi,tifRgbvi,tifGci,#tifSipi,
                tifNgrdvi,tifMexg,tifMgvri,tifCive,
                tifDem3m, tifSlope3m, tifAspect3m,
                tifSol3m,
                glcmNir) 
# rename
names(b2data) <- c('blue','green','red','nir',
                   'Dvi','Evi','Gari','Gli','Gndvi','Mnli','Msavi','Nli','Ndvi','mNdgi','Osavi','Rdvi','Savi',
                   'Tdvi','Vari','Varvi','Arvi','Wdrvi','Grvi','Sr','Bgi','Exgi','Vdvi','Rgbvi','Gci',#'Sipi',
                   'Ngrdvi','Mexg','Mgvri','Cive',
                   'dem','slope','aspect',
                   'silt','sand','clay','shc','ph','om','bd','wp','fc',
                   'Tmean','Tvar','Tcorr') 

# drop unneeded layers based on removal criteria
b2data <- dropLayer(b2data, sort(c(7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,26,27,30,31,32,33,36,38,45)))

# normalize the raster
# calculate min value per layer
min_val <- raster::cellStats(b2data, stat = 'min', na.rm = TRUE) 

# calculate max value per layer
max_val <- raster::cellStats(b2data, stat = 'max', na.rm = TRUE) 

# calculate range value per layer
range_val <- max_val - min_val 

# min-max normalize 
n2data <- scale(b2data, center = min_val, scale = range_val) 
n2data

#### load points for training/testing ####
# define the shapefile name and directory
shpName <- "20241118_mr_hs_kl_pnt.shp"
# shpName <- "augmented_field_data_6.shp" # if using the augmented dataset
shpPath <- file.path(WorkDir, shpName)

# load the point shapefile
shpPoints <- st_read(shpPath)

# if using the augmented field data file, run this line to fix the column name
# names(shpPoints)[names(shpPoints) == "class"] <- "Class"

# add an 'id' column to shpPoints
shpPoints <- shpPoints %>%
  dplyr::mutate(id = seq_along(Class))

# select relevant columns, rename for consistency
dfrPoints <- shpPoints %>%
  st_drop_geometry() %>%
  dplyr::select(id, Class)

# # If using ground truth data, this is where it can be integrated for model development
# groundTruthPath <- file.path(ShapeDir, "HS_field_data_UTM_2024-06-0307.shp")
# ground_truth_data <- st_read(groundTruthPath)
# 
# # match CRS to another dataset (e.g., tifDem3m)
# crs_s1data <- st_crs(s1data)
# ground_truth_data <- st_transform(ground_truth_data, crs_s1data)
# 
# # split the 'Comment' column at the comma and keep only the first part (the class)
# # i had originally inputted as 'class, description' and we want just class in this instance
# ground_truth_data$Comment <- sapply(strsplit(as.character(ground_truth_data$Comment), ","), `[`, 1)
# 
# # rename the 'Comment' column to 'Class' for consistency
# names(ground_truth_data)[names(ground_truth_data) == "Comment"] <- "Class"
# 
# # Split the 'Comment' column and extract the class, then rename
# ground_truth_data <- ground_truth_data %>%
#   #dplyr::mutate(class = sapply(strsplit(as.character(Comment), ","), `[`, 1)) %>%
#   dplyr::filter(Class %in% 1:6) %>% # keep only valid classes (1-6)
#   dplyr::mutate(id = seq_along(Class))%>%
#   #dplyr::select(-Comment) %>% # drop the original 'Comment' column
#   na.omit()
# 
# # reassign class 6 to join class 1, as class definitions changed post-fieldwork
# ground_truth_data <- ground_truth_data %>%
#   dplyr::filter(Class %in% 1:5) %>%  # keep only valid classes (1-5)
#   dplyr::mutate(
#     Class = ifelse(Class == 6, 1, Class),  # reassign Class 6 to 1
#     id = seq_along(Class)  # Create unique IDs
#   ) %>%
#   na.omit()

# # keep only relevant columns, as before
# gfrPoints <- ground_truth_data %>%
#   st_drop_geometry() %>%
#   dplyr::select(id, Class)
# 
# # extract band values to points
# samples <- raster::extract(x = n2data,    
#                            y = ground_truth_data, 
#                            method = 'simple',
#                            df = TRUE) 
# 
# # your samples layer must have a column for each image in the raster stack, a column for the land cover class that point represents, an X and Y column
# samples <- samples %>% 
#   dplyr::left_join(gfrPoints, by = c('ID' = "id")) %>%  
#   drop_na()
# head(samples); tail(samples)
# names(samples)
# samplesbkp <- samples # make a back up

# extract band values to points
samples <- raster::extract(x = n2data,    
                           y = shpPoints, 
                           method = 'simple',
                           df = TRUE)  

# your samples layer must have a column for each image in the raster stack, a column for the land cover class that point represents, an X and Y column
samples <- samples %>% 
  dplyr::left_join(dfrPoints, by = c('ID' = "id")) %>%  
  drop_na()
head(samples); tail(samples)
names(samples)
samplesbkp <- samples # make a back up

# zero and near-zero variance feature variables
nearZeroVar(samples, saveMetrics = TRUE)

# remove correlated  variables
# run correlation again if needed - skip this step the first time to find out which layers are highly correlated
#samples <- samplesbkp %>% 
#  dplyr::select(-'Tdvi', -'Savi', -'Sr', -'Msavi', -'Wdrvi', -'Osavi', -'Vdvi', -'Gli', -'Cive', -'Exgi',
#                -'Arvi', -'Mnli', -'Nli', -'Rdvi', -'Ngrdvi', -'Gari', -'Grvi', -'fc', -'bd', -'Vari',
#                -'Mgvri', -'mNdgi', -'Varvi', -'Mexg', -'Gndvi', -'aspect', -'sand')
names(samplesbkp); names(samples)
dplyr::glimpse(samples)

## exploratory data analysis      
dplyr::select(samples, -`ID`,-`Class`) %>% 
  skimr::skim() 

# if creating sample subsets, use this chunk first before moving on to creating training/testing split
# sampleIndex <- createDataPartition(samples$Class,
#                                    times = 1,     # times specifies how many splits to perform
#                                    p = 0.0379,       # p designates the split - 70/30
#                                    list = FALSE)  # n = 3413
# length(sampleIndex)
# sampleSet <- samples[sampleIndex,]
# prop.table(table(sampleSet$class))
# samples <- sampleSet

# split the data frame into 70-30 by class: training set and test set 
# partition based on the proportion from the response variable
set.seed(380843)
trainIndex <- createDataPartition(samples$Class, 
                                  times = 1,     # times specifies how many splits to perform
                                  p = 0.7,       # p designates the split - 70/30
                                  list = FALSE)  # n = 3413
length(trainIndex)
trainingSet <- samples[trainIndex,]
testSet <- samples[-trainIndex,]
head(trainIndex, n = 10); head(trainingSet, n = 5); head(testSet, n = 5)
names(samples); names(trainingSet); names(testSet)

# correlated feature variables
base_cor <-  cor(trainingSet[,(2:(ncol(trainingSet)-1))],
                 use = "complete.obs",  # "everything" "all.obs" "complete.obs" "na.or.complete" "pairwise.complete.obs".
                 method = "spearman")    # "pearson", "kendall", "spearman"

# save correlation matix
basecor <- abs(base_cor)
basecor[basecor < 0.9] = NA
write.table(basecor, file = paste(baseName, "baseCorrelation.csv", sep ="_"), row.names = TRUE, col.names = TRUE, sep = ',')

# check for any extreme correlations close to 1 in magnitude
extreme_cor <- sum(abs(base_cor[upper.tri(base_cor)]) > .85)
extreme_cor

# we assess a summary of the correlation values:
summary(base_cor[upper.tri(base_cor)])

#### create rf model and accuracy assessment ####
# set up a resampling method in the model training process
tc <- trainControl(method = "repeatedcv",            # repeated cross-validation of the training data: "cv"
                   number = 20,                      # number of folds
                   repeats = 15,                     # number of repeats
                   summaryFunction = defaultSummary, # function to compute performance metrics (twoClassSummary | multiClassSummary)
                   search = "grid",                  # "random" or grid" uses the default grid search routine
                   allowParallel = TRUE,             # allow use of multiple cores if specified in training
                   savePredictions = "final",        # "final" saves the predictions for the optimal tuning parameters.
                   verboseIter = TRUE)               # view the training iterations

# generate a grid search of candidate hyper-parameter values for inclusion into the models training
rf.grid <- expand.grid(.mtry = c(1:(ncol(trainingSet) - 2))) 

# doSNOW package to enable caret to train in parallel. 
cl <- makeCluster(32, type = "SOCK")
# register cluster 
registerDoSNOW(cl)

# begin training the models
rf_model <- caret::train(x = trainingSet[,(2:(ncol(trainingSet)-1))],              # OR: x = trainingSet[,(2:(ncol(trainingSet)-1))], # 'x =' as a data frame of our input columns | as.factor(`class`) ~ . -`ID`
                         y = as.factor(as.integer(as.factor(trainingSet$Class))),  # OR: y = as.factor(as.integer(as.factor(trainingSet$class))), # 'y =' as a vector of our target column | data = trainingSet
                         method = "rf",                                            # "ranger"
                         metric = "Accuracy",                                      # "ROC" | "Accuracy"
                         trainControl = tc,                                    
                         ntree = 600,                                              # Number of trees to grow. the default value is 500. A higher number of trees give you better performance
                         tuneGrid = rf.grid,
                         na.action = na.omit,                                      # na.pass
                         allowParallel = TRUE)                                     # importance = TRUE, Time difference: 1.73382 mins

# shutdown cluster
stopCluster(cl)

# figure and shows the relationship between the resampled performance values and the number of PLS components
plot(rf_model, main = "The best mtry numbers on model’s accuracy produced by grid search") 
getTrainPerf(rf_model)
rf_model 
rf_model$finalModel
rf_model$results

# get the confusion matrix
rf_mmtx <- rf_model$finalModel$confusion 
rf_mmtx <- rf_mmtx[,1:5]

# calculate the user’s accuracies 
usersacc <- diag(rf_mmtx) / rowSums(rf_mmtx) * 100
print(usersacc)

# calculate the producer’s accuracies
prodsacc <- diag(rf_mmtx) / colSums(rf_mmtx) * 100
print(prodsacc)

# calculate the overall accuracy 
overacc <- sum(diag(rf_mmtx)) / sum(rf_mmtx) * 100
print(overacc)

# mtry: Number of variables randomly selected as testing conditions at each split of decision trees. default value is sqr(col).
varimp_rf <- varImp(rf_model)

#### apply predictions to the test portion of the dataset ####
# make predictions on the test set
rf_predictions <- predict(rf_model, 
                          newdata = testSet)
head(rf_predictions)

# the option type = "prob" can be used to compute class probabilities from the model
rf_probs <- predict(rf_model, newdata = testSet, type = "prob")
head(rf_probs)

# apply the random forest model to the full raster
library(prettymapr)

rf_prediction = raster::predict(object = n2data, 
                                model = rf_model, 
                                fun = predict, 
                                na.rm = TRUE)

plot(rf_prediction,
     col = c('#4d4d4d','#6baed6','#ffffcc','#c2e699','#006837'), 
     main = "vegetation classification",
     addnortharrow(pos = "topright", scale = 0.5, pad = c(0.1, 0.1)))

# caret contains a function to compute the confusion matrix and associated statistics for the model fit:
confusionMatrix(data = rf_predictions, reference = as.factor(testSet$Class), mode = 'everything')

# Create a confusion matrix
rf_test_mtx <- table(Predicted = rf_predictions, Actual = testSet$Class) # Replace 'Class' with the actual column name

# calculate User's Accuracy
usersacc_test <- diag(rf_test_mtx) / rowSums(rf_test_mtx) * 100
print(usersacc_test)

# calculate Producer's Accuracy
prodsacc_test <- diag(rf_test_mtx) / colSums(rf_test_mtx) * 100
print(prodsacc_test)

# calculate Overall Accuracy
overacc_test <- sum(diag(rf_test_mtx)) / sum(rf_test_mtx) * 100
print(overacc_test)

# save model
saveRDS(rf_model, 
        file = paste(baseName, "model_rf.RDS", sep ="_"))        


