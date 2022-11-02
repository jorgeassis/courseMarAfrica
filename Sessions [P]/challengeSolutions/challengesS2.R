# ----------------------------------------------------
# ----------------------------------------------------
# Challenge 2.1
# The main objective of this challenge is to plot a shapefile with the distribution of a seagrass species. The plot must include the landmass shapefile of the world.

library(raster)

data1 <- "Data/vectorShapefiles/seagrass/cymodoceaNodosa.shp"
data2 <- "Data/vectorShapefiles/globalLandmass/world.shp"

?shapefile

records <- shapefile(data1)
world <- shapefile(data2)

plot(world)
plot(records)

plot(world)
plot(records, add=TRUE,col="red")

# ----------------------------------------------------
# ----------------------------------------------------
# Challenge 2.2
# The main objective of this challenge is to crop a global map to the Azores islands region and plot it.

library(raster)

data <- "Data/vectorShapefiles/globalLandmass/world.shp"

world <- shapefile(data)

# Crop with extent

regionAzores <- extent(-36, -18, 34, 42)

azores <- crop(world,regionAzores)
plot(azores)

# ----------------------------------------------------
# ----------------------------------------------------
# Challenge 2.3
# The main objective of this challenge is to make two plots of the global maximum and minimum SST.

library(raster)

SST <- raster("Data/rasterLayers/Climate/Present/OceanTemperature Surface Pred Mean.tif")

plot(SST, main ="Sea Surface Temperatures")

# ----------------------------------------------------
# ----------------------------------------------------
# Challenge 2.4
# The main objective of this challenge is to plot the estimated warming of the ocean.

library(raster)

presentSST <- raster("Data/rasterLayers/Climate/Present/OceanTemperature Surface Pred Mean.tif")
futureSST <- raster("Data/rasterLayers/Climate/RCP85/OceanTemperature Surface Pred Mean.tif")

differenceSST <- futureSST - presentSST
plot(differenceSST, main="Global warming" )

# ----------------------------------------------------
# ----------------------------------------------------
# Challenge 2.5
# The main objective of this challenge is to plot the estimated warming of the ocean in West Africa.

library(raster)

presentSST <- raster("Data/rasterLayers/Climate/Present/OceanTemperature Surface Pred Mean.tif")
futureSST <- raster("Data/rasterLayers/Climate/RCP85/OceanTemperature Surface Pred Mean.tif")

differenceSST <- futureSST - presentSST

westAfricaExtent <- extent(-30, 22, -37.5, 36)
differenceSSTWAfrica <- crop(differenceSST,westAfricaExtent)
plot(differenceSSTWAfrica, main="West Africa warming" )

# ----------------------------------------------------
# ----------------------------------------------------
# Challenge 2.6
# The main objective of this individual assignment is to **extract** the depth range used (values) by a mediterranean coral and plot it as an histogram.

data1 <- "Data/rasterLayers/BathymetryDepthMean.tif"
data2 <- "Data/dataBases/Paramuricea_clavata.csv"

bathymetry <- raster(data1)
occurrences <- read.csv(data2,sep=";")

colnames(occurrences)

depthsUsed <- extract(bathymetry,occurrences[,c("Lon","Lat")])
hist(depthsUsed,breaks=100)
