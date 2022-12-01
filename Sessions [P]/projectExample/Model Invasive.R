## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
## Recipe 2
# The main objetive of this Recipe is to fit a BRT **model** to biodiversity data predict distributions for the native and invasive regions

## -----------------------
## -----------------------

# Set the working directory to where the data is located
setwd(".../projectExample")

# Load main functions
source("sourceFunctions.R")

## -----------------------
# 01. open final records

# load clean occurrence data with two columns only for Lon and Lat (follow Recipe 1)
records <- read.csv("Data/mycleanRecords.csv", sep = ";")

## -----------------------
# 02. open and crop landmass

# Read polygon defining global landmasses
world <- shapefile("Data/globalLandmass/world.shp")

# Define my extent
myExtent <- extent(-20,47.5,10,60)

# Crop polygon to my region
myRegion <- crop(world,myExtent)
plot(myRegion,col="gray",border="gray")

# plot my records
points(records, pch=20, col="Black")

# Save image file to PNG or PDF [Plots panel -> Export]

## -----------------------
# 02. open and crop environmental layers

# Load layers
environmentPresent <- stackDirectory('Data/rastersPresent/')
plot(environmentPresent)

# mask to coastLine [if needed!!] and plot
maskCoastal <- raster("Data/CoastLine.tif")
plot(maskCoastal, col="black")

environmentPresent <- mask(environmentPresent,maskCoastal) 
plot(environmentPresent)

# crop and plot layers
environmentPresent <- crop(environmentPresent,myExtent)
plot(environmentPresent)

# Save image file to PNG or PDF [Plots panel -> Export]

## -----------------------
## -----------------------
# Inside function

myRasterLayers <- environmentPresent
records <- records

options(warn=-1)
spobj1 <- records
spobj1 <- spobj1[which(!is.na(spobj1[,1])),] 
spobj1 <- spobj1[which(!is.na(spobj1[,2])),] 
spobj2 <- subset(myRasterLayers,1)
spobj2 <- crop(spobj2,extent(c(min(spobj1[,1])-1,max(spobj1[,1])+1,min(spobj1[,2])-1,max(spobj1[,2])+1)))

toCorrect <- which(is.na(raster::extract(spobj2,spobj1[c(1,2)])))
overLand <- 0
dist <- 25

if( length(toCorrect) > 0) {
  
  spobj2Cells <- Which(!is.na(spobj2),cells=T)
  corrected <- xyFromCell(spobj2,spobj2Cells)
  
  for(i in 1:length(toCorrect)) {
    
    dists <- spDistsN1(as.matrix(corrected),as.matrix(spobj1[,c(1,2)])[toCorrect[i],],longlat = TRUE)
    
    if( min(dists) <= dist){ 
      closest <- which.min(spDistsN1(as.matrix(corrected),as.matrix(spobj1[,c(1,2)])[toCorrect[i],],longlat = TRUE))
      spobj1[toCorrect[i],1] <- corrected[closest,1]
      spobj1[toCorrect[i],2] <- corrected[closest,2]
    }
    
    if(min(dists) > dist){ 
      spobj1[toCorrect[i],1] <- NA
      spobj1[toCorrect[i],2] <- NA
      overLand <- overLand + 1
    }
  }
  
}

options(warn=0)

records <- spobj1[which(!is.na(spobj1[,1])),] 

# generate pseudo absences
pseudoAbs <- pseudoAbsences(myRasterLayers,records,n=1000)
plot(pseudoAbs)
points(records, col="red")

## -----------------------
# 04. fit a model with best hyperparameters

# extract environmental values and make a data.frame with PA information

p <- records
a <- pseudoAbs
m <- subset(myRasterLayers,1)

p.i <- raster::extract(m,p)
p <- p[which(!is.na(p.i)),]

p <- xyFromCell(m,cellFromXY(m,p))
p <- unique(p)

a.i <- raster::extract(m,a)
a <- a[which(!is.na(a.i)),]

a <- xyFromCell(m,cellFromXY(m,a))
a <- unique(a)

modelData <- prepareSWD(species = "Model species", p = p, a = a, env = myRasterLayers)

# Generate cross validation folds
folds <- getBlocks(modelData)

# fit a BRT model with cross-validation and 
model <- train("BRT", modelData, folds = folds)

# given a set of possible hyperparameter values for BRT
h <- list(interaction.depth = c(1,2,3,4) , shrinkage = c(0.01,0.001) )

# test all possible combinations of hyperparameter values
exp1 <- gridSearch(model, hypers = h, metric = "auc")

# fit a BRT model to the dataset with the best hyperparameter values
model <- train("BRT", modelData, folds = folds , interaction.depth=exp1@results[which.max(exp1@results$test_AUC),"interaction.depth"], shrinkage=exp1@results[which.max(exp1@results$test_AUC),"shrinkage"] )
aucModel <- getAUC(model, test = TRUE)
viModel <- varImp(model, permut = 5)
plotVarImp(viModel)

# inspect response curves
names(myRasterLayers)
plotResponse(model, var = "OceanTemperature.Benthic.Mean.Pred.LtMax", type = "logistic", only_presence = FALSE, marginal = FALSE, rug = FALSE, color="Black")
plotResponse(model, var = "DissolvedMolecularOxygen.Benthic.Mean.Pred.Mean", type = "logistic", only_presence = FALSE, marginal = FALSE, rug = FALSE, color="Black")
plotResponse(model, var = "OceanTemperature.Benthic.Mean.Pred.LtMin", type = "logistic", only_presence = FALSE, marginal = FALSE, rug = FALSE, color="Black")

## -----------------------
# 06. predict to produce maps

# predict with BRT to raster stack
mapNative <- predict(model, environmentPresent, type=c("logistic"))
plot(mapNative)

# Save image file to PNG or PDF [Plots panel -> Export]

## -----------------------
# 08. predict to an invasive region

# Load layers
environmentPresentInvasive <- stackDirectory('Data/rastersPresent/')
plot(environmentPresentInvasive)

# mask to coastLine [if needed!!] and plot
maskCoastal <- raster("Data/CoastLine.tif")
plot(maskCoastal, col="black")

environmentPresentInvasive <- mask(environmentPresentInvasive,maskCoastal)
plot(environmentPresentInvasive)

# crop and plot layers with a new extent
myNewExtent <- extent(-137,-100,14,46)
environmentPresentInvasive <- crop(environmentPresentInvasive,myNewExtent)
plot(environmentPresentInvasive)

# predict with BRT to raster stack
mapInvasive <- predict(model, environmentPresentInvasive, type=c("logistic"))
plot(mapInvasive)

# Save image file to PNG or PDF [Plots panel -> Export]

par(mfrow=c(1,3))
plot(mapPresent)
plot(mapRCP26 - mapPresent)
plot(mapRCP85 - mapPresent)

dev.off()

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------