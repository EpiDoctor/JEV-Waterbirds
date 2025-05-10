library(sf)      
library(terra)
library(tidyverse)
library(tmap)
library(phytools)
library(FD)
library(ape)
source("trova.r")
source('melodic.r')
source("doublerda.R")
source("gawdis.R")

# First run script "C:/Users/thego/Dropbox/JEV/Data/manuscript.aus.com.jev.hmsc.data.prep.revised.21.02.2025.R" to set up species and trait data

#############################
# Functional traits analyses
#############################

# Waterbird species pool-weighted means (PWMs) for each site, first log transforming occurrence data:

birdCWM.1 <- functcomp(trait.data.jsdm, log(as.matrix(y.abund.jsdm) + 1))
head(birdCWM.1, 14)
tail(birdCWM.1, 14)

# Merge functional PWMs with sf spatial polygon df:

birdCWM.1$plot_id <- as.numeric(row.names(birdCWM.1))
head(birdCWM.1)
tail(birdCWM.1)
dim(birdCWM.1)

abundance <- merge(x = abundance, 
                          y = birdCWM.1, 
                          by.x = "ID", 
                          by.y = "plot_id", 
                          all.x = T)

dim(abundance)
head(abundance)
tail(abundance)

phylosig(consensus.tree, trait.data.jsdm$HWI, test = T, nsim = 1000, method = "K")
phylosig(consensus.tree, trait.data.jsdm$HWI, test = T, nsim = 1000, method = "lambda")

trait.hwi <- trait.data.jsdm$HWI
names(trait.hwi) <- row.names(trait.data.jsdm)
head(trait.hwi, 5)
head(trait.data.jsdm, 5)
tail(trait.hwi, 5)
tail(trait.data, 5)

hwi.cont.map <- contMap(consensus.tree, trait.hwi, type = "fan", fsize = .5, plot = F)

################################
# Functional diversity analyses
################################

# First compute dissimilarity matrix:
colnames(trait.data.jsdm)

# with population density in:
dissim.tr <- gawdis(trait.data.jsdm, w.type = "equal", groups = c(1, 1, 1, 2, 2, 2, 3, 4, 5, 6), fuzzy = TRUE)
dissim.tr

# with population density out:
dissim.tr <- gawdis(trait.data.jsdm, w.type = "equal", groups = c(1, 1, 1, 2, 2, 2, 3, 4, 5), fuzzy = TRUE)
dissim.tr

mpd.bird <- mpd(log(as.matrix(y.abund.jsdm) + 1), as.matrix(dissim.tr), abundance.weighted = T)

mpd.df <- data.frame(mpd.bird)

mpd.df$plot_id <- as.numeric(row.names(mpd.df))
head(mpd.df)
tail(mpd.df)
dim(mpd.df)

abundance <- merge(x = abundance, 
                   y = mpd.df, 
                   by.x = "ID", 
                   by.y = "plot_id", 
                   all.x = T)

head(abundance)
tail(abundance)
dim(abundance)

#############################################
# Convert functional trait fields to rasters
#############################################

colnames(abundance)

f.diet.inv <- terra::rasterize(vect(abundance), field = "Diet_invert", predictors$ardea.alba, touches = F)
writeRaster(f.diet.inv, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.diet.inv.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)
f.diet.fish <- terra::rasterize(vect(abundance), field = "Diet_fish", predictors$ardea.alba, touches = F)
writeRaster(f.diet.fish, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.diet.fish.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)
f.diet.plant <- terra::rasterize(vect(abundance), field = "Diet_plant", predictors$ardea.alba, touches = F)
writeRaster(f.diet.plant, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.diet.plant.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)

f.water.surface <- terra::rasterize(vect(abundance), field = "Water_surface", predictors$ardea.alba, touches = F)
writeRaster(f.water.surface, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.water.surface.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)
f.below.surface <- terra::rasterize(vect(abundance), field = "Below_surface", predictors$ardea.alba, touches = F)
writeRaster(f.below.surface, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.below.surface.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)
f.ground <- terra::rasterize(vect(abundance), field = "Ground", predictors$ardea.alba, touches = F)
writeRaster(f.ground, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.ground.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)

f.hwi <- terra::rasterize(vect(abundance), field = "HWI", predictors$ardea.alba, touches = F)
writeRaster(f.hwi, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.hwi.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)
f.body.mass.lg <- terra::rasterize(vect(abundance), field = "Body_mass_ln", predictors$ardea.alba, touches = F)
writeRaster(f.body.mass.lg, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.body.mass.lg.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)
f.egg.mass <- terra::rasterize(vect(abundance), field = "Egg_mass", predictors$ardea.alba, touches = F)
writeRaster(f.egg.mass, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.egg.mass.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)

f.mean.density <- terra::rasterize(vect(abundance), field = "Pop_density", predictors$ardea.alba, touches = F)
writeRaster(f.mean.density, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.div.mean.density.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)

mpd.bird.ras <- terra::rasterize(vect(abundance), field = "mpd.bird", predictors$ardea.alba, touches = F)
writeRaster(mpd.bird.ras, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.mpd.bird.ras.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)

# create new raster and tif for mpd computed without population density:
mpd.bird.ras.sans.density <- terra::rasterize(vect(abundance), field = "mpd.bird", predictors$ardea.alba, touches = F)
writeRaster(mpd.bird.ras.sans.density, filename="C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific host func div/SDM/5arcmin/func.mpd.bird.ras.sans.density.5.arcmin.ip.sdm.tif", filetype = "GTiff", overwrite = T)
