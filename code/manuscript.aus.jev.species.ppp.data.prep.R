library(sf)
library(terra)   
library(tidyverse)
library(tmap)
library(spatstat)
library(viridis)

# Australia shapefile
states <- read_sf("C:/Users/thego/Dropbox/JEV/Data/aust_cd66states.shp")
plot(states[0])
head(states)

states <- states %>%  
  st_set_crs("EPSG:4326")
states

east <- states %>% 
  filter(STE == 4 |
           STE == 3 |
           STE == 2 |
           STE == 1) 

east
plot(east[0])

###################
## JEV point data
###################

# First, read in WOAH data from followup Australia 2022 outbreak reports:

jev.1 <- read_csv("C:/Users/thego/Dropbox/JEV/Data/FUR_154421_05082022.csv")
dim(jev.1)
head(jev.1)
jev.1 <- jev.1 %>% 
  dplyr::select(outbreak_ID, 
                latitude, 
                longitude)
dim(jev.1)
head(jev.1)

jev.2 <- read_csv("C:/Users/thego/Dropbox/JEV/Data/FUR_154673_05082022.csv")
dim(jev.2)
head(jev.2)
jev.2 <- jev.2 %>% 
  dplyr::select(outbreak_ID, 
                latitude, 
                longitude)
dim(jev.2)
head(jev.2)

jev.3 <- read_csv("C:/Users/thego/Dropbox/JEV/Data/FUR_154770_05082022.csv")
dim(jev.3)
head(jev.3)
jev.3 <- jev.3 %>% 
  dplyr::select(outbreak_ID, 
                latitude, 
                longitude)
dim(jev.3)
head(jev.3)

jev.4 <- read_csv("C:/Users/thego/Dropbox/JEV/Data/FUR_154952_05082022.csv")
dim(jev.4)
head(jev.4)
jev.4 <- jev.4 %>% 
  dplyr::select(outbreak_ID, 
                latitude, 
                longitude)
dim(jev.4)
head(jev.4)

# Let's combine these into a single file - there will of course be overlaps, but we can easily remove the duplicates after:

jev <- bind_rows(jev.1, jev.2, jev.3, jev.4)
jev <- jev %>% 
  rename(eventID = outbreak_ID)

dim(jev)
head(jev)

# First, let's set up spatial points:

summary(jev$latitude)
summary(jev$longitude)

jev_sf <- jev %>% 
  st_as_sf(coords = c("longitude", "latitude")) %>% # set coordinates
  st_set_crs("EPSG:4326") # set geographic CRS

# Now we remove those duplicates:

jev_sf <- st_difference(jev_sf)
jev_sf

tm_shape(east) + 
  tmap_options(check.and.fix = TRUE) +
  tm_polygons() +
  tm_shape(jev_sf) + 
  tm_dots(col = "red")

###########################
# Independent testing data
###########################

# testing data derived from additional commercial piggery outbreaks, and NSW Health reported mosquito surveillance

jev.test <- read_csv("C:/Users/thego/Dropbox/JEV/Data/external_validation_data.csv")
dim(jev.test)
head(jev.test)

jev.test_sf <- jev.test %>% 
  st_as_sf(coords = c("longitude", "latitude")) %>% # set coordinates
  st_set_crs("EPSG:4326") # set geographic CRS

tm_shape(east) + 
  tmap_options(check.and.fix = TRUE) +
  tm_polygons() +
  tm_shape(jev_sf) + 
  tm_dots(col = "red") +
  tm_shape(jev.test_sf) + 
  tm_dots(col = "blue")

###########
## Rasters
###########

# APL pig and piggery densities:

pig.names <- c("pig.dens", 
               "mean.pig.dens", 
               "piggery.dens")
pig.list <- list.files(path = "C:/Users/thego/Dropbox/JEV/Data/APL Rasters", full.names = T)
pig.list

predictors.pigs <- rast(pig.list)

names(predictors.pigs) <- pig.names

predictors.pigs

# bird species habitat suitability rasters - estimated from SDMs:

bird.names <- c("amaurornis.isabellina", "amaurornis.olivacea", "amaurornis.phoenicurus", "anas.castanea", "anas.gibberifrons", 
                "anas.gracilis", "anas.platyrhynchos", "anas.superciliosa", "anastomus.oscitans", "anser.anser", 
                "anseranas.semipalmata", "ardea.alba", "ardea.cinerea", "ardea.pacifica", "ardea.purpurea", "ardea.sumatrana", 
                "ardeola.bacchus", "ardeola.speciosa", "aythya.australis", "aythya.fuligula", "biziura.lobata", "botaurus.poiciloptilus",
                "bubulcus.ibis", "butorides.striata", "cairina.moschata", "cereopsis.novaehollandiae", "chenonetta.jubata",
                "ciconia.stormi", "cygnus.atratus", "dendrocygna.arcuata", "dendrocygna.eytoni", "dendrocygna.guttata",
                "dendrocygna.javanica", "dupetor.flavicollis", "egretta.eulophotes", "egretta.garzetta", "egretta.intermedia", 
                "egretta.novaehollandiae", "egretta.picata", "egretta.sacra", "ephippiorhynchus.asiaticus", "fulica.atra", "gallicrex.cinerea",
                "gallinula.chloropus", "gallinula.mortierii", "gallinula.tenebrosa", "gallinula.ventralis", "gallirallus.castaneoventris",
                "gallirallus.philippensis", "gallirallus.striatus", "gallirallus.torquatus", "grus.antigone", "grus.rubicunda",
                "ixobrychus.cinnamomeus", "ixobrychus.dubius", "ixobrychus.sinensis", "leptoptilos.javanicus", "lewinia.pectoralis",
                "malacorhynchus.membranaceus", "microcarbo.melanoleucos", "microcarbo.niger", "mycteria.cinerea", "mycteria.leucocephala",
                "nettapus.coromandelianus", "nettapus.pulchellus", "nycticorax.caledonicus", "nycticorax.nycticorax", "oxyura.australis",
                "pelecanus.conspicillatus", "phalacrocorax.carbo", "phalacrocorax.fuscescens", "phalacrocorax.sulcirostris", 
                "phalacrocorax.varius", "platalea.flavipes", "platalea.regia", "plegadis.falcinellus", "porphyrio.indicus", "porphyrio.melanotus",
                "porzana.cinerea", "porzana.fluminea", "porzana.fusca", "porzana.pusilla", "porzana.tabuensis", "radjah.radjah", "rallina.tricolor",
                "spatula.clypeata", "spatula.querquedula", "spatula.rhynchotis", "stictonetta.naevosa", "tadorna.tadornoides", 
                "threskiornis.molucca", "threskiornis.spinicollis")

bird.names

bird.list <- list.files(path = "C:/Users/thego/Dropbox/JEV/Data/Indo-Pacific region sdms/abundance/5arcmin", full.names = T)
bird.list
predictors.birds <- rast(bird.list)
names(predictors.birds) <- bird.names
predictors.birds

# We'll crop and mask rasters as needed to the states under investigation:

predictors.birds <- crop(predictors.birds, east)
predictors.birds <- mask(predictors.birds, vect(east))

predictors.pigs <- mask(predictors.pigs, vect(east))

# Resample the following rasters to the bird rasters:

predictors.pigs <- resample(x = predictors.pigs, 
                     y = predictors.birds$ardea.alba, 
                     method = "bilinear")

# shift data for pig density so logarithm can be calculated for offset in ppms:

predictors.pigs <- predictors.pigs + 1

# full feature stack to assess all correlation
# Create raster stack to create analytic dataframe later, and to explore correlations among variables:

preds.names <- c("mean.pig.dens",
                 "piggery.dens",
                 "bird.rich")

predictors <- rast(list(predictors.pigs$mean.pig.dens,
                        predictors.pigs$piggery.dens,
                        bird.richness))

names(predictors) <- preds.names
predictors

############################
## Sample background points
############################

# JEV detection is expected to be greater in higher density piggeries
# So we sample background points proportional to mean pig density 

# use spatSample to sample proportional to that which is expected to bias reporting: pig density

set.seed(12)
bg.pts <- spatSample(predictors$mean.pig.dens, 
                     size = 10000, 
                     method = "weights", 
                     xy = T)

bg.pts = st_as_sf(bg.pts, coords = c("x", "y"), crs = 4326)
plot(predictors$mean.pig.dens)
points(bg.pts, 
       pch = 19, 
       cex = .1, 
       col = "red")


##############
## Spatstat
##############

# Spatstat requires it's own data format for points, so we can convert directly from sf:

# create point process from sf
jev.ppp <- as.ppp(
  X = st_coordinates(jev_sf), 
  W = as.owin(st_bbox(east)))

# to weight background points using points weighted by pig density we need to reset quadrature using quadscheme:

jev.bckg.ppp <- as.ppp(
  X = st_coordinates(bg.pts),
  W = as.owin(st_bbox(east))) 

jev.q <- quadscheme(data = jev.ppp, jev.bckg.ppp)

# finally, create point process for testing data:

dim(jev.test_sf)
jev.test.ppp <- as.ppp(
  X = st_coordinates(jev.test_sf), 
  W = as.owin(st_bbox(east)))

###########################
## Rasters to pixel images
###########################

# Finally, spatstat requires individual pixel images for each raster so we need to convert SpatRasters to pixel images:

# define function to convert:

as.im.SpatRaster1 <- function(X) {
  X <- X[[1]]
  rs <- terra::res(X)
  e <- as.vector(terra::ext(X))
  out <- list(
    v = as.matrix(X, wide=TRUE)[nrow(X):1, ],
    dim = dim(X)[1:2],
    xrange = e[1:2],
    yrange = e[3:4],
    xstep = rs[1],
    ystep = rs[2],
    xcol = e[1] + (1:ncol(X)) * rs[1] + 0.5 * rs[1],
    yrow = e[4] - (nrow(X):1) * rs[2] + 0.5 * rs[2],
    type = "real",
    units  = list(singular=units(X), plural=units(X), multiplier=1)
  )
  attr(out$units, "class") <- "unitname"
  attr(out, "class") <- "im"
  out
}

ras_list <- as.list(predictors.birds) # Convert SpatRaster stack to list
pixel_list <- list() # initialise the list of pixel images
pixel_list <- lapply(ras_list,  as.im.SpatRaster1) # apply pixel im conversion function
names(pixel_list) <- bird.names

piggery.dens.img <- as.im.SpatRaster1(predictors$piggery.dens)

