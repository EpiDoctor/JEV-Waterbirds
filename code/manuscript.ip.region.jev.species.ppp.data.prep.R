library(sf)
library(terra)   
library(tidyverse)
library(tmap)
library(spatstat)
library(viridis)

# World shapefile
region <- read_sf("C:/Users/thego/Dropbox/JEV/Data/WB_countries_Admin0_10m.shp")
head(region)

# Select countries for Central Indo-Pacific region

region <- region %>% 
  filter(NAME_EN == "Indonesia" |
           NAME_EN == "East Timor" |
           NAME_EN == "Brunei" |
           NAME_EN == "Papua New Guinea"|
           NAME_EN == "Australia" |
           NAME_EN == "New Caledonia" |
           NAME_EN == "Philippines" |
           NAME_EN == "Solomon Islands"|
           NAME_EN == "Vanuatu" |
           NAME_EN == "Malaysia") 
region
plot(region[0])

# need to specify region without Australia for different pig density for background point sampling:

ip.pigs.region <- region %>% 
  filter(NAME_EN == "Indonesia" |
           NAME_EN == "East Timor" |
           NAME_EN == "Brunei" |
           NAME_EN == "Papua New Guinea"|
           NAME_EN == "New Caledonia" |
           NAME_EN == "Philippines" |
           NAME_EN == "Solomon Islands"|
           NAME_EN == "Vanuatu" |
           NAME_EN == "Malaysia") 
ip.pigs.region
plot(ip.pigs.region[0])

# Australia shapefile
states <- read_sf("C:/Users/thego/Dropbox/JEV/Data/aust_cd66states.shp")
plot(states[0])
head(states)

states <- states %>%  
  st_set_crs("EPSG:4326")
states

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

tm_shape(region) + 
  tmap_options(check.and.fix = TRUE) +
  tm_polygons() +
  tm_shape(jev_sf) + 
  tm_dots(col = "red")

# Now we add data from additional detection from around the central IP region:

jev.ip <- read_csv("C:/Users/thego/Dropbox/JEV/Data/IP region JEV detection.csv")
dim(jev.ip)
head(jev.ip)
jev.ip <- jev.ip %>% 
  dplyr::select(eventID, 
                latitude, 
                longitude)
dim(jev.ip)
head(jev.ip)

jev_ip <- jev.ip %>% 
  st_as_sf(coords = c("longitude", "latitude")) %>% # set coordinates
  st_set_crs("EPSG:4326") # set geographic CRS

dim(jev_ip)
head(jev_ip)

# Now we combine the Australian and IP region JEV data: 
jev_sf <- rbind(jev_sf, jev_ip)

dim(jev_sf)
head(jev_sf)

tm_shape(region) + 
  tmap_options(check.and.fix = TRUE) +
  tm_polygons() +
  tm_shape(jev_sf) + 
  tm_dots(col = "red")

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

ip.pigs <- rast("C:/Users/thego/Dropbox/JEV/Data/Glb_Pigs_CC2006_AD.tif")
ip.pigs

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

ip.pigs <- crop(ip.pigs, ip.pigs.region) 

ip.pigs <- mask(ip.pigs, vect(ip.pigs.region))
predictors.pigs <- mask(predictors.pigs, vect(states))

# Now we aggregate pigs to 5 arc minutes as needed:

ip.pigs <- aggregate(ip.pigs, 
                     fact = 10, 
                     fun = sum, 
                     na.rm = TRUE)

# We'll now merge the APL pig density for Australia with GLB pig density for IP region, 
# which we'll use for background point selection:

pig.dens <- sprc(predictors.pigs$mean.pig.dens,
                 ip.pigs)

pig.dens <- merge(pig.dens)

# Resample the following rasters to the bird rasters:

pig.dens <- resample(x = pig.dens, 
                     y = predictors.birds$amaurornis.isabellina, 
                     method = "bilinear")

############################
## Sample background points
############################

# JEV detection is expected to be greater in higher density piggeries
# So we sample background points proportional to mean pig density 

# use spatSample to sample proportional to that which is expected to bias reporting: pig density

set.seed(12)
bg.pts <- spatSample(pig.dens, 
                     size = 20000, 
                     method = "weights", 
                     xy = T)

bg.pts = st_as_sf(bg.pts, coords = c("x", "y"), crs = 4326)
plot(pig.dens)
points(bg.pts, 
       pch = 19, 
       cex = .1, 
       col = "blue")


##############
## Spatstat
##############

# Spatstat requires it's own data format for points, so we can convert directly from sf:

# create point process from sf
jev.ppp <- as.ppp(
  X = st_coordinates(jev_sf), 
  W = as.owin(st_bbox(region)))

# to weight background points using points weighted by pig density we need to reset quadrature using quadscheme:

jev.bckg.ppp <- as.ppp(
  X = st_coordinates(bg.pts),
  W = as.owin(st_bbox(region))) 

jev.q <- quadscheme(data = jev.ppp, jev.bckg.ppp)

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


