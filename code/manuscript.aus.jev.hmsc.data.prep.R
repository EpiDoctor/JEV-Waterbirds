library(sf)      
library(terra)
library(tidyverse)
library(tmap)
library(Hmsc)
library(ape)
library(phytools)
library(missForest)

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

# First, read in WOAH data from followup reports:

jev.1 <- read_csv("C:/Users/thego/Dropbox/JEV/Data/FUR_154421_05082022.csv")
dim(jev.1)
head(jev.1)
jev.1 <- jev.1 %>% 
  dplyr::select(outbreak_reference, 
                outbreak_ID, 
                start_date, 
                latitude, 
                longitude, 
                location,
                affected_population)
dim(jev.1)
head(jev.1)

jev.2 <- read_csv("C:/Users/thego/Dropbox/JEV/Data/FUR_154673_05082022.csv")
dim(jev.2)
head(jev.2)
jev.2 <- jev.2 %>% 
  dplyr::select(outbreak_reference, 
                outbreak_ID, 
                start_date, 
                latitude, 
                longitude, 
                location,
                affected_population)
dim(jev.2)
head(jev.2)

jev.3 <- read_csv("C:/Users/thego/Dropbox/JEV/Data/FUR_154770_05082022.csv")
dim(jev.3)
head(jev.3)
jev.3 <- jev.3 %>% 
  dplyr::select(outbreak_reference, 
                outbreak_ID, 
                start_date, 
                latitude, 
                longitude, 
                location,
                affected_population)
dim(jev.3)
head(jev.3)

jev.4 <- read_csv("C:/Users/thego/Dropbox/JEV/Data/FUR_154952_05082022.csv")
dim(jev.4)
head(jev.4)
jev.4 <- jev.4 %>% 
  dplyr::select(outbreak_reference, 
                outbreak_ID, 
                start_date, 
                latitude, 
                longitude, 
                location,
                affected_population)
dim(jev.4)
head(jev.4)

# Let's combine these into a single file - there will of course be overlaps, but we can easily remove the 
# duplicates after:

jev <- bind_rows(jev.1, jev.2, jev.3, jev.4)
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

###########
## Rasters
###########

# Distance to waterways:

waterway.dist <- rast("C:/Users/thego/Dropbox/Raster/Distance to major waterways/Distance to major waterways Australia/aus_osm_dst_waterway_100m_2016.tif")

# Land cover and landscape metrics based on IUCN Level 1 classes:

circle.grass <- rast("C:/Users/thego/Dropbox/JEV/Data/Aus landscape metrics/metric.class.circle.400.5.arcmin.aus.tif")
arable <- rast("C:/Users/thego/Dropbox/Raster/IUCN habitat classification data/Level 1 classes/iucn_habitatclassification_fraction_lvl2__1401_Arable Land__ver003.tif")

# Copernicus Global Surface Water Project - water seasonality - number of months water present per pixel in 2021

water.season <- rast("C:/Users/thego/Dropbox/Raster/Copernicus - Global Surface Water Project/Seasonality 2021/Australia/Processed/Presence/copernicus.water.seasonality.2021.30.seconds.tif")

# hydrological flow accumulation and flood hazard data:

flow.acc <- rast("C:/Users/thego/Dropbox/Raster/HydroSHEDS/Hydro flow accumulation/Australia/v1.4/hyd_au_acc_30s.tif")

# piggery data

pig.names <- c("pig.dens", "mean.pig.dens", "piggery.dens")
pig.list <- list.files(path = "C:/Users/thego/Dropbox/JEV/Data/APL Rasters", full.names = T)
pig.list

predictors.pigs <- rast(pig.list)

names(predictors.pigs) <- pig.names
predictors.pigs

# Only presence is recorded for land class area metrics for each land class in this raster stack so we need to 
# convert absence of each class to 0 so absence is represented: 
arable <- subst(arable, NA, 0)
circle.grass <- subst(circle.grass, NA, 0)

# We'll crop and mask these rasters to the states under investigation:
arable <- crop(arable, east)
flow.acc <- crop(flow.acc, east)
waterway.dist <- crop(waterway.dist, east)
predictors.pigs <- crop(predictors.pigs, east)

circle.grass <- mask(circle.grass, vect(east))
arable <- mask(arable, vect(east))
waterway.dist <- mask(waterway.dist, vect(east))
flow.acc <- mask(flow.acc, vect(east))
predictors.pigs <- mask(predictors.pigs, vect(east))

# Now we aggregate to 5 arc minutes:

waterway.dist <- aggregate(waterway.dist, 
                              fact = 100, 
                              fun = mean, 
                              na.rm = TRUE)

# Water present less than 25% of the year (seasonal/temp):
water.season.temp.presence.25 <- water.season
water.season.temp.presence.25[water.season.temp.presence.25 > 0 & water.season.temp.presence.25 < 4] <- 1
water.season.temp.presence.25[water.season.temp.presence.25 >= 4] <- 0

water.season.temp.presence.25 <- aggregate(water.season.temp.presence.25, 
                                           fact = 10, 
                                           fun = sum, 
                                           na.rm = TRUE)

flow.acc <- aggregate(flow.acc, 
                      fact = 10, 
                      fun = mean, 
                      na.rm = TRUE)

arable <- aggregate(arable,
                    fact = 10, 
                    fun = mean, 
                    na.rm = TRUE)


# bird species habitat suitability rasters - estimated from SDMs
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

predictors.birds <- crop(predictors.birds, east)
predictors.birds <- mask(predictors.birds, vect(east))

# resample to birds

arable <- resample(x = arable,
                   y = predictors.birds$ardea.alba, 
                   method = "bilinear")

circle.grass <- resample(x = circle.grass,
                         y = predictors.birds$ardea.alba, 
                         method = "bilinear")

flow.acc <- resample(x = flow.acc,
                     y = predictors.birds$ardea.alba, 
                     method = "bilinear")

water.season.temp.presence.25 <- resample(x = water.season.temp.presence.25,
                                          y = predictors.birds$ardea.alba, 
                                          method = "bilinear")

waterway.dist <- resample(x = waterway.dist,
                          y = predictors.birds$ardea.alba, 
                          method = "bilinear")

predictors.pigs <- resample(x = predictors.pigs, 
                            y = predictors.birds$ardea.alba, 
                            method = "bilinear")

# Create raster stack to create analytic dataframe later, and to explore correlations among variables:

predictors <- rast(list(predictors.birds,
                        arable, 
                        circle.grass, 
                        flow.acc, 
                        water.season.temp.presence.25, 
                        waterway.dist,
                        predictors.pigs$piggery.dens))

predictors

#######################
## Create grid polygon
#######################

# NOTE: the object, "abundance", generated below is NOT intending to indicate that this is a species abundance data object
# as would be derived from community data, which these data DO NOT comprise.

abundance <- as.polygons(predictors, dissolve = F, values = T, aggregate = F) |>
  st_as_sf()

abundance$ID <- 1:nrow(abundance)
head(abundance)

# Now we'll add JEV as an additional environmental variable:

abundance$jev_count <- lengths(st_intersects(abundance, jev_sf))
head(abundance)
summary(abundance$jev_count)

#############
# Trait data
#############

# Elton host trait data:

traits <- read_csv("BirdFuncDat.csv")
dim(traits)
head(traits)
colnames(traits)
colnames(traits) <- gsub("-", ".", colnames(traits))

traits <- traits %>%  
  dplyr::select(Scientific,
                Diet.Inv, 
                Diet.Vend, 
                Diet.Vect, 
                Diet.Vfish, 
                Diet.Vunk, 
                Diet.Scav, 
                Diet.Fruit, 
                Diet.Nect, 
                Diet.Seed, 
                Diet.PlantO, 
                Diet.5Cat, 
                ForStrat.watbelowsurf, 
                ForStrat.wataroundsurf, 
                ForStrat.ground, 
                ForStrat.SpecLevel,
                BodyMass.Value)
dim(traits)
head(traits)                                                
summary(traits)

# Now trim birds dataset to only waterbirds of Australia -  these include only species extant in Australia:

aus.species <- read_csv("Australia water host families for traits analyses.csv")
dim(aus.species)
colnames(aus.species)
head(aus.species, 15)
tail(aus.species, 15)

traits <- merge(x = traits, 
                y = aus.species, 
                by.x = "Scientific", 
                by.y = "species.name.resolved", 
                all = FALSE)
dim(traits)
head(traits)

# Bird density data from TetraDENSITY:

bird.density <- read_csv("TetraDENSITY.csv")
dim(bird.density)
head(bird.density)

bird.density$binomial <- paste(bird.density$Genus, 
                               bird.density$Species)

bird.density.avg <- bird.density %>%
  group_by(binomial) %>%
  summarize(mean_density = mean(Density, na.rm = TRUE)) %>% 
  dplyr::select(binomial, 
                mean_density)

traits <- merge(x = traits, 
                y = bird.density.avg, 
                by.x = "Scientific", 
                by.y = "binomial", 
                all.x = T)
dim(traits)
head(traits)

# AVONET

avonet <- read_csv("AVONET_Raw_Data.csv")
avonet
dim(avonet)
colnames(avonet)
colnames(avonet)[22] <- "HWI"

avonet.avg <- avonet %>%
  group_by(Species1_BirdLife) %>%
  summarize(mean_HWI = mean(HWI, na.rm = TRUE)) %>% 
  dplyr::select(Species1_BirdLife, 
                mean_HWI)

head(avonet.avg)
summary(avonet.avg$mean_HWI)

traits <- merge(x = traits, 
                y = avonet.avg, 
                by.x = "Scientific", 
                by.y = "Species1_BirdLife", 
                all.x = T)
dim(traits)
head(traits, 14)

# Lislevand bird traits (Ecological Archives E087-047):

life.hx.traits <- read_csv("avian_ssd_jan07.csv")
dim(life.hx.traits)
head(life.hx.traits)
life.hx.traits[life.hx.traits == "-999" | life.hx.traits == "-999.0" | life.hx.traits == "-999.00"] <- NA
head(life.hx.traits)
print(life.hx.traits$Species_name)
colnames(life.hx.traits)
summary(life.hx.traits$Clutch_size)
summary(life.hx.traits$Egg_mass)
summary(life.hx.traits$Mating_System)

life.hx.traits <- life.hx.traits %>%
  dplyr::select(Species_name, 
                Clutch_size, 
                Egg_mass)

traits <- merge(x = traits, 
                y = life.hx.traits, 
                by.x = "Scientific", 
                by.y = "Species_name", 
                all.x = T)
dim(traits)
head(traits)
print(traits$Scientific)
traits <- traits[-12,]
dim(traits)
head(traits, 14)
traits$Scientific
summary(traits)
summary(traits$mean_density)
summary(traits$Egg_mass)
summary(traits$mean_HWI)

##########################
# Missing data imputation
##########################

# Random forest imputation:

colnames(traits)

traits$Diet.Inv <- as.numeric(traits$Diet.Inv)
traits$Diet.Vend <- as.numeric(traits$Diet.Vend)
traits$Diet.Vect <- as.numeric(traits$Diet.Vect)
traits$Diet.Vfish <- as.numeric(traits$Diet.Vfish)
traits$Diet.Vunk <- as.numeric(traits$Diet.Vunk)
traits$Diet.Scav <- as.numeric(traits$Diet.Scav)
traits$Diet.5Cat <- as.factor(traits$Diet.5Cat)
traits$ForStrat.watbelowsurf <- as.numeric(traits$ForStrat.watbelowsurf)
traits$ForStrat.wataroundsurf <- as.numeric(traits$ForStrat.wataroundsurf)
traits$ForStrat.ground <- as.numeric(traits$ForStrat.ground)
traits$wos <- as.numeric(traits$wos)
traits$BodyMass.Value.lg <- log(traits$BodyMass.Value)

# Imputation using missForest
colnames(traits)
# imputation done excluding first column of species names:
set.seed(10)
traits.imp <- missForest(traits[,c(2:19,21:25)])

# now we need to add species names back to imputed traits dataframe:
head(traits.imp$ximp)
head(traits)
tail(traits.imp$ximp)
tail(traits)

traits.imp$ximp$Scientific <- traits$Scientific
head(traits.imp$ximp)
head(traits)
tail(traits.imp$ximp)
tail(traits)

# Now let's compare distribution of imputed traits to traits with missing data:
dim(traits.imp$ximp)
summary(traits$mean_HWI)
summary(traits.imp$ximp$mean_HWI)
summary(traits$Egg_mass)
summary(traits.imp$ximp$Egg_mass)
summary(traits$mean_density)
summary(traits.imp$ximp$mean_density)

plot(density(traits.imp$ximp$mean_HWI), 
     lty = 3, lwd = 2, col = "red", 
     main = "Hand-wing index", 
     xlab = "(Wing length-Secondary length/Wing length)*100")
lines(density(traits$mean_HWI, na.rm = T), lwd = 2)
plot(density(traits.imp$ximp$Egg_mass), 
     lty = 3, lwd = 2, col = "red", 
     main = "Egg mass", 
     xlab = "grams")
lines(density(traits$Egg_mass, na.rm = T), lwd = 2)
plot(density(traits.imp$ximp$mean_density), 
     lty = 3, lwd = 2, col = "red", 
     main = "Population density", 
     xlab = "n")
lines(density(traits$mean_density, na.rm = T), lwd = 2)

############
# Phylogeny
############

tree <- read.nexus("C:/Users/thego/Dropbox/JEV/Data/Bird reservoir data/VertLife tree/Australia all waterbird species/Australia waterbird phylo tree with Hackett backbone.nex")
tree

set.seed(12)
random.tree <- sample(tree,size=1)[[1]]
random.tree

set.seed(12)
consensus.tree <- ls.consensus(tree, 
                               start = random.tree)

consensus.tree

plot(consensus.tree, 
     cex = .75,
     show.tip.label = T)

###################
# HMSC data set up
###################

colnames(abundance)
aus.birds <- c("anas.castanea", "anas.gracilis", "anas.platyrhynchos", "anas.superciliosa", "anser.anser", 
               "anseranas.semipalmata", "ardea.alba", "ardea.pacifica", "ardea.sumatrana", 
               "aythya.australis", "aythya.fuligula", "biziura.lobata", "botaurus.poiciloptilus",
               "bubulcus.ibis", "butorides.striata", "cairina.moschata", "cereopsis.novaehollandiae", "chenonetta.jubata",
               "cygnus.atratus", "dendrocygna.arcuata", "dendrocygna.eytoni", "dendrocygna.guttata",
               "dupetor.flavicollis", "egretta.garzetta", "egretta.intermedia", 
               "egretta.novaehollandiae", "egretta.picata", "egretta.sacra", "ephippiorhynchus.asiaticus", "fulica.atra",
               "gallinula.mortierii", "gallinula.tenebrosa", "gallinula.ventralis", "gallirallus.castaneoventris",
               "gallirallus.philippensis", "grus.antigone", "grus.rubicunda",
               "ixobrychus.dubius", "lewinia.pectoralis", "malacorhynchus.membranaceus", "microcarbo.melanoleucos",
               "nettapus.coromandelianus", "nettapus.pulchellus", "nycticorax.caledonicus", "oxyura.australis",
               "pelecanus.conspicillatus", "phalacrocorax.carbo", "phalacrocorax.fuscescens", "phalacrocorax.sulcirostris", 
               "phalacrocorax.varius", "platalea.flavipes", "platalea.regia", "plegadis.falcinellus", "porphyrio.melanotus",
               "porzana.cinerea", "porzana.fluminea", "porzana.pusilla", "porzana.tabuensis", "radjah.radjah", "rallina.tricolor",
               "spatula.clypeata", "spatula.querquedula", "spatula.rhynchotis", "stictonetta.naevosa", "tadorna.tadornoides", 
               "threskiornis.molucca", "threskiornis.spinicollis")

y.abund.jsdm <- as.data.frame(abundance)[,c(1:92)]
y.abund.jsdm <- y.abund.jsdm %>%
  dplyr::select(aus.birds)

colnames(y.abund.jsdm) <- gsub(".", "_", fixed = T, colnames(y.abund.jsdm))
colnames(y.abund.jsdm) <- str_to_title(colnames(y.abund.jsdm))
colnames(y.abund.jsdm)
y.abund.jsdm <- y.abund.jsdm[ , order(names(y.abund.jsdm))]
colnames(y.abund.jsdm)

# resolve binomials with VertLife tree:

colnames(y.abund.jsdm)[7] <- "Casmerodius_albus"
colnames(y.abund.jsdm)[23] <- "Ixobrychus_flavicollis"
colnames(y.abund.jsdm)[25] <- "Mesophoyx_intermedia"
colnames(y.abund.jsdm)[27] <- "Ardea_picata"
colnames(y.abund.jsdm)[34] <- "Eulabeornis_castaneoventris"
colnames(y.abund.jsdm)[38] <- "Ixobrychus_minutus"
colnames(y.abund.jsdm)[41] <- "Phalacrocorax_melanoleucos"
colnames(y.abund.jsdm)[54] <- "Porphyrio_porphyrio"
colnames(y.abund.jsdm)[59] <- "Tadorna_radjah"
colnames(y.abund.jsdm)[61] <- "Anas_clypeata"
colnames(y.abund.jsdm)[62] <- "Anas_querquedula"
colnames(y.abund.jsdm)[63] <- "Anas_rhynchotis"

y.abund.jsdm <- y.abund.jsdm[ , order(names(y.abund.jsdm))]
colnames(y.abund.jsdm)

colnames(abundance)
x.data.jsdm <- as.data.frame(abundance)[,c(93:98,101)] 
x.data.jsdm$jev_pres <- x.data.jsdm$jev_count
x.data.jsdm$jev_pres[x.data.jsdm$jev_pres > 0] <- 1
x.data.jsdm$jev_pres <- as.factor(x.data.jsdm$jev_pres)
colnames(x.data.jsdm)
x.data.jsdm <- x.data.jsdm[,c(-7)]
xy.latlong <- st_coordinates(st_centroid(abundance))
study.design <- data.frame(sample = as.factor(abundance$ID))
rownames(xy.latlong) = study.design[ ,1]
head(y.abund.jsdm)
env.names <- c("Croplands", 
               "RCC_grassland", 
               "Flow_accumulation", 
               "Transient_wetlands", 
               "Dist_to_waterways",
               "Piggeries",
               "JEV_present")

colnames(x.data.jsdm) <- env.names
head(x.data.jsdm)
head(xy.latlong)

# Trait data

colnames(traits.imp$ximp)
head(traits.imp$ximp, 14)
trait.data.jsdm <- traits.imp$ximp
trait.data.jsdm <- trait.data.jsdm %>%
  dplyr::select(Diet.Inv,
                Diet.Vfish,
                Diet.PlantO,
                ForStrat.watbelowsurf,
                BodyMass.Value.lg,
                wos,
                mean_HWI,
                mean_density,
                Egg_mass, 
                Scientific)
  
trait.data.jsdm$Scientific <- gsub(" ", "_", trait.data.jsdm$Scientific)
head(trait.data.jsdm)
dim(trait.data.jsdm)
rownames(trait.data.jsdm) <- trait.data.jsdm$Scientific
colnames(trait.data.jsdm)

trait.data.jsdm <- trait.data.jsdm[,c(1:9)]

# need to remove population density because imputation was suboptimal for this trait:
trait.data.jsdm <- trait.data.jsdm[,-8]

# As above with environmental data, let's rename vars for better visualisations:
colnames(trait.data.jsdm)

# with population density out:
trait.names <- c("Diet_inverts", 
                 "Diet_fish",
                 "Diet_plants",
                 "Below_water",
                 "Body_mass_ln",
                 "Studiedness",
                 "HWI", 
                 "Egg_mass")

colnames(trait.data.jsdm) <- trait.names
head(trait.data.jsdm)

# Finally we'll remove Gallinula mortierii since this species is limited to Tasmania and is not included in any of the 
# landscape units for this analysis:

colnames(y.abund.jsdm)
y.abund.jsdm <- y.abund.jsdm[,-33]
rownames(trait.data.jsdm)
trait.data.jsdm <- trait.data.jsdm[-33,]
consensus.tree <- drop.tip(consensus.tree, "Gallinula_mortierii")
