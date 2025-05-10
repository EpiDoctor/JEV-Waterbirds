library(sf)
library(terra)   
library(tidyverse)
library(tmap)
library(spatstat)
library(usdm)
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

###########################
## Kernel density estimate
###########################

jev.intens <- density.ppp(jev.ppp, 
                          sigma = 2, 
                          diggle = T)

jev.intens.ras <- rast(jev.intens)
crs(jev.intens.ras)  <- "EPSG:4326"

jev.intens.ras
jev.intens.ras <- mask(jev.intens.ras, vect(region))
plot(jev.intens.ras, 
     col = rev(viridis(10, option = "A")))

###########################
## Estimate the K-function
###########################

set.seed(4)
kf.jev.annual.env <- envelope(jev.ppp, Kest)
plot(kf.jev.annual.env)

#############
# Null model
#############

null <- ppm(jev.q, ~1) 
summary(null)

AIC(null)

##################
## Bivariate PPMs
##################

# Functional trait models with reduced background data capturing locations with bird communities only

jev.1 <- ppm(jev.q, ~ trait.diet.fish,, 
             covariates = list(trait.diet.fish = trait.diet.fish.img))

jev.1
AIC(jev.1)

jev.2 <- ppm(jev.q, ~ trait.diet.inv,
             covariates = list(trait.diet.inv = trait.diet.inv.img))

jev.2
AIC(jev.2)

jev.3 <- ppm(jev.q, ~ trait.diet.plant,
             covariates = list(trait.diet.plant = trait.diet.plant.img))
             
jev.3
AIC(jev.3)

jev.4 <- ppm(jev.q, ~ trait.body.mass.lg,
             covariates = list(trait.body.mass.lg = trait.body.mass.lg.img))

jev.4
AIC(jev.4)

jev.5 <- ppm(jev.q, ~ trait.egg.mass, 
             covariates = list(trait.egg.mass = trait.egg.mass.img))

jev.5
AIC(jev.5)

jev.6 <- ppm(jev.q, ~ trait.hwi,
             covariates = list(trait.hwi = trait.hwi.img))

jev.6
AIC(jev.6)

jev.7 <- ppm(jev.q, ~ trait.water.surface, 
             covariates = list(trait.water.surface = trait.water.surface.img))

jev.7
AIC(jev.7)

jev.8 <- ppm(jev.q, ~ trait.below.water.surface, 
             covariates = list(trait.below.water.surface = trait.below.water.surface.img))

jev.8
AIC(jev.8)

jev.9 <- ppm(jev.q, ~ trait.ground,
             covariates = list(trait.ground = trait.ground.img))

jev.9
AIC(jev.9)

# population density no longer included because of suboptimal imputation of this trait
jev.10 <- ppm(jev.q, ~ trait.mean.density,
              covariates = list(trait.mean.density = trait.mean.density.img))

jev.10
AIC(jev.10)

# MPD computed without population density because of suboptimal imputation of this trait
jev.11 <- ppm(jev.q, ~ trait.mpd.sans.dens, 
              covariates = list(trait.mpd.sans.dens = trait.mpd.sans.dens.img))
jev.11
AIC(jev.11)

jev.12 <- ppm(jev.q, ~ trait.mpd.sans.dens + I(trait.mpd.sans.dens^2), 
              covariates = list(trait.mpd.sans.dens = trait.mpd.sans.dens.img))
jev.12
AIC(jev.12)

jev.13 <- ppm(jev.q, ~ bird.rich, 
              covariates = list(bird.rich = bird.rich.img))
jev.13
AIC(jev.13)

jev.14 <- ppm(jev.q, ~ bird.rich + I(bird.rich^2), 
              covariates = list(bird.rich = bird.rich.img))
jev.14
AIC(jev.14)

# Let's look at the correlation and VIF for all the features that were bivariately associated with JEV occurrences so we know where collinearity could be a
# problem:

predictors.multi.ppms.func <- rast(list(predictors$trait.hwi,
                                        predictors$trait.for.below.surface,
                                        predictors$trait.egg.mass,
                                        predictors$trait.diet.plant,
                                        predictors$trait.func.div.mpd.sans.dens,
                                        predictors$bird.rich)
                                   )

predictors.multi.ppms.func

pairs(predictors.multi.ppms.func)

# Next extract data at presence and background points to examine VIF:

bg.data <- terra::extract(predictors.multi.ppms.func, vect(bg.pts))
head(bg.data)
dim(bg.data)

pres.data <- terra::extract(predictors.multi.ppms.func, vect(jev_sf))
head(pres.data)
dim(pres.data)

# VIF for outbreak locations:

set.seed(4)
vif.pres <- vifstep(pres.data)
vif.pres

# VIF for background locations:

set.seed(4)
vif.abs <- vifstep(bg.data)
vif.abs

###############################
## inhomogeneous multiple PPMs
###############################

##################
# diversity model
##################

jev.div <- ppm(jev.q, ~ trait.mpd.sans.dens + bird.rich + I(bird.rich^2), 
               covariates = list(trait.mpd.sans.dens = trait.mpd.sans.dens.img,
                                 bird.rich = bird.rich.img))

jev.div
exp(coef(jev.div))
exp(confint(jev.div))
AIC(jev.div)

# stepwise model
jev.step.1 <- step(jev.div)

jev.step.1
exp(coef(jev.step.1))
exp(confint(jev.step.1))
AIC(jev.step.1)

##############
# trait model
##############

jev.trait <- ppm(jev.q, ~ trait.hwi + trait.below.water + trait.diet.plant,
                 covariates = list(trait.hwi = trait.hwi.img,
                                   trait.below.water = trait.below.water.surface.img,
                                   trait.diet.plant = trait.diet.plant.img))
                

jev.trait
exp(coef(jev.trait))
exp(confint(jev.trait))
AIC(jev.trait)

# stepwise model
jev.step.2 <- step(jev.trait)

jev.step.2
exp(coef(jev.step.2))
exp(confint(jev.step.2))
AIC(jev.step.2)

#############
# Full model
#############

jev.full <- ppm(jev.q, ~ trait.hwi + trait.diet.plant + trait.mpd.sans.dens + bird.rich, 
                           covariates = list(trait.hwi = trait.hwi.img,
                                             trait.below.water = trait.below.water.surface.img,
                                             trait.diet.plant = trait.diet.plant.img,
                                             trait.mpd.sans.dens = trait.mpd.sans.dens.img,
                                             bird.rich = bird.rich.img))

jev.full
exp(coef(jev.full))
exp(confint(jev.full))
AIC(jev.full)

# stepwise model
jev.step.3 <- step(jev.full)

jev.step.3
exp(coef(jev.step.3))
exp(confint(jev.step.3))
AIC(jev.step.3)

##############
# final model
##############

final.model <- ppm(jev.q, ~ trait.hwi + trait.diet.plant + bird.rich, 
                   covariates = list(trait.hwi = trait.hwi.img,
                                     trait.diet.plant = trait.diet.plant.img,
                                     bird.rich = bird.rich.img))

final.model
exp(coef(final.model))
exp(confint(final.model))
AIC(final.model)
