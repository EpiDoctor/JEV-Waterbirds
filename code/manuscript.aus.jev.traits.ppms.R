library(sf)
library(terra)   
library(tidyverse)
library(tmap)
library(spatstat)
library(usdm)
library(viridis)

###########################
## Kernel density estimate
###########################

jev.intens <- density.ppp(jev.ppp, 
                          sigma = 2, 
                          diggle = T)

jev.intens.ras <- rast(jev.intens)
crs(jev.intens.ras)  <- "EPSG:4326"

jev.intens.ras
jev.intens.ras <- mask(jev.intens.ras, vect(east))
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

null <- ppm(jev.q, ~1 + offset(log(pig.dens)), covariates = list(pig.dens = piggery.dens.img)) 
summary(null)

AIC(null)

##################
## Bivariate PPMs
##################

# Functional trait models with reduced background data capturing locations with ardeid communities only

jev.1 <- ppm(jev.q, ~ trait.diet.fish + offset(log(pig.dens)), covariates = list(trait.diet.fish = trait.diet.fish.img, 
                                                                                          pig.dens = piggery.dens.img))
jev.1
AIC(jev.1)

jev.2 <- ppm(jev.q, ~ trait.diet.inv + offset(log(pig.dens)), covariates = list(trait.diet.inv = trait.diet.inv.img, 
                                                                                 pig.dens = piggery.dens.img))
jev.2
AIC(jev.2)

jev.3 <- ppm(jev.q, ~ trait.diet.plant + offset(log(pig.dens)), covariates = list(trait.diet.plant = trait.diet.plant.img, 
                                                                                pig.dens = piggery.dens.img))
jev.3
AIC(jev.3)

jev.4 <- ppm(jev.q, ~ trait.body.mass.lg + offset(log(pig.dens)), covariates = list(trait.body.mass.lg = trait.body.mass.lg.img, 
                                                                               pig.dens = piggery.dens.img))
jev.4
AIC(jev.4)

jev.5 <- ppm(jev.q, ~ trait.egg.mass + offset(log(pig.dens)), covariates = list(trait.egg.mass = trait.egg.mass.img, 
                                                                               pig.dens = piggery.dens.img))
jev.5
AIC(jev.5)

jev.6 <- ppm(jev.q, ~ trait.hwi + offset(log(pig.dens)), covariates = list(trait.hwi = trait.hwi.img, 
                                                                               pig.dens = piggery.dens.img))
jev.6
AIC(jev.6)

jev.7 <- ppm(jev.q, ~ trait.water.surface + offset(log(pig.dens)), covariates = list(trait.water.surface = trait.water.surface.img, 
                                                                                          pig.dens = piggery.dens.img))
jev.7
AIC(jev.7)

jev.8 <- ppm(jev.q, ~ trait.below.water.surface + offset(log(pig.dens)), covariates = list(trait.below.water.surface = trait.below.water.surface.img, 
                                                                                     pig.dens = piggery.dens.img))
jev.8
AIC(jev.8)

jev.9 <- ppm(jev.q, ~ trait.ground + offset(log(pig.dens)), covariates = list(trait.ground = trait.ground.img, 
                                                                                           pig.dens = piggery.dens.img))
jev.9
AIC(jev.9)

# population density no longer included because of suboptimal imputation of this trait
jev.10 <- ppm(jev.q, ~ trait.mean.density + offset(log(pig.dens)), covariates = list(trait.mean.density = trait.mean.density.img, 
                                                                                     pig.dens = piggery.dens.img))
jev.10
AIC(jev.10)

# MPD computed without population density because of suboptimal imputation of this trait
jev.11 <- ppm(jev.q, ~ trait.mpd.sans.dens + offset(log(pig.dens)), covariates = list(trait.mpd.sans.dens = trait.mpd.sans.dens.img, 
                                                                             pig.dens = piggery.dens.img))
jev.11
AIC(jev.11)

jev.12 <- ppm(jev.q, ~ trait.mpd.sans.dens + I(trait.mpd.sans.dens^2) + offset(log(pig.dens)), covariates = list(trait.mpd.sans.dens = trait.mpd.sans.dens.img, 
                                                                            pig.dens = piggery.dens.img))
jev.12
AIC(jev.12)

jev.13 <- ppm(jev.q, ~ bird.rich + offset(log(pig.dens)), covariates = list(bird.rich = bird.rich.img, 
                                                                            pig.dens = piggery.dens.img))
jev.13
AIC(jev.13)

jev.14 <- ppm(jev.q, ~ bird.rich + I(bird.rich^2) + offset(log(pig.dens)), covariates = list(bird.rich = bird.rich.img, 
                                                                            pig.dens = piggery.dens.img))
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

tiff("Predictor correlation Aus community metrics func div.tiff", width = 12, height = 12, units = 'in', res = 300)
pairs(predictors.multi.ppms.func)
dev.off()

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

jev.div <- ppm(jev.q, ~ trait.mpd.sans.dens + bird.rich + I(bird.rich^2) + offset(log(pig.dens)), 
                  covariates = list(trait.mpd.sans.dens = trait.mpd.sans.dens.img,
                                    bird.rich = bird.rich.img,
                                    pig.dens = piggery.dens.img))

jev.div
exp(coef(jev.div))
exp(confint(jev.div))
AIC(jev.div)

p.jev.div <- predict.ppm(jev.div, covariates = list(trait.mpd.sans.dens = trait.mpd.sans.dens.img,
                                                    bird.rich = bird.rich.img,
                                                    pig.dens = piggery.dens.img), 
                            locations = piggery.dens.img)
# Test model performance:

set.seed(2)
auc.jev.div <- auc(jev.test.ppp, p.jev.div)
auc.jev.div

# stepwise model
jev.step.1 <- step(jev.div)

jev.step.1
exp(coef(jev.step.1))
exp(confint(jev.step.1))
AIC(jev.step.1)

##############
# trait model
##############

jev.trait <- ppm(jev.q, ~ trait.hwi + trait.below.water.surface + offset(log(pig.dens)), 
                     covariates = list(trait.hwi = trait.hwi.img,
                                       trait.below.water.surface = trait.below.water.surface.img,
                                       pig.dens = piggery.dens.img))

jev.trait
exp(coef(jev.trait))
exp(confint(jev.trait))
AIC(jev.trait)

p.jev.trait <- predict.ppm(jev.trait, covariates = list(trait.hwi = trait.hwi.img,
                                                        trait.below.water.surface = trait.below.water.surface.img,
                                                        pig.dens = piggery.dens.img), 
                         locations = piggery.dens.img)
# Test model performance:

set.seed(2)
auc.jev.trait <- auc(jev.test.ppp, p.jev.trait)
auc.jev.trait

# stepwise model
jev.step.2 <- step(jev.trait)

jev.step.2
exp(coef(jev.step.2))
exp(confint(jev.step.2))
AIC(jev.step.2)

#############
# Full model
#############

jev.full <- ppm(jev.q, ~ trait.hwi + trait.below.water.surface + bird.rich + I(bird.rich^2) + offset(log(pig.dens)), 
                  covariates = list(trait.hwi = trait.hwi.img,
                                    trait.below.water.surface = trait.below.water.surface.img,
                                    bird.rich = bird.rich.img,
                                    pig.dens = piggery.dens.img))

jev.full
exp(coef(jev.full))
exp(confint(jev.full))
AIC(jev.full)

p.jev.full <- predict.ppm(jev.full, covariates = list(trait.hwi = trait.hwi.img,
                                                      trait.below.water.surface = trait.below.water.surface.img,
                                                      bird.rich = bird.rich.img,
                                                      pig.dens = piggery.dens.img), 
                            locations = piggery.dens.img)
# Test model performance:

set.seed(2)
auc.jev.full <- auc(jev.test.ppp, p.jev.full)
auc.jev.full

# stepwise model
jev.step.3 <- step(jev.full)

jev.step.3
exp(coef(jev.step.3))
exp(confint(jev.step.3))
AIC(jev.step.3)

# final model 2

final.model <- ppm(jev.q, ~ trait.hwi + bird.rich + I(bird.rich^2) + offset(log(pig.dens)), 
                     covariates = list(trait.hwi = trait.hwi.img,
                                       bird.rich = bird.rich.img,
                                       pig.dens = piggery.dens.img))

final.model
exp(coef(final.model))
exp(confint(final.model))
AIC(final.model)

# for estimates and performance
p.final.model <- predict.ppm(final.model, covariates = list(trait.hwi = trait.hwi.img,
                                                                bird.rich = bird.rich.img,
                                                                pig.dens = piggery.dens.img), 
                               locations = piggery.dens.img)
# Test model performance:

set.seed(2)
auc.final.model <- auc(jev.test.ppp, p.final.model)
auc.final.model

# for confidence limits
p.final.model.ci <- predict.ppm(final.model, covariates = list(trait.hwi = trait.hwi.img,
                                                               bird.rich = bird.rich.img,
                                                               pig.dens = piggery.dens.img), 
                                  locations = piggery.dens.img, interval = "confidence")


# Sensitivity analysis for model performance - testing performance against the lower confidence limit of the intensity:

set.seed(2)
auc.final.model.sn.lower <- auc(jev.test.ppp, p.final.model.ci$`2.5%`)
auc.final.model.sn.lower

set.seed(2)
auc.final.model.sn.upper <- auc(jev.test.ppp, p.final.model.ci$`97.5%`)
auc.final.model.sn.upper

# Convert pixel images to rasters for visualisations

pfinal.model.ras <- rast(p.final.model)
pfinal.model.lower <- rast(p.final.model.ci$`2.5%`)
pfinal.model.upper <- rast(p.final.model.ci$`97.5%`)

# create deciles

qnt.pFull.ras <- global(pfinal.model.ras, fun = quantile, probs = seq(0, 1, 0.10), na.rm = T)
m <- c(qnt.pFull.ras[1,1:2], 1,
       qnt.pFull.ras[1,2:3], 2,
       qnt.pFull.ras[1,3:4], 3,
       qnt.pFull.ras[1,4:5], 4,
       qnt.pFull.ras[1,5:6], 5,
       qnt.pFull.ras[1,6:7], 6,
       qnt.pFull.ras[1,7:8], 7,
       qnt.pFull.ras[1,8:9], 8,
       qnt.pFull.ras[1,9:10], 9,
       qnt.pFull.ras[1,10:11], 10)
cat.mat <- matrix(as.numeric(m), ncol=3, byrow=TRUE)
pFull.ras.cat <- classify(pfinal.model.ras, cat.mat, include.lowest=TRUE)
qnt.pFull.lower <- global(pfinal.model.lower, fun = quantile, probs = seq(0, 1, 0.10), na.rm = T)
m.l <- c(qnt.pFull.lower[1,1:2], 1,
       qnt.pFull.lower[1,2:3], 2,
       qnt.pFull.lower[1,3:4], 3,
       qnt.pFull.lower[1,4:5], 4,
       qnt.pFull.lower[1,5:6], 5,
       qnt.pFull.lower[1,6:7], 6,
       qnt.pFull.lower[1,7:8], 7,
       qnt.pFull.lower[1,8:9], 8,
       qnt.pFull.lower[1,9:10], 9,
       qnt.pFull.lower[1,10:11], 10)
cat.mat.l <- matrix(as.numeric(m.l), ncol=3, byrow=TRUE)
pFull.lower.cat <- classify(pfinal.model.lower, cat.mat.l, include.lowest=TRUE)
qnt.pFull.upper <- global(pfinal.model.upper, fun = quantile, probs = seq(0, 1, 0.10), na.rm = T)
m.u <- c(qnt.pFull.upper[1,1:2], 1,
         qnt.pFull.upper[1,2:3], 2,
         qnt.pFull.upper[1,3:4], 3,
         qnt.pFull.upper[1,4:5], 4,
         qnt.pFull.upper[1,5:6], 5,
         qnt.pFull.upper[1,6:7], 6,
         qnt.pFull.upper[1,7:8], 7,
         qnt.pFull.upper[1,8:9], 8,
         qnt.pFull.upper[1,9:10], 9,
         qnt.pFull.upper[1,10:11], 10)
cat.mat.u <- matrix(as.numeric(m.u), ncol=3, byrow=TRUE)
pFull.upper.cat <- classify(pfinal.model.upper, cat.mat.u, include.lowest=TRUE)

plot(pFull.lower.cat, col = rev(viridis(10, option = "A")))
plot(pFull.ras.cat, col = rev(viridis(10, option = "A")))
plot(pFull.upper.cat, col = rev(viridis(10, option = "A")))
