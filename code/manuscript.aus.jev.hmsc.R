library(sf)      
library(terra)
library(tidyverse)
library(tmap)
library(Hmsc)
library(ape)
library(phytools)

#######
# HMSC
#######

# Now we need to set up the data structures for Hmsc before we can implement the models:

# Let's look at correlation between traits and between environmental covariates:

head(trait.data.jsdm)
cor(trait.data.jsdm)
cor(x.data.jsdm[,1:5])

# We're modelling how species pools associated with landscape features that have been found to influence JEV outbreaks, so
# we'll create 2 model structures - one for JEV-associated landscape features and the second for the JEV-occurrence 
# landscapes directly:

# set up model data then set up model

# Model 1: Landscape feature-based

# formula for JEV-associated landscape features model
X.formula.1 <- ~ Croplands + RCC_grassland + Flow_accumulation + Transient_wetlands + Dist_to_waterways

# formula for traits
trait.formula <- ~ Diet_inverts + Below_water + Body_mass_ln + HWI + Studiedness

model.land = Hmsc(Y = y.abund.jsdm,
                 XData = x.data.jsdm,  XFormula = X.formula.1,
                 TrData = trait.data.jsdm, TrFormula = trait.formula,
                 phyloTree = consensus.tree,
                 distr="probit", 
                 studyDesign = study.design
)

model.land

# Model 1: Landscape feature-based

nChains <- 4
thin <- 30 
samples <- 1000 
transient <- 200*thin 
verbose <- 200*thin 

set.seed(4)
model.1 <- sampleMcmc(model.land, 
                      thin = thin, 
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      verbose = verbose)

model.1

# convert the posterior to a coda object to examine model estimates:
model.1.post <- convertToCodaObject(model.1)

# examine convergence
gelman.diag(model.1.post$Beta)$psrf

# examine parameter estimates:
round(summary(model.1.post$Beta[[1]])$quantiles, 3)
round(summary(model.1.post$Gamma[[1]])$quantiles, 3)

# Let's look at the extent to which the species occurrences depend on the environmental covariates

post.beta.1 <- getPostEstimate(model.1, parName = "Beta")
plotBeta(model.1, post = post.beta.1, param = "Sign", plotTree = T, supportLevel = 0.95, split = 0.4, spNamesNumbers = c(T, F))

# Now let's look at the extent to which species niches depend on species traits

post.gamma.1 <- getPostEstimate(model.1, parName = "Gamma")
plotGamma(model.1, post = post.gamma.1, param = "Support")

# evaluate model fit:
# predicted values:
preds.1 <- computePredictedValues(model.1)
fit.model.1 <- evaluateModelFit(hM = model.1, predY = preds.1)

fit.model.1$TjurR2
summary(fit.model.1$TjurR2)

# Data 2: JEV-based

# formula for JEV occurrence only model
X.formula.2 <- ~ JEV_present

# formula for traits
trait.formula <- ~ Diet_inverts + Below_water + Body_mass_ln + HWI + Studiedness

model.jev = Hmsc(Y = y.abund.jsdm,
                 XData = x.data.jsdm,  XFormula = X.formula.2,
                 TrData = trait.data.jsdm, TrFormula = trait.formula,
                 phyloTree = consensus.tree,
                 distr="probit",
                 studyDesign = study.design
)

model.jev

# Model 2: JEV-based

nChains <- 4
thin <- 4 # model converges with less thinned than model 1 above
samples <- 1000 
transient <- 200*thin 
verbose <- 200*thin 

set.seed(4)
model.2 <- sampleMcmc(model.jev, 
                      thin = thin, 
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      verbose = verbose)

model.2

# convert the posterior to a coda object to examine model estimates:
model.2.post <- convertToCodaObject(model.2)

# examine convergence
gelman.diag(model.2.post$Beta)$psrf

# examine parameter estimates:
round(summary(model.2.post$Beta[[1]])$quantiles, 3)
round(summary(model.2.post$Gamma[[1]])$quantiles, 3)

# Let's look at the extent to which the species occurrences depend on the environmental covariates

post.beta.2 <- getPostEstimate(model.2, parName = "Beta")
plotBeta(model.2, post = post.beta.2, param = "Sign", plotTree = T, supportLevel = 0.95, split = 0.4, spNamesNumbers = c(T, F))

# Now let's look at the extent to which species "niches" depend on species traits

post.gamma.2 <- getPostEstimate(model.2, parName = "Gamma")
plotGamma(model.2, post = post.gamma.2, param = "Sign", supportLevel = 0.95)

# evaluate model fit:
# predicted values:

preds.2 <- computePredictedValues(model.2)
fit.model.2 <- evaluateModelFit(hM = model.2, predY = preds.2)

fit.model.2$TjurR2
summary(fit.model.2$TjurR2)

