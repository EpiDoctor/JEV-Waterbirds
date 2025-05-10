library(tidyverse)
library(sf)
library(tmap)
library(terra)
library(viridis)
library(spatstat)

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

################################################################################
## Bivariate PPMs for "species pool presence" (SPP) and landscape scale of 5 km
################################################################################

# Ardeid models

jev.1 <- ppm(jev.q, ~ ardea.alba.spp + offset(log(pig.dens)), 
             covariates = list(ardea.alba.spp = pixel_list[["ardea.alba"]], pig.dens = piggery.dens.img))

jev.1
AIC(jev.1)

jev.2 <- ppm(jev.q, ~ ardea.pacifica.spp + offset(log(pig.dens)), 
             covariates = list(ardea.pacifica.spp = pixel_list[["ardea.pacifica"]], pig.dens = piggery.dens.img))

jev.2
AIC(jev.2)

jev.3 <- ppm(jev.q, ~ ardea.sumatrana.spp + offset(log(pig.dens)), 
             covariates = list(ardea.sumatrana.spp = pixel_list[["ardea.sumatrana"]], pig.dens = piggery.dens.img))

jev.3
AIC(jev.3)

jev.4 <- ppm(jev.q, ~ botaurus.poiciloptilus.spp + offset(log(pig.dens)), 
             covariates = list(botaurus.poiciloptilus.spp = pixel_list[["botaurus.poiciloptilus"]], pig.dens = piggery.dens.img))
jev.4
AIC(jev.4)

jev.5 <- ppm(jev.q, ~ bubulcus.coromandus.spp + offset(log(pig.dens)), 
             covariates = list(bubulcus.coromandus.spp = pixel_list[["bubulcus.ibis"]], pig.dens = piggery.dens.img))
jev.5
AIC(jev.5)

jev.6 <- ppm(jev.q, ~ butorides.striata.spp + offset(log(pig.dens)), 
             covariates = list(butorides.striata.spp = pixel_list[["butorides.striata"]], pig.dens = piggery.dens.img))
jev.6
AIC(jev.6)

jev.7 <- ppm(jev.q, ~ dupetor.flavicollis.spp + offset(log(pig.dens)), 
             covariates = list(dupetor.flavicollis.spp = pixel_list[["dupetor.flavicollis"]], pig.dens = piggery.dens.img))
jev.7
AIC(jev.7)

jev.8 <- ppm(jev.q, ~ egretta.garzetta.spp + offset(log(pig.dens)), 
             covariates = list(egretta.garzetta.spp = pixel_list[["egretta.garzetta"]], pig.dens = piggery.dens.img))
jev.8
AIC(jev.8)

jev.9 <- ppm(jev.q, ~ egretta.intermedia.spp + offset(log(pig.dens)), 
             covariates = list(egretta.intermedia.spp = pixel_list[["egretta.intermedia"]], pig.dens = piggery.dens.img))
jev.9
AIC(jev.9)

jev.10 <- ppm(jev.q, ~ egretta.novaehollandiae.spp + offset(log(pig.dens)), 
              covariates = list(egretta.novaehollandiae.spp = pixel_list[["egretta.novaehollandiae"]], pig.dens = piggery.dens.img))
jev.10
AIC(jev.10)

jev.11 <- ppm(jev.q, ~ egretta.picata.spp + offset(log(pig.dens)), 
              covariates = list(egretta.picata.spp = pixel_list[["egretta.picata"]], pig.dens = piggery.dens.img))
jev.11
AIC(jev.11)

jev.12 <- ppm(jev.q, ~ egretta.sacra.spp + offset(log(pig.dens)), 
              covariates = list(egretta.sacra.spp = pixel_list[["egretta.sacra"]], pig.dens = piggery.dens.img))
jev.12
AIC(jev.12)

jev.13 <- ppm(jev.q, ~ ixobrychus.dubius.spp + offset(log(pig.dens)), 
              covariates = list(ixobrychus.dubius.spp = pixel_list[["ixobrychus.dubius"]], pig.dens = piggery.dens.img))
jev.13
AIC(jev.13)

jev.14 <- ppm(jev.q, ~ nycticorax.caledonicus.spp + offset(log(pig.dens)), 
              covariates = list(nycticorax.caledonicus.spp = pixel_list[["nycticorax.caledonicus"]], pig.dens = piggery.dens.img))
jev.14
AIC(jev.14)

# Anatid models

jev.15 <- ppm(jev.q, ~ anas.castanea.spp + offset(log(pig.dens)), 
              covariates = list(anas.castanea.spp = pixel_list[["anas.castanea"]], pig.dens = piggery.dens.img))
jev.15
AIC(jev.15)

jev.16 <- ppm(jev.q, ~ anas.gracilis.spp + offset(log(pig.dens)), 
              covariates = list(anas.gracilis.spp = pixel_list[["anas.gracilis"]], pig.dens = piggery.dens.img))
jev.16
AIC(jev.16)

jev.17 <- ppm(jev.q, ~ anas.platyrhynchos.spp + offset(log(pig.dens)), 
              covariates = list(anas.platyrhynchos.spp = pixel_list[["anas.platyrhynchos"]], pig.dens = piggery.dens.img))
jev.17
AIC(jev.17)

jev.18 <- ppm(jev.q, ~ anas.superciliosa.spp + offset(log(pig.dens)), 
              covariates = list(anas.superciliosa.spp = pixel_list[["anas.superciliosa"]], pig.dens = piggery.dens.img))
jev.18
AIC(jev.18)

jev.19 <- ppm(jev.q, ~ anser.anser.spp + offset(log(pig.dens)), 
              covariates = list(anser.anser.spp = pixel_list[["anser.anser"]], pig.dens = piggery.dens.img))
jev.19
AIC(jev.19)

jev.20 <- ppm(jev.q, ~ anseranas.semipalmata.spp + offset(log(pig.dens)), 
              covariates = list(anseranas.semipalmata.spp = pixel_list[["anseranas.semipalmata"]], pig.dens = piggery.dens.img))
jev.20
AIC(jev.20)

jev.21 <- ppm(jev.q, ~ aythya.australis.spp + offset(log(pig.dens)), 
              covariates = list(aythya.australis.spp = pixel_list[["aythya.australis"]], pig.dens = piggery.dens.img))
jev.21
AIC(jev.21)

jev.22 <- ppm(jev.q, ~ aythya.fuligula.spp + offset(log(pig.dens)), 
              covariates = list(aythya.fuligula.spp = pixel_list[["aythya.fuligula"]], pig.dens = piggery.dens.img))
jev.22
AIC(jev.22)

jev.23 <- ppm(jev.q, ~ biziura.lobata.spp + offset(log(pig.dens)), 
              covariates = list(biziura.lobata.spp = pixel_list[["biziura.lobata"]], pig.dens = piggery.dens.img))
jev.23
AIC(jev.23)

jev.24 <- ppm(jev.q, ~ cairina.moschata.spp + offset(log(pig.dens)), 
              covariates = list(cairina.moschata.spp = pixel_list[["cairina.moschata"]], pig.dens = piggery.dens.img))
jev.24
AIC(jev.24)

jev.25 <- ppm(jev.q, ~ cereopsis.novaehollandiae.spp + offset(log(pig.dens)), 
              covariates = list(cereopsis.novaehollandiae.spp = pixel_list[["cereopsis.novaehollandiae"]], pig.dens = piggery.dens.img))
jev.25
AIC(jev.25)

jev.26 <- ppm(jev.q, ~ chenonetta.jubata.spp + offset(log(pig.dens)), 
              covariates = list(chenonetta.jubata.spp = pixel_list[["chenonetta.jubata"]], pig.dens = piggery.dens.img))
jev.26
AIC(jev.26)

jev.27 <- ppm(jev.q, ~ cygnus.atratus.spp + offset(log(pig.dens)), 
              covariates = list(cygnus.atratus.spp = pixel_list[["cygnus.atratus"]], pig.dens = piggery.dens.img))
jev.27
AIC(jev.27)

jev.28 <- ppm(jev.q, ~ dendrocygna.arcuata.spp + offset(log(pig.dens)), 
              covariates = list(dendrocygna.arcuata.spp = pixel_list[["dendrocygna.arcuata"]], pig.dens = piggery.dens.img))
jev.28
AIC(jev.28)

jev.29 <- ppm(jev.q, ~ dendrocygna.eytoni.spp + offset(log(pig.dens)), 
              covariates = list(dendrocygna.eytoni.spp = pixel_list[["dendrocygna.eytoni"]], pig.dens = piggery.dens.img))
jev.29
AIC(jev.29)

jev.30 <- ppm(jev.q, ~ dendrocygna.guttata.spp + offset(log(pig.dens)), 
              covariates = list(dendrocygna.guttata.spp = pixel_list[["dendrocygna.guttata"]], pig.dens = piggery.dens.img))
jev.30
AIC(jev.30)

jev.31 <- ppm(jev.q, ~ malacorhynchus.membranaceus.spp + offset(log(pig.dens)), 
              covariates = list(malacorhynchus.membranaceus.spp = pixel_list[["malacorhynchus.membranaceus"]], pig.dens = piggery.dens.img))
jev.31
AIC(jev.31)

jev.32 <- ppm(jev.q, ~ nettapus.coromandelianus.spp + offset(log(pig.dens)), 
              covariates = list(nettapus.coromandelianus.spp = pixel_list[["nettapus.coromandelianus"]], pig.dens = piggery.dens.img))
jev.32
AIC(jev.32)

jev.33 <- ppm(jev.q, ~ nettapus.pulchellus.spp + offset(log(pig.dens)), 
              covariates = list(nettapus.pulchellus.spp = pixel_list[["nettapus.pulchellus"]], pig.dens = piggery.dens.img))
jev.33
AIC(jev.33)

jev.34 <- ppm(jev.q, ~ oxyura.australis.spp + offset(log(pig.dens)), 
              covariates = list(oxyura.australis.spp = pixel_list[["oxyura.australis"]], pig.dens = piggery.dens.img))
jev.34
AIC(jev.34)

jev.35 <- ppm(jev.q, ~ radjah.radjah.spp + offset(log(pig.dens)), 
              covariates = list(radjah.radjah.spp = pixel_list[["radjah.radjah"]], pig.dens = piggery.dens.img))
jev.35
AIC(jev.35)

jev.36 <- ppm(jev.q, ~ spatula.clypeata.spp + offset(log(pig.dens)), 
              covariates = list(spatula.clypeata.spp = pixel_list[["spatula.clypeata"]], pig.dens = piggery.dens.img))
jev.36
AIC(jev.36)

jev.37 <- ppm(jev.q, ~ spatula.querquedula.spp + offset(log(pig.dens)), 
              covariates = list(spatula.querquedula.spp = pixel_list[["spatula.querquedula"]], pig.dens = piggery.dens.img))
jev.37
AIC(jev.37)

jev.38 <- ppm(jev.q, ~ spatula.rhynchotis.spp + offset(log(pig.dens)), 
              covariates = list(spatula.rhynchotis.spp = pixel_list[["spatula.rhynchotis"]], pig.dens = piggery.dens.img))
jev.38
AIC(jev.38)

jev.39 <- ppm(jev.q, ~ stictonetta.naevosa.spp + offset(log(pig.dens)), 
              covariates = list(stictonetta.naevosa.spp = pixel_list[["stictonetta.naevosa"]], pig.dens = piggery.dens.img))
jev.39
AIC(jev.39)

jev.40 <- ppm(jev.q, ~ tadorna.tadornoides.spp + offset(log(pig.dens)), 
              covariates = list(tadorna.tadornoides.spp = pixel_list[["tadorna.tadornoides"]], pig.dens = piggery.dens.img))
jev.40
AIC(jev.40)

# rallids

jev.41 <- ppm(jev.q, ~ fulica.atra.spp + offset(log(pig.dens)), 
              covariates = list(fulica.atra.spp = pixel_list[["fulica.atra"]], pig.dens = piggery.dens.img))
jev.41
AIC(jev.41)

jev.42 <- ppm(jev.q, ~ gallinula.tenebrosa.spp + offset(log(pig.dens)), 
              covariates = list(gallinula.tenebrosa.spp = pixel_list[["gallinula.tenebrosa"]], pig.dens = piggery.dens.img))
jev.42
AIC(jev.42)

jev.43 <- ppm(jev.q, ~ gallinula.ventralis.spp + offset(log(pig.dens)), 
              covariates = list(gallinula.ventralis.spp = pixel_list[["gallinula.ventralis"]], pig.dens = piggery.dens.img))
jev.43
AIC(jev.43)

jev.44 <- ppm(jev.q, ~ gallirallus.castaneoventris.spp + offset(log(pig.dens)), 
              covariates = list(gallirallus.castaneoventris.spp = pixel_list[["gallirallus.castaneoventris"]], pig.dens = piggery.dens.img))
jev.44
AIC(jev.44)

jev.45 <- ppm(jev.q, ~ gallirallus.philippensis.spp + offset(log(pig.dens)), 
              covariates = list(gallirallus.philippensis.spp = pixel_list[["gallirallus.philippensis"]], pig.dens = piggery.dens.img))
jev.45
AIC(jev.45)

jev.46 <- ppm(jev.q, ~ lewinia.pectoralis.spp + offset(log(pig.dens)), 
              covariates = list(lewinia.pectoralis.spp = pixel_list[["lewinia.pectoralis"]], pig.dens = piggery.dens.img))
jev.46
AIC(jev.46)

jev.47 <- ppm(jev.q, ~ porphyrio.melanotus.spp + offset(log(pig.dens)), 
              covariates = list(porphyrio.melanotus.spp = pixel_list[["porphyrio.melanotus"]], pig.dens = piggery.dens.img))
jev.47
AIC(jev.47)

jev.48 <- ppm(jev.q, ~ porzana.cinerea.spp + offset(log(pig.dens)), 
              covariates = list(porzana.cinerea.spp = pixel_list[["porzana.cinerea"]], pig.dens = piggery.dens.img))
jev.48
AIC(jev.48)

jev.49 <- ppm(jev.q, ~ porzana.fluminea.spp + offset(log(pig.dens)), 
              covariates = list(porzana.fluminea.spp = pixel_list[["porzana.fluminea"]], pig.dens = piggery.dens.img))
jev.49
AIC(jev.49)

jev.50 <- ppm(jev.q, ~ porzana.pusilla.spp + offset(log(pig.dens)), 
              covariates = list(porzana.pusilla.spp = pixel_list[["porzana.pusilla"]], pig.dens = piggery.dens.img))
jev.50
AIC(jev.50)

jev.51 <- ppm(jev.q, ~ porzana.tabuensis.spp + offset(log(pig.dens)), 
              covariates = list(porzana.tabuensis.spp = pixel_list[["porzana.tabuensis"]], pig.dens = piggery.dens.img))
jev.51
AIC(jev.51)

jev.52 <- ppm(jev.q, ~ rallina.tricolor.spp + offset(log(pig.dens)), 
              covariates = list(rallina.tricolor.spp = pixel_list[["rallina.tricolor"]], pig.dens = piggery.dens.img))
jev.52
AIC(jev.52)

# phalacrocoracids

jev.53 <- ppm(jev.q, ~ microcarbo.melanoleucos.spp + offset(log(pig.dens)), 
              covariates = list(microcarbo.melanoleucos.spp = pixel_list[["microcarbo.melanoleucos"]], pig.dens = piggery.dens.img))
jev.53
AIC(jev.53)

jev.54 <- ppm(jev.q, ~ phalacrocorax.carbo.spp + offset(log(pig.dens)), 
              covariates = list(phalacrocorax.carbo.spp = pixel_list[["phalacrocorax.carbo"]], pig.dens = piggery.dens.img))
jev.54
AIC(jev.54)

jev.55 <- ppm(jev.q, ~ phalacrocorax.fuscescens.spp + offset(log(pig.dens)), 
              covariates = list(phalacrocorax.fuscescens.spp = pixel_list[["phalacrocorax.fuscescens"]], pig.dens = piggery.dens.img))
jev.55
AIC(jev.55)

jev.56 <- ppm(jev.q, ~ phalacrocorax.sulcirostris.spp + offset(log(pig.dens)), 
              covariates = list(phalacrocorax.sulcirostris.spp = pixel_list[["phalacrocorax.sulcirostris"]], pig.dens = piggery.dens.img))
jev.56
AIC(jev.56)

jev.57 <- ppm(jev.q, ~ phalacrocorax.varius.spp + offset(log(pig.dens)), 
              covariates = list(phalacrocorax.varius.spp = pixel_list[["phalacrocorax.varius"]], pig.dens = piggery.dens.img))
jev.57
AIC(jev.57)

# threskis

jev.58 <- ppm(jev.q, ~ platalea.flavipes.spp + offset(log(pig.dens)), 
              covariates = list(platalea.flavipes.spp = pixel_list[["platalea.flavipes"]], pig.dens = piggery.dens.img))
jev.58
AIC(jev.58)

jev.59 <- ppm(jev.q, ~ platalea.regia.spp + offset(log(pig.dens)), 
              covariates = list(platalea.regia.spp = pixel_list[["platalea.regia"]], pig.dens = piggery.dens.img))
jev.59
AIC(jev.59)

jev.60 <- ppm(jev.q, ~ plegadis.falcinellus.spp + offset(log(pig.dens)), 
              covariates = list(plegadis.falcinellus.spp = pixel_list[["plegadis.falcinellus"]], pig.dens = piggery.dens.img))
jev.60
AIC(jev.60)

jev.61 <- ppm(jev.q, ~ threskiornis.molucca.spp + offset(log(pig.dens)), 
              covariates = list(threskiornis.molucca.spp = pixel_list[["threskiornis.molucca"]], pig.dens = piggery.dens.img))
jev.61
AIC(jev.61)

jev.62 <- ppm(jev.q, ~ threskiornis.spinicollis.spp + offset(log(pig.dens)), 
              covariates = list(threskiornis.spinicollis.spp = pixel_list[["threskiornis.spinicollis"]], pig.dens = piggery.dens.img))
jev.62
AIC(jev.62)

# others

jev.63 <- ppm(jev.q, ~ ephippiorhynchus.asiaticus.spp + offset(log(pig.dens)), 
              covariates = list(ephippiorhynchus.asiaticus.spp = pixel_list[["ephippiorhynchus.asiaticus"]], pig.dens = piggery.dens.img))
jev.63
AIC(jev.63)

jev.64 <- ppm(jev.q, ~ grus.antigone.spp + offset(log(pig.dens)), 
              covariates = list(grus.antigone.spp = pixel_list[["grus.antigone"]], pig.dens = piggery.dens.img))
jev.64
AIC(jev.64)

jev.65 <- ppm(jev.q, ~ grus.rubicunda.spp + offset(log(pig.dens)), 
              covariates = list(grus.rubicunda.spp = pixel_list[["grus.rubicunda"]], pig.dens = piggery.dens.img))
jev.65
AIC(jev.65)


jev.66 <- ppm(jev.q, ~ pelecanus.conspicillatus.spp + offset(log(pig.dens)), 
              covariates = list(pelecanus.conspicillatus.spp = pixel_list[["pelecanus.conspicillatus"]], pig.dens = piggery.dens.img))
jev.66
AIC(jev.66)
