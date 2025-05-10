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

null <- ppm(jev.q, ~1) #+ offset(log(pig.dens)), covariates = list(pig.dens = piggery.dens.img)) 
summary(null)

AIC(null)

##################
## Bivariate PPMs
##################

# Waterbird models

jev.1 <- ppm(jev.q, ~ ardea.alba.spp, 
             covariates = list(ardea.alba.spp = pixel_list[["ardea.alba"]]))

jev.1
AIC(jev.1)

jev.2 <- ppm(jev.q, ~ ardea.pacifica.spp,
             covariates = list(ardea.pacifica.spp = pixel_list[["ardea.pacifica"]]))

jev.2
AIC(jev.2)

jev.3 <- ppm(jev.q, ~ ardea.sumatrana.spp,
             covariates = list(ardea.sumatrana.spp = pixel_list[["ardea.sumatrana"]]))

jev.3
AIC(jev.3)

jev.4 <- ppm(jev.q, ~ botaurus.poiciloptilus.spp,
             covariates = list(botaurus.poiciloptilus.spp = pixel_list[["botaurus.poiciloptilus"]]))

jev.4
AIC(jev.4)

jev.5 <- ppm(jev.q, ~ bubulcus.coromandus.spp,
             covariates = list(bubulcus.coromandus.spp = pixel_list[["bubulcus.ibis"]]))

jev.5
AIC(jev.5)

jev.6 <- ppm(jev.q, ~ butorides.striata.spp, 
             covariates = list(butorides.striata.spp = pixel_list[["butorides.striata"]]))

jev.6
AIC(jev.6)

jev.7 <- ppm(jev.q, ~ dupetor.flavicollis.spp, 
             covariates = list(dupetor.flavicollis.spp = pixel_list[["dupetor.flavicollis"]]))

jev.7
AIC(jev.7)

jev.8 <- ppm(jev.q, ~ egretta.garzetta.spp,
             covariates = list(egretta.garzetta.spp = pixel_list[["egretta.garzetta"]]))

jev.8
AIC(jev.8)

jev.9 <- ppm(jev.q, ~ egretta.intermedia.spp, 
             covariates = list(egretta.intermedia.spp = pixel_list[["egretta.intermedia"]]))

jev.9
AIC(jev.9)

jev.10 <- ppm(jev.q, ~ egretta.novaehollandiae.spp,
              covariates = list(egretta.novaehollandiae.spp = pixel_list[["egretta.novaehollandiae"]]))

jev.10
AIC(jev.10)

jev.11 <- ppm(jev.q, ~ egretta.picata.spp,
              covariates = list(egretta.picata.spp = pixel_list[["egretta.picata"]]))

jev.11
AIC(jev.11)

jev.12 <- ppm(jev.q, ~ egretta.sacra.spp,
              covariates = list(egretta.sacra.spp = pixel_list[["egretta.sacra"]]))

jev.12
AIC(jev.12)

jev.13 <- ppm(jev.q, ~ ixobrychus.dubius.spp,
              covariates = list(ixobrychus.dubius.spp = pixel_list[["ixobrychus.dubius"]]))

jev.13
AIC(jev.13)

jev.14 <- ppm(jev.q, ~ nycticorax.caledonicus.spp,
              covariates = list(nycticorax.caledonicus.spp = pixel_list[["nycticorax.caledonicus"]]))

jev.14
AIC(jev.14)

# out of Oz

jev.i <- ppm(jev.q, ~ ardea.cinerea.spp,
              covariates = list(ardea.cinerea.spp = pixel_list[["ardea.cinerea"]]))
jev.i
AIC(jev.i)

jev.ii <- ppm(jev.q, ~ ardea.purpurea.spp,
             covariates = list(ardea.purpurea.spp = pixel_list[["ardea.purpurea"]]))
jev.ii
AIC(jev.ii)

jev.iii <- ppm(jev.q, ~ ardeola.bacchus.spp,
             covariates = list(ardeola.bacchus.spp = pixel_list[["ardeola.bacchus"]]))
jev.iii
AIC(jev.iii)

jev.iv <- ppm(jev.q, ~ ardeola.speciosa.spp,
             covariates = list(ardeola.speciosa.spp = pixel_list[["ardeola.speciosa"]]))
jev.iv
AIC(jev.iv)

jev.v <- ppm(jev.q, ~ egretta.eulophotes.spp,
             covariates = list(egretta.eulophotes.spp = pixel_list[["egretta.eulophotes"]]))
jev.v
AIC(jev.v)

jev.vi <- ppm(jev.q, ~ ixobrychus.cinnamomeus.spp,
             covariates = list(ixobrychus.cinnamomeus.spp = pixel_list[["ixobrychus.cinnamomeus"]]))
jev.vi
AIC(jev.vi)

jev.vii <- ppm(jev.q, ~ ixobrychus.sinensis.spp,
             covariates = list(ixobrychus.sinensis.spp = pixel_list[["ixobrychus.sinensis"]]))
jev.vii
AIC(jev.vii)

jev.viii <- ppm(jev.q, ~ nycticorax.nycticorax.spp,
             covariates = list(nycticorax.nycticorax.spp = pixel_list[["nycticorax.nycticorax"]]))
jev.viii
AIC(jev.viii)

# Anatid models

jev.16 <- ppm(jev.q, ~ anas.castanea.spp,
              covariates = list(anas.castanea.spp = pixel_list[["anas.castanea"]]))
                                                                              
jev.16
AIC(jev.16)

jev.17 <- ppm(jev.q, ~ anas.gracilis.spp, 
              covariates = list(anas.gracilis.spp = pixel_list[["anas.gracilis"]]))
                                                                  
jev.17
AIC(jev.17)

jev.18 <- ppm(jev.q, ~ anas.platyrhynchos.spp,
              covariates = list(anas.platyrhynchos.spp = pixel_list[["anas.platyrhynchos"]]))
                                                                     
jev.18
AIC(jev.18)

jev.19 <- ppm(jev.q, ~ anas.superciliosa.spp,
              covariates = list(anas.superciliosa.spp = pixel_list[["anas.superciliosa"]]))

jev.19
AIC(jev.19)

jev.20 <- ppm(jev.q, ~ anser.anser.spp,
              covariates = list(anser.anser.spp = pixel_list[["anser.anser"]]))

jev.20
AIC(jev.20)

jev.21 <- ppm(jev.q, ~ anseranas.semipalmata.spp,
              covariates = list(anseranas.semipalmata.spp = pixel_list[["anseranas.semipalmata"]]))

jev.21
AIC(jev.21)

jev.22 <- ppm(jev.q, ~ aythya.australis.spp,
              covariates = list(aythya.australis.spp = pixel_list[["aythya.australis"]]))

jev.22
AIC(jev.22)

jev.23 <- ppm(jev.q, ~ aythya.fuligula.spp,
              covariates = list(aythya.fuligula.spp = pixel_list[["aythya.fuligula"]]))

jev.23
AIC(jev.23)

jev.24 <- ppm(jev.q, ~ biziura.lobata.spp, 
              covariates = list(biziura.lobata.spp = pixel_list[["biziura.lobata"]]))

jev.24
AIC(jev.24)

jev.25 <- ppm(jev.q, ~ cairina.moschata.spp, 
              covariates = list(cairina.moschata.spp = pixel_list[["cairina.moschata"]]))

jev.25
AIC(jev.25)

jev.26 <- ppm(jev.q, ~ cereopsis.novaehollandiae.spp,
              covariates = list(cereopsis.novaehollandiae.spp = pixel_list[["cereopsis.novaehollandiae"]]))

jev.26
AIC(jev.26)

jev.27 <- ppm(jev.q, ~ chenonetta.jubata.spp,
              covariates = list(chenonetta.jubata.spp = pixel_list[["chenonetta.jubata"]]))

jev.27
AIC(jev.27)

jev.28 <- ppm(jev.q, ~ cygnus.atratus.spp,
              covariates = list(cygnus.atratus.spp = pixel_list[["cygnus.atratus"]]))

jev.28
AIC(jev.28)

jev.29 <- ppm(jev.q, ~ dendrocygna.arcuata.spp, 
              covariates = list(dendrocygna.arcuata.spp = pixel_list[["dendrocygna.arcuata"]]))

jev.29
AIC(jev.29)

jev.30 <- ppm(jev.q, ~ dendrocygna.eytoni.spp,
              covariates = list(dendrocygna.eytoni.spp = pixel_list[["dendrocygna.eytoni"]]))

jev.30
AIC(jev.30)

jev.31 <- ppm(jev.q, ~ dendrocygna.guttata.spp,
              covariates = list(dendrocygna.guttata.spp = pixel_list[["dendrocygna.guttata"]]))

jev.31
AIC(jev.31)

jev.32 <- ppm(jev.q, ~ dendrocygna.javanica.spp,
              covariates = list(dendrocygna.javanica.spp = pixel_list[["dendrocygna.javanica"]]))

jev.32
AIC(jev.32)

jev.33 <- ppm(jev.q, ~ malacorhynchus.membranaceus.spp,
              covariates = list(malacorhynchus.membranaceus.spp = pixel_list[["malacorhynchus.membranaceus"]]))

jev.33
AIC(jev.33)

jev.34 <- ppm(jev.q, ~ nettapus.coromandelianus.spp,
              covariates = list(nettapus.coromandelianus.spp = pixel_list[["nettapus.coromandelianus"]])) 

jev.34
AIC(jev.34)

jev.35 <- ppm(jev.q, ~ nettapus.pulchellus.spp,
              covariates = list(nettapus.pulchellus.spp = pixel_list[["nettapus.pulchellus"]]))

jev.35
AIC(jev.35)

jev.36 <- ppm(jev.q, ~ oxyura.australis.spp,
              covariates = list(oxyura.australis.spp = pixel_list[["oxyura.australis"]]))

jev.36
AIC(jev.36)

jev.37 <- ppm(jev.q, ~ radjah.radjah.spp,
              covariates = list(radjah.radjah.spp = pixel_list[["radjah.radjah"]])) 

jev.37
AIC(jev.37)

jev.38 <- ppm(jev.q, ~ spatula.clypeata.spp,
              covariates = list(spatula.clypeata.spp = pixel_list[["spatula.clypeata"]])) 

jev.38
AIC(jev.38)

jev.39 <- ppm(jev.q, ~ spatula.querquedula.spp,
              covariates = list(spatula.querquedula.spp = pixel_list[["spatula.querquedula"]])) 

jev.39
AIC(jev.39)

jev.40 <- ppm(jev.q, ~ spatula.rhynchotis.spp,
              covariates = list(spatula.rhynchotis.spp = pixel_list[["spatula.rhynchotis"]])) 

jev.40
AIC(jev.40)

jev.41 <- ppm(jev.q, ~ stictonetta.naevosa.spp,
              covariates = list(stictonetta.naevosa.spp = pixel_list[["stictonetta.naevosa"]])) 

jev.41
AIC(jev.41)

jev.42 <- ppm(jev.q, ~ tadorna.tadornoides.spp,
              covariates = list(tadorna.tadornoides.spp = pixel_list[["tadorna.tadornoides"]])) 

jev.42
AIC(jev.42)

# out of Oz

jev.ix <- ppm(jev.q, ~ anas.gibberifrons.spp,
                covariates = list(anas.gibberifrons.spp = pixel_list[["anas.gibberifrons"]]))
jev.ix
AIC(jev.ix)

jev.x <- ppm(jev.q, ~ anastomus.oscitans.spp,
                covariates = list(anastomus.oscitans.spp = pixel_list[["anastomus.oscitans"]]))
jev.x
AIC(jev.x)

# rallids

jev.43 <- ppm(jev.q, ~ amaurornis.olivacea.spp,
              covariates = list(amaurornis.olivacea.spp = pixel_list[["amaurornis.olivacea"]]))

jev.43
AIC(jev.43)

jev.44 <- ppm(jev.q, ~ fulica.atra.spp,
              covariates = list(fulica.atra.spp = pixel_list[["fulica.atra"]]))

jev.44
AIC(jev.44)

jev.45 <- ppm(jev.q, ~ gallinula.mortierii.spp,
              covariates = list(gallinula.mortierii.spp = pixel_list[["gallinula.mortierii"]]))

jev.45
AIC(jev.45)

jev.46 <- ppm(jev.q, ~ gallinula.tenebrosa.spp,
              covariates = list(gallinula.tenebrosa.spp = pixel_list[["gallinula.tenebrosa"]]))

jev.46
AIC(jev.46)

jev.47 <- ppm(jev.q, ~ gallinula.ventralis.spp,
              covariates = list(gallinula.ventralis.spp = pixel_list[["gallinula.ventralis"]]))

jev.47
AIC(jev.47)

jev.48 <- ppm(jev.q, ~ gallirallus.castaneoventris.spp,
              covariates = list(gallirallus.castaneoventris.spp = pixel_list[["gallirallus.castaneoventris"]])) 

jev.48
AIC(jev.48)

jev.49 <- ppm(jev.q, ~ gallirallus.philippensis.spp,
              covariates = list(gallirallus.philippensis.spp = pixel_list[["gallirallus.philippensis"]])) 

jev.49
AIC(jev.49)

jev.50 <- ppm(jev.q, ~ lewinia.pectoralis.spp,
              covariates = list(lewinia.pectoralis.spp = pixel_list[["lewinia.pectoralis"]]))

jev.50
AIC(jev.50)

jev.51 <- ppm(jev.q, ~ porphyrio.melanotus.spp, 
              covariates = list(porphyrio.melanotus.spp = pixel_list[["porphyrio.melanotus"]]))

jev.51
AIC(jev.51)

jev.52 <- ppm(jev.q, ~ porzana.cinerea.spp,
              covariates = list(porzana.cinerea.spp = pixel_list[["porzana.cinerea"]]))

jev.52
AIC(jev.52)

jev.53 <- ppm(jev.q, ~ porzana.fluminea.spp,
              covariates = list(porzana.fluminea.spp = pixel_list[["porzana.fluminea"]]))

jev.53
AIC(jev.53)

jev.54 <- ppm(jev.q, ~ porzana.pusilla.spp,
              covariates = list(porzana.pusilla.spp = pixel_list[["porzana.pusilla"]]))

jev.54
AIC(jev.54)

jev.55 <- ppm(jev.q, ~ porzana.tabuensis.spp,
              covariates = list(porzana.tabuensis.spp = pixel_list[["porzana.tabuensis"]]))

jev.55
AIC(jev.55)

jev.56 <- ppm(jev.q, ~ rallina.tricolor.spp,
              covariates = list(rallina.tricolor.spp = pixel_list[["rallina.tricolor"]]))

jev.56
AIC(jev.56)

# out of Oz

jev.xi <- ppm(jev.q, ~ amaurornis.isabellina.spp,
                covariates = list(amaurornis.isabellina.spp = pixel_list[["amaurornis.isabellina"]]))
jev.xi
AIC(jev.xi)

jev.xii <- ppm(jev.q, ~ amaurornis.phoenicurus.spp,
                covariates = list(amaurornis.phoenicurus.spp = pixel_list[["amaurornis.phoenicurus"]]))
jev.xii
AIC(jev.xii)

jev.xiii <- ppm(jev.q, ~ gallicrex.cinerea.spp,
                covariates = list(gallicrex.cinerea.spp = pixel_list[["gallicrex.cinerea"]]))
jev.xiii
AIC(jev.xiii)

jev.xiv <- ppm(jev.q, ~ gallinula.chloropus.spp,
                covariates = list(gallinula.chloropus.spp = pixel_list[["gallinula.chloropus"]]))
jev.xiv
AIC(jev.xiv)

jev.xv <- ppm(jev.q, ~ gallirallus.striatus.spp,
                covariates = list(gallirallus.striatus.spp = pixel_list[["gallirallus.striatus"]]))
jev.xv
AIC(jev.xv)

jev.xvi <- ppm(jev.q, ~ gallirallus.torquatus.spp,
                covariates = list(gallirallus.torquatus.spp = pixel_list[["gallirallus.torquatus"]]))
jev.xvi
AIC(jev.xvi)

jev.xvii <- ppm(jev.q, ~ porphyrio.indicus.spp,
                covariates = list(porphyrio.indicus.spp = pixel_list[["porphyrio.indicus"]]))
jev.xvii
AIC(jev.xvii)

jev.xviii <- ppm(jev.q, ~ porzana.fusca.spp,
                covariates = list(porzana.fusca.spp = pixel_list[["porzana.fusca"]]))
jev.xviii
AIC(jev.xviii)

# phalacrocoracids

jev.57 <- ppm(jev.q, ~ microcarbo.melanoleucos.spp, 
              covariates = list(microcarbo.melanoleucos.spp = pixel_list[["microcarbo.melanoleucos"]])) 

jev.57
AIC(jev.57)

jev.58 <- ppm(jev.q, ~ phalacrocorax.carbo.spp,
              covariates = list(phalacrocorax.carbo.spp = pixel_list[["phalacrocorax.carbo"]]))

jev.58
AIC(jev.58)

jev.59 <- ppm(jev.q, ~ phalacrocorax.fuscescens.spp,
              covariates = list(phalacrocorax.fuscescens.spp = pixel_list[["phalacrocorax.fuscescens"]]))

jev.59
AIC(jev.59)

jev.60 <- ppm(jev.q, ~ phalacrocorax.sulcirostris.spp,
              covariates = list(phalacrocorax.sulcirostris.spp = pixel_list[["phalacrocorax.sulcirostris"]])) 

jev.60
AIC(jev.60)

jev.61 <- ppm(jev.q, ~ phalacrocorax.varius.spp,
              covariates = list(phalacrocorax.varius.spp = pixel_list[["phalacrocorax.varius"]])) 

jev.61
AIC(jev.61)

# out of Oz

jev.xix <- ppm(jev.q, ~ microcarbo.niger.spp,
                covariates = list(microcarbo.niger.spp = pixel_list[["microcarbo.niger"]]))
jev.xix
AIC(jev.xix)

# threskis

jev.62 <- ppm(jev.q, ~ platalea.flavipes.spp,
              covariates = list(platalea.flavipes.spp = pixel_list[["platalea.flavipes"]])) 

jev.62
AIC(jev.62)

jev.63 <- ppm(jev.q, ~ platalea.regia.spp,
              covariates = list(platalea.regia.spp = pixel_list[["platalea.regia"]])) 

jev.63
AIC(jev.63)

jev.64 <- ppm(jev.q, ~ plegadis.falcinellus.spp,
              covariates = list(plegadis.falcinellus.spp = pixel_list[["plegadis.falcinellus"]]))

jev.64
AIC(jev.64)

jev.65 <- ppm(jev.q, ~ threskiornis.molucca.spp,
              covariates = list(threskiornis.molucca.spp = pixel_list[["threskiornis.molucca"]]))

jev.65
AIC(jev.65)

jev.66 <- ppm(jev.q, ~ threskiornis.spinicollis.spp,
              covariates = list(threskiornis.spinicollis.spp = pixel_list[["threskiornis.spinicollis"]]))

jev.66
AIC(jev.66)

# others

jev.67 <- ppm(jev.q, ~ ephippiorhynchus.asiaticus.spp,
              covariates = list(ephippiorhynchus.asiaticus.spp = pixel_list[["ephippiorhynchus.asiaticus"]]))

jev.67
AIC(jev.67)

jev.68 <- ppm(jev.q, ~ grus.antigone.spp,
              covariates = list(grus.antigone.spp = pixel_list[["grus.antigone"]])) 

jev.68
AIC(jev.68)

jev.69 <- ppm(jev.q, ~ grus.rubicunda.spp,
              covariates = list(grus.rubicunda.spp = pixel_list[["grus.rubicunda"]]))

jev.69
AIC(jev.69)

jev.70 <- ppm(jev.q, ~ pelecanus.conspicillatus.spp, 
              covariates = list(pelecanus.conspicillatus.spp = pixel_list[["pelecanus.conspicillatus"]])) 

jev.70
AIC(jev.70)

# out of Oz

jev.xx <- ppm(jev.q, ~ leptoptilos.javanicus.spp,
                covariates = list(leptoptilos.javanicus.spp = pixel_list[["leptoptilos.javanicus"]]))
jev.xx
AIC(jev.xx)

jev.xxi <- ppm(jev.q, ~ mycteria.cinerea.spp,
                covariates = list(mycteria.cinerea.spp = pixel_list[["mycteria.cinerea"]]))
jev.xxi
AIC(jev.xxi)

jev.xxii <- ppm(jev.q, ~ mycteria.leucocephala.spp,
                covariates = list(mycteria.leucocephala.spp = pixel_list[["mycteria.leucocephala"]]))
jev.xxii
AIC(jev.xxii)

jev.xxiii <- ppm(jev.q, ~ ciconia.stormi.spp,
                covariates = list(ciconia.stormi.spp = pixel_list[["ciconia.stormi"]]))
jev.xxiii
AIC(jev.xxiii)
