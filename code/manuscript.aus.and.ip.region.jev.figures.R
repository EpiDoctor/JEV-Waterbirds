library(viridis)

##########
# Figures
##########

#Figure 1

tiff("Figure 1 revised 21.02.2025.tiff", width = 12, height = 8, units = "in", res = 600, compression = c("lzw"))
par(mfrow=c(1,2))

plot(predictors$bird.rich, 
     main = "Waterbird richness", 
     cex = 1.25, 
     col = rev(viridis(200, option = "A")))
plot(jev_sf[0], 
     pch = 19, 
     col = "blue", 
     cex = .35, 
     add = T)

plot(predictors$trait.hwi, 
     main = "PWM Hand-wing index", 
     cex = 1.25, 
     col = rev(viridis(200, option = "A")))
plot(jev_sf[0], 
     pch = 19, 
     col = "skyblue", 
     cex = .35, 
     add = T)

dev.off()

# figure 2

tiff("Figure 2 revised 21.02.2025.tiff", width = 10, height = 4, units = "in", res = 600, compression = c("lzw"))
par(mfrow=c(1,3))

plot(predictors$bird.rich, 
     main = "Waterbird richness", 
     cex = 1.25, 
     col = rev(viridis(200, option = "A")))
plot(jev_sf[0], 
     pch = 19, 
     col = "skyblue", 
     cex = .35, 
     add = T)

plot(predictors$trait.hwi, 
     main = "PWM Hand-wing index", 
     cex = 1.25, 
     col = rev(viridis(200, option = "A")))
plot(jev_sf[0], 
     pch = 19, 
     col = "skyblue", 
     cex = .35, 
     add = T)

plot(predictors$trait.diet.plant, 
     main = "PWM % plant diet", 
     cex = 1.25, 
     col = rev(viridis(200, option = "A")))
plot(jev_sf[0], 
     pch = 19, 
     col = "blue", 
     cex = .35, 
     add = T)

dev.off()

# figure 3

tiff("Figure 3 revised 21.02.2025.tiff", width = 12, height = 4, units = "in", res = 600, compression = c("lzw"))
par(mfrow=c(1,3), omi = c(.1, .1, .1, .2), mai = c(.25, 0.25, 0.25, 0.25))

plot(pFull.lower.cat, 
     main = "Lower 95% C.L.",  
     col = rev(viridis(10, option = "A")), 
     cex.axis = 1.5, 
     cex.main = 1.5,
     legend = F) 

plot(pFull.ras.cat, 
     main = "Species pool landscape targets", 
     col = rev(viridis(10, option = "A")), 
     cex.axis = 1.5, 
     cex.main = 1.5,
     legend = F)

plot(pFull.upper.cat, 
     main = "Upper 95% C.L.", 
     col = rev(viridis(10, option = "A")), 
     cex.axis = 1.5, 
     cex.main = 1.5) 

dev.off()

# figure 4

tiff("Figure 4 revised 21.02.2025.tiff", width = 12, height = 12, units = "in", res = 600, compression = c("lzw"))
par(mfrow=c(1,1), omi = c(.1, .1, .1, .2))

plotBeta(model.1, 
         post = post.beta.1, 
         param = "Sign", 
         plotTree = T, 
         supportLevel = 0.95, 
         split = 0.5, 
         spNamesNumbers = c(T, F),
         cex = c(0.95, 0.85, 0.8),
         covOrder = "Vector",
         covVector = c(2, 3, 4, 5, 6),
         covNamesNumbers = c(F, T),
         colors = colorRampPalette(c("#400F73FF", "#FCFDBFFF", "#E85362FF")),
         newplot = F)

dev.off()

# figure 5

tiff("Figure 5 revised 21.02.2025.tiff", width = 12, height = 12, units = "in", res = 600, compression = c("lzw"))
par(mfrow=c(1,1), omi = c(.1, .1, .1, .2))

plotBeta(model.2, 
         post = post.beta.2, 
         param = "Sign", 
         plotTree = T, 
         supportLevel = 0.95, 
         split = 0.5, 
         spNamesNumbers = c(T, F),
         cex = c(0.95, 0.85, 0.8),
         covNamesNumbers = c(F, T),
         colors = colorRampPalette(c("#400F73FF", "#FCFDBFFF", "#E85362FF")),
         newplot = T
)

dev.off()

# figure 6

tiff("Figure 6 revised 21.02.2025.tiff", width = 12, height = 8, units = "in", res = 600, compression = c("lzw"))
par(mfrow=c(1, 2), omi = c(.1, .1, .1, .2))

plotGamma(model.1, 
          post = post.gamma.1, 
          param = "Sign",
          supportLevel = 0.95,
          cex = c(0.7, 0.75, 0.8),
          covOrder = "Vector",
          covVector = c(2, 3, 4, 5, 6),
          trOrder = "Vector",
          trVector = c(2, 3, 4, 5, 6),
          trNamesNumbers = c(T, F),
          covNamesNumbers = c(F, T),
          colors = colorRampPalette(c("#400F73FF", "#FCFDBFFF", "#E85362FF")), # #07071DFF
          newplot = F
)

plotGamma(model.2, 
          post = post.gamma.2, 
          param = "Sign",
          supportLevel = 0.95,
          cex = c(0.7, 0.75, 0.8),
          trOrder = "Vector",
          trVector = c(2, 3, 4, 5, 6),
          trNamesNumbers = c(T, F),
          covNamesNumbers = c(F, T),
          colors = colorRampPalette(c("#400F73FF", "#FCFDBFFF", "#E85362FF")), # #07071DFF
          newplot = F
)

dev.off()

# figure 7

tree <- keep.tip(consensus.tree, intersect(consensus.tree$tip.label, colnames(y.abund.jsdm))) 
trait.hwi <- trait.data.jsdm %>%
  dplyr::select(HWI)
head(trait.hwi, 5)
hwi.cont.map <- contMap(tree, trait.hwi, fsize = .5)
hwi.cont.map <- setMap(hwi.cont.map, rev(viridis(7, option = "A")))

tiff("Figure 7 revised 21.02.2025.tiff", width = 12, height = 12, units = "in", res = 600, compression = c("lzw"))

plot(hwi.cont.map, fsize = c(0.8,0.8), leg.txt="Hand wing index", rev(viridis(7, option = "A")))

dev.off()

# Figure S1

tiff("Figure S1.tiff", width = 12, height = 8, units = "in", res = 600, compression = c("lzw"))
par(mfrow=c(2,3))

plot(predictors$waterway.dist, 
     main = "Distance to waterways (km)", 
     cex = 1.25, 
     col = viridis(200, option = "D"))
plot(jev_sf[0], 
     pch = 19, 
     col = "salmon", 
     cex = .5, 
     add = T)

plot(flow.cat, 
     main = "Surface water flow accumulation (deciles)", 
     cex = 1.25, 
     col = rev(viridis(200, option = "D")))
plot(jev_sf[0], 
     pch = 19, 
     col = "salmon", 
     cex = .5, 
     add = T)

plot(predictors$value, 
     main = "Temporary wetlands (%)", 
     cex = 1.25, 
     col = rev(viridis(200, option = "D")))
plot(jev_sf[0], 
     pch = 19, 
     col = "salmon", 
     cex = .5, 
     add = T)

plot(predictors$`iucn_habitatclassification_fraction_lvl2__1401_Arable Land__ver003`, 
     main = "Cultivated land (%)", 
     cex = 1.25, 
     col = rev(viridis(200, option = "A")))
plot(jev_sf[0], 
     pch = 19, 
     col = "skyblue", 
     cex = .5, 
     add = T)

plot(predictors$class_circle_mn_400, 
     main = "Mean grassland RCC", 
     cex = 1.25, 
     col = rev(viridis(200, option = "A")))
plot(jev_sf[0], 
     pch = 19, 
     col = "skyblue", 
     cex = .5, 
     add = T)

plot(predictors$ardeid.richness, 
     main = "Ardeid richness", 
     cex = 1.25, 
     col = rev(viridis(200, option = "A")))
plot(jev_sf[0], 
     pch = 19, 
     col = "blue", 
     cex = .5, 
     add = T)

dev.off()

# Figure S2

tiff("Figure S2 revised b.tiff", width = 10, height = 4, units = "in", res = 600, compression = c("lzw"))
par(mfrow=c(1,3))

plot(density(traits.imp$ximp$mean_HWI), 
     lty = 3, lwd = 2, col = "red",
     ylim = c(0,0.145),
     main = "Hand-wing index", 
     xlab = "(Wing length - Secondary length/Wing length)*100")
lines(density(traits$mean_HWI, na.rm = T), lwd = 2)
plot(density(traits.imp$ximp$Egg_mass), 
     lty = 3, lwd = 2, col = "red", 
     ylim = c(0,0.145),
     main = "Egg mass", 
     xlab = "grams")
lines(density(traits$Egg_mass, na.rm = T), lwd = 2)
plot(density(traits.imp$ximp$mean_density), 
     lty = 3, lwd = 2, col = "red", 
     ylim = c(0,0.145),
     main = "Population density", 
     xlab = expression("birds/km"^2))
lines(density(traits$mean_density, na.rm = T), lwd = 2)

dev.off()