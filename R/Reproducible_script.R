

###################################################################################
####################################################################################
## Reproducible script for Cortell-Nicolau, A., Kandler, A., 'name of the paper', Journal of Archaeological Method and Theory
####################################################################################
####################################################################################

## Load library and data
#devtools::install_github("acortell3/Simumorph") 
library(Simumorph)

## Load necessary objects (produced with morphospace.R)
geo_out <- readRDS("../Utilities/geo_out.rds") ## Observed shapes
morpho_pars <- readRDS("../Utilities/morpho_pars.rds") ## Observed parameters
amp_pha_mat <- readRDS("../Utilities/amp_pha_mat.rds") ## Amplitude and phase matrix
amp_pha_cov <- readRDS("../Utilities/amp_pha_cov.rds") ## Covariance matrix

## Load utilities
sims <- 100
npts <- 120
	
mat <- matrix(c(seq(1,100), rep(101,20)), nrow = 12, byrow = T) ## For plots
seed <- 123

###################################################################################
## Figure 3. Stability covariance matrix
png("../Figures/Fig_3.png", res = 50, height = 800, width = 800)
valid_morphospace(x = amp_pha_mat, k = 4, rnd = 100, msize = 2, lsize = 1.6, trans = 0.7)
dev.off()
###################################################################################

###################################################################################
## Figure 4a. Case study 1
set.seed(seed)

AtoAa <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = 1, method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 0.2, e = 0.05, f = 100)

## Plot
png("../Figures/Fig_4a.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(AtoAa$Shapes[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
	polygon(AtoAa$Shapes[[i]], col = adjustcolor("seagreen", alpha = 0.85))
}
plot(unlist(AtoAa$P.distances), type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "seagreen3", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()
###################################################################################

###################################################################################
## Figure 4b. Case study 1 (dynamic)
set.seed(seed)

dyn_e <- data.frame("time" = c(1,10,40,70,85,95),
		    "e" = c(0.02,0.03,0.25,0.2,0.15,0.1))

AtoAb <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = 1, method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 0.2, f = 100, dynamic_e = dyn_e)

## Plot
png("../Figures/SI_figures/Fig_4b.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(AtoAb$Shapes[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
	polygon(AtoAb$Shapes[[i]], col = adjustcolor("seagreen", alpha = 0.85))
}
plot(unlist(AtoAb$P.distances), type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "seagreen3", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()
###################################################################################

###################################################################################
## Figure 5. Case study 2a
set.seed(seed)

AtoBa <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = which(rownames(amp_pha_mat)=="G18_m_G18"), method = "AtoB", sim = sims, npts = npts, only.shapes = F, a = 0.5, e = 0.5, f = 100, max.attempts = 500, speedAtoB = 0)

## Plot
png("../Figures/Fig_5a.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(AtoBa$Shapes[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
	polygon(AtoBa$Shapes[[i]], col = adjustcolor("sienna4", alpha = 0.85))
}
plot(unlist(AtoBa$P.distances), type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "sienna4", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()
###################################################################################

###################################################################################
## Figure 5. Case study 2b (higher dist)
set.seed(seed)

AtoBb <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = which(rownames(amp_pha_mat)=="G18_m_G18"), method = "AtoB", sim = 500, npts = npts, only.shapes = F, a = 0.5, e = 0.5, f = 100, max.attempts = 500, speedAtoB = 0.0035)

samp <- seq(1,500,5)

## Plot
png("../Figures/SI_figures/Fig_5b.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(AtoBb$Shapes[[samp[i]]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",samp[i]))
	polygon(AtoBb$Shapes[[samp[i]]], col = adjustcolor("sienna4", alpha = 0.85))
}
plot(unlist(AtoBb$P.distances), type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "sienna4", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()
###################################################################################

###################################################################################
## Figure 5. Case study 2c (lower dist)
set.seed(seed)

AtoBc <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = which(rownames(amp_pha_mat)=="G18_m_G18"), method = "AtoB", sim = sims, npts = npts, only.shapes = F, a = 0.5, e = 0.5, f = 100, max.attempts = 500, speedAtoB = -0.002)

## Plot
png("../Figures/SI_figures/Fig_5c.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(AtoBc$Shapes[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
	polygon(AtoBc$Shapes[[i]], col = adjustcolor("sienna4", alpha = 0.85))
}
plot(unlist(AtoBc$P.distances), type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "sienna4", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()
###################################################################################
###################################################################################
## Figure 6. Case study 3
set.seed(seed)

### Let's try with the full morphospace
targets_multi <- c(1:nrow(amp_pha_mat))
target_names <- rownames(amp_pha_mat)
c_a <- 3
c_f <- 50


AtoMult_G1 <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G1.1_m_G1.1"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)

## Create vector with minimum procrustes distances and their 
pdist_G1 <- data.frame("target" = character(),
		       "distance" = numeric())

for (i in 1:sims){
	pdist_G1[i,1] <- target_names[which.min(AtoMult_G1$P.distances[,i])]
	pdist_G1[i,2] <- min(AtoMult_G1$P.distances[,i])
}


## Plot
png("../Figures/Fig_6a.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(AtoMult_G1$Shapes[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i," / S = ", pdist_G1[i,1]))
	polygon(AtoMult_G1$Shapes[[i]], col = adjustcolor("slategray3", alpha = 0.85))
}
plot(pdist_G1$distance, type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "slategrey", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()

## Go for G2
AtoMult_G2 <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G2_m_G2"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)

## Create vector with minimum procrustes distances and their 
pdist_G2 <- data.frame("target" = character(),
		       "distance" = numeric())

for (i in 1:sims){
	pdist_G2[i,1] <- target_names[which.min(AtoMult_G2$P.distances[,i])]
	pdist_G2[i,2] <- min(AtoMult_G2$P.distances[,i])
}


## Plot
png("../Figures/SI_figures/Fig_6b.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(AtoMult_G2$Shapes[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i," / S = ", pdist_G2[i,1]))
	polygon(AtoMult_G2$Shapes[[i]], col = adjustcolor("slategray3", alpha = 0.85))
}
plot(pdist_G2$distance, type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "slategrey", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()


### Now G9
AtoMult_G9 <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)

## Create vector with minimum procrustes distances and their 
pdist_G9 <- data.frame("target" = character(),
		       "distance" = numeric())

for (i in 1:sims){
	pdist_G9[i,1] <- target_names[which.min(AtoMult_G9$P.distances[,i])]
	pdist_G9[i,2] <- min(AtoMult_G9$P.distances[,i])
}


## Plot
png("../Figures/SI_figures/Fig_6c.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(AtoMult_G9$Shapes[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i," / S = ", pdist_G9[i,1]))
	polygon(AtoMult_G9$Shapes[[i]], col = adjustcolor("slategray3", alpha = 0.85))
}
plot(pdist_G9$distance, type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "slategrey", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()


### And G18
AtoMult_G18 <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G18_m_G18"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)

## Create vector with minimum procrustes distances and their 
pdist_G18 <- data.frame("target" = character(),
	        	"distance" = numeric())

for (i in 1:sims){
	pdist_G18[i,1] <- target_names[which.min(AtoMult_G18$P.distances[,i])]
	pdist_G18[i,2] <- min(AtoMult_G18$P.distances[,i])
}


## Plot
png("../Figures/SI_figures/Fig_6d.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(AtoMult_G18$Shapes[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i," / S = ", pdist_G18[i,1]))
	polygon(AtoMult_G18$Shapes[[i]], col = adjustcolor("slategray3", alpha = 0.85))
}
plot(pdist_G18$distance, type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "slategrey", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()
###################################################################################

###################################################################################
## Figure 7. Case study 4
###################################################################################
set.seed(seed)

Free <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = nrow(amp_pha_mat), method = "Free", sim = sims, npts = npts, only.shapes = F, a = 0.2, e = 0.05, f = 100, max.attempts = 500)

## Plot
png("../Figures/Fig_7.png", res = 50, height = 1500, width = 1500)
layout(mat)
for (i in 1:sims){
	plot(Free$Shapes[[i]], type = "l", lwd = 1.5, cex.main = 2, xlab = "", ylab = "", bty = "n", main = paste0("t = ",i))
	polygon(Free$Shapes[[i]], col = adjustcolor("lightsalmon4", alpha = 0.85))
}
plot(unlist(Free$P.distances), type = "l", lwd = 3, cex.main = 2, cex.lab = 2.2, mgp = c(1.8,0.7,0), col = "lightsalmon4", main = "Procrustes distances", xlab = "Time", ylab = "Dist", bty = "n")
dev.off()






