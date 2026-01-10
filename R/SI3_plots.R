
#####################################################################################
######## The following SI script includes the plots for the additional SI figures for <paper name>. It explores the variation in different ways of the four implemented methods
#####################################################################################

################ EXTENDED SIMULATIONS AtoA
## 100 most distant shapes to initial shapes. Need Procrustes distances, time of max divergence, and shapes themselves
SI_AtoA <- readRDS("../SI3_results/SI_AtoA_res.rds")
AtoA_mdi <- SI_AtoA[[1]]
AtoA_md <- SI_AtoA[[2]]
AtoA_shapes <- SI_AtoA[[3]]

## Prepare objects for plots
## df
dens_AtoA_mdi <- data.frame("x" = density(AtoA_mdi)$x,
		 	    "y" = density(AtoA_mdi)$y)

dens_AtoA_md <- data.frame("x" = density(AtoA_md)$x,
		 	   "y" = density(AtoA_md)$y)

## cis
AtoA_mdi_quant <- quantile(AtoA_mdi, probs = c(0.025, 0.5, 0.975))
AtoA_md_quant <- quantile(AtoA_md, probs = c(0.025, 0.5, 0.975))

dens_AtoA_mdi_ci <- dens_AtoA_mdi[dens_AtoA_mdi$x <= AtoA_mdi_quant[3] & dens_AtoA_mdi$x >= AtoA_mdi_quant[1],] 
dens_AtoA_md_ci <- dens_AtoA_md[dens_AtoA_md$x <= AtoA_md_quant[3] & dens_AtoA_md$x >= AtoA_md_quant[1],] 


dens_AtoA_mdi_ci$y[1] <- 0
dens_AtoA_mdi_ci$y[nrow(dens_AtoA_mdi_ci)] <- 0

dens_AtoA_md_ci$y[1] <- 0
dens_AtoA_md_ci$y[nrow(dens_AtoA_md_ci)] <- 0

png("../SI_figures/SI3/Fig_1_SI3.png", res = 100, height = 1000, width = 1000)
par(mfrow = c(1,2))
plot(x=dens_AtoA_mdi$x, y=dens_AtoA_mdi$y, xlab = "Timestep of max divergence", col = "blue", lwd = 1.5, main = "Max divergence (AtoA)", type = "l", ylab = "value")
polygon(x=dens_AtoA_mdi$x, y=dens_AtoA_mdi$y, col = "lightblue1", border = NA)
polygon(x=dens_AtoA_mdi_ci$x, y=dens_AtoA_mdi_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_AtoA_mdi$x, y=dens_AtoA_mdi$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(AtoA_mdi), lty = 2, col = "navyblue")
abline(v = median(AtoA_mdi), lty = 2, col = "red")
plot(x=dens_AtoA_md$x, y = dens_AtoA_md$y, xlab = "Max Procrustes distance", col = "blue", lwd = 1.5, main = "Max Procrustes distance (AtoA)", type = "l", ylab = "value")
polygon(x=dens_AtoA_md$x, y=dens_AtoA_md$y, col = "lightblue1", border = NA)
polygon(x=dens_AtoA_md_ci$x, y=dens_AtoA_md_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_AtoA_md$x, y=dens_AtoA_md$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(AtoA_md), lty = 2, col = "navyblue")
abline(v = median(AtoA_md), lty = 2, col = "red")

dev.off()

################ EXTENDED SIMULATIONS AtoB
## When do we reach maximum convergence? Need Procrustes distances, time of max convergence, and shapes themselves

SI_AtoB <- readRDS("../SI3_results/SI_AtoB_res.rds")
AtoB_mdi <- SI_AtoB[[1]]
AtoB_md <- SI_AtoB[[2]]
AtoB_shapes <- SI_AtoB[[3]]

## Prepare objects for plots
## df
dens_AtoB_mdi <- data.frame("x" = density(AtoB_mdi)$x,
		 	    "y" = density(AtoB_mdi)$y)

dens_AtoB_md <- data.frame("x" = density(AtoB_md)$x,
		 	   "y" = density(AtoB_md)$y)

## cis
AtoB_mdi_quant <- quantile(AtoB_mdi, probs = c(0.025, 0.5, 0.975))
AtoB_md_quant <- quantile(AtoB_md, probs = c(0.025, 0.5, 0.975))

dens_AtoB_mdi_ci <- dens_AtoB_mdi[dens_AtoB_mdi$x <= AtoB_mdi_quant[3] & dens_AtoB_mdi$x >= AtoB_mdi_quant[1],] 
dens_AtoB_md_ci <- dens_AtoB_md[dens_AtoB_md$x <= AtoB_md_quant[3] & dens_AtoB_md$x >= AtoB_md_quant[1],] 


dens_AtoB_mdi_ci$y[1] <- 0
dens_AtoB_mdi_ci$y[nrow(dens_AtoB_mdi_ci)] <- 0

dens_AtoB_md_ci$y[1] <- 0
dens_AtoB_md_ci$y[nrow(dens_AtoB_md_ci)] <- 0

png("../SI_figures/SI3/Fig_2_SI3.png", res = 100, height = 1000, width = 1000)
par(mfrow = c(1,2))
plot(x=dens_AtoB_mdi$x, y=dens_AtoB_mdi$y, xlab = "Timestep of max convergence", col = "blue", lwd = 1.5, main = "Max convergence (AtoB)", type = "l", ylab = "value")
polygon(x=dens_AtoB_mdi$x, y=dens_AtoB_mdi$y, col = "lightblue1", border = NA)
polygon(x=dens_AtoB_mdi_ci$x, y=dens_AtoB_mdi_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_AtoB_mdi$x, y=dens_AtoB_mdi$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(AtoB_mdi), lty = 2, col = "navyblue")
abline(v = median(AtoB_mdi), lty = 2, col = "red")
plot(x=dens_AtoB_md$x, y = dens_AtoB_md$y, xlab = "Max convergence distance", col = "blue", lwd = 1.5, main = "Min Procrustes distance (AtoB)", type = "l", ylab = "value")
polygon(x=dens_AtoB_md$x, y=dens_AtoB_md$y, col = "lightblue1", border = NA)
polygon(x=dens_AtoB_md_ci$x, y=dens_AtoB_md_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_AtoB_md$x, y=dens_AtoB_md$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(AtoB_md), lty = 2, col = "navyblue")
abline(v = median(AtoB_md), lty = 2, col = "red")

dev.off()

################ EXTENDED SIMULATIONS AtoMult

library(wordcloud)

## When is it a different shape? How many shapes does it go through? And to which shapes does it change? Comparison with and without initial shape in target pool

## G1a

SI_AtoMulta_G1 <- readRDS("../SI3_results/SI_AtoMulta_G1_res.rds")
Multa_G1_fd <- sapply(SI_AtoMulta_G1, function(x) min(which(x != x[1]))) ## First element different from 1
Multa_G1_dc <- sapply(SI_AtoMulta_G1, function(x) length(unique(x))) ## How many different types 
Multa_G1_vc <- as.data.frame(table(unlist(SI_AtoMulta_G1))) ## Count of how many appearances of each
colnames(Multa_G1_vc) <- c("Type", "Appearances")

## Prepare for density plot of first element different from 1
## df
dens_Multa_G1_fd <- data.frame("x" = density(Multa_G1_fd)$x,
		 	       "y" = density(Multa_G1_fd)$y)

## cis
Multa_G1_fd_quant <- quantile(Multa_G1_fd, probs = c(0.025, 0.5, 0.975))

dens_Multa_G1_fd_ci <- dens_Multa_G1_fd[dens_Multa_G1_fd$x <= Multa_G1_fd_quant[3] & dens_Multa_G1_fd$x >= Multa_G1_fd_quant[1],] 

dens_Multa_G1_fd_ci$y[1] <- 0
dens_Multa_G1_fd_ci$y[nrow(dens_Multa_G1_fd_ci)] <- 0

## Prepare index for barplot. Select five with most appearances
xlab_n <- rep("",nrow(Multa_G1_vc))
sel_types <- order(Multa_G1_vc$Appearances, decreasing = TRUE)[1:5]
xlab_n[sel_types] <- as.character(Multa_G1_vc$Type[sel_types])

png("../SI_figures/SI3/Fig_3_SI3.png", res = 140, height = 1200, width = 1200)
layout(matrix(c(1,2,3,3),ncol = 2, byrow = T))
plot(x=dens_Multa_G1_fd$x, y=dens_Multa_G1_fd$y, xlab = "Timestep to first different type", col = "blue", lwd = 1.5, main = "First different type (AtoMulta_G1)", type = "l", ylab = "value")
polygon(x=dens_Multa_G1_fd$x, y=dens_Multa_G1_fd$y, col = "lightblue1", border = NA)
polygon(x=dens_Multa_G1_fd_ci$x, y=dens_Multa_G1_fd_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Multa_G1_fd$x, y=dens_Multa_G1_fd$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(Multa_G1_fd), lty = 2, col = "navyblue")
abline(v = median(Multa_G1_fd), lty = 2, col = "red")

wordcloud(words = as.character(Multa_G1_vc$Type), freq = Multa_G1_vc$Appearances, scale = c(2,0.5), colors = brewer.pal(8, "Dark2"), random.order = T)
title("Type frequency")

barplot(Multa_G1_vc$Appearances, names.arg = xlab_n, las = 2, col = "gold4", main = "Types present in the simulation", cex.names = 0.6)
dev.off()

## G2a

SI_AtoMulta_G2 <- readRDS("../SI3_results/SI_AtoMulta_G2_res.rds")
Multa_G2_fd <- sapply(SI_AtoMulta_G2, function(x) min(which(x != x[1]))) ## First element different from 1
Multa_G2_dc <- sapply(SI_AtoMulta_G2, function(x) length(unique(x))) ## How many different types 
Multa_G2_vc <- as.data.frame(table(unlist(SI_AtoMulta_G2))) ## Count of how many appearances of each
colnames(Multa_G2_vc) <- c("Type", "Appearances")

## Prepare for density plot of first element different from 1
## df
dens_Multa_G2_fd <- data.frame("x" = density(Multa_G2_fd)$x,
		 	       "y" = density(Multa_G2_fd)$y)

## cis
Multa_G2_fd_quant <- quantile(Multa_G2_fd, probs = c(0.025, 0.5, 0.975))

dens_Multa_G2_fd_ci <- dens_Multa_G2_fd[dens_Multa_G2_fd$x <= Multa_G2_fd_quant[3] & dens_Multa_G2_fd$x >= Multa_G2_fd_quant[1],] 

dens_Multa_G2_fd_ci$y[1] <- 0
dens_Multa_G2_fd_ci$y[nrow(dens_Multa_G2_fd_ci)] <- 0

## Prepare index for barplot. Select five with most appearances
xlab_n <- rep("",nrow(Multa_G2_vc))
sel_types <- order(Multa_G2_vc$Appearances, decreasing = TRUE)[1:5]
xlab_n[sel_types] <- as.character(Multa_G2_vc$Type[sel_types])

png("../SI_figures/SI3/Fig_4_SI3.png", res = 140, height = 1200, width = 1200)
layout(matrix(c(1,2,3,3),ncol = 2, byrow = T))
plot(x=dens_Multa_G2_fd$x, y=dens_Multa_G2_fd$y, xlab = "Timestep to first different type", col = "blue", lwd = 1.5, main = "First different type (AtoMulta G2)", type = "l", ylab = "value")
polygon(x=dens_Multa_G2_fd$x, y=dens_Multa_G2_fd$y, col = "lightblue1", border = NA)
polygon(x=dens_Multa_G2_fd_ci$x, y=dens_Multa_G2_fd_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Multa_G2_fd$x, y=dens_Multa_G2_fd$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(Multa_G2_fd), lty = 2, col = "navyblue")
abline(v = median(Multa_G2_fd), lty = 2, col = "red")

wordcloud(words = as.character(Multa_G2_vc$Type), freq = Multa_G2_vc$Appearances, scale = c(2,0.5), colors = brewer.pal(8, "Dark2"), random.order = T)
title("Type frequency")

barplot(Multa_G2_vc$Appearances, names.arg = xlab_n, las = 2, col = "gold4", main = "Types present in the simulation", cex.names = 0.6)
dev.off()

## G9a

SI_AtoMulta_G9 <- readRDS("../SI3_results/SI_AtoMulta_G9_res.rds")
Multa_G9_fd <- sapply(SI_AtoMulta_G9, function(x) min(which(x != x[1]))) ## First element different from 1
Multa_G9_fd[which(Multa_G9_fd == Inf)] <- NA ## For the cases where it doesn't change 
Multa_G9_fd <- na.omit(Multa_G9_fd)
Multa_G9_dc <- sapply(SI_AtoMulta_G9, function(x) length(unique(x))) ## How many different types 
Multa_G9_vc <- as.data.frame(table(unlist(SI_AtoMulta_G9))) ## Count of how many appearances of each
colnames(Multa_G9_vc) <- c("Type", "Appearances")

## Prepare for density plot of first element different from 1
## df
dens_Multa_G9_fd <- data.frame("x" = density(Multa_G9_fd)$x,
		 	       "y" = density(Multa_G9_fd)$y)

## cis
Multa_G9_fd_quant <- quantile(Multa_G9_fd, probs = c(0.025, 0.5, 0.975))

dens_Multa_G9_fd_ci <- dens_Multa_G9_fd[dens_Multa_G9_fd$x <= Multa_G9_fd_quant[3] & dens_Multa_G9_fd$x >= Multa_G9_fd_quant[1],] 

dens_Multa_G9_fd_ci$y[1] <- 0
dens_Multa_G9_fd_ci$y[nrow(dens_Multa_G9_fd_ci)] <- 0

## Prepare index for barplot. Select five with most appearances
xlab_n <- rep("",nrow(Multa_G9_vc))
sel_types <- order(Multa_G9_vc$Appearances, decreasing = TRUE)[1:5]
xlab_n[sel_types] <- as.character(Multa_G9_vc$Type[sel_types])

png("../SI_figures/SI3/Fig_5_SI3.png", res = 140, height = 1200, width = 1200)
layout(matrix(c(1,2,3,3),ncol = 2, byrow = T))
plot(x=dens_Multa_G9_fd$x, y=dens_Multa_G9_fd$y, xlab = "Timestep to first different type", col = "blue", lwd = 1.5, main = "First different type (AtoMulta G9)", type = "l", ylab = "value")
polygon(x=dens_Multa_G9_fd$x, y=dens_Multa_G9_fd$y, col = "lightblue1", border = NA)
polygon(x=dens_Multa_G9_fd_ci$x, y=dens_Multa_G9_fd_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Multa_G9_fd$x, y=dens_Multa_G9_fd$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(Multa_G9_fd), lty = 2, col = "navyblue")
abline(v = median(Multa_G9_fd), lty = 2, col = "red")
text(x = 70, y = 0.08, paste0(100-length(Multa_G9_fd)," sims did not change"))

wordcloud(words = as.character(Multa_G9_vc$Type), freq = Multa_G9_vc$Appearances, scale = c(2,0.5), colors = brewer.pal(8, "Dark2"), random.order = T)
title("Type frequency")

barplot(Multa_G9_vc$Appearances, names.arg = xlab_n, las = 2, col = "gold4", main = "Types present in the simulation", cex.names = 0.6)
dev.off()

## G18a

#### Comment the code and not use it bc it stays within one single shape for several simulations

#SI_AtoMulta_G18 <- readRDS("../SI3_results/SI_AtoMulta_G18_res.rds")
#Multa_G18_fd <- sapply(SI_AtoMulta_G18, function(x) min(which(x != x[1]))) ## First element different from 1
#Multa_G18_fd[which(Multa_G18_fd == Inf)] <- NA ## For the cases where it doesn't change 
#Multa_G18_fd <- na.omit(Multa_G18_fd)
#Multa_G18_dc <- sapply(SI_AtoMulta_G18, function(x) length(unique(x))) ## How many different types 
#Multa_G18_vc <- as.data.frame(table(unlist(SI_AtoMulta_G18))) ## Count of how many appearances of each
#colnames(Multa_G18_vc) <- c("Type", "Appearances")

## Prepare for density plot of first element different from 1
## df
#dens_Multa_G18_fd <- data.frame("x" = density(Multa_G18_fd)$x,
#		 	       "y" = density(Multa_G18_fd)$y)

## cis
#Multa_G18_fd_quant <- quantile(Multa_G18_fd, probs = c(0.025, 0.5, 0.975))

#dens_Multa_G18_fd_ci <- dens_Multa_G18_fd[dens_Multa_G18_fd$x <= Multa_G18_fd_quant[3] & dens_Multa_G18_fd$x >= Multa_G18_fd_quant[1],] 

#dens_Multa_G18_fd_ci$y[1] <- 0
#dens_Multa_G18_fd_ci$y[nrow(dens_Multa_G18_fd_ci)] <- 0

## Prepare index for barplot. Select five with most appearances
#xlab_n <- rep("",nrow(Multa_G18_vc))
#sel_types <- order(Multa_G18_vc$Appearances, decreasing = TRUE)[1:5]
#xlab_n[sel_types] <- as.character(Multa_G18_vc$Type[sel_types])

#layout(matrix(c(1,2,3,3),ncol = 2, byrow = T))
#plot(x=dens_Multa_G18_fd$x, y=dens_Multa_G18_fd$y, xlab = "Timestep to first different type", col = "blue", lwd = 1.5, main = "First different type", type = "l", ylab = "value")
#polygon(x=dens_Multa_G18_fd$x, y=dens_Multa_G18_fd$y, col = "lightblue1", border = NA)
#polygon(x=dens_Multa_G18_fd_ci$x, y=dens_Multa_G18_fd_ci$y, col = "lightblue3", border = NA)
#polygon(x=dens_Multa_G18_fd$x, y=dens_Multa_G18_fd$y, col = adjustcolor("blue", alpha = 0), border = "black")
#lines(x = rep(mean(dens_Multa_G18_fd$x),2), y = c(0,dens_Multa_G18_fd$y[which.min(abs(dens_Multa_G18_fd$x - mean(dens_Multa_G18_fd$x)))]), col = "navyblue", lty = 1.5)
#text(x = 70, y = 0.08, paste0(100-length(Multa_G18_fd)," sims did not change"))

#par(mar = c(1,1,2,1))
#wordcloud(words = as.character(Multa_G18_vc$Type), freq = Multa_G18_vc$Appearances, scale = c(2,0.5), colors = brewer.pal(8, "Dark2"), random.order = F)
#title("Type frequency")

#par(mar = c(5.1,4.1,4.1,2.1))
#barplot(Multa_G18_vc$Appearances, names.arg = xlab_n, las = 2, col = "gold4", main = "Types present in the simulation", cex.names = 0.6)


######### Now without the initial shape within the target pool


## G1b

SI_AtoMultb_G1 <- readRDS("../SI3_results/SI_AtoMultb_G1_res.rds")
Multb_G1_fd <- sapply(SI_AtoMultb_G1, function(x) min(which(x != x[1]))) ## First element different from 1
Multb_G1_dc <- sapply(SI_AtoMultb_G1, function(x) length(unique(x))) ## How many different types 
Multb_G1_vc <- as.data.frame(table(unlist(SI_AtoMultb_G1))) ## Count of how many appearances of each
colnames(Multb_G1_vc) <- c("Type", "Appearances")

## Prepare for density plot of first element different from 1
## df
dens_Multb_G1_fd <- data.frame("x" = density(Multb_G1_fd)$x,
		 	       "y" = density(Multb_G1_fd)$y)

## cis
Multb_G1_fd_quant <- quantile(Multb_G1_fd, probs = c(0.025, 0.5, 0.975))

dens_Multb_G1_fd_ci <- dens_Multb_G1_fd[dens_Multb_G1_fd$x <= Multb_G1_fd_quant[3] & dens_Multb_G1_fd$x >= Multb_G1_fd_quant[1],] 

dens_Multb_G1_fd_ci$y[1] <- 0
dens_Multb_G1_fd_ci$y[nrow(dens_Multb_G1_fd_ci)] <- 0

## Prepare index for barplot. Select five with most appearances
xlab_n <- rep("",nrow(Multb_G1_vc))
sel_types <- order(Multb_G1_vc$Appearances, decreasing = TRUE)[1:5]
xlab_n[sel_types] <- as.character(Multb_G1_vc$Type[sel_types])

png("../SI_figures/SI3/Fig_6_SI3.png", res = 140, height = 1200, width = 1200)
layout(matrix(c(1,2,3,3),ncol = 2, byrow = T))
plot(x=dens_Multb_G1_fd$x, y=dens_Multb_G1_fd$y, xlab = "Timestep to first different type", col = "blue", lwd = 1.5, main = "First different type (AtoMultb G1)", type = "l", ylab = "value")
polygon(x=dens_Multb_G1_fd$x, y=dens_Multb_G1_fd$y, col = "lightblue1", border = NA)
polygon(x=dens_Multb_G1_fd_ci$x, y=dens_Multb_G1_fd_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Multb_G1_fd$x, y=dens_Multb_G1_fd$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(Multb_G1_fd), lty = 2, col = "navyblue")
abline(v = median(Multb_G1_fd), lty = 2, col = "red")

wordcloud(words = as.character(Multb_G1_vc$Type), freq = Multb_G1_vc$Appearances, scale = c(2,0.5), colors = brewer.pal(8, "Dark2"), random.order = T)
title("Type frequency")

barplot(Multb_G1_vc$Appearances, names.arg = xlab_n, las = 2, col = "gold4", main = "Types present in the simulation", cex.names = 0.6)
dev.off()

## G2b

SI_AtoMultb_G2 <- readRDS("../SI3_results/SI_AtoMultb_G2_res.rds")
Multb_G2_fd <- sapply(SI_AtoMultb_G2, function(x) min(which(x != x[1]))) ## First element different from 1
Multb_G2_dc <- sapply(SI_AtoMultb_G2, function(x) length(unique(x))) ## How many different types 
Multb_G2_vc <- as.data.frame(table(unlist(SI_AtoMultb_G2))) ## Count of how many appearances of each
colnames(Multb_G2_vc) <- c("Type", "Appearances")

## Prepare for density plot of first element different from 1
## df
dens_Multb_G2_fd <- data.frame("x" = density(Multb_G2_fd)$x,
		 	       "y" = density(Multb_G2_fd)$y)

## cis
Multb_G2_fd_quant <- quantile(Multb_G2_fd, probs = c(0.025, 0.5, 0.975))

dens_Multb_G2_fd_ci <- dens_Multb_G2_fd[dens_Multb_G2_fd$x <= Multb_G2_fd_quant[3] & dens_Multb_G2_fd$x >= Multb_G2_fd_quant[1],] 

dens_Multb_G2_fd_ci$y[1] <- 0
dens_Multb_G2_fd_ci$y[nrow(dens_Multb_G2_fd_ci)] <- 0

## Prepare index for barplot. Select five with most appearances
xlab_n <- rep("",nrow(Multb_G2_vc))
sel_types <- order(Multb_G2_vc$Appearances, decreasing = TRUE)[1:5]
xlab_n[sel_types] <- as.character(Multb_G2_vc$Type[sel_types])

png("../SI_figures/SI3/Fig_7_SI3.png", res = 130, height = 1200, width = 1200)
layout(matrix(c(1,2,3,3),ncol = 2, byrow = T))
plot(x=dens_Multb_G2_fd$x, y=dens_Multb_G2_fd$y, xlab = "Timestep to first different type", col = "blue", lwd = 1.5, main = "First different type (AtoMultb G2)", type = "l", ylab = "value")
polygon(x=dens_Multb_G2_fd$x, y=dens_Multb_G2_fd$y, col = "lightblue1", border = NA)
polygon(x=dens_Multb_G2_fd_ci$x, y=dens_Multb_G2_fd_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Multb_G2_fd$x, y=dens_Multb_G2_fd$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(Multb_G2_fd), lty = 2, col = "navyblue")
abline(v = median(Multb_G2_fd), lty = 2, col = "red")

wordcloud(words = as.character(Multb_G2_vc$Type), freq = Multb_G2_vc$Appearances, scale = c(2,0.5), colors = brewer.pal(8, "Dark2"), random.order = T)
title("Type frequency")

barplot(Multb_G2_vc$Appearances, names.arg = xlab_n, las = 2, col = "gold4", main = "Types present in the simulation", cex.names = 0.6)
dev.off()

## G9b

SI_AtoMultb_G9 <- readRDS("../SI3_results/SI_AtoMultb_G9_res.rds")
Multb_G9_fd <- sapply(SI_AtoMultb_G9, function(x) min(which(x != x[1]))) ## First element different from 1
Multb_G9_dc <- sapply(SI_AtoMultb_G9, function(x) length(unique(x))) ## How many different types 
Multb_G9_vc <- as.data.frame(table(unlist(SI_AtoMultb_G9))) ## Count of how many appearances of each
colnames(Multb_G9_vc) <- c("Type", "Appearances")

## Prepare for density plot of first element different from 1
## df
dens_Multb_G9_fd <- data.frame("x" = density(Multb_G9_fd)$x,
		 	       "y" = density(Multb_G9_fd)$y)

## cis
Multb_G9_fd_quant <- quantile(Multb_G9_fd, probs = c(0.025, 0.5, 0.975))

dens_Multb_G9_fd_ci <- dens_Multb_G9_fd[dens_Multb_G9_fd$x <= Multb_G9_fd_quant[3] & dens_Multb_G9_fd$x >= Multb_G9_fd_quant[1],] 

dens_Multb_G9_fd_ci$y[1] <- 0
dens_Multb_G9_fd_ci$y[nrow(dens_Multb_G9_fd_ci)] <- 0

## Prepare index for barplot. Select five with most appearances
xlab_n <- rep("",nrow(Multb_G9_vc))
sel_types <- order(Multb_G9_vc$Appearances, decreasing = TRUE)[1:5]
xlab_n[sel_types] <- as.character(Multb_G9_vc$Type[sel_types])

png("../SI_figures/SI3/Fig_8_SI3.png", res = 140, height = 1200, width = 1200)
layout(matrix(c(1,2,3,3),ncol = 2, byrow = T))
plot(x=dens_Multb_G9_fd$x, y=dens_Multb_G9_fd$y, xlab = "Timestep to first different type", col = "blue", lwd = 1.5, main = "First different type (AtoMultb G9)", type = "l", ylab = "value")
polygon(x=dens_Multb_G9_fd$x, y=dens_Multb_G9_fd$y, col = "lightblue1", border = NA)
polygon(x=dens_Multb_G9_fd_ci$x, y=dens_Multb_G9_fd_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Multb_G9_fd$x, y=dens_Multb_G9_fd$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(Multb_G9_fd), lty = 2, col = "navyblue")
abline(v = median(Multb_G9_fd), lty = 2, col = "red")

wordcloud(words = as.character(Multb_G9_vc$Type), freq = Multb_G9_vc$Appearances, scale = c(2,0.5), colors = brewer.pal(8, "Dark2"), random.order = T)
title("Type frequency")

barplot(Multb_G9_vc$Appearances, names.arg = xlab_n, las = 2, col = "gold4", main = "Types present in the simulation", cex.names = 0.4)
dev.off()

## G18b

SI_AtoMultb_G18 <- readRDS("../SI3_results/SI_AtoMultb_G18_res.rds")
Multb_G18_fd <- sapply(SI_AtoMultb_G18, function(x) min(which(x != x[1]))) ## First element different from 1
Multb_G18_dc <- sapply(SI_AtoMultb_G18, function(x) length(unique(x))) ## How many different types 
Multb_G18_vc <- as.data.frame(table(unlist(SI_AtoMultb_G18))) ## Count of how many appearances of each
colnames(Multb_G18_vc) <- c("Type", "Appearances")

## Prepare for density plot of first element different from 1
## df
dens_Multb_G18_fd <- data.frame("x" = density(Multb_G18_fd)$x,
		 	       "y" = density(Multb_G18_fd)$y)

## cis
Multb_G18_fd_quant <- quantile(Multb_G18_fd, probs = c(0.025, 0.5, 0.975))

dens_Multb_G18_fd_ci <- dens_Multb_G18_fd[dens_Multb_G18_fd$x <= Multb_G18_fd_quant[3] & dens_Multb_G18_fd$x >= Multb_G18_fd_quant[1],] 

dens_Multb_G18_fd_ci$y[1] <- 0
dens_Multb_G18_fd_ci$y[nrow(dens_Multb_G18_fd_ci)] <- 0

## Prepare index for barplot. Select five with most appearances
xlab_n <- rep("",nrow(Multb_G18_vc))
sel_types <- order(Multb_G18_vc$Appearances, decreasing = TRUE)[1:5]
xlab_n[sel_types] <- as.character(Multb_G18_vc$Type[sel_types])

png("../SI_figures/SI3/Fig_9_SI3.png", res = 140, height = 1200, width = 1200)
layout(matrix(c(1,2,3,3),ncol = 2, byrow = T))
plot(x=dens_Multb_G18_fd$x, y=dens_Multb_G18_fd$y, xlab = "Timestep to first different type", col = "blue", lwd = 1.5, main = "First different type (AtoMultb G18)", type = "l", ylab = "value")
polygon(x=dens_Multb_G18_fd$x, y=dens_Multb_G18_fd$y, col = "lightblue1", border = NA)
polygon(x=dens_Multb_G18_fd_ci$x, y=dens_Multb_G18_fd_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Multb_G18_fd$x, y=dens_Multb_G18_fd$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(Multb_G18_fd), lty = 2, col = "navyblue")
abline(v = median(Multb_G18_fd), lty = 2, col = "red")

wordcloud(words = as.character(Multb_G18_vc$Type), freq = Multb_G18_vc$Appearances, scale = c(2,0.5), colors = brewer.pal(8, "Dark2"), random.order = T)
title("Type frequency")

barplot(Multb_G18_vc$Appearances, names.arg = xlab_n, las = 2, col = "gold4", main = "Types present in the simulation", cex.names = 0.4)
dev.off()


################ EXTENDED SIMULATIONS AtoFree
## How much does it deviate? Maximum distance to the whole of the morphospace and to the initial shape

SI_Free <- readRDS("../SI3_results/SI_Free_res.rds")


## Prepare objects for plots
## df
dens_Free_max <- data.frame("x" = density(SI_Free$Dists_to_morph[,2])$x,
		 	    "y" = density(SI_Free$Dists_to_morph[,2])$y)

dens_Free_min <- data.frame("x" = density(SI_Free$Dists_to_morph[,5])$x,
		 	    "y" = density(SI_Free$Dists_to_morph[,5])$y)

dens_Free_max_pd <- data.frame("x" = density(SI_Free$Dists_to_morph[,1])$x,
		 	      "y" = density(SI_Free$Dists_to_morph[,1])$y)

dens_Free_min_pd <- data.frame("x" = density(SI_Free$Dists_to_morph[,4])$x,
		 	      "y" = density(SI_Free$Dists_to_morph[,4])$y)

## cis
Free_max_quant <- quantile(SI_Free$Dists_to_morph[,2], probs = c(0.025, 0.5, 0.975))
Free_min_quant <- quantile(SI_Free$Dists_to_morph[,5], probs = c(0.025, 0.5, 0.975))

Free_max_pd_quant <- quantile(SI_Free$Dists_to_morph[,1], probs = c(0.025, 0.5, 0.975))
Free_min_pd_quant <- quantile(SI_Free$Dists_to_morph[,4], probs = c(0.025, 0.5, 0.975))

dens_Free_max_ci <- dens_Free_max[dens_Free_max$x <= Free_max_quant[3] & dens_Free_max$x >= Free_max_quant[1],] 
dens_Free_min_ci <- dens_Free_min[dens_Free_min$x <= Free_min_quant[3] & dens_Free_min$x >= Free_min_quant[1],] 

dens_Free_max_pd_ci <- dens_Free_max_pd[dens_Free_max_pd$x <= Free_max_pd_quant[3] & dens_Free_max_pd$x >= Free_max_pd_quant[1],] 
dens_Free_min_pd_ci <- dens_Free_min_pd[dens_Free_min_pd$x <= Free_min_pd_quant[3] & dens_Free_min_pd$x >= Free_min_pd_quant[1],] 

dens_Free_max_ci$y[1] <- 0
dens_Free_max_ci$y[nrow(dens_Free_max_ci)] <- 0

dens_Free_min_ci$y[1] <- 0
dens_Free_min_ci$y[nrow(dens_Free_min_ci)] <- 0

dens_Free_max_pd_ci$y[1] <- 0
dens_Free_max_pd_ci$y[nrow(dens_Free_max_pd_ci)] <- 0

dens_Free_min_pd_ci$y[1] <- 0
dens_Free_min_pd_ci$y[nrow(dens_Free_min_pd_ci)] <- 0

png("../SI_figures/SI3/Fig_10_SI3.png", res = 140, height = 1200, width = 1200)
par(mfrow = c(2,2))

## Divergence times

plot(x=dens_Free_max$x, y=dens_Free_max$y, xlab = "Timestep of max divergence", col = "blue", lwd = 1.5, main = "Max divergence Free simulation", type = "l", ylab = "value")
polygon(x=dens_Free_max$x, y=dens_Free_max$y, col = "lightblue1", border = NA)
polygon(x=dens_Free_max_ci$x, y=dens_Free_max_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Free_max$x, y=dens_Free_max$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(SI_Free$Dists_to_morph[,2]), lty = 2, col = "navyblue")
abline(v = median(SI_Free$Dists_to_morph[,2]), lty = 2, col = "red")
plot(x=dens_Free_min$x, y = dens_Free_min$y, xlab = "Timestep of first divergence", col = "blue", lwd = 1.5, main = "First divergence Free simulation", type = "l", ylab = "value")
polygon(x=dens_Free_min$x, y=dens_Free_min$y, col = "lightblue1", border = NA)
polygon(x=dens_Free_min_ci$x, y=dens_Free_min_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Free_min$x, y=dens_Free_min$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(SI_Free$Dists_to_morph[,5]), lty = 2, col = "navyblue")
abline(v = median(SI_Free$Dists_to_morph[,5]), lty = 2, col = "red")

## Procrustes distances

plot(x=dens_Free_max_pd$x, y=dens_Free_max_pd$y, xlab = "Maximum Procrustes distance", col = "blue", lwd = 1.5, main = "Max distance to morphospace Free", type = "l", ylab = "value")
polygon(x=dens_Free_max_pd$x, y=dens_Free_max_pd$y, col = "lightblue1", border = NA)
polygon(x=dens_Free_max_pd_ci$x, y=dens_Free_max_pd_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Free_max_pd$x, y=dens_Free_max_pd$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(SI_Free$Dists_to_morph[,1]), lty = 2, col = "navyblue")
abline(v = median(SI_Free$Dists_to_morph[,1]), lty = 2, col = "red")
plot(x=dens_Free_min_pd$x, y = dens_Free_min_pd$y, xlab = "Minimum Procrustes distance", col = "blue", lwd = 1.5, main = "Min distance to morphospace Free", type = "l", ylab = "value")
polygon(x=dens_Free_min_pd$x, y=dens_Free_min_pd$y, col = "lightblue1", border = NA)
polygon(x=dens_Free_min_pd_ci$x, y=dens_Free_min_pd_ci$y, col = "lightblue3", border = NA)
polygon(x=dens_Free_min_pd$x, y=dens_Free_min_pd$y, col = adjustcolor("blue", alpha = 0), border = "black")
abline(v = mean(SI_Free$Dists_to_morph[,4]), lty = 2, col = "navyblue")
abline(v = median(SI_Free$Dists_to_morph[,4]), lty = 2, col = "red")

dev.off()
