
######################################################
### SIMULATION OF MORPHOMETRIC EVOLUTION OVER SINGLE SHAPE 
#######################################################


######################################################
#### Load libraries and utilities
#######################################################

## Libraries
library(Momocs) ## For GMM functions
library(vegan) ## For manual Procrustes
library(tmvtnorm) ## For truncated multivariate normal
library(sf) ## to check internal intersections

## Load necessary functions
source("trig_functions.R")

## Load necessary objects (produced with morphospace.R)
geo_out <- readRDS("../Utilities/geo_out.rds") ## Observed shapes
morpho_pars <- readRDS("../Utilities/morpho_pars.rds") ## Observed parameters
amp_pha_mat <- readRDS("../Utilities/amp_pha_mat.rds") ## Amplitude and phase matrix
amp_pha_cov <- readRDS("../Utilities/amp_pha_cov.rds") ## Covariance matrix


##### SIMULATION ON PLAUSIBLE SHAPES


## Other utilities
set.seed(1)
n.harm <- 7 # Number of harmonics
init <- 337 ## Which geometric will start the simulation (G17.1)
npts <- 120 ## Number of point to rebuild the outline
phase_thres <- 0.2 ## Limit to allow phase variation
sim <- 100 # number of simulations
proc_thres <- 0.1 # Threshold for distance
amp_pha <- amp_pha_mat[init,] ## Starting shape
sim_pars <- list() ## list to store the parameters of the simulation

# Thresholds for truncated multivariate normals
upper_thres <- apply(amp_pha_mat,2,max)
lower_thres <- apply(amp_pha_mat,2,min)

## Recalibrate phases to smooth jumps
lower_thres[15:28] <- unlist(as.data.frame(amp_pha_mat[init,15:28])) - phase_thres
upper_thres[15:28] <- unlist(as.data.frame(amp_pha_mat[init,15:28])) + phase_thres
 	
## Bound phases to [-pi,pi]	
lower_thres[15:28] <- pmin(pmax(lower_thres[15:28],-pi),pi)
upper_thres[15:28] <- pmin(pmax(upper_thres[15:28],-pi),pi)

## Create reference object to compare distances (will need later)
coeff_list <- list("an" = morpho_pars[26][1:7],
		   "bn" = morpho_pars[26][8:14],
		   "cn" = morpho_pars[26][15:21],
		   "dn" = morpho_pars[26][22:28])

ref_shape <- efourier_i(coeff_list, nb.h = n.harm, nb.pts = npts) # So far first shape
ref_shape <- coo_sample(ref_shape, n = npts)

## Build objects to store values
## Amplitudes and phases
Axmat <- matrix(nrow = n.harm, ncol = sim)
Aymat <- matrix(nrow = n.harm, ncol = sim)
Phyxmat <- matrix(nrow = n.harm, ncol = sim)
Phyymat <- matrix(nrow = n.harm, ncol = sim)

## Distances
proc_distances <- rep(NA,sim)

## Parameters
par_sim <- matrix(nrow = sim, ncol = n.harm*4)
colnames(par_sim) <- c(paste0("a",seq(1,n.harm)),paste0("b",seq(1,n.harm)),
		       paste0("c",seq(1,n.harm)),paste0("d",seq(1,n.harm)))

## Store simulation figures
simuls <- list()
index <- 1


png("../Figures/Fig_1a.png", res = 50, height = 1500, width = 1500)
par(mfrow = c(10,10))
#for (i in 1:sim){
while (index <= sim){
	## Set lower and upper threshold of the phases with a small kernel from the observed phases to avoid large jumps
	lower_thres[15:28] <- unlist(as.data.frame(amp_pha[15:28])) - phase_thres
	upper_thres[15:28] <- unlist(as.data.frame(amp_pha[15:28])) + phase_thres
 	
	## Bound to [-pi,pi]	
	lower_thres[15:28] <- pmin(pmax(lower_thres[15:28],-pi),pi)
	upper_thres[15:28] <- pmin(pmax(upper_thres[15:28],-pi),pi)

	cand_vals <- rtmvnorm(1,mean = unlist(as.data.frame(amp_pha)), sigma = amp_pha_cov/100, lower = lower_thres, upper = upper_thres)
	
	Ax <- cand_vals[c(1,3,5,7,9,11,13)]
	Ay <- cand_vals[c(2,4,6,8,10,12,14)]
	Phyx <- cand_vals[c(15,17,19,21,23,25,27)]
	Phyy <- cand_vals[c(16,18,20,22,24,26,28)]

	## cos and sin inverted to comply with Momocs
	mat_par_vals <- data.frame("an" = Ax*sin(Phyx),
				   "bn" = Ax*cos(Phyx),
			     	   "cn" = Ay*sin(Phyy),
			     	   "dn" = Ay*cos(Phyy))	
	
	## Truncating amplitudes there's no need for weighting
	mat_par_vals_l <- list("an" = mat_par_vals[,1],
			     "bn" = mat_par_vals[,2],
			     "cn" = mat_par_vals[,3],
			     "dn" = mat_par_vals[,4])

	## Build candidate_f
	candidate_f <- efourier_i(mat_par_vals_l, nb.h = n.harm, nb.pts = npts) # So far first shape
	candidate_f_check <- rbind(candidate_f,candidate_f[1,])## Close polygon for validity
	candidate_f_check <- st_polygon(list(candidate_f_check))
	
	if(st_is_valid(candidate_f_check)){
		## Manually center, scale and rotate shapes for distance computation	
		candidate_s <- coo_sample(coo_center(coo_scale(candidate_f)), n  = npts)
	
		# Center
		ref_centered <- scale(ref_shape, scale = FALSE, center = TRUE)
		cand_centered <- scale(candidate_s, scale = FALSE, center = TRUE)
	
		# Scale to unit size
		ref_scaled <- ref_centered / sqrt(sum(ref_centered^2))
		cand_scaled <- cand_centered / sqrt(sum(cand_centered^2))
	
		# Rotate candidate for best match
		svd_result <- svd(t(ref_scaled) %*% cand_scaled)
		rotation_matrix <- svd_result$v %*% t(svd_result$u)
		cand_aligned <- cand_scaled %*% rotation_matrix

		## Compute Procrustes distance
		proc_dist <- sqrt(procrustes(ref_scaled,cand_aligned,scale = TRUE)$ss)
		
		if (proc_dist <= proc_thres){
			## Plot without rotation for visualisation
			plot(cand_scaled, type = "l", lwd = 1.5, asp = 1, cex.main = 2,  main = paste0("t = ",index), xlab = "", ylab = "")
			polygon(cand_scaled, col = "lightsalmon3")

			## Store values for posterior plotting and assessment
			simuls[[index]] <- cand_aligned
			par_sim[index,1:7] <- mat_par_vals[,1]
			par_sim[index,8:14] <- mat_par_vals[,2]
			par_sim[index,15:21] <- mat_par_vals[,3]
			par_sim[index,22:28] <- mat_par_vals[,4]
			amp_pha <- cand_vals
		
			print(index) ## To see progress
			index <- index + 1
		}
	}
}

dev.off()


###### RELOAD VALUES FOR SIMULATION ON 'NOT' PLAUSIBLE SHAPES

phase_thres <- 1 ## Less control in phase variation to allow wiggles

# Thresholds for truncated multivariate normals
upper_thres <- apply(amp_pha_mat,2,max)
lower_thres <- apply(amp_pha_mat,2,min)

## Recalibrate phases to smooth jumps
lower_thres[15:28] <- unlist(as.data.frame(amp_pha_mat[init,15:28])) - phase_thres
upper_thres[15:28] <- unlist(as.data.frame(amp_pha_mat[init,15:28])) + phase_thres
 	
## Bound phases to [-pi,pi]	
lower_thres[15:28] <- pmin(pmax(lower_thres[15:28],-pi),pi)
upper_thres[15:28] <- pmin(pmax(upper_thres[15:28],-pi),pi)

## Create reference object to compare distances (will need later)
coeff_list <- list("an" = morpho_pars[26][1:7],
		   "bn" = morpho_pars[26][8:14],
		   "cn" = morpho_pars[26][15:21],
		   "dn" = morpho_pars[26][22:28])

ref_shape <- efourier_i(coeff_list, nb.h = n.harm, nb.pts = npts) # So far first shape
ref_shape <- coo_sample(ref_shape, n = npts)

## Build objects to store values
## Amplitudes and phases
Axmat <- matrix(nrow = n.harm, ncol = sim)
Aymat <- matrix(nrow = n.harm, ncol = sim)
Phyxmat <- matrix(nrow = n.harm, ncol = sim)
Phyymat <- matrix(nrow = n.harm, ncol = sim)

## Distances
proc_distances <- rep(NA,sim)

## Parameters
par_sim_np <- matrix(nrow = sim, ncol = n.harm*4)
colnames(par_sim_np) <- c(paste0("a",seq(1,n.harm)),paste0("b",seq(1,n.harm)),
		       paste0("c",seq(1,n.harm)),paste0("d",seq(1,n.harm)))

## Store simulation figures
simuls <- list()
index <- 1


png("../Figures/Fig_1b.png", res = 50, height = 1500, width = 1500)
par(mfrow = c(10,10))
#for (i in 1:sim){
while (index <= sim){
	## Set lower and upper threshold of the phases with a small kernel from the observed phases to avoid large jumps
	lower_thres[15:28] <- unlist(as.data.frame(amp_pha[15:28])) - phase_thres
	upper_thres[15:28] <- unlist(as.data.frame(amp_pha[15:28])) + phase_thres
 	
	## Bound to [-pi,pi]	
	lower_thres[15:28] <- pmin(pmax(lower_thres[15:28],-pi),pi)
	upper_thres[15:28] <- pmin(pmax(upper_thres[15:28],-pi),pi)

	cand_vals <- rtmvnorm(1,mean = unlist(as.data.frame(amp_pha)), sigma = amp_pha_cov, lower = lower_thres, upper = upper_thres)
	
	Ax <- cand_vals[c(1,3,5,7,9,11,13)]
	Ay <- cand_vals[c(2,4,6,8,10,12,14)]
	Phyx <- cand_vals[c(15,17,19,21,23,25,27)]
	Phyy <- cand_vals[c(16,18,20,22,24,26,28)]

	## cos and sin inverted to comply with Momocs
	mat_par_vals <- data.frame("an" = Ax*sin(Phyx),
				   "bn" = Ax*cos(Phyx),
			     	   "cn" = Ay*sin(Phyy),
			     	   "dn" = Ay*cos(Phyy))	
	
	## Truncating amplitudes there's no need for weighting
	mat_par_vals_l <- list("an" = mat_par_vals[,1],
			     "bn" = mat_par_vals[,2],
			     "cn" = mat_par_vals[,3],
			     "dn" = mat_par_vals[,4])

	## Simluation is totall free (includes wiggles)
	## Build candidate_f
	candidate_f <- efourier_i(mat_par_vals_l, nb.h = n.harm, nb.pts = npts) # So far first shape
			## Plot without rotation for visualisation
			plot(candidate_f, type = "l", lwd = 1.5, asp = 1,  cex.main = 2, main = paste0("t = ",index), xlab = "", ylab = "")
			polygon(candidate_f, col = "lightsalmon3")
			## Store values for posterior plotting and assessment
			par_sim_np[index,1:7] <- mat_par_vals[,1]
			par_sim_np[index,8:14] <- mat_par_vals[,2]
			par_sim_np[index,15:21] <- mat_par_vals[,3]
			par_sim_np[index,22:28] <- mat_par_vals[,4]
			amp_pha <- cand_vals
		
			print(index) ## To see progress
			index <- index + 1
		#}
	#}
}

dev.off()


# PRODUCE PCAs

plausible_pca <- prcomp(par_sim, scale. = TRUE)

png("../Figures/Fig_1c.png", res = 300, height = 1800, width = 1800)
plot(plausible_pca$x[,1], plausible_pca$x[,2], xlab = "PC1", ylab = "PC2", col = "indianred4", main = "PCA on plausible shapes", pch = 16)
dev.off()


not_plausible_pca <- prcomp(par_sim_np, scale. = TRUE)

png("../Figures/Fig_1d.png", res = 300, height = 1800, width = 1800)
plot(not_plausible_pca$x[,1], not_plausible_pca$x[,2], xlab = "PC1", ylab = "PC2", col = "indianred4", main = "PCA on not plausible shapes", pch = 16)
dev.off()

