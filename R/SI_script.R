#####################################################################################
######## The following SI script includes extended simulations to analyse the model in the paper for <paper name>. It explores the variation in different ways of the four implemented methods
#####################################################################################

## Load libraries and data
library(Simumorph)
library(parallel)

## Load necessary objects (produced with morphospace.R)
geo_out <- readRDS("../Utilities/geo_out.rds") ## Observed shapes
morpho_pars <- readRDS("../Utilities/morpho_pars.rds") ## Observed parameters
amp_pha_mat <- readRDS("../Utilities/amp_pha_mat.rds") ## Amplitude and phase matrix
amp_pha_cov <- readRDS("../Utilities/amp_pha_cov.rds") ## Covariance matrix

## Utilities
sims <- 100
npts <- 120
n_iter <- 100
ncores <- 10
seed <- 123

#a############### EXTENDED SIMULATIONS AtoA
## 100 most distant shapes to initial shapes. Need Procrustes distances, time of max divergence, and shapes themselves

SI_AtoA_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 0.2, e = 0.03, f = 100)}

RNGkind("L'Ecuyer-CMRG") ### To preserve internal seeds
set.seed(seed)

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoA <- mclapply(1:n_iter, SI_AtoA_fun, mc.cores = ncores)

## Extract needed statistics
AtoA_max_div_index <- vapply(SI_AtoA, function(x) {which.max(x$P.distance)}, integer(1))
AtoA_max_div <- vapply(SI_AtoA, function(x) (max(x$P.distance)), numeric(1))
AtoA_max_shapes <- lapply(seq_along(SI_AtoA), function(i) {SI_AtoA[[i]]$Shapes[[AtoA_max_div_index[i]]]})

## Save results
SI_AtoA_res <- list("indexes" = AtoA_max_div_index,
		    "max_dist" = AtoA_max_div,
		    "Shapes" = AtoA_max_shapes)
saveRDS(SI_AtoA_res, "../SI_results/SI_AtoA_res.rds")

################ EXTENDED SIMULATIONS AtoB
## When do we reach maximum convergence? Need Procrustes distances, time of max convergence, and shapes themselves

SI_AtoB_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat)=="G2_m_G2"), target = which(rownames(amp_pha_mat)=="G18_m_G18"), method = "AtoB", sim = 500, npts = npts, only.shapes = F, a = 0.5, e = 0.5, f = 100, max.attempts = 500, speedAtoB = 0.0035)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoB <- mclapply(1:n_iter, SI_AtoB_fun, mc.cores = ncores)

AtoB_min_div_index <- vapply(SI_AtoB, function(x) {which.min(x$P.distance)}, integer(1))
AtoB_min_div <- vapply(SI_AtoB, function(x) (min(x$P.distance)), numeric(1))
AtoB_min_shapes <- lapply(seq_along(SI_AtoB), function(i) {SI_AtoB[[i]]$Shapes[[AtoB_min_div_index[i]]]})

## Save results
SI_AtoB_res <- list("indexes" = AtoB_min_div_index,
		    "min_dist" = AtoB_min_div,
		    "Shapes" = AtoB_min_shapes)
saveRDS(SI_AtoB_res, "../SI_results/SI_AtoB_res.rds")

################ EXTENDED SIMULATIONS AtoMult
## When is it a different shape? How many shapes does it go through? And to which shapes does it change? Comparison with and without initial shape in target pool

## We use the same types as in the paper

#### G1

#### AtoMulta. Include initial shape in target pool
targets_multi <- c(1:nrow(amp_pha_mat))
target_names <- rownames(amp_pha_mat)
c_a <- 3
c_f <- 50

SI_AtoMulta_G1_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G1.1_m_G1.1"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoMulta_G1 <- mclapply(1:n_iter, SI_AtoMulta_G1_fun, mc.cores = ncores)

SI_AtoMulta_G1_shapes <- as.data.frame(matrix(NA,nrow = n_iter, ncol = sims))

for (i in 1:n_iter){
	index_min <- apply(SI_AtoMulta_G1[[i]]$P.distance,2,which.min)
	SI_AtoMulta_G1_shapes[i,] <- target_names[unlist(index_min)]
}

## Save results
saveRDS(SI_AtoMulta_G1_shapes, "../SI_results/SI_AtoMulta_G1_res.rds")

#### G2

SI_AtoMulta_G2_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G2_m_G2"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoMulta_G2 <- mclapply(1:n_iter, SI_AtoMulta_G2_fun, mc.cores = ncores)

SI_AtoMulta_G2_shapes <- as.data.frame(matrix(NA,nrow = n_iter, ncol = sims))

for (i in 1:n_iter){
	index_min <- apply(SI_AtoMulta_G2[[i]]$P.distance,2,which.min)
	SI_AtoMulta_G2_shapes[i,] <- target_names[unlist(index_min)]
}

## Save results
saveRDS(SI_AtoMulta_G2_shapes, "../SI_results/SI_AtoMulta_G2_res.rds")

#### G9

SI_AtoMulta_G9_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoMulta_G9 <- mclapply(1:n_iter, SI_AtoMulta_G9_fun, mc.cores = ncores)

SI_AtoMulta_G9_shapes <- as.data.frame(matrix(NA,nrow = n_iter, ncol = sims))

for (i in 1:n_iter){
	index_min <- apply(SI_AtoMulta_G9[[i]]$P.distance,2,which.min)
	SI_AtoMulta_G9_shapes[i,] <- target_names[unlist(index_min)]
}

## Save results
saveRDS(SI_AtoMulta_G9_shapes, "../SI_results/SI_AtoMulta_G9_res.rds")

#### G18
## Increased e so the sim does not crash
SI_AtoMulta_G18_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G18_m_G18"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.06, f = c_f)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoMulta_G18 <- mclapply(1:n_iter, SI_AtoMulta_G18_fun, mc.cores = ncores)

SI_AtoMulta_G18_shapes <- as.data.frame(matrix(NA,nrow = n_iter, ncol = sims))

for (i in 1:n_iter){
	index_min <- apply(SI_AtoMulta_G18[[i]]$P.distance,2,which.min)
	SI_AtoMulta_G18_shapes[i,] <- target_names[unlist(index_min)]
}

## Save results
saveRDS(SI_AtoMulta_G18_shapes, "../SI_results/SI_AtoMulta_G18_res.rds")

### AtoMultb. Don't include initial shape in target pool

#### G1

target_names <- rownames(amp_pha_mat)[!grepl("G1.1",rownames(amp_pha_mat))]
targets_multi <- grep("G1.1", rownames(amp_pha_mat), invert = TRUE)

c_a <- 3
c_f <- 50


SI_AtoMultb_G1_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G1.1_m_G1.1"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoMultb_G1 <- mclapply(1:n_iter, SI_AtoMultb_G1_fun, mc.cores = ncores)

SI_AtoMultb_G1_shapes <- as.data.frame(matrix(NA,nrow = n_iter, ncol = sims))

for (i in 1:n_iter){
	index_min <- apply(SI_AtoMultb_G1[[i]]$P.distance,2,which.min)
	SI_AtoMultb_G1_shapes[i,] <- target_names[unlist(index_min)]
}

## Save results
saveRDS(SI_AtoMultb_G1_shapes, "../SI_results/SI_AtoMultb_G1_res.rds")

#### G2
target_names <- rownames(amp_pha_mat)[!grepl("G2",rownames(amp_pha_mat))]
targets_multi <- grep("G2",rownames(amp_pha_mat), invert = TRUE)

SI_AtoMultb_G2_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G2_m_G2"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoMultb_G2 <- mclapply(1:n_iter, SI_AtoMultb_G2_fun, mc.cores = ncores)

SI_AtoMultb_G2_shapes <- as.data.frame(matrix(NA,nrow = n_iter, ncol = sims))

for (i in 1:n_iter){
	index_min <- apply(SI_AtoMultb_G2[[i]]$P.distance,2,which.min)
	SI_AtoMultb_G2_shapes[i,] <- target_names[unlist(index_min)]
}

## Save results
saveRDS(SI_AtoMultb_G2_shapes, "../SI_results/SI_AtoMultb_G2_res.rds")

#### G9
target_names <- rownames(amp_pha_mat)[!grepl("G9",rownames(amp_pha_mat))]
targets_multi <- grep("G9",rownames(amp_pha_mat), invert = TRUE)

## We need to do it with dynamic e because otherwise the algorithm is not able to escape from the G18 shape
dyn_e <- data.frame("time" = c(1,10,20,30,40,50),
		    "e" = c(0.14,0.13,0.12,0.11,0.1,0.09))

SI_AtoMultb_G9_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, dynamic_e = dyn_e, f = c_f)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoMultb_G9 <- mclapply(1:n_iter, SI_AtoMultb_G9_fun, mc.cores = ncores)

SI_AtoMultb_G9_shapes <- as.data.frame(matrix(NA,nrow = n_iter, ncol = sims))

for (i in 1:n_iter){
	index_min <- apply(SI_AtoMultb_G9[[i]]$P.distance,2,which.min)
	SI_AtoMultb_G9_shapes[i,] <- target_names[unlist(index_min)]
}

## Save results
saveRDS(SI_AtoMultb_G9_shapes, "../SI_results/SI_AtoMultb_G9_res.rds")

#### G18
target_names <- rownames(amp_pha_mat)[!grepl("G18",rownames(amp_pha_mat))]
targets_multi <- grep("G18",rownames(amp_pha_mat), invert = TRUE)

## We need to do it with dynamic e because otherwise the algorithm is not able to escape from the G18 shape
dyn_e <- data.frame("time" = c(1,10,20,30,40,75,90),
		    "e" = c(0.16,0.15,0.14,0.13,0.12,0.11,0.10))

SI_AtoMultb_G18_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G18_m_G18"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, dynamic_e = dyn_e, f = c_f)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoMultb_G18 <- mclapply(1:n_iter, SI_AtoMultb_G18_fun, mc.cores = ncores)

SI_AtoMultb_G18_shapes <- as.data.frame(matrix(NA,nrow = n_iter, ncol = sims))

for (i in 1:n_iter){
	index_min <- apply(SI_AtoMultb_G18[[i]]$P.distance,2,which.min)
	SI_AtoMultb_G18_shapes[i,] <- target_names[unlist(index_min)]
}

## Save results
saveRDS(SI_AtoMultb_G18_shapes, "../SI_results/SI_AtoMultb_G18_res.rds")

################ EXTENDED SIMULATIONS AtoFree
## How much does it deviate? Maximum distance to the whole of the morphospace and to the initial shape

SI_Free_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = nrow(amp_pha_mat), method = "Free", sim = sims, npts = npts, only.shapes = F, a = 0.2, e = 0.05, f = 100, max.attempts = 500)}

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_Free <- mclapply(1:n_iter, SI_Free_fun, mc.cores = ncores)

## Maximum distance to the initial shape
Free_max_div_index <- vapply(SI_Free, function(x) {which.max(x$P.distance)}, integer(1))
Free_max_div <- vapply(SI_Free, function(x) (max(x$P.distance)), numeric(1))
Free_max_shapes <- lapply(seq_along(SI_Free), function(i) {SI_Free[[i]]$Shapes[[Free_max_div_index[i]]]})

## Distances to the full morphospace

## Build target morphospace
tar_morph <- list()

for (i in 1:nrow(amp_pha_mat)){
	tar_morph[[i]] <- build_s(unlist(amp_pha_mat[i,]),fou.pars = F, npts = npts)
}

## Df to store distances
dists_to_morph <- data.frame("Max_distance" = rep(NA,length(SI_Free)),
			     "Which_max_distance" = rep(NA,length(SI_Free)),
			     "Type_max_distance" = rep(NA,length(SI_Free)),
			     "Min_distance" = rep(NA,length(SI_Free)),
			     "Which_min_distance" = rep(NA,length(SI_Free)),
			     "Type_min_distance" =  rep(NA,length(SI_Free)))


for (j in 1:length(SI_Free)){
	dists <- t(vapply(seq_len(sims), function(i) {proc_dist(SI_Free[[j]]$Shapes[[i]],tar_morph, multi = T)},numeric(length(tar_morph))))

	## The rationale from maximum distance is: Compute the minimum distance for each row (the closest target shape to the simulated shape) and, from there, which is the maximum distance. That is, the maximum distance to the most similar shape
	dists_to_morph[j,1] <- max(apply(dists,1,min))
	dists_to_morph[j,2] <- which(dists == dists_to_morph[j,1], arr.ind = T)[1]
	dists_to_morph[j,3] <- target_names[which.min(dists[dists_to_morph[j,2],])]

	dists_to_morph[j,4] <- min(dists)
	dists_to_morph[j,5] <- which(dists == min(dists), arr.ind = T)[1]
	dists_to_morph[j,6] <- target_names[which.min(dists[dists_to_morph[j,5],])]

}

#  Save results
SI_Free_res <- list("indexes" = Free_max_div_index,
		    "max_dist" = Free_max_div,
		    "Shapes" = Free_max_shapes,
		    "Dists_to_morph" = dists_to_morph)

saveRDS(SI_Free_res, "../SI_results/SI_Free_res.rds")
