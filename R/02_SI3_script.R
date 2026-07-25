####################################################################################
###################################################################################
###################################################################################
#a# Scripts for Cortell-Nicolau, A., Kandler, A., 'Simulating shape variation in material culture: The Simumorph R-package', Journal of Archaeological Method and Theory
####################################################################################

####################################################################################
######## The following SI script includes extended simulations to analyse the model. It explores the variation in different ways of the four implemented methods
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

################ ASSESS RESTRICTIONS

set.seed(seed)

## Unrestricted
SI_AtoA_unres <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 4, e = 1, f = 1, int.allowed = TRUE)
saveRDS(SI_AtoA_unres, "../SI3_results/SI_AtoA_unres.rds")

## Avoid self-crossing
SI_AtoA_nocross <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 4, e = 1, f = 1)
saveRDS(SI_AtoA_nocross, "../SI3_results/SI_AtoA_nocross.rds")

## Include s only
SI_AtoA_s <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 4, e = 1, f = 100)
saveRDS(SI_AtoA_s, "../SI3_results/SI_AtoA_s.rds")

## Include epsilon only
## Commented because it doesn't produce results. alpha and s are too large for the algorithm to find e
#SI_AtoA_epsilon <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 4, e = 0.03, f = 1)
#saveRDS(SI_AtoA_epsilon, "../SI3_results/SI_AtoA_epsilon.rds")

## Include alpha only
## Same as above
#SI_AtoA_alpha <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 0.2, e = 1, f = 1)
#saveRDS(SI_AtoA_alpha, "../SI3_results/SI_AtoA_alpha.rds")

## Include s and epsilon
## Works, but veeery slow!
SI_AtoA_epsilon <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 4, e = 0.03, f = 100)
saveRDS(SI_AtoA_epsilon, "../SI3_results/SI_AtoA_s_epsilon.rds")

## Include s and alpha
SI_AtoA_alpha <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 0.2, e = 1, f = 100)
saveRDS(SI_AtoA_alpha, "../SI3_results/SI_AtoA_s_alpha.rds")

## Include alpha and epsilon
## Again, no restults. As above
#SI_AtoA_alpha <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 0.2, e = 0.03, f = 1)
#saveRDS(SI_AtoA_alpha, "../SI3_results/SI_AtoA_alpha_epsilon.rds")

## Include all of them
SI_AtoA_alpha <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G9_m_G9"), target = which(rownames(amp_pha_mat) == "G9_m_G9"), method = "AtoA", sim = sims, npts = npts, only.shapes = F, a = 0.2, e = 0.03, f = 100)
saveRDS(SI_AtoA_alpha, "../SI3_results/SI_AtoA_all.rds")


################ EXTENDED SIMULATIONS AtoA
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
saveRDS(SI_AtoA_res, "../SI3_results/SI_AtoA_res.rds")

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
saveRDS(SI_AtoB_res, "../SI3_results/SI_AtoB_res.rds")

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
saveRDS(SI_AtoMulta_G1_shapes, "../SI3_results/SI_AtoMulta_G1_res.rds")

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
saveRDS(SI_AtoMulta_G2_shapes, "../SI3_results/SI_AtoMulta_G2_res.rds")

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
saveRDS(SI_AtoMulta_G9_shapes, "../SI3_results/SI_AtoMulta_G9_res.rds")

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
saveRDS(SI_AtoMulta_G18_shapes, "../SI3_results/SI_AtoMulta_G18_res.rds")

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
saveRDS(SI_AtoMultb_G1_shapes, "../SI3_results/SI_AtoMultb_G1_res.rds")

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
saveRDS(SI_AtoMultb_G2_shapes, "../SI3_results/SI_AtoMultb_G2_res.rds")

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
saveRDS(SI_AtoMultb_G9_shapes, "../SI3_results/SI_AtoMultb_G9_res.rds")

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
saveRDS(SI_AtoMultb_G18_shapes, "../SI3_results/SI_AtoMultb_G18_res.rds")

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

saveRDS(SI_Free_res, "../SI3_results/SI_Free_res.rds")


###################################################################################
###################################################################################

######## SENSITIVITY ANALYSIS

####################################################################################
####################################################################################

## Load library and data (already loaded above)
#library(Simumorph)
#library(parallel)

## Load necessary objects (produced with morphospace.R)
#geo_out <- readRDS("../Utilities/geo_out.rds") ## Observed shapes
#morpho_pars <- readRDS("../Utilities/morpho_pars.rds") ## Observed parameters
#amp_pha_mat <- readRDS("../Utilities/amp_pha_mat.rds") ## Amplitude and phase matrix
#amp_pha_cov <- readRDS("../Utilities/amp_pha_cov.rds") ## Covariance matrix

## Load utilities
sims <- 5
sims_in <- 100
npts <- 120

## Wrapper to repeat simulation with different seed if max attempts is reached within simumorph
safe_simumorph <- function(seed, ..., max_attempts = 50){
	attempt <- 1
	repeat {
		## Seed per entry
		current_seed <- seed + attempt
		set.seed(current_seed)
		
		## Make sure the iteration goes through
		temp <- try(simumorph(...),silent = TRUE)
		
		## When it works
		if (!inherits(temp, "try-error")){
			return(list(result = temp,seed_used = current_seed,attempts = attempt))
		}
		
		## In case it stals for too long
		if (attempt >= max_attempts){
			stop(paste("Maximum retry attempts reached for seed",seed))
		}
		attempt <- attempt + 1
	}
}

###########################
###########################
##### AtoA
###########################
###########################
store_names <- c("alpha","epsilon","s","maxPd","tmaxPd","minPd","tminPd","seed","attempts")
AtoA_list_temp <- vector("list", length(sims))
a_min <- 0.1
a_max <- pi
e_val <- 0.03
s_min <- 50
s_max <- 170

## Simulate
AtoA_par_sens <- function(i){
	## Generate reproducible seed for this iteration
	sim_seed <- sample.int(.Machine$integer.max, 1)
        set.seed(sim_seed)

        a_val <- runif(1, a_min, a_max)
        s_val <- runif(1, s_min, s_max)
	
        temp <- safe_simumorph(seed = sim_seed,x = amp_pha_cov,m.space = amp_pha_mat,init = which(rownames(amp_pha_mat) == "G9_m_G9"),target = which(rownames(amp_pha_mat) == "G9_m_G9"),method = "AtoA",sim = sims_in,npts = npts,only.shapes = FALSE,a = a_val,e = e_val,f = s_val)
	    
        AtoA_sens <- temp$result
        AtoA_temp <- c(a_val,e_val,s_val,max(AtoA_sens$P.distances),which.max(AtoA_sens$P.distances),min(AtoA_sens$P.distances),which.min(AtoA_sens$P.distances),temp$seed_used,temp$attempts)
    	names(AtoA_temp) <- store_names

	P_dists <- temp$result$P.distances

        AtoA_list_temp[[i]] <- AtoA_sens
        names(AtoA_list_temp) <- paste0("AtoA_sim_",i)
	return(list(df = AtoA_temp,full = AtoA_sens, P_dists = P_dists))
}

## Execute
AtoA_res <- mclapply(1:sims,AtoA_par_sens,mc.cores = detectCores() - 2,mc.set.seed = TRUE)

## Organise
AtoA_sens_df <- do.call(rbind,lapply(AtoA_res, function(x) x$df))
#AtoA_sens_full <- lapply(AtoA_res, function(x) x$full)
AtoA_sens_pdist <- do.call(rbind,lapply(AtoA_res, function(x) x$P_dists))

## Save
saveRDS(AtoA_sens_df,"../SI3_results/AtoA_sens_df.rds")
#saveRDS(AtoA_sens_full,"../SI3_results/AtoA_sens_full.rds")
saveRDS(AtoA_sens_pdist,"../SI3_results/AtoA_sens_pdist.rds")

###########################
###########################
##### AtoB
###########################
###########################

e_val <- 0.5
d_min <- -0.02
d_max <- 0.035

store_names2 <- c("alpha","epsilon","s","maxPd","tmaxPd","minPd","tminPd","delta","seed","attempts")
list_temp <- vector("list",length(sims))
sims_in <- 100

## Simulate
AtoB_par_sens <- function(i){
	sim_seed <- sample.int(.Machine$integer.max, 1)
	set.seed(sim_seed)
	
	a_val <- runif(1, a_min, a_max)
	s_val <- runif(1, s_min, s_max)
	d_val <- runif(1, d_min, d_max)
	
	temp <- safe_simumorph(seed = sim_seed,x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat)=="G2_m_G2"), target = which(rownames(amp_pha_mat)=="G18_m_G18"), method = "AtoB", sim = sims_in, npts = npts, only.shapes = F, a = a_val, e = e_val, f = s_val, max.attempts = 500, speedAtoB = d_val)
	
	AtoB_sens <- temp$result
	df_temp <- c(a_val,e_val,s_val,max(AtoB_sens$P.distances),which.max(AtoB_sens$P.distances),min(AtoB_sens$P.distances),which.min(AtoB_sens$P.distances),d_val,temp$seed_used,temp$attempts)
	names(df_temp) <- store_names2
	
	list_temp[[i]] <- AtoB_sens
	names(list_temp) <- paste0("AtoB_sim_s_", s_val, "_", i)
	return(list(df = df_temp, full = list_temp))
}

## Execute
AtoB_res <- mclapply(1:sims,AtoB_par_sens,mc.cores = detectCores() - 2)

## Organise
AtoB_sens_df <- do.call(rbind,lapply(AtoB_res, function(x) x$df))
#AtoB_sens_full <- lapply(AtoB_res, function(x) x$full)

## Save
saveRDS(AtoB_sens_df,"../SI3_results/AtoB_sens_df.rds")
#saveRDS(AtoB_sens_full,"../SI3_results/AtoB_sens_full.rds")

###########################
###########################
##### AtoMult
###########################
###########################
AtoMult_sens_full_temp <- list()

e_val <- 0.05
s_min <- 50
s_max <- 170

list_temp <- vector("list",length(sims))
targets_multi <- c(1:nrow(amp_pha_mat))

## Simulate
AtoMult_par_sens <- function(i){
	## Generate reprodusible seed for this iteration
	sim_seed <- sample.int(.Machine$integer.max, 1)
	set.seed(sim_seed)
	
	a_val <- runif(1, a_min, a_max)
	s_val <- runif(1, s_min, s_max)
	

	temp <- safe_simumorph(seed = sim_seed, x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G1.1_m_G1.1"), target = targets_multi, method = "AtoMult", sim = sims_in, npts = npts, only.shapes = F, a = a_val, e = e_val, f = s_val)
		
	AtoMult_sens <- temp$result
	## Arrange for single target
	mult_pd <- AtoMult_sens$P.distances
	wmax <- which(mult_pd == max(mult_pd), arr.ind = TRUE)[1,2]	
	wmin <- which(mult_pd == min(mult_pd), arr.ind = TRUE)[1,2]	
	df_temp <- c(a_val,e_val,s_val,max(mult_pd),wmax,min(mult_pd),wmin,temp$seed_used,temp$attempts)
    	names(df_temp) <- store_names
	
	list_temp[[i]] <- AtoMult_sens
	names(list_temp) <- paste0("AtoMult_sim_s_", s_val, "_", i)
	return(list(df = df_temp, full = list_temp))
}

## Execute
AtoMult_res <- mclapply(1:sims,AtoMult_par_sens,mc.cores = detectCores() - 2)

## Organise
AtoMult_sens_df <- do.call(rbind,lapply(AtoMult_res, function(x) x$df))
#AtoMult_sens_full <- lapply(AtoMult_res, function(x) x$full)

## Save
saveRDS(AtoMult_sens_df,"../SI3_results/AtoMult_sens_df.rds")
#saveRDS(AtoMult_sens_full,"../SI3_results/AtoMult_sens_full.rds")

###########################
###########################
##### Free
###########################
###########################
sims_in <- 500
list_temp <- vector("list",length(sims))

## Simulate
Free_par_sens <- function(i){
	## Generate reproducible seed for this iteration
	sim_seed <- sample.int(.Machine$integer.max, 1)
	set.seed(sim_seed)
	
	a_val <- runif(1, a_min, a_max)
	s_val <- runif(1, s_min, s_max)
	
	temp <- safe_simumorph(seed = sim_seed,x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = nrow(amp_pha_mat), method = "Free", sim = sims_in, npts = npts, only.shapes = F, a = a_val, e = e_val, f = s_val, max.attempts = 500)
		
	Free_sens <- temp$result
	df_temp <- c(a_val,e_val,s_val,max(Free_sens$P.distances),which.max(Free_sens$P.distances),min(Free_sens$P.distances),which.min(Free_sens$P.distances),temp$seed_used,temp$attempts)
	names(df_temp) <- store_names

	P_dists <- temp$result$P.distances
	list_temp[[i]] <- Free_sens
	names(list_temp) <- paste0("Free_sim_s_", s_val, "_", i)
	return(list(df = df_temp, full = list_temp, P_dists = P_dists))
}

## Execute
Free_res <- mclapply(1:sims,Free_par_sens,mc.cores = detectCores() - 2)

## Organise
Free_sens_df <- do.call(rbind,lapply(Free_res, function(x) x$df))
#Free_sens_full <- lapply(Free_res, function(x) x$full)
Free_sens_pdist <- t(vapply(Free_res, function(x) {x$P_dists}, numeric(500)))

## Save
saveRDS(Free_sens_df,"../SI3_results/Free_sens_df.rds")
#saveRDS(Free_sens_full,"../SI3_results/Free_sens_full.rds")
saveRDS(Free_sens_pdist,"../SI3_results/Free_sens_pdist.rds")





