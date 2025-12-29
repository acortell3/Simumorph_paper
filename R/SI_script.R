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
n_iter <- 10
ncores <- 10
seed <- 123


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
saveRDS(SI_AtoA_res, "../SI_results/SI_AtoA_res.rds")

################ EXTENDED SIMULATIONS AtoB
## When do we reach maximum convergence? Need Procrustes distances, time of max convergence, and shapes themselves

SI_AtoB_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat)=="G2_m_G2"), target = which(rownames(amp_pha_mat)=="G18_m_G18"), method = "AtoB", sim = 500, npts = npts, only.shapes = F, a = 0.5, e = 0.5, f = 100, max.attempts = 500, speedAtoB = 0.0035)}

set.seed(123)

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoB <- mclapply(1:n_iter, SI_AtoB_fun, mc.cores = ncores)

AtoB_min_div_index <- vapply(SI_AtoB, function(x) {which.min(x$P.distance)}, integer(1))
AtoB_min_div <- vapply(SI_AtoB, function(x) (min(x$P.distance)), numeric(1))
AtoB_min_shapes <- lapply(seq_along(SI_AtoB), function(i) {SI_AtoB[[i]]$Shapes[[AtoB_min_div_index[i]]]})

## Save results
SI_AtoB_res <- list("indexes" = AtoB_max_div_index,
		    "max_dist" = AtoB_max_div,
		    "Shapes" = AtoB_max_shapes)
saveRDS(SI_AtoB_res, "../SI_results/SI_AtoB_res.rds")

################ EXTENDED SIMULATIONS AtoMult
## When is it a different shape? How many shapes does it go through? And to which shapes does it change? Comparison with and without initial shape in target pool

#### AtoMulta. Include initial shape in target pool
targets_multi <- c(1:nrow(amp_pha_mat))
target_names <- rownames(amp_pha_mat)
c_a <- 3
c_f <- 50


SI_AtoMulta_G1_fun <- function(i) {simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G1.1_m_G1.1"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)}

set.seed(123)

## !!!!!! WARNING!!!!!!!! Windows users, not that mclapply only works on Linux. You'll have to adapt this to foreach
SI_AtoMulta_G1 <- mclapply(1:n_iter, SI_AtoMulta_G1_fun, mc.cores = ncores)

SI_AtoMulta_G1_stats <- as.data.frame(matrix(NA, nrow = n_iter, ncol = sims))

## Retreive names of shapes
for (j in 1:n_iter){

	#AtoMult_G1 <- SI_AtoMulta_G1[[j]]
	## Create vector with minimum procrustes distances and their 
	#pdist_G1 <- data.frame("target" = character(),
	#		       "distance" = numeric())
	
	for (i in 1:sims){
		SI_AtoMulta_G1_stats[j,] <- target_names[which.min(SI_AtoMulta_G1[[j]]$P.distances[2:length(target_names),i])]
		#pdist_G1[i,1] <- target_names[which.min(AtoMult_G1a$P.distances[,i])]
		#pdist_G1[i,2] <- min(AtoMult_G1$P.distances[,i])
	}

	#SI_AtoMulta_G1_stats[[j]] <- pdist_G1
}


AtoMulta_min_div_index <- vapply(SI_AtoMulta, function(x) {which.min(x$P.distance)}, integer(1))
AtoMulta_min_div <- vapply(SI_AtoMulta, function(x) (min(x$P.distance)), numeric(1))
AtoMulta_min_shapes <- lapply(seq_along(SI_AtoMulta), function(i) {SI_AtoMulta[[i]]$Shapes[[AtoMulta_min_div_index[i]]]})

## Save results
SI_AtoMulta_res <- list("indexes" = AtoMulta_max_div_index,
			"max_dist" = AtoMulta_max_div,
			"Shapes" = AtoMulta_max_shapes)

saveRDS(SI_AtoMulta_res, "../SI_results/SI_AtoMulta_res.rds")
### AtoMultb. Don't include initial shape in target pool
### Let's try with the full morphospace
target_names <- rownames(amp_pha_mat)[!grepl("G1.1",rownames(amp_pha_mat))]
targets_multi <- grep("G1.1", rownames(amp_pha_mat), invert = TRUE)

c_a <- 3
c_f <- 50


SI_AtoMult_without_init <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G1.1_m_G1.1"), target = targets_multi, method = "AtoMult", sim = sims, npts = npts, only.shapes = F, a = c_a, e = 0.05, f = c_f)


################ EXTENDED SIMULATIONS AtoFree
## How much does it deviate? Maximum distance to the whole of the morphospace and to the initial shape

SI_Free <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = nrow(amp_pha_mat), method = "Free", sim = sims, npts = npts, only.shapes = F, a = 0.2, e = 0.05, f = 100, max.attempts = 500)


