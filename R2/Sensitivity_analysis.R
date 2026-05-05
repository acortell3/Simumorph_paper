

###################################################################################
###################################################################################
## Sensitivity analysis for Cortell-Nicolau, A., Kandler, A., 'Simulating shape variation in material culture: The Simumorph R-package', Journal of Archaeological Method and Theory
####################################################################################
####################################################################################

## Load library and data
#devtools::install_github("acortell3/Simumorph") 
library(Simumorph)
library(parallel)

## Load necessary objects (produced with morphospace.R)
geo_out <- readRDS("../Utilities/geo_out.rds") ## Observed shapes
morpho_pars <- readRDS("../Utilities/morpho_pars.rds") ## Observed parameters
amp_pha_mat <- readRDS("../Utilities/amp_pha_mat.rds") ## Amplitude and phase matrix
amp_pha_cov <- readRDS("../Utilities/amp_pha_cov.rds") ## Covariance matrix

## Load utilities
sims <- 1000
sims_in <- 100
npts <- 120
	
mat <- matrix(c(seq(1,100), rep(101,20)), nrow = 12, byrow = T) ## For plots
#seed <- 123

### AtoA
## What is the maximum P_d achieved depending on the parameters?
## When is that distance achieved?
## What is the minimum distance achieved depending on the parameters?
## When is that distance achieved?

## Set utilities
AtoA_sens_df <- data.frame("alpha" = numeric(),
			   "epsilon" = numeric(),
			   "s" = numeric(),
			   "maxPd" = numeric(),
			   "tmaxPd" = numeric(),
			   "minPd" = numeric(),
			   "tminPd" = numeric(),
			   "seed" = numeric())
AtoA_sens_full <- list()

AtoA_sens_full_temp <- list()

a_min <- 0.1
a_max <- 0.5
e_min <- 0.03
e_max <- 0.07
f_vals <- seq(50,170,30)

AtoA_par_sens <- function(i){
	set.seed(i)
	
	a_val <- runif(1, a_min, a_max)
	e_val <- runif(1, e_min, e_max)
	
	df_temp <- data.frame(alpha = numeric(length(f_vals)),
			      epsilon = numeric(length(f_vals)),
			      s = numeric(length(f_vals)),
			      maxPd = numeric(length(f_vals)),
			      tmaxPd = numeric(length(f_vals)),
			      minPd = numeric(length(f_vals)),
			      tminPd = numeric(length(f_vals)),
			      seed = numeric(length(f_vals)))
	
	list_temp <- vector("list", length(f_vals))
	
	for (j in seq_along(f_vals)){

		AtoA_sens <- simumorph(x = amp_pha_cov,m.space = amp_pha_mat,init = which(rownames(amp_pha_mat) == "G9_m_G9"),target = which(rownames(amp_pha_mat) == "G9_m_G9"),method = "AtoA",sim = sims_in,npts = npts,only.shapes = FALSE,a = a_val,e = e_val,f = f_vals[j])
		
		df_temp[j, ] <- c(a_val,e_val,f_vals[j],max(AtoA_sens$P.distances),which.max(AtoA_sens$P.distances),min(AtoA_sens$P.distances),which.min(AtoA_sens$P.distances),i)
		
		list_temp[[j]] <- AtoA_sens
		names(list_temp)[j] <- paste0("AtoA_sim_s_", f_vals[j], "_", i)
	}
	return(list(df = df_temp, full = list_temp))
}

AtoA_res <- mclapply(1:sims,AtoA_par_sens,mc.cores = detectCores() - 2)
AtoA_sens_df <- do.call(rbind,lapply(AtoA_res, function(x) x$df))
AtoA_sens_full <- unlist(lapply(AtoA_res, function(x) x$full),recursive = FALSE)

saveRDS(AtoA_sens_df,"../Sens_results/AtoA_sens_df.rds")
saveRDS(AtoA_sens_full,"../Sens_results/AtoA_sens_full.rds")

### AtoB
## What is the maximum P_d achieved depending on the parameters?
## When is that distance achieved?
## What is the minimum distance achieved depending on the parameters?
## When is that distance achieved?

## Set utilities
AtoB_sens_df <- data.frame("alpha" = numeric(),
			   "epsilon" = numeric(),
			   "s" = numeric(),
			   "maxPd" = numeric(),
			   "tmaxPd" = numeric(),
			   "minPd" = numeric(),
			   "tminPd" = numeric(),
			   "seed" = numeric())
AtoB_sens_full <- list()

AtoB_sens_full_temp <- list()

a_min <- 0.3
a_max <- 1
e_min <- 0.5
e_max <- 1

AtoB_par_sens <- function(i){
	set.seed(i)
	
	a_val <- runif(1, a_min, a_max)
	e_val <- runif(1, e_min, e_max)
	
	df_temp <- data.frame(alpha = numeric(length(f_vals)),
			      epsilon = numeric(length(f_vals)),
			      s = numeric(length(f_vals)),
			      maxPd = numeric(length(f_vals)),
			      tmaxPd = numeric(length(f_vals)),
			      minPd = numeric(length(f_vals)),
			      tminPd = numeric(length(f_vals)),
			      seed = numeric(length(f_vals)))
	
	list_temp <- vector("list", length(f_vals))
	
	for (j in seq_along(f_vals)){
	
		AtoB_sens <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat)=="G2_m_G2"), target = which(rownames(amp_pha_mat)=="G18_m_G18"), method = "AtoB", sim = sims_in, npts = npts, only.shapes = F, a = a_val, e = e_val, f = f_vals[j], max.attempts = 500, speedAtoB = 0)
		
		df_temp[j, ] <- c(a_val,e_val,f_vals[j],max(AtoB_sens$P.distances),which.max(AtoB_sens$P.distances),min(AtoB_sens$P.distances),which.min(AtoB_sens$P.distances),i)
		
		list_temp[[j]] <- AtoB_sens
		names(list_temp)[j] <- paste0("AtoB_sim_s_", f_vals[j], "_", i)
	}
	return(list(df = df_temp, full = list_temp))
}

AtoB_res <- mclapply(1:sims,AtoB_par_sens,mc.cores = detectCores() - 2)
AtoB_sens_df <- do.call(rbind,lapply(AtoB_res, function(x) x$df))
AtoB_sens_full <- unlist(lapply(AtoB_res, function(x) x$full),recursive = FALSE)

saveRDS(AtoB_sens_df,"../Sens_results/AtoB_sens_df.rds")
saveRDS(AtoB_sens_full,"../Sens_results/AtoB_sens_full.rds")

### AtoMult
## What is the maximum P_d achieved depending on the parameters?
## When is that distance achieved?
## What is the minimum distance achieved depending on the parameters?
## When is that distance achieved?

## Set utilities
AtoMult_sens_df <- data.frame("alpha" = numeric(),
			   "epsilon" = numeric(),
			   "s" = numeric(),
			   "maxPd" = numeric(),
			   "tmaxPd" = numeric(),
			   "minPd" = numeric(),
			   "tminPd" = numeric(),
			   "seed" = numeric())
AtoMult_sens_full <- list()

AtoMult_sens_full_temp <- list()

a_min <- 1.5
a_max <- 3.5
e_min <- 0.03
e_max <- 0.09

targets_multi <- c(1:nrow(amp_pha_mat))

AtoMult_par_sens <- function(i){
	set.seed(i)
	
	a_val <- runif(1, a_min, a_max)
	e_val <- runif(1, e_min, e_max)
	
	df_temp <- data.frame(alpha = numeric(length(f_vals)),
			      epsilon = numeric(length(f_vals)),
			      s = numeric(length(f_vals)),
			      maxPd = numeric(length(f_vals)),
			      tmaxPd = numeric(length(f_vals)),
			      minPd = numeric(length(f_vals)),
			      tminPd = numeric(length(f_vals)),
			      seed = numeric(length(f_vals)))
	
	list_temp <- vector("list", length(f_vals))
	
	for (j in seq_along(f_vals)){
	
		AtoMult_sens <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat) == "G1.1_m_G1.1"), target = targets_multi, method = "AtoMult", sim = sims_in, npts = npts, only.shapes = F, a = a_val, e = e_val, f = f_vals[j])
		
		## Arrange for single target
		mult_pd <- AtoMult_sens$P.distances
		wmax <- which(mult_pd == max(mult_pd), arr.ind = TRUE)[1,2]	
		wmin <- which(mult_pd == min(mult_pd), arr.ind = TRUE)[1,2]	
		df_temp[j, ] <- c(a_val,e_val,f_vals[j],max(mult_pd),wmax,min(mult_pd),wmin,i)
		
		list_temp[[j]] <- AtoMult_sens
		names(list_temp)[j] <- paste0("AtoMult_sim_s_", f_vals[j], "_", i)
	}
	return(list(df = df_temp, full = list_temp))
}

AtoMult_res <- mclapply(1:sims,AtoMult_par_sens,mc.cores = detectCores() - 2)
AtoMult_sens_df <- do.call(rbind,lapply(AtoMult_res, function(x) x$df))
AtoMult_sens_full <- unlist(lapply(AtoMult_res, function(x) x$full),recursive = FALSE)

saveRDS(AtoMult_sens_df,"../Sens_results/AtoMult_sens_df.rds")
saveRDS(AtoMult_sens_full,"../Sens_results/AtoMult_sens_full.rds")


### Free
## What is the maximum P_d achieved depending on the parameters?
## When is that distance achieved?
## What is the minimum distance achieved depending on the parameters?
## When is that distance achieved?

## Set utilities
Free_sens_df <- data.frame("alpha" = numeric(),
			   "epsilon" = numeric(),
			   "s" = numeric(),
			   "maxPd" = numeric(),
			   "tmaxPd" = numeric(),
			   "minPd" = numeric(),
			   "tminPd" = numeric(),
			   "seed" = numeric())
Free_sens_full <- list()

Free_sens_full_temp <- list()

a_min <- 1.5
a_max <- 3.5
e_min <- 0.03
e_max <- 0.09

targets_multi <- c(1:nrow(amp_pha_mat))

Free_par_sens <- function(i){
	set.seed(i)
	
	a_val <- runif(1, a_min, a_max)
	e_val <- runif(1, e_min, e_max)
	
	df_temp <- data.frame(alpha = numeric(length(f_vals)),
			      epsilon = numeric(length(f_vals)),
			      s = numeric(length(f_vals)),
			      maxPd = numeric(length(f_vals)),
			      tmaxPd = numeric(length(f_vals)),
			      minPd = numeric(length(f_vals)),
			      tminPd = numeric(length(f_vals)),
			      seed = numeric(length(f_vals)))
	
	list_temp <- vector("list", length(f_vals))
	
	for (j in seq_along(f_vals)){
		
		Free_sens <- simumorph(x = amp_pha_cov, m.space = amp_pha_mat, init = 1, target = nrow(amp_pha_mat), method = "Free", sim = sims_in, npts = npts, only.shapes = F, a = a_val, e = e_val, f = f_vals[j], max.attempts = 500)
		
		df_temp[j, ] <- c(a_val,e_val,f_vals[j],max(Free_sens$P.distances),which.max(Free_sens$P.distances),min(Free_sens$P.distances),which.min(Free_sens$P.distances),i)
		
		list_temp[[j]] <- Free_sens
		names(list_temp)[j] <- paste0("Free_sim_s_", f_vals[j], "_", i)
	}
	return(list(df = df_temp, full = list_temp))
}

Free_res <- mclapply(1:sims,Free_par_sens,mc.cores = detectCores() - 2)
Free_sens_df <- do.call(rbind,lapply(Free_res, function(x) x$df))
Free_sens_full <- unlist(lapply(Free_res, function(x) x$full),recursive = FALSE)

saveRDS(Free_sens_df,"../Sens_results/Free_sens_df.rds")
saveRDS(Free_sens_full,"../Sens_results/Free_sens_full.rds")


