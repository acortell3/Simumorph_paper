

###################################################################################
###################################################################################

####### SENSITIVITY ANALYSIS

####################################################################################
####################################################################################

## Load library and data (already loaded above)
library(Simumorph)
library(parallel)

## Load necessary objects (produced with morphospace.R)
geo_out <- readRDS("../Utilities/geo_out.rds") ## Observed shapes
morpho_pars <- readRDS("../Utilities/morpho_pars.rds") ## Observed parameters
amp_pha_mat <- readRDS("../Utilities/amp_pha_mat.rds") ## Amplitude and phase matrix
amp_pha_cov <- readRDS("../Utilities/amp_pha_cov.rds") ## Covariance matrix

## Load utilities
#sima <- 100
sims <- 100
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
#a_min <- 0.1
#a_max <- pi
e_val <- 0.03
#s_min <- 50
#s_max <- 170
a_vals <- seq(0.6,1.4,0.2)
s_vals <- seq(60,140,20)

## Simulate
#AtoA_par_sens <- function(i){
#	for (j in 1:5){
	## Generate reproducible seed for this iteration
#	sim_seed <- sample.int(.Machine$integer.max, 1)
#        set.seed(sim_seed)

        #a_val <- runif(1, a_min, a_max)
        #s_val <- runif(1, s_min, s_max)
#	a_val <- aval[i]
#	s_val <- sval[j]

#       temp <- safe_simumorph(seed = sim_seed,x = amp_pha_cov,m.space = amp_pha_mat,init = which(rownames(amp_pha_mat) == "G9_m_G9"),target = which(rownames(amp_pha_mat) == "G9_m_G9"),method = "AtoA",sim = sims_in,npts = npts,only.shapes = FALSE,a = a_val,e = e_val,f = s_val)
	    
#        AtoA_sens <- temp$result
#        AtoA_temp <- c(a_val,e_val,s_val,max(AtoA_sens$P.distances),which.max(AtoA_sens$P.distances),min(AtoA_sens$P.distances),which.min(AtoA_sens$P.distances),temp$seed_used,temp$attempts)
#    	names(AtoA_temp) <- store_names

#	P_dists <- temp$result$P.distances

#        AtoA_list_temp[[i]] <- AtoA_sens
#        names(AtoA_list_temp) <- paste0("AtoA_sim_",i)
	
#	return(list(df = AtoA_temp,full = AtoA_sens, P_dists = P_dists))
#
#	}
 #}


pars <- expand.grid(a = a_vals,s = s_vals,rep = 1:sims)

AtoA_par_sens <- function(i){

    sim_seed <- sample.int(.Machine$integer.max, 1)

    a_val <- pars$a[i]
    s_val <- pars$s[i]

    temp <- safe_simumorph(seed = sim_seed,x = amp_pha_cov,m.space = amp_pha_mat,init = which(rownames(amp_pha_mat) == "G9_m_G9"),target = which(rownames(amp_pha_mat) == "G9_m_G9"),method = "AtoA",sim = sims_in,npts = npts,only.shapes = FALSE,a = a_val,e = e_val,f = s_val)

    AtoA_sens <- temp$result

    AtoA_temp <- c(a_val,e_val,s_val,max(AtoA_sens$P.distances),which.max(AtoA_sens$P.distances),min(AtoA_sens$P.distances),which.min(AtoA_sens$P.distances),temp$seed_used,temp$attempts)

    names(AtoA_temp) <- store_names

    list(df = AtoA_temp,full = AtoA_sens,P_dists = AtoA_sens$P.distances)
}

## Execute
AtoA_res <- mclapply(1:nrow(pars),AtoA_par_sens,mc.cores = detectCores() - 2,mc.set.seed = TRUE)

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
a_val <- 0.8
s_val <- 100
#d_min <- -0.02
#d_max <- 0.035
d_vals <- seq(-0.004,0.004,0.001)
pars <- expand.grid(d = d_vals,rep = 1:sims)
store_names2 <- c("alpha","epsilon","s","maxPd","tmaxPd","minPd","tminPd","delta","seed","attempts")
list_temp <- vector("list",length(sims))
sims_in <- 500

## Simulate
AtoB_par_sens <- function(i){
	sim_seed <- sample.int(.Machine$integer.max, 1)
	set.seed(sim_seed)
	
	d_val <- pars$d[i]
	
	temp <- safe_simumorph(seed = sim_seed,x = amp_pha_cov, m.space = amp_pha_mat, init = which(rownames(amp_pha_mat)=="G2_m_G2"), target = which(rownames(amp_pha_mat)=="G18_m_G18"), method = "AtoB", sim = sims_in, npts = npts, only.shapes = F, a = a_val, e = e_val, f = s_val, max.attempts = 500, speedAtoB = d_val)
	
	AtoB_sens <- temp$result
	df_temp <- c(a_val,e_val,s_val,max(AtoB_sens$P.distances),which.max(AtoB_sens$P.distances),min(AtoB_sens$P.distances),which.min(AtoB_sens$P.distances),d_val,temp$seed_used,temp$attempts)
	names(df_temp) <- store_names2
	
	list_temp[[i]] <- AtoB_sens
	names(list_temp) <- paste0("AtoB_sim_s_", s_val, "_", i)
	return(list(df = df_temp, full = list_temp))
}

## Execute
AtoB_res <- mclapply(1:nrow(pars),AtoB_par_sens,mc.cores = detectCores() - 2)

## Organise
AtoB_sens_df <- do.call(rbind,lapply(AtoB_res, function(x) x$df))
#AtoB_sens_full <- lapply(AtoB_res, function(x) x$full)
AtoB_sens_pdist <- do.call(rbind,lapply(AtoB_res, function(x) x$full[[length(x$full)]]$P.distances))

## Save
saveRDS(AtoB_sens_df,"../SI3_results/AtoB_sens_df.rds")
saveRDS(AtoB_sens_pdist,"../SI3_results/AtoB_sens_pdist.rds")
#saveRDS(AtoB_sens_full,"../SI3_results/AtoB_sens_full.rds")

#saveRDS(AtoA_sens_full,"../SI3_results/AtoA_sens_full.rds")
###########################
###########################
##### Free
###########################
###########################
sims_in <- 500
list_temp <- vector("list",length(sims))

pars <- expand.grid(a = a_vals,s = s_vals,rep = 1:sims)

## Simulate
Free_par_sens <- function(i){
	## Generate reproducible seed for this iteration
	sim_seed <- sample.int(.Machine$integer.max, 1)
	set.seed(sim_seed)
	
	a_val <- pars$a[i]
	s_val <- pars$s[i]
	#a_val <- runif(1, a_min, a_max)
	#s_val <- runif(1, s_min, s_max)
	
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
Free_res <- mclapply(1:nrow(pars),Free_par_sens,mc.cores = detectCores() - 2)

## Organise
Free_sens_df <- do.call(rbind,lapply(Free_res, function(x) x$df))
#Free_sens_full <- lapply(Free_res, function(x) x$full)
Free_sens_pdist <- t(vapply(Free_res, function(x) {x$P_dists}, numeric(500)))

## Save
saveRDS(Free_sens_df,"../SI3_results/Free_sens_df.rds")
#saveRDS(Free_sens_full,"../SI3_results/Free_sens_full.rds")
saveRDS(Free_sens_pdist,"../SI3_results/Free_sens_pdist.rds")




