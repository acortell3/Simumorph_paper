
######################################################
### CREATION OF THE MORPHOSPACE
#######################################################

#######################################################
#### Load the different outlines (produced in inkscape)
######################################################
## Load outlines for GMM analysis and perform PCA to extract morphometric value
library(Momocs) ## For GMM functions
library(matrixcalc) ## Check if the covariance matrix is positive definite

## Load necessary functions
source("trig_functions.R")

## Some utilities
n.harm <- 7 # Number of harmonics
npts <- 120

types <- c("G1.1", "G1.2", "G2", "G3.1", "G3.2", "G4.1",
	   "G4.2", "G5.1", "G5.2", "G6", "G7.1", "G7.2",
	   "G9", "G10", "G11", "G12.1", "G12.2", "G13.1",
	   "G13.2", "G14.1", "G14.2", "G15.1", "G15.2",
	   "G16.1", "G16.2", "G17.1", "G17.2", "G18")
geo_out <- list()
#morpho_pars <- list()

## Extract outlines
dir <- getwd()
setwd("../Parametric_space_outlines/")
set.seed(1)
geo_list <- c(1:28)
geo_list <- paste0(geo_list,".jpg")
geo_out <- Out(import_jpg(geo_list, auto.notcentered = T), fac = types)
setwd(dir)

## Extract parameters from outlines. To avoid flipphing shapes, align before through the x axis (coo_alignxax) and set normalisation to false, as explained in the documentation (it should be controlled for by coo_center(coo_scale(), but for precaution
morpho_pars <- efourier(coo_center(coo_scale(coo_alignxax(geo_out))),nb.h = n.harm, norm = F, start = T)

## Check that shapes are indeed consistently aligned
#to_plot_m <- list()

#for (i in 1:28){
#	morpho_par_list <- list("an" = morpho_pars_vals[i,1:7],
#				"bn" = morpho_pars_vals[i,8:14],
#		      		"cn" = morpho_pars_vals[i,15:21],
#		      		"dn" = morpho_pars_vals[i,22:28],
#				"ao" = 0,
#				"co" = 0)
#	to_plot_m[[i]] <- efourier_i(ef=morpho_par_list, nb.h = n.harm, nb.pts = npts)
#}

#plot(to_plot_m[[1]], type = "l", xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
#for (i in 2:28){
#	lines(to_plot_m[[i]])
#}

## Expand the morphospace. for a symmetric matrix of G1-G28, G_ij is the mean between G_i and G_j
expanded_ms <- list()
mat_index <- 1 ## Need this to populate the matrix
for (i in 1:length(types)){
	for (j in i:length(types)){
		shape <- apply(cbind(morpho_pars[i],morpho_pars[j]),1,mean)
		expanded_ms[[mat_index]] <- shape
		mat_index <- mat_index + 1
	}
}

## Extract parameters without constants
par_vals <- sapply(expanded_ms,function(x) unlist(x))[c(0:length(types)),]

#######################################################
## 2.2 Create covariance matrix for amplitude and phase
#######################################################
## Matrix to store amplitude and phase
amp_pha_mat <- data.frame("Ax1" = rep(0,length(expanded_ms)), "Ay1" = rep(0,length(expanded_ms)), "Ax2" = rep(0,length(expanded_ms)), "Ay2" = rep(0,length(expanded_ms)), "Ax3" = rep(0,length(expanded_ms)), "Ay3" = rep(0,length(expanded_ms)), "Ax4" = rep(0,length(expanded_ms)),"Ay4" = rep(0,length(expanded_ms)),"Ax5" = rep(0,length(expanded_ms)), "Ay5" = rep(0,length(expanded_ms)), "Ax6" = rep(0,length(expanded_ms)), "Ay6" = rep(0,length(expanded_ms)), "Ax7" = rep(0,length(expanded_ms)), "Ay7" = rep(0,length(expanded_ms)),"Phix1" = rep(0,length(expanded_ms)), "Phiy1" = rep(0,length(expanded_ms)), "Phix2" = rep(0,length(expanded_ms)), "Phiy2" = rep(0,length(expanded_ms)), "Phix3" = rep(0,length(expanded_ms)), "Phiy3" = rep(0,length(expanded_ms)), "Phix4" = rep(0,length(expanded_ms)), "Phiy4" = rep(0,length(expanded_ms)), "Phix5" = rep(0,length(expanded_ms)), "Phiy5" = rep(0,length(expanded_ms)), "Phix6" = rep(0,length(expanded_ms)), "Phiy6" = rep(0,length(expanded_ms)), "Phix7" = rep(0,length(expanded_ms)), "Phiy7" = rep(0,length(expanded_ms)))

## Sort matrix per harmonic
harm1 <- t(par_vals[c(1,8,15,22),])
harm2 <- t(par_vals[c(2,9,16,23),])
harm3 <- t(par_vals[c(3,10,17,24),])
harm4 <- t(par_vals[c(4,11,18,25),])
harm5 <- t(par_vals[c(5,12,19,26),])
harm6 <- t(par_vals[c(6,13,20,27),])
harm7 <- t(par_vals[c(7,14,21,28),])

## Convert to list for looping
harms <- list(harm1,harm2,harm3,harm4,harm5,harm6,harm7)

## Populate with amplitudes
index <- 1 ## I need this for correspondence between list and matrix elements
for (i in 1:length(harms)){
	amp_pha_mat[,index] <- amplitude(harms[[i]], coordinate = "x")
	index <- index + 1
	amp_pha_mat[,index] <- amplitude(harms[[i]], coordinate = "y")
	index <- index + 1	
}

## Populate with phases
index <- 15 ## I need this for correspondence between list and matrix elements
for (i in 1:length(harms)){
	amp_pha_mat[,index] <- phase(harms[[i]], coordinate = "x")
	index <- index + 1
	amp_pha_mat[,index] <- phase(harms[[i]], coordinate = "y")
	index <- index + 1	
}

## Create rownames for amp_pha_mat
names_apm <- rep(NA,nrow(amp_pha_mat))

index <- 1
for (i in 1:length(types)){
	for (j in i:length(types)){
		names_apm[index] <- paste0(types[i],"_m_",types[j])
		index <- index + 1
	}
}

rownames(amp_pha_mat) <- names_apm

#### Plot the whole morphospace (to check for unplausible shapes)
## Vector to select shapes (done before now because needed for colouring)
not_selected <- c(20,23,25,47,48,50,51,52,73,76,78,101,122,123,125,127,145,146,148,150,170,188,189,191,193,208,211,213,227,230,232,246,248,250,262,263,265,267,278,288,298,300,307,310,321,323,332,337,344,347,348,353,356,367,369,370,372,375,377,382,388,390,395,398) ## After initial visual inspection of the extended morphospace
col_vec <- rep("darkslategrey",nrow(amp_pha_mat))
col_vec[not_selected] <- "khaki3"

## For indexing
chunks <- list("chunk1" = c(1,100),
	       "chunk2" = c(101,200),
	       "chunk3" = c(201,300),
	       "chunk4" = c(301,400),
	       "chunk5" = c(401,406))
index <- 1 # For colouring

for (j in 1:length(chunks)){
	png(paste0("../Figures/SI_figures/Morphospace_",chunks[[j]][1],"_",chunks[[j]][2],".png"), res = 50, height = 1500, width = 1500)
	par(mfrow = c(10,10))
	for (i in chunks[[j]][1]:chunks[[j]][2]){
		## Prepare shape
		amp_pha <- amp_pha_mat[i,] 
		amp_pha <- unlist(as.data.frame(amp_pha))
		
		## Rebuild parameters from harmonic/phase values
		Ax <- amp_pha[c(1,3,5,7,9,11,13)]
		Ay <- amp_pha[c(2,4,6,8,10,12,14)]
		Phyx <- amp_pha[c(15,17,19,21,23,25,27)]
		Phyy <- amp_pha[c(16,18,20,22,24,26,28)]
		
		## cos and sin inverted to comply with Momocs
		mat_par_vals <- data.frame("an" = Ax*sin(Phyx),
		        		   "bn" = Ax*cos(Phyx),
		     		      	   "cn" = Ay*sin(Phyy),
		     	   	      	   "dn" = Ay*cos(Phyy))
		
		mat_par_vals_l <- list("an" = mat_par_vals[,1],
			               "bn" = mat_par_vals[,2],
			               "cn" = mat_par_vals[,3],
		        	       "dn" = mat_par_vals[,4])
		
		amp_pha_shape <- efourier_i(mat_par_vals_l, nb.h = n.harm, nb.pts = npts)
		amp_pha_shape <- coo_sample(amp_pha_shape, n = npts)

		plot(amp_pha_shape, type = "n", asp = 1, xlab = "", ylab = "", main = rownames(amp_pha_mat)[i])
		polygon(amp_pha_shape, col = col_vec[index], border = "black")	
		index <- index + 1
	}

	dev.off()
}

## Remove not accepted shapes
amp_pha_mat <- amp_pha_mat[-not_selected,]

## Create covariance matrix
amp_pha_cov <- cov(amp_pha_mat)

## Can check whether it's positive definite (it is)
# is.positive.definite(amp_pha_cov)
######################################################
# 2.3 Save the objects produced, to use later in the simulation
###################################################

## Save observed shapes
saveRDS(geo_out,"../Utilities/geo_out.rds")

## Save observed parameters
saveRDS(morpho_pars,"../Utilities/morpho_pars.rds")

## Save amplitude and phases matrix
saveRDS(amp_pha_mat,"../Utilities/amp_pha_mat.rds")

## Save covariance matrix
saveRDS(amp_pha_cov,"../Utilities/amp_pha_cov.rds")




