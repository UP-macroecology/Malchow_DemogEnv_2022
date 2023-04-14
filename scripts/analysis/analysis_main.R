
###---------------------------------------------------------------------------------------
##   GENERATE OUTPUTS AND PLOTS FROM CALIBRATED DEMOGRAPHY-ENVIRONMENT RELATIONSHIPS
###---------------------------------------------------------------------------------------

# Main script for analysing the posterior samples and for plotting results; loads the functions stored in the other files in the same folder

# Contents:
# 1.) Create and store combined MCMCs, extract quantiles of marginals 
# 2.) Compose demography-environment relationships from calibrated parameters
# 3.) Plot demography-environment relationships
# 4.) Plot growth-rate maps
# 5.) Make posterior predictions & Calculate c-Index & Plot abundance time series
# 6.) Climate attribution: a) Total abundance, b) Intrinsic growth rate change


# include packages

library(raster)
library(RangeShiftR)
library(coda)
library(bayesplot)
library(BayesianTools)
library(Hmisc)
library(xtable)
library(ggplot2)
library(pROC)


# read functions from other files in this folder

source("scripts/analysis/combine_chains.R")
source("scripts/analysis/docu_colors.R")
cols_specID <- rev(cols_qualitative)
cols_specID <- rev(cols_qualitative[c(8:10,1:7)])
cols_specID[c(1,2,4)] <- c(NA,cols_specID[1],"#20b56a") #"#20b56a",NA,"#52198a"


# species names
species_names <- c("Chaffinch", 
				   "Bullfinch",
				   "Crested tit",
				   "Eurasian woodcreeper",
				   "Eurasian nuthatch",
				   "Goldcrest",
				   "Dunnock",
				   "Common linnet",
				   "Ring ouzel",
				   "Alpine accentor")

# choose type of graphics device: latex (TRUE) or pdf (FALSE)?
latex_device <- FALSE

# parameter names
par_names_tex <- c("Dens-dep. $b^{-1}$","Fec. $\\rho_0$","Juv. surv. $s_j$","Ad. surv. $s_a$","Emig. prob. $e_1$","Disp. dist. $\\overline d$","Dispersion $\\nu$","\\beta_{fc,p1}","\\beta_{fc,p2}","\\beta_{fc,t1}","\\beta_{fc,t2}","\\beta_{jS,p1}","\\beta_{jS,p2}","\\beta_{jS,t1}","\\beta_{jSt2}","\\beta_{aS,p1}","\\beta_{aS,p2}","\\beta_{aS,t1}")
par_names     <- c("DensDep",           "Fecund",        "juvSurv",         "adSurv",         "EmigProb",         "DispDist",                 "GPsize",           "Fcb.pr1",      "Fcb.pr2",      "Fcb.tm1",      "Fcb.tm2",      "jSb.pr1",      "jSb.pr2",      "jSb.tm1",      "jSb.tm2",     "aSb.pr1",      "aSb.pr2",      "aSb.tm"       )



###--- 1.) Create and store combined MCMCs

# Load and combine single MCMCs, generate some summary statistics, and store everything to new objects

## a) get info about the MCMC to identify the correct file names to read:

samplr <- "CA_mcmc" # likelihood type
iter <- 6e+04       # number of iterations in chain
specs <- spec0 <- 2:10       # species IDs
spec0[9] <- 0               # species IDs as used in batch numbers
batches <- matrix(rep((500+spec0),each=3) + (4:6)*10, nrow=3)					#--- Changed batch nr. 
thinningrate <- 100  # rate of thinning

sel_specs <- c(2,3,4,5,7,8,9,10) # -> 6 not converged



## b) Create combined MCMC files: 

# loop over columns in batches; for each batchID read single MCMC from file and combine MCMCs of all batchIDs in a column; then save to file

if(FALSE){ # deactivate the loop once combined MCMCs are generated				#- ran once
	for(specID in sel_specs ) {
	#for(specID in 1:ncol(batches) ) {
		
		simnames <- batches[,which(specs==specID)]
		simnames <- simnames[!is.na(simnames)]
		print(specID)
		
		# read the three independent MCMCs and create combined MCMC sampler list
		sampler_list <- combineMCMCs(samplr,specID,simnames,iter)
		
		# generate plots
		if(TRUE){
			# open graphics device
			if(latex_device){ tikz(file = paste0(path_doc,"CA_comb_Spec",specID,"_it",iter,"_",simnames[1],"-",simnames[length(simnames)],".tex"), width = 5.76, height = 3.84, standAlone = FALSE)
			}else{ pdf(file = paste0("results/CA_out/CA_comb_Spec",specID,"_it",iter,"_",simnames[1],"-",simnames[length(simnames)],".pdf"),width = 12, height = 8) }
			# plot trace and marginal posteriors
			plot(sampler_list, trace = !latex_device, density = TRUE, smooth = FALSE)
			dev.off()
		}
		
		# save combined chain and info in R object, then remove
		save(sampler_list, #sampler_info, posterior_df, correlations, 
			 file = paste0("results/CA_out/CA_comb_Spec",specID,"_it",iter,"_",simnames[1],"-",simnames[length(simnames)],".Rdata"))
		
		rm(sampler_list, simnames, specID) #, sampler_info, posterior_df, correlations
	}
}


## c) Extract posterior quantiles

calc_quantiles <- c(0.10,0.20,0.25,0.50,0.75,0.80,0.90,0.025,0.975) 						#--- added additional quantiles here

if(FALSE){																		#--- ran once
	quants <- lapply(sel_specs, function(specID){ # 1:ncol(batches)
		
		simnames <- batches[,which(specs==specID)]
		simnames <- simnames[!is.na(simnames)]
	
		# read the combined MCMC
		load( file = paste0("results/CA_out/CA_comb_Spec",specID,"_it",iter,"_",simnames[1],"-",simnames[length(simnames)],".Rdata") )
		
		# convert to BayesionTools object
		BT_smplst <- convertCoda(sampler_list, names = par_names)
		
		# discard burn-in (= first half of chain)
		info_sample <- getSample(BT_smplst, parametersOnly = T, coda = T, start = 1+iter/2) #-- changed to 2 as we have 3x60.000 iterations
		npar <- ncol(info_sample[[1]])
		
		# calculate quartiles
		spec_quartiles <- sapply(1:npar, function(col){
			# extract samples of one parameter only and unlist to combine sub-chains
			par_smpl <- unlist(info_sample[,col])
			# calc quartiles
			quantile(par_smpl, probs = calc_quantiles )
		})
		colnames(spec_quartiles) <- colnames(info_sample[[1]])
		
		# remove loaded object
		rm(sampler_list)
		
		# return
		spec_quartiles
	})
	
	# assign species names
	names(quants) <- species_names[sel_specs]
	
	# save quantiles
	save(quants, file = paste0("results/analysis/species_post_quantiles.Rdata"))
	
}else{
	# load quantiles
	load(paste0("results/analysis/species_post_quantiles.Rdata"))
}


###--- 2.) Compose demography-environment relationships from calibrated parameters for each species


## a) Fecundity

loc_dens <- 1

DER_Fecund <- function(Fec0, Fc.pr1, Fc.pr2, Fc.tm1, Fc.tm2, DDp, Tbr=0, Pbr=0, Hab=100){ 
	exp(
		log(Fec0) + 
		Fc.pr1*Pbr + 
		Fc.pr2*Pbr^2 + 
		Fc.tm1*Tbr + 
		Fc.tm2*Tbr^2 ) *
	# density-dependence
	exp( - 1 * eval(loc_dens) / DDp / Hab )
}


## b) Juvenile survival

DER_juvSurv <- function(jSv0, jS.pr1, jS.pr2, jS.tm1, jS.tm2, Tat=0, Pwn=0){
	plogis(
		jS.pr1*Pwn + 
		jS.pr2*Pwn^2 + 
		jS.tm1*Tat + 
		jS.tm2*Tat^2 ,
		location = log((1-jSv0)/jSv0))
}


## c) Adult survival

DER_adSurv <- function(aSv0, aS.pr1, aS.pr2, aS.tm1, Twn=0, Pwn=0){
	plogis(
		aS.pr1*Pwn + 
		aS.pr2*Pwn^2 + 
		aS.tm1*Twn ,
		location = log((1-aSv0)/aSv0))
}




if(FALSE){																		#--- ran once with CCcorr TRUE and once with CCcorr FALSE

# Load species climate space quantiles
load("results/analysis/species_preds_quantiles.Rdata")							#-- copied from archive folder

use_quantiles <- calc_quantiles[c(1,3,4,5,7,8,9)]                               #--- added additional quantiles here
use_climquants <- paste0(calc_quantiles[c(1,3,4,5,7)]*100,"%")					#--- we use 0.1 and 0.9 
post_samplesize <- 396															#--- sample size to 396
pred_resolution <- 0.01

CCcorr <- FALSE 																
if(CCcorr){
	load("data/climatemaps/aggrVariables_zTrafod_CCcorrected.Rdata")
	climvar_chg <- clim_means_lm["year",]*20
}

	# loop over species
	sp_DERs <- lapply(sel_specs, function(specID) {
		
		name <- names(spec_preds_quants)[specID]
		cat(c(specID,name,"\n"))
		
		simnames <- batches[,which(specs==specID)]
		simnames <- simnames[!is.na(simnames)]
		
		# read the combined MCMC
		load( file = paste0("results/CA_out/CA_comb_Spec",specID,"_it",iter,"_",simnames[1],"-",simnames[length(simnames)],".Rdata") )
		
		# convert to BayesionTools object
		BT_smplst <- convertCoda(sampler_list, names = par_names)
		
		# remove loaded object again
		rm(sampler_list)
		
		# select needed parameters, discard burn-in (= first half of chain) & draw posterior sample
		parSel <- c(1:4,8:18)
		post_sample <- getSample(BT_smplst, parametersOnly = T, coda = F, start = 1+iter/2, whichParameters = parSel, numSamples = post_samplesize) #-- changed from 3 to 2 as we now have 3*60.000 iterations
		# predictor quantiles
		clim_quants <- spec_preds_quants[[specID]]
		
		# generate predictor axes
		if(CCcorr){
			x_Tbr <- c(clim_quants["50%","breeding_tmp"], clim_quants["50%","breeding_tmp"] + climvar_chg[1] )
			x_Pbr <- c(clim_quants["50%","breeding_pr"],  clim_quants["50%","breeding_pr"]  + climvar_chg[2] )
			x_Tat <- c(clim_quants["50%","postbreed_tmp"],clim_quants["50%","postbreed_tmp"]+ climvar_chg[3] )
			x_Twn <- c(clim_quants["50%","winter_tmin"],  clim_quants["50%","winter_tmin"]  + climvar_chg[4] )
			x_Pwn <- c(clim_quants["50%","winter_pr"],    clim_quants["50%","winter_pr"]    + climvar_chg[5] )
		}else{
			x_Tbr <- seq(clim_quants["10%","breeding_tmp"], clim_quants["90%","breeding_tmp"], pred_resolution)
			x_Pbr <- seq(clim_quants["10%","breeding_pr"],  clim_quants["90%","breeding_pr"],  pred_resolution)
			x_Tat <- seq(clim_quants["10%","postbreed_tmp"],clim_quants["90%","postbreed_tmp"],pred_resolution)
			x_Twn <- seq(clim_quants["10%","winter_tmin"],  clim_quants["90%","winter_tmin"],  pred_resolution)
			x_Pwn <- seq(clim_quants["10%","winter_pr"],    clim_quants["90%","winter_pr"],    pred_resolution)
		}
		x_Hab <- seq(clim_quants["10%","hab_suit"],     clim_quants["90%","hab_suit"],     pred_resolution)
		
	
		
		### Loop though posterior samples and calculate species demog.env.rels
		
		## a) Fecundity
		
		DER_Fc_Tbr <- lapply(use_climquants, function(climQ){
			
			post_pred <- sapply(1:post_samplesize, function(s){
				
				smp <- post_sample[s,]
				
				DER_Fecund(Fec0 = smp["Fecund"],
						   Fc.pr1 = smp["Fcb.pr1"],
						   Fc.pr2 = smp["Fcb.pr2"],
						   Fc.tm1 = smp["Fcb.tm1"],
						   Fc.tm2 = smp["Fcb.tm2"],
						   DDp  = smp["DensDep"],
						   Tbr = x_Tbr,
						   Pbr = clim_quants[climQ,"breeding_pr"],
						   Hab = clim_quants["50%","hab_suit"])
			})
			apply(post_pred, 1, quantile, probs = use_quantiles )
		})
		names(DER_Fc_Tbr) <- paste0("pred2_",use_climquants)
		
		
		DER_Fc_Pbr <- lapply(use_climquants, function(climQ){
			
			post_pred <- sapply(1:post_samplesize, function(s){
				
				smp <- post_sample[s,]
				
				DER_Fecund(Fec0 = smp["Fecund"],
						   Fc.pr1 = smp["Fcb.pr1"],
						   Fc.pr2 = smp["Fcb.pr2"],
						   Fc.tm1 = smp["Fcb.tm1"],
						   Fc.tm2 = smp["Fcb.tm2"],
						   DDp  = smp["DensDep"],
						   Tbr = clim_quants[climQ,"breeding_tmp"],
						   Pbr = x_Pbr,
						   Hab = clim_quants["50%","hab_suit"])
			})
			apply(post_pred, 1, quantile, probs = use_quantiles )
		})
		names(DER_Fc_Pbr) <- paste0("pred2_",use_climquants)
		
		
		
		## b) Juvenile survival
		
		DER_jS_Tat <- lapply(use_climquants, function(climQ){
			
			post_pred <- sapply(1:post_samplesize, function(s){
				
				smp <- post_sample[s,]
				
				DER_juvSurv(jSv0 = smp["juvSurv"], 
							jS.pr1 = smp["jSb.pr1"],
							jS.pr2 = smp["jSb.pr2"],
							jS.tm1 = smp["jSb.tm1"],
							jS.tm2 = smp["jSb.tm2"],
							Tat = x_Tat, 
							Pwn = clim_quants[climQ,"winter_pr"])
			})
			apply(post_pred, 1, quantile, probs = use_quantiles )
		})
		names(DER_jS_Tat) <- paste0("pred2_",use_climquants)
		
		
		DER_jS_Pwn <- lapply(use_climquants, function(climQ){
			
			post_pred <- sapply(1:post_samplesize, function(s){
				
				smp <- post_sample[s,]
				
				DER_juvSurv(jSv0 = smp["juvSurv"],
							jS.pr1 = smp["jSb.pr1"],
							jS.pr2 = smp["jSb.pr2"],
							jS.tm1 = smp["jSb.tm1"],
							jS.tm2 = smp["jSb.tm2"],
							Tat = clim_quants[climQ,"postbreed_tmp"],
							Pwn = x_Pwn)
			})
			apply(post_pred, 1, quantile, probs = use_quantiles )
		})
		names(DER_jS_Pwn) <- paste0("pred2_",use_climquants)
		
		
		
		## c) Adult survival
		
		DER_aS_Twn <- lapply(use_climquants, function(climQ){
			
			post_pred <- sapply(1:post_samplesize, function(s){
				
				smp <- post_sample[s,]
				
				DER_adSurv(aSv0 = smp["adSurv"],
						   aS.pr1 = smp["aSb.pr1"],
						   aS.pr2 = smp["aSb.pr2"],
						   aS.tm1 = smp["aSb.tm"],
						   Twn = x_Twn,
						   Pwn = clim_quants[climQ,"winter_pr"])
			})
			apply(post_pred, 1, quantile, probs = use_quantiles )
		})
		names(DER_aS_Twn) <- paste0("pred2_",use_climquants)
		
		
		DER_aS_Pwn <- lapply(use_climquants, function(climQ){
			
			post_pred <- sapply(1:post_samplesize, function(s){
				
				smp <- post_sample[s,]
				
				DER_adSurv(aSv0 = smp["adSurv"],
						   aS.pr1 = smp["aSb.pr1"],
						   aS.pr2 = smp["aSb.pr2"],
						   aS.tm1 = smp["aSb.tm"],
						   Twn = clim_quants[climQ,"winter_tmin"],
						   Pwn = x_Pwn)
			})
			apply(post_pred, 1, quantile, probs = use_quantiles )
		})
		names(DER_aS_Pwn) <- paste0("pred2_",use_climquants)
		
		
		# Return
		list(DER_Fc_Tbr = DER_Fc_Tbr,
			 DER_Fc_Pbr = DER_Fc_Pbr,
			 DER_jS_Tat = DER_jS_Tat,
			 DER_jS_Pwn = DER_jS_Pwn,
			 DER_aS_Twn = DER_aS_Twn,
			 DER_aS_Pwn = DER_aS_Pwn,
			 preds = list(x_Tbr = x_Tbr,
			 			  x_Pbr = x_Pbr,
			 			  x_Tat = x_Tat,
			 			  x_Twn = x_Twn,
			 			  x_Pwn = x_Pwn,
			 			  x_Hab = ifelse(CCcorr,NA,x_Hab))
			 )
	
	# End of species loop	
	})
	names(sp_DERs) <- species_names[sel_specs]
	
	
	# save DERs
	save(sp_DERs, file = paste0("results/analysis/sp_DERs_largerCI",ifelse(CCcorr,"_CCcorr",""),".Rdata"))

}else{
	# load DERs
	load(paste0("results/analysis/sp_DERs_largerCI.Rdata"))
}





###--- 3.) Plot demography-environment relationships

# if(TRUE){																		#--- ran once
if(FALSE){
	# use_quantiles <- calc_quantiles[c(1,4,7)]									#--- not used here -> quantile to be used in plot are selected with variable DER_quants later
	use_climquants <- paste0(calc_quantiles[c(1,4,7)]*100,"%")					#--- only used for species specific plots
		
	# Plotting options
	
	col_shade <- c(t_col(cols_qualitative[10],percent = 50),t_col(cols_qualitative[3],percent = 50))
	col_shade[3] <- col_shade[1]
	
	col_line <- c(makeColDarker(cols_qualitative[10]),makeColDarker(cols_qualitative[3]))
	col_line[3] <- col_line[1]
	
	xlimits <- list(x_Tbr = c(-1.6,1.1), x_Pbr = c(-1.2,2.2), x_Tat = c(-1.7,1.1), x_Twn = c(-1.7,1.1), x_Pwn = c(-1.1,2.1) , x_Hab = c(0,100) )
	fec_max <- 7
	linetype <- c(2,1,2)
	
	
	##-- 3.1.) Per species, all DERs individually
	
	# Loop over species															#-- did not change anything as plots are not used in manuscript
	
	for(specID in sel_specs){
		
		sp_name <- species_names[specID]
		cat(c(specID, sp_name, "\n"))
		
		# get stored demogr.-env. rel.s
		this_DERs <- getElement(sp_DERs, sp_name)
	
		
		##-- a) Fecundity
		
		# - over breeding season precipitation
		png(file = paste0("results/analysis/demogenv_rels/posterior_preds/Fec_Pbr_",specID,".png"), width = 720, height = 480)
			plot(NULL, xlim = xlimits$x_Pbr , ylim = c(0,fec_max) ,  main = sp_name, xlab = "Breeding season precipitation" , ylab = "Fecundity") 
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_Fc_Pbr, clqu)
				polygon(x = c(this_DERs$preds$x_Pbr,rev(this_DERs$preds$x_Pbr)), y = c(this_resp["25%",],rev(this_resp["75%",])),  col = col_shade[i], border = NA)
			}
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_Fc_Pbr, clqu)
				lines(x = this_DERs$preds$x_Pbr, y = this_resp["50%",],  col = col_line[i], lwd = 2, lty = linetype[i] )
			}
		dev.off()
		
		# - over breeding season temperature
		png(file = paste0("results/analysis/demogenv_rels/posterior_preds/Fec_Tbr_",specID,".png"), width = 720, height = 480)
			plot(NULL, xlim = xlimits$x_Tbr , ylim = c(0,fec_max) ,  main = sp_name, xlab = "Breeding season temperature" , ylab = "Fecundity") 
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_Fc_Tbr, clqu)
				polygon(x = c(this_DERs$preds$x_Tbr,rev(this_DERs$preds$x_Tbr)), y = c(this_resp["25%",],rev(this_resp["75%",])),  col = col_shade[i], border = NA)
			}
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_Fc_Tbr, clqu)
				lines(x = this_DERs$preds$x_Tbr, y = this_resp["50%",],  col = col_line[i], lwd = 2, lty = linetype[i] )
			}
		dev.off()
		
		
		##-- b) Juvenile survival
		
		# - over winter precipitation
		png(file = paste0("results/analysis/demogenv_rels/posterior_preds/jSv_Pwn_",specID,".png"), width = 720, height = 480)
			plot(NULL, xlim = xlimits$x_Pwn , ylim = c(0,1) ,  main = sp_name, xlab = "Winter precipitation" , ylab = "Juvenile survival") 
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_jS_Pwn, clqu)
				polygon(x = c(this_DERs$preds$x_Pwn,rev(this_DERs$preds$x_Pwn)), y = c(this_resp["25%",],rev(this_resp["75%",])),  col = col_shade[i], border = NA)
			}
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_jS_Pwn, clqu)
				lines(x = this_DERs$preds$x_Pwn, y = this_resp["50%",],  col = col_line[i], lwd = 2, lty = linetype[i] )
			}
		dev.off()
		
		# - over autumn temperature
		png(file = paste0("results/analysis/demogenv_rels/posterior_preds/jSv_Tat_",specID,".png"), width = 720, height = 480)
			plot(NULL, xlim = xlimits$x_Tat , ylim = c(0,1) ,  main = sp_name, xlab = "Autumn temperature" , ylab = "Juvenile survival") 
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_jS_Tat, clqu)
				polygon(x = c(this_DERs$preds$x_Tat,rev(this_DERs$preds$x_Tat)), y = c(this_resp["25%",],rev(this_resp["75%",])),  col = col_shade[i], border = NA)
			}
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_jS_Tat, clqu)
				lines(x = this_DERs$preds$x_Tat, y = this_resp["50%",],  col = col_line[i], lwd = 2, lty = linetype[i] )
			}
		dev.off()
		
		
		
		##-- c) Adult survival 
		
		# - over winter precipitation
		png(file = paste0("results/analysis/demogenv_rels/posterior_preds/aSv_Pwn_",specID,".png"), width = 720, height = 480)
			plot(NULL, xlim = xlimits$x_Pwn , ylim = c(0,1) ,  main = sp_name, xlab = "Winter precipitation" , ylab = "Adult survival") 
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_aS_Pwn, clqu)
				polygon(x = c(this_DERs$preds$x_Pwn,rev(this_DERs$preds$x_Pwn)), y = c(this_resp["25%",],rev(this_resp["75%",])),  col = col_shade[i], border = NA)
			}
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_aS_Pwn, clqu)
				lines(x = this_DERs$preds$x_Pwn, y = this_resp["50%",],  col = col_line[i], lwd = 2, lty = linetype[i] )
			}
		dev.off()
		
		# - over winter min-temperature
		png(file = paste0("results/analysis/demogenv_rels/posterior_preds/aSv_Twn_",specID,".png"), width = 720, height = 480)
			plot(NULL, xlim = xlimits$x_Twn , ylim = c(0,1) ,  main = sp_name, xlab = "Winter temperature" , ylab = "Adult survival") 
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_aS_Twn, clqu)
				polygon(x = c(this_DERs$preds$x_Twn,rev(this_DERs$preds$x_Twn)), y = c(this_resp["25%",],rev(this_resp["75%",])),  col = col_shade[i], border = NA)
			}
			for(i in c(1,3,2) ){
				clqu <- paste0("pred2_",use_climquants[i])
				this_resp <- getElement(this_DERs$DER_aS_Twn, clqu)
				lines(x = this_DERs$preds$x_Twn, y = this_resp["50%",],  col = col_line[i], lwd = 2, lty = linetype[i] )
			}
		dev.off()
	
	# end species loop
	}
	
	##-- 3.2.) Per demog.-env. relation, all species in one plot
	
	DER_quants <- c("2.5%","10%", "90%", "97.5%")	                            #-- adapted index when generating plots (once 80% CI and once 95% CI)
	
	cols <- cols_specID
	col_shade <- sapply(cols,t_col,percent = 50)
	col_shade1 <- sapply(cols,t_col,percent = 25)                               #-- added lighter shade for pultipanel plot
	col_line <- sapply(cols,makeColDarker)
	
	# use only median for second predictor (not three different levels as in the individual plots)
	clqu <- "pred2_50%"
	fec_max <- 4
	
	# load climate means and SDs to back-z-transform 
	load(file = "data/climatemaps/zTrafo.Rdata")
	
	zBacktrafo <- function(x, var){ 
		out <- (x*envVars_tot_stdev[var])+envVars_tot_means[var]
		if(substr(var, 3, 3) == "T") out <- out-273.15
		if(substr(var, 3, 3) == "P") out <- out/100
		out
	}
	
	names(envVars_tot_means) <- names(envVars_tot_stdev) <- c("x_Tbr", "x_Pbr", "x_Tat", "x_Twn", "x_Pwn") #-- needed to make adaptions to run the code (names were not correct)
	
	# back-trafo limits
	xlimits <- list(x_Tbr = zBacktrafo(c(-1.6,1.1), "x_Tbr"), 
					x_Pbr = zBacktrafo(c(-1.2,1.7), "x_Pbr"), 
					x_Tat = zBacktrafo(c(-1.7,1.1), "x_Tat"), 
					x_Twn = zBacktrafo(c(-1.7,1.1), "x_Twn"), 
					x_Pwn = zBacktrafo(c(-1.1,2.1), "x_Pwn"), 
					x_Hab = c(0,100) )
	
	names(envVars_tot_means) <- names(envVars_tot_stdev) <- names(xlimits)[-6]
	
	
	## MULTIPLOT
	
	# tikz(file = "results/analysis/demogenv_rels/posterior_preds/DERs_all_noZ_largeCI.tex", width = 5.76, height = 4.84, standAlone = FALSE) #, width = 3.84, height = 3.22
	# pdf(file = "results/analysis/demogenv_rels/posterior_preds/DERs_all_noZ_largeCI_95.pdf", width = 5.76, height = 4.84) 
	pdf(file = "results/analysis/demogenv_rels/posterior_preds/DERs_all_noZ_largeCI_80.pdf", width = 5.76, height = 4.84) #-- adapt name according to CI; used grapic device for first figures
	{
		
		par(mfcol=c(3,2), cex = 0.9, mgp = c(1.7,0.6,0), mar=c(3.4,3.2,0.1,0.1))
		
		## Column 1
		
		# - Fecundity over breeding season temperature
		plot(NULL, xlim = xlimits$x_Tbr , ylim = c(0,fec_max) ,  main = "", xlab = "Breeding season temperature [$^\\circ$ C]" , ylab = "Fecundity $\\rho$", yaxp = c(0,fec_max,2)) 
		for(specID in sel_specs){
			print(specID)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_Fc_Tbr, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Tbr, "x_Tbr")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID] #-- adapted quaniles
					, border = NA)
			lines(x = this_pred , y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
		}
		
		# - Juvenile survival over autumn temperature
		plot(NULL, xlim = xlimits$x_Tat , ylim = c(0,1) ,  main = "", xlab = "Autumn temperature [$^\\circ$ C]" , ylab = "Juvenile survival $s_j$", yaxp = c(0,1,2))
		for(specID in sel_specs){
			print(specID)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_jS_Tat, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Tat, "x_Tat")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID], border = NA)
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
		}
		
		# - Adult survival over winter min-temperature
		plot(NULL, xlim = xlimits$x_Twn , ylim = c(0,1) ,  main = "", xlab = "Winter minimum temperature [$^\\circ$ C]" , ylab = "Adult survival $s_a$", yaxp = c(0,1,2))
		for(specID in sel_specs){
			print(specID)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_aS_Twn, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Twn, "x_Twn")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID], border = NA) #-- adapted quaniles
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
		}
		
		## Column 2
		par(mar=c(3.4,0.2,0.1,0.2)) #-- adapted margins
		
		# - Fecundity over breeding season precipitation
		plot(NULL, xlim = xlimits$x_Pbr , ylim = c(0,fec_max) ,  main = "", xlab = "Breeding season precipitation [mm]", yaxt = "n") 
		axis(1, at = NULL, labels = FALSE)
		for(specID in sel_specs){
			print(specID)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_Fc_Pbr, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Pbr, "x_Pbr")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID]  #-- adapted quaniles
					, border = NA)
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
		}
		
		# - Juvenile survival over winter precipitation
		plot(NULL, xlim = xlimits$x_Pwn , ylim = c(0,1) ,  main = "", xlab = "Winter precipitation [mm]", yaxt = "n")
		for(specID in sel_specs){
			print(specID)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_jS_Pwn, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Pwn, "x_Pwn")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID], border = NA)
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
		}
		
		# - Adult survival over winter precipitation
		plot(NULL, xlim = xlimits$x_Pwn , ylim = c(0,1) ,  main = "", xlab = "Winter precipitation [mm]", yaxt = "n")
		for(specID in sel_specs){
			print(specID)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_aS_Pwn, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Pwn, "x_Pwn")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID], border = NA) #-- adapted quaniles
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
		}
		
		# end multiplot
	}
	dev.off()
	
	# generate multiplot 6x8
	# tikz(file = "results/analysis/demogenv_rels/posterior_preds/DERs_all_noZ_largeCI.tex", width = 5.76, height = 4.84, standAlone = FALSE) #, width = 3.84, height = 3.22
	species_codes <- c("Bf",
					   "Ct",
					   "Tc",
					   "Nh",
					   "Du",
					   "Li",
					   "Ro",
					   "Aa")
	
	pdf(file = "results/analysis/demogenv_rels/posterior_preds/DERs_all_noZ_largeCI_80+95_6x8.pdf", width = 5.76, height = 4.84*2) # currently this graphic device is used
	{
		# Layout of plot
		layout.matrix <- matrix(c(
			55,56,57,58,59,60,61,62,
			1, 7, 13, 19, 25, 31, 37, 43,
			rep(49,8),
			2, 8, 14, 20, 26, 32, 38, 44,
			rep(50,8),
			3, 9, 15, 21, 27, 33, 39, 45,
			rep(51,8),
			4, 10, 16, 22, 28, 34, 40, 46, 
			rep(52,8),
			5, 11, 17, 23, 29, 35, 41, 47, 
			rep(53,8),
			6, 12, 18, 24, 30, 36, 42, 48, 
			rep(54,8)
		), nrow = 13, ncol = 8, byrow = T)
		# ), nrow = 12, ncol = 8, byrow = T)
		
		layout(mat = layout.matrix,
			   heights = c(1,4,1,4,1,4,1,4,1,4,1,4,1), # Heights of the two rows
			   # heights = c(4,1,4,1,4,1,4,1,4,1,4,1), # Heights of the two rows
			   widths = c(8,5,5,5,5,5,5,5)) # Widths of the two columns

		col_shade1 <- sapply(cols,t_col,percent = 75)
		
		# par(mfcol=c(6, length(sel_specs)), cex = 0.9, mgp = c(1.7,0.6,0), mar=c(2,3.4,0.2,0.1))
		par(cex = 0.8, cex.axis=0.7, mgp = c(1.7,0.6,0), mar=c(2,3.4,0.2,0.1))
		for (specID in sel_specs){
			## row 1	
			if(specID==sel_specs[1]) par(mar=c(1.5,3.4,0.2,0.1))
			if(specID==sel_specs[8]) par(mar=c(1.5,0.2,0.2,0.4)) 
			if(specID %in% sel_specs[-c(1,8)]) par(mar=c(1.5,0.2,0.2,0.1))
			# - Fecundity over breeding season temperature
			plot(NULL, xlim = xlimits$x_Tbr , ylim = c(0,fec_max) ,  main = "", xlab="", ylab = ifelse(specID==sel_specs[1],"Fecundity $\\rho$",""), yaxp = c(0,fec_max,2), yaxt='n', xaxt='n') 
			print(specID)
			if(specID==sel_specs[1]) axis(2, at = NULL, labels = TRUE)
			axis(1, at = c(4,8,12), labels = TRUE)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_Fc_Tbr, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Tbr, "x_Tbr")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[1],],rev(this_resp[DER_quants[4],])), col = col_shade1[specID]
					, border = NA)
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID]
					, border = NA)
			lines(x = this_pred , y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
			
			## row 2
			# - Fecundity over breeding season precipitation
			plot(NULL, xlim = xlimits$x_Pbr , ylim = c(0,fec_max) ,  main = "", xlab="", ylab = ifelse(specID==sel_specs[1],"Fecundity $\\rho$",""), yaxp = c(0,fec_max,2), yaxt='n') 
			print(specID)
			if(specID==sel_specs[1]) axis(2, at = NULL, labels = TRUE)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_Fc_Pbr, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Pbr, "x_Pbr")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[1],],rev(this_resp[DER_quants[4],])), col = col_shade1[specID]
					, border = NA)
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID]
					, border = NA)
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
			
			## row 3
			# - Juvenile survival over autumn temperature
			plot(NULL, xlim = xlimits$x_Tat , ylim = c(0,1) ,  main = "", xlab="", ylab = ifelse(specID==sel_specs[1],"Juvenile survival $s_j$",""), yaxp = c(0,1,2), xaxt='n', yaxt='n')
			print(specID)
			if(specID==sel_specs[1]) axis(2, at = NULL, labels = TRUE)
			axis(1, at = c(0,4,8), labels = TRUE)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_jS_Tat, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Tat, "x_Tat")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[1],],rev(this_resp[DER_quants[4],])), col = col_shade1[specID], border = NA)
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID], border = NA)
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
			
			## row 4
			# if(specID==sel_specs[8]) par(mar=c(4,3.4,0.2,0.1)) else par(mar=c(1,3.4,0.2,0.1))		
			# - Juvenile survival over winter precipitation
			plot(NULL, xlim = xlimits$x_Pwn , ylim = c(0,1) ,  main = "", xlab="", ylab = ifelse(specID==sel_specs[1],"Juvenile survival $s_j$",""), yaxp = c(0,1,2), yaxt='n')
			print(specID)
			if(specID==sel_specs[1]) axis(2, at = NULL, labels = TRUE)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_jS_Pwn, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Pwn, "x_Pwn")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[1],],rev(this_resp[DER_quants[4],])), col = col_shade1[specID], border = NA)
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID], border = NA)
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
			
			## row 5
			# if(specID==sel_specs[8]) par(mar=c(4,3.4,0.2,0.1)) else par(mar=c(1,3.4,0.2,0.1))
			# - Adult survival over winter min-temperature
			plot(NULL, xlim = xlimits$x_Twn , ylim = c(0,1) ,  main = "", xlab="", ylab = ifelse(specID==sel_specs[1],"Adult survival $s_a$",""), yaxp = c(0,1,2), xaxt='n', yaxt='n')
			print(specID)
			if(specID==sel_specs[1]) axis(2, at = NULL, labels = TRUE)
			axis(1, at = c(-15,-10,-5), labels = c(-15,"",-5))
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_aS_Twn, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Twn, "x_Twn")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[1],],rev(this_resp[DER_quants[4],])), col = col_shade1[specID], border = NA)
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID], border = NA)
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
			
			## row 6
			# if(specID==sel_specs[8]) par(mar=c(4,3.4,0.2,0.1)) else par(mar=c(1,3.4,0.2,0.1))
			# - Adult survival over winter precipitation
			plot(NULL, xlim = xlimits$x_Pwn , ylim = c(0,1) ,  main = "", xlab="", ylab = ifelse(specID==sel_specs[1],"Adult survival $s_a$",""), yaxp = c(0,1,2), yaxt='n')
			print(specID)
			if(specID==sel_specs[1]) axis(2, at = NULL, labels = TRUE)
			sp_name <- species_names[specID]
			this_DERs <- getElement(sp_DERs, sp_name)
			this_resp <- getElement(this_DERs$DER_aS_Pwn, clqu)
			this_pred <- zBacktrafo(this_DERs$preds$x_Pwn, "x_Pwn")
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[1],],rev(this_resp[DER_quants[4],])), col = col_shade1[specID], border = NA)
			polygon(x = c(this_pred,rev(this_pred)), y = c(this_resp[DER_quants[2],],rev(this_resp[DER_quants[3],])), col = col_shade[specID], border = NA)
			lines(x = this_pred, y = this_resp["50%",],  col = col_line[specID], lwd = 2, lty = 1)
		}
		par(mar=c(0.1,0.1,0.1,0.1), cex=0.9)
		# axis labels
		plot.new()
		text(x=0.5,y=0.5, label="Breeding season temperature [$^\\circ$ C]")
		plot.new()
		text(x=0.5,y=0.5, label="Breeding season precipitation [mm]")
		plot.new()
		text(x=0.5,y=0.5, label="Autumn temperature [$^\\circ$ C]")
		plot.new()
		text(x=0.5,y=0.5, label="Winter precipitation [mm]")
		plot.new()
		text(x=0.5,y=0.5, label="Winter minimum temperature [$^\\circ$ C]")
		plot.new()
		text(x=0.5,y=0.5, label="Winter precipitation [mm]")
		par(mar=c(0.1,0.1,0.1,0.1), cex=0.9)
		# # species names
		plot.new()
		text(x=0.8,y=0.5, label=species_codes[1])
		plot.new()
		text(x=0.5,y=0.5, label=species_codes[2])
		plot.new()
		text(x=0.5,y=0.5, label=species_codes[3])
		plot.new()
		text(x=0.5,y=0.5, label=species_codes[4])
		plot.new()
		text(x=0.5,y=0.5, label=species_codes[5])
		plot.new()
		text(x=0.5,y=0.5, label=species_codes[6])
		plot.new()
		text(x=0.5,y=0.5, label=species_codes[7])
		plot.new()
		text(x=0.5,y=0.5, label=species_codes[8])
		
		
		# end multiplot
	}
	dev.off()
	

}


###--- 4.) Plot growth-rate maps


if(FALSE){                                                                      #-- Ran once + plot_ratesMaps.R

# read maps of climate and habitat suitability
load("data/climatemaps/aggrVariables_zTrafod.Rdata")
load("data/habitatmaps/habitatmaps_masked.Rdata")

# Swiss border for masking
bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')

# Make template stack for output
empty_map <- yearly_envStacks_standard$year_1999$breeding_tmp
values(empty_map)[!is.na(values(empty_map))] <- 0

growthr_stack_empty <- stack(replicate(4,empty_map))
names(growthr_stack_empty) <- c("fec","jSv","aSv","rInf")

# Object to store time series of spatial mean of growth rate
growth_TS <- matrix(NA, ncol = 10, nrow = 21)

# Loop over species
for(specID in sel_specs){
	
	sp_name_eng <- species_names[specID]
	sp_name_lat <- names(habitatmaps_1km$corine2006)[specID]
	cat(c(specID, sp_name_eng, "\n"))
	
	# get stored posterior medians
	this_median <- getElement(quants,sp_name_eng)["50%",]
	
	# loop over years
	for(year in 1999:2018){
		print(year)
		
		growthr_stack <- growthr_stack_empty
		
		this_habitat <- habitatmaps_1km[[which(names(habitatmaps_1km) == paste0("corine",year))]]
		
		pred_maps <- stack(yearly_envStacks_standard[[(year-1998)]], hab = this_habitat[[which(names(this_habitat) == sp_name_lat )]] )
		names(pred_maps)[6] <- "habitat"
		
		
		# Fecundity
		values(growthr_stack$fec) <- DER_Fecund(Fec0   = this_median["Fecund"],
												Fc.pr1 = this_median["Fcb.pr1"],
												Fc.pr2 = this_median["Fcb.pr2"],
												Fc.tm1 = this_median["Fcb.tm1"],
												Fc.tm2 = this_median["Fcb.tm2"],
												DDp  = this_median["DensDep"],
												Tbr = values(pred_maps$breeding_tmp),
												Pbr = values(pred_maps$breeding_pr),
												Hab = values(pred_maps$habitat))
		
		# Juvenile survival
		values(growthr_stack$jSv) <- DER_juvSurv(jSv0   = this_median["juvSurv"],
												 jS.pr1 = this_median["jSb.pr1"],
												 jS.pr2 = this_median["jSb.pr2"],
												 jS.tm1 = this_median["jSb.tm1"],
												 jS.tm2 = this_median["jSb.tm2"],
												 Tat = values(pred_maps$postbreed_tmp),
												 Pwn = values(pred_maps$winter_pr))
		
		# Adult survival
		values(growthr_stack$aSv) <- DER_adSurv(aSv0   = this_median["adSurv"],
												aS.pr1 = this_median["aSb.pr1"],
												aS.pr2 = this_median["aSb.pr2"],
												aS.tm1 = this_median["aSb.tm"],
												Twn = values(pred_maps$winter_tmin),
												Pwn = values(pred_maps$winter_pr))
		
		# Calculate growth rate
		values(growthr_stack$rInf) <- values(growthr_stack$aSv)/2 + sqrt(((values(growthr_stack$aSv))^2)/4 + (values(growthr_stack$jSv) * values(growthr_stack$fec)) )
		
		# mask to Swiss border
		growthr_stack <- raster::mask(growthr_stack, bg)
		
		# get mean growth rate over Switzerland
		growth_TS[(year-1998),specID] <- mean(values(growthr_stack$rInf),na.rm=T)
		
		# Plot all
			#png(file = paste0("results/analysis/rates_maps/spec",specID,"/allRates_year",year,"_spec",specID,".png"), width = 960, height = 640)
			#	print(spplot(growthr_stack, main = sp_name_eng, zlim = c(0,2)))
			#dev.off()
		
			# only growth rate
			#my.settings <- list(
			#	strip.background = list(col="gray90")
			#	#,strip.border=list(col="transparent")
			#	#,layout = c(1,2)
			#)
			#col_palette_grain <- 11
			#brk <- seq(-0.4,2.4,0.2)
			#cols <- pal_diverging(length(brk))
			
			#png(file = paste0("results/analysis/rates_maps/spec",specID,"/growth0_year",year,"_spec",specID,".png"), width = 720, height = 480)
			#	print(spplot(growthr_stack$rInf, main = sp_name_eng, col.regions = cols, zlim = c(0,2.4), at = brk, colorkey = list(space = "right", width = .75), par.settings = my.settings, par.strip.text = list(cex=1, lines=1)))
			#dev.off()
		
		# Save stack, too
		save(growthr_stack, file = paste0("results/analysis/rates_maps/stacks/stack_spec",specID,"_year",year,".Rdata"))
	
	#end year loop
	}
# end loop over species
}
# growth_TS = list("50%"=growth_TS_med,"75%"=growth_TS_upp,"25%"=growth_TS_low)
# save(growth_TS, file = "results/analysis/rates_maps/spatial_mean/growthrate_TS.Rdata")
}


###--- 5.) Make posterior predictions & Calculate c-Index, plot times series

###--- 5.1.) Make posterior predictions

if(FALSE){                                                                      #-- ran once with CCcorr TRUE and once with CCcorr FALSE 
	rm(postpred_counts)
	# load CC corrected climate data if needed
	CCcorr <- FALSE
	if(CCcorr){
		load(envir = .GlobalEnv, file = "data/climatemaps/aggrVariables_zTrafod_CCcorrected.Rdata")
		# winter precipitation of year 2019 is missing. replace with that of 2018 because last years survival doesn't count into abundance anyway
		clim_CCcorrd[[21]] <- addLayer(clim_CCcorrd[[21]],clim_CCcorrd[[20]]$winter_pr)
	}
	
	# get MHB data
	load("data/species_counts/MHBabunds_1999-2021.Rdata")
	# store MHB cell numbers
	MHBcells <- species_counts[[1]]$cell
	
	# Swiss border for masking
	bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')
	
	# number of draws from the posterior
	post_samplesize <- 396
	# postpred_counts <- lapply(sel_specs, function(specID){ # BUG: lapply() function is not working correctly: reshaping seem not to work, pred_mhbcells is incorrect!
	                                                         # used for loop instead
															# calculated occupancy probability directly together with mean and var

	for(specID in sel_specs) {	
		sp_name <- species_names[specID]
		sp_MHB <- species_counts[[specID]][,1:22]
		
		cat(c(specID, sp_name, "\n"))
	
		# remove old data
		rm(post_pred_total, post_pred_mean, post_pred_var, s) # made sure, that old variables are removed; specifically 's' was masked by .GlobalEnv (2nd species in loop onwards)
		
		#- a) Prepare RS setup
		
		# load one of the calibration results files to get simulation setup that was used for calibration
		filename <- paste0("results/CA_out/CA_mcmc_Spec",specID,"_Batch_55",ifelse(specID==10,0,specID),"_it60000.RData") #-- used batch 55 as basis
		attach(filename, name = "CAsetup")
		
		# replace climate data with CC corrected one
		if(CCcorr) assign("yearly_envStacks_standard", value = clim_CCcorrd, pos = "CAsetup")
		
		# set number of RS replicates to 1
		s@simul@Replicates <- 1
		
		# reshape counts to long format
		sp_MHB <- reshape(sp_MHB, direction = "long", varying = list(2:22), timevar = "Year", idvar = "cell", v.names = "Count")
	
		# prepare object to store time series of total abundance for all posterior draws
		post_pred_total <- array(NA, dim = c(post_samplesize, s@simul@Years-s@simul@OutStartPop))
		dimnames(post_pred_total)[[2]] <- paste0("year_", (1+s@simul@OutStartPop):s@simul@Years)
		
	
		#- b) Run simulations
		
		# read the combined MCMC
		simnames <- batches[,which(specs==specID)]
		simnames <- simnames[!is.na(simnames)]
		filename <- paste0("results/CA_out/CA_comb_Spec",specID,"_it",iter,"_",simnames[1],"-",simnames[length(simnames)],".Rdata")
		attach(filename, name = "currChain")
		
		# get sample from joint posterior
		BT_smplst <- convertCoda(get("sampler_list",pos="currChain"), names = par_names)
		post_sample <- getSample(BT_smplst, parametersOnly = T, coda = F, start = 1+iter/2, numSamples = post_samplesize ) #-- changed 3 to 2 as we now hav 3*60.000 iterations
		detach("currChain", character.only = TRUE)
	
		# loop over samples
		for (smpl in 1:post_samplesize){
			if(smpl %% 10 == 1) print(paste("sample #",smpl))
			s@simul@Replicates <- 1												#-- reset value to 1
			# generate prediction -> returns raster stack
			smpl_pred <- runRSmodel(post_sample[smpl,])
			
			# crop to Switzerland to exclude population in buffer area from the total
			smpl_pred <- raster::mask(smpl_pred, bg)
			
			# get total abundance for this fold (this gets stored for each posterior sample)
			post_pred_total[smpl,] <- sapply(1:ncol(post_pred_total), FUN = function(lyr_ix){ 
				if(nlayers(smpl_pred)>=lyr_ix) sum(values(smpl_pred[[lyr_ix]]), na.rm=T )
			})
			
			# extract counts in MHB cells
			smpl_pred <- raster::extract(smpl_pred, MHBcells, df = TRUE)
			smpl_pred$ID <- MHBcells
			
			# add all predicted counts and their squares together to calculate mean & variance at the end
			if(smpl==1) {
				post_pred_mean <- post_pred_var  <- post_pred_occ <- smpl_pred  #-- added calculation of occupancy probability
				post_pred_var[,-1]  <- (smpl_pred[,-1])^2
				post_pred_occ [,-1] <- as.integer(!(!smpl_pred[,-1]))
			}else{
				post_pred_mean[,-1] <- post_pred_mean[,-1] + smpl_pred[,-1]
				post_pred_var[,-1] <- post_pred_var[,-1] + (smpl_pred[,-1])^2
				post_pred_occ[,-1] <- post_pred_occ [,-1] + as.integer(!(!smpl_pred[,-1]))
			}
		}
		
		# detach namespaces again
		detach("CAsetup", character.only = TRUE)
		
		# calculate final values of mean & variance in MHB cells for this fold
		post_pred_mean[,-1] <- post_pred_mean[,-1]/post_samplesize
		post_pred_occ[,-1] <- post_pred_occ[,-1]/post_samplesize
		post_pred_var[,-1] <- (post_pred_var[,-1]/post_samplesize) - (post_pred_mean[,-1])^2
		
		##- c) Process results
		
		# reshape data frames to long format
		post_pred_mean <- reshape(post_pred_mean, direction = "long", varying = list(2:22), timevar = "Year", idvar = "ID", v.names = "Abund")
		post_pred_var  <- reshape(post_pred_var , direction = "long", varying = list(2:22), timevar = "Year", idvar = "ID", v.names = "varAbund")
		post_pred_occ <- reshape(post_pred_occ, direction = "long", varying = list(2:22), timevar = "Year", idvar = "ID", v.names = "Occ")
		
		if(identical(post_pred_mean$ID,post_pred_var$ID)){post_pred <- cbind(post_pred_mean,post_pred_var[3],post_pred_occ[3])
		}else{warning("Format error in prediction output")}
		
		print(head(post_pred))
		
		# projections to MHB cells:
		pred_mhbcells <- merge.data.frame(sp_MHB,post_pred,by.x = c("cell","Year"), by.y = c("ID","Year"))
		
		# store Rdata for each species
		postpred_counts_specs <- list(pred_mhbcells = pred_mhbcells,            #-- generate data for each species and store it for each species (as now we have a for loop and not lapply)
			 pred_total    = post_pred_total)
		
		# names(postpred_counts_specs) <- species_names[specID]
		
		#- Save all results to file
		save(postpred_counts_specs, file = paste0("results/analysis/postpred_counts_spec_",specID,ifelse(CCcorr,"_CCcorr",""),".RData"))
		
		rm(postpred_counts_specs)
		
		# build return object
		# postpred_counts <- c(postpred_counts,list(pred_mhbcells = pred_mhbcells,
		# 	 pred_total    = post_pred_total))
		# list(pred_mhbcells = pred_mhbcells,
		# 	pred_total    = post_pred_total)
		
	# end loop over species
	}
	
	# combine lists
	postpred_counts <- lapply(sel_specs, function(specID){
		load(paste0("results/analysis/postpred_counts_spec_",specID,ifelse(CCcorr,"_CCcorr",""),".RData"))
		postpred <- postpred_counts_specs
		rm(postpred_counts_specs)
		postpred
		})
	names(postpred_counts) <- species_names[sel_specs]
		
	#- Save all results to file
	save(postpred_counts, file = paste0("results/analysis/postpred_counts",ifelse(CCcorr,"_CCcorr",""),".RData"))
}else{
	# load predictions from file
	load("results/analysis/postpred_counts.RData")
	# load("results/analysis/postpred_counts_CCcorr.RData")
}
# 

###--- 5.1.) Plot total abundance time series & Calculate c-Index              #-- also calculate mean and std of total abundance
                                                                               #-- ploted relative total abundance and diff. to BBI in one plot

if(TRUE){
	calc_quantiles <- c(0.025,0.05,0.10,0.20,0.25,0.50,0.75,0.80,0.90,0.95,0.975) #-- calc_quantiles was already used before with different order
	plot_quantiles <- paste0(calc_quantiles[c(1,6,11)]*100,"%")                    #-- larger CI only used for temporary plotting
	years <- 1999:2019
	
	cols <- cols_specID
	col_shade <- sapply(cols,t_col,percent = 50)
	col_line <- sapply(cols,makeColDarker)
	
	# Plot time series of relative total abundance for all species in one plot
	# tikz(file = "results/analysis/total_abund/plot_timeseries/totAbund_allSpec.tex", width = 3.84, height = 2.6, standAlone = FALSE) #, width = 5.76, height = 4.84
	pdf(file = "results/analysis/total_abund/plot_timeseries/totAbund_allSpec.pdf", width = 3.84, height = 2.6) #, width = 5.76, height = 4.84 #-- for now used pdf graphic device
	{
		par( cex = 0.9, mgp = c(2,0.6,0), mar=c(3.4,3.2,0.1,0.1) )
		plot(NULL, xlim = range(years)+c(0,0.3), ylim = c(0.5,2.3), xlab = "Year" , ylab = "Relative total abundance") 
		for(specID in rev(sel_specs)){
			sp_name <- species_names[specID]
			cat(c(specID, sp_name, "\n"))
			
			# get counts
			sp_pred <- postpred_counts[[which(sel_specs==specID)]]
			
			# calculate counts quantiles
			pred_total_quants <- apply(sp_pred$pred_total, 2, quantile, probs = calc_quantiles )
			
			# normalize all quantiles to median counts of year 1999 (which is the fifth simulated year)
			pred_total_quants <- pred_total_quants / pred_total_quants["50%","year_5"]
			
			# plot
			polygon(x = c(years,rev(years)), y = c(pred_total_quants["10%",],rev(pred_total_quants["90%",])),  col = col_shade[specID], border = NA) #-- 80% CI
			lines(x = years, y = pred_total_quants["50%",],  col = col_line[specID], lwd = 2, type = "b", pch=19, cex=0.5)
		}
	}
	dev.off()
	
	# Plot time series of difference of simulated relative total abundance to Breeding bird index
	
	Art_names <- read.csv("MHB_data/key2birdnames.csv") #-- adapted path
	ArtID <- c(5550, 5480, 3830, 3940, 3910, 4820, 4900, 5370, 4230, 4910)
	BBI <- read.csv("MHB_data/MHB 1999-2019 vorläufige Ergebnisse Basis 1999.csv") #-- adapted path
	BBI_sel <- subset(BBI, Species %in% ArtID)
	
	# tikz(file = "results/analysis/total_abund/plot_timeseries/dBBI_allSpec.tex", width = 3.84, height = 2.6, standAlone = FALSE) #, width = 5.76, height = 4.84
	pdf(file = "results/analysis/total_abund/plot_timeseries/dBBI_allSpec.pdf", width = 3.84, height = 2.6) #, width = 5.76, height = 4.84 #-- for now used pdf graphic device
	{
		par( cex = 0.9, mgp = c(2,0.6,0), mar=c(3.4,3.2,0.1,0.1) )
		plot(NULL, xlim = range(years)+c(0,0.3), ylim = c(-0.9,0.9), xlab = "Year" , ylab = "Difference of relative abundance to BBI") 
		for(specID in rev(sel_specs)){
			sp_name <- species_names[specID]
			cat(c(specID, sp_name, "\n"))
			
			# get counts
			sp_pred <- postpred_counts[[which(sel_specs==specID)]]
			
			# calculate counts quantiles
			pred_total_quants <- apply(sp_pred$pred_total, 2, quantile, probs = calc_quantiles )
			
			# normalize all quantiles to median counts of year 1999 (which is the fifth simulated year)
			pred_total_quants <- pred_total_quants / pred_total_quants["50%","year_5"]
			
			# now, get Breeding bird index
			sp_BBI <- subset(BBI_sel, Species == ArtID[specID], select = Model)
			
			# plot difference of simulated relative abundance and BBI
			polygon(x = c(years,rev(years)), y = c(pred_total_quants["10%",]-sp_BBI$Model,rev(pred_total_quants["90%",]-sp_BBI$Model)),  col = col_shade[specID], border = NA) #-- now 80% CI
			lines(x = years, y = pred_total_quants["50%",]-sp_BBI$Model,  col = col_line[specID], lwd = 2, type = "b", pch=19, cex=0.5 )
		}
	}
	dev.off()
	
	# Plot both time series in one plot
	tikz(file = "results/analysis/total_abund/plot_timeseries/totAbund_dBBI_allSpec_80CI.tex", width = 6.85, height = 3, standAlone = FALSE) #, width = 5.76, height = 4.84
	# pdf(file = "results/analysis/total_abund/plot_timeseries/totAbund_dBBI_allSpec_95CI.pdf", width = 6.85, height = 3) #, width = 5.76, height = 4.84
	{
		# Layout of plot
		layout.matrix <- matrix(c(1,3,2,4
		), nrow = 1, ncol = 4, byrow = T)

		
		layout(mat = layout.matrix,
			   heights = c(1,1,1,1), # Heights of the two rows
			   widths = c(9,2,9,1)) # Widths of the two columns
		
		par(cex = 0.8,  mgp = c(1.7,0.6,0), mar=c(3.5,3,0.2,0.1))
		
		##- plot relative total abundance
		plot(NULL, xlim = range(years)+c(0,1.3), ylim = c(0.5,2.5), xlab = "Year" , ylab = "Relative total abundance") 
		plot.text <- data.frame(x1=c(NA,rep(2020,4),NA,rep(2020,4)),
								y1=c(NA,0.8,1.03,1.55,1.4,NA,1.09,2.5,0.95,0.87), 
								x2=c(NA,rep(2020,4),NA,rep(2020,4)),
								y2=c(NA,-0.05,-0.15,-0.65,0.40,NA,-0.25,1.0,0.10,0.2), 
								label=c("","BF", "CT","TC","NH","","DU","LI","RO","AA"))
		for(specID in rev(sel_specs)){
			sp_name <- species_names[specID]
			cat(c(specID, sp_name, "\n"))
			
			# get counts
			sp_pred <- postpred_counts[[which(sel_specs==specID)]]
			
			# calculate counts quantiles
			pred_total_quants <- apply(sp_pred$pred_total, 2, quantile, probs = calc_quantiles )
			
			# normalize all quantiles to median counts of year 1999 (which is the fifth simulated year)
			pred_total_quants <- pred_total_quants / pred_total_quants["50%","year_5"]
			
			# plot
			polygon(x = c(years,rev(years)), y = c(pred_total_quants["10%",],rev(pred_total_quants["90%",])),  col = col_shade[specID], border = NA)
			lines(x = years, y = pred_total_quants["50%",],  col = col_line[specID], lwd = 2, type = "b", pch=19, cex=0.5)
			text(x=plot.text$x1[specID],y=plot.text$y1[specID], label=plot.text$label[specID],cex=0.8)
		}
		
		##- plot difference of simulated relative total abundance to Breeding bird index
		plot(NULL, xlim = range(years)+c(0,1.3), ylim = c(-0.9,1.2), xlab = "Year" , ylab = "Deviance of rel. abundance from BBI") 
		for(specID in rev(sel_specs)){
			sp_name <- species_names[specID]
			cat(c(specID, sp_name, "\n"))
			
			# get counts
			sp_pred <- postpred_counts[[which(sel_specs==specID)]]
			
			# calculate counts quantiles
			pred_total_quants <- apply(sp_pred$pred_total, 2, quantile, probs = calc_quantiles )
			
			# normalize all quantiles to median counts of year 1999 (which is the fifth simulated year)
			pred_total_quants <- pred_total_quants / pred_total_quants["50%","year_5"]
			
			# now, get Breeding bird index
			sp_BBI <- subset(BBI_sel, Species == ArtID[specID], select = Model)
			
			# plot difference of simulated relative abundance and BBI
			polygon(x = c(years,rev(years)), y = c(pred_total_quants["10%",]-sp_BBI$Model,rev(pred_total_quants["90%",]-sp_BBI$Model)),  col = col_shade[specID], border = NA)
			lines(x = years, y = pred_total_quants["50%",]-sp_BBI$Model,  col = col_line[specID], lwd = 2, type = "b", pch=19, cex=0.5 )
			text(x=plot.text$x2[specID],y=plot.text$y2[specID], label=plot.text$label[specID],cex=0.8)
		}
		
		##- plot visual help for interpretation
		
		# for total abundance
		par(mar=c(4,0.1,0.2,2.0), cex=0.7)
		# for total abundance
		plot(NULL, xlim = c(0,10), ylim = c(0.5,2.5), axes = F, xaxt='n', xlab="")
		arrows(x0=3, y0=1, x1 = 3, y1 = 2.3, length = 0.1, angle = 30,
			   code = 2,lwd=1.5)
		arrows(x0=3, y0=1, x1 = 3, y1 = 0.5, length = 0.1, angle = 30,
			   code = 2,lwd=1.5)
		lines(x=c(0,6), y=c(1,1),lwd=1.5)
		text(x=8,y=1.25, label="increase",srt=90)
		text(x=8,y=0.75, label="decrease",srt=90)
		
		par(mar=c(4,0.1,0.2,0.4), cex=0.7)
		# for BBI difference
		plot(NULL, xlim = c(0,10), ylim = c(-0.9,1.2), axes = F, xaxt='n', xlab="")
		lines(x=c(0,6), y=c(0,0),lwd=1.5)
		arrows(x0=3, y0=0, x1 = 3, y1 = 0.5, length = 0.1, angle = 30,
			   code = 2,lwd=1.5)
		arrows(x0=3, y0=0, x1 = 3, y1 = -0.5, length = 0.1, angle = 30,
			   code = 2,lwd=1.5)
		text(x=8,y=0.5, label="overestimated",srt=90)
		text(x=8,y=-0.5, label="underestimated",srt=90)
		
	}
	dev.off()
	
	##- additional: plot mean and standard deviation
	
	# Plot time series of relative total abundance for all species in one plot -> mean and standard deviation
	# tikz(file = "results/analysis/total_abund/plot_timeseries/totAbund_allSpec_mean.tex", width = 3.84, height = 2.6, standAlone = FALSE) #, width = 5.76, height = 4.84
	pdf(file = "results/analysis/total_abund/plot_timeseries/totAbund_allSpec_mean.pdf", width = 3.84, height = 2.6) #, width = 5.76, height = 4.84
	{
		par( cex = 0.9, mgp = c(2,0.6,0), mar=c(3.4,3.2,0.1,0.1) )
		plot(NULL, xlim = range(years)+c(0,0.3), ylim = c(0.5,2.3), xlab = "Year" , ylab = "Relative total abundance") 
		for(specID in rev(sel_specs)){
			sp_name <- species_names[specID]
			cat(c(specID, sp_name, "\n"))
			
			# get counts
			sp_pred <- postpred_counts[[which(sel_specs==specID)]]
			
			# calculate counts mean and std
			mean_std <- function(x) c(mean(x),sd(x))
			
			pred_total_mean <- apply(sp_pred$pred_total, 2, mean_std)
			
			# normalize all quantiles to median counts of year 1999 (which is the fifth simulated year)
			pred_total_mean <- pred_total_mean / pred_total_mean[1,"year_5"]
			
			# plot
			polygon(x = c(years,rev(years)), y = c(pred_total_mean[1,]-pred_total_mean[2,],rev(pred_total_mean[1,]+pred_total_mean[2,])),  col = col_shade[specID], border = NA)
			lines(x = years, y = pred_total_mean[1,],  col = col_line[specID], lwd = 2, type = "b", pch=19, cex=0.5)
		}
	}
	dev.off()
	
	# Plot time series of  relative total abundance for all species in one plot --> mean and standard deviation
	# tikz(file = "results/analysis/total_abund/plot_timeseries/totAbund_dBBI_allSpec_mean.tex", width = 3.84, height = 2.6, standAlone = FALSE) #, width = 5.76, height = 4.84
	pdf(file = "results/analysis/total_abund/plot_timeseries/dBBI_allSpec_mean.pdf", width = 3.84, height = 2.6) #, width = 5.76, height = 4.84
	{
		par( cex = 0.9, mgp = c(2,0.6,0), mar=c(3.4,3.2,0.1,0.1) )
		plot(NULL, xlim = range(years)+c(0,0.3), ylim = c(-0.9,0.9), xlab = "Year" , ylab = "Difference of relative abundance to BBI") 
		for(specID in rev(sel_specs)){
			sp_name <- species_names[specID]
			cat(c(specID, sp_name, "\n"))
			
			# get counts
			sp_pred <- postpred_counts[[which(sel_specs==specID)]]
			
			# calculate counts mean and std
			mean_std <- function(x) c(mean(x),sd(x))
			
			pred_total_mean <- apply(sp_pred$pred_total, 2, mean_std)
			
			# normalize all to mean counts of year 1999 (which is the fifth simulated year)
			pred_total_mean <- pred_total_mean / pred_total_mean[1,"year_5"]
			
			# now, get Breeding bird index
			sp_BBI <- subset(BBI_sel, Species == ArtID[specID], select = Model)
			
			# plot difference of simulated relative abundance and BBI
			polygon(x = c(years,rev(years)), y = c(pred_total_mean[1,]-pred_total_mean[2,]-sp_BBI$Model,rev(pred_total_mean[1,]+pred_total_mean[2,]-sp_BBI$Model)),  col = col_shade[specID], border = NA)
			lines(x = years, y = pred_total_mean[1,]-sp_BBI$Model,  col = col_line[specID], lwd = 2, type = "b", pch=19, cex=0.5)
		}
	}
	dev.off()
	
	# Plot both series in one plot
	
	tikz(file = "results/analysis/total_abund/plot_timeseries/totAbund_dBBI_allSpec_mean.tex", width = 6.85, height = 3, standAlone = FALSE) #, width = 5.76, height = 4.84
	# pdf(file = "results/analysis/total_abund/plot_timeseries/totAbund_dBBI_allSpec_mean.pdf", width = 6.85, height = 3.84) #, width = 5.76, height = 4.84
	{
		# Layout of plot
		layout.matrix <- matrix(c(1,3,2,4
		), nrow = 1, ncol = 4, byrow = T)
		
		
		layout(mat = layout.matrix,
			   heights = c(1,1,1,1), # Heights of the two rows
			   widths = c(9,2,9,1)) # Widths of the two columns
		
		par(cex = 0.8,  mgp = c(1.7,0.6,0), mar=c(3.5,3,0.2,0.1))
		
		##- plot relative total abundance
		plot(NULL, xlim = range(years)+c(0,1.3), ylim = c(0.5,2.5), xlab = "Year" , ylab = "Relative total abundance") 
		plot.text <- data.frame(x1=c(NA,rep(2020,4),NA,rep(2020,4)),
								y1=c(NA,0.8,1.03,1.55,1.4,NA,1.09,2.5,0.95,0.87), 
								x2=c(NA,rep(2020,4),NA,rep(2020,4)),
								y2=c(NA,-0.05,-0.15,-0.65,0.40,NA,-0.25,1.0,0.10,0.2), 
								label=c("","BF", "CT","TC","NH","","DU","LI","RO","AA"))
		for(specID in rev(sel_specs)){
			sp_name <- species_names[specID]
			cat(c(specID, sp_name, "\n"))
			
			# get counts
			sp_pred <- postpred_counts[[which(sel_specs==specID)]]
			
			# calculate counts mean and std
			mean_std <- function(x) c(mean(x),sd(x))
			
			pred_total_mean <- apply(sp_pred$pred_total, 2, mean_std)
			
			# normalize all quantiles to median counts of year 1999 (which is the fifth simulated year)
			pred_total_mean <- pred_total_mean / pred_total_mean[1,"year_5"]
			
			# plot
			polygon(x = c(years,rev(years)), y = c(pred_total_mean[1,]-pred_total_mean[2,],rev(pred_total_mean[1,]+pred_total_mean[2,])),  col = col_shade[specID], border = NA)
			lines(x = years, y = pred_total_mean[1,],  col = col_line[specID], lwd = 2, type = "b", pch=19, cex=0.5)
			text(x=plot.text$x1[specID],y=plot.text$y1[specID], label=plot.text$label[specID],cex=0.8)
		}
		
		##- plot difference of simulated relative total abundance to Breeding bird index
		plot(NULL, xlim = range(years)+c(0,1.3), ylim = c(-0.9,1.2), xlab = "Year" , ylab = "Deviance of rel. abundance from BBI") 
		for(specID in rev(sel_specs)){
			sp_name <- species_names[specID]
			cat(c(specID, sp_name, "\n"))
			
			# get counts
			sp_pred <- postpred_counts[[which(sel_specs==specID)]]
			
			# calculate counts mean and std
			mean_std <- function(x) c(mean(x),sd(x))
			
			pred_total_mean <- apply(sp_pred$pred_total, 2, mean_std)
			
			# normalize all to mean counts of year 1999 (which is the fifth simulated year)
			pred_total_mean <- pred_total_mean / pred_total_mean[1,"year_5"]
			
			# now, get Breeding bird index
			sp_BBI <- subset(BBI_sel, Species == ArtID[specID], select = Model)
			
			# plot difference of simulated relative abundance and BBI
			polygon(x = c(years,rev(years)), y = c(pred_total_mean[1,]-pred_total_mean[2,]-sp_BBI$Model,rev(pred_total_mean[1,]+pred_total_mean[2,]-sp_BBI$Model)),  col = col_shade[specID], border = NA)
			lines(x = years, y = pred_total_mean[1,]-sp_BBI$Model,  col = col_line[specID], lwd = 2, type = "b", pch=19, cex=0.5)
			text(x=plot.text$x2[specID],y=plot.text$y2[specID], label=plot.text$label[specID],cex=0.8)
		}
		
		##- plot visual help for interpretation
		
		# for total abundance
		par(mar=c(4,0.1,0.2,2.0), cex=0.7)
		# for total abundance
		plot(NULL, xlim = c(0,10), ylim = c(0.5,2.3), axes = F, xaxt='n', xlab="")
		arrows(x0=3, y0=1, x1 = 3, y1 = 2.3, length = 0.1, angle = 30,
			   code = 2,lwd=1.5)
		arrows(x0=3, y0=1, x1 = 3, y1 = 0.5, length = 0.1, angle = 30,
			   code = 2,lwd=1.5)
		lines(x=c(0,6), y=c(1,1),lwd=1.5)
		text(x=8,y=1.25, label="increase",srt=90)
		text(x=8,y=0.75, label="decrease",srt=90)
		
		par(mar=c(4,0.1,0.2,0.4), cex=0.7)
		# for BBI difference
		plot(NULL, xlim = c(0,10), ylim = c(-0.9,0.9), axes = F, xaxt='n', xlab="")
		lines(x=c(0,6), y=c(0,0),lwd=1.5)
		arrows(x0=3, y0=0, x1 = 3, y1 = 0.5, length = 0.1, angle = 30,
			   code = 2,lwd=1.5)
		arrows(x0=3, y0=0, x1 = 3, y1 = -0.5, length = 0.1, angle = 30,
			   code = 2,lwd=1.5)
		text(x=8,y=0.5, label="overestimated",srt=90)
		text(x=8,y=-0.5, label="underestimated",srt=90)
		
	}
	dev.off()
	
	
	# Now, for each single species
	
	col_shade <- c(t_col(cols_qualitative[10],percent = 50),t_col(cols_qualitative[3],percent = 50))
	col_line <- c(makeColDarker(cols_qualitative[10]),makeColDarker(cols_qualitative[3]))
	
	## -- calculate spatial and temporal c-index (based on abund. in MHB) 
	# Loop over species 
	postpred_cindex <- sapply(sel_specs, function(specID){
		
		sp_name <- species_names[specID]
		cat(c(specID, sp_name, "\n"))
		
		sp_pred <- postpred_counts[[which(sel_specs==specID)]]
		
		
		##-- a) Total abundance of single species
		
		# calculate counts quantiles
		pred_total_quants <- apply(sp_pred$pred_total, 2, quantile, probs = calc_quantiles )
		
		# plot time series of total abundance
		#png(file = paste0("results/analysis/total_abund/plot_timeseries/totAbund_",specID,".png"), width = 720, height = 480)
			plot(NULL, xlim = range(years), ylim = range(pred_total_quants[plot_quantiles,]), main = sp_name, xlab = "Year" , ylab = "Total abundance") 
			polygon(x = c(years,rev(years)), y = c(pred_total_quants["25%",],rev(pred_total_quants["75%",])),  col = col_shade[2], border = NA)
			lines(x = years, y = pred_total_quants["50%",],  col = col_line[2], lwd = 2, type = "b" )
		#dev.off()
		
		
		##-- b) spatiotemporal c-index
		cix_mhbcells <- with(sp_pred$pred_mhbcells, rcorr.cens(Abund, Count))   
		
		
		##-- c) spatial AUC and c-index
		# # get median of dispersion parameter
		# NB_size <- quants[[which(names(quants)==sp_name)]]["50%","GPsize"]
		# # generate occupancy probability
		# sp_pred$pred_occprob <- dnbinom( 0, mu = sp_pred$pred_mhbcells$Abund, size = NB_size)
		# flatten MHB to P/A
		# sp_pred$pred_mhbcells$mhb_PA <- as.integer(as.logical(sp_pred$pred_mhbcells$Count))
		# calculate AUC
		# occ_roc <- roc(sp_pred$mhb_PA, sp_pred$pred_occprob, na.rm = TRUE, smooth = FALSE, ci = TRUE, plot = FALSE, ret = c("roc", "coords"))
		# occ_roc <- roc(sp_pred$pred_mhbcells$mhb_PA, sp_pred$pred_mhbcells$Occ, na.rm = TRUE, smooth = FALSE, ci = TRUE, plot = FALSE, ret = c("roc", "coords"))
		
		# ##-- d) R^2 // coefficient of determination
		# # variance of reference data
		# tss <- sum( (sp_pred$pred_mhbcells$Count - mean(sp_pred$pred_mhbcells$Count, na.rm=T))^2, na.rm=T )
		# # variance of residuals
		# rss <- sum( (sp_pred$pred_mhbcells$Abund - sp_pred$pred_mhbcells$Count)^2, na.rm=T )
		# # R^2
		# Rsqr <- 1-rss/tss
		
		
		##-- e) RSME // root mean squared error
		# rsme <- sqrt( mean( (sp_pred$pred_mhbcells$Abund - sp_pred$pred_mhbcells$Count)^2, na.rm=T ) ) #-- TODO output of rsme
		
	
		# Return
		cix_mhbcells
		 #occ_roc
		 #Rsqr
		 #rsme
	})
	colnames(postpred_cindex) <- species_names[sel_specs]
	
	# Plot spatiotemporal c-index
	postpred_cindex <- data.frame(t(postpred_cindex[c("C Index","S.D."),] ))
	postpred_cindex <- cbind(Species = rownames(postpred_cindex),postpred_cindex)
	names(postpred_cindex) <- c("Species","C_Index","SD")
	
	cix_ggplot <- 
		ggplot(postpred_cindex,
			   aes(
			   	x = C_Index,
			   	y = Species
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = C_Index+2*SD, xmin = C_Index-2*SD, height = .2)) +
		xlim(NA, 1) +
		theme_bw() +
		theme(panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank(),
			  panel.background = element_rect(colour = "black"),
			  axis.text = element_text(colour = "black"),
			  panel.grid = element_line(colour = "gray") ) +
		#coord_flip() + 
		scale_y_discrete(name ="", limits=rev(rownames(postpred_cindex))) +
		labs(x = "C-index", y="")

	# export as tikz
	tikz(file = "results/analysis/cIndex/plot_cix_spatialtemporal.tex", width = 3.2, height = 1.9, standAlone = FALSE)
	plot(cix_ggplot)
	dev.off()
	
	## -- calculate spatial (mean yearly AUC based on abundance) c-index and spatial (mean yearly AUC based on pres/abs and occ. prob) AUC (each based on MHB)
	# Loop over species 
	if(exists("postpred_cindex_spatial")) rm(postpred_cindex_spatial)
	for(specID in sel_specs){
		sp_name <- species_names[specID]
		cat(c(specID, sp_name, "\n"))
		
		sp_pred <- postpred_counts[[which(sel_specs==specID)]]
		
		##-- a) Total abundance of single species
		
		# calculate counts quantiles
		pred_total_quants <- apply(sp_pred$pred_total, 2, quantile, probs = calc_quantiles)
		
		##-- b) AUC spatial
		
		# flatten MHB to P/A
		sp_pred$pred_mhbcells$mhb_PA <- as.integer(as.logical(sp_pred$pred_mhbcells$Count))
		
		# calculate spatial AUC and c-index per year
		spatial_auc_cindex <- sapply(1:21, function(yr){
			pred_mhbcells <- sp_pred$pred_mhbcells[sp_pred$pred_mhbcells$Year==yr,]
			
			# compare abundances
			cix_mhbcells <- with(pred_mhbcells, rcorr.cens(Abund, Count))

			# compare presence/absence
			auc_mhbcells <- with(pred_mhbcells, rcorr.cens(Occ, mhb_PA))
			# auc_mhbcells <- roc(pred_mhbcells$mhb_PA, pred_mhbcells$Occ, na.rm = TRUE, smooth = FALSE, ci = TRUE, plot = TRUE, auc=TRUE, ret = c("roc", "coords"))
			c(cindex = cix_mhbcells,
				  auc = auc_mhbcells)
		})
		
		##-- c) calculate mean c-index (based on abundance) and auc (based on presence/absence)
		mean_auc_cindex  <- apply(spatial_auc_cindex, 1, mean)
		if(!exists("postpred_cindex_spatial")) postpred_cindex_spatial <- mean_auc_cindex
			else postpred_cindex_spatial <- cbind(postpred_cindex_spatial,as.vector(mean_auc_cindex))

	}	
	colnames(postpred_cindex_spatial) <- species_names[sel_specs]
	
	# Plot c-index and auc
	postpred_cindex_spatial <- data.frame(t(postpred_cindex_spatial[c("cindex.C Index","cindex.S.D.","auc.C Index","auc.S.D."),] ))
	postpred_cindex_spatial <- cbind(Species = rownames(postpred_cindex_spatial),postpred_cindex_spatial)
	names(postpred_cindex_spatial) <- c("Species","c.index.C_Index","c.index.SD","auc.C_Index","auc.SD")
	
	cix_ggplot <- 
		ggplot(postpred_cindex_spatial,
			   aes(
			   	x = c.index.C_Index,
			   	y = Species
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = c.index.C_Index+2*c.index.SD, xmin = c.index.C_Index-2*c.index.SD, height = .2)) +
		xlim(NA, 1) +
		theme_bw() +
		theme(panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank(),
			  panel.background = element_rect(colour = "black"),
			  axis.text = element_text(colour = "black"),
			  panel.grid = element_line(colour = "gray") ) +
		#coord_flip() + 
		scale_y_discrete(name ="", limits=rev(rownames(postpred_cindex))) +
		labs(x = "C-index", y="")
	auc_ggplot <- 
		ggplot(postpred_cindex_spatial,
			   aes(
			   	x = auc.C_Index,
			   	y = Species
			   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = auc.C_Index+2*auc.SD, xmin = auc.C_Index-2*auc.SD, height = .2)) +
		xlim(NA, 1) +
		theme_bw() +
		theme(panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank(),
			  panel.background = element_rect(colour = "black"),
			  axis.text = element_text(colour = "black"),
			  panel.grid = element_line(colour = "gray") ) +
		#coord_flip() + 
		scale_y_discrete(name ="", limits=rev(rownames(postpred_cindex))) +
		labs(x = "AUC", y="")
	
	# export as tikz
	tikz(file = "results/analysis/cIndex/plot_cix_spatial.tex", width = 3.2, height = 1.9, standAlone = FALSE)
		plot(cix_ggplot)
	dev.off()
	tikz(file = "results/analysis/cIndex/plot_auc_spatial.tex", width = 3.2, height = 1.9, standAlone = FALSE)
	plot(auc_ggplot)
	dev.off()
	
	
	## calculate temporal AUC (based BBI data)
	
	Art_names <- read.csv("MHB_data/key2birdnames.csv")
	ArtID <- c(5550, 5480, 3830, 3940, 3910, 4820, 4900, 5370, 4230, 4910)
	BBI <- read.csv("MHB_data/MHB 1999-2019 vorläufige Ergebnisse Basis 1999.csv")
	BBI_sel <- subset(BBI, Species %in% ArtID)
	
	if(exists("postpred_cindex_BBI")) rm(postpred_cindex_BBI)
	for(specID in rev(sel_specs)){
		sp_name <- species_names[specID]
		cat(c(specID, sp_name, "\n"))
		
		# get counts
		sp_pred <- postpred_counts[[which(sel_specs==specID)]]
		
		# calculate counts quantiles
		pred_total_quants <- apply(sp_pred$pred_total, 2, quantile, probs = calc_quantiles )
		
		# normalize all quantiles to median counts of year 1999 (which is the fifth simulated year)
		pred_total_quants <- pred_total_quants / pred_total_quants["50%","year_5"]
		
		# now, get Breeding bird index
		sp_BBI <- subset(BBI_sel, Species == ArtID[specID], select = c(Time,Model))
		
		sp_BBI_pred <- cbind(sp_BBI, "pred"=pred_total_quants["50%",])
		# calculate c-index
		(cix_BBI <- rcorr.cens(sp_BBI_pred$pred, sp_BBI_pred$Model))
		if(!exists("postpred_cindex_BBI")) postpred_cindex_BBI <- cix_BBI
		else postpred_cindex_BBI <- cbind(postpred_cindex_BBI,as.vector(cix_BBI))
	}
	colnames(postpred_cindex_BBI) <- species_names[sel_specs]

	postpred_cindex_BBI<- data.frame(t(postpred_cindex_BBI[c("C Index","S.D."),] ))
	postpred_cindex_BBI <- cbind(Species = rownames(postpred_cindex_spatial),postpred_cindex_BBI)
	names(postpred_cindex_BBI) <- c("Species","C_Index","SD")
	
	# combine all indices
	postpred_cindex_spatial_auc <- postpred_cindex_spatial[,-c(2:3)]
	names(postpred_cindex_spatial_auc) <- names(postpred_cindex)[-4]
	postpred_cindex_spatial_auc$index <- "spatial_auc"
	postpred_cindex_spatial_ci <- postpred_cindex_spatial[,-c(4:5)]
	names(postpred_cindex_spatial_ci) <- names(postpred_cindex)[-4]
	postpred_cindex_spatial_ci$index <- "spatial_cindex"
	postpred_cindex$index <- "spatialtemporal_cindex_MHB"
	postpred_cindex_BBI$index <- "temporal_cindex_BBI"
	cindex <- rbind(postpred_cindex,postpred_cindex_spatial_auc,postpred_cindex_spatial_ci, postpred_cindex_BBI)
	write.table(cindex,file="results/analysis/tables/cindices.txt", sep="\t")
	save(cindex,file="results/analysis/tables/cindices.Rdata")
	# plot all indices
	
	ggplot(cindex,
		   aes(
		   	x = C_Index,
		   	y = index
		   )) +
		geom_point() +
		geom_errorbarh(aes(xmax = C_Index+2* SD, xmin =  C_Index-2* SD, height = .2)) +
		facet_grid(.~Species) +
		# xlim(NA, 1) +
		theme_bw() +
		theme(panel.grid.major.x = element_blank(),
			  panel.grid.minor.x = element_blank(),
			  panel.background = element_rect(colour = "black"),
			  axis.text = element_text(colour = "black"),
			  panel.grid = element_line(colour = "gray") ) +
		# coord_flip() +
		labs(x = "Index", y="")
}


###--- 6.) Climate attribution

###--- 6.a) Intrinsic growth rate sensitivity to climate change
if(FALSE){
	# load(paste0("results/analysis/sp_DERs_CCcorr.Rdata")) # Annes values
	load(paste0("results/analysis/sp_DERs_largerCI_CCcorr.Rdata")) # Jettes values

	growthrate <- function(fec,jSv,aSv){ aSv/2 + sqrt( (aSv^2)/4 + (jSv * fec)) }
	
	# Loop over species
	growthrate_chg <- lapply(sel_specs, function(specID){
		
		sp_name <- species_names[specID]
		cat(c(specID, sp_name, "\n"))
		
		# get demographic rates values at median predictor and CC corrected
		sp_DER_chg <- getElement(sp_DERs,sp_name)
		fec0 = sp_DER_chg$DER_Fc_Tbr$`pred2_50%`["50%",1]
		jSv0 = sp_DER_chg$DER_jS_Tat$`pred2_50%`["50%",1]
		aSv0 = sp_DER_chg$DER_aS_Twn$`pred2_50%`["50%",1]
		
		# calculate base growth rate
		grate0 <- growthrate(fec = fec0, jSv = jSv0, aSv = aSv0)
		names(grate0) <- NULL
		
		# extract growth rate changes
		grate_CC  <- c(Tbr = growthrate(fec = sp_DER_chg$DER_Fc_Tbr$`pred2_50%`["50%",2], jSv = jSv0, aSv = aSv0),
					   Pbr = growthrate(fec = sp_DER_chg$DER_Fc_Pbr$`pred2_50%`["50%",2], jSv = jSv0, aSv = aSv0),
					   Tat = growthrate(fec = fec0, jSv = sp_DER_chg$DER_jS_Tat$`pred2_50%`["50%",2], aSv = aSv0),
					   Twn = growthrate(fec = fec0, jSv = jSv0, aSv = sp_DER_chg$DER_aS_Twn$`pred2_50%`["50%",2]),
					   Pwn = growthrate(fec = fec0, jSv = sp_DER_chg$DER_jS_Pwn$`pred2_50%`["50%",2], aSv = sp_DER_chg$DER_aS_Pwn$`pred2_50%`["50%",2]))
		names(grate_CC) <- c("Tbr","Pbr","Tat","Twn","Pwn")
		
		# return
		list(medGR  = grate0,
			 CCcorr = grate_CC,
			 absChg = grate_CC-grate0,
			 relChg = grate_CC/grate0)
	})
	names(growthrate_chg) <- species_names[sel_specs]
	
	# Make latex table
	growthrate_relChg <- data.frame(t(sapply(growthrate_chg, function(sp){sp$relChg})))
	
	colnames(growthrate_relChg) <- c("$T_{br}$","$P_{br}$","$T_{at}$","$T_{wn}$","$P_{wn}$")
	rownames(growthrate_relChg)[rownames(growthrate_relChg)=="Eurasian woodcreeper"] <- "E. treecreeper"
	rownames(growthrate_relChg)[rownames(growthrate_relChg)=="Eurasian nuthatch"] <- "E. nuthatch"
	
	# set options
	options(xtable.floating = FALSE)
	options(xtable.booktabs = TRUE)
	options(xtable.sanitize.text.function = function(x){x})
	
	# print table
	print(xtable(growthrate_relChg, 
				 type = "latex",
				 align = "lccccc",
				 floating = FALSE),
		  # file = paste0("results/analysis/tables/growthrate_relChg.tex"))
		file = paste0("results/analysis/tables/growthrate_relChg.txt"))
	
	# extract medium growthrate
	(growthrate_med <- t(data.frame(t(sapply(growthrate_chg, function(sp){sp$medGR})))))
	
	colnames(growthrate_med) <- "r_${br}$"
	
	# combinge data
	growthrate_combined <- cbind(growthrate_med, growthrate_relChg)
	rownames(growthrate_combined)<-rownames(growthrate_relChg)
	rownames(growthrate_combined)[rownames(growthrate_combined)=="Eurasian woodcreeper"] <- "E. treecreeper"
	rownames(growthrate_combined)[rownames(growthrate_combined)=="Eurasian nuthatch"] <- "E. nuthatch"
	
	# set options
	options(xtable.floating = FALSE)
	options(xtable.booktabs = TRUE)
	options(xtable.sanitize.text.function = function(x){x})
	
	# print table
	print(xtable(growthrate_combined, 
				 type = "latex",
				 align = "lcccccc",
				 floating = FALSE),
		  # file = paste0("results/analysis/tables/growthrate_relChg.tex"))
		  file = paste0("results/analysis/tables/growthrate_combined.txt"))
}



###--- 6.b) Total abundance
# -> compare mean total abundance within last three years (2017-19)

if(FALSE){
	# Read in yearly stacks of CC-corrected climate variables
	attach("results/analysis/postpred_counts_CCcorr.RData", name = "CCcorr")
	postpred_counts_CCcorr <- get("postpred_counts", pos="CCcorr")
	detach("CCcorr", character.only = TRUE)
	
	# Loop over species ##- TODO: not sure if calculation is correct? (when to summarize quantiles and calculate relative vaues?)
	attr_abund <- sapply(sel_specs, function(specID){
		
		sp_name <- species_names[specID]
		cat(c(specID, sp_name, "\n"))
		
		# get abundances for historic climate & calculate mean abundance over last three years
		sp_pred <- postpred_counts[[which(sel_specs==specID)]]$pred_total
		sp_pred <- apply(sp_pred[,19:21], 2, quantile, probs = calc_quantiles )
		sp_pred <- rowMeans(sp_pred)
		
		# do the same for CC-corrected climate
		sp_pred_CCcorr <- postpred_counts_CCcorr[[which(sel_specs==specID)]]$pred_total
		sp_pred_CCcorr <- apply(sp_pred_CCcorr[,19:21], 2, quantile, probs = calc_quantiles )
		sp_pred_CCcorr <- rowMeans(sp_pred_CCcorr)
		
		# return ratio of both
		sp_pred/sp_pred_CCcorr
	})
	colnames(attr_abund) <- species_names[sel_specs]
	
	round(attr_abund[c("2.5%","50%","97.5%"),],2)

}
