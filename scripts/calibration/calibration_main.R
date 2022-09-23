#!/usr/bin/env Rscript

###-----------------------------------------------------------------
##   RangeShiftR: calibration main script
###-----------------------------------------------------------------

# This is the top-level script to run the Bayesian calibration with RangeShiftR.
# 
# It starts by including the required libraries, functions and data, sets various 
# parameters that control the type of simulation, runs a sensitivity analysis, and
# the calibration with optimizer or MCMC.


## libraries ------------------------------------------------------------------

library(raster)
library(RangeShiftR)

## directories ----------------------------------------------------------------

RSdir <- "model/"

## helper functions ------------------------------------

# these helper functions are mostly used to organise the code sections,
# as they rely on parameters from the calling environment

source("scripts/calibration/readandpreparedata.R")
source("scripts/calibration/createRSsim.R")
source("scripts/calibration/createInitInds.R")
source("scripts/calibration/runRSmodel.R")
source("scripts/calibration/CAtarget.R")


# get command line arguments
shell_args = commandArgs(trailingOnly=TRUE)
if(length(shell_args)==3) {
	species_nr <- as.integer(shell_args[1])
	analysis_nr <- as.integer(shell_args[2])
	batch_prefix <- as.integer(shell_args[3])
	arg_iterations <- as.integer(shell_args[4])
}else{
	print("No cl args given!\n")
	species_nr <- 8
	analysis_nr <- 4
	batch_prefix <- 11
	arg_iterations <- 6e4
}	# default params




###-- 1.) Set simulation & calibration parameters -----------------------------------------------------------
{
	#--- calibration parameters:
	set_calib <- list(likdist = "GamPois", sampler = "DEzs", iter = arg_iterations, parallel = 3, InfVal = -1e6) 
	# set switches to choose which analyses to run
	switch_analyses <- c("SA_oat" = FALSE, "SA_morris" = FALSE, "SA_bt" = FALSE,
						 "CA_mcmc" = FALSE, "CA_opt" = FALSE)
	switch_analyses[analysis_nr] <- TRUE
	
	#---- data type:
	simulatedData <- FALSE
	aggregate_space <- c(NA,20,25,30,25,15,25,30,25,30)[species_nr] # spatially aggregate observations?
	aggregate_time <- 0
	
	#--- simulation parameters:
	batchnr <- 10*batch_prefix+ifelse(species_nr<10,species_nr,0)
	nr_replicates <- 8
	nrYrMHB <- 23 # nr of years in MHB data set (1999-2021)
	nrYrPred <- 21  # nr of years with climate layers (1999-2019)
	nrYrBurnin <- 4 # (1995-1998) start in Atlas period (1993-1996), as this was used for initialisation
	years_ca <- 1999:2019
	nrYears <- length(years_ca)+nrYrBurnin # years in output stack: 0...(Years-1), in Pop output: 0...Years
	OutStartPop <- nrYrBurnin
	OutIntPop <- 1
	beta_limit <- 5
	beta_sd <- 1
	maxDemogRates <- c(Fecund = 8, juvSurv = 1, adSurv = 1)
	sigma <- 0.01 # 1/(15^2)  # value to replace zero counts with
}



###-- 2.) Read and prepare data -------------------------------------------------------------

{
	# get species name
	species_names <- c("Fringilla.coelebs", 
					   "Pyrrhula.pyrrhula",
					   "Lophophanes.cristatus",
					   "Certhia.familiaris",
					   "Sitta.europaea",
					   "Regulus.regulus",
					   "Prunella.modularis",
					   "Linaria.cannabina",
					   "Turdus.torquatus",
					   "Prunella.collaris")
	
	species_name <- species_names[species_nr]
	print(species_name)

	# get species and environmental data
	readandpreparedata()
	
	# set parameter ranges, reference params & SDs for priors
	bl  <- beta_limit
	CA_params  <- data.frame(name = c( "DensDep",     "Fecund"   , "juvSurv" ,  "adSurv" , "EmigProb",   "DispDist"   ,"GPsize","Fcb.pr1","Fcb.pr2","Fcb.tm1","Fcb.tm2","jSb.pr1","jSb.pr2","jSb.tm1","jSb.tm2","aSb.pr1","aSb.pr2","aSb.tm"),
							 min  = c(   1.0e-4 ,       0.01     ,    0.01   ,     0.01  ,    0.01   ,       1        ,   1.00 ,    -bl  ,    -bl  ,    -bl  ,    -bl  ,    -bl  ,    -bl  ,    -bl  ,    -bl  ,    -bl  ,    -bl  ,    -bl ),
							 def  = c(   3.0e-2 ,trt[,'fec_mean'],trt[,'jSv'],trt[,'aSv'],trt[,'emg'],trt[,'dpd']/1000,  50.00 ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0  ,    0.0 ),
							 sd   = c(   2.0e-2 , trt[,'fec_sd'] ,    0.20   ,     0.20  ,    0.20   ,       0.5      , 150.00 , beta_sd , beta_sd , beta_sd , beta_sd , beta_sd , beta_sd , beta_sd , beta_sd , beta_sd , beta_sd , beta_sd ),
							 max  = c(   1.0e-0 ,      10.00     ,    0.99   ,     0.99  ,    0.99   ,      20        , 500.00 ,     bl  ,     bl  ,     bl  ,     bl  ,     bl  ,     bl  ,     bl  ,     bl  ,     bl  ,     bl  ,     bl ))
	
	# select parameters for calibration
	par_sel <- c(1:nrow(CA_params))
	# index of variance parameter
	GPsize_parix <- which(par_sel==which(CA_params$name == "GPsize"))
}



###-- 3.) Create RS simulation -----------------------------------------------------

# reference parameters (to be used when not chosen for calibration)
ref_par <- CA_params[["def"]]
names(ref_par) <- CA_params$name

# create RSsim
s <- createRSsim()



###-- 4.) Test runs -----------------------------------------------------

#params <- ref_par[par_sel]
#out <- runRSmodel(params)

#CAtarget(CA_params$def[par_sel])



###-- 5.) Sensitivity analysis

calib_name_core <- paste0(names(which(switch_analyses)),"_Spec",species_nr,"_Batch_",batchnr)


##----- a) Local sensitivity analysis

if(switch_analyses["SA_oat"]){
	MHB_counts <- createArtificialMHBdata(ref_par)
	#par_sel <- c(7:9)
	#GPsize_parix <- which(par_sel==which(CA_params$name == "GPsize"))
	oneFactorSA <- sensitivityOAT(ll_points = 9, ll_replics = 9, plot = TRUE, minVal = set_calib$InfVal)
	pdf(file = paste0("results/SA_out/SAplot_",calib_name_core,".pdf"), width = 10, height = 7)
	plotSA_OAT(oneFactorSA)
	dev.off()
}


##----- b) Global sensitivity analysis with Morris screening

if(switch_analyses["SA_morris"]){
	require(sensitivity)
	require(BayesianTools)

	MHB_counts <- createArtificialMHBdata(ref_par)

	CAtarget_Pll <- generateParallelExecuter(CAtarget, parallel = ifelse(set_calib$parallel>1,set_calib$parallel,FALSE) ) # parallelize with BT

	morrisResult <- morris(model = CAtarget_Pll$parallelFun,
						   factors = CA_params$name[par_sel],
						   r = 15,
						   design = list(type = "oat", levels = 5, grid.jump = 3),
						   binf = CA_params$min[par_sel],
						   bsup = CA_params$max[par_sel],
						   scale = T) # scales calculated sensitivity in units of given interval, so that sensitivity is given w.r.t. its uncertainty
	#plot(morrisResult)
}

if(!any(switch_analyses[c("CA_mcmc","CA_opt")])) save.image(paste0("results/SA_out/",calib_name_core,".RData"))




###-- 6.) Calibration

##----- a) Calibration with BayesianTools MCMC

if(switch_analyses["CA_mcmc"]){
	require(BayesianTools)

	prior <- createTruncatedNormalPrior(mean = CA_params$def[par_sel],
										sd   = CA_params$sd[par_sel],
										lower = CA_params$min[par_sel],
										upper = CA_params$max[par_sel])
	#prior <- createUniformPrior(best = CA_params$def[par_sel],
	#                            lower = CA_params$min[par_sel],
	#                            upper = CA_params$max[par_sel])
	#CA_params$def[par_sel] <- prior$sampler(1)

	setup <- createBayesianSetup(CAtarget,
								 prior = prior,
								 names = CA_params$name[par_sel],
								 parallel = ifelse(set_calib$parallel>1,set_calib$parallel,FALSE),
								 plotBest = CA_params$def[par_sel],
								 plotLower = CA_params$min[par_sel],
								 plotUpper = CA_params$max[par_sel])

	if(switch_analyses["SA_bt"]){
		pdf(file = paste0("results/SA_out/SA_bt_",calib_name_core,"_",batchnr,".pdf"), width = 12, height = 10)
		sen <- plotSensitivity(setup)
		dev.off()
	}


	if(newMCMC <- TRUE){
		settings = list(iterations = set_calib$iter,
						#f = ifelse(set_calib$sampler=="DEzs",1.7,2.38/sqrt(2*length(par_sel))),
						#eps.add = 0.01,
						startValue = ifelse(set_calib$parallel>1, set_calib$parallel, 3),
						parallel = (set_calib$parallel>1)
		)

	}else{
		# restart DE sampler
		print("restart MCMC:\n")
		# load previous MCMC
		attach(paste0("results/CA_out/",paste0(names(which(switch_analyses)),"_Spec",species_nr,"_Batch_",batchnr-100),"_it",set_calib$iter,".RData"), name = "prevChain")
		previousMCMC <- MCMCout
		detach("prevChain", character.only = TRUE)
		nrChains <- set_calib$parallel
		samples_x <- getSample(previousMCMC, start = 1, thin = 1, coda = FALSE)
		nrPar <- ncol(samples_x)
		settings <- list(iterations = set_calib$iter,
						 Z = samples_x,
						 startValue =  samples_x[(nrow(samples_x)-(nrChains-1)):nrow(samples_x), ],
						 parallel = (set_calib$parallel>1))
	}
	
	MCMCout <- runMCMC(bayesianSetup = setup, sampler = set_calib$sampler, settings = settings)

	summary(MCMCout)
	save.image(paste0("results/CA_out/",calib_name_core,"_it",set_calib$iter,".RData"))

	pdf(file = paste0("results/CA_out/",calib_name_core,"_it",set_calib$iter,".pdf"), width = 12, height = 10)
	tracePlot(MCMCout, start=1)
	marginalPlot(MCMCout, start=1)
	dev.off()

}


##----- b) Calibration with Optimizer DEoptim

if(switch_analyses["CA_opt"]){
	require(DEoptim)
	##-- not here yet -- #
}
