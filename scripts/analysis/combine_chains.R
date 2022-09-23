
###-----------------------------------------------------------------
##   SUMMARISE AND COMBINE INDEPENDENT MCMCS
###-----------------------------------------------------------------

# Function to read the calibration outputs for repeat runs and create combined MCMC sampler list

combineMCMCs <- function(sampler,species,batchnrs,chain_iter){
	require(BayesianTools)
	
	# construct file names from info given via parameters
	filename  <- paste0("results/CA_out/",sampler,"_Spec",species,"_Batch_",(batchnrs-100),"_it",chain_iter,".RData")
	filename2 <- paste0("results/CA_out/",sampler,"_Spec",species,"_Batch_",batchnrs      ,"_it",chain_iter,".RData")
	
	# select parameters
	parSel <- (1:18) #[-(5:6)]
	
	# get object 'MCMCout' (that contains the separate chain) from each file and store in a list
	MCMC_chains <- lapply(1:length(batchnrs), function(sim) {
		# Chain 1
		name <- paste0("run_",batchnrs[sim]-100)
		attach(filename[sim], name = name)
		mcmcChain <- get("MCMCout", name)
		detach(name, character.only = TRUE)
		mcmcChain1 <- getSample(mcmcChain, coda = T, whichParameters = parSel)
		rm(mcmcChain)
		# Chain 2
		name <- paste0("run_",batchnrs[sim])
		attach(filename2[sim], name = name)
		mcmcChain <- get("MCMCout", name)
		detach(name, character.only = TRUE)
		mcmcChain2 <- getSample(mcmcChain, coda = T, whichParameters = parSel)
		rm(mcmcChain)
		# concatenate chains
		mcmcChainC <- lapply(1:length(mcmcChain1), function(c){mcmc(rbind(mcmcChain1[[c]],mcmcChain2[[c]]))})
		return(mcmcChainC)
	})
	
	# combine individual chains and return combined MCMC
	return(mcmc.list(unlist(MCMC_chains, recursive = FALSE)))
	#return(createMcmcSamplerList(MCMC_chains))
}
