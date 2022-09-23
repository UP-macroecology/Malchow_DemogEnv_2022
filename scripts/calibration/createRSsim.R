
###-----------------------------------------------------------------
##   CREATE RANGESHIFTER SIMULATION
###-----------------------------------------------------------------


# The function createRSsim() initialises the RangeShifter simulation.

# The model structure for the IBM is specified here. All calibration parameters 
# are set to their reference values. The function uses variables from
# the parent (global) environment to set the parameters.

createRSsim <- function() {
	
	##--- simulation params
	
	simul <- Simulation(Years = nrYears,
						Replicates = nr_replicates,
						Simulation = 1,
						OutStartPop = OutStartPop, 
						OutIntPop = OutIntPop,
						OutIntRange = 0,
						ReturnPopRaster = T,
						CreatePopFile = F )
	
	
	##--- landscape params + demographic scalings
	
	# create a placeholder for demographic scaling layers
	empty_demogLayer <- habitat_matrices[[1]]
	empty_demogLayer[!is.na(empty_demogLayer)] <- 100
	# determine array dimensions
	dimX <<- ncol(habitat_matrices[[1]])
	dimY <<- nrow(habitat_matrices[[1]])
	nlayers <- 3
	dim_array = c(dimY,dimX,nlayers)
	# create array of dynamic scalings with empty placeholder layers
	tmp_scaleLayers <- rep(list(array(rep(empty_demogLayer,nlayers), dim = dim_array)),length(habitat_matrices))
	
	# landscape params
	land <- ImportedLandscape(LandscapeFile = habitat_matrices,
							  demogScaleLayers = tmp_scaleLayers, # list with one array per year change
							  OriginCoords = LLC_coords,
							  DynamicLandYears = c(0,(nrYrBurnin+1):(nrYrPred+nrYrBurnin-1)), # years = land@DynamicLandYears+1995
							  HabPercent = TRUE,
							  Resolution = Resolution,
							  K_or_DensDep = ref_par["DensDep"] )
	
	##--- demographic params
	
	nr_stages <<- 2
	
	demog <- Demography(ReproductionType = 0,
						StageStruct = StageStructure(Stages = nr_stages,
													 TransMatrix = matrix(c(           0             , maxDemogRates["Fecund"],
													 								  maxDemogRates["juvSurv"] , maxDemogRates["adSurv"]),
													 					 nrow = nr_stages,
													 					 byrow = T),
													 #MaxAge = 12,
													 FecDensDep = TRUE,
													 FecStageWtsMatrix = matrix(c( 0,   0,
													 							  0, 1.0),
													 						   nrow = nr_stages, byrow = T),
													 FecLayer = c(NA,1), # matrix(c(NA,1)),
													 #DevLayer = matrix(c(NA,NA)),
													 SurvLayer = c(2,3) #matrix(c(2,3))
						)
	)
	
	##--- dispersal params
	
	disp <- Dispersal(Emigration = Emigration(StageDep = TRUE,
											  EmigProb = matrix(c(0:(nr_stages-1),ref_par["EmigProb"],0), nrow = nr_stages )),
					  Transfer   = DispersalKernel(Distances = matrix(1000*ref_par["DispDist"])),
					  Settlement = Settlement(Settle = 3))
	
	
	##--- initialisation params
	
	init <- Initialise(InitType = 2, 
					   InitIndsFile = "NULL", # "initinds.txt", # List init
					   InitIndsList = replicate(simul@Replicates, expr = createInitInds(NULL), simplify = FALSE )
	)
	
	##--- param master
	
	s <- RSsim(batchnum = batchnr, land = land, demog = demog, dispersal = disp, init = init, simul = simul)
	return(s)
}
