
###-----------------------------------------------------------------
##   READ AND PREPARE DATA FOR SIMULATION AND CALIBRATION
###-----------------------------------------------------------------


# The function readandpreparedata() reads the data necessary for the simulation and 
# calibration and assigns them to the parent (global) environment.

# First, the prepared data for the creation of the initialisation scenarios is read.
# Then, the habitat maps are read from ASCII files and turned into matrices for use in RangeShiftR.
# Next, the observations are read and processes according to optional spatial and temporal aggregation.
# Finally, the indices of the observed cells and years are recorded for the comparison of 
# observed and simulated data in the calibration step.


readandpreparedata <- function(){
	
	
	#--- 1.) get species traits
	
	load(file = "data/species_traits/species_traits.Rdata")
	trt <<- species_traits[species_traits$species_name==eval(species_name),]
	rm(species_traits)
	
	
	#--- 2.) get initialisation data
	
	init_poisMeans <<- read.csv(paste0("data/initpop/init_PoissonMeans/",species_name,".csv"))
	
	
	#--- 3.) get MHB data
	
	load(file = "data/species_counts/MHB_aggregated_counts.Rdata")
	sample_cells_abd <<- MHB_counts_sp[[eval(species_name)]]$sample_cells_abd
	MHB_counts <<- MHB_counts_sp[[eval(species_name)]]$MHB_counts
	rm(MHB_counts_sp)
	
	# years to be compared
	mhb_years <<- (0:(nrYrMHB-1))+nrYrBurnin   # MHB-sampled years to be compared with simulation
	simul_years <<- seq(from = OutStartPop, by = OutIntPop, to = nrYears-1)  # recorded simulated years
	ix_years <<- which(simul_years %in% mhb_years)   # indices of sampled years within vector of simulated years
	
	
	#--- 4.) get habitat maps
	
	habitat_names <- paste0("data/habitatmaps/RS_maps/habitat_",species_name,"_",sprintf("%02d",years_ca),".asc")
	habitatmap <- raster(habitat_names[1])
	habitat_matrices <- lapply(habitat_names, 
							   FUN = function(filename){
							   	matrix(raster(filename), 
							   		   ncol = ncol(habitatmap), 
							   		   nrow = nrow(habitatmap), byrow = TRUE)})
	names(habitat_matrices) <- paste0("year_",sprintf("%02d",years_ca))

	Resolution <<- res(habitatmap)[1]
	LLC_coords <<- c(xmin(habitatmap),ymin(habitatmap))
	nRow_map   <<- nrow(habitatmap)
	assign("habitat_matrices", habitat_matrices, envir = .GlobalEnv)
	
	
	#--- 5.) get climate maps
	
	load(envir = .GlobalEnv, file = "data/climatemaps/aggrVariables_zTrafod.Rdata")
	# winter precipitation of year 2019 is missing. replace with that of 2018 because last years survival doesn't count into abundance anyway
	yearly_envStacks_standard[[21]] <<- addLayer(yearly_envStacks_standard[[21]],yearly_envStacks_standard[[20]]$winter_pr)
	
}
