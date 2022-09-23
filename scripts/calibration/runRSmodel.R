
###-----------------------------------------------------------------
##   RUN RANGESHIFTER SIMULATION
###-----------------------------------------------------------------


# The function runRSmodel() runs the RangeShifter simulation for a given set of parameters.

# It modifies the simulation parameters to be calibrated, validates the new parameter settings,
# draws a new seed, creates a new initialisation scenario, calculates the demographic scaling
# layers from the climate layers and the demography-environment parameters, runs the simulation,
# and returns the output as a raster stack.


runRSmodel <- function(params){
	
	##-- 0.) get set of params
	
	p <- ref_par
	p[par_sel] <- params
	
	
	##-- 1.) reset parameters
	
	s@land@K_or_DensDep <- p["DensDep"]
	
	s@dispersal@Emigration@EmigProb[1,2] <- p["EmigProb"]
	
	s@dispersal@Transfer@Distances <- 1000*matrix(p["DispDist"])
	
	s@init@InitIndsList <- replicate(s@simul@Replicates, 
									 expr = createInitInds(init_poisMeans),
									 simplify = FALSE )
	
	
	##-- 2.) calculate demography layers
	
	# means of demographic rates
	meandemog <- c( p['Fecund'] , p['juvSurv'] , p['adSurv'] )
	location_mu <- c(log(meandemog['Fecund']), log((1-meandemog[2:3])/meandemog[2:3]))
	
	# check number of layers in demog. scalings / habitat maps and climate layers
	if(length(s@land@demogScaleLayers)!=length(yearly_envStacks_standard)) stop("differing number of years in habitat maps and climate maps")
	
	# calculate absolute values of demographic traits
	for(year in 1:length(s@land@demogScaleLayers)) {
		#print(year)
		climate <- yearly_envStacks_standard[[year]]
		# fecundity -> log link
		#s@land@demogScaleLayers[[year]][,,1] <- 100 * matrix(  exp( location_mu['Fecund'] + p['Fcb.pr1']*climate$breeding_pr + p['Fcb.pr2']*(climate$breeding_pr)^2 + p['Fcb.tm1']*climate$breeding_tmp + p['Fcb.tm2']*(climate$breeding_tmp)^2 ) / maxDemogRates['Fecund'] , ncol = dimX, nrow = dimY, byrow = TRUE)
		s@land@demogScaleLayers[[year]][,,1] <-  100 * matrix( min( exp( location_mu['Fecund'] + p['Fcb.pr1']*climate$breeding_pr + p['Fcb.pr2']*(climate$breeding_pr)^2 + p['Fcb.tm1']*climate$breeding_tmp + p['Fcb.tm2']*(climate$breeding_tmp)^2 ) / maxDemogRates['Fecund'], 1 ) , ncol = dimX, nrow = dimY, byrow = TRUE)
		# juvenile survival -> logit link
		s@land@demogScaleLayers[[year]][,,2] <- 100 * plogis( matrix( p['jSb.pr1']*climate$winter_pr + p['jSb.pr2']*(climate$winter_pr)^2 + p['jSb.tm1']*climate$postbreed_tmp + p['jSb.tm2']*(climate$postbreed_tmp)^2, ncol = dimX, nrow = dimY, byrow = TRUE), location = location_mu['juvSurv'] )
		# adult survival -> logit link
		s@land@demogScaleLayers[[year]][,,3] <- 100 * plogis( matrix( p['aSb.pr1']*climate$winter_pr + p['aSb.pr2']*(climate$winter_pr)^2 + p['aSb.tm']*climate$winter_pr, ncol = dimX, nrow = dimY, byrow = TRUE) , location = location_mu['juvSurv'] )
	}
	
	
	##-- 3.) Set seed and validate RS master
	
	s@control@seed <- runif(1, 1000000, 9999999)
	# assign modified parameter master to global environment
	#assign("s", s, envir = .GlobalEnv)
	
	validateRS <- FALSE
	try(validateRS <- validateRSparams(s), silent = TRUE)
	if(validateRS)
		capture.output(
			out <- RunRS(s, dirpath = RSdir) #,type = c("output", "message")
		)
	else out <- NA
	
	return(out)
}
