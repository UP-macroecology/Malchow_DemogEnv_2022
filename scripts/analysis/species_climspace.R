
###---------------------------------------------------------------------------------------
##   RECORD CLIMATE SPACE COVERED BY EACH SPECIES
###---------------------------------------------------------------------------------------


library(raster)


## 1.) Species occurrences

# Read in species counts
load("data/species_counts/MHBabunds_1999-2021.Rdata")
MHBcells <- species_counts[[1]]$cell

# Coerce to presence/absence
species_counts <- lapply(species_counts, function(sp){
	sp[-1] <- sp[-1]>0
	sp[is.na(sp)] <- FALSE
	sp
})


## 2.) Climate variables

# Read in yearly stacks of climate variables
load("data/climatemaps/aggrVariables_zTrafod.Rdata")

# Extract climate variables at MHB sites
clim <- lapply(yearly_envStacks_standard, extract, y=MHBcells)


## 3.) Habitat suitabilities

# Read in yearly stacks habitat suitabilites for each species
load("data/habitatmaps/habitatmaps_masked.Rdata")

# record years in each data set
hab_years  <- as.integer(unlist(strsplit( names(habitatmaps_1km), "corine"))[seq(2,2* length(habitatmaps_1km),2)])
clim_years <- as.integer(unlist(strsplit( names(clim), "year_"))[seq(2,2* length(clim),2)])

spec_habitat <- 
	lapply( names(habitatmaps_1km[[1]]), function(spec){
		print(spec)
		yearmaps <- lapply( habitatmaps_1km, '[[', spec)
		yearlists <- lapply(yearmaps, extract, y=MHBcells)
		yearlists[order(hab_years[which(hab_years %in% clim_years)])]
	})
names(spec_habitat) <- names(habitatmaps_1km[[1]])




## 4.) Collect values of all variables at all presences for each year

clim_presences <- lapply(1:length(species_counts), function(spec){
	
	sp_name <- names(species_counts)[spec]
	sp_counts <- species_counts[[spec]]
	
	# get species habitat
	sp_hab <- spec_habitat[[sp_name]]
		
	# prepare empty data frame
	out1 <- cbind( year = NULL, clim[[1]][FALSE,] )
	out2 <- NULL

	for(yearID in 1:length(clim)){
		
		# collect climate predictors
		if(yearID<21) out1 <- rbind(out1,cbind( yearID, clim[[yearID]][sp_counts[,yearID+1],]))
		else out1 <- rbind(out1, cbind( yearID, clim[[yearID]][sp_counts[,yearID+1],], winter_pr=NA))
		
		# collect habitat suitability
		out2 <- c(out2,sp_hab[[yearID]][sp_counts[,yearID+1]])
	}
	
	cbind(out1, hab_suit = out2)
})
names(clim_presences) <- names(species_counts)

# Save as object
save(clim_presences, file = "data/climatemaps/aggrVar_atPresences.Rdata")



## 5.) Calculate quantiles

calc_quantiles <- c(0.10,0.20,0.25,0.50,0.75,0.80,0.90)


spec_preds_quants <- lapply(1:length(species_counts), function(spec){
	
	sp_name <- names(species_counts)[spec]
	sp_preds <- clim_presences[[spec]]
	
	# calculate quartiles
	apply(sp_preds, 2, quantile, probs = calc_quantiles, na.rm = TRUE )

})
names(spec_preds_quants) <- names(species_counts)

# save quantiles
save(spec_preds_quants, file = paste0("results/analysis/species_preds_quantiles.Rdata"))


