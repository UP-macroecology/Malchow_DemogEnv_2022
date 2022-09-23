
###-----------------------------------------------------------------
##   CREATE INITIAL INDIVIDUALS LIST
###-----------------------------------------------------------------


# The function createInitInds() creates a list of individuals for initialisation as returns it as a data frame.

# This object is needed for the initialisation option 'From initial individuals list file' (InitType = 2).
# Setting 'InitIndsFile = "NULL"' lets RangeShiftR expect a data frame instead of a file to read.
# This data frame is then given through the parameter 'InitIndsList'.
# This function is meant to be called from the calibration_main.R script, in which the required objects 
# that contain the initialisation data are loaded.


createInitInds <- function(pois_mean) {

	if(is.null(pois_mean)){
		return(data.frame(Year = integer(), Species = integer(), X = integer(), Y = integer(), Ninds = integer(), Age = integer(), Stage = integer()) )
	}

	# draw from Poisson
	pois_mean$pois_draw <- rpois(nrow(pois_mean), lambda = pois_mean$pois_mean)
	#pois_mean$habitat <- extract(habitatmaps, cellFromRowCol(habitatmaps, row=pois_mean$row, col=pois_mean$col))
	
	# discard zero draws
	pois_mean <- pois_mean[pois_mean$pois_draw>0,]
	
	# prepare initialisation data frame
	# note that x,y coordinates start counting at 0, and that y-coordinate is reversed:
	ind_table <- data.frame(Year = 0, Species = 0, X = pois_mean[['col']]-1, Y = nRow_map-pois_mean[['row']], Ninds = pois_mean[['pois_draw']], Age = 1, Stage = 1)

	return(ind_table)
}

# 
# for(species_nr in c(5,8,9,10)){ #1:10
# 	
# 	#species_nr <- 8
# 	
# 	mhb_data <- cbind(rowColFromCell(ref_map,mhb_data$cell),mhb_data)
# 	ertr <- mhb_data[,c(1:3,4)]
# 	names(ertr) <- c('row_mhb','col_mhb','cell_mhb','count_mhb')
# 	ertr <- cbind(ertr,habitat_mhb=extract(ref_map,mhb_data$cell))
# 	names(init_poisMeans) <- c('row_ini','col_ini','pois_mean','habitat_ini')
# 	ertrr <- merge(ertr,init_poisMeans, by.x = c('row_mhb','col_mhb'), by.y = c('row_ini','col_ini'))
# 	head(ertrr)
# 	
# 	rty <- ertrr[ertrr$habitat_mhb==0,]
# 	ab <- sum(rty$count_mhb,na.rm=T)/nrow(rty)
# 	ba <- sum(ertrr$count_mhb,na.rm=T)/nrow(ertrr)
# 	print( paste(ab,",",ba) )
# 	print(  paste(ab/ba,",",nrow(rty)/nrow(ertrr)) )
# 	print(  sum(rty$count_mhb,na.rm=T)/sum(ertrr$count_mhb,na.rm=T)  )
# 	
# 	values(ref_map)[mhb_data[mhb_data$X1999>0,'cell']] <- 300
# 	values(ref_map)[rty[rty$count_mhb>0,'cell_mhb']] <- -200
# 	plot(ref_map, main = species_name, col = hcl.colors(20, "Plasma", rev=T))
# 
# }
