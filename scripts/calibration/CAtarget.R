
###-----------------------------------------------------------------
##   CALIBRATION TARGET
###-----------------------------------------------------------------


# The function CAtarget() defines the likelihood, our calibration target.


CAtarget <- function (params) {
	
	gof <- -Inf
	out = runRSmodel(params)
	if (class(out) == "RasterStack") {
		out <- stackApply(out, indices = as.integer(sapply(strsplit(names(out), "year"), "[[", 2)) + 1 - nrYrBurnin, fun = mean)
		sim_data <- raster::extract(out[[ix_years[ix_years %in% seq(nlayers(out))]]], sample_cells_abd$cell, df = TRUE)
		sim_data$ID <- sample_cells_abd$blockID
		#names(sim_data)[names(sim_data) == "ID"] <- "blockID"
		notsim_years <- sum(!ix_years %in% seq(nlayers(out)))
		if (notsim_years) sim_data <- cbind(sim_data, data.frame(matrix(0, ncol = notsim_years, nrow = nrow(sample_cells_abd))))
		sim_data <- aggregate(sim_data, by = list(sim_data$ID), FUN = mean )[-c(1,2)]
		sim_data[sim_data == 0] <- sigma
		GPsize <- params[GPsize_parix]
		gof <- sum(dnbinom(x = round(as.matrix(MHB_counts)), # round() when temporal aggregation 
						   mu = as.matrix(sim_data), size = GPsize, log = TRUE), 
				   na.rm = TRUE)
	}
	else return(set_calib$InfVal)
	return(gof)
}
