
load(file = "results/analysis/rates_maps/spatial_mean/growthrate_TS.Rdata")

# Plot time series of mean growth rate for all species in one plot

years <- 1999:2018

cols <- cols_specID
col_shade <- sapply(cols,t_col,percent = 50)
col_line <- sapply(cols,makeColDarker)

tikz(file = "results/analysis/rates_maps/spatial_mean/meanGrowthRate_allSpec.tex", width = 3.84, height = 2.6, standAlone = FALSE) #, width = 5.76, height = 4.84
{
	par( cex = 0.9, mgp = c(2,0.6,0), mar=c(3.4,3.2,0.1,0.1) )
	plot(NULL, xlim = range(years), ylim = c(0.4,1.2), xlab = "Year" , ylab = "Mean growth rate") 
	for(specID in rev(sel_specs)){
		if(specID!=10){
			polygon(x = c(years,rev(years)), y = c(growth_TS[["75%"]][-21,specID],rev(growth_TS[["25%"]][-21,specID])),  col = col_shade[specID], border = NA)
			lines(x = years, y = growth_TS[["50%"]][-21,specID],  col = col_line[specID], lwd = 2, type = "b" )
		}
	}
}
dev.off()
