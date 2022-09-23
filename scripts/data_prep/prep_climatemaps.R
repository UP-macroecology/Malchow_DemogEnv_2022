
###-----------------------------------------------------------------
##   PREPARE HABITATMAPS
###-----------------------------------------------------------------


# This script creates the climate maps for all years (1999-2019), based on the Chelsa V2.1 climatologies.

# The aggregated climatic variables are the following:

# breeding season (April – July):
# 	mean daily air temperature
# 	mean precipitation
 
# after breeding season (September - November):	
# 	mean daily air temperature
 	
# winter (December - February):
# 	minimum air temperature
# 	mean precipitation
 


### 0.) Load packages and data 

library(raster)

# Retrieve CORINE data 
mt_path <- NULL  # path to downloaded Chelsa v2.1 monthly layers, cropped to extent of Switzerland

# Get Swiss border and apply 12km-buffer for masking
bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')
bg_12k <- buffer(bg,12090)


### 1.) Generate target climatic variables per year (in absolute values)

# Loop over years
years <- 1999:2019

if(FALSE){ # comment out loop, instead load saved results
	yearly_envStacks <- sapply(years, function(year){
		
		## Calculate climate layers
		print(year)
		# a) Daily mean temperature
		stk_tas <- brick(paste0(mt_path,'CHELSA_tas_',year,'_V.2.1.tif'))
		
		averageed_tas <- stackApply(stk_tas, indices = c(1,1,1,2,2,2,2,1,3,3,3,1), fun=mean)
		breeding_tas <- averageed_tas[[2]]
		postbreed_tas <- averageed_tas[[3]]
		
		# b) minimum air temperature
		stk_tasmin1 <- brick(paste0(mt_path,'CHELSA_tasmin_',year,'_V.2.1.tif'))
		if(year<2019){
			stk_tasmin2 <- brick(paste0(mt_path,'CHELSA_tasmin_',year+1,'_V.2.1.tif'))
			
			winter_tasmin <- mean(stack(stk_tasmin1[[12]],stk_tasmin2[[1:2]]))
		}else{
			winter_tasmin <- stk_tasmin1[[12]]
		}
		
		# c) precipitation
		stk_pr1 <- brick(paste0(mt_path,'CHELSA_pr_',year,  '_V.2.1.tif'))
		if(year<2019){
			stk_pr2 <- brick(paste0(mt_path,'CHELSA_pr_',year+1,'_V.2.1.tif'))
			
			breeding_pr <- mean(stk_pr1[[4:7]])
			winter_pr <- mean(stack(stk_pr1[[12]],stk_pr2[[1:2]]))
		}else{ # year 2019
			breeding_pr <- mean(stk_pr1[[4:6]])
			winter_pr <- NA
		}
		
		## Make stack, re-transform temperatures (are given by Chelsa in units of 1/10 Kelvin), and mask to border
		year_stack <- stack(breeding_tas/10, breeding_pr, postbreed_tas/10, winter_tasmin/10, winter_pr)
		year_stack <- mask(year_stack, bg_12k)
		if(year<2019){names(year_stack) <- c("breeding_tmp", "breeding_pr", "postbreed_tmp", "winter_tmin", "winter_pr")
		}else{names(year_stack) <- c("breeding_tmp", "breeding_pr", "postbreed_tmp", "winter_tmin")}
		
		return(year_stack)
	})
	names(yearly_envStacks) <- paste0("year_",years)
	#save(yearly_envStacks, file = "data/climatemaps/aggrVariables.Rdata")
}else{
	load(file = "data/climatemaps/aggrVariables.Rdata")
}



### 2.) Standardize values of climatic variables

## 2.1.) Calculate mean and variance over all years

# calculate means and square means over space for each climatic variable, fill up NA for precipitation variables in year 2019

envVars_yearly_means <- sapply(yearly_envStacks, function(year_stack){
	tmp <- sapply(1:nlayers(year_stack), function(lyr){ mean(values(year_stack[[lyr]]), na.rm=T) } )
	names(tmp) <- names(year_stack)
	return(tmp)
})
envVars_yearly_means$year_2019 <- c(envVars_yearly_means$year_2019, winter_pr=NA)
#envVars_yearly_means$year_2019 <- envVars_yearly_means$year_2019[names(envVars_yearly_means$year_2018)]

envVars_yearly_squaremeans <- sapply(yearly_envStacks, function(year_stack){
	tmp <- sapply(1:nlayers(year_stack), function(lyr){ mean(values(year_stack[[lyr]])^2, na.rm=T) } )
	names(tmp) <- names(year_stack)
	return(tmp)
})
envVars_yearly_squaremeans$year_2019 <- c(envVars_yearly_squaremeans$year_2019, winter_pr=NA)
#envVars_yearly_squaremeans$year_2019 <- envVars_yearly_squaremeans$year_2019[names(envVars_yearly_squaremeans$year_2018)]


# from this, calculate overall means and square means (spatial and temporal), calculate variance and SD

envVars_tot_means <- rowMeans(as.data.frame(envVars_yearly_means), na.rm = T)
envVars_tot_squaremeans <- rowMeans(as.data.frame(envVars_yearly_squaremeans), na.rm = T)

envVars_tot_variances <- envVars_tot_squaremeans - (envVars_tot_means)^2

envVars_tot_stdev <- sqrt(envVars_tot_variances)

save(envVars_tot_means, envVars_tot_stdev, file = "data/climatemaps/zTrafo.Rdata")

## 2.2.) Standardise climate variables using z-transform

yearly_envStacks_standard <- yearly_envStacks
for (year in 1:length(yearly_envStacks)) {
	for (var in 1:nlayers(yearly_envStacks[[year]])) {
		var_name <- names(yearly_envStacks[[year]][[var]])
		yearly_envStacks_standard[[year]][[var]] <- (yearly_envStacks[[year]][[var]]-envVars_tot_means[var_name])/envVars_tot_stdev[var_name]
	}
}


### 3.) Save 

# as Rdata object
save(yearly_envStacks_standard, file = "data/climatemaps/aggrVariables_zTrafod.Rdata")

# as ASCII rasters
for (year in 1:length(yearly_envStacks_standard)) {
	for (var in 1:nlayers(yearly_envStacks_standard[[year]])) {
		var_name <- names(yearly_envStacks[[year]][[var]])
		writeRaster(round(yearly_envStacks_standard[[year]][[var]], digits = 3),
					format="ascii",
					filename = paste0("data/climatemaps/RS_maps/climate_",var_name,"_",years[year],".asc"),
					NAflag = -9L,
					#datatype='INT1S',
					overwrite = F )
	}
}

# also save mean and SDs used for z-transform as Rdata object to be able to re-transform
save(envVars_tot_means,envVars_tot_stdev, file = "data/climatemaps/zTransformPars.Rdata")



##-- 4.) Take out linear trend in climate predictors

# Read in yearly stacks of climate variables
#load("data/climatemaps/aggrVariables_zTrafod.Rdata")

# Calculate spatial mean for all years and each variable
clim_means <- sapply(yearly_envStacks_standard, function(yr_stack){
	out <- sapply(1:nlayers(yr_stack), function(var){mean(values(yr_stack[[var]]), na.rm=T) })
	if(length(out)==4) out[5] <- NA
	out
})

# Add year
clim_means <- rbind(1999:2019,clim_means)
# Make data frame
clim_means <- data.frame(t(clim_means))
names(clim_means) <- c("year",names(yearly_envStacks_standard[[1]]))

# Plot
#png(file = paste0("data/climatemaps/aggrVariables_zTrafod_means.png"), width = 720, height = 480)
	plot(NULL, xlim = range(clim_means[,1]), ylim = range(clim_means[,2:6], na.rm=T), main = "Climate predictors", xlab = "Year", ylab = "Spatial mean of z-Transform")
	for(i in 2:6) lines(clim_means[,1],clim_means[,i],col=i,lwd=2)
	legend("topleft", legend = names(clim_means)[2:6], col = 2:6, lwd=2)
#dev.off()

# Linear models
clim_means_lm <- sapply(2:6, function(var){with(clim_means,lm(as.formula(paste0(names(clim_means)[var],"~year"))))$coefficients})
colnames(clim_means_lm) <- names(yearly_envStacks_standard[[1]])

# Go through years & climate predictors and correct the linear trend
clim_CCcorrd <- yearly_envStacks_standard

for(yr in seq(length(clim_CCcorrd)) ){
	for(var in seq(nlayers(clim_CCcorrd[[1]])) ){
		if(yr!=21 || var!=5) values(clim_CCcorrd[[yr]][[var]]) <- values(clim_CCcorrd[[yr]][[var]])-(clim_means_lm["year",var]*yr)
	}
}

# Save CC-corrected climates as Rdata object
save(clim_CCcorrd,clim_means_lm, file = "data/climatemaps/aggrVariables_zTrafod_CCcorrected.Rdata")
