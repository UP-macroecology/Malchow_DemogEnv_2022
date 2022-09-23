
###-----------------------------------------------------------------
##   PREPARE HABITATMAPS
###-----------------------------------------------------------------


# This script creates the habitat maps for all species and years (1999-2020), based on their habitat preferences and CORINE land cover data.


### 0.) Load packages and data 

library(raster)

# load species habitat preferences
# (use data_prep/prep_speciestraits.R to create this object)
load("data/habitatmaps/species_preferences/habitat_traits.Rdata")

# load CORINE land cover layers of 2018
# (use data_prep/prep_CORINElayers.R to create this object)
load("data/habitatmaps/landcover/corine_layers.Rdata") 


### 1.) Exclude land cover classes with less than 1% cover in Switzerland

# number of populated cells within each layer
nr_cells <- lapply(corine_layers, function(corine_yr){sapply(1:nlayers(corine_yr), function(i){sum(values(corine_yr[[i]]),na.rm=T)})})
# total number of cells
total_cells <- ncell(corine_layers[[4]][[44]])
# percentage cover of each layer
perc_cover <- lapply(nr_cells, '/', total_cells)
# select those land cover classes with at least 1% cover over all of Switzerland
clc_over1perc <- lapply(perc_cover,  function(yr){which(yr>0.01 )})
# check if same lc classes for all CORINE years:
table(unlist(clc_over1perc)) # yes, all appear exactly 4 times -> take one of the years
# get CLC names
clc_over1perc <- names(corine_layers[[1]])[clc_over1perc[[1]]]


### 2.) Mapping preference habitat classes to land cover classes

# create empty data frame to hold preferred CORINE land cover classes per species
CORINE_classes <- data.frame(habitat_traits$Species, matrix(0,nrow=nrow(habitat_traits),ncol=length(clc_over1perc)))
names(CORINE_classes) <- c("Species",clc_over1perc)

# mapping 
CORINE_classes$CLC_112 <- habitat_traits$Human.settlements    # Discontinuous urban fabric
#CORINE_classes$CLC_211    # Non-irrigated arable land
#CORINE_classes$CLC_231    # Pastures
#CORINE_classes$CLC_242    # Complex cultivation patterns
CORINE_classes$CLC_311 <- habitat_traits$Deciduous.forest	 # Broad-leaved forest
CORINE_classes$CLC_312 <- habitat_traits$Coniferous.forest    # Coniferous forest
CORINE_classes$CLC_313 <- habitat_traits$Deciduous.forest    # Mixed forest  --  (because all species in Deciduous.forest are also in Coniferous.forest)
CORINE_classes$CLC_321 <- habitat_traits$Shrub + habitat_traits$Woodland    # Natural grasslands  --  (because all species in Woodland are also in Shrub and vice versa)
CORINE_classes$CLC_322 <- habitat_traits$Mountain.meadows + habitat_traits$Grassland     # Moors and heathland
CORINE_classes$CLC_324 <- habitat_traits$Shrub + habitat_traits$Woodland    # Transitional woodland-shrub
CORINE_classes$CLC_332 <- habitat_traits$Rocks    # Bare rocks
CORINE_classes$CLC_333 <- habitat_traits$Mountain.meadows + habitat_traits$Grassland   # Sparsely vegetated areas  --  (+Grassland because this only concerns Turdus torquatus and it has a 0 for Mountain.meadows)
#CORINE_classes$CLC_335   # Glaciers and perpetual snow
#CORINE_classes$CLC_512   # Water bodies

# after comparison with maps from Vogelwarte.ch, adjust habitat preference for Regulus regulus by adding Mixed forest (CLC_313)
CORINE_classes[CORINE_classes$Species=="Regulus regulus","CLC_313"] <- 1
# ... and for Linaria cannabina by adding Pastures (CLC_231)
CORINE_classes[CORINE_classes$Species=="Linaria cannabina","CLC_231"] <- 1

# save(CORINE_classes, file = "data/habitatmaps/species_preferences/landcover_species.Rdata")


### 2.) Create habitat maps

# create empty stack
empty_stack <- stack( replicate(nrow(habitat_traits), corine_layers[[1]][[44]]))
names(empty_stack) <- habitat_traits$Species

# create habitat maps by collecting presences of habitat in cells
habitatmaps <- lapply(corine_layers, function(corine_year){
	this_stack <- empty_stack
	for( sp in 1:nlayers(this_stack) ){
		for (clc in colnames(CORINE_classes)[-1]) {
			if(CORINE_classes[sp,clc]){
				values(this_stack[[sp]]) <- values(this_stack[[sp]]) + values(corine_year[[clc]])
			}
		}
	}
	return(this_stack)
})

# sum up number of 1-ha habitat cells per 1-km^2
habitatmaps_1km <- lapply(habitatmaps, function(habitat_year){
	aggregate(habitat_year, 10, sum)
})

# Re-project to Swiss coordinates - use background mask for this
bg <- raster('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')
habitatmaps_1km <- lapply(habitatmaps_1km, projectRaster, bg, method = "ngb")
habitatmaps_1km <- lapply(habitatmaps_1km, resample, bg)

spplot(habitatmaps_1km[[4]])

# save maps at this point, still (almost) rectangular, i.e. not masked to Swiss border
	#save(habitatmaps_1km, file = "data/habitatmaps/habitatmaps_nomask.Rdata")
	load(file = "data/habitatmaps/habitatmaps_nomask.Rdata")


### 3.) Generate RangeShifter maps

# now, apply Swiss border plus a buffer of 12km
bg_12k <- buffer(bg, 12090) # need to add a little because of the reprojection
habitatmaps_1km <- lapply(habitatmaps_1km, mask, bg_12k)

# show difference in habitat amount between 2000 and 2018
spplot(habitatmaps_1km[[4]]-habitatmaps_1km[[1]],col.regions = hcl.colors(30, palette = "Purple-Green", rev = FALSE), at = seq(-104,104, length.out = 30))


# interpolate to yearly time steps

# 1. interval 1999-2005
hab_diff <- habitatmaps_1km$corine2006 - habitatmaps_1km$corine2000
habitatmaps_1km$corine1999 <- habitatmaps_1km$corine2000-(1/6*hab_diff)
habitatmaps_1km$corine2001 <- habitatmaps_1km$corine2000+(1/6*hab_diff)
habitatmaps_1km$corine2002 <- habitatmaps_1km$corine2000+(2/6*hab_diff)
habitatmaps_1km$corine2003 <- habitatmaps_1km$corine2000+(3/6*hab_diff)
habitatmaps_1km$corine2004 <- habitatmaps_1km$corine2000+(4/6*hab_diff)
habitatmaps_1km$corine2005 <- habitatmaps_1km$corine2000+(5/6*hab_diff)
# 2. interval 2007-2011
hab_diff <- habitatmaps_1km$corine2012 - habitatmaps_1km$corine2006
habitatmaps_1km$corine2007 <- habitatmaps_1km$corine2006+(1/6*hab_diff)
habitatmaps_1km$corine2008 <- habitatmaps_1km$corine2006+(2/6*hab_diff)
habitatmaps_1km$corine2009 <- habitatmaps_1km$corine2006+(3/6*hab_diff)
habitatmaps_1km$corine2010 <- habitatmaps_1km$corine2006+(4/6*hab_diff)
habitatmaps_1km$corine2011 <- habitatmaps_1km$corine2006+(5/6*hab_diff)
# 3. interval 2013-2017
hab_diff <- habitatmaps_1km$corine2018 - habitatmaps_1km$corine2012
habitatmaps_1km$corine2013 <- habitatmaps_1km$corine2012+(1/6*hab_diff)
habitatmaps_1km$corine2014 <- habitatmaps_1km$corine2012+(2/6*hab_diff)
habitatmaps_1km$corine2015 <- habitatmaps_1km$corine2012+(3/6*hab_diff)
habitatmaps_1km$corine2016 <- habitatmaps_1km$corine2012+(4/6*hab_diff)
habitatmaps_1km$corine2017 <- habitatmaps_1km$corine2012+(5/6*hab_diff)
# 4. interval 2019-2020
hab_diff <- habitatmaps_1km$corine2018 - habitatmaps_1km$corine2000
habitatmaps_1km$corine2019 <- habitatmaps_1km$corine2018+(1/18*hab_diff)
habitatmaps_1km$corine2020 <- habitatmaps_1km$corine2018+(2/18*hab_diff)

# for the extrapolated years (1999, 2019, 2020) restrict habitat values to [0,100]
extrapolated_years <- c("corine1999","corine2019","corine2020")
for(year in extrapolated_years){
	for(sp in names(habitatmaps_1km$corine2000)) {
			neg_cells <- which(values(habitatmaps_1km[[year]][[sp]]) < 0 )
			values(habitatmaps_1km[[year]][[sp]])[neg_cells] <- 0
			pos_cells <- which(values(habitatmaps_1km[[year]][[sp]]) > 100 )
			values(habitatmaps_1km[[year]][[sp]])[pos_cells] <- 100
	}
}

# time series of habitat amount per species
habitat_sums <- lapply(habitatmaps_1km, function(year_stack){
	sapply(1:nlayers(year_stack), function(yr_layer){ sum(values(year_stack[[yr_layer]]), na.rm=T)})
})
habitat_sums <- as.data.frame(habitat_sums)
habitat_sums <- habitat_sums[,sort(names(habitat_sums))]
row.names(habitat_sums) <- habitat_traits$Species

years <- 1999:2020
par(mar=c(4,4,0.5,0.5))
plot(NULL,xlim=range(years),ylim=range(habitat_sums),xlab="year",ylab="habitat amount [ha]")
for(sp in 1:nrow(habitat_sums)) {
	lines(years,habitat_sums[sp,], type="b",col=sp)
	if(sp==2) xpos <- 2006.6 else if(sp==3) xpos <- 2010.8 else if(sp==5) xpos <- 2018.5 else xpos <- 2015
	if(sp%in%c(2:5)) yoff <- -40000 else yoff <- 44000
	text(xpos,habitat_sums[sp,16]+yoff, row.names(habitat_sums)[sp], col=sp)
}


# save finished maps as Rdata object
 #save(habitatmaps_1km, file = "data/habitatmaps/habitatmaps_masked.Rdata")
load(file = "data/habitatmaps/habitatmaps_masked.Rdata")

# save as ASCII rasters
lapply(names(habitatmaps_1km), function(year){
	sapply(names(habitatmaps_1km$corine2000), function(sp){
		writeRaster((round(habitatmaps_1km[[year]][[sp]], digits = 0)),
					format="ascii",
					filename = paste0("data/habitatmaps/RS_maps/habitat_",sp,"_",strsplit(year,"corine")[[1]][2],".asc"),
					NAflag = -9L,
					datatype='INT1S',
					overwrite = T )
	})
})

