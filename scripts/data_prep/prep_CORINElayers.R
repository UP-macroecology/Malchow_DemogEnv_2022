
###-----------------------------------------------------------------
##   PREPARE HABITATMAPS
###-----------------------------------------------------------------


# This script reads the CORINE land cover data from the years 2000, 2006, 2012, 2018; crops them to Swiss extend, prepares a single layer
# for each land cover class and year, reprojects these to Swiss coordinates, and saves them as a stack.
# These land cover maps are used to create habitat maps from the species' habitat preferences in the script prep_habitatmaps.R.


library(raster)

# Retrieve CORINE data from COPERNICUS

dl_path <- NULL # set download path here

corine2018 <- raster(paste0(dl_path,'/CLC2018/DATA/U2018_CLC2018_V2020_20u1.tif'))
corine2012 <- raster(paste0(dl_path,'/CLC2012/DATA/U2018_CLC2012_V2020_20u1.tif'))
corine2006 <- raster(paste0(dl_path,'/CLC2006/DATA/U2012_CLC2006_V2020_20u1.tif'))
corine2000 <- raster(paste0(dl_path,'/CLC2000/DATA/U2006_CLC2000_V2020_20u1.tif'))

corine_stack <- stack(corine2000,corine2006,corine2012,corine2018)
names(corine_stack) <- c("corine2000","corine2006","corine2012","corine2018")

# Crop
ext <- c(4000000, 4400000, 2500000,2750000)
corine_stack <- crop(corine_stack, ext)
rm(corine2000,corine2006,corine2012,corine2018)
spplot(corine_stack)

# Create stack of empty layers to store values in 
corine_emptylayer <- corine_stack[[1]]
values(corine_emptylayer) <- 0
corine_emptystack <- stack(replicate(44,corine_emptylayer)) # 44 = nr of LC classes in CORINE
names(corine_emptystack) <- c('CLC_111', 'CLC_112', 'CLC_121', 'CLC_122', 'CLC_123', 'CLC_124', 'CLC_131', 'CLC_132', 'CLC_133', 'CLC_141', 'CLC_142', 'CLC_211', 'CLC_212', 'CLC_213', 'CLC_221', 'CLC_222', 'CLC_223', 'CLC_231', 'CLC_241', 'CLC_242', 'CLC_243', 'CLC_244', 'CLC_311', 'CLC_312', 'CLC_313', 'CLC_321', 'CLC_322', 'CLC_323', 'CLC_324', 'CLC_331', 'CLC_332', 'CLC_333', 'CLC_334', 'CLC_335', 'CLC_411', 'CLC_412', 'CLC_421', 'CLC_422', 'CLC_423', 'CLC_511', 'CLC_512', 'CLC_521', 'CLC_522', 'CLC_523')

# Make one stack of layers of land cover classes per CORINE year
corine_layers <- lapply(1:nlayers(corine_stack), function(yr){
	this_stack <- corine_emptystack
	for (i in 1:44) values(this_stack[[i]]) <- (values(corine_stack[[yr]])==i)
	return(this_stack)
})
names(corine_layers) <- c("corine2000","corine2006","corine2012","corine2018")

# Save to Rdata object
save(corine_layers, file = "data/habitatmaps/landcover/corine_layers.Rdata")
