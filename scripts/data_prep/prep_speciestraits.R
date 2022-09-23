
###-----------------------------------------------------------------
##   PREPARE SPECIES TRAITS
###-----------------------------------------------------------------


# This script collects and stores the required demographic and dispersal traits for each species.

# Specifically, these are the survival probabilities of juveniles and adults, the fecundity, the emigration rate and mean dispersal distance.
# The values are used as prior means in the Bayesian calibration.
# The demographic traits are taken from:

# Storchová, Lenka, and David Hořák. "Life‐history characteristics of European birds." Global Ecology and Biogeography 27.4 (2018): 400-406.

# The dispersal traits are derived from exponential dispersal kernel means estimated from capture-recapture data. Because RangeShiftR requires
# that the mean dispersal distance is larger than the cell size (here 1000m), we calculate the emigration rate as the proportion of the
# probability mass of the kernel that exceeds this limit, and the mean dispersal distance as the mean of this part of the kernel.


### 0.) Species names

# names of our target species
species_names <- c("Certhia familiaris",
				   "Fringilla coelebs",
				   "Linaria cannabina",
				   "Lophophanes cristatus",
				   "Prunella collaris",
				   "Prunella modularis",
				   "Pyrrhula pyrrhula",
				   "Regulus regulus",
				   "Sitta europaea",
				   "Turdus torquatus")


### 1.) Demographic traits

# read table of traits downloaded from Storchová & Hořák (2018) SI
traits_df <- read.table("data/habitatmaps/species_preferences/Storchova_&_Horak_2018/Life-history characteristics of European birds.txt", sep = "\t", header = TRUE)

head(traits_df)


# the columns to be extracted (descriptions in the line comments are taken from the README of the SI)
traits_names <- c("Species",
				  "Clutch_MIN",              # Minimum clutch size [number of eggs]
				  "Clutch_MAX",              # Maximum clutch size [number of eggs]
				  "Clutch_MEAN",             # Mean clutch size [number of eggs]
				  "Broods.per.year",         # Mean number of broods per breeding season (replaced broods are not included) [number of broods]
				  "Nest.type",               # Type of nest [G = ground, on ground directly; H = hole, in tree, bank, ground, crevice; OA = open-arboreal, cup in bush, tree, on cliff ledge; CA = closed-arboreal; GC = ground close, nest in tussock very close to ground but not directly on ground, hidden in surrounded vegetation]
				  #"Mating system",          # Type of mating system [M = monogamous, PG = polygynous, PA = polyandrous, PM = promiscuous] # -> selected species all monogamous 
				  "Incubation.period",       # Mean length of eggs' incubation [days]
				  "Nestling.period",         # Average age of young when leaving nest [days]
				  "Fledging.period",         # Average age of young when fledging [days]
				  "Age.of.independence",     # Average age when young totally separate off parents - young are not fed or protected by parents [days]
				  "Age.of.first.breeding",   # Average age of the first breeding [years]
				  "Life.span",               # Maximum life span recorded in wild [years]
				  "Post.fledging.mortality", # Mean mortality of young in the first year of their life [%]
				  "Mortality.of.adults",     # Mean annual mortality of adults [%]
				  "Association.outside.the.breeding.season",  # Association of adults outside breeding season [GR = gregarious; PA = in pairs; SO = solitary]
				  "Territoriality",          # Defense of a territory (defended area occupied exclusively by a single bird, pair or larger social unit) [yes/no]
				  "Sedentary",               # Species lives in the same area in both breeding and non-breeding season [yes/no]
				  "Short.distance.migrant",  # Species migrates within the range of the Western Palearctic in non-breeding season [yes/no]
				  # Species occupies the following habitat in breeding area [yes/no]:
				  "Deciduous.forest", "Coniferous.forest", "Woodland", "Shrub", "Savanna", "Tundra", "Grassland", "Mountain.meadows", "Reed", "Swamps", "Desert", "Freshwater", "Marine", "Rocks", "Human.settlements")

# subset dataframe for selected species and traits
(species_traits <- traits_df[traits_df$Species %in% species_names, traits_names[c(1:5,13:14)]])

# Calculate fecundity as clutch size * #broods * nestling survival; we assume a nestling survival for all species of 50%.
# We further account for females only and assume a sex ratio of 1:1

nestling_survival <- 0.50
female_ratio <- 0.50

species_traits$fec_min  <- species_traits$Clutch_MIN  * species_traits$Broods.per.year * nestling_survival * female_ratio
species_traits$fec_max  <- species_traits$Clutch_MAX  * species_traits$Broods.per.year * nestling_survival * female_ratio
species_traits$fec_mean <- species_traits$Clutch_MEAN * species_traits$Broods.per.year * nestling_survival * female_ratio
species_traits$fec_sd   <- ((species_traits$fec_max - species_traits$fec_mean) + (species_traits$fec_mean - species_traits$fec_min)) / 4 # assuming the whole range respresents 2*sigma


# Post-fledging survival in accounted for in the juvenile survival probability.
# For species with no estimate of post-fledging survival, assume a mean prior of 50%:

species_traits[is.na(species_traits$Post.fledging.mortality),'Post.fledging.mortality'] <- 50

species_traits$jSv <- 1-species_traits$Post.fledging.mortality/100


# Lastly adult survival:
# For species with no estimate of adult survival, we assume a mean prior of 50%:

species_traits[is.na(species_traits$Mortality.of.adults),'Mortality.of.adults'] <- 50

species_traits$aSv <- 1-species_traits$Mortality.of.adults/100



### 2.) Dispersal traits

dispersal_df1 <- read.table("data/species_traits/Guillermo/Anne_dispersal_parameters_sp.csv",  sep = ",", header = TRUE)
dispersal_df2 <- read.table("data/species_traits/Guillermo/Anne_dispersal_parameters_sp2.csv", sep = ",", header = TRUE)

# subset data frames
dispersal_df1 <- subset(dispersal_df1, type == "average" & function_id == "Exponential", select = c(species, mean_parameter1, sd_parameter1))
dispersal_df2 <- subset(dispersal_df2, type == "average" & function_id == "Exponential", select = c(species, mean_parameter1, sd_parameter1))

# exclude robin / Erithacus rubecula
dispersal_df1 <- dispersal_df1[dispersal_df1$species!="Erithacus rubecula",]

# join data frames
dispersal_df <- rbind(dispersal_df1, dispersal_df2)

# harmonise species names
dispersal_df[dispersal_df$species=="Parus cristatus",'species'] <- "Lophophanes cristatus"
dispersal_df[dispersal_df$species=="Carduelis cannabina",'species'] <- "Linaria cannabina"

# transform from units of kilometre to metre
dispersal_df[,c('mean_parameter1','sd_parameter1')] <- 1000*dispersal_df[,c('mean_parameter1','sd_parameter1')]

# make rows for the species for which we have no dispersal data
missing_species <- species_traits$Species[!species_traits$Species %in% dispersal_df$species]
# assume mean and SD dispersal distance as the mean and SD of the mean dispersal distance of all given species
dispersal_df <- rbind(dispersal_df,	data.frame( species = missing_species,
												matrix(c(rep(mean(dispersal_df$mean_parameter1),length(missing_species)),
													     rep(sd(dispersal_df$mean_parameter1),length(missing_species))   ),
													   ncol = ncol(dispersal_df)-1,
													   dimnames = list(rows = NULL, cols = names(dispersal_df)[-1]))))
 

# Now, calculate RangeShifter parameters; emigration probability and mean dispersal distance
# use half cell size to account for random placement within cell 
cellsize <- 1000

dispersal_df$emg   <- pexp(cellsize/2, rate = 1/dispersal_df$mean_parameter1, lower.tail = FALSE )
dispersal_df$emg_u <- pexp(cellsize/2, rate = 1/(dispersal_df$mean_parameter1+dispersal_df$sd_parameter1), lower.tail = FALSE )
dispersal_df$emg_l <- pexp(cellsize/2, rate = 1/(dispersal_df$mean_parameter1-dispersal_df$sd_parameter1), lower.tail = FALSE )

dispersal_df$dpd   <- 1/dexp(cellsize/2, rate = 1/dispersal_df$mean_parameter1)
dispersal_df$dpd_u <- 1/dexp(cellsize/2, rate = 1/(dispersal_df$mean_parameter1+dispersal_df$sd_parameter1))
dispersal_df$dpd_l <- 1/dexp(cellsize/2, rate = 1/(dispersal_df$mean_parameter1-dispersal_df$sd_parameter1))




### 3.) Combine to one object and save

# merge data frames os demographic and dispersal traits
species_traits <- merge(species_traits, dispersal_df, by.x = "Species", by.y = "species")

# create names with points as delimiter because this is how they will be read in the calibration script
species_traits$species_name <- sapply(strsplit(species_traits$Species, " "), paste, collapse = ".")


# save to .Rdata object
save(species_traits, file = "data/species_traits/species_traits.Rdata")
