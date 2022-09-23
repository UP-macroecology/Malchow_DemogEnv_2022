
###-----------------------------------------------------------------
##   PREPARE LIST OF SPECIES HABITAT PREFERENCES
###-----------------------------------------------------------------


# This script extracts species traits from the publication:
#
# Storchová, Lenka, and David Hořák. "Life‐history characteristics of European birds." Global Ecology and Biogeography 27.4 (2018): 400-406.
#
# for our target bird species. The collected traits describe habitat preferences and are used to generate habitatmaps from CORINE data in the script
# prep_habitatmaps.R


### 1.) Subsetting species traits dataframe from Storchová & Hořák (2018)

# read table of traits downloaded from Storchová & Hořák (2018) SI
traits_df <- read.table("data/habitatmaps/species_preferences/Storchova_&_Horak_2018/Life-history characteristics of European birds.txt", sep = "\t", header = TRUE)

head(traits_df)

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
(species_traits <- traits_df[traits_df$Species %in% species_names, traits_names])


### 2.) Focus on habitat traits

# select only habitat preferences, roughly order species by these preferences, and save as object 
habitat_traits <- species_traits[c(2,4,5,1,9,8,7,3,10,6),c(1,33,19:22,25,26,32)]

# save to Rdata object
save(habitat_traits, file = "data/habitatmaps/species_preferences/habitat_traits.Rdata")
