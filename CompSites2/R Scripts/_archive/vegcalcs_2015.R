# Revision 2 of veg community analysis
# 29/June/2021

# libraries
library(tidyverse)
library(BiodiversityR)

# DATA IMPORT & PREP

# Path relative to working directory. Run getwd() to see your current working directory. It should print your path to "CompSites".
# If not, ensure you are working with the CompSites R project provided in the CompSites folder: "CompSites.Rproj".
# Note that "." here represents the current working directory.

veg <- read.csv("./FieldData/2015/csv/15-003-B.csv", fileEncoding = "UTF-8-BOM") # Modify filepath per site.
veg <- select(veg, 1:8) # workaround for some import issues of 2015 data


# Prep the data
veg$PERCENT_COVER <- as.numeric(veg$PERCENT_COVER) # ensure numeric cover data
veg <- subset(veg, COMMUNITY == 1 | is.na(COMMUNITY)) # select only community 1. Include "NA" for single community sites.

# transform data for richness calculations
veg.wide <- subset(veg, select = c(-COMMUNITY, -SITE_ID, -MAX_LH, -ORIGIN, -Scientific)) %>%
  spread(Common, PERCENT_COVER) # long to wide format
veg.wide[is.na(veg.wide)] <- 0 # empty cells are zero
veg.wide <- veg.wide[-1] # remove extra columns
veg.wide <- veg.wide %>% select(-any_of(c("WOOD", "MUD", "LITTER", "ROCK", "ALGAE",
                                          "mud", "mud ", "wood", "wood ", "litter", "litter ",
                                          "Ground", "ground", "Ground ", "ground ", "bare ground",
                                          "rock", "rock ", "rocks", "gravel", "sand", "sand ",
                                          "log", "log ", "Log", "Log ", "leaf litter",
                                          "Bare ground at stems", "Bare ground from above",
                                          "bare ground at stems", "bare ground from above"))) # exclude non-veg entries

# generate species lists
species <- subset(veg$Common, veg$ORIGIN == "N" | veg$ORIGIN == "E" | veg$ORIGIN == "I" | veg$ORIGIN == "U") %>% # unique species list
  unique() 
species.nat <- subset(veg$Common, veg$ORIGIN == "N") %>% # unique native species
  unique()
species.inv <- subset(veg$Common, veg$ORIGIN == "I") %>% # unique invasive species
  unique()
species.exo <- subset(veg$Common, veg$ORIGIN == "E") %>% # unique exotic species
  unique()
species.unk <- subset(veg$Common, veg$ORIGIN == "U") %>% # unique unknown species
  unique()


# native and invasive subsets
veg.wide.nat <- veg.wide %>% select(any_of(species.nat)) # select only native species
veg.wide.inv <- veg.wide %>% select(any_of(species.inv)) # select only invasive species
veg.wide.exo <- veg.wide %>% select(any_of(species.exo)) # select any non-native species
veg.wide.unk <- veg.wide %>% select(any_of(species.unk)) # select any non-native species


# CALCULATIONS

#Percent cover summary
PC_mean <- sapply(veg.wide, mean, na.rm = TRUE)
PC_sd <- sapply(veg.wide, sd, na.rm = TRUE)

# mean height of tallest Carex lyngbyei
lyngbyeHeight <- mean(veg$MAX_LH, na.rm = TRUE)
lyngbyesd <- sd(veg$MAX_LH, na.rm = TRUE)

# richness (native and total)
richness <- specnumber(veg.wide)
richness.nat <- specnumber(veg.wide.nat)

# shannon-weiner diversity index (native)
shannon <- diversity(veg.wide.nat, index = "shannon")

# simpson's diversity index (native)
simpson <- diversity(veg.wide.nat, index = "simpson")

# relative abundance
rel_ab <- function(origin, total) {
  originSums <- rowSums(origin)
  totalSums <- rowSums(total)
  return(mean(originSums/totalSums, na.rm = TRUE))
}

rel_ab_sd <- function(origin, total) {
  originSums <- rowSums(origin)
  totalSums <- rowSums(total)
  return(sd(originSums/totalSums, na.rm = TRUE))
}

natives <- rel_ab(veg.wide.nat, veg.wide)
nativesd <- rel_ab_sd(veg.wide.nat, veg.wide)

invasives <- rel_ab(veg.wide.inv, veg.wide)
invasivesd <- rel_ab_sd(veg.wide.inv, veg.wide)

exotics <- rel_ab(veg.wide.exo, veg.wide)
exoticsd <- rel_ab_sd(veg.wide.exo, veg.wide)

unknowns <- rel_ab(veg.wide.unk, veg.wide)
unknownsd <- rel_ab_sd(veg.wide.unk, veg.wide)


# RESULTS (modify filepath)

result <- data.frame(lyngbyeHeight,
                     lyngbyesd,
                     mean(richness),
                     mean(richness.nat),
                     mean(simpson),
                     sd(simpson),
                     mean(shannon),
                     sd(shannon),
                     natives,
                     nativesd,
                     exotics,
                     exoticsd,
                     invasives,
                     invasivesd,
                     unknowns,
                     unknownsd)

PC_result <- data.frame(PC_mean, PC_sd)

write.csv(result, "./Results/2015/15-003-B-results.csv") # veg analysis results
write.csv(species, "./Results/2015/15-003-B-species.csv") # unique species list
write.csv(PC_result, "./Results/2015/15-003-B-percentcover.csv") # summary of percent cover for each species


