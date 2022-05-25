#library("ggplot2")
library("tidyverse")
#library("car")
#library("visreg") #for model visualisation
#library("robustHD") #for model visualisation
#library("sjPlot") #for model visualisation
#library("cowplot") #for panel plots
#library("sjmisc") # for model output table
#library("sjlabelled") # for model output table

#LOADING MASTER DATA .CSV 
MASTERDATA <- read.csv("~/Documents/R/CompSites2/FieldData/MODEL1/SiteData_Master.csv") 

#Ensuring columns are in correct format
MASTERDATA$SAMPLE_YEAR <- as.factor(MASTERDATA$SAMPLE_YEAR)

#Transforming mudflat data to 0-1 scale for binomial model
MASTERDATA$MUDFLAT_BIN = (MASTERDATA$PRCNT_MUDFLAT/100)


##Creating Subsets 
#Fraser River Subset (removing Pitt, Serpentine and Nicomekyl Rivers)
FRESITES <- MASTERDATA %>%
  filter(RIVER == "Fraser") 
#Fraser Comp Site Subset (no natural/reference sites) 
FRECOMPSITES <- FRESITES %>%
  filter(REFERENCE == "NO") 

###RESEARCH QUESTION #1: What factors lead to marsh recession? 

#####MODEL#####
#MODEL 1: Percent Mudflat (recessed marsh)
#I did not use a mixed effects model for this model, with the rationale that the random effects used in our other models (site, sample year) are not relevant
#One could argue that "Year" could be used as a random effect, but much of the data used in these site-based models was acquired in 2021 using remote sensing

#Standard linear regression (used in DUC 2021 report)
#using a linear model with bound (proportional) data is problematic, also the skewed distribution (~half of sites had no recession)
RECMODEL1 <-lm(PRCNT_MUDFLAT ~ (INLAND  + SHEAR_BOOM + SLOUGH + OFFSHORE_STRUCTURE + AGE + AREA_MAPPED + DIST_UPRIVER_KM + ARM + PRCNT_EDGE2 + ELEV_ADJ), data = FRECOMPSITES)
summary(RECMODEL1)

#Potential alternative: quasibinomial model
#dependent variable had to be between 0-1, so we are using a transformed %mudflat column (MUDFLAT_BIN) - see line 17
RECMODEL2 <-glm(MUDFLAT_BIN ~ (INLAND  + SHEAR_BOOM + SLOUGH + OFFSHORE_STRUCTURE + AGE + AREA_MAPPED + DIST_UPRIVER_KM + ARM + PRCNT_EDGE2 + ELEV_ADJ), family = quasibinomial, data = FRECOMPSITES)
summary(RECMODEL2)
