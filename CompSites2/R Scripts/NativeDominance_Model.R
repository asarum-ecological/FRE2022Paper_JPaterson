library("ggplot2")
library("tidyverse")
library("lme4")

#LOADING MASTER DATA .CSV 
#note that models 2-4 use the same data table
MASTERDATA <- read.csv("~/Documents/R/CompSites2/FieldData/MODELS2-4/PLOTDATA_MASTER.csv") 

#Ensuring columns are in correct format
MASTERDATA$REFERENCE <- as.factor(MASTERDATA$REFERENCE)
MASTERDATA$SITE <- as.factor(MASTERDATA$SITE)
MASTERDATA$RC_Invasive <- as.numeric(MASTERDATA$RC_Invasive)

##Creating Subsets 
#Fraser River Subset (removing Pitt, Serpentine and Nicomekyl Rivers)
FRESITES <- MASTERDATA %>%
  filter(RIVER == "Fraser") 
#Fraser Comp Site Subset (no natural/reference sites) 
FRECOMPSITES <- FRESITES %>%
  filter(REFERENCE == "No") 

###RESEARCH QUESTION #2: Which factors influence the dominance of native species in created marshes?
#MODEL 2:Native Dominance
#elevation*distance upriver is under the assumption that elevation-related stresses are most pronounced at estuary mouth
#sample year was originally a random effect, but was moved to fixed effect because there were only 2 sample years

#Standard linear regression (used in DUC 2021 report)
#using a linear model with bound (proportional) data  may be problematic
NATDOM_MODEL1 <- lmer(RC_Native~(INLAND + ARM + SAMPLING_AGE + KM_UPRIVER*ELEV_ADJ + PROX_CHAN + SAMPLE_YEAR) + (1|SITE),data = FRECOMPSITES)
summary(NATDOM_MODEL1)
