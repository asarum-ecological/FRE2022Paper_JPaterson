library("ggplot2")
library("tidyverse")
library("lme4")


#LOADING MASTER DATA .CSV 
MASTERDATA <- read.csv("~/Documents/R/CompSites2/FieldData/MODELS2-4/PLOTDATA_MASTER.csv") 

#Ensuring columns are in correct format
MASTERDATA$REFERENCE <- as.factor(MASTERDATA$REFERENCE)
MASTERDATA$SITE <- as.factor(MASTERDATA$SITE)


##Creating Subsets 
#Fraser River Subset (removing Pitt, Serpentine and Nicomekyl Rivers)
FRESITES <- MASTERDATA %>%
  filter(RIVER == "Fraser") 

###RESEARCH QUESTION #3: Which factors influence plant community diversity in created and natural tidal marshes?
#MODEL 3: Native Richness
#elevation*distance upriver is under the assumption that elevation-related stresses are most pronounced at estuary mouth
#sample year was originally a random effect, but was moved to fixed effect because there were only 2 sample years

#Standard linear regression (used in DUC 2021 report)
#using a linear model with count(richness) data may be problematic
NATRIC_MODEL1 <- lmer(NAT_RICH~(INLAND + ARM + REFERENCE + PROX_CHAN + KM_UPRIVER*ELEV_ADJ + SAMPLE_YEAR) + (1|SITE),data = FRESITES)
summary(NATRIC_MODEL1)

#Potential alternative: poisson model
#variables likely need to be rescaled
NATRIC_MODEL2 <- glmer(NAT_RICH~(INLAND + ARM + REFERENCE + PROX_CHAN + KM_UPRIVER*ELEV_ADJ + SAMPLE_YEAR) + (1|SITE), family = poisson, data = FRESITES)
summary(NATRIC_MODEL2)
