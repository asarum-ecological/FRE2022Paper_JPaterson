## ---------------------------
##
## Script name: NativeRichness_Model.R 
##
## Purpose of script: Model species richness of native species in vegetation plots
##
## Author: Daniel Stewart & James E Paterson
##
## Date Created:
##
## 
## Email: 
##
## ---------------------------

##### Load libraries & data ---------------------------------------------------

library(ggplot2)
library(dplyr)
library(magrittr)
library(lme4)


#LOADING MASTER DATA .CSV 
MASTERDATA <- read.csv("CompSites2/FieldData/MODELS2-4/PLOTDATA_MASTER.csv") 

#Ensuring columns are in correct format
MASTERDATA$REFERENCE <- as.factor(MASTERDATA$REFERENCE)
MASTERDATA$SITE <- as.factor(MASTERDATA$SITE)


##Creating Subsets 
#Fraser River Subset (removing Pitt, Serpentine and Nicomekyl Rivers)
FRESITES <- MASTERDATA %>%
  filter(RIVER == "Fraser") 

# Re-scale numeric predictor variables
fresites_scale <- FRESITES %>%
  mutate(sampling_age_scale = scale(SAMPLING_AGE),
         km_upriver_scale = scale(KM_UPRIVER),
         elev_adj_scale = scale(ELEV_ADJ),
         prox_chan_scale = scale(PROX_CHAN),
         sample_year_group = as.factor(SAMPLE_YEAR))


##### Fit GLMER explaining native species richness ----------------------------

###RESEARCH QUESTION #3: Which factors influence plant community diversity in created and natural tidal marshes?
#MODEL 3: Native Richness
#elevation*distance upriver is under the assumption that elevation-related stresses are most pronounced at estuary mouth
#sample year was originally a random effect, but was moved to fixed effect because there were only 2 sample years


# Fit poisson glmer
natrich_glmer <- glmer(NAT_RICH~INLAND + ARM + REFERENCE + prox_chan_scale + km_upriver_scale*elev_adj_scale + sample_year_group + (1|SITE),
                       family = poisson,
                       data = fresites_scale)

# Diagnostics
par(mfrow = c(2, 2))
hist(resid(natrich_glmer),breaks = 50)
qqnorm(resid(natrich_glmer))
qqline(resid(natrich_glmer))
plot(resid(natrich_glmer)~fitted(natrich_glmer))
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(natrich_glmer), residuals(natrich_glmer)))
hist(fresites_scale$NAT_RICH, breaks = 50)

# Summary
summary(natrich_glmer)

# Test main effects
car::Anova(natrich_glmer)
