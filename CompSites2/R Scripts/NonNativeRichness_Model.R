## ---------------------------
##
## Script name: NonNativeRichness_Model.R 
##
## Purpose of script: Model species richness of non-native species in vegetation plots
##
## Author: Daniel Stewart & James E Paterson
##
## Date Created:
##
## 
## Email:  daniel.stewart@asarum.org, j_paterson@ducks.ca
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

# Scale numeric predictor variables
fresites_scale <- FRESITES %>%
  mutate(sampling_age_scale = scale(SAMPLING_AGE),
         km_upriver_scale = scale(KM_UPRIVER),
         elev_adj_scale = scale(ELEV_ADJ),
         prox_chan_scale = scale(PROX_CHAN),
         sample_year_group = as.factor(SAMPLE_YEAR))


##### Model non-native species richness ---------------------------------------

###RESEARCH QUESTION #3: Which factors influence plant community diversity in created and natural tidal marshes?
#MODEL 3: Non-Native Richness
#elevation*distance upriver is under the assumption that elevation-related stresses are most pronounced at estuary mouth
#sample year was originally a random effect, but was moved to fixed effect because there were only 2 sample years


##### Fit model for species richness of non-native species in plots -----------

# Fit poisson glmer
nonatrich_glmer <- glmer(NN_RICH~INLAND + ARM + REFERENCE + prox_chan_scale + km_upriver_scale*elev_adj_scale + sample_year_group + (1|SITE),
                       family = poisson,
                       data = fresites_scale)

# Diagnostics
par(mfrow = c(2, 2))
hist(resid(nonatrich_glmer ),breaks = 50)
qqnorm(resid(nonatrich_glmer ))
qqline(resid(nonatrich_glmer ))
plot(resid(nonatrich_glmer )~fitted(nonatrich_glmer ))
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(nonatrich_glmer ), residuals(nonatrich_glmer )))
hist(fresites_scale$NN_RICH, breaks = 50)

# Summary
summary(nonatrich_glmer)

# Test main effects
car::Anova(nonatrich_glmer)

#Coefficient Plot
plot_model(nonatrich_glmer, show.values = TRUE, value.offset = .3, title = "Non-Native Richness", ci.lvl = .95,sort.est = TRUE,
           axis.lim = c(0.55,1.5),
           axis.labels = c('Sample Year','Elevation:Distance Upriver','Closed Embayment [Yes]','Reference Site [Yes]','Channel Proximity','Arm [North]','Elevation','Distance Upriver'))


