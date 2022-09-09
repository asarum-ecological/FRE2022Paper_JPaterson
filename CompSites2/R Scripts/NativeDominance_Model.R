## ---------------------------
##
## Script name:NativeDominance_Model.R 
##
## Purpose of script: Model proportion of plots covered in native vegetation
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
#note that models 2-4 use the same data table
MASTERDATA <- read.csv("CompSites2/FieldData/MODELS2-4/PLOTDATA_MASTER.csv") 

#Ensuring columns are in correct format
MASTERDATA$REFERENCE <- as.factor(MASTERDATA$REFERENCE)
MASTERDATA$SITE <- as.factor(MASTERDATA$SITE)
MASTERDATA$RC_Invasive <- as.numeric(MASTERDATA$RC_Invasive)

##Creating Subsets 

#old data for comp site-only analysis
# FRECOMPSITES <- MASTERDATA %>%
  # #Fraser River Subset (removing Pitt, Serpentine and Nicomekyl Rivers)
  # filter(RIVER == "Fraser",
  #        #Fraser Comp Site Subset (no natural/reference sites) 
  #        REFERENCE == "No") 

FREALLSITES <- MASTERDATA %>%
  #Fraser River Subset (removing Pitt, Serpentine and Nicomekyl Rivers)
  filter(RIVER == "Fraser")

# Scale numeric variables
#old code: comp sites only
# #fre_scale <- FRECOMPSITES %>%
#   mutate(sampling_age_scale = scale(SAMPLING_AGE),
#          km_upriver_scale = scale(KM_UPRIVER),
#          elev_adj_scale = scale(ELEV_ADJ),
#          prox_chan_scale = scale(PROX_CHAN),
#          sample_year_group = as.factor(SAMPLE_YEAR),
#          native_cover = RC_Native/100)

# Scale numeric variables
fre_scale <- FREALLSITES %>%
  mutate(km_upriver_scale = scale(KM_UPRIVER),
         elev_adj_scale = scale(ELEV_ADJ),
         prox_chan_scale = scale(PROX_CHAN),
         sample_year_group = as.factor(SAMPLE_YEAR),
         native_cover = RC_Native/100)

##### Fit models for native vegetation dominance ------------------------------

###RESEARCH QUESTION #2: Which factors influence the dominance of native species in created marshes?
#MODEL 2:Native Dominance
#elevation*distance upriver is under the assumption that elevation-related stresses are most pronounced at estuary mouth
#sample year was originally a random effect, but was moved to fixed effect because there were only 2 sample years

# Add two-stage
# 1 for difference between reference vs created (remove sampling age & INLAND)
# 2. test for an age effect within re

# GLMER with a binomial response (proportion)
natdom_glmer <- glmer(native_cover~INLAND + ARM + REFERENCE + km_upriver_scale*elev_adj_scale + prox_chan_scale + sample_year_group + (1|SITE),
                     data = fre_scale,
                     family = "binomial",
                     control=glmerControl(optimizer="bobyqa", 
                                          optCtrl = list(maxfun = 100000)))

# #old code for comp site only model
# natdom_glmer <- glmer(native_cover~INLAND + ARM + sampling_age_scale + km_upriver_scale*elev_adj_scale + prox_chan_scale + sample_year_group + (1|SITE),
#                       data = fre_scale,
#                       family = "binomial",
#                       control=glmerControl(optimizer="bobyqa", 
#                                            optCtrl = list(maxfun = 100000)))


# Model Diagnostics
par(mfrow = c(2, 2))
hist(resid(natdom_glmer),breaks = 50)
qqnorm(resid(natdom_glmer))
qqline(resid(natdom_glmer))
plot(resid(natdom_glmer)~fitted(natdom_glmer))
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(natdom_glmer), residuals(natdom_glmer)))
hist(fre_scale$native_cover, breaks = 50)


# Summary
summary(natdom_glmer)

# Main effect tests
car::Anova(natdom_glmer)
