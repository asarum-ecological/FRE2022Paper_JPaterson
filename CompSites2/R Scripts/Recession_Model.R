## ---------------------------
##
## Script name: Recession_Model.R 
##
## Purpose of script: Model marsh recession at FRE
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
library(gamlss)

#LOADING MASTER DATA .CSV 
MASTERDATA <- read.csv("CompSites2/FieldData/MODEL1/SiteData_Master.csv") 

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

##### Fit model for marsh recession -------------------------------------------

###RESEARCH QUESTION #1: What factors lead to marsh recession? 

#####MODEL#####
#MODEL 1: Percent Mudflat (recessed marsh)
#I did not use a mixed effects model for this model, with the rationale that the random effects used in our other models (site, sample year) are not relevant
#One could argue that "Year" could be used as a random effect, but much of the data used in these site-based models was acquired in 2021 using remote sensing

# Add integer version & generic "erosion protection"
FRECOMPSITES <- FRECOMPSITES %>%
  mutate(recessed = as.integer(PRCNT_MUDFLAT),
         erosion_protection = as.factor(ifelse(SHEAR_BOOM == "Present"|
                                       #SLOUGH == "Yes"|
                                       OFFSHORE_STRUCTURE == "Present",
                                     "y", "n")),
         sampling_age_scale = scale(SAMPLING_AGE),
         km_upriver_scale = scale(DIST_UPRIVER_KM),
         elev_adj_scale = scale(ELEV_ADJ)
         ) 

# Is response variable zero-inflated?
hist(FRECOMPSITES$recessed, breaks = 25)


# gamlss can fit Beta inflated distribution (BEINF), which is appropriate for inflated fixed-boundary responses (0,1) including boundary cases
m1 <- gamlss(MUDFLAT_BIN ~ erosion_protection + SLOUGH + AGE + km_upriver_scale*elev_adj_scale + ARM, 
             sigma.formula=~1, #erosion_protection + SLOUGH + AGE + km_upriver_scale*elev_adj_scale + ARM, 
             nu.formula=~erosion_protection + SLOUGH + AGE + km_upriver_scale*elev_adj_scale + ARM , #erosion_protection + AGE  + DIST_UPRIVER_KM + ARM + PRCNT_EDGE2 + ELEV_ADJ, 
             tau.formula=~1, 
             family = BEINF, 
             data = na.omit(FRECOMPSITES))
# Need to better understand the parameters in this model, see: http://www.gamlss.com/wp-content/uploads/2013/01/gamlss-manual.pdf


# Diagnostics
plot(m1)

# Summary
summary(m1)


##### Predicted effects on marsh recession ------------------------------------

# Two main effects are:
# 1. Distance up-river,
# 2. Elevation

# Predicting effect of distance upriver
predict_recession = data.frame(erosion_protection = "n",
                               SLOUGH = "No",
                               AGE = 0,
                               km_upriver_scale = seq(from = min(FRECOMPSITES$km_upriver_scale),
                                                      to = max(FRECOMPSITES$km_upriver_scale),
                                                      by = 0.1),
                               elev_adj_scale = 0,
                               ARM = "North") %>%
  mutate(predicted = predict(m1, newdata = ., type = "response"))

# Plot
ggplot(data = predict_recession, aes(x = km_upriver_scale, y = predicted))+
  stat_smooth(col = "black") +
  geom_point(data = FRECOMPSITES, aes(x = km_upriver_scale, y = MUDFLAT_BIN)) +
  labs(x = "Distance upriver (scaled)", y = "Marsh recession (proportion)") + 
  theme_classic()

# Predicting effect of distance upriver
predict_recession_elev = data.frame(erosion_protection = "n",
                               SLOUGH = "No",
                               AGE = 0,
                               km_upriver_scale = 0,
                               elev_adj_scale = seq(from = min(FRECOMPSITES$elev_adj_scale),
                                                    to = max(FRECOMPSITES$elev_adj_scale),
                                                    by = 0.1),
                               ARM = "North") %>%
  mutate(predicted = predict(m1, newdata = ., type = "response"))

# Plot
ggplot(data = predict_recession_elev, aes(x = elev_adj_scale, y = predicted))+
  stat_smooth(col = "black") +
  geom_point(data = FRECOMPSITES, aes(x = elev_adj_scale, y = MUDFLAT_BIN)) +
  labs(x = "Elevation (scaled)", y = "Marsh recession (proportion)") + 
  theme_classic()
