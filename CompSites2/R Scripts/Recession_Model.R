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
## Email:  daniel.stewart@asarum.org, j_paterson@ducks.ca
##
## ---------------------------


##### Load libraries & data ---------------------------------------------------

library(ggplot2)
library(dplyr)
library(magrittr)
library(gamlss)
library(cowplot)

#LOADING MASTER DATA .CSV 
MASTERDATA <- read.csv("CompSites2/FieldData/MODEL1/SiteData_Master.csv")

# Data filters, scaling, updating formats
FRECOMPSITES <- MASTERDATA %>%
  filter(# Only include Fraser River sites, and remove Pitt, Serpentine and Nicomekyl Rivers
         RIVER == "Fraser",
         # Remove reference sites and only include created wetlands
         REFERENCE == "NO") %>%
  mutate(# Transform mudflat data to 0-1 scale for binomial model
         MUDFLAT_BIN = PRCNT_MUDFLAT/100,
         # Make sure year is a factor
         SAMPLE_YEAR = as.factor(SAMPLE_YEAR),
         # Add integer version of response & generic "erosion protection" to reduce number of predictor variables
         recessed = as.integer(PRCNT_MUDFLAT),
         erosion_protection = as.factor(ifelse(SHEAR_BOOM == "Present"|
                                                 #SLOUGH == "Yes"|
                                                 OFFSHORE_STRUCTURE == "Present",
                                               "y", "n")),
         # Scaling variables on very different scales
         # Note "as.numeric" added to remove extra attributes added by "scale()" that cause errors in predict
         sampling_age_scale = as.numeric(scale(SAMPLING_AGE)),
         km_upriver_scale = as.numeric(scale(DIST_UPRIVER_KM)),
         elev_adj_scale = as.numeric(scale(ELEV_ADJ)),
         percent_edge_scale = as.numeric(scale(PRCNT_EDGE2)),
         area_mapped_scale = as.numeric(scale(AREA_MAPPED)),
         age_scale = as.numeric(scale(AGE))
  )  %>%
  # Remove variables with NA values (to avoid error in gamlss, 
  # does not remove rows, only unused columns in this analyis)
  dplyr::select(!c(AREA_INTENDED, GRAZING, CARELYN_MH:RC_UnknownSD))

# Is response variable zero-inflated?
hist(FRECOMPSITES$recessed, breaks = 25)

###RESEARCH QUESTION #1: What factors lead to marsh recession? 

#####MODEL#####
#MODEL 1: Percent Mudflat (recessed marsh)
#I did not use a mixed effects model for this model, with the rationale that the random effects used in our other models (site, sample year) are not relevant
#One could argue that "Year" could be used as a random effect, but much of the data used in these site-based models was acquired in 2021 via imagery analysis


##### Fit GAMLSS Beta Zero Inflated Model for Marsh Recession -------------------------------------------

# Analyses using gamlss to fit a Beta zero inflated distribution model, which accounts for zero inflation
# but requires responses < 1.

# gamlss can fit Beta inflated distribution (BEINF), which is appropriate for inflated fixed-boundary responses (0,1) including boundary cases
# Fit recession response for mu (first equation) and nu. Mu affects cases >0 & <1, nu affects cases = 0, sigma affects precision
gamlss_bezi <- gamlss(MUDFLAT_BIN ~ erosion_protection + age_scale + km_upriver_scale*elev_adj_scale + ARM + percent_edge_scale + area_mapped_scale, 
             sigma.formula=~1,
             nu.formula=~erosion_protection + age_scale + km_upriver_scale*elev_adj_scale + ARM + percent_edge_scale + area_mapped_scale, #erosion_protection + AGE  + DIST_UPRIVER_KM + ARM + PRCNT_EDGE2 + ELEV_ADJ, 
             family = BEZI, 
             data = FRECOMPSITES %>%
               mutate(MUDFLAT_BIN = ifelse(MUDFLAT_BIN == 1,
                                            0.999,
                                            MUDFLAT_BIN))
             )

# More on the parameters in this model, see: http://www.gamlss.com/wp-content/uploads/2013/01/gamlss-manual.pdf

# Diagnostics
plot(gamlss_bezi)

# Summary
summary(gamlss_bezi)


##### Predicted effects on marsh recession ------------------------------------

# Three main effects of interest are:
# 1. Distance up-river,
# 2. Elevation
# 3. Project Age

# Predicting effect of distance upriver
predict_recession <- data.frame(erosion_protection = "n",
                               age_scale = 0,
                               km_upriver_scale = seq(from = min(FRECOMPSITES$km_upriver_scale),
                                                      to = max(FRECOMPSITES$km_upriver_scale),
                                                      by = 0.01),
                               elev_adj_scale = 0,
                               percent_edge_scale = mean(FRECOMPSITES$percent_edge_scale),
                               area_mapped_scale = mean(FRECOMPSITES$area_mapped_scale),
                               ARM = "North") %>%
  mutate(predicted_bezi = predict(gamlss_bezi, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         km_upriver =  (sd(FRECOMPSITES$DIST_UPRIVER_KM)*km_upriver_scale)+ mean(FRECOMPSITES$DIST_UPRIVER_KM))


# Distance upriver/recession Plot 
Plot2A <- ggplot(data = predict_recession, aes(x = km_upriver, y = predicted_bezi))+
  stat_smooth(col = "black", method = "loess") + 
  geom_point(data = FRECOMPSITES, aes(x = DIST_UPRIVER_KM, y = MUDFLAT_BIN)) +
  labs(x = "Distance upriver (km)", y = "Marsh recession (proportion)") + 
  annotate("text", x = 0, y = 1.0, label = "  (a)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 


# Predicting effect of elevation
predict_recession_elev = data.frame(erosion_protection = "n",
                               age_scale = 0,
                               km_upriver_scale = 0,
                               elev_adj_scale = seq(from = min(FRECOMPSITES$elev_adj_scale),
                                                    to = max(FRECOMPSITES$elev_adj_scale),
                                                    by = 0.01),
                               percent_edge_scale = mean(FRECOMPSITES$percent_edge_scale),
                               area_mapped_scale = mean(FRECOMPSITES$area_mapped_scale),
                               ARM = "North") %>%
  mutate(predicted_bezi = predict(gamlss_bezi, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         elev_adj =  (sd(FRECOMPSITES$ELEV_ADJ)*elev_adj_scale)+ mean(FRECOMPSITES$ELEV_ADJ))


# Elevation/recession plot
Plot2B <- ggplot(data = predict_recession_elev, aes(x = elev_adj, y = predicted_bezi))+
  stat_smooth(col = "black",method = "loess") +
  geom_point(data = FRECOMPSITES, aes(x = ELEV_ADJ, y = MUDFLAT_BIN)) +
  labs(x = "Elevation (m)", y ="") + 
  annotate("text", x = 0, y = 1.0, label = " (b)") +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

# Predicting effect of project age
predict_recession_age = data.frame(erosion_protection = "n",
                                    km_upriver_scale = 0,
                                    elev_adj_scale = 0,
                                    age_scale = seq(from = min(FRECOMPSITES$age_scale),
                                                   to = max(FRECOMPSITES$age_scale),
                                                   by = 0.1),
                                    percent_edge_scale = mean(FRECOMPSITES$percent_edge_scale),
                                    area_mapped_scale = mean(FRECOMPSITES$area_mapped_scale),
                                    ARM = "North") %>%
  mutate(predicted_bezi = predict(gamlss_bezi, newdata = ., type = "response"),
         age = (sd(FRECOMPSITES$AGE)*age_scale)+ mean(FRECOMPSITES$AGE))

# Age/recession plot
Plot2C <- ggplot(data = predict_recession_age, aes(x = age, y = predicted_bezi))+
  stat_smooth(col = "black",method = "loess") +
  geom_point(data = FRECOMPSITES, aes(x = AGE, y = MUDFLAT_BIN)) +
  labs(x = "Project age (years)", y = "") + 
  xlim(6,40) +
  annotate("text", x = 6, y = 1.0, label = "  (c)") +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#figure for paper (uses "cowplot" package)
cowplot::plot_grid(Plot2A,Plot2B,Plot2C, ncol = 3)

#model summary table html that can be copied into MS
sjPlot::tab_model(gamlss_bezi)

