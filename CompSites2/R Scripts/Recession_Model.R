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

#Add integer version & generic "erosion protection"
FRECOMPSITES <- FRECOMPSITES %>%
  mutate(recessed = as.integer(PRCNT_MUDFLAT),
         erosion_protection = as.factor(ifelse(SHEAR_BOOM == "Present"|
                                                 #SLOUGH == "Yes"|
                                                 OFFSHORE_STRUCTURE == "Present",
                                               "y", "n")),
         # Scaling variables on very different scales
         # Note "as.numeric" added to remove extra attributes added by "scale()" that cause errors in predict.glm
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


##### Fit model for marsh recession with binomial glm -------------------------------------------

# One approach is to use binomial GLM with proportion of marsh that is recessed (open mudflat) as the response variable
recession_glm <- glm(MUDFLAT_BIN ~ erosion_protection + age_scale + km_upriver_scale*elev_adj_scale + ARM + percent_edge_scale + area_mapped_scale, 
             family = "binomial", 
             data = FRECOMPSITES)

# Need to better understand the parameters in this model, see: http://www.gamlss.com/wp-content/uploads/2013/01/gamlss-manual.pdf

# Diagnostics
par(mfrow = c(2,2))
plot(recession_glm)
par(mfrow = c(1,1))
# Fit not great; maybe caused by zero-inflation
# Is approximately okay, though

# Summary
summary(recession_glm)

# Main effects test
car::Anova(recession_glm, type = 3)

# Is the strong positive effect of distance upriver present if removing the few points far upriver?
recession_glm_reduced <- glm(MUDFLAT_BIN ~ erosion_protection + age_scale + km_upriver_scale*elev_adj_scale + ARM + percent_edge_scale + area_mapped_scale, 
                     family = "binomial", 
                     data = FRECOMPSITES %>%
                       filter(km_upriver_scale <= 1.96 &
                                km_upriver_scale >= -1.96))

# Diagnostics
par(mfrow = c(2,2))
plot(recession_glm_reduced)
par(mfrow = c(1,1))


# Summary
summary(recession_glm_reduced)

# Main effects test
car::Anova(recession_glm_reduced, type = 3)


##### Predicted effects on marsh recession with GLM ------------------------------------

# Two main effects are:
# 1. Distance up-river,
# 2. Elevation

# Predicting effect of distance upriver
predict_recession_glm <- data.frame(erosion_protection = "n",
                               age_scale = 0,
                               km_upriver_scale = seq(from = min(FRECOMPSITES$km_upriver_scale),
                                                      to = max(FRECOMPSITES$km_upriver_scale),
                                                      by = 0.1),
                               elev_adj_scale = 0,
                               percent_edge_scale = mean(FRECOMPSITES$percent_edge_scale),
                               area_mapped_scale = mean(FRECOMPSITES$area_mapped_scale),
                               ARM = "North") %>%
  mutate(predicted = predict(recession_glm, newdata = ., 
                             type = "response"),
         predicted_se = predict(recession_glm, newdata = ., 
                                type = "response", se = TRUE)$se.fit)

# Plot
ggplot(data = predict_recession_glm, aes(x = km_upriver_scale, y = predicted))+
  geom_ribbon(aes(ymin = predicted-predicted_se,
                  ymax = predicted + predicted_se),
              fill = "lightgray") +
  stat_smooth(col = "black") +
  geom_point(data = FRECOMPSITES, aes(x = km_upriver_scale, y = MUDFLAT_BIN)) +
  labs(x = "Distance upriver (scaled)", y = "Marsh recession (proportion)") + 
  theme_classic()


# Predicting effect of distance upriver with reduced model
predict_recession_glm_reduced <- data.frame(erosion_protection = "n",
                                    age_scale = 0,
                                    km_upriver_scale = seq(from = min(FRECOMPSITES$km_upriver_scale),
                                                           to = max(FRECOMPSITES$km_upriver_scale),
                                                           by = 0.1),
                                    elev_adj_scale = 0,
                                    percent_edge_scale = mean(FRECOMPSITES$percent_edge_scale),
                                    area_mapped_scale = mean(FRECOMPSITES$area_mapped_scale),
                                    ARM = "North") %>%
  mutate(predicted = predict(recession_glm_reduced, newdata = ., 
                             type = "response"),
         predicted_se = predict(recession_glm_reduced, newdata = ., 
                                type = "response", se = TRUE)$se.fit)

# Plot
ggplot(data = predict_recession_glm_reduced, aes(x = km_upriver_scale, y = predicted))+
  geom_ribbon(data = predict_recession_glm, aes(x = km_upriver_scale,
                                                ymin = predicted-predicted_se,
                                                ymax = predicted + predicted_se),
              fill = "darkgray") +
  geom_ribbon(aes(ymin = predicted-predicted_se,
                  ymax = predicted + predicted_se),
              fill = "lightgray") +
  stat_smooth(col = "black") +
  stat_smooth(data = predict_recession_glm, aes(x = km_upriver_scale, y = predicted),
              col = "red") +
  geom_point(data = FRECOMPSITES, aes(x = km_upriver_scale, y = MUDFLAT_BIN)) +
  labs(x = "Distance upriver (scaled)", y = "Marsh recession (proportion)") + 
  theme_classic()

# Predicting effect of distance upriver
predict_recession_elev_glm <- data.frame(erosion_protection = "n",
                                    age_scale = 0,
                                    km_upriver_scale = 0,
                                    elev_adj_scale = seq(from = min(FRECOMPSITES$elev_adj_scale),
                                                         to = max(FRECOMPSITES$elev_adj_scale),
                                                         by = 0.1),
                                    percent_edge_scale = mean(FRECOMPSITES$percent_edge_scale),
                                    area_mapped_scale = mean(FRECOMPSITES$area_mapped_scale),
                                    ARM = "North") %>%
  mutate(predicted = predict(recession_glm, newdata = ., type = "response"),
         predicted_se = predict(recession_glm, newdata = ., 
                                type = "response", se = TRUE)$se.fit)

# Plot
ggplot(data = predict_recession_elev_glm, aes(x = elev_adj_scale, y = predicted))+
  geom_ribbon(aes(ymin = predicted-predicted_se,
                  ymax = predicted + predicted_se),
              fill = "lightgray") +
  stat_smooth(col = "black") +
  geom_point(data = FRECOMPSITES, aes(x = elev_adj_scale, y = MUDFLAT_BIN)) +
  labs(x = "Elevation (scaled)", y = "Marsh recession (proportion)") + 
  theme_classic()

##### Fit model for marsh recession -------------------------------------------

# ALternate analyses using gamlss to fit a Beta-inflated distribution model, which accounts for zero inflation
# and accommodates fixed-boundary responses (e.g. [0,1] as in here).

# gamlss can fit Beta inflated distribution (BEINF), which is appropriate for inflated fixed-boundary responses (0,1) including boundary cases
# Fit same variable response for mu (first equation) and nu. Mu affects cases >0 p <1, nu affects po and p1
m1 <- gamlss(MUDFLAT_BIN ~ erosion_protection + age_scale + km_upriver_scale*elev_adj_scale + ARM + percent_edge_scale + area_mapped_scale, 
             sigma.formula=~1,
             nu.formula=~erosion_protection + age_scale + km_upriver_scale*elev_adj_scale + ARM + percent_edge_scale + area_mapped_scale, #erosion_protection + AGE  + DIST_UPRIVER_KM + ARM + PRCNT_EDGE2 + ELEV_ADJ, 
             tau.formula=~1, 
             family = BEINF, 
             data = FRECOMPSITES)

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
                               age_scale = 0,
                               km_upriver_scale = seq(from = min(FRECOMPSITES$km_upriver_scale),
                                                      to = max(FRECOMPSITES$km_upriver_scale),
                                                      by = 0.1),
                               elev_adj_scale = 0,
                               percent_edge_scale = mean(FRECOMPSITES$percent_edge_scale),
                               area_mapped_scale = mean(FRECOMPSITES$area_mapped_scale),
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
                               age_scale = 0,
                               km_upriver_scale = 0,
                               elev_adj_scale = seq(from = min(FRECOMPSITES$elev_adj_scale),
                                                    to = max(FRECOMPSITES$elev_adj_scale),
                                                    by = 0.1),
                               percent_edge_scale = mean(FRECOMPSITES$percent_edge_scale),
                               area_mapped_scale = mean(FRECOMPSITES$area_mapped_scale),
                               ARM = "North") %>%
  mutate(predicted = predict(m1, newdata = ., type = "response"))

# Plot
ggplot(data = predict_recession_elev, aes(x = elev_adj_scale, y = predicted))+
  stat_smooth(col = "black") +
  geom_point(data = FRECOMPSITES, aes(x = elev_adj_scale, y = MUDFLAT_BIN)) +
  labs(x = "Elevation (scaled)", y = "Marsh recession (proportion)") + 
  theme_classic()


# Trying reduced gamlss
FRECOMPSITES_reduced <- FRECOMPSITES %>%
  filter(km_upriver_scale <= 1.96)
m2 <- gamlss(MUDFLAT_BIN ~ erosion_protection + age_scale + km_upriver_scale*elev_adj_scale + ARM + percent_edge_scale + area_mapped_scale, 
             sigma.formula=~1,
             nu.formula=~erosion_protection + age_scale + km_upriver_scale*elev_adj_scale + ARM + percent_edge_scale + area_mapped_scale, #erosion_protection + AGE  + DIST_UPRIVER_KM + ARM + PRCNT_EDGE2 + ELEV_ADJ, 
             tau.formula=~1, 
             family = BEINF, 
             data = FRECOMPSITES_reduced)

# Diagnostics
plot(m2)

# Summary
summary(m2)

# Predicting effect of distance upriver
predict_recession_m2 = data.frame(erosion_protection = "n",
                               age_scale = 0,
                               km_upriver_scale = seq(from = min(FRECOMPSITES$km_upriver_scale),
                                                      to = max(FRECOMPSITES$km_upriver_scale),
                                                      by = 0.1),
                               elev_adj_scale = 0,
                               percent_edge_scale = mean(FRECOMPSITES$percent_edge_scale),
                               area_mapped_scale = mean(FRECOMPSITES$area_mapped_scale),
                               ARM = "North") %>%
  mutate(predicted = predict(m2, newdata = ., type = "response"))

# Plot
ggplot(data = predict_recession_m2, aes(x = km_upriver_scale, y = predicted))+
  stat_smooth(col = "black") +
  geom_point(data = FRECOMPSITES, aes(x = km_upriver_scale, y = MUDFLAT_BIN)) +
  labs(x = "Distance upriver (scaled)", y = "Marsh recession (proportion)") + 
  theme_classic()
