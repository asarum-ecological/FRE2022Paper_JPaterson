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
MuMIn::r.squaredGLMM(nonatrich_glmer)

# Summary
summary(nonatrich_glmer)

# Test main effects
car::Anova(nonatrich_glmer)

#Coefficient Plot
Figure4c = plot_model(nonatrich_glmer, show.values = TRUE, value.offset = .3, title = "Non-Native Richness", ci.lvl = .95,sort.est = TRUE,
           axis.lim = c(0.55,1.5),
           axis.labels = c('Sample Year','Elevation:Distance Upriver','Closed Embayment [Yes]','Reference Site [Yes]','Channel Proximity','Arm [North]','Elevation','Distance Upriver'))


###Plot Generation for Figure 6###----------------------------------------------------------------------

#work in progress!!

#Predicting effect of elevation on native dominance
predict_natdom_elev_glm <- data.frame(INLAND = "Yes",
                                      SITE = "02-001",
                                      REFERENCE = "No",
                                      prox_chan_scale = 0,
                                      km_upriver_scale = 0,
                                      sample_year_group = "2015",
                                      elev_adj_scale = seq(from = min(fre_scale$elev_adj_scale),
                                                           to = max(fre_scale$elev_adj_scale),
                                                           by = 0.1),
                                      ARM = "North") %>%
  mutate(predicted_glm = predict(natdom_glmer, newdata = ., type = "response"))

#Elevation/dominance plot
Fig6a <- ggplot(data = predict_natdom_elev_glm, aes(x = elev_adj_scale, y = predicted_glm))+
  stat_smooth(col = "black") +
  geom_point(data = fre_scale, aes(x = elev_adj_scale, y = native_cover),alpha = 0.3) +
  labs(x = "Elevation (scaled)", y = "Native Dominance (proportion)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#Predicting effect of distance upriver on native dominance
predict_natdom_distupr_glm <- data.frame(INLAND = "Yes",
                                         SITE = "02-001",
                                         REFERENCE = "No",
                                         prox_chan_scale = 0,
                                         elev_adj_scale = 0,
                                         sample_year_group = "2015",
                                         km_upriver_scale = seq(from = min(fre_scale$km_upriver_scale),
                                                                to = max(fre_scale$km_upriver_scale),
                                                                by = 0.1),
                                         ARM = "North") %>%
  mutate(predicted_glm = predict(natdom_glmer, newdata = ., type = "response"))

#Dist upriver/dominance plot
Fig6b <- ggplot(data = predict_natdom_distupr_glm, aes(x = km_upriver_scale, y = predicted_glm))+
  stat_smooth(col = "black") +
  geom_point(data = fre_scale, aes(x = km_upriver_scale, y = native_cover),alpha = 0.3) +
  labs(x = "Distance Upriver (scaled)", y = "") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#Predicting effect of reference on native dominance
predict_natdom_ref_glm <- data.frame(INLAND = "Yes",
                                     SITE = "02-001",
                                     prox_chan_scale = 0,
                                     km_upriver_scale = 0,
                                     elev_adj_scale = 0,
                                     sample_year_group = "2015",
                                     REFERENCE = fre_scale$REFERENCE,
                                     ARM = "North") %>%
  mutate(predicted_glm = predict(natdom_glmer, newdata = ., type = "response"))

#Reference/Dominance Plot
Fig6c <- ggplot(data = predict_natdom_ref_glm, aes(x = REFERENCE, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = REFERENCE, y = native_cover)) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = native_cover),alpha = 0.3) +
  labs(x = "Reference Site", y = "Native Dominance (proportion)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#Predicting effect of embayment on native dominance
predict_natdom_emb_glm <- data.frame(REFERENCE = "No",
                                     SITE = "02-001",
                                     prox_chan_scale = 0,
                                     km_upriver_scale = 0,
                                     elev_adj_scale = 0,
                                     sample_year_group = "2015",
                                     INLAND = fre_scale$INLAND,
                                     ARM = "North") %>%
  mutate(predicted_glm = predict(natdom_glmer, newdata = ., type = "response"))

#Embayment/Dominance Plot
Fig6d <- ggplot(data = predict_natdom_emb_glm, aes(x = INLAND, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = INLAND, y = native_cover)) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = native_cover),alpha = 0.3) +
  labs(x = "Closed Embayment", y = "") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#formation of panel figure using cowplot
cowplot::plot_grid(Fig3a,Fig3b, Fig3c, Fig3d, ncol = 2)


