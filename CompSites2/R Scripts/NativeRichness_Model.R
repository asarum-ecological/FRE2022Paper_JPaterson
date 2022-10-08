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
MuMIn::r.squaredGLMM(natrich_glmer)

# Summary
summary(natrich_glmer)

# Test main effects
car::Anova(natrich_glmer)

#Coefficient Plot
set_theme(base = theme_classic()) #To remove the background color and the grids
Figure5a <- plot_model(natrich_glmer, show.values = TRUE, value.offset = .3, title = "Native Richness", ci.lvl = .95,sort.est = TRUE,
           axis.lim = c(0.55,1.5),
           axis.labels = c('Closed Embayment [Yes]','Sample Year','Elevation:Distance Upriver','Arm [North]','Reference Site [Yes]','Channel Proximity','Elevation','Distance Upriver'))





#######Plot generation for Figure 5########--------------------------------------------------------
# Three main effects of interest are:
# 1. Distance up-river,
# 2. Elevation
# 3. Reference Site
# 4. Closed Embayment 

#1. Predicting effect of distance upriver on native dominance--------------------
predict_natrich_distupr_glm <- data.frame(INLAND = "Yes",
                                         SITE = "02-001",
                                         REFERENCE = "No",
                                         prox_chan_scale = 0,
                                         elev_adj_scale = 0,
                                         sample_year_group = "2015",
                                         km_upriver_scale = seq(from = min(fre_scale$km_upriver_scale),
                                                                to = max(fre_scale$km_upriver_scale),
                                                                by = 0.1),
                                         ARM = "North") %>%
  mutate(predicted_glm = predict(natrich_glmer, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         DISTUPR_ADJ =  (sd(fre_scale$KM_UPRIVER)*km_upriver_scale)+ mean(fre_scale$KM_UPRIVER))


#Distance Upriver/Native Dominance plot
Fig5a <- ggplot(data = predict_natrich_distupr_glm, aes(x = DISTUPR_ADJ, y = predicted_glm))+
  stat_smooth(col = "black",method = "loess") +
  geom_point(data = fre_scale, aes(x = KM_UPRIVER, y = NAT_RICH),alpha = 0.2) +
  labs(x = "Distance Upriver (km)", y = "Native Richness") + 
  annotate("text", x = 0, y = 13.0, label = "  (a)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#2. Predicting effect of elevation on native dominance---------------------------
predict_natrich_elev_glm <- data.frame(INLAND = "Yes",
                                      SITE = "02-001",
                                      REFERENCE = "No",
                                      prox_chan_scale = 0,
                                      km_upriver_scale = 0,
                                      sample_year_group = "2015",
                                      elev_adj_scale = seq(from = min(fre_scale$elev_adj_scale),
                                                           to = max(fre_scale$elev_adj_scale),
                                                           by = 0.1),
                                      ARM = "North") %>%
  mutate(predicted_glm = predict(natrich_glmer, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         ELEV_ADJ =  (sd(fre_scale$ELEV_ADJ)*elev_adj_scale)+ mean(fre_scale$ELEV_ADJ))


#Elevation/Native Dominance plot
Fig5b <- ggplot(data = predict_natrich_elev_glm, aes(x = ELEV_ADJ, y = predicted_glm))+
  stat_smooth(col = "black",method = "loess") +
  geom_point(data = fre_scale, aes(x = ELEV_ADJ, y = NAT_RICH),alpha = 0.2) +
  labs(x = "Elevation (m)", y = "") + 
  annotate("text", x = -1, y = 13, label = "  (b)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#3. Effect of reference site on native dominance---------------------------
#Predicting effect of reference on native dominance
predict_natrich_ref_glm <- data.frame(INLAND = "Yes",
                                     SITE = "02-001",
                                     prox_chan_scale = 0,
                                     km_upriver_scale = 0,
                                     elev_adj_scale = 0,
                                     sample_year_group = "2015",
                                     REFERENCE = fre_scale$REFERENCE,
                                     ARM = "North") %>%
  mutate(predicted_glm = predict(natrich_glmer, newdata = ., type = "response"))

#Reference/Dominance Plot
Fig5c <- ggplot(data = predict_natrich_ref_glm, aes(x = REFERENCE, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = REFERENCE, y = NAT_RICH)) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = NAT_RICH),alpha = 0.2) +
  labs(x = "Reference Site", y = "") + 
  annotate("text", x = .5, y = 13, label = "  (c)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#4. Effect of embayment on native dominance---------------------------

#Predicting effect of embayment on native dominance
predict_natrich_ref_glm <- data.frame(REFERENCE = "No",
                                     SITE = "02-001",
                                     prox_chan_scale = 0,
                                     km_upriver_scale = 0,
                                     elev_adj_scale = 0,
                                     sample_year_group = "2015",
                                     INLAND = fre_scale$INLAND,
                                     ARM = "North") %>%
  mutate(predicted_glm = predict(natrich_glmer, newdata = ., type = "response"))

#Embayment/Dominance Plot
Fig5d <- ggplot(data = predict_natrich_ref_glm, aes(x = INLAND, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = INLAND, y = NAT_RICH)) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = NAT_RICH),alpha = 0.2) +
  labs(x = "Closed Embayment", y = "") + 
  annotate("text", x = .5, y = 13, label = "  (d)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 



#formation of panel figure using cowplot
NATPANEL <- cowplot::plot_grid(Fig5a,Fig5b, Fig5c, Fig5d, ncol = 4, nrow =1)


#produce model summary table html that can be copied into MS
sjPlot::tab_model(natrich_glmer)




