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
library(sjPlot)

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
MuMIn::r.squaredGLMM(natdom_glmer)

# Summary
summary(natdom_glmer)


# Main effect tests
car::Anova(natdom_glmer)


#Coefficient Plot
set_theme(base = theme_classic()) #To remove the background color and the grids
Figure4a = plot_model(natdom_glmer, show.values = TRUE, value.offset = .3, title = "Native Dominance", ci.lvl = .95,sort.est = TRUE,
           axis.lim = c(0.35,2),
           axis.labels = c('Embayed [Yes]','Dist. Upr.','Chan. Prox.','Ref. Site [Yes]','Elev.','Elev.:Dist. Upr.','Year','Arm [North]'))

###### Predicted effects on native dominance  ------------------------------------

# Three main effects of interest are:
# 1. Distance up-river,
# 2. Elevation
# 3. Reference Site
# 4. Closed Embayment 


#1. Predicting effect of distance upriver on native dominance--------------------
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
  mutate(predicted_glm = predict(natdom_glmer, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         DISTUPR_ADJ =  (sd(fre_scale$KM_UPRIVER)*km_upriver_scale)+ mean(fre_scale$KM_UPRIVER))


#Distance Upriver/Native Dominance plot
Fig3a <- ggplot(data = predict_natdom_distupr_glm, aes(x = DISTUPR_ADJ, y = predicted_glm))+
  stat_smooth(col = "black",method = "loess") +
  geom_point(data = fre_scale, aes(x = KM_UPRIVER, y = native_cover),alpha = 0.2) +
  labs(x = "Distance Upriver (km)", y = "Native Dominance (proportion)") + 
  annotate("text", x = -4, y = 1.0, label = "  (a)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#2. Predicting effect of elevation on native dominance---------------------------
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
  mutate(predicted_glm = predict(natdom_glmer, newdata = ., type = "response"),
# Add unscaled (original) variable by multiplying scaled value by STD and adding mean
ELEV_ADJ =  (sd(fre_scale$ELEV_ADJ)*elev_adj_scale)+ mean(fre_scale$ELEV_ADJ))


#Elevation/Native Dominance plot
Fig3b <- ggplot(data = predict_natdom_elev_glm, aes(x = ELEV_ADJ, y = predicted_glm))+
  stat_smooth(col = "black",method = "loess") +
  geom_point(data = fre_scale, aes(x = ELEV_ADJ, y = native_cover),alpha = 0.2) +
  labs(x = "Elevation (m)", y = "") + 
  annotate("text", x = -1.15, y = 1.0, label = "  (b)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#3. Effect of reference site on native dominance---------------------------
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
Fig3c <- ggplot(data = predict_natdom_ref_glm, aes(x = REFERENCE, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = REFERENCE, y = native_cover)) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = native_cover),alpha = 0.2) +
  labs(x = "Reference Site", y = "Native Dominance (proportion)") + 
  annotate("text", x = .5, y = 1.0, label = "  (c)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#3. Effect of reference site on native dominance---------------------------
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
Fig3c <- ggplot(data = predict_natdom_ref_glm, aes(x = REFERENCE, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = REFERENCE, y = native_cover)) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = native_cover),alpha = 0.2,outlier.shape = NA) +
  labs(x = "Reference Site", y = "Native Dominance (proportion)") + 
  annotate("text", x = .5, y = 1.0, label = "(c)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#4. Effect of embayment on native dominance---------------------------

#Predicting effect of embayment on native dominance
predict_natdom_ref_glm <- data.frame(REFERENCE = "No",
                                      SITE = "02-001",
                                      prox_chan_scale = 0,
                                      km_upriver_scale = 0,
                                      elev_adj_scale = 0,
                                      sample_year_group = "2015",
                                      INLAND = fre_scale$INLAND,
                                      ARM = "North") %>%
  mutate(predicted_glm = predict(natdom_glmer, newdata = ., type = "response"))

#Embayment/Dominance Plot
Fig3d <- ggplot(data = predict_natdom_ref_glm, aes(x = INLAND, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = INLAND, y = native_cover)) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = native_cover),alpha = 0.2,outlier.shape = NA) +
  labs(x = "Closed Embayment", y = "") + 
  annotate("text", x = .5, y = 1.0, label = "(d)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 


#formation of panel figure using cowplot
cowplot::plot_grid(Fig3a,Fig3b, Fig3c, Fig3d, ncol = 2)


#produce model summary table html that can be copied into MS
sjPlot::tab_model(natdom_glmer)


cowplot::plot_grid(Figure4a,Figure4b,Figure4c,ncol = 3)


###OLD CODE-----------------------------------------------------------------------
fre_scale$ELEVATION_2tile <- ntile(fre_scale$elev_adj_scale, 2)
fre_scale$ELEVATION_3tile <- ntile(fre_scale$elev_adj_scale, 3)
x <- fre_scale$elev_adj_scale

fre_scale$elev_adj_scale_group <-
  case_when(x > mean(x)+sd(x) ~ "high",
            x < mean(x)+sd(x) & x > mean(x)-sd(x) ~ "average",
            x < mean(x)-sd(x) ~ "low")

count(fre_scale,fre_scale$elev_adj_scale_group)
fre_scale$elev_adj_scale_group <- factor(fre_scale$elev_adj_scale_group, levels = c("high", "average", "low"))

HIGHSITES = subset(fre_scale, elev_adj_scale_group == "high")
MEDIUMSITES = subset(fre_scale, elev_adj_scale_group == "average")
LOWSITES = subset(fre_scale, elev_adj_scale_group == "low")
