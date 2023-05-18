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
           axis.labels = c('Year','Elev.:Dist. Upr.','Embayed [Yes]','Ref. Site [Yes]','Chan. Prox.','Arm [North]','Elev.','Dist. Upr.'))


#######Plot generation for Figure 5########--------------------------------------------------------
# Four main effects of interest are:
# 1. Distance up-river (with elevation as an interaction)
# 2. Elevation (with distance upriver as an interaction)
# 3. Reference Site
# 4. Closed Embayment 

#1. Effect of distance upriver on non-native richness with elevation as an interacting effect 

quantiles <- quantile(fre_scale$elev_adj_scale) #calculating quantiles
fre_scale$ELEVATION_2tile <- ntile(fre_scale$elev_adj_scale, 2)
fre_scale$ELEVATION_3tile <- ntile(fre_scale$elev_adj_scale, 3)
x <- fre_scale$elev_adj_scale

fre_scale$elev_adj_scale_group <- #categorizing by elevation class (not really necessary, but used for legend)
  case_when(x >= 0.6388005 ~ "upper",
            x < 0.6388005 & x > -0.6091675 ~ "mid",
            x <= -0.6091675 ~ "lower")
fre_scale$elev_adj_scale_group <- factor(fre_scale$elev_adj_scale_group, levels = c("lower","mid","upper")) #reorder categories for legend

#upper Elevation: Predicting effect of distance upriver on non-native richness, classified by elevation class
predict_nonatrich_distupr_glm_h <- data.frame(INLAND = "Yes",
                                              SITE = "02-001",
                                              REFERENCE = "No",
                                              prox_chan_scale = 0,
                                              elev_adj_scale = 0.6388005, #based on quantile calculation (line 189)
                                              sample_year_group = "2015",
                                              km_upriver_scale = seq(from = min(fre_scale$km_upriver_scale),
                                                                     to = max(fre_scale$km_upriver_scale),
                                                                     by = 0.1),
                                              ARM = "North") %>%
  mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         DISTUPR_ADJ =  (sd(fre_scale$KM_UPRIVER)*km_upriver_scale)+ mean(fre_scale$KM_UPRIVER))

#mid Elevation: Predicting effect of distance upriver on non-native richness, classified by elevation class
predict_nonatrich_distupr_glm_m <- data.frame(INLAND = "Yes",
                                              SITE = "02-001",
                                              REFERENCE = "No",
                                              prox_chan_scale = 0,
                                              elev_adj_scale = 0.1039572,
                                              sample_year_group = "2015",
                                              km_upriver_scale = seq(from = min(fre_scale$km_upriver_scale),
                                                                     to = max(fre_scale$km_upriver_scale),
                                                                     by = 0.1),
                                              ARM = "North") %>%
  mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         DISTUPR_ADJ =  (sd(fre_scale$KM_UPRIVER)*km_upriver_scale)+ mean(fre_scale$KM_UPRIVER))

predict_nonatrich_distupr_glm_l <- data.frame(INLAND = "Yes",
                                              SITE = "02-001",
                                              REFERENCE = "No",
                                              prox_chan_scale = 0,
                                              elev_adj_scale = -0.6091675,
                                              sample_year_group = "2015",
                                              km_upriver_scale = seq(from = min(fre_scale$km_upriver_scale),
                                                                     to = max(fre_scale$km_upriver_scale),
                                                                     by = 0.1),
                                              ARM = "North") %>%
  mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         DISTUPR_ADJ =  (sd(fre_scale$KM_UPRIVER)*km_upriver_scale)+ mean(fre_scale$KM_UPRIVER))


Fig5k <- ggplot(data = predict_nonatrich_distupr_glm_l, aes(x = DISTUPR_ADJ, y = predicted_glm)) +
  stat_smooth(col = "#5e3c99",method = "loess", data = predict_nonatrich_distupr_glm_l,size=.6) +
  stat_smooth(col = "#fdb863", method = "loess", data = predict_nonatrich_distupr_glm_m,size=.6) +
  stat_smooth(col = "#e66101",method = "loess", data = predict_nonatrich_distupr_glm_h,size=.6) +
  geom_point(data = fre_scale, aes(x = KM_UPRIVER, y = NN_RICH),alpha = 0.09) +
  labs(x = "Distance Upriver (km)", y = "Non-Native Richness/plot") + 
  ylim(0,13) +
  annotate("text", x = -1.15, y = 13.0, label = "  (e)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11))


#2. Effect of elevation on non-native richness, with distance upriver as an interacting effect

quantiles <- quantile(fre_scale$km_upriver_scale) #calculating quantiles
fre_scale$km_upriver_2tile <- ntile(fre_scale$km_upriver_scale, 2)
fre_scale$km_upriver_3tile <- ntile(fre_scale$km_upriver_scale, 3)
x <- fre_scale$km_upriver_scale
fre_scale$km_upriver_scale_group <- #quantile groups for legend 
  case_when(x >= 0.7214262 ~ "upper",
            x < 0.7214262 & x > -0.7559121 ~ "mid",
            x <= -0.7559121 ~ "lower")
fre_scale$km_upriver_scale_group <- factor(fre_scale$km_upriver_scale_group, levels = c("lower","mid","upper")) #determining order of legend

#upper Distance Upriver: Predicting effect of elevation on non-native richness, classified by distance upriver class
predict_nonatrich_elev_glm_h <- data.frame(INLAND = "Yes",
                                           SITE = "02-001",
                                           REFERENCE = "No",
                                           prox_chan_scale = 0,
                                           km_upriver_scale = 0.7214262, #based on quantile calculation (line 189)
                                           sample_year_group = "2015",
                                           elev_adj_scale = seq(from = min(fre_scale$elev_adj_scale),
                                                                to = max(fre_scale$elev_adj_scale),
                                                                by = 0.1),
                                           ARM = "North") %>%
  mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         ELEV_ADJ =  (sd(fre_scale$ELEV_ADJ)*elev_adj_scale)+ mean(fre_scale$ELEV_ADJ))


#mid Elevation: Predicting effect of elevation on non-native richness, classified by distance upriver class
predict_nonatrich_elev_glm_m <- data.frame(INLAND = "Yes",
                                           SITE = "02-001",
                                           REFERENCE = "No",
                                           prox_chan_scale = 0,
                                           km_upriver_scale = -0.2172876,
                                           sample_year_group = "2015",
                                           elev_adj_scale = seq(from = min(fre_scale$elev_adj_scale),
                                                                to = max(fre_scale$elev_adj_scale),
                                                                by = 0.1),
                                           ARM = "North") %>%
  mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         ELEV_ADJ =  (sd(fre_scale$ELEV_ADJ)*elev_adj_scale)+ mean(fre_scale$ELEV_ADJ))

#lower Elevation: Predicting effect of elevation on non-native richness, classified by distance upriver class
predict_nonatrich_elev_glm_l <- data.frame(INLAND = "Yes",
                                           SITE = "02-001",
                                           REFERENCE = "No",
                                           prox_chan_scale = 0,
                                           km_upriver_scale = -0.7559121,
                                           sample_year_group = "2015",
                                           elev_adj_scale = seq(from = min(fre_scale$km_upriver_scale),
                                                                to = max(fre_scale$km_upriver_scale),
                                                                by = 0.1),
                                           ARM = "North") %>%
  mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"),
         # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
         ELEV_ADJ =  (sd(fre_scale$ELEV_ADJ)*elev_adj_scale)+ mean(fre_scale$ELEV_ADJ))


Fig5l <- ggplot(data = predict_nonatrich_elev_glm_l, aes(x = ELEV_ADJ, y = predicted_glm)) +
  stat_smooth(col = "#5e3c99",method = "loess", data = predict_nonatrich_elev_glm_l,size=.6) +
  stat_smooth(col = "#fdb863", method = "loess", data = predict_nonatrich_elev_glm_m,size=.6) +
  stat_smooth(col = "#e66101",method = "loess", data = predict_nonatrich_elev_glm_h,size=.6) +
  geom_point(data = fre_scale, aes(x = ELEV_ADJ, y = NN_RICH),alpha = 0.09) +
  labs(x = "Elevation (m)", y = "") + 
  annotate("text", x = -1.15, y = 13.0, label = "  (f)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11))


#3. Effect of reference site on native dominance---------------------------
#Predicting effect of reference on native dominance
predict_nonatrich_ref_glm <- data.frame(INLAND = "Yes",
                                      SITE = "02-001",
                                      prox_chan_scale = 0,
                                      km_upriver_scale = 0,
                                      elev_adj_scale = 0,
                                      sample_year_group = "2015",
                                      REFERENCE = fre_scale$REFERENCE,
                                      ARM = "North") %>%
  mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"))

#Reference/Dominance Plot
Fig5i <- ggplot(data = predict_nonatrich_ref_glm, aes(x = REFERENCE, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = REFERENCE, y = NN_RICH),outlier.shape = NA) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = NN_RICH),alpha = 0.09) +
  labs(x = "Reference Site", y = "") + 
  annotate("text", x = .5, y = 13, label = "(h)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#4. Effect of embayment on native dominance---------------------------

#Predicting effect of embayment on native dominance
predict_nonatrich_emb_glm <- data.frame(REFERENCE = "No",
                                      SITE = "02-001",
                                      prox_chan_scale = 0,
                                      km_upriver_scale = 0,
                                      elev_adj_scale = 0,
                                      sample_year_group = "2015",
                                      INLAND = fre_scale$INLAND,
                                      ARM = "North") %>%
  mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"))

#Embayment/Dominance Plot
Fig5j <- ggplot(data = predict_nonatrich_emb_glm, aes(x = INLAND, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = INLAND, y = NN_RICH),outlier.shape = NA) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = NN_RICH),alpha = 0.09) +
  labs(x = "Closed Embayment", y = "") + 
  annotate("text", x = .5, y = 13, label = "  (g)") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 



#formation of panel figure using cowplot
NONATPANEL <-cowplot::plot_grid(Fig5k,Fig5l,Fig5j,Fig5i, nrow = 1, ncol = 4)
LEGEND <- cowplot::plot_grid(NULL,legend,ncol = 5)

cowplot::plot_grid(NATPANEL,NONATPANEL,LEGEND,nrow = 3, rel_heights = c(1,1,.1))

#produce model summary table html that can be copied into MS
sjPlot::tab_model(nonatrich_glmer)

###OLD CODE-------------

# #1. Predicting effect of distance upriver on non-native richness--------------------
# predict_nonatrich_distupr_glm <- data.frame(INLAND = "Yes",
#                                             SITE = "02-001",
#                                             REFERENCE = "No",
#                                             prox_chan_scale = 0,
#                                             elev_adj_scale = 0,
#                                             sample_year_group = "2015",
#                                             km_upriver_scale = seq(from = min(fre_scale$km_upriver_scale),
#                                                                    to = max(fre_scale$km_upriver_scale),
#                                                                    by = 0.1),
#                                             ARM = "North") %>%
#   mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"),
#          # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
#          DISTUPR_ADJ =  (sd(fre_scale$KM_UPRIVER)*km_upriver_scale)+ mean(fre_scale$KM_UPRIVER))
# 
# 
# #Distance Upriver/Native Dominance plot
# Fig5g <- ggplot(data = predict_nonatrich_distupr_glm, aes(x = DISTUPR_ADJ, y = predicted_glm))+
#   stat_smooth(col = "black",method = "loess") +
#   geom_point(data = fre_scale, aes(x = KM_UPRIVER, y = NN_RICH),alpha = 0.09) +
#   labs(x = "Distance Upriver (km)", y = "Non-Native Richness") + 
#   annotate("text", x = 0, y = 13.0, label = "  (e)") + 
#   theme_classic() +
#   theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 
# 
# #2. Predicting effect of elevation on native dominance---------------------------
# predict_nonatrich_elev_glm <- data.frame(INLAND = "Yes",
#                                          SITE = "02-001",
#                                          REFERENCE = "No",
#                                          prox_chan_scale = 0,
#                                          km_upriver_scale = 0,
#                                          sample_year_group = "2015",
#                                          elev_adj_scale = seq(from = min(fre_scale$elev_adj_scale),
#                                                               to = max(fre_scale$elev_adj_scale),
#                                                               by = 0.1),
#                                          ARM = "North") %>%
#   mutate(predicted_glm = predict(nonatrich_glmer, newdata = ., type = "response"),
#          # Add unscaled (original) variable by multiplying scaled value by STD and adding mean
#          ELEV_ADJ =  (sd(fre_scale$ELEV_ADJ)*elev_adj_scale)+ mean(fre_scale$ELEV_ADJ))
# 
# 
# #Elevation/Native Dominance plot
# Fig5h <- ggplot(data = predict_nonatrich_elev_glm, aes(x = ELEV_ADJ, y = predicted_glm))+
#   stat_smooth(col = "black",method = "loess") +
#   geom_point(data = fre_scale, aes(x = ELEV_ADJ, y = NN_RICH),alpha = 0.09) +
#   labs(x = "Elevation (m)", y = "") + 
#   annotate("text", x = -1, y = 13, label = "  (f)") + 
#   theme_classic() +
#   theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 