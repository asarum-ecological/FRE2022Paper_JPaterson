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





#######Plot generation for Figure 6########--------------------------------------------------------

#work in progress!!

#Predicting effect of distance upriver and elevation on native richness
#interaction of elevation an distance upriver
#first have to calculate mean, and mean +/- sd for visualisation
fre_scale$elev_adj_scale_2tile <- ntile(fre_scale$elev_adj_scale, 2)
fre_scale$elev_adj_scale_3tile <- ntile(fre_scale$elev_adj_scale, 3)
x <- fre_scale$elev_adj_scale

fre_scale$elev_adj_scale_group <-
  case_when(x > mean(x)+sd(x) ~ "high",
            x < mean(x)+sd(x) & x > mean(x)-sd(x) ~ "average",
            x < mean(x)-sd(x) ~ "low")

count(fre_scale,fre_scale$elev_adj_scale_group)
fre_scale$elev_adj_scale_group <- factor(fre_scale$elev_adj_scale_group, levels = c("high", "average", "low"))


predict_natrich_elevdist_glm <- data.frame(INLAND = "Yes",
                                       SITE = "02-001",
                                       REFERENCE = "No",
                                       prox_chan_scale = 0,
                                       elev_adj_scale = 0,
                                       sample_year_group = "2015",
                                       km_upriver_scale = seq(from = min(fre_scale$elev_adj_scale),
                                                             to = max(fre_scale$elev_adj_scale),
                                                             by = 0.1),
                                       ARM = "North") %>%
  mutate(predicted_glm = predict(natrich_glmer, newdata = ., type = "response"))

#Channel Prox/Nat Rich plot
Fig6a <- ggplot(data = predict_natrich_elevdist_glm , aes(x = km_upriver_scale, y = predicted_glm))+
  geom_point(data = fre_scale, aes(x = km_upriver_scale, y = NAT_RICH,group= elev_adj_scale_group, colour = elev_adj_scale_group),alpha = 0.3) +
  geom_smooth() +
  labs(x = "Distance Upriver (scaled)", y = "Native Richness/Plot") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#plotting interaction effects
visreg::visreg(natrich_glmer,"km_upriver_scale", by = "elev_adj_scale", overlay=TRUE,partial = FALSE, gg=TRUE) + 
  theme_bw()+
  xlab("Distance Upriver (km)") + ylab("Relative % Cover Native") +
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#Predicting effect of elevation on native dominance
predict_natrich_prox_glm <- data.frame(INLAND = "Yes",
                                      SITE = "02-001",
                                      REFERENCE = "No",
                                      elev_adj_scale = 0,
                                      km_upriver_scale = 0,
                                      sample_year_group = "2015",
                                      prox_chan_scale = seq(from = min(fre_scale$elev_adj_scale),
                                                           to = max(fre_scale$elev_adj_scale),
                                                           by = 0.1),
                                      ARM = "North") %>%
  mutate(predicted_glm = predict(natrich_glmer, newdata = ., type = "response"))

#Channel Prox/Nat Rich plot
Fig6b <- ggplot(data = predict_natrich_prox_glm, aes(x = prox_chan_scale, y = predicted_glm))+
  stat_smooth(col = "black") +
  geom_point(data = fre_scale, aes(x = prox_chan_scale, y = NAT_RICH),alpha = 0.3) +
  labs(x = "Channel Proximity (scaled)", y = "Native Richness/Plot") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

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
Fig6c <- ggplot(data = predict_natrich_ref_glm, aes(x = REFERENCE, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = REFERENCE, y = NAT_RICH)) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = NAT_RICH),alpha = 0.3) +
  labs(x = "Reference Site", y = "Native Richness/Plot") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#Predicting effect of embayment on native dominance
predict_natrich_emb_glm <- data.frame(REFERENCE = "No",
                                     SITE = "02-001",
                                     prox_chan_scale = 0,
                                     km_upriver_scale = 0,
                                     elev_adj_scale = 0,
                                     sample_year_group = "2015",
                                     INLAND = fre_scale$INLAND,
                                     ARM = "North") %>%
  mutate(predicted_glm = predict(natrich_glmer, newdata = ., type = "response"))

#Embayment/Dominance Plot
Fig6d <- ggplot(data = predict_natrich_emb_glm, aes(x = INLAND, y = predicted_glm))+
  geom_boxplot(data = fre_scale, aes(x = INLAND, y = NAT_RICH)) +
  geom_jitter(data = fre_scale, aes(x = REFERENCE, y = NAT_RICH),alpha = 0.3) +
  labs(x = "Closed Embayment", y = "") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11)) 

#formation of panel figure using cowplot
cowplot::plot_grid(Fig3a,Fig3b, Fig3c, Fig3d, ncol = 2)

#produce model summary table html that can be copied into MS
sjPlot::tab_model(natrich_glmer)




