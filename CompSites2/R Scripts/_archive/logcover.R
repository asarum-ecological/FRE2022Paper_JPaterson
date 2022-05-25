# Calculate percent cover of substrates using veg survey data
#
#

library(tidyverse)

# IMPORT AND CALCULATIONS

setwd("./FieldData/2021")
files <- list.files(pattern = ".csv")

debris <- vector()

for(i in 1:length(files)){
  site <- read.csv(files[i], fileEncoding = "UTF-8-BOM")
  site <- subset(site, site$COMMUNITY == 1)
  site$PERCENT_COVER <- as.numeric(site$PERCENT_COVER)
  
  siteWide <- subset(site, select = c(-COMMENTS, -COMMUNITY, -Site_Number, -MAX_LH_CM, -ORIGIN)) %>%
    spread(SPECIES_CODE, PERCENT_COVER)
  siteWide[is.na(siteWide)] <- 0

  
  if(any(names(siteWide) == "WOOD")){
    meanWood <- mean(siteWide$WOOD)
  } else {
    meanWood <- 0
  }
  
  if(any(names(siteWide) == "LOG")){
    meanLogs <- mean(siteWide$LOG)
  } else {
    meanLogs <- 0
  }
  
  meanDebris <- meanWood + meanLogs
  debris <- c(debris, meanDebris)
}


# RESULTS

result <- data.frame(SITE_ID = files,
                     DEBRIS = debris/100)

setwd("~/GitHub/CompSites")
write.csv(result, "./Results/2021/LogCover.csv")


