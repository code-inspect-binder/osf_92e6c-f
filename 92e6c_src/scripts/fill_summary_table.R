# This file takes outputs from other source files and combines them into a single dataset for use in producing results/figures reported in Kraft et al.

# Packages and functions
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(gtable)
library(grid)
library(cowplot)
library(gmodels)

confidence_interval <- function(vector, interval) {
  
  vector <- vector[which(!is.na(vector))]
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}


# Load various data files
returns <- read_excel("ape_return_rates.xlsx")  # great ape data
all <- read.csv("summary_data.csv", stringsAsFactors = F)  # species/population values for mass/TEE/body fat/etc
comp <- read_excel("comparative_return_rates.xlsx", sheet=1) # cross cultural data
had <- read.csv("hadza_all.csv", stringsAsFactors = F)  # hadza data
tsimane <- read.csv("tsimane_all.csv", stringsAsFactors = F)  # tsimane data

# set age at which to extract values for humans. Age 40 chosen to reflect a near-peak productive adult.
p_age <- 40

# add columns to main dataframe for variables calculated from posteriors 
all$low.time_forage_and_moving <- NA
all$hi.time_forage_and_moving <- NA
all$Rg <- NA
all$Rg_low <- NA
all$Rg_hi <- NA
all$Rn <- NA
all$Rn_low <- NA
all$Rn_hi <- NA
all$EE <- NA
all$EE_low <- NA
all$EE_hi <- NA


## HADZA
all$prod[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$Ea[which(had$age == p_age & had$male == 1)]
all$prod[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$Ea[which(had$age == p_age & had$male == 0)]

all$hi.prod[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$Ea.hi[which(had$age == p_age & had$male == 1)]
all$hi.prod[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$Ea.hi[which(had$age == p_age & had$male == 0)]

all$low.prod[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$Ea.lo[which(had$age == p_age & had$male == 1)]
all$low.prod[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$Ea.lo[which(had$age == p_age & had$male == 0)]

all$cost[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$cost[which(had$age == p_age & had$male == 1)]
all$cost[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$cost[which(had$age == p_age & had$male == 0)]

all$hi.cost[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$hi.cost[which(had$age == p_age & had$male == 1)]
all$hi.cost[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$hi.cost[which(had$age == p_age & had$male == 0)]

all$low.cost[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$low.cost[which(had$age == p_age & had$male == 1)]
all$low.cost[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$low.cost[which(had$age == p_age & had$male == 0)]

all$time_forage_and_moving[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$hrs.worked[which(had$age == p_age & had$male == 1)]
all$time_forage_and_moving[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$hrs.worked[which(had$age == p_age & had$male == 0)]

all$hi.time_forage_and_moving[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$hi.hrs.worked[which(had$age == p_age & had$male == 1)]
all$hi.time_forage_and_moving[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$hi.hrs.worked[which(had$age == p_age & had$male == 0)]

all$low.time_forage_and_moving[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$low.hrs.worked[which(had$age == p_age & had$male == 1)]
all$low.time_forage_and_moving[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$low.hrs.worked[which(had$age == p_age & had$male == 0)]

all$Rg[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$Rg[which(had$age == p_age & had$male == 1)]
all$Rg[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$Rg[which(had$age == p_age & had$male == 0)]

all$Rg_low[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$Rg_low[which(had$age == p_age & had$male == 1)]
all$Rg_low[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$Rg_low[which(had$age == p_age & had$male == 0)]

all$Rg_hi[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$Rg_hi[which(had$age == p_age & had$male == 1)]
all$Rg_hi[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$Rg_hi[which(had$age == p_age & had$male == 0)]

all$Rn[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$Rn[which(had$age == p_age & had$male == 1)]
all$Rn[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$Rn[which(had$age == p_age & had$male == 0)]

all$Rn_low[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$Rn_low[which(had$age == p_age & had$male == 1)]
all$Rn_low[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$Rn_low[which(had$age == p_age & had$male == 0)]

all$Rn_hi[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$Rn_hi[which(had$age == p_age & had$male == 1)]
all$Rn_hi[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$Rn_hi[which(had$age == p_age & had$male == 0)]

all$EE[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$EE[which(had$age == p_age & had$male == 1)]
all$EE[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$EE[which(had$age == p_age & had$male == 0)]

all$EE_low[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$EE_low[which(had$age == p_age & had$male == 1)]
all$EE_low[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$EE_low[which(had$age == p_age & had$male == 0)]

all$EE_hi[which(all$Species == "Human (Hadza)" & all$Sex == "male")] <- had$EE_hi[which(had$age == p_age & had$male == 1)]
all$EE_hi[which(all$Species == "Human (Hadza)" & all$Sex == "female")] <- had$EE_hi[which(had$age == p_age & had$male == 0)]

## TSIMANE
all$prod[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$Ea[which(tsimane$age == p_age & tsimane$male == 1)]
all$prod[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$Ea[which(tsimane$age == p_age & tsimane$male == 0)]

all$low.prod[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$Ea.lo[which(tsimane$age == p_age & tsimane$male == 1)]
all$low.prod[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$Ea.lo[which(tsimane$age == p_age & tsimane$male == 0)]

all$hi.prod[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$Ea.hi[which(tsimane$age == p_age & tsimane$male == 1)]
all$hi.prod[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$Ea.hi[which(tsimane$age == p_age & tsimane$male == 0)]

all$time_forage_and_moving[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$hrs.worked[which(tsimane$age == p_age & tsimane$male == 1)]
all$time_forage_and_moving[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$hrs.worked[which(tsimane$age == p_age & tsimane$male == 0)]

all$low.time_forage_and_moving[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$low.hrs.worked[which(tsimane$age == p_age & tsimane$male == 1)]
all$low.time_forage_and_moving[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$low.hrs.worked[which(tsimane$age == p_age & tsimane$male == 0)]

all$hi.time_forage_and_moving[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$hi.hrs.worked[which(tsimane$age == p_age & tsimane$male == 1)]
all$hi.time_forage_and_moving[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$hi.hrs.worked[which(tsimane$age == p_age & tsimane$male == 0)]

all$cost[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$cost[which(tsimane$age == p_age & tsimane$male == 1)]
all$cost[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$cost[which(tsimane$age == p_age & tsimane$male == 0)]

all$low.cost[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$low.cost[which(tsimane$age == p_age & tsimane$male == 1)]
all$low.cost[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$low.cost[which(tsimane$age == p_age & tsimane$male == 0)]

all$hi.cost[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$hi.cost[which(tsimane$age == p_age & tsimane$male == 1)]
all$hi.cost[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$hi.cost[which(tsimane$age == p_age & tsimane$male == 0)]

all$Rg[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$Rg[which(tsimane$age == p_age & tsimane$male == 1)]
all$Rg[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$Rg[which(tsimane$age == p_age & tsimane$male == 0)]

all$Rg_low[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$Rg_low[which(tsimane$age == p_age & tsimane$male == 1)]
all$Rg_low[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$Rg_low[which(tsimane$age == p_age & tsimane$male == 0)]

all$Rg_hi[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$Rg_hi[which(tsimane$age == p_age & tsimane$male == 1)]
all$Rg_hi[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$Rg_hi[which(tsimane$age == p_age & tsimane$male == 0)]

all$Rn[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$Rn[which(tsimane$age == p_age & tsimane$male == 1)]
all$Rn[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$Rn[which(tsimane$age == p_age & tsimane$male == 0)]

all$Rn_low[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$Rn_low[which(tsimane$age == p_age & tsimane$male == 1)]
all$Rn_low[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$Rn_low[which(tsimane$age == p_age & tsimane$male == 0)]

all$Rn_hi[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$Rn_hi[which(tsimane$age == p_age & tsimane$male == 1)]
all$Rn_hi[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$Rn_hi[which(tsimane$age == p_age & tsimane$male == 0)]

all$EE[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$EE[which(tsimane$age == p_age & tsimane$male == 1)]
all$EE[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$EE[which(tsimane$age == p_age & tsimane$male == 0)]

all$EE_low[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$EE_low[which(tsimane$age == p_age & tsimane$male == 1)]
all$EE_low[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$EE_low[which(tsimane$age == p_age & tsimane$male == 0)]

all$EE_hi[which(all$Species == "Human (Tsimane)" & all$Sex == "male")] <- tsimane$EE_hi[which(tsimane$age == p_age & tsimane$male == 1)]
all$EE_hi[which(all$Species == "Human (Tsimane)" & all$Sex == "female")] <- tsimane$EE_hi[which(tsimane$age == p_age & tsimane$male == 0)]


## Non-human great apes
# Errors combine in quadrature when we are combining estimates of confidence intervals, see here for details (http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf)
apes <- returns %>% 
  mutate(time_forage_and_moving = time_spent_feeding + time_spent_moving) %>%
  group_by(Species, Sex) %>%
  summarise(low.time_forage_and_moving = confidence_interval(time_forage_and_moving, 0.95)[1]/60,
            hi.time_forage_and_moving = confidence_interval(time_forage_and_moving, 0.95)[2]/60,
            time_forage = mean(time_spent_feeding, na.rm=T)/60, 
            time_forage_and_moving = mean(time_forage_and_moving, na.rm=T)/60) %>%
  left_join(select(all, Species, Sex, cost, hi.cost, low.cost, prod, low.prod, hi.prod)) %>% 
  mutate(Rg = prod/time_forage_and_moving, 
         Rn = (prod-cost)/time_forage_and_moving, 
                      EE = prod/cost,
         Rg_low = Rg - (Rg*sqrt((abs(prod-low.prod)/prod)^2 + (abs(time_forage_and_moving-low.time_forage_and_moving)/time_forage_and_moving)^2)),
         Rg_hi = Rg + (Rg*sqrt((abs(prod-hi.prod)/prod)^2 + (abs(time_forage_and_moving-hi.time_forage_and_moving)/time_forage_and_moving)^2)),
         Rn_low = Rn - (Rn*sqrt((abs(prod-low.prod)/prod)^2 + (abs(time_forage_and_moving-low.time_forage_and_moving)/time_forage_and_moving)^2)),
         Rn_hi = Rn + (Rn*sqrt((abs(prod-low.prod)/prod)^2 + (abs(time_forage_and_moving-low.time_forage_and_moving)/time_forage_and_moving)^2)),
         #for apes that lack EE cost error, use only time error
          EE_low = EE - (EE*sqrt((abs(prod-low.prod)/prod)^2 + (abs(cost-low.cost)/cost)^2)),
          EE_hi = EE + (EE*sqrt((abs(prod-hi.prod)/prod)^2 + (abs(cost-hi.cost)/cost)^2))
         )


# populate great ape cells in main dataframe
for(i in which(all$Species %in% apes$Species)){
  
  NA.cells <- which(is.na(all[i,]) & names(all) %in% names(apes))
  
  all[i, NA.cells] <- apes[which(apes$Species == all$Species[i] & apes$Sex == all$Sex[i]), match(names(all)[NA.cells], names(apes))]
  
}

## Calculate fat free mass (% fat estimates come from sources noted in data file)
all$FFM <- all$mass-(all$mass*all$body_fat)

# for E_i: (TEE-cost)/(mass^0.75), we use the high and low costs if we have them to generate intervals
# note this doesn't use any error propagation. 
all$E_i <- (all$TEE-all$cost)/(all$FFM^0.75)
all$E_i_low <- (all$TEE - all$hi.cost)/(all$FFM^0.75)
all$E_i_high <- (all$TEE - all$low.cost)/(all$FFM^0.75)


save.image(file = "analysis_data.RData")
