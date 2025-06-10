############## Convert time into energy and combine with production  ####################

# required packages
library(dplyr)
library(rstan)
library(rethinking)
library(shinystan)
library(tidyr)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)

set.seed(122)  # since there is some sampling to conform posteriors of different length, ensure consistency

p.hrs.female <- readRDS("female_hours_final.Rds")
p.hrs.male <- readRDS("male_hours_final.Rds")
ss <- readRDS("tsimane_energetic_summaries.Rds")  # Database of summarized energetic cost values of different Tsimane activities based on field respirometry. nkh abbreviation = net kcal/hr (energetic costs for individual of average mass, net of baseline costs of standing/sitting at rest)

## Info to match column numbers to activities

## For males
# 1 = Chop
# 2 = Clear
# 3 = Eat: 0.5*60 (kcal/hr: Value from Montgomery and Johnson 1977)
# 4 = Firewood/water collection
# 5 = Hoeing
# 6 = Food processing
# 7 = Tool manufacture
# 8 = Hunting/other subsistence activities characterized by normal paced walking
# 9 = Rice pounding (processing)
# 10 = Fishing
# 11 = Harvest
# 12 = Yucca digging (harvest)
# 13 = Non-subsistence activities
male.acts <- c("CHOP", "CLEAR", "EAT", "FWC", "HOE", "light_processing", "MANU", "NW", "RICE", "SW", "WR", "YUCCA", "nonsubsistence" )

## For females
# 1 = Chop
# 2 = Clear
# 3 = Eat:
# 4 = Firewood/water collection
# 5 = Food processing
# 6 = Tool manufacture
# 7 = Hunting/other subsistence activities characterized by normal paced walking
# 8 = Rice pounding (processing)
# 9 = Fishing
# 10 = Harvest
# 11 = Yucca digging (harvest)
# 12 = Non-subsistence activities
female.acts <- c("CHOP", "CLEAR", "EAT", "FWC", "light_processing", "MANU", "NW", "RICE", "SW", "WR", "YUCCA", "nonsubsistence" )

age_seq <- seq(15, 75, by=1)

#################################
##### CALCULATE MALE COSTS ########
###################################

# get the 7 activity mean costs (Note: non-subsistence is ignored)
male.costs <- filter(ss, Male == 1)
mc <- male.costs$mu.nkh[match(male.acts, male.costs$act.code)]

mc_all <- mc
mc_all[which(male.acts %in% c("nonsubsistence"))] <- 0 # add a zero at the end for nonsubsistence costs


# matrix multiply so that columns multiply by elements of vector
pe.list <- lapply(p.hrs.male, function(x) x %*% diag(mc_all))

#### Make stacked figure of costs by age
stack <- lapply(pe.list, colMeans, na.rm=T)
stacks <- as.data.frame(do.call(rbind, stack))
names(stacks) <- male.acts
stacks$age <- seq(15, 75, 1)

stacks$garden <- stacks$CHOP + stacks$CLEAR + stacks$HOE
stacks$process <- stacks$light_processing + stacks$RICE + stacks$EAT
stacks$harvest <- stacks$WR + stacks$YUCCA

stack_plot <- stacks %>% dplyr::select(-nonsubsistence, -CHOP, -CLEAR, -HOE, -light_processing, -RICE, -WR, -YUCCA) %>% gather(key="act", value="cost", FWC, MANU, NW, SW, garden, process, harvest, EAT)
stack_plot$act <- factor(stack_plot$act, levels = c("FWC", "garden", "harvest", "EAT", "MANU", "process", "NW",  "SW"))
male_cost <- ggplot(stack_plot, aes(x=age, y=cost, fill=act)) + geom_area() +
  scale_fill_manual(name="Activity", labels=c("Firewood and water collection", "Gardening", "Harvest", "Eating", "Manufacture", "Food processing", "Hunting/other foraging",  "Fishing"), values = brewer.pal(8, "Set1")) + labs(x="Age (years)", y="Energy cost (kcal/day)") +
  theme_classic(base_size=16) + ggtitle("Tsimane: Men") + lims(y=c(0,750))

# sum across rows to get subsistence costs at given age
costs_male <- lapply(pe.list, rowSums, na.rm=T)

# combine the list into a dataframe (columns are ages, rows are posterior samples)
costs_male <- do.call(cbind, costs_male)
age.foraging.costs.male <- apply(costs_male, 2, median)
age.foraging.costs.male.PI <- apply(costs_male, 2, HPDI, prob=0.95)

# male hours spent foraging
p.hrs.male.sub <- lapply(p.hrs.male, function(x) x[, -which(male.acts %in% c("nonsubsistence"))])
male.hrs.worked <- apply(do.call(cbind, lapply(p.hrs.male.sub, rowSums, na.rm=T)), 2, median)
male.hrs.worked.PI <- apply(do.call(cbind, lapply(p.hrs.male.sub, rowSums, na.rm=T)), 2, HPDI)

tot.male_all <- data.frame(age = seq(15, 75, 1), male = 1, cost = age.foraging.costs.male, low.cost = age.foraging.costs.male.PI[1,], hi.cost = age.foraging.costs.male.PI[2,], hrs.worked = male.hrs.worked, low.hrs.worked = male.hrs.worked.PI[1,], hi.hrs.worked = male.hrs.worked.PI[2,])



######################################
###### CALCULATE FEMALE COSTS #######
#####################################

# get the 7 activity mean costs (Note: non-subsistence is ignored)
female.costs <- filter(ss, Male == 0)
fc <- female.costs$mu.nkh[match(female.acts, female.costs$act.code)]

# At this point, we can save various configurations, including different 
# activities in subsistence costs by zeroing out those that we want to exclude
fc_all <- fc
fc_all[which(female.acts %in% c("nonsubsistence"))] <- 0 # add a zero at the end for nonsubsistence costs

# matrix multiply so that columns multiply by elements of vector
pe.list <- lapply(p.hrs.female, function(x) x %*% diag(fc_all))


#### Make stacked figure of costs by age
stack <- lapply(pe.list, colMeans, na.rm=T)
stacks <- as.data.frame(do.call(rbind, stack))
names(stacks) <- female.acts
stacks$age <- seq(15,75,1)

stacks$garden <- stacks$CHOP + stacks$CLEAR
stacks$process <- stacks$light_processing + stacks$RICE + stacks$EAT
stacks$harvest <- stacks$WR + stacks$YUCCA

stack_plot <- stacks %>% dplyr::select(-nonsubsistence, -CHOP, -CLEAR, -light_processing, -RICE, -WR, -YUCCA) %>% gather(key="act", value="cost", FWC, MANU, NW, SW, garden, process, harvest, EAT)
stack_plot$act <- factor(stack_plot$act, levels = c("FWC", "garden", "harvest", "EAT", "MANU", "process", "NW",  "SW"))
female_cost <- ggplot(stack_plot, aes(x=age, y=cost, fill=act)) + geom_area() +
  scale_fill_manual(name="Activity", labels=c("Firewood and water collection", "Gardening", "Harvest", "Eating", "Manufacture", "Food processing", "Hunting/other foraging",  "Fishing"), values = brewer.pal(8, "Set1")) +
  labs(x="Age (years)", y="Energy cost (kcal/day)") + theme_classic(base_size=16) +
  ggtitle("Tsimane: Women") + lims(y=c(0, 750))


# sum across rows to get subsistence costs at given age
costs_female <- lapply(pe.list, rowSums, na.rm=T)

# combine the list into a dataframe (columns are ages, rows are posterior samples)
costs_female <- do.call(cbind, costs_female)
age.foraging.costs.female <- apply(costs_female, 2, median)
age.foraging.costs.female.PI <- apply(costs_female, 2, HPDI, prob=0.95)

# female hours spent foraging
p.hrs.female.sub <- lapply(p.hrs.female, function(x) x[, -which(female.acts %in% c("nonsubsistence"))])
female.hrs.worked <- apply(do.call(cbind, lapply(p.hrs.female.sub, rowSums, na.rm=T)), 2, median)
female.hrs.worked.PI <- apply(do.call(cbind, lapply(p.hrs.female.sub, rowSums, na.rm=T)), 2, HPDI)

tot.female_all <- data.frame(age = seq(15,75,1), male=0, cost = age.foraging.costs.female, low.cost = age.foraging.costs.female.PI[1,], hi.cost = age.foraging.costs.female.PI[2,], hrs.worked = female.hrs.worked, low.hrs.worked = female.hrs.worked.PI[1,], hi.hrs.worked = female.hrs.worked.PI[2,])


tsimane <- bind_rows(tot.male_all, tot.female_all)


##############################################
##### combine with Tsimane production data #####
##############################################
library(brms)
tsimane_age_mean <- 19.09   # these are needed to back transform between age z-scores used in the model and natural scale 
tsimane_age_sd <- 18.10
tsimane_returns_model <- readRDS("tsimane_returns_model_gamma_spline.Rds")
tsimane_returns_model <- readRDS("tsimane_returns_model_gamma_spline_final.Rds")

# extract full posterior predictions for each level
tsi_dat <- fitted(tsimane_returns_model, 
            newdata= expand.grid(age_z = (age_seq-tsimane_age_mean)/tsimane_age_sd,  Male=c("1", "0")),
                      summary=F)  # note that the use of fitted vs predict means that the residual (observation-level) variance is not included


# Add Ea to dataset
tsimane$Ea <- apply(tsi_dat, 2, median)
tsimane$Ea.lo <- apply(tsi_dat, 2, HPDI, prob=0.95)[1,]
tsimane$Ea.hi <- apply(tsi_dat, 2, HPDI, prob=0.95)[2,]

# Calculate net Ea
# conform posterior first (due to different number of posterior samples)
costs_female_conform <- apply(costs_female, 2, function(x) sample(x, size=dim(costs_male)[1], replace = T))

net_Ea_mat <- tsi_dat - cbind(costs_male, costs_female_conform)

tsimane$net_EA <- apply(net_Ea_mat, 2, median)
tsimane$net_Ea_low <- apply(net_Ea_mat, 2, HPDI, prob=0.95)[1,]
tsimane$net_Ea_hi <- apply(net_Ea_mat, 2, HPDI, prob=0.95)[2,]

# Calculate efficiency (F)
EE_mat <- tsi_dat/cbind(costs_male, costs_female_conform)

tsimane$EE <- apply(EE_mat, 2, median)
tsimane$EE_low <- apply(EE_mat, 2, HPDI, prob=0.95)[1,]
tsimane$EE_hi <- apply(EE_mat, 2, HPDI, prob=0.95)[2,]

# Calculate Rg
male.time <- do.call(cbind, lapply(p.hrs.male.sub, rowSums, na.rm=T))
female.time <- do.call(cbind, lapply(p.hrs.female.sub, rowSums, na.rm=T))

female.time.conform <- apply(female.time, 2, function(x) sample(x, size=dim(tsi_dat)[1], replace = T))

Rg_mat <- tsi_dat/cbind(male.time, female.time.conform) 
tsimane$Rg <- apply(Rg_mat, 2, median)
tsimane$Rg_low <- apply(Rg_mat, 2, HPDI, prob=0.95)[1,]
tsimane$Rg_hi <- apply(Rg_mat, 2, HPDI, prob=0.95)[2,]

# Calculate Rn
Rn_mat <- net_Ea_mat/cbind(male.time, female.time.conform)

tsimane$Rn <- apply(Rn_mat, 2, median)
tsimane$Rn_low <- apply(Rn_mat, 2, HPDI, prob=0.95)[1,]
tsimane$Rn_hi <- apply(Rn_mat, 2, HPDI, prob=0.95)[2,]


# Save file
write.csv(tsimane, file = "tsimane_all.csv", row.names = F)
