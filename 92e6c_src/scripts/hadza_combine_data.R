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

set.seed(145)  # since there is some sampling to conform posteriors of different length, ensure consistency

load("hadza_female_hours_final_05_29_2020.RData")
load("hadza_male_hours_final_05_29_2020.RData")

## Info to match column numbers to activities

## For males
# 1 = Firewood/water collection
# 2 = Eat
# 3 = Food processing
# 4 = Tool manufacture
# 5 = non-subsistence
# 6 = out-of-camp
# 7 = pounding baobab

# column 6 is out of camp, so we are going to multiply that by the proportion time spent doing various out of camp things from  focal follows. We then cbind that back to the other columns, and rename the columns appropriately.
# breakdown of time out of camp:
# walking = 0.632
# other = 0.224
# digging = 0 (excluded from vector here)
# running = 0.001 (excluded from vector here, this amount added to walking)
# chopping = 0.049
# resting = 0.095

props_male <- c(0.632, 0.224, 0.049, 0.095)
p.hrs.male2 <- lapply(p.hrs.male, function(x) cbind(x[,c(1:5,7)], t(props_male %*% t(x[,6]))))

male.acts <- c("carry_firewood_water", "eat", "light_processing", "manufacture", "nonsubsistence", "pounding_baobab", 'walking', 'other_out_of_camp', 'chopping', 'resting')


###### Convert to energy costs #####
## Note that in all calculations below we generate net costs (the cost above and beyond standing or sitting at rest)

## Some outside data sources were used to fill in gaps:
# carry_firewood_water: (3.2-1.39)*60 -- Uses "Other" category in Montgomery & Johnson 1977 which includes getting firewood, fetching water, labor expended in exchange for trade good, and travel to trading posts. 
# eat: (1.9-1.39)*60 (Montgomery & Johnson 1977)
# light_processing: (1.9 - 1.39)*60  (numbers come from Montgomery & Johnson 1977 for Machiguenga. See table 5 - we use value for "Food Preparation" and subtract resting (sitting) to get net cost)
# manufacture: (2.7 - 1.15)*60 -- numbers come from Montgomery & Johnson 1977 for Machiguenga. Here we use the value for bow and arrow manufacture because that is what represents most of tsimane male manufacture for men. We use the value from table 5 for "Manufacture" for women. We subtract resting (sitting) to get net cost.
# nonsubsistence: 0
# pounding_baobab: 2.25 kcal/min (This value is assumed to be equal to the metabolic cost of digging)
# walking:  0.53 kcal/kg/km * 50.9 kg / 1000 m/km * 1.19 m/s * 60 s/min * 60 min/hr  (note that this starts with a cost of transport value in units of kcal/kg/km, which is then multiplied by average Hadza male mass and walking speed to get costs per unit time)
# other_out_of_camp: (1.9 - 1.39)*60 -- notes from focal follows indicate that this "Other out of camp" activity category mostly included processing and relatively low expenditure activities, so the value above for processing is assumed
# chopping: 5.05 kcal/min
# resting (out of camp): this is essentially non-subsistence, so cost = 0
# We also estimate that men climb 10 meters per day on average. Metabolic cost = 0.015 kcal/kg/meter, mass = 54.4 kg
male.cost.vector <- c((3.2-1.39)*60, (1.9-1.39)*60, (1.9 - 1.39)*60, (2.7 - 1.15)*60, 0, 2.25*60, 0.53*50.9*1.19*60 * 60 / 1000, (1.9 - 1.39)*60, 5.05*60, 0)

## For females (order of activity categories)
# 1 = Firewood/water collection
# 2 = Eat
# 3 = Food processing
# 4 = Tool manufacture
# 5 = non-subsistence
# 6 = out-of-camp
# 7 = pounding baobab

# column 6 is out of camp, so we are going to multiply that by the proportion time spent doing various out of camp things from focal follows. We then cbind that back to the other columns, and rename the columns appropriately.
# breakdown of time out of camp:
# walking = 0.306
# other = 0.206
# digging = 0.386
# running = 0.000 (excluded from vector here)
# chopping = 0.000 (excluded from vector here)
# resting = 0.101 
props_female <- c(0.306, 0.206, 0.386, 0.101)
p.hrs.female2 <- lapply(p.hrs.female, function(x) cbind(x[,c(1:5,7)], t(props_female %*% t(x[,6]))))

female.acts <- c("carry_firewood_water", "eat", "light_processing", "manufacture", "nonsubsistence", "pounding_baobab", 'walking', 'other_out_of_camp', 'digging', 'resting')

## Convert to energy costs
## Note that in all calculations below we generate net costs (the cost above and beyond standing or sitting at rest)

# carry_firewood_water: (1.8-1.39)*60 -- Uses "Other" category in Montgomery & Johnson 1977 which includes getting firewood, fetching water, labor expended in exchange for trade good, and travel to trading posts. 
# eat: (1.3-1.15)*60 (Montgomery & Johnson 1977)
# light_processing: (1.9 - 1.15)*60   ## numbers come from Montgomery & Johnson 1977 for Machiguenga. See table 5 - we use value for "Food Preparation" and subtract resting (sitting) to get net cost.
# manufacture: (1.7 - 1.15)*60 -- numbers come from Montgomery & Johnson 1977 for Machiguenga. We use the value from table 5 for "Manufacture" for women. We subtract resting (sitting) to get net cost. 
# nonsubsistence: 0
# pounding_baobab: 2.25 kcal/min (This value is assumed to be equal to the metabolic cost of digging)
# walking: 0.53 kcal/kg/km * 43.4 kg / 1000 m/km * 1.03 m/s * 60 s/min * 60 min/hr (note that this starts with a cost of transport value in units of kcal/kg/km, which is then multiplied by average Hadza female mass and walking speed to get costs per unit time)
# other_out_of_camp: (1.9 - 1.15)*60 -- notes from focal follows indicate that this "Other out of camp" activity category mostly included processing and relatively low expenditure activities, so the value above for processing is assumed
# digging: 2.25 kcal/min
# resting (out of camp): this is essentially non-subsistence, so cost = 0
female.cost.vector <- c((1.8-1.39)*60, (1.3-1.15)*60, (1.9 - 1.15)*60, (1.7 - 1.15)*60, 0, 2.25*60, 0.53 * 43.4 * 1.03 * 60 * 60 / 1000 , (1.9 - 1.15)*60, 2.25*60, 0)

age_seq <- seq(15, 75, 1)


##### CALCULATE  COSTS

# At this point, we can save various configurations, including different 
# activities in subsistence costs by zeroing out those that we want to exclude

# matrix multiply so that columns multiply by elements of vector
pe.list.hadza_male <- lapply(p.hrs.male2, function(x) x %*% diag(male.cost.vector))
pe.list.hadza_male <- lapply(pe.list.hadza_male, function(x) cbind(x, 0.015*50.9*10)) # Add in fixed climbing cost of 10 m/day: 0.015 kcal/kg/meter * 50.9 kg * 10 meter/day
pe.list.hadza_female <- lapply(p.hrs.female2, function(x) x %*% diag(female.cost.vector))


#### Make stacked figure of costs by age
stack <- lapply(pe.list.hadza_male, colMeans, na.rm=T)
stacks <- as.data.frame(do.call(rbind, stack))
names(stacks) <- c(male.acts, "climbing")
stacks$age <- age_seq

stacks$Processing <- stacks$light_processing + stacks$other_out_of_camp + stacks$pounding_baobab

stack_plot_male <- stacks %>% select(-nonsubsistence, -resting, -light_processing, -other_out_of_camp, -pounding_baobab) %>% gather(key="act", value="cost", carry_firewood_water:climbing, Processing)
male_stack_hadza <- ggplot(stack_plot_male, aes(x=age, y=cost, fill=act)) + geom_area() + scale_fill_manual(name="Activity", labels=c("Firewood and water collection", "Chopping", "Climbing", "Eating", "Manufacture", "Food processing", "Walking (during foraging)"), values = brewer.pal(7, "Set1")) + labs(x="Age (years)", y="Energy cost (kcal/day)") + theme_classic(base_size=16) +
  ggtitle("Hadza: Men") + lims(y=c(0, 750))

stack <- lapply(pe.list.hadza_female, colMeans, na.rm=T)
stacks <- as.data.frame(do.call(rbind, stack))
names(stacks) <- female.acts
stacks$age <- age_seq

stacks$Processing <- stacks$light_processing + stacks$other_out_of_camp + stacks$pounding_baobab

stack_plot_female <- stacks %>% select(-nonsubsistence, -resting, -light_processing, -other_out_of_camp, -pounding_baobab) %>% gather(key="act", value="cost", carry_firewood_water:digging, Processing)
fhcols <- brewer.pal(8, "Set1")[c(1,2,4,5,6,7)]
fhcols[2] <- "#33CCFF"
female_stack_hadza <- ggplot(stack_plot_female, aes(x=age, y=cost, fill=act)) + geom_area() + scale_fill_manual(name="Activity", labels=c("Firewood and water collection", "Digging", "Eating", "Manufacture", "Food processing", "Walking (during foraging)"), values = fhcols) + labs(x="Age (years)", y="Energy cost (kcal/day)") + theme_classic(base_size=16) +
  ggtitle("Hadza: Women") + lims(y=c(0, 750))




# final male
# sum across rows to get subsistence costs at given age
costs_male <- lapply(pe.list.hadza_male, rowSums, na.rm=T)

# combine the list into a dataframe (columns are ages, rows are posterior samples)
costs_male <- do.call(cbind, costs_male)
age.foraging.costs.male <- apply(costs_male, 2, median)
age.foraging.costs.male.PI <- apply(costs_male, 2, HPDI, prob=0.95)

# male hours spent working
p.hrs.male.sub <- lapply(p.hrs.male2, function(x) x[, -which(male.acts %in% c("nonsubsistence", "resting"))])
male.hrs.worked <- apply(do.call(cbind, lapply(p.hrs.male.sub, rowSums, na.rm=T)), 2, median)
male.hrs.worked.PI <- apply(do.call(cbind, lapply(p.hrs.male.sub, rowSums, na.rm=T)), 2, HPDI)

tot.male_all <- data.frame(age = age_seq, male = 1, cost = age.foraging.costs.male, low.cost = age.foraging.costs.male.PI[1,], hi.cost = age.foraging.costs.male.PI[2,], hrs.worked = male.hrs.worked, low.hrs.worked = male.hrs.worked.PI[1,], hi.hrs.worked = male.hrs.worked.PI[2,])



# female final
# sum across rows to get subsistence costs at given age
costs_female <- lapply(pe.list.hadza_female, rowSums, na.rm=T)

# combine the list into a dataframe (columns are ages, rows are posterior samples)
costs_female <- do.call(cbind, costs_female)
age.foraging.costs.female <- apply(costs_female, 2, median)
age.foraging.costs.female.PI <- apply(costs_female, 2, HPDI, prob=0.95)


# female hours spent foraging
p.hrs.female.sub <- lapply(p.hrs.female2, function(x) x[, -which(female.acts %in% c("nonsubsistence", "resting"))])
female.hrs.worked <- apply(do.call(cbind, lapply(p.hrs.female.sub, rowSums, na.rm=T)), 2, median)
female.hrs.worked.PI <- apply(do.call(cbind, lapply(p.hrs.female.sub, rowSums, na.rm=T)), 2, HPDI)

tot.female_all <- data.frame(age = age_seq, male=0, cost = age.foraging.costs.female, low.cost = age.foraging.costs.female.PI[1,], hi.cost = age.foraging.costs.female.PI[2,], hrs.worked = female.hrs.worked, low.hrs.worked = female.hrs.worked.PI[1,], hi.hrs.worked = female.hrs.worked.PI[2,])

# combine male/female
total_hadza <- bind_rows(tot.male_all, tot.female_all)

##############################################
##### combine with Hadza production data #####
##############################################
library(brms)
hadza_age_mean <- 30.94   # these are needed to back transform between age z-scores used in the model and natural scale 
hadza_age_sd <- 18.79
load("C:/Users/tskra/Dropbox/energetics_tsimane_hadza/manuscript/submission2_science/revise_and_resubmit/published data/fitted_models/hadza_lognormal_returns_model.RData")
load("/Users/thomaskraft/Dropbox/energetics_tsimane_hadza/manuscript/submission2_science/revise_and_resubmit/published data/fitted_models/hadza_lognormal_returns_model.RData")

# extract full posterior predictions for each level
had_rr_dat <- fitted(hadza_lognormal_returns_model, 
                     newdata= expand.grid(age_z = (age_seq-hadza_age_mean)/hadza_age_sd,  sex=c("M", "F")),
                     re_formula=NA, summary=F)  # note that the use of fitted vs predict means that the residual (observation-level) variance is not included


## Add in out of camp production: We start with published observations that adult hadza men get ~306 extra 
## kcal/hr when foraging out of camp than is returned, and hadza women get ~70 kcal/hr more than is returned
# men: 306 kcal/hr 
# women: 70 kcal/hr
## We assume these to be peak values for adults, so we scale these down according to our empirical estimates 
## by age here as follows:
out_of_camp_scaling <- apply(had_rr_dat, 2, mean)
hadza_male_out_of_camp_Ea_per_hr <- 306*(out_of_camp_scaling[1:61]/max(out_of_camp_scaling[1:61]))  # men
hadza_female_out_of_camp_Ea_per_hr <- 70*(out_of_camp_scaling[62:122]/max(out_of_camp_scaling[62:122])) # women

# multiply by the time spent foraging (columns) by vector of per hour out of camp returns using matrix multiplication
male.time.out <- do.call(cbind, lapply(p.hrs.male, "[", , c(6))) # (recall that that column 6 in this matrix corresponds to out of camp time)
hadza_male_out_of_camp_kcal_per_day <- male.time.out %*% diag(hadza_male_out_of_camp_Ea_per_hr)

female.time.out <- do.call(cbind, lapply(p.hrs.female, "[", , c(6))) # (recall that that column 6 in this matrix corresponds to out of camp time)
hadza_female_out_of_camp_kcal_per_day <- female.time.out %*% diag(hadza_female_out_of_camp_Ea_per_hr)

# add out of camp returns to those brought back to camp
# we must first conform the matrices because they the models had different numbers of posterior samples
num_rows <- dim(had_rr_dat)[1]
hadza_male_out_of_camp_kcal_per_day_conform <- apply(hadza_male_out_of_camp_kcal_per_day, 2, function(x) sample(x, size=num_rows, replace = T))
hadza_female_out_of_camp_kcal_per_day_conform <- apply(hadza_female_out_of_camp_kcal_per_day, 2, function(x) sample(x, size=num_rows, replace = T))

# add together
had_rr_dat[,1:length(age_seq)] <- had_rr_dat[,1:length(age_seq)] + hadza_male_out_of_camp_kcal_per_day_conform
had_rr_dat[,(length(age_seq)+1):(2*length(age_seq))] <- had_rr_dat[,(length(age_seq)+1):(2*length(age_seq))] + hadza_female_out_of_camp_kcal_per_day_conform

# Add Ea to dataset
total_hadza$Ea <- apply(had_rr_dat, 2, median)
total_hadza$Ea.lo <- apply(had_rr_dat, 2, HPDI, prob=0.95)[1,]
total_hadza$Ea.hi <- apply(had_rr_dat, 2, HPDI, prob=0.95)[2,]

# Calculate net Ea
male.costs.conform <- apply(costs_male, 2, function(x) sample(x, size=num_rows, replace = T))
female.costs.conform <- apply(costs_female, 2, function(x) sample(x, size=num_rows, replace = T))

net_Ea_mat <- had_rr_dat - cbind(male.costs.conform, female.costs.conform)

total_hadza$net_EA <- apply(net_Ea_mat, 2, median)
total_hadza$net_Ea_low <- apply(net_Ea_mat, 2, HPDI, prob=0.95)[1,]
total_hadza$net_Ea_hi <- apply(net_Ea_mat, 2, HPDI, prob=0.95)[2,]

# Calculate efficiency (F)
EE_mat <- had_rr_dat/cbind(male.costs.conform, female.costs.conform)

total_hadza$EE <- apply(EE_mat, 2, median)
total_hadza$EE_low <- apply(EE_mat, 2, HPDI, prob=0.95)[1,]
total_hadza$EE_hi <- apply(EE_mat, 2, HPDI, prob=0.95)[2,]

# Calculate Rg
male.time <- do.call(cbind, lapply(p.hrs.male.sub, rowSums, na.rm=T))
female.time <- do.call(cbind, lapply(p.hrs.female.sub, rowSums, na.rm=T))

male.time.conform <- apply(male.time, 2, function(x) sample(x, size=num_rows, replace = T))
female.time.conform <- apply(female.time, 2, function(x) sample(x, size=num_rows, replace = T))

Rg_mat <- had_rr_dat/cbind(male.time.conform, female.time.conform) 
total_hadza$Rg <- apply(Rg_mat, 2, median)
total_hadza$Rg_low <- apply(Rg_mat, 2, HPDI, prob=0.95)[1,]
total_hadza$Rg_hi <- apply(Rg_mat, 2, HPDI, prob=0.95)[2,]

# Calculate Rn
Rn_mat <- net_Ea_mat/cbind(male.time.conform, female.time.conform)

total_hadza$Rn <- apply(Rn_mat, 2, median)
total_hadza$Rn_low <- apply(Rn_mat, 2, HPDI, prob=0.95)[1,]
total_hadza$Rn_hi <- apply(Rn_mat, 2, HPDI, prob=0.95)[2,]



# Save file
write.csv(total_hadza, file = "hadza_all.csv", row.names = F)
