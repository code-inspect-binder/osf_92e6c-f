# Analysis of hadza food production by age/sex

# required packages
library(sjPlot)
library(lubridate)
library(tidyverse)
library(brms)


## Load data: data file "hfp" contains individual-level foraging return data (kcal) at the daily level along with the date, unique camp id, age, sex, and unique person id. 
hfp <- read.csv("food_pro_data.csv", as.is=T)

h <- hfp[hfp$age > 5,]
h$age_z <- (h$age-mean(h$age, na.rm=T))/sd(h$age, na.rm=T)
h$month <- month(ymd(h$date))

# basic plot of data
ggplot(h, aes(x=age, y=log(kcal+1), col=sex)) + geom_jitter(alpha=0.1) + stat_smooth()

# prior predictive checks 
# PPC's are useful here for defining weakly informative priors for polynomial (cubic) lognormal model, with the central goal of constraining implausibly high values given the lognormal distribution. We know from domain knowledge deriving from studies of savannah foragers that most outcomes will range between 0-20,000 kcal/day, with some very high returns mixed in from large game. As can be seen, these priors remain quite weak but still allow some values beyond what is likely in reality. Ultimately this constrains the parameter space to prevent poor convergence but model outcomes are driven primarily by the data.
fit_test <-
  brm(
    bf(kcal~age_z*sex + I(age_z^2)*sex + I(age_z^3)*sex + (1|person_id) + (1|camp) + (1|month), 
       hu ~age_z*sex + I(age_z^2)*sex + I(age_z^3)*sex + (1|person_id) + (1|camp) + (1|month)), 
    data=h,
    family=hurdle_lognormal,
    sample_prior = 'only',
    warmup = 500,
    iter = 1000,
    cores = 4,
    control = list(adapt_delta = 0.95, max_treedepth = 10),
    prior = c(set_prior('normal(0, 2)', coef = "age_z"),
              set_prior('normal(0, 3)', coef = "age_z:sexM"),
              set_prior('normal(0, 3)', coef = "Iage_zE2"),
              set_prior('normal(0, 3)', coef = "Iage_zE3"),
              set_prior('normal(0, 5)', coef = "sexM"),
              set_prior('normal(0, 5)', coef = "sexM:Iage_zE2"),
              set_prior('normal(0, 5)', coef = "sexM:Iage_zE3"),
              set_prior('student_t(3, 0, 2)', class="sd"),
              set_prior('normal(8, 2)', class = "Intercept"),
              set_prior('normal(0,3)', class="b", dpar="hu")))
fit_test_cd <- conditional_effects(fit_test, effects=c("age_z:sex")) 
plot(fit_test_cd)[[1]] + ggplot2::lims(x=c(-1,2), y=c(0,20000))


## Zero-adjusted Gamma or Lognormal model perform well, see Koster & McElreath 2014
## Run model
hadza_lognormal_returns_model <- brm(bf(kcal~age_z*sex + I(age_z^2)*sex + I(age_z^3)*sex + (1|person_id) + (1|camp) + (1|month), 
                     hu ~age_z*sex + I(age_z^2)*sex + I(age_z^3)*sex + (1|person_id) + (1|camp) + (1|month)),
                                           data=h,
                                           cores = 4, chains=4, iter = 3000, warmup = 1500, 
                                           family=hurdle_lognormal, 
                                           control = list(adapt_delta = 0.95),
                  prior = c(set_prior('normal(0, 2)', coef = "age_z"),
                            set_prior('normal(0, 3)', coef = "age_z:sexM"),
                            set_prior('normal(0, 3)', coef = "Iage_zE2"),
                            set_prior('normal(0, 3)', coef = "Iage_zE3"),
                            set_prior('normal(0, 5)', coef = "sexM"),
                            set_prior('normal(0, 5)', coef = "sexM:Iage_zE2"),
                            set_prior('normal(0, 5)', coef = "sexM:Iage_zE3"),
                            set_prior('student_t(3, 0, 2)', class="sd"),
                            set_prior('normal(8, 2)', class = "Intercept"),
                            set_prior('normal(0,3)', class="b", dpar="hu")))
# convergence of chains, Rhat, and effective sample size all look good

# alternative gamma model, 
had_mod1_gamma <- brm(bf(kcal~age_z*sex + I(age_z^2)*sex + I(age_z^3)*sex + (1|person_id) + (1|camp) + (1|month), 
                     hu ~age_z*sex + I(age_z^2)*sex + I(age_z^3)*sex + (1|person_id) + (1|camp) + (1|month)),
                  data=h,
                  cores = 4, chains=4, iter = 1000, warmup = 500, 
                  family=hurdle_gamma, 
                  control = list(adapt_delta = 0.95),
                  prior = c(set_prior('normal(0, 1)', coef = "age_z"),
                            set_prior('normal(0, 1)', coef = "age_z:sexM"),
                            set_prior('normal(0, 1)', coef = "Iage_zE2"),
                            set_prior('normal(0, 1)', coef = "Iage_zE3"),
                            set_prior('normal(0, 3)', coef = "sexM"),
                            set_prior('normal(0, 1)', coef = "sexM:Iage_zE2"),
                            set_prior('normal(0, 1)', coef = "sexM:Iage_zE3"),
                            set_prior('student_t(3, 0, 2)', class="sd"),
                            set_prior('normal(8, 2)', class = "Intercept"),
                            set_prior('normal(0,2)', class="b", dpar="hu")))



save(hadza_lognormal_returns_model, file="../hadza_lognormal_returns_model.RData") # hurdle lognormal version
save(had_mod1_gamma, file="../had_mod1_gamma.RData")  #hurdle gamma version

# Plot the marginal/conditional effects for predicted values and probability of zero-days
p <- conditional_effects(hadza_lognormal_returns_model, effects=c("age_z:sex"))  # includes the hurdle component
p2 <- conditional_effects(hadza_lognormal_returns_model, effects=c("age_z:sex"), dpar="hu") # to see the hurdle component (Pr zero-days)

plot(p)[[1]] + scale_x_continuous(breaks = c(-1, 0, 1, 2), labels= round(c(-1, 0, 1, 2)*sd(h$age, na.rm=T) + mean(h$age, na.rm=T))) + ggplot2::labs(x="Age", y="kcal/day") +ggplot2::lims(y=c(0,6500)) + theme_classic()
plot(p2)[[1]] + ggplot2::lims(y=c(0,1)) + labs(x="Age", y="Probability of zero day\n(hurdle component of model)") + theme_classic() + scale_x_continuous(breaks = c(-1, 0, 1, 2), labels= round(c(-1, 0, 1, 2)*sd(h$age, na.rm=T) + mean(h$age, na.rm=T))) 


# Posterior predictive check
pp_check(hadza_lognormal_returns_model, nsamples = 50) + scale_x_continuous(trans='log1p')
pp_check(had_mod1_gamma, nsamples = 50) + scale_x_continuous(trans='log1p')

pp_check(hadza_lognormal_returns_model, type = "stat", stat="prop_zero", nsamples = 50)
pp_check(had_mod1_gamma, type = "stat", stat="prop_zero", nsamples = 50)


# Compare IC
hadza_lognormal_returns_model <- add_criterion(hadza_lognormal_returns_model, criterion = c("loo", "waic"))
had_mod1_gamma <- add_criterion(had_mod1_gamma, criterion = c("loo", "waic"))
loo_compare(hadza_lognormal_returns_model, had_mod1_gamma)  # both models demonstrate reasonable fit, but the lognormal is clearly preferred on the basis of IC comparison. From a theoretical standpoint the lognormal is also generally heavier tailed by some metrics than the gamma, which is likely to be useful to capture some extreme values from men hunting large game. 


# extract predictions with error estimates from model
age_seq <- (seq(10, 75, 1)-mean(h$age, na.rm=T))/sd(h$age, na.rm=T)
had_rr_dat <- fitted(hadza_lognormal_returns_model, newdata= expand.grid(age_z = age_seq,  sex=c("M", "F")), re_formula=NA)  # note that the use of fitted vs predict means that the residual (observation-level) variance is not included
had_rr_dat <- cbind(had_rr_dat, expand.grid(age_z = age_seq,  sex=c("M", "F")), age=seq(10, 75, 1))
names(had_rr_dat) <- c("Ea", "Ea_se", "Ea.lo", "Ea.hi", "age_z", "sex", "age")

# save data to be used in other scripts
write.csv(had_rr_dat, file="../hadza_return_data_lnorm.csv", row.names = F)