library(dplyr)
library(rstan)
library(rethinking)
library(shinystan)
library(tidyr)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

############ MALE TSIMANE #################

# #######################################################################################################
# ## The first step toward plotting predicted probabilities of our fixed effects is to extract the
# ## posterior samples using the extract.samples function in the STAN package.
# ## Because of the memory requirements for this operation, it was run on the HPC and the 
# ## posterior samples were exported directly.
# ## The following command was used to extract samples:

# post <- extract.samples(model_fit)

# #######################################################################################################

load("./post_male_final.RData")

## We can visually inspect the individual-level random effects by creating a new object, v_est.
## The pairs function in the base R package allows us to visualize the pairwise correlation of the
## random effects. Based on the analysis of Koster/McElreath, a logical thing to do here would be
## to examine the model output using the precis function and then examine those with strong correlations. 
## Subsequently, the dens function allows us to see the posterior distribution of the correlation
## between those categories.

# male tsimane activities in order they appear in columns of posterior
male_activs <- c("CHOP","CLEAR","EAT","FWC","HOE","light_processing","MANU","NW","RICE","SW","WR","YUCCA","z")  
v_est <- apply(post$v_id,2:3,median)
pairs(v_est)
plot(v_est[,5],v_est[,7])
dens(post$Rho_id[,5,7])

# mean and SD of male ages for z-score conversion
sd_age_male <- 15.32
mean_age_male <- 34.26

#######################################################################################################
## Multinomial version of link function for the Tsimane (male) model
## Much like the link function in the rethinking package, the following function can be used to
## generate predictions for a customized data frame. Conventionally, researchers allow one of the
## predictors to vary while holding other predictors at a constant value (often the sample mean).
## The script begins by defining a sequence length that should correspond to the dimensions of
## data frame that is created to generate the predictions. We use 100 in this case, but other
## values are possible.
## We define additional quantities at the beginning of the function, and the function derives
## the values from the posterior. For example, K is the number of response categories.
## The quantity, ns, is the number of samples.
## The function relies on an intermediate step, the creation of ptemp, which is tied to the likelihood
## function. Essentially, the function works by taking the supplied values for each of the i rows
## in the newly created data frame, then multipying those values by the value of the parameter in each
## sample from the posterior. In other words, for each of the samples from the posterior, the 
## parameters vary, which are used to generate predicted probabilities for each of the supplied
## combinations of values in the new data frame that we will create. We can subsequently 
## use the distribution of predicted values to understand the effets of the values we supply.
## It is common for researchers to average over the random effects. In other words, predicted
## probabilies are based only on the fixed effects, not accounting for the random effects.
## To accomplish that, our function requires users to supply a 0 for the random effects. By contrast,
## one can include random effects by supplying the integer that corresponds to the unit of interest
## in the higher-level classification. The varying intercept for that unit will then be incorporated
## as part of the predicted probability.
#######################################################################################################

seq.length <- length(seq(15, 75, 1))

link.mn <- function( data ) {
  K <- dim(post$v_id)[3] + 1
  ns <- dim(post$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post$a[,k] + 
          post$bA[,k] * data$age_z[i] + 
          post$bQ[,k] * data$age_zq[i] + 
          post$bT[,k] * data$time_z[i] + 
          post$bTQ[,k] * data$time_zq[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post$v_id[,data$id[i],k]
        if ( data$com_id[i]>0 ) ptemp <- ptemp + post$v_com[,data$com_id[i],k]
        if ( data$month_id[i]>0 ) ptemp <- ptemp + post$v_month[,data$month_id[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}



## The following generates a vector of ages which we use to generate predicted probabilities.
## We use the seq function to generate a sequence across the range of ages in the empirical data.
## Note that the length of this vector is seq.length, as defined above.
## Set the target age range to generate 100 points between ages 15 and 75
low_age <- (15 - mean_age_male)/sd_age_male 
high_age <- (75 - mean_age_male)/sd_age_male

age_seq <- seq (from= low_age, to= high_age, length.out = seq.length)

## We create a data frame from age_seq and its second order polynomial, holding most of predictors
## at their sample mean. In the original Koster+McElreath paper, they specified 8:00 am for the time
## of day in order to effect a more sensible alignment with the empirical data. In this case, we 
## hold that parameter at 0, which therefore corresponds with the mean time of day in our dataset
## for behavioral observations (0.532 = ~12:45).
## Also, as we noted earlier in the notes about the multinomial link function, we set "id" to 
## zero in order to average over the random effects.
pred_dat <- data.frame(
  id = 0 ,
  age_z = age_seq,
  age_zq = age_seq^2,
  time_z = 0,
  time_zq = 0,
  com_id = 0,
  month_id = 0
)

## The following is an object with 100 predicted probabilities (as in seq.length) corresponding
## to the different values of standardized ages that were supplied in the pred_dat data frame,
## as repeated for each of the K response categories.
p <- link.mn ( pred_dat )

## The following calculates the mean of the predicted samples.
p_mean <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )

## Alternative way to do this. I find the syntax more intuitive than the alternative, so could help others see what is done.
# do.call(cbind, lapply(p, function(i) apply(i, 2, mean))) 

###############################
# Calculating predicted probabilities at a single time of day cannot tell
# us about how much time an individual of a given age/sex/whatever class 
# is expected to spend on a particular behavior over the course of a "normal"
# day. This is particularly the case if you had more complex functional forms
# for the time parameter, and behaviors that are bi- or even tri-modal over the
# course of a day.

# The following code implements a procedure that at each draw from the 
# posterior integrates across the time of day from 7 am to 7 pm (this is
# when our samples were collected, and when the behaviors of interest occur).
# Thus, the output of this new link function is a composite "total time" spent per day
# across the range specified (could be whole day if measures spanned it).

# Nothing fancy is done here, but an integration is performed within the sapply() call 
# across the day and night interval specified in the function call, using the fitted model and
# our input data at which we want to predict. 

link.day.hours <- function( data, post, morning, night ) {
  K <- dim(post$v_id)[3] + 1
  ns <- dim(post$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- length(age_seq)
  
  p <- list()
  
  for ( i in 1:n ) {
    
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        sapply( 1:ns, function(j) {    
          integrate(
            function(x) 
exp(post$a[j,k] + post$bA[j,k]*data$age_z[i] + post$bQ[j,k]*data$age_zq[i] + post$bT[j,k]*x+ post$bTQ[j,k]*x^2)/
              ( 1 + (
exp(post$a[j,1]+ post$bA[j,1]*data$age_z[i] + post$bQ[j,1]*data$age_zq[i]+ post$bT[j,1]*x+ post$bTQ[j,1]*x^2) +
exp(post$a[j,2]+ post$bA[j,2]*data$age_z[i] + post$bQ[j,2]*data$age_zq[i]+ post$bT[j,2]*x+ post$bTQ[j,2]*x^2) +
exp(post$a[j,3]+ post$bA[j,3]*data$age_z[i] + post$bQ[j,3]*data$age_zq[i]+ post$bT[j,3]*x+ post$bTQ[j,3]*x^2) +
exp(post$a[j,4]+ post$bA[j,4]*data$age_z[i] + post$bQ[j,4]*data$age_zq[i]+ post$bT[j,4]*x+ post$bTQ[j,4]*x^2) +
exp(post$a[j,5]+ post$bA[j,5]*data$age_z[i] + post$bQ[j,5]*data$age_zq[i]+ post$bT[j,5]*x+ post$bTQ[j,5]*x^2) +
exp(post$a[j,6]+ post$bA[j,6]*data$age_z[i] + post$bQ[j,6]*data$age_zq[i]+ post$bT[j,6]*x+ post$bTQ[j,6]*x^2) +
exp(post$a[j,7]+ post$bA[j,7]*data$age_z[i] + post$bQ[j,7]*data$age_zq[i]+ post$bT[j,7]*x+ post$bTQ[j,7]*x^2) +
exp(post$a[j,8]+ post$bA[j,8]*data$age_z[i] + post$bQ[j,8]*data$age_zq[i]+ post$bT[j,8]*x+ post$bTQ[j,8]*x^2) +
exp(post$a[j,9]+ post$bA[j,9]*data$age_z[i] + post$bQ[j,9]*data$age_zq[i]+ post$bT[j,9]*x+ post$bTQ[j,9]*x^2) +
exp(post$a[j,10]+ post$bA[j,10]*data$age_z[i] + post$bQ[j,10]*data$age_zq[i]+ post$bT[j,10]*x+ post$bTQ[j,10]*x^2) +
exp(post$a[j,11]+ post$bA[j,11]*data$age_z[i] + post$bQ[j,11]*data$age_zq[i]+ post$bT[j,11]*x+ post$bTQ[j,11]*x^2) +
exp(post$a[j,12]+ post$bA[j,12]*data$age_z[i] + post$bQ[j,12]*data$age_zq[i]+ post$bT[j,12]*x+ post$bTQ[j,12]*x^2) 
              )),lower=morning, upper=night)$value * (12/(night-morning))
        })
        
        # could add back in if needed, but not needed unless predicting at specific level of random effect
        # if ( data$id[i]>0 ) ptemp <- ptemp + post$v_id[,data$id[i],k]
        # if ( data$com_id[i]>0 ) ptemp <- ptemp + post$v_com[,data$com_id[i],k]
        # if ( data$month_id[i]>0 ) ptemp <- ptemp + post$v_month[,data$month_id[i],k]
      } else {
        sapply( 1:ns, function(j) { 
          integrate(
            function(x) 1/( 1 + (
exp(post$a[j,1]+ post$bA[j,1]*data$age_z[i] + post$bQ[j,1]*data$age_zq[i]+ post$bT[j,1]*x+ post$bTQ[j,1]*x^2) +
exp(post$a[j,2]+ post$bA[j,2]*data$age_z[i] + post$bQ[j,2]*data$age_zq[i]+ post$bT[j,2]*x+ post$bTQ[j,2]*x^2) +
exp(post$a[j,3]+ post$bA[j,3]*data$age_z[i] + post$bQ[j,3]*data$age_zq[i]+ post$bT[j,3]*x+ post$bTQ[j,3]*x^2) +
exp(post$a[j,4]+ post$bA[j,4]*data$age_z[i] + post$bQ[j,4]*data$age_zq[i]+ post$bT[j,4]*x+ post$bTQ[j,4]*x^2) +
exp(post$a[j,5]+ post$bA[j,5]*data$age_z[i] + post$bQ[j,5]*data$age_zq[i]+ post$bT[j,5]*x+ post$bTQ[j,5]*x^2) +
exp(post$a[j,6]+ post$bA[j,6]*data$age_z[i] + post$bQ[j,6]*data$age_zq[i]+ post$bT[j,6]*x+ post$bTQ[j,6]*x^2) +
exp(post$a[j,7]+ post$bA[j,7]*data$age_z[i] + post$bQ[j,7]*data$age_zq[i]+ post$bT[j,7]*x+ post$bTQ[j,7]*x^2) +
exp(post$a[j,8]+ post$bA[j,8]*data$age_z[i] + post$bQ[j,8]*data$age_zq[i]+ post$bT[j,8]*x+ post$bTQ[j,8]*x^2) +
exp(post$a[j,9]+ post$bA[j,9]*data$age_z[i] + post$bQ[j,9]*data$age_zq[i]+ post$bT[j,9]*x+ post$bTQ[j,9]*x^2) +
exp(post$a[j,10]+ post$bA[j,10]*data$age_z[i] + post$bQ[j,10]*data$age_zq[i]+ post$bT[j,10]*x+ post$bTQ[j,10]*x^2) +
exp(post$a[j,11]+ post$bA[j,11]*data$age_z[i] + post$bQ[j,11]*data$age_zq[i]+ post$bT[j,11]*x+ post$bTQ[j,11]*x^2) +
exp(post$a[j,12]+ post$bA[j,12]*data$age_z[i] + post$bQ[j,12]*data$age_zq[i]+ post$bT[j,12]*x+ post$bTQ[j,12]*x^2) 
            )),lower=morning, upper=night)$value * (12/(night-morning))
        })
      }
      #return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    #for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}

mean_time_prop <- 0.5348547
sd_time_prop <- 0.1473534

p.hrs.male <- link.day.hours(data=pred_dat, post=post, morning = ((7/24)-mean_time_prop)/sd_time_prop, night = ((19/24)-mean_time_prop)/sd_time_prop)


saveRDS(p.hrs.male, file="male_hours_final.Rds")


rm(list=ls())


###########################################
####### FOR FEMALE TSIMANE ################
###########################################

# #######################################################################################################
# ## The first step toward plotting predicted probabilities of our fixed effects is to extract the
# ## posterior samples using the extract.samples function in the STAN package.
# ## Because of the memory requirements for this operation, it was run on the cluster and the 
# ## posterior samples were exported directly.
# ## The following command was used to extract samples:

# post <- extract.samples(model_fit)

# #######################################################################################################

load("./post_female_final")


## We can visually inspect the individual-level random effects by creating a new object, v_est.
## The pairs function in the base R package allows us to visualize the pairwise correlation of the
## random effects. Based on the analysis of Koster/McElreath, a logical thing to do here would be
## to examine the model output using the precis function and then examine those with strong correlations. 
## Subsequently, the dens function allows us to see the posterior distribution of the correlation
## between those categories.

# female tsimane activities in order they appear in columns of posterior
female_activs <- c("CHOP","CLEAR","EAT","FWC","light_processing","MANU","NW","RICE","SW","WR","YUCCA","z")  
v_est <- apply(post_f$v_id,2:3,median)
pairs(v_est)
plot(v_est[,5],v_est[,7])
dens(post_f$Rho_id[,5,7])

# mean and SD of female ages for z-score conversion
sd_age_female <- 14.89
mean_age_female <- 32.66

#######################################################################################################
## Multinomial version of link function for the Tsimane (female) model
## Much like the link function in the rethinking package, the following function can be used to
## generate predictions for a customized data frame. Conventionally, researchers allow one of the
## predictors to vary while holding other predictors at a constant value (often the sample mean).
## The script begins by defining a sequence length that should correspond to the dimensions of
## data frame that is created to generate the predictions. We use 100 in this case, but other
## values are possible.
## We define additional quantities at the beginning of the function, and the function derives
## the values from the posterior. For example, K is the number of response categories.
## The quantity, ns, is the number of samples.
## The function relies on an intermediate step, the creation of ptemp, which is tied to the likelihood
## function. Essentially, the function works by taking the supplied values for each of the i rows
## in the newly created data frame, then multipying those values by the value of the parameter in each
## sample from the posterior. In other words, for each of the samples from the posterior, the 
## parameters vary, which are used to generate predicted probabilities for each of the supplied
## combinations of values in the new data frame that we will create. We can subsequently 
## use the distribution of predicted values to understand the effets of the values we supply.
## It is common for researchers to average over the random effects. In other words, predicted
## probabilies are based only on the fixed effects, not accounting for the random effects.
## To accomplish that, our function requires users to supply a 0 for the random effects. By contrast,
## one can include random effects by supplying the integer that corresponds to the unit of interest
## in the higher-level classification. The varying intercept for that unit will then be incorporated
## as part of the predicted probability.
#######################################################################################################

seq.length <- length(seq(15, 75, 1))

link.mn <- function( data ) {
  K <- dim(post_f$v_id)[3] + 1
  ns <- dim(post_f$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- seq.length
  
  softmax2 <- function(x) {
    x <- max(x) - x
    exp(-x)/sum(exp(-x))
  }
  
  p <- list()
  
  for ( i in 1:n ) {
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        ptemp <- post_f$a[,k] + 
          post_f$bA[,k] * data$age_z[i] + 
          post_f$bQ[,k] * data$age_zq[i] + 
          post_f$bT[,k] * data$time_z[i] + 
          post_f$bTQ[,k] * data$time_zq[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post_f$v_id[,data$id[i],k]
        if ( data$com_id[i]>0 ) ptemp <- ptemp + post_f$v_com[,data$com_id[i],k]
        if ( data$month_id[i]>0 ) ptemp <- ptemp + post_f$v_month[,data$month_id[i],k]
      } else {
        ptemp <- rep(0,ns)
      }
      return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}


## The following generates a vector of ages which we use to generate predicted probabilities.
## We use the seq function to generate a sequence across the range of ages in the empirical data.
## Note that the length of this vector is seq.length, as defined above.
## Set the target age range to generate 100 points between ages 15 and 75
low_age <- (15 - mean_age_female)/sd_age_female 
high_age <- (75 - mean_age_female)/sd_age_female

age_seq <- seq (from= low_age, to= high_age, length.out = seq.length)

## We create a data frame from age_seq and its second order polynomial, holding most of predictors
## at their sample mean. In the original Koster+McElreath paper, they specified 8:00 am for the time
## of day in order to effect a more sensible alignment with the empirical data. In this case, we 
## hold that parameter at 0, which therefore corresponds with the mean time of day in our dataset
## for behavioral observations (0.532 = ~12:45).
## Also, as we noted earlier in the notes about the multinomial link function, we set "id" to 
## zero in order to average over the random effects.
pred_dat <- data.frame(
  id = 0 ,
  age_z = age_seq,
  age_zq = age_seq^2,
  time_z = 0,
  time_zq = 0,
  com_id = 0,
  month_id = 0
)

## The following is an object with 100 predicted probabilities (as in seq.length) corresponding
## to the different values of standardized ages that were supplied in the pred_dat data frame,
## as repeated for each of the K response categories.
p <- link.mn ( pred_dat )

## The following calculates the mean of the predicted samples.
p_mean <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )

## Alternative way to do this. I find the syntax more intuitive than the alternative, so could help others see what is done.
# do.call(cbind, lapply(p, function(i) apply(i, 2, mean))) 

###############################
# Calculating predicted probabilities at a single time of day cannot tell
# us about how much time an individual of a given age/sex/whatever class 
# is expected to spend on a particular behavior over the course of a "normal"
# day. This is particularly the case if you had more complex functional forms
# for the time parameter, and behaviors that are bi- or even tri-modal over the
# course of a day.

# The following code implements a procedure that at each draw from the 
# posterior integrates across the time of day from 7 am to 7 pm (this is
# when our samples were collected, and when the behaviors of interest occur).
# Thus, the output of this new link function is a composite "total time" spent per day
# across the range specified (could be whole day if measures spanned it).

# Nothing fancy is done here, but an integration is performed within the sapply() call 
# across the day and night interval specified in the function call, using the fitted model and
# our input data at which we want to predict.

# note that this has 1 less column (number of activity categories) than male version
link.day.hours <- function( data, post, morning, night ) {
  K <- dim(post_f$v_id)[3] + 1
  ns <- dim(post_f$v_id)[1]
  if ( missing(data) ) stop( "BOOM: Need data argument" )
  n <- length(age_seq)
  
  p <- list()
  
  for ( i in 1:n ) {
    
    p[[i]] <- sapply( 1:K , function(k) {
      if ( k < K ) {
        sapply( 1:ns, function(j) {    
          integrate(
            function(x) 
              exp(post_f$a[j,k] + post_f$bA[j,k]*data$age_z[i] + post_f$bQ[j,k]*data$age_zq[i] + post_f$bT[j,k]*x+ post_f$bTQ[j,k]*x^2)/
              ( 1 + (
                exp(post_f$a[j,1]+ post_f$bA[j,1]*data$age_z[i] + post_f$bQ[j,1]*data$age_zq[i]+ post_f$bT[j,1]*x+ post_f$bTQ[j,1]*x^2) +
                  exp(post_f$a[j,2]+ post_f$bA[j,2]*data$age_z[i] + post_f$bQ[j,2]*data$age_zq[i]+ post_f$bT[j,2]*x+ post_f$bTQ[j,2]*x^2) +
                  exp(post_f$a[j,3]+ post_f$bA[j,3]*data$age_z[i] + post_f$bQ[j,3]*data$age_zq[i]+ post_f$bT[j,3]*x+ post_f$bTQ[j,3]*x^2) +
                  exp(post_f$a[j,4]+ post_f$bA[j,4]*data$age_z[i] + post_f$bQ[j,4]*data$age_zq[i]+ post_f$bT[j,4]*x+ post_f$bTQ[j,4]*x^2) +
                  exp(post_f$a[j,5]+ post_f$bA[j,5]*data$age_z[i] + post_f$bQ[j,5]*data$age_zq[i]+ post_f$bT[j,5]*x+ post_f$bTQ[j,5]*x^2) +
                  exp(post_f$a[j,6]+ post_f$bA[j,6]*data$age_z[i] + post_f$bQ[j,6]*data$age_zq[i]+ post_f$bT[j,6]*x+ post_f$bTQ[j,6]*x^2) +
                  exp(post_f$a[j,7]+ post_f$bA[j,7]*data$age_z[i] + post_f$bQ[j,7]*data$age_zq[i]+ post_f$bT[j,7]*x+ post_f$bTQ[j,7]*x^2) +
                  exp(post_f$a[j,8]+ post_f$bA[j,8]*data$age_z[i] + post_f$bQ[j,8]*data$age_zq[i]+ post_f$bT[j,8]*x+ post_f$bTQ[j,8]*x^2) +
                  exp(post_f$a[j,9]+ post_f$bA[j,9]*data$age_z[i] + post_f$bQ[j,9]*data$age_zq[i]+ post_f$bT[j,9]*x+ post_f$bTQ[j,9]*x^2) +
                  exp(post_f$a[j,10]+ post_f$bA[j,10]*data$age_z[i] + post_f$bQ[j,10]*data$age_zq[i]+ post_f$bT[j,10]*x+ post_f$bTQ[j,10]*x^2) +
                  exp(post_f$a[j,11]+ post_f$bA[j,11]*data$age_z[i] + post_f$bQ[j,11]*data$age_zq[i]+ post_f$bT[j,11]*x+ post_f$bTQ[j,11]*x^2) 
              )),lower=morning, upper=night)$value * (12/(night-morning))
        })
        
        # could add back in if needed, but not needed unless predicting at specific level of random effect
        # if ( data$id[i]>0 ) ptemp <- ptemp + post_f$v_id[,data$id[i],k]
        # if ( data$com_id[i]>0 ) ptemp <- ptemp + post_f$v_com[,data$com_id[i],k]
        # if ( data$month_id[i]>0 ) ptemp <- ptemp + post_f$v_month[,data$month_id[i],k]
      } else {
        sapply( 1:ns, function(j) { 
          integrate(
            function(x) 1/( 1 + (
              exp(post_f$a[j,1]+ post_f$bA[j,1]*data$age_z[i] + post_f$bQ[j,1]*data$age_zq[i]+ post_f$bT[j,1]*x+ post_f$bTQ[j,1]*x^2) +
                exp(post_f$a[j,2]+ post_f$bA[j,2]*data$age_z[i] + post_f$bQ[j,2]*data$age_zq[i]+ post_f$bT[j,2]*x+ post_f$bTQ[j,2]*x^2) +
                exp(post_f$a[j,3]+ post_f$bA[j,3]*data$age_z[i] + post_f$bQ[j,3]*data$age_zq[i]+ post_f$bT[j,3]*x+ post_f$bTQ[j,3]*x^2) +
                exp(post_f$a[j,4]+ post_f$bA[j,4]*data$age_z[i] + post_f$bQ[j,4]*data$age_zq[i]+ post_f$bT[j,4]*x+ post_f$bTQ[j,4]*x^2) +
                exp(post_f$a[j,5]+ post_f$bA[j,5]*data$age_z[i] + post_f$bQ[j,5]*data$age_zq[i]+ post_f$bT[j,5]*x+ post_f$bTQ[j,5]*x^2) +
                exp(post_f$a[j,6]+ post_f$bA[j,6]*data$age_z[i] + post_f$bQ[j,6]*data$age_zq[i]+ post_f$bT[j,6]*x+ post_f$bTQ[j,6]*x^2) +
                exp(post_f$a[j,7]+ post_f$bA[j,7]*data$age_z[i] + post_f$bQ[j,7]*data$age_zq[i]+ post_f$bT[j,7]*x+ post_f$bTQ[j,7]*x^2) +
                exp(post_f$a[j,8]+ post_f$bA[j,8]*data$age_z[i] + post_f$bQ[j,8]*data$age_zq[i]+ post_f$bT[j,8]*x+ post_f$bTQ[j,8]*x^2) +
                exp(post_f$a[j,9]+ post_f$bA[j,9]*data$age_z[i] + post_f$bQ[j,9]*data$age_zq[i]+ post_f$bT[j,9]*x+ post_f$bTQ[j,9]*x^2) +
                exp(post_f$a[j,10]+ post_f$bA[j,10]*data$age_z[i] + post_f$bQ[j,10]*data$age_zq[i]+ post_f$bT[j,10]*x+ post_f$bTQ[j,10]*x^2) +
                exp(post_f$a[j,11]+ post_f$bA[j,11]*data$age_z[i] + post_f$bQ[j,11]*data$age_zq[i]+ post_f$bT[j,11]*x+ post_f$bTQ[j,11]*x^2) 
            )),lower=morning, upper=night)$value * (12/(night-morning))
        })
      }
      #return(ptemp)
    })
    ## The values are converted to probabilities using the softmax function
    ## which ensures that the predicted values across categories sum to
    ## 100% probabilities.
    #for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
  }
  return(p)
}


mean_time_prop <- 0.5315371
sd_time_prop <- 0.1467263

p.hrs.female <- link.day.hours(data=pred_dat, post=post, morning = ((7/24)-mean_time_prop)/sd_time_prop, night = ((19/24)-mean_time_prop)/sd_time_prop)



saveRDS(p.hrs.female, file="female_hours_final.Rds")


rm(list=ls())