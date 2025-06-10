# This file takes the output of the multinomial model fit in Stan to Hadza scan sampling, and further processes it before application.

# required packages
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


# ##############################################################################
# ## The first step toward plotting predicted probabilities of our fixed effects is to extract the
# ## posterior samples using the extract.samples function in the STAN package.
# ## Because of the memory requirements for this operation, it was run on the HPC and the 
# ## posterior samples were exported directly.
# ###############################################################################

#############################
####### FOR MALES ##########
############################
load("./a_male_hadza_final_5_29_2020.RData") # loads "a_male" into environment (this is the dataframe used in the multinomial model of behavior)
load("./post_male_hadza_final_5_29_2020.RData")  # loads "post_male" into environment (this holds posterior samples extracted from fitted multinomial model of behavior using >> extract.samples(fitted_mod)). 

# rename inputs
post <- post_male
a <- a_male

# mean and SD of male ages for z-score conversion
sd_age_male <- 18.53
mean_age_male <- 35.78

# sequence of age inputs within range of model: (15, 75) years old
seq.length <- length(seq(15, 75, 1))

# Define link function with softmax transformation
# The following function is taken from Koster and McElreath (2017)
# Modifications made to accommodate current model structure (e.g. addition of com_id, month)
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

# Calculate predictions for average time spent on activities over the whole day.
low_age <- (15 - mean_age_male)/sd_age_male 
high_age <- (75 - mean_age_male)/sd_age_male 

age_seq <- (seq(15, 75, 1)- mean_age_male)/sd_age_male 

pred_dat <- data.frame(
  id = 0 ,
  age_z = age_seq,
  age_zq = age_seq^2,
  time_z = 0,
  time_zq = 0,
  com_id = 0,
  month_id = 0
)


# Generate predictions using link function
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
                  exp(post$a[j,6]+ post$bA[j,6]*data$age_z[i] + post$bQ[j,6]*data$age_zq[i]+ post$bT[j,6]*x+ post$bTQ[j,6]*x^2)
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
                exp(post$a[j,6]+ post$bA[j,6]*data$age_z[i] + post$bQ[j,6]*data$age_zq[i]+ post$bT[j,6]*x+ post$bTQ[j,6]*x^2)
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

mean_time_prop <- 0.5338086
sd_time_prop <- 0.1469649

p.hrs.male <- link.day.hours(pred_dat, post=post, morning = ((7/24)-mean_time_prop)/sd_time_prop, night = ((19/24)-mean_time_prop)/sd_time_prop)


save(p.hrs.male, file="hadza_male_hours_final_05_29_2020.RData")


###########################################################################################
####### FOR FEMALES ########################################################################
#############################################################################################

load("/Users/thomaskraft/Dropbox/energetics_tsimane_hadza/data/a_female_hadza_final_05_29_2020.RData")
load("/Users/thomaskraft/Dropbox/energetics_tsimane_hadza/data/post_female_hadza_final_05_29_2020.RData")

# mean and SD of female ages for z-score conversion
sd_age_female <- 16.40
mean_age_female <- 37.31

# activity vector (order is relevant, keep as is): these correspond to the levels in the posterior
activs <- c("carry_firewood_water", "eat", "light_processing", "manufacture", "z", "out_of_camp", "pounding_baobab")   

seq.length <- length(seq(15,75, 1))

link.mn <- function( data ) {
  K <- dim(post_female$v_id)[3] + 1
  ns <- dim(post_female$v_id)[1]
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
        ptemp <- post_female$a[,k] + 
          post_female$bA[,k] * data$age_z[i] + 
          post_female$bQ[,k] * data$age_zq[i] + 
          post_female$bT[,k] * data$time_z[i] + 
          post_female$bTQ[,k] * data$time_zq[i]
        if ( data$id[i]>0 ) ptemp <- ptemp + post_female$v_id[,data$id[i],k]
        if ( data$com_id[i]>0 ) ptemp <- ptemp + post_female$v_com[,data$com_id[i],k]
        if ( data$month_id[i]>0 ) ptemp <- ptemp + post_female$v_month[,data$month_id[i],k]
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


low_age <- (15 - mean_age_female)/sd_age_female
high_age <- (75 - mean_age_female)/sd_age_female

# everything will be calculated over the "adult" interval from 15 to 75
age_seq <- (seq(15, 75, 1)- mean_age_female)/sd_age_female


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


p <- link.mn ( pred_dat )

## The following calculates the mean of the predicted samples.
p_mean <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )


pred_dat <- data.frame(
  id = 0 ,
  age_z = age_seq,
  age_zq = age_seq^2,
  time_z = 0,
  time_zq = 0,
  com_id = 0,
  month_id = 0
)

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
                  exp(post$a[j,6]+ post$bA[j,6]*data$age_z[i] + post$bQ[j,6]*data$age_zq[i]+ post$bT[j,6]*x+ post$bTQ[j,6]*x^2)
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
                exp(post$a[j,6]+ post$bA[j,6]*data$age_z[i] + post$bQ[j,6]*data$age_zq[i]+ post$bT[j,6]*x+ post$bTQ[j,6]*x^2)
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

p.hrs.female <- link.day.hours(pred_dat, post=post_female, morning = ((7/24)-mean_time_prop)/sd_time_prop, night = ((19/24)-mean_time_prop)/sd_time_prop)

save(p.hrs.female, file="hadza_female_hours_final_05_29_2020.RData")
