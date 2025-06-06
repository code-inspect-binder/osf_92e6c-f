## This script modifies the multinomial multilevel behavior model from Koster 
## and McElreath (2017) "The Multinomial Analysis of Behavior," from their 
## supplemental file. It is modified to work with the energetics and time 
## allocation data for Hadza of Tanzania 

## Note that full credit for the development of this model lies with Jeremy 
## Koster and Richard McElreath, and not the authors of the present paper.

## This script uses functions from several R packages and their dependencies:
library(tidyr)
library(dplyr)
library(rstan)
library(rethinking)
library(chron)


## Prepare the data
source(file = "scan_prepare.R") #loads data and makes ready for analysis

## each row of data file "a" contains a behavioral observation that corresponds with an energetic category in addition to information on date, time of day, community, and focal attributes such as age/sex

# Model for either male or female
#a <- hm  # male
a <- hf  #female

# reduce to adults only
a <- a[which(a$age >= 15), ]

# make response variable a factor
a$work <- as.factor(a$work)



##The following are notes from the original authors that remain relevant here:

## Later in this script, we will use the coerce_index function to turn our
## response variable, currently a factor, into an integer that serves as the
## outcome variable(s) in the STAN model. By changing the "reference" level to "z",
## we ensure that the reference level is last in alphabetical order and therefore
## becomes the reference category by default in integer format. There are 
## many possible alternatives to this approach that would equally ensure that
## the desired reference level is the last among K categories in the response.
levels(a$work)[levels(a$work)=="nonsubsistence"] <- "z"

## Time as record was in 00:00 format whereas the "times" function in the 
## chron package requires 00:00:00 format. The following nested code
## first appends the seconds to the times, then converts them into times,
## then converts them to a vector of numbers that specific times as the
## proportion of a 24-hour day that has elapsed at the observed time. For
## example, noon yields a value of 0.5 whereas an observation as 6:00 PM
## yields a value of 0.75.

a$time <- strftime(a$date_time, format="%H:%M:%S")
a$time_prop <- as.numeric(times(a$time))

# need factor formats for coerce_index below
a$id <- as.factor(a$id)
a$camp <- as.factor(a$camp)


##The following are notes from the original authors that remain relevant here:

## What follows is the preparation of data needed in the STAN model. In general,
## STAN requires categorical variables (including for higher-level random effects)
## to be formatted as sequential integers. It therefore lacks the convenience of
## packages such as lme4, which makes the necessary conversions of categorical 
## variables without effort by the applied researcher. From the rethinking package,
## the coerce_index function is designed to convert vectors into the integer format
## required by STAN.
## In general, we employ z-score transformation of our continuous variables. This
## is primarily to make the sampling more efficient. Note that standardization 
## is based on sample means and standard deviations, where the sample is all 
## observations in the dataset.
## Transformed data that correspond to the observational data are appended to
## the "d' data frame. We also create index variables (N, K, N_id, N_month)
## that will be necessary for the STAN model code. These latter variables index
## the number of observations, the number of response categories, and the number
## of units in the higher-level classifications (for the random effects).
N <- nrow(a) ## Number of observations
a$y <- coerce_index (a$work) ## Renaming response variable
K <- max(a$y) ## Number of response categories
a$id <- coerce_index (a$id) ## Index of observed individuals
N_id <- max(a$id) ## Number of observed individuals
a$com_id <- coerce_index (a$camp) ## Index of observed individuals' communities
N_com <- max(a$com_id) ## Number of communities
a$month_id <- coerce_index (a$month) ## Index of months in which observations occurred
N_month <- max(a$month) ## Number of months (12)
a$age_z <- (a$age-mean(a$age))/sd(a$age) ## Standardized age of observed individuals
a$age_zq <-a$age_z^2 ## Quadratic transformation of standardized age
a$time_z <- (a$time_prop - mean(a$time_prop))/sd(a$time_prop) ## Standardized time of day
a$time_zq <- a$time_z^2 ## Quadratic transformation of standardized time of day

#######################################################################################################
## Data lists
## In the call to STAN, there is a "data' argument, which calls for a list of data to be used.
## We find it helpful to prepare lists in advance of the call to STAN, and what follows is
## data lists for each of the four models that we implement in this script.
## Note that we use suffixes to distinguish our models:
## (i) is the suffix for models with only individual-level random effects
## (ihm) is the suffix for models with random effects for individuals (i), household (h), and month (m)
## (iF) has individual-level random effects and the fixed effects
## (ihmF) has all three aforementioned random effects and the fixed effects

## For energetics project:
## (icmF) has random effects for individuals (i), community (c), and month (m), and all fixed effects

## Fixed effects include age, age squared, time of day, and time of day squared

#######################################################################################################

dat_list_icmF <- list(
  K = K,
  N = N,
  N_id = N_id,
  N_com = N_com,
  N_month = N_month,
  y = a$y,
  id = a$id,
  com_id = a$com_id,
  month_id = a$month_id,
  age_z = a$age_z,
  age_zq = a$age_zq,
  time_z = a$time_z,
  time_zq = a$time_zq
)

#######################################################################################################
## The call to STAN also requires model code. Using the same suffixes as above, we generate the model
## code for the four models. We begin by stating the data, distinguishing between integers and "real"
## continuous variables. Those variables originating in the "d" data frame are appended with [N] in
## brackets to denote that there is a data point corresponding to each observation in our dataset.
## We note that STAN model code can be somewhat idiosyncratic and the order in which variables are
## defined in a block of code can be somewhat rigid. In the "generated quantities" block, for instance,
## all matrices must be defined prior to the generation of the correlations (i.e., the rhos for the
## random effects).
## Further annotations in the model code follow the dual slashes "//"
#######################################################################################################

model_code_icmF <- "
data{
int N;
int N_id;
int N_com;
int N_month;
int y[N];
int id[N];
int com_id[N];
int month_id[N];
real age_z[N];
real age_zq[N];
real time_z[N];
real time_zq[N];
int K; // number of categories
}
parameters{
real a[K-1];                  // intercepts for each behavior
real bA[K-1];				  // fixed effect for age
real bQ[K-1];				// fixed effect for age squared
real bT[K-1];				// fixed effect for time of day
real bTQ[K-1];				// fixed effect for time of day squared
matrix[K-1,N_id] z_id;      		// matrix of indiv-level standardized random effects
vector<lower=0>[K-1] sigma_id;  	 // stddev of indiv-level random effects
cholesky_factor_corr[K-1] L_Rho_id;  // correlation matrix of indiv-level random effects
matrix[K-1,N_com] z_com;     		 // matrix of community-level standardized random effects
vector<lower=0>[K-1] sigma_com;   // stddev of community-level random effects
cholesky_factor_corr[K-1] L_Rho_com;         // correlation matrix of community-level random effects
matrix[K-1,N_month] z_month;      // matrix of month-level standardized random effects
vector<lower=0>[K-1] sigma_month;   // stddev of month-level random effects
cholesky_factor_corr[K-1] L_Rho_month;         // correlation matrix of month-level random effects   
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_com,K-1] v_com;
matrix[N_month,K-1] v_month;
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_com = (diag_pre_multiply(sigma_com,L_Rho_com) * z_com)';
v_month = (diag_pre_multiply(sigma_month,L_Rho_month) * z_month)';
}
model{
// priors
a ~ normal(0,14);
bA ~ normal(0,1);
bQ ~ normal(0,1);
bT ~ normal(0,1);
bTQ ~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);
to_vector(z_com) ~ normal(0,1);
sigma_com ~ exponential(1);
L_Rho_com ~ lkj_corr_cholesky(2);
to_vector(z_month) ~ normal(0,1);
sigma_month ~ exponential(1);
L_Rho_month ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] <- a[k] + bA[k] * age_z[i] + bQ[k] * age_zq[i] + bT[k] * time_z[i] + bTQ[k] * time_zq[i] + v_id[id[i],k] + v_com[com_id[i],k] + v_month[month_id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_com;
matrix[K-1,K-1] Rho_month;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_com = L_Rho_com * L_Rho_com';
Rho_month = L_Rho_month * L_Rho_month';

for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] <- a[k] + bA[k] * age_z[i] + bQ[k] * age_zq[i] + bT[k] * time_z[i] + bTQ[k] * time_zq[i] + v_id[id[i],k] + v_com[com_id[i],k] + v_month[month_id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

#######################################################################################################
## What follows is the call to STAN for each of the four models. Prior to the call, we define
## (1) starting values for the fixed effects, (2) the variance of the random effects for each of the
## K-1 responses, (3) the variance-covariance matrix of the random effects for the K-1 responses, and
## (4) the unit-by-unit matrix of standardized random effects.
## In each case, we opt for uninformative starting values, which are then indexed by the init object
## that is supplied to the model call.
## We also supply values for the number of chains -- see McElreath's Statistical Rethinking for
## guidance on the number of chains and other considerations in calls to STAN, including the
## adapt_delta argument, which dictates the proportion of the proposed samples that must be within
## acceptable bounds. Raising this number toward 1 has benefits of refinining the estimator during
## warmup, but the tradeoff is that warmup requires more time as a result.
## Note that we retain the use of suffixes to distinguish between models for these auxiliary objects
## and quantities.
#######################################################################################################

################################################
## Model icmF

start_icmF <- list(
  a = rep(0,K-1),
  bA = rep(0,K-1),
  bQ = rep(0,K-1),
  bT = rep(0,K-1),
  bTQ = rep(0,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_com = rep(1,K-1),
  L_Rho_com = diag(K-1),
  z_com = matrix(0,nrow=K-1,ncol=N_com),
  sigma_month = rep(1,K-1),
  L_Rho_month = diag(K-1),
  z_month = matrix(0,nrow=K-1,ncol=N_month)
)

n_chains_icmF <- 8 
init_icmF <- list()
for ( i in 1:n_chains_icmF ) init_icmF[[i]] <- start_icmF

mfit_hadza_female <- stan( model_code=model_code_icmF , data=dat_list_icmF , chains=n_chains_icmF , cores= n_chains_icmF , warmup=500, iter=1000, init=init_icmF , control = list(adapt_delta = 0.95))

save.image("rstan_hadza_female_05_29_2020.RData")


# do some processing on cluster (crashes weak laptop)
# these outputs get used in subsequent analysis file
post_female <- extract.samples(mfit_hadza_female)
a_female <- a
save(post_female, file="post_female_hadza_final_05_29_2020.RData")
save(a_female, file= "a_female_hadza_final_05_29_2020.RData")


#diagnostics 
check_hmc_diagnostics(mfit_hadza_female)