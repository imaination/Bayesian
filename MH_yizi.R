### copy owned by yizi
### June 8.
### Authors: Ashley, Yizi, Jack
### MetropolisHastrings.R
#
# This script is simulating numbers for a linear regression model.
# We have three parameters, slopt, intercept, and standard deviation
# of error. 
# 
# likelihood distribution: normal
# prior distributions: slope: normal, intercept: normal, sd: unif
# proposal distributions: normal

set.seed(54321)
library(coda)
library(BayesianTools)

# real parameters
trueSlope = 5
trueIntercept = 2
trueSd = 10
sampleSize = 422

# range of x
x = (-(sampleSize-1)/2):((sampleSize-1)/2) 
# linear regression model w/ normally distributed error
y = trueSlope * x + trueIntercept + rnorm(n=sampleSize, mean=0, sd=trueSd)

likelihood_func <- function(param){
  # get likelihood function
  slope = param[1]
  intercept = param[2]
  sd = param[3]
  y_hat = slope*x + intercept
  # print(y_hat)
  singlelikelihood = dnorm(y, mean=y_hat, sd=sd, log=T)
  which(is.nan(singlelikelihood))
  # print(singlelikelihood)
  nlikelihood = sum(singlelikelihood)
  #print(nlikelihood)
  nlikelihood
}

prior_func <- function(param){
  # prior 
  slope = param[1]
  intercept = param[2]
  sd = param[3]
  slope_prior = dunif(slope, min=0, max=10, log=T)
  intercept_prior = dnorm(intercept, mean=3, sd=2, log=T)
  sd_prior = dnorm(sd, sd=5, log=T)
  return(slope_prior+intercept_prior+sd_prior)
  
}

posterior_func <- function(param){
  # posterior
  return(likelihood_func(param) + prior_func(param))
}

proposal_func <- function(param){
  # proposal function
  return(rnorm(3, mean=param, sd = c(0.1, 0.4, 0.9)))
}

metropolis_hastings_func <- function(start_value, iter, burnIn){
  # actual function
  chain = array(dim = c(iter+1,3))
  chain[1,] = start_value
  for (i in 1:iter){
    proposal = proposal_func(chain[i,])
    prob = exp(posterior_func(proposal) - posterior_func(chain[i,]))
    if (runif(1) < prob){
      chain[i+1,] = proposal
    }
    else{
      chain[i+1,] = chain[i,]
    }
  }
  chain = chain[-(1:burnIn),]
  return(chain)
}

# set initial point
start_value = c(3, 2, 1)
S = 50000
burnIn = 10000
chain = metropolis_hastings_func(start_value, S, burnIn)
 
# analysis
acceptRatio = 1-mean(duplicated(chain))   
summary(chain)
par(mar=c(2,2,2,2))
plot(mcmc(chain))

# check partial correlatation btw parameters
correlationPlot(chain)

# check convergence 
chain2 = metropolis_hastings_func(start_value, S, burnIn)
combinedchains = mcmc.list(mcmc(chain), mcmc(chain2))
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)

  

