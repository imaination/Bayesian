### MH algorithm w/ analysis 
### June 2.
### Authors: Ashly, Bella, Jack
### MetropolisHastrings.R
#
# This script is simulating numbers for a linear regression model.
# We have three parameters, slopt, intercept, and standard deviation
# of error. 
# 
library(coda)

trueA = 5
trueB = 0
trueSd = 10
sampleSize = 31

x = (-(sampleSize-1)/2):((sampleSize-1)/2)
y =  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

likelihood = function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  
  pred = a*x + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

prior = function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  # NOTE:
  # a and b can have the same distribution (could be normal)
  # sigma2 of error could have inv-chisq distribution
  aprior = dunif(a, min=0, max=10, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+sdprior)
}

proposalfunction = function(param){
  return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(likelihood(proposal)+ prior(proposal) - likelihood(chain[i,])- prior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue = c(4,2,8)
chain = run_metropolis_MCMC(startvalue, 10000)
result = mcmc(chain)

# analysis
summary(result)
plot(result)

# check partial correlatation btw parameters
library(BayesianTools)
correlationPlot(chain)

# check convergence 
chain2 = run_metropolis_MCMC(startvalue, 10000)
result2 = mcmc(chain2)
combinedchains = mcmc.list(result, result2)
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)
