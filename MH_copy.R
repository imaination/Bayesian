### copy owned by yizi
### June 2.
### Authors: Ashly, Bella, Jack
### MetropolisHastrings.R
#
# This script is simulating numbers for a linear regression model.
# We have three parameters, slopt, intercept, and standard deviation
# of error. 
# 
# likelihood distribution: poisson distirbution, lambda=1
# prior distributions: slope: gamma(3,7), intercept: unif, sd: exponential(5)
# proposal distributions: Normal((1,2,3), (1,4,9))

library("coda")

likelihood_func <- function(param, x, y){
# get likelihood function
	slope = param[1]
	intercept = param[2]
	sd = param[3]
	y_hat = slope*x + intercept
	# print(y_hat)
	singlelikelihood = dnorm(y, mean=y_hat, sd=sd, log=T)
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
	slope_prior = dgamma(slope, shape=3, rate=7, log=T)
	intercept_prior = dunif(intercept, min=0, max=10, log=T)
	sd_prior = dexp(sd, rate=5, log=T)
	return(slope_prior+intercept_prior+sd_prior)
	
}

posterior_func <- function(param, x, y){
# posterior
  l = likelihood_func(param, x, y)
  p = prior_func(param)
	return(l + p)
}

proposal_func <- function(param){
# proposal function
	return(rnorm(3, mean=param, sd = c(0.1, 0.4, 0.9)))
}

metropolis_Hastings_func <- function(start_value, iter, x, y){
# actual function
	chain = array(dim = c(iter+1,3))
	chain[1,] = start_value
	for (i in 1:iter){
		proposal = proposal_func(chain[i,])
		prob = exp(posterior_func(proposal, x, y) - posterior_func(chain[i,], x, y))
		if (runif(1) < prob){
			chain[i+1,] = proposal
		}
		else{
			chain[i+1,] = chain[i,]
		}
	}
	result = mcmc(chain)
	return(result)
}

# real parameters
trueSlope = 5
trueIntercept = 2
trueSd = 10
sampleSize = 422

# range of x
x = (-(sampleSize-1)/2):((sampleSize-1)/2) 
# linear regression model w/ normally distributed error
y = trueSlope * x + trueIntercept + rnorm(n=sampleSize, mean=0, sd=trueSd)

# set initial point
start_value = c(3, 2, 1)
chain = metropolis_Hastings_func(start_value, 20000, x, y)

burnIn = 10000
acceptRatio = 1-mean(duplicated(chain[-(1:burnIn),]))    ### REDEFINE ?

# analysis
summary(chain)
plot(mcmc(chain))
