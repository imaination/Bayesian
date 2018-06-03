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

# real parameters
trueSlope = 5
trueIntercept = 2
trueSd = 10
smapleSize = 422

# range of x
x = (-(sampleSize-1)/2):((sampleSize-1)/2) 
# linear regression model w/ normally distributed error
y = trueSlope * x + trueIntercept + rnorm(n=sampleSize, mean=0, sd=trueSD)

likelihood_func <- function(param){
# get likelihood function
	slope = param[1]
	intercept = param[2]
	sd = paramp[3]
	y_hat = slope*x + intercept
	singlelikelihood = dpois(y, lamda = 1)
	nlikelihood = prod(singlelikelihood)
	return nlikelihood
}

prior_func <- function(param){
# prior 
	slope = param[1]
	intercept = param[2]
	sd = param[3]
	slope_prior = dgamma(slope, shape=3, rate=7)
	intercept_prior = dunif(intercept, min=0, max=10)
	sd_prior = dexp(sd, rate=5)
	return(slope_prior*intercept_prior*sd_prior)
	
}

prosterior_func <- function(param){
# posterior
	return(likelihood_func(param)*prior_func(param))
}

proposal_func <- function(param){
# proposal function
	return(rnorm(3, mean=param, sd = c(0.1, 0.4, 0.9)))
}

metropolis_Hastings_func <- function(start_value, iter){
# actual function
	chain = matrix(c(), nrow=iter+1, ncol=3)
	chain[1,] = start_value
	for (i in 1:iter){
		proposal = proposal_func(chain[i,])
		prob = posterior_func(proposal) / posterior(chain[i,])
		if (runif(1) < prob){
			chain[i+1,] = proposal
		}
		else{
			chain[i+1,] = chain[i,]
		}
	}
	return(mcmc(chain))
}

# real parameters
trueSlope = 5
trueIntercept = 2
trueSd = 10
smapleSize = 422

# range of x
x = (-(sampleSize-1)/2):((sampleSize-1)/2) 
# linear regression model w/ normally distributed error
y = trueSlope * x + trueIntercept + rnorm(n=sampleSize, mean=0, sd=trueSD)

# set initial point
start_value = c(3, 2, 1)
chain = metropolis_Hastings_func(start_value, 20000)

burnIn = 10000
acceptRatio = 1-mean(duplicated(chain[-(1:burnIn),]))    ### REDEFINE ?

# analysis
summary(chain)
plot(chain)
