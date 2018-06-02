# A simple Metropolis-Hastings MCMC in R
# https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/

## (1) setting true values
# Creating test data
trueA = 5
trueB = 0
trueSd = 10 
sampleSize = 312

# around x=0
x = (-(sampleSize-1)/2):((sampleSize-1)/2)
# linear model with normally distributed errors mean as 0, sd as trueSd
y = trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

# test code for 'true plot'
# plot(x,y, main="Test Data")

## (2) likelihood functions
likelihood = function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  pred = a*x + b 
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = TRUE)
  sumll = sum(singlelikelihoods)
  sumll
}
slopevalues = function(x) {return(likelihood(c(x, trueB, trueSd)))}

# test code
slopelikelihoods <- lapply(seq(3, 7, by=.05), slopevalues )
# plot (seq(3, 7, by=.05), slopelikelihoods, type="l", 
#       xlab = "values of slope parameter a", ylab = "Log likelihood")


## (3) Priors and Posteriors
# Prior distribution
prior <- function(param){
    a = param[1]
    b = param[2]
    sd = param[3]
    aprior = dunif(a, min=0, max=10, log = T)
    bprior = dnorm(b, sd = 5, log = T)
    sdprior = dunif(sd, min=0, max=30, log = T)
    return(aprior+bprior+sdprior)
}

posterior <- function(param){
   return (likelihood(param) + prior(param))
}


## (4) Proposal function and actual simulation
proposalfunction <- function(param){
    return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
}
 
run_metropolis_MCMC <- function(startvalue, iterations){
    chain = array(dim = c(iterations+1,3))
    chain[1,] = startvalue
    for (i in 1:iterations){
        proposal = proposalfunction(chain[i,])
         
        probab = exp(posterior(proposal) - posterior(chain[i,]))
        if (runif(1) < probab){
            chain[i+1,] = proposal
        }else{
            chain[i+1,] = chain[i,]
        }
    }
    return(chain)
}
 
# initial point
# we can set multiple initial points and see if they converge
startvalue = c(4,0,10) 
chain = run_metropolis_MCMC(startvalue, 20000)
 
burnIn = 10000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))


par(mfrow = c(2,3))
hist(chain[-(1:burnIn),1],nclass=30, , main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = trueA, col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = trueB, col="red" )
hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]) )
abline(v = trueSd, col="red" )
plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = trueA, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = trueB, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSd, col="red" )
 
# for comparison:
summary(lm(y~x))

