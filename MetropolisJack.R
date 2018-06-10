# Metropolis Sampler for two parameters
library(MASS)


# set true parameters
true_theta = c(0.1, 0.7, 0.2)
alpha = c(2,2,2)
x = rmultinom(20, size=7, prob=true_theta) 
init = apply(x,1,sum)/ sum(x)

posteriorGenerator = function(param) {
# likelihood + prior
  likelihood = sum(log(param^x)) # '^' is an element-wised operation here
  prior = sum(log(param^alpha)) # dirichlet
  likelihood + prior 
}

proposalGenerator = function(param){
# Using bivariate normal distribution 
# as our jumping distribution
  mvrnorm(n=1, param, matrix(c(1,0,0,1), 2, 2))

}

run_metropolis = function(initial, burn, thin, size) {
  chain = array(dim = c(size+1, 4))
  chain[1,1:3] = initial
  for (i in 1:size) {
    proposal = proposalGenerator(chain[i,1:2])
    logratio = posteriorGenerator(c(proposal, 1-proposal)) - posteriorGenerator(chain[i,1:3]) # log
	ratio = exp(logratio)
    if(runif < ratio) {
	  chain[i+1,1:3] = c(proposal, 1-sum(c(proposal,1-proposal)))
	  chain[i+1,4] = posteriorGenerator(proposal)
	} else{
	  chain[i+1,1:3] = chain[i, 1:3]
	  chain[i+1,4] = posteriorGenerator(chain[i,1:3])
	}
  }
  return(chain)
}

test=run_metropolis(init, 1000,100,1000000)


