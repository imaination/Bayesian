library(MCMCpack)
library(BayesianTools)


# Data
true_theta <- c(0.1,0.7,0.2)
alpha <- c(2,2,2)
x <- rmultinom(1,size=7,prob = true_theta)    # FIX THIS

# Posterior
posterior_func <- function(param){
  single_likelihood <- dmultinom(x, prob=param)
  single_prior <- ddirichlet(param, alpha)
  nlikelihood <- sum(log(single_likelihood))
  nprior <- sum(log(single_prior))
  return(nlikelihood + nprior)
}

# Proposal
scale_factor = 1000  # Scale variance of proposal function
proposal_func <- function(param, scale_factor){
  param <- rdirichlet(1, param * scale_factor)
  return(param)
}

# Metropolis MCMC algorithm
run_metropolis_MCMC <- function(start_value, burnin, thin, iter, s = scale_factor){
  theta_mat = array(dim = c(iter+1,3))
  theta_mat[1,1:3] = start_value
  for(i in 1:iter){
    proposal = proposal_func(theta_mat[i,1:3], s)
    prob = exp(posterior_func(proposal) - posterior_func(theta_mat[i,1:3]))
    if (runif(1) < prob){
      theta_mat[i+1,1:3] = proposal
    }else{
      theta_mat[i+1,1:3] = theta_mat[i,1:3]
    }
    
  }
  sample1<-theta_mat[,1][seq(burnin+1,iter,thin)]
  sample2<-theta_mat[,2][seq(burnin+1,iter,thin)]
  
  chain<-matrix(nrow=length(sample1), ncol=2)
  chain[,1]<-sample1
  chain[,2]<-sample2
  
  return(chain) 
}

start_value = c(0.2, 0.4, 0.4)
chain = run_metropolis_MCMC(start_value, 1000,100,100000)

# analysis
summary(chain)
par(mar=c(2,2,2,2))
plot(mcmc(chain))

# check partial correlatation btw parameters
correlationPlot(chain)

# check convergence 
chain2<-run_metropolis_MCMC(c(0.1,0.2,0.7), 1000, 100,100000)
combinedchains = mcmc.list(mcmc(chain), mcmc(chain2))
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)
