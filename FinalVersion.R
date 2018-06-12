library(MASS)
library(coda)
library(MCMCpack)
library(BayesianTools)

set.seed(54321)

true_theta<-c(0.1,0.7,0.2)
alpha<-c(2,2,2)
x<-rmultinom(20,size=7,prob = true_theta)   # FIX THIS

gibbs_func<-function(start_value,burnin,thin,iter){
  update<-start_value
  theta_mat<-mat.or.vec(2,iter)
  for(i in 1:iter){
    theta1<-(1-update[2])*(rbeta(1,shape1 = alpha[1]+sum(x[1,]),shape2 = alpha[3]+sum(x[3,])))
    theta2<-(1-theta1)*(rbeta(1,shape1 = alpha[2]+sum(x[2,]),shape2 = alpha[3]+sum(x[3,])))
    update<-c(theta1,theta2)
    theta_mat[,i]<-update
  }
  sample1<-theta_mat[1,][seq(burnin+1,iter,thin)]
  sample2<-theta_mat[2,][seq(burnin+1,iter,thin)]
  
  chain<-matrix(nrow=length(sample1), ncol=2)
  chain[,1]<-sample1
  chain[,2]<-sample2
  
  return(chain)
  
}

# Posterior
posterior_func <- function(param){
  single_likelihood <- apply(x, 2, function(t) dmultinom(t, prob=param))
  single_prior <- ddirichlet(param, alpha)
  nlikelihood <- sum(log(single_likelihood))
  nprior <- sum(log(single_prior))
  return(nlikelihood + nprior)
  #posterior = sum(apply(post_alpha, 2, function(t) ddirichlet(param, t)))
  # posterior = ddirichlet(param, post_alpha)
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

start_value = apply(x, 1, sum) / 140

start = Sys.time()
gibbs_chain<-gibbs_func(start_value[1:2],1000,100,100000)
cat('gibbs time spent: ', Sys.time() - start, '\n')

start = Sys.time()
metro_chain = run_metropolis_MCMC(start_value, 1000,100,100000)
cat("metropolis time spent: ", Sys.time() - start, '\n')

# analysis
summary(gibbs_chain)
par(mar=c(2,2,2,2))

# jpeg('gibbs_chain.jpg')
plot(mcmc(gibbs_chain))
# dev.off()

# jpeg('metropolis_chain.jpg')
plot(mcmc(metro_chain))
# dev.off()

# # check partial correlatation btw parameters
# correlationPlot(chain)
# 
# # check convergence
# chain2<-gibbs(c(sum(x[1,])/140,sum(x[2,])/140),1000,100,1000000)
combinedchains = mcmc.list(mcmc(gibbs_chain), mcmc(metro_chain))
# jpeg('combinedChains.jpg')
plot(combinedchains)
# dev.off()

# jpeg('diagnostic.jpg')
# gelman.diag(combinedchains)
# dev.off()

# jpeg("gelman.jpg")
gelman.plot(combinedchains)
# dev.off()
