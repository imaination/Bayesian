library(coda)
library(MCMCpack)
library(BayesianTools)

true_theta<-c(0.1,0.7,0.2)
param<-c(2,2,2)
x<-rmultinom(20,size=7,prob = true_theta)

gibbs_func<-function(start_value,burnin,thin,iter){
  update<-start_value
  theta_mat<-mat.or.vec(2,iter)
  for(i in 1:iter){
    theta1<-(1-update[2])*(rbeta(1,shape1 = param[1]+sum(x[1,]),shape2 = param[3]+sum(x[3,])))
    theta2<-(1-theta1)*(rbeta(1,shape1 = param[2]+sum(x[2,]),shape2 = param[3]+sum(x[3,])))
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

chain<-gibbs_func(c(sum(x[1,])/140,sum(x[2,])/140),1000,100,1000000)

# analysis
summary(chain)
par(mar=c(2,2,2,2))
plot(mcmc(chain))

# check partial correlatation btw parameters
correlationPlot(chain)

# check convergence 
chain2<-gibbs(c(sum(x[1,])/140,sum(x[2,])/140),1000,100,1000000)
combinedchains = mcmc.list(mcmc(chain), mcmc(chain2))
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)
