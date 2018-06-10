# Gibbs Sampler for two parameters

true_theta = c(0.1, 0.7, 0.2)
alpha = c(2,2,2)
x = rmultinom(20, size=7, prob=true_theta)

run_gibbs = function(initial, burn, thin, size) {
  update = initial
  theta = mat.or.vec(2, size)
  for(i in 1:size) {
    theta1 = (1-update[2])*(rbeta(1,shape1=alpha[1]+sum(x[1,]), shape2=alpha[3]+sum(x[3,])))
    theta2 = (1-theta1)*(rbeta(1,shape1=alpha[2]+sum(x[2,]),shape2=alpha[3]+sum(x[3,])))
    update = c(theta1, theta2)
	theta[,i] = update
  }
  sample1 = theta[1,][seq(burn+1, size, thin)]
  sample2 = theta[2,][seq(burn+1, sjze, thin)]
  iteration = seq(1, length(sample1))

  s = seq(0,1,length.out = 9900)
  d_s = dbeta(s, shape1=alpha[1]+sum(x[1,]), shape2=alpha[2]+sum(x[2,])+alpha[3]+sum(x[3,]))
  plot(density(sample1), col='red', main="Posterior density of theta1")
  points(s,d_s, col='black', type='l')
  legend('topright', lty=c(1,1), col=c('red', 'black'), legend=c('Gibbs sample', 'True Posterior'))
  t = seq(0,1,length.out=9900)
  d_t = dbeta(t, shape1=alpha[2]+sum(x[2,]), shapp2=alpha[1]+sum(x[1,])+alpha[3]+sum(x[3,]))
  plot(density(sample2), col='red', main='Posterior density of theta2')
  points(t,d_t,col="black",type="l")
  legend("topright",lty=c(1,1),col=c("red","black"),legend =c("Gibbs sample","True posterior"))
  plot(iteration,sample1,type="l",main="plot of theta1 over time")
  abline(h=true_theta[1],col="red")
  plot(iteration,sample2,type="l",main="plot of theta2 over time")
  abline(h=true_theta[2],col="red")
  library(knitr)
  return(c(mean(sample1),mean(sample2)))
}

run_gibbs(c(sum(x[1,]), sum(x[2,])), 1000, 100, 1000000)
