# Gibbs sampler from:http://hedibert.org/wp-content/uploads/2013/12/lab1-realtime.R.txt
# Classical and Bayesian approaches to the 
# standard normal linear regression model

# simulating data from a simple linear regression
n = 100
e = rnorm(n)
x = runif(n)
y = 2+5*x+e

# scatterplot of x versus y
plot(x,y)

# How do I run a simple OLS regression R?
?lm

# simple linear regression via OLS
lm(y~x)


# Save the linear regression results
lab1reg = lm(y~x)

names(lab1reg)

lab1reg$coefficients

lab1reg$residuals

# Plot the residuals
plot(lab1reg$residuals)
abline(h=0,lty=2)

# Plot fitted line
plot(x,y,pch=16)
lines(x,lab1reg$fitted,col=2,lwd=4)

# Create design matrix
X = cbind(1,x)

# Computing estimate 
bhat = solve(t(X)%*%X)%*%t(X)%*%y

Se = t(y-X%*%bhat)%*%(y-X%*%bhat)

#Alternative way of computing Se
error = y-X%*%bhat
Se = sum(error^2)


# Drawing the distribution of betahat2
iXX = solve(t(X)%*%X)

V = Se[1,1]*iXX

sdbeta2 = sqrt(V[2,2])

bs = seq(bhat[2]-5*sdbeta2,
         bhat[2]+5*sdbeta2,
         length=1000)

plot(bs,dnorm(bs,bhat[2],sdbeta2))

# Now performing Bayesian inference
# via Gibbs sampler
# ---------------------------------

# Prior hyperparameters
b0    = matrix(0,2,1)
B0    = diag(100,2)
n0    = 6
sig20 = 1

# initial value for beta
b = matrix(0,2,1)

# Gibbs sampler in action
# -----------------------
set.seed(2935957)
niter = 10000
draws = matrix(0,niter,3)
for (iter in 1:niter){
  # Sampling from full conditional of sigma2
  part1   = (n0+n)/2
  sumres2 = sum((y-X%*%b)^2)
  part2   = (n0*sig20+sumres2)/2
  sig2    = 1/rgamma(1,part1,part2)
  
  # Sampling from full conditional of beta
  B1 = solve(solve(B0)+t(X)%*%X/sig2)
  b1 = B1%*%(solve(B0)%*%b0+t(X)%*%y/sig2)
  L  = t(chol(B1))
  b  = b1 + L%*%rnorm(2)
  
  # Storing the draws
  draws[iter,] = c(b[1],b[2],sig2)
}

pdf(file="mcmcdraws.pdf",width=10,height=15)
par(mfrow=c(3,2))
plot(draws[10:niter,1],xlab="Interations",ylab="",main="beta1",type="l")
hist(draws[10:niter,1],xlab="",main="",prob=TRUE)
plot(draws[10:niter,2],xlab="Interations",ylab="",main="beta2",type="l")
hist(draws[10:niter,2],xlab="",main="",prob=TRUE)
plot(draws[10:niter,3],xlab="Interations",ylab="",main="sigma2",type="l")
hist(draws[10:niter,3],xlab="",main="",prob=TRUE)
dev.off()

write(t(draws),"mcmcdraws.out",ncol=3)

