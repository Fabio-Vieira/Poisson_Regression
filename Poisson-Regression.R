#Bayesian Poisson Regression, with prior Student's t and Metropolis Hastings
#with random walk log adaptive proposal
###########################################################################
####Simulating regression model

N <- 20
beta <- c(0.5,2)
x <- cbind(rep(1,N),runif(N,min=0,max=1))
lambda <- exp(x%*%beta)
y <- rpois(N,lambda = lambda)
cbind(y,x)

hist(y)

###########################################################################
#Posterior distribution of Betas with prior Student's with df = 1, as in
#http://www.stat.columbia.edu/~gelman/research/published/priors11.pdf

PostBeta <- function(Beta, Y, X){
  N <- length(Y)
  p <- length(Beta)
  Sigma <- solve(diag(1,p,p))
  #log of the core of a multivariate student's t, with df = 1
  prior <- -(1+p)/2 * log(1 + t(beta)%*%Sigma%*%beta)
  #log of poisson likelihood
  like <- sum(-exp(X%*%Beta) + X%*%Beta * Y)
  return(prior + like) #bayes rule
}

#Metropolis Hastings algorithm with log adaptive proposal, as described in
#http://www.personal.psu.edu/bas59/papers/adaptive_techreport.pdf

updateBeta <- function(Y,X,Beta,Sm,Sigma,t,BetaMean){
  accept <- NULL 
  count <- 0
  c <- .8
  ggamma <- 1/t^c
  proposal <- Beta + sqrt(Sm*Sigma) * rnorm(1)
  prob <- min(1,exp(PostBeta(proposal,Y,X) - PostBeta(Beta,Y,X)))
  
  if(runif(1) < prob){
    accept <- proposal
    count <- 1
  } else {
    accept <- Beta
  }
  #Adapting the proposal
  lSm <- log(Sm)+ggamma*(prob-.234)
  Sm <- exp(lSm)	
  Sigma <- Sigma+ggamma*(((Beta-BetaMean)^2)-Sigma)
  BetaMean <- BetaMean+ggamma*(Beta-BetaMean)
  return(list(accept,count,Sm,Sigma,BetaMean))
}

#################################################################################
#MCMC 
Niter <- 20000
Y <- y
X <- x
Beta.out <- array(NA, dim = c(Niter,dim(X)[2]))
Count <- array(0, dim = Niter)

#Initial values
Beta.out[1,] <- beta
Sm <- 2.4^2 #Described in Shaby and Wells
Sigma <- 1
BetaMean <- 1

for(i in 2:Niter){
  Updating <- updateBeta(Y,X,Beta.out[i-1,],Sm,Sigma,i,BetaMean)
  Beta.out[i,] <- Updating[[1]]
  Count[i] <- Updating[[2]]
  Sm <- Updating[[3]]
  Sigma <- Updating[[4]]
  BetaMean <- Updating[[5]]
  print(i)
}

#Checking convergence
plot(Beta.out[,1],type='l')
abline(h=beta[1],col='red')
hist(Beta.out[,1])
abline(v=quantile(Beta.out[,1],probs=0.025),lty=2,col='blue')
abline(v=quantile(Beta.out[,1],probs=0.975),lty=2,col='blue')
abline(v=beta[1],col='red')
mean(Beta.out[,1])

plot(Beta.out[,2],type='l')
abline(h=beta[2],col='red')
hist(Beta.out[,2])
abline(v=quantile(Beta.out[,2],probs=0.025),lty=2,col='blue')
abline(v=quantile(Beta.out[,2],probs=0.975),lty=2,col='blue')
abline(v=beta[2],col='red')
mean(Beta.out[,2])
