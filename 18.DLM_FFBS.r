# This is a Gibbs sampler implementation for the analysis of first-order DLMs.
# The model is
#      y_t   = theta_t + epsilon_t, epsilon_t ~ N(0,sigma2),
#    theta_t = theta_{t-1} + nu_t,  nu_t ~ N(0,sigma2 * W).
# In this parametrization, W is a signal-to-noise ratio parameter. 
#
# The priors are the following:
# theta_1 ~ N(a_1, R_1)
# sigma2 ~ IG(nsig/2, nssig/2)
# W ~ IG(nW/2, nsW/2)
# 
# The latent variables theta_{1:T} are simulated jointly using FFBS.
# sigma2 and W are simulated from their full condiotinal distributions 
# that are inverse gamma distributions.

# By M. A. R. Ferreira (April 2012, updated April 2015).

require(dlm)
require(MCMCpack)
theta.simulation <- function(Tbig,y,sigma2,W,a1,R1)
{ 
  # Define the auxiliary vectors:
  a <- rep(0,Tbig)
  R <- rep(0,Tbig)
  f <- rep(0,Tbig)
  Q <- rep(0,Tbig)
  A <- rep(0,Tbig)
  e <- rep(0,Tbig)
  m <- rep(0,Tbig)
  C <- rep(0,Tbig)
  # Kalman filter:
  a[1] <- a1
  R[1] <- R1
  f[1] <- a[1]
  Q[1] <- R[1]+sigma2
  A[1] <- R[1] / Q[1]
  e[1] <- y[1] - f[1]
  C[1] <- 1.0 / (1.0/R[1] + 1.0/sigma2)          # Equivalent to C[1] <- R[1] - A[1]^2 * Q[1]  and 
  m[1] <- C[1] * (a[1]/R[1] + y[1]/sigma2)       # m[1] <- a[1] + A[1] * e[1], but numerically more stable
  for(t in 2:Tbig)
  {
    a[t] <- m[t-1]
    R[t] <- C[t-1] + W * sigma2
    f[t] <- a[t]
    Q[t] <- R[t]+sigma2
    A[t] <- R[t] / Q[t]
    e[t] <- y[t] - f[t]
    m[t] <- a[t] + A[t] * e[t]
    C[t] <- R[t] - A[t]^2 * Q[t]
  }

  # Now, the backward sampler:
  theta[Tbig] <- rnorm(1,m[Tbig],sqrt(C[Tbig]))
  for(t in (Tbig-1):1)
  {
    H <- 1.0 / ( 1.0/(W * sigma2) + 1.0/C[t] )
    h <- H * ( theta[t+1]/(W * sigma2) + m[t]/C[t] )
    theta[t] <- rnorm(1,h,sqrt(H))
  }
  theta
}

sigma2.simulation <- function(Tbig,y,theta,W,n,ns)
{ # This is a Gibbs step.
  a <- n + 2*Tbig - 1
  b <- ns + sum((y-theta)^2) + sum((theta[2:Tbig]-theta[1:(Tbig-1)])^2)/W 
  sigma2 <- 1.0 / rgamma(1,a/2,b/2)
  sigma2
}

W.simulation <- function(Tbig,y,theta,sigma2,n,ns)
{ # This is a Gibbs step.
  a <- n + Tbig - 1
  b <- ns + sum((theta[2:Tbig]-theta[1:(Tbig-1)])^2) / sigma2
  W <- 1.0 / rgamma(1,a/2,b/2)
  W
}
require(MCMCpack)
drawIGpost <- function(x, a=0.001, b=0.001) {
  return(rinvgamma(1, a+length(x)/2, b+sum(x^2)/2))
}


# Initializing the parameters, latent and observation vectors.
sigma2 <- 40
sigmat2 <- 20
sigmad2 <- 0.2
W <- 5
Tbig <- 200
theta_org <- rep(0,Tbig)
delta <- rep(0,Tbig)
y <- rep(0,Tbig)

# Simulating the data
theta_org[1] <- rnorm(1,0,sqrt(sigmat2))
delta[1] <- rnorm(1,0,sqrt(sigmad2))
y[1] <- theta[1] + rnorm(1,0,sqrt(sigma2))
for(t in 2:Tbig)
{
  theta_org[t] <- theta_org[t-1] + delta[t-1]+rnorm(1,0,sqrt(sigmat2))
  delta[t] <- delta[t-1] + rnorm(1,0,sqrt(sigmad2))
  y[t] <- theta_org[t] + rnorm(1,0,sqrt(sigma2))
}

plot(y,type="l",lwd=2)
lines(theta_org,lwd=2,col="red")

# Gibbs sampler
G <- 5000

sigma2G <- rep(0,G)
WG <- rep(0,G)
thetaG <- matrix(0,nrow=2,ncol=Tbig)
thetaG1 <- matrix(0,nrow=G,ncol=Tbig)
thetaG2 <- matrix(0,nrow=G,ncol=Tbig)

sigma2G[1] <- 1
WG[1] <- 1
thetaG[1,] <- rep(0,Tbig)
thetaG2[1,]<- rep(0,Tbig)

a1 <- 0
R1 <- 10^20
nsig <- 0.001
nssig <- 0.001
nW <- 0.001
nsW <- 0.001
#nW  <- 6.2   # A priori W has the first 3 moments.
#nsW <- 0.42  # A a priori E[W] = 0.1
# nW  <- 0
# nsW <- 0

# for(g in 2:G)
# {
#   thetaG <- theta.simulation(Tbig,y,sigma2G[g-1],WG[g-1],a1,R1)
#   thetaG1[g,] <- thetaG[1,]
#   thetaG2[g,]<- thetaG[2,]
#   
#   sigma2G[g] <- sigma2.simulation(Tbig,y,thetaG1[g,],WG[g-1],nsig,nssig)
#   
#   WG[g] <-    W.simulation(Tbig,y,thetaG1[g,],sigma2G[g],nW,nsW)
#   WG2[g] <-    W.simulation(Tbig,y,thetaG2[g,],sigma2G[g],nW,nsW)
# }
V=5;
W1=20;
W2=2;

theta1.reps <- matrix(0,nrow=G,ncol=Tbig+1)
theta2.reps <- matrix(0,nrow=G,ncol=Tbig+1)
V.reps <- rep(0,G)
W1.reps <- rep(0,G)
W2.reps <- rep(0,G)

for (i in 1:G) {
  cat(i,"\n")
  # Sample states
  mod <- dlmModPoly(2, dV=V, dW=c(W1,W2))
  filt <- dlmFilter(y, mod)
  theta <- dlmBSample(filt)
  # Sample V and W
  V <- drawIGpost(y-theta[1:length(y),1])
  W1 <- drawIGpost((theta[2:Tbig,1]-theta[1:(Tbig-1),1]-theta[1:(Tbig-1),2]))
  W2 <- drawIGpost((theta[2:Tbig,2]-theta[1:(Tbig-1),2]))
   # V <-  sigma2.simulation(Tbig,y,theta[,1],W1,nsig,nssig)
   # W1 <- W.simulation(Tbig,y,theta[,1],V,nW,nsW)
   # W2 <- W.simulation(Tbig,y,theta[,2],V,nW,nsW)
  
  # Save iterations
  V.reps[i] <- V
  W1.reps[i] <- W1
  W2.reps[i] <- W2
  theta1.reps[i,] = theta[,1]
  theta2.reps[i,] = theta[,2]
}


par(mfrow=c(2,2))
plot(V.reps,type="l")
plot(W1.reps,type="l")
hist(V.reps[1000:G])
hist(W1.reps[1000:G])



plot(theta1.reps[,1],type="l")
plot(theta1.reps[,10],type="l")
plot(theta1.reps[,20],type="l")
plot(theta1.reps[,30],type="l")
plot(theta1.reps[,40],type="l")
plot(theta1.reps[,50],type="l")

acf(theta1.reps[,1])
acf(theta1.reps[,10])
acf(theta1.reps[,20])
acf(theta1.reps[,30])
acf(theta1.reps[,40])
acf(theta1.reps[,50])
acf(theta1.reps[,60])
acf(theta1.reps[,70])
acf(theta1.reps[,80])
acf(theta1.reps[,90])
acf(theta1.reps[,100])

thetamean <- apply(theta1.reps[1000:G,],2,mean)
thetavar <- apply(theta1.reps[1000:G,],2,var)

ll = min(c(min(thetamean-2*sqrt(thetavar)),y))
ul = max(c(max(thetamean+2*sqrt(thetavar)),y))
zero = rep(0,length(y))

par(mfrow=c(1,1))
plot(thetamean,lty=1,ylim=c(ll,ul),type="l",lwd=2)
lines(theta_org,col="red",lwd=2)
lines(thetamean+1.96*sqrt(thetavar),lty=2,lwd=2)
lines(thetamean-1.96*sqrt(thetavar),lty=2,lwd=2)
lines(zero,lty=3)
points(y,lwd=2)

