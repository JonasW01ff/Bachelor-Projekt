#install.packages("numDeriv")
#install.packages("comprehenr")
#install.packages("Hmisc")
#install.packages("matlib")
library("numDeriv")
library("matlib")
#library("comprehenr")
#library("Hmisc")

logLikeFunc <- function(X,kappa,theta,sigma2){
  "
  Tager daglig data og retunere loglikelihood funktionsværdien,
  dvs. $\Delta$ t er antaget at være 1, og vi måler i dage.
  
  ----
  param X
  ::Type list af numerics eller én numeric
    En liste eller ét tal med observation(erne/en)
  
  param kappa
  ::Type numeric
    kappa parameteren som et tal
  
  param theta
  ::Type numeric
    theta parameteren som et tal
  
  param sigma2
  ::Type numeric
    sigma^2 parameteren som et tal
  
  return
  ::Type list af numerics eller én numeric
    funktionsvaerdi(erne/en)
  
  "
  dt <- 1 # Bare en FYI
  a <- exp(-kappa*dt)
  b <- theta*(1-exp(-kappa*dt))
  v <- sigma2*(1-exp(-2*kappa*dt))
  logLike <- 0
  for (i in 2:length(X)){
    logLike <- logLike + -1/2*log(v*2*pi)-(X[i]-X[i-1]*a-b)^2/(2*v)
  }
  return(logLike)
}

KvadratiskTeststoerelse <- function(logLikeFunc,X,kappa,theta,sigma2){
  "
  Givet en loglikelihood funktion og observationerne, findes
  den kvadratiske score test størelse
  
  ----
  param X
  ::Type list af numerics eller én numeric
    En liste med observationerne
  
  param kappa
  ::Type numeric
    kappa parameteren som et tal
  
  param theta
  ::Type numeric
    theta parameteren som et tal
  
  param sigma2
  ::Type numeric
    sigma^2 parameteren som et tal
  
  param logLikeFunc
  ::Type R funktion
    En R funktion der svare til din loglikelihood funktion
  
  return
  ::Type numeric
    Teststørelsen $Q_n$ j.f. Lauritzen2021 afsnit 5.4.2.
  "
  S <- mean(numDeriv::grad(logLikeFunc,X,kappa,theta,sigma2))
  Ii <- matlib::inv(mean(numDeriv::hessian(logLikeFunc,X,kappa,theta,sigma2)))
  teststoerelse <- n*S*Ii*t(S)
  return(teststoerelse)
}