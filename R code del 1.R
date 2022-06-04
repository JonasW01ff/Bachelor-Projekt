#install.packages("numDeriv", dependencies=TRUE)
#install.packages("matlib", dependencies=TRUE)

#library("numDeriv")
library("MASS")
#library("scatterplot3d")
library("quantmod")
library("ggplot2")
library("reshape2")
library("matlib")

logLikeFunc <- function(X,kappa,theta,sigma2){
  "
  Tager daglig data og retunere loglikelihood funktionsv?rdien,
  dvs. $Delta$ t er antaget at v?re 1, og vi m?ler i dage.
  
  ----
  param X
  ::Type list af numerics eller ?n numeric
    En liste eller ?t tal med observation(erne/en)
  
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
  ::Type list af numerics eller ?n numeric
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

KvadratiskTeststoerelse <- function(X){
  "
  Givet en loglikelihood funktion og observationerne, findes
  den kvadratiske score teststoerelse konfidens regioner i form
  af en punkt sky
  
  ----
  param X
  ::Type list af numerics eller ?n numeric
    En liste med observationerne
  
  param logLikeFunc
  ::Type R funktion
    En R funktion der svare til din loglikelihood funktion
  
  return
  ::Type 6 vectors
    kappa'er <= taerskel, kappa'er over taerskel, ligeså med theta'er og derefter sigma^2'er
  "
  Kappas = c()
  kappas = c()
  Thetas = c()
  thetas = c()
  Sigmas2 = c()
  sigmas2 = c()
  taerskel = qchisq(.05, df=3)
  for (kappa in seq(-1,1,0.03)){
    print(kappa)
    for (theta in seq(-1,1,0.1)){
      for (sigma2 in seq(0.01,1,0.1)){
        #print(c(kappa,theta,sigma2))
        if(is.nan(logLikeFunc(X,kappa,theta,sigma2))){next}
        logLikeFunc_data <- function(arg){
          return(logLikeFunc(X,arg[1],arg[2],arg[3]))
        }
        S <- c(NaN,NaN,NaN)
        S <- numDeriv::grad(logLikeFunc_data,c(kappa,theta,sigma2))
        #for (i in 1:3){
        #  vars <- c(kappa,theta,sigma2)
        #  func <- function(vari) {
        #    arg <- c(kappa,theta,sigma2)
        #    arg[i] <- vari
        #    return(logLikeFunc(X,arg[1],arg[2],arg[3]))}
        #  S[i] <- mean(numDeriv::grad(func,vars[i], method="Richardson")) 
        #}
        I <- matrix(NaN,nrow=3,ncol=3)
        I <- -numDeriv::hessian(logLikeFunc_data,c(kappa,theta,sigma2))
        #for (i in 1:3){
        #  for (j in 1:3){
        #    vars <- c(kappa,theta,sigma2)
        #    func_base <- function(vari_base){
        #      func <- function(vari) {
        #        arg <- c(kappa,theta,sigma2)
        #        arg[i] <- vari
        #        arg[j] <- vari_base
        #        return(logLikeFunc(X,arg[1],arg[2],arg[3]))}
        #      return(numDeriv::grad(func,vars[i], method="Richardson"))
        #    }
            #I[i,j] <- -mean(numDeriv::grad(func_base,vars[j], method="Richardson"))
        #    logLikeFunc_data <- function(kappa,theta,sigma2){
        #      return(logLikeFunc(X,kappa,theta,sigma2))
        #    }
        #    I <- -mean(numDeriv::hessian(logLikeFunc_data,kappa=kappa,theta=theta,sigma2=sigma2))
        #  }
        #}

        Ii <- MASS::ginv(I)
        teststoerelse <- length(X)^(-2)*t(S)%*%Ii%*%S
        if (teststoerelse <= taerskel){
          # add params til Cap
          Kappas <- c(Kappas,kappa)
          Thetas <- c(Thetas,theta)
          Sigmas2 <- c(Sigmas2, sigma2)

        }
        else{
          # add params til low
          kappas <- c(kappas,kappa)
          thetas <- c(thetas,theta)
          sigmas2 <- c(sigmas2, sigma2)
        }
        
      }
    }
  }
  return(cbind(Kappas,Thetas,Sigmas2))
}

asympMLE95 <- function(data,hat_kappa,hat_theta,hat_sigma2,dt=1,HESS=FALSE){
  kappa_ <- hat_kappa
  theta_ <- hat_theta
  sigma2_ <- hat_sigma2
  
  logLikeFunc_data <- function(arg){
    return(logLikeFunc(data,arg[1],arg[2],arg[3]))
  }
  
  
  ell_grad <- function(data,theta,kappa,sigma2){
    x_1 <- head(data,-1)
    x_0 <- tail(data,-1)
    kappa_ <- kappa
    theta_ <- theta
    sigma2_ <- sigma2
    s2_      <- function(){return( sigma2_*(1-exp(-2*kappa_*dt))/(2*kappa_))}
    mu_      <- function(x_1){return( x_1*exp(-kappa_*dt)+theta_*(1-exp(-kappa_*dt)))}
    mu_theta <- function(){return( 1-exp(-kappa_*dt))}
    mu_kappa <- function(x_1){return( -x_1/dt*exp(-kappa_*dt)+theta_/dt*exp(-kappa_*dt))}
    s2_kappa <- function(){return( (sigma2_/dt*exp(-2*kappa_*dt)*kappa_-2*sigma2_*(1-exp(-2*kappa_*dt)))/(4*kappa_^2))}
    s2_sigma <- function(){return( sqrt(sigma2_)/kappa_*(1-exp(-2*kappa_*dt)))}
    phi      <- function(y,mu,s2){return( exp(-(y-mu)^2/(2*s2))/sqrt(s2*2*pi))}
    phi_mu   <- function(y,mu,s2){return( phi(y,mu,s2)*(y-mu)/s2)}
    phi_s2   <- function(y,mu,s2){return( phi(y,mu,s2)*((y-mu)^2/(2*s2^2)-1/(2*s2)))}
    out = c(NaN,NaN,NaN)
    out[1] <- sum((phi_mu(x_0,mu_(x_1),s2_())*mu_kappa(x_1) + phi_s2(x_0,mu_(x_1),s2_())*s2_kappa()   )/phi(x_0,mu_(x_1),s2_())) # ell_kappa
    out[2] <- sum(phi_mu(x_0,mu_(x_1),s2_())*mu_theta()/phi(x_0,mu_(x_1),s2_())) # ell_theta
    out[3] <- sum(phi_s2(x_0,mu_(x_1),s2_())*s2_sigma()/phi(x_0,mu_(x_1),s2_())) # ell_sigma
    return( out)
  }
  ell_kappa_data <- function(arg){return(ell_grad(data,arg[1],arg[2],arg[3])[1])}
  ell_theta_data <- function(arg){return(ell_grad(data,arg[1],arg[2],arg[3])[2])}
  ell_sigma_data <- function(arg){return(ell_grad(data,arg[1],arg[2],arg[3])[3])}
  

  if (HESS==FALSE){
    I_row1 <- -numDeriv::grad(ell_kappa_data,c(kappa_,theta_,sigma2_))
    I_row2 <- -numDeriv::grad(ell_theta_data,c(kappa_,theta_,sigma2_))
    I_row3 <- -numDeriv::grad(ell_sigma_data,c(kappa_,theta_,sigma2_))
    I <- matrix(c(I_row1,I_row2,I_row3),nrow=3,byrow=TRUE)
    Iinv <- inv(I)
  }
  else{
    I <- -numDeriv::hessian(logLikeFunc_data,c(kappa_,theta_,sigma2_))
    Iinv <- MASS::ginv(I)
  }
  kappa_ <- hat_kappa
  theta_ <- hat_theta
  sigma2_ <- hat_sigma2
  konfUP <- c(hat_kappa, hat_theta, hat_sigma2)+ sqrt(diag(Iinv))*qnorm(0.975)
  konfDOWN <- c(hat_kappa, hat_theta, hat_sigma2)- sqrt(diag(Iinv))*qnorm(0.975)
  output <- matrix(c(c("","Kappa","Theta","Sigma^2"),c("Upper asymp. 95%-konf.",konfUP),c("Lower asymp. 95%-konf.",konfDOWN)),nrow=3, byrow=TRUE)
  
  return(output)
}

MLEcheck <- function(kappa_,theta_, sigma2_,n,n_simul){
  b_ = theta_*(1-exp(-kappa_))
  a_ = exp(-kappa_)
  v_ <- sigma2_*(1-exp(-2*kappa_))/(2*kappa_)
  hat_kappas <- c()
  hat_thetas <- c()
  hat_sigmas2<- c()
  for (simulation in 1:n_simul){
    #print(c("Simulation:",simulation))
    Xt <- 0
    Xtu <- Xt # last X
    X_1 <- c(Xt)
    X <- c()
    for (u in 1:n){
      u <- 1
      middel <- exp(-kappa_*u)*Xtu+theta_*(1-exp(-kappa_*u))
      varians <- sigma2_*(1-exp(-2*kappa_*u))/(2*kappa_)
      Xtu <- rnorm(1,middel,varians^0.5)
      X <- c(X,Xtu)
      X_1 <- c(X_1,Xtu)
    }
    data <- X_1
    
    X_1 <- head(X_1,-1)
    x0 <- sum(X)
    x1 <- sum(X_1)
    x01 <- sum(X*X_1)
    x11 <- sum(X_1*X_1)
    
    hat_b <- sum(X-X_1*(x01/x11) )/sum(1-X_1*(x1/x11) )
    hat_a <- (x01-x1*hat_b)/x11
    hat_v <- mean((X-X_1*hat_a-hat_b)^2)
    
    hat_kappa <- -log(hat_a)
    hat_theta <- hat_b/(1-exp(-hat_kappa))
    hat_sigma2 <- 2*hat_kappa*hat_v/(1-exp(-2*hat_kappa))
    
    hat_kappas <- c(hat_kappas, hat_kappa)
    hat_thetas <- c(hat_thetas, hat_theta)
    hat_sigmas2<- c(hat_sigmas2,hat_sigma2)
  }
  print(asympMLE95(data, hat_kappa, hat_theta, hat_sigma2))
  hat_kappas <- na.omit(hat_kappas)
  hat_thetas <- na.omit(hat_thetas)
  hat_sigmas2 <- na.omit(hat_sigmas2)
  output <- matrix(c( c("","2.5%","97.5%"),
  c("simul. 95% Kappa in:",quantile(hat_kappas,0.025),quantile(hat_kappas,0.975)),
  c("simul.95% Theta in:",quantile(hat_thetas,0.025),quantile(hat_thetas,0.975)),
  c("simul.95% Sigma^2 in:",quantile(hat_sigmas2,0.025),quantile(hat_sigmas2,0.975))),nrow=4,byrow=TRUE)
  print(output)


  kappa_ <- hat_kappa
  theta_ <- hat_theta
  sigma2_ <- hat_sigma2
  
}
# Når der kun køres 5
MLEcheck(0.1,0.1,0.1,5,100)
MLEcheck(0.1,0.1,0.1,100,100)
MLEcheck(0.05,0.2,0.4,5,100)
MLEcheck(0.05,0.2,0.4,100,100)
MLEcheck(0.05,0.2,0.4,1000,100)
MLEcheck(0.05,0.2,0.4,10000,100)
#MLEcheck(0.05,0.2,0.4,100000,100)

MLE95 <- function(X){
  data <- X
  X_1 <- head(X,-1)
  X <- tail(X,-1)
  x0 <- sum(X)
  x1 <- sum(X_1)
  x01 <- sum(X*X_1)
  x11 <- sum(X_1*X_1)
  
  hat_b <- sum(X-X_1*(x01/x11) )/sum(1-X_1*(x1/x11) )
  hat_a <- (x01-x1*hat_b)/x11
  hat_v <- mean((X-X_1*hat_a-hat_b)^2)
  
  hat_kappa <- -log(hat_a)
  hat_theta <- hat_b/(1-exp(-hat_kappa))
  hat_sigma2 <- 2*hat_kappa*hat_v/(1-exp(-2*hat_kappa))
  print(c("MLE kappa:",hat_kappa))
  print(c("MLE theta:",hat_theta))
  print(c("MLE sigma2:",hat_sigma2))
  
  kappa_ <- hat_kappa
  theta_ <- hat_theta
  sigma2_ <- hat_sigma2
  b_ = theta_*(1-exp(-kappa_))
  a_ = exp(-kappa_)
  v_ <- sigma2_*(1-exp(-2*kappa_))/(2*kappa_)
  n <- length(X)
  hat_kappas <- c()
  hat_thetas <- c()
  hat_sigmas2<- c()
  for (simulation in 1:100){
    #print(c("Simulation:",simulation))
    Xt <- 0
    Xtu <- Xt # last X
    X_1 <- c(Xt)
    X <- c()
    for (u in 1:n){
      u <- 1
      middel <- exp(-kappa_*u)*Xtu+theta_*(1-exp(-kappa_*u))
      varians <- sigma2_*(1-exp(-2*kappa_*u))/(2*kappa_)
      Xtu <- rnorm(1,middel,varians^0.5)
      X <- c(X,Xtu)
      X_1 <- c(X_1,Xtu)
    }
    
    X_1 <- head(X_1,-1)
    x0 <- sum(X)
    x1 <- sum(X_1)
    x01 <- sum(X*X_1)
    x11 <- sum(X_1*X_1)
    
    hat_b <- sum(X-X_1*(x01/x11) )/sum(1-X_1*(x1/x11) )
    hat_a <- (x01-x1*hat_b)/x11
    hat_v <- mean((X-X_1*hat_a-hat_b)^2)
    
    hat_kappa <- -log(hat_a)
    hat_theta <- hat_b/(1-exp(-hat_kappa))
    hat_sigma2 <- 2*hat_kappa*hat_v/(1-exp(-2*hat_kappa))
    
    hat_kappas <- c(hat_kappas, hat_kappa)
    hat_thetas <- c(hat_thetas, hat_theta)
    hat_sigmas2<- c(hat_sigmas2,hat_sigma2)
  }
  hat_kappas <- na.omit(hat_kappas)
  hat_thetas <- na.omit(hat_thetas)
  hat_sigmas2 <- na.omit(hat_sigmas2)
  print(c("95% Kappa in:",quantile(hat_kappas,0.025),quantile(hat_kappas,0.975)))
  print(c("95% Theta in:",quantile(hat_thetas,0.025),quantile(hat_thetas,0.975)))
  print(c("95% Sigma^2 in:",quantile(hat_sigmas2,0.025),quantile(hat_sigmas2,0.975)))
  
  logLikeFunc_data <- function(arg){
    return(logLikeFunc(data,arg[1],arg[2],arg[3]))
  }
}

SDEcreator <- function(kappa,theta,sigma2,n,Xt){
  X <- c(Xt)
  for (i in 1:(n-1)){
    dWt <- rnorm(1,0,1)
    Xdt <- kappa*(theta-Xt)+sigma2*dWt
    Xt <- Xt+Xdt
    X <- c(X,Xt)
  }
  return(X)
}

# <- KvadratiskTeststoerelse(X)
#scatterplot3d(punktsky)
#print("kappa, confidens interval: ")
#print(c(min(punktsky[,1]),max(punktsky[,1])))
#print("theta, confidens interval: ")
#print(c(min(punktsky[,2]),max(punktsky[,2])))
#print("sigma2, confidens interval: ")
#print(c(min(punktsky[,3]),max(punktsky[,3])))

# VIX data
print("VIK-data")
getSymbols("^VIX",src="yahoo")
X<-as.vector(VIX[,6])#head(VIX[,6],100)
par(mfrow=c(1,2))
qqnorm(X, pch = 1, frame = FALSE, main="Normal QQ-plot af VIX")
qqline(X, col = "steelblue", lwd = 2)
qqnorm(log(X), pch = 1, frame = FALSE, main="Normal QQ-plot af ln(VIX)")
qqline(log(X), col = "steelblue", lwd = 2)
MLE95(X)
MLE95(log(X))

X_1 <- head(X,-1)
X <- tail(X,-1)
x0 <- sum(X)
x1 <- sum(X_1)
x01 <- sum(X*X_1)
x11 <- sum(X_1*X_1)

hat_b <- sum(X-X_1*(x01/x11) )/sum(1-X_1*(x1/x11) )
hat_a <- (x01-x1*hat_b)/x11
hat_v <- mean((X-X_1*hat_a-hat_b)^2)

hat_kappa <- -log(hat_a)
hat_theta <- hat_b/(1-exp(-hat_kappa))
hat_sigma2 <- 2*hat_kappa*hat_v/(1-exp(-2*hat_kappa))


linSDE <- SDEcreator(hat_kappa,hat_theta,hat_sigma2,length(VIX[,6]),as.numeric( VIX[1,6]))
#dim(linGBM) <- c(length(linGBM),1)
data <- data.frame(linSDE, VIX[,6])
data$day <- rownames(data)

ggplot(melt(data, id="day")) + geom_line(aes(x = day, y=value, group=variable, colour=variable)) + xlab("Dag")+ylab("VIX")+   theme(legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                    legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                    legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                      legend.title = element_text(size=20), #change legend title font size
                                                                                                                                    legend.text = element_text(size=20), #change legend text font size
                                                                                                                                    axis.text=element_text(size=20),
                                                                                                                                    axis.title=element_text(size=20))


X<-as.vector(VIX[,6])
X <- log(X)
X_1 <- head(X,-1)
X <- tail(X,-1)
x0 <- sum(X)
x1 <- sum(X_1)
x01 <- sum(X*X_1)
x11 <- sum(X_1*X_1)

hat_b <- sum(X-X_1*(x01/x11) )/sum(1-X_1*(x1/x11) )
hat_a <- (x01-x1*hat_b)/x11
hat_v <- mean((X-X_1*hat_a-hat_b)^2)

hat_kappa <- -log(hat_a)
hat_theta <- hat_b/(1-exp(-hat_kappa))
hat_sigma2 <- 2*hat_kappa*hat_v/(1-exp(-2*hat_kappa))
loglinSDE <- SDEcreator(hat_kappa,hat_theta,hat_sigma2,length(VIX[,6]),log(as.numeric( VIX[1,6])))
data <- data.frame(loglinSDE, log(VIX[,6]))
names(data) <- c("loglinSDE","logVIX_Adjusted")
data$day <- rownames(data)
ggplot(melt(data, id="day")) + geom_path(aes(x = day, y=value, group=variable, colour=variable))+ xlab("Dag")+ylab("VIX")+   theme(legend.key.size = unit(1, 'cm'), #change legend key size
                                      legend.key.height = unit(1, 'cm'), #change legend key height
                                      legend.key.width = unit(1, 'cm'), #change legend key width
                                      legend.title = element_text(size=20), #change legend title font size
                                      legend.text = element_text(size=20), #change legend text font size
                                      axis.text=element_text(size=20),
                                      axis.title=element_text(size=20))


######
## Empirical hedging
######
print("S&P 500-data")
getSymbols("^GSPC",src="yahoo")
X_SPX<-as.vector(GSPC[,6])
data <- cbind(VIX[,6]/as.numeric(VIX[1,6])*100,GSPC[,6]/as.numeric(GSPC[1,6])*100)
data <- data.frame(data)
data["day"] = index(VIX)
par(mfrow=c(1,1))
data <- melt(data ,  id.vars = 'day', variable.name = 'indeks')
ggplot(data, aes(day,value)) + geom_line(aes(colour = indeks))



