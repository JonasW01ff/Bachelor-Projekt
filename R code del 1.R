#install.packages("numDeriv", dependencies=TRUE)
#install.packages("matlib", dependencies=TRUE)

library("numDeriv")
library("runner")
library("MASS")
#library("scatterplot3d")
library("quantmod")
library("ggplot2")
library("reshape2")
library("matlib")

##########
##  Empirical analysis of a linear SDE
##########
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

  #dt <- 1 # Bare en FYI
  a <- exp(-kappa*dt)
  b <- theta*(1-exp(-kappa*dt))
  v <- sigma2*(1-exp(-2*kappa*dt))
  #logLike <- sum(-1/2*log(v*2*pi)-(X[2:length(X)]-X[1:length(X)-1]*a-b)^2/(2*v))
  #for (i in 2:length(X)){
  #  logLike <- logLike + -1/2*log(v*2*pi)-(X[i]-X[i-1]*a-b)^2/(2*v)
  #}
  return(sum(-1/2*log(v*2*pi)-(X[2:length(X)]-X[1:length(X)-1]*a-b)^2/(2*v)))
}

KvadratiskTeststoerelse <- function(X){
  "
  [Bruges ikke længere!]
  Givet en loglikelihood funktion og observationerne, findes
  den kvadratiske score teststoerelse konfidens regioner i form
  af en punkt sky
  
  ----
  param X
  ::Type list af numerics eller ?n numeric
    En liste med observationerne
  
  return
  ::Type 3 vectors
    Punktsky med variable der gav en kvadratisk testscore mindre end eller lig med tærsklen, i.e. chi^2 kvatilen.
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
        I <- matrix(NaN,nrow=3,ncol=3)
        I <- -numDeriv::hessian(logLikeFunc_data,c(kappa,theta,sigma2))


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
  "
  Bruger de manuelt udledte første afledte, og derefter ved numerisk approximation
  findes de anden afledte som gradient af de første afledte.
  Dette bliver til fisher informationen of at assymptotisk 95% konf. interval
  kan findes
  
  
  ------
  :param data:
    Type vector of floats
      vector med data der følger GBM
      
  :param hat_kappa:
    Type float
      MLE for kappa
  
  :param hat_theta:
    Type float
      MLE for theta
  
  :param hat_sigma:
    Type float
      MLE for sigma^2
      
  :param dt:
    Type float
      den ækvidistante distance mellem obervationener, i.e. data
  
  :param HESS:
    Type bool
      Hvis TRUE, så bruges numderiv hessian functionen i stedet for grad på de førsteafledte for at finde fisher informationen.
      Hvis FALSE, så bruges numderiv grad på de først afledte til at finde fischer informationen
  
  :return:
    Type Matrix
      3x4 matrix med de assymptotiske konfidens intervaller for de 3 MLE'er
  
  "
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
    I <- (I+t(I))/2 # giver noget symmetri bedre approx. for anden afledte
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
  "
  Finder MLE konfidensinterval ved simulation, og ved asymptotiske resultater
  og tjekker dem op mod hindanden ud fra givne MLE'er
  
  
  ------
  :param data:
    Type vector of floats
      vector med data der følger GBM
      
  :param hat_kappa:
    Type float
      MLE for kappa
  
  :param hat_theta:
    Type float
      MLE for theta
  
  :param hat_sigma:
    Type float
      MLE for sigma^2
      
  
  :return:
    Type print
      3x4 matrix med de assymptotiske konfidens intervaller for de 3 MLE'er
      4x3 matrix med de simulleret konfidensintervaller for de 3 MLE'er
  
  "
  b_ = theta_*(1-exp(-kappa_))
  a_ = exp(-kappa_)
  v_ <- sigma2_*(1-exp(-2*kappa_))/(2*kappa_)
  hat_kappas <- rep(NaN,n_simul)
  hat_thetas <- rep(NaN,n_simul)
  hat_sigmas2<- rep(NaN,n_simul)
  u <- 1
  varians <- sigma2_*(1-exp(-2*kappa_*u))/(2*kappa_)
  sd      <- varians^0.5
  konstant_led <- theta_*(1-exp(-kappa_*u))
  konstant_fakor <- exp(-kappa_*u)
  for (simulation in 1:n_simul){
    #print(c("Simulation:",simulation))
    Xt <- 0
    #Xtu <- Xt # last X
    X_1 <- rep(NaN,n+1)
    X_1[1] <- Xt
    X_1_pre <- rnorm(n,konstant_led,sd)
    for (u_dummy in 1:n){
      #middel <- konstant_fakor*Xtu+konstant_led#theta_*(1-exp(-kappa_*u))
      X_1[u_dummy+1] <- X_1_pre[u_dummy]+konstant_fakor*X_1[u_dummy]
    }
    data <- X_1
    X <- tail(X_1,-1)
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
    
    hat_kappas[simulation] <- hat_kappa
    hat_thetas[simulation] <- hat_theta
    hat_sigmas2[simulation]<- hat_sigma2
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
MLEcheck(0.05,0.2,0.4,100000,100)
MLEcheck(0.05,0.2,0.4,1000000,100)

MLE95 <- function(X, n_simul=100){
  "
  Finder MLE og konfidensinterval ved simulation,
  og printer dem
  
  
  ------
  :param data:
    Type vector of floats
      vector med data der følger GBM
      
  :return:
    Type print
      print 3 gange med MLE'er
      4x3 matrix med de simulleret konfidensintervaller for de 3 MLE'er
  
  "
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
  hat_kappas <- rep(NaN,n_simul)
  hat_thetas <- rep(NaN,n_simul)
  hat_sigmas2<- rep(NaN,n_simul)
  u <- 1
  varians <- sigma2_*(1-exp(-2*kappa_*u))/(2*kappa_)
  sd      <- varians^0.5
  konstant_led <- theta_*(1-exp(-kappa_*u))
  konstant_fakor <- exp(-kappa_*u)
  for (simulation in 1:n_simul){
    #print(c("Simulation:",simulation))
    Xt <- 0
    #Xtu <- Xt # last X
    X_1 <- rep(NaN,n+1)
    X_1[1] <- Xt
    X_1_pre <- rnorm(n,konstant_led,sd)
    for (u_dummy in 1:n){
      #middel <- konstant_fakor*Xtu+konstant_led#theta_*(1-exp(-kappa_*u))
      X_1[u_dummy+1] <- X_1_pre[u_dummy]+konstant_fakor*X_1[u_dummy]
    }
    X <- tail(X_1,-1)
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
    
    hat_kappas[simulation] <- hat_kappa
    hat_thetas[simulation] <- hat_theta
    hat_sigmas2[simulation]<- hat_sigma2
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
  "
  Simulere en sti ud fra givet kappa theta og sigma^2
  
  
  ------
  :param hat_kappa:
    Type float
      MLE for kappa
  
  :param hat_theta:
    Type float
      MLE for theta
  
  :param hat_sigma:
    Type float
      MLE for sigma^2
  
  :param n:
    Type float
      antal observationer der skal simuleres
  
  :param Xt:
    Type float
      Start værdi for den GBM.
      
  :return:
    Type vector
      En vector med den GBM sti.
  
  "
  X <- rep(NaN,n)
  X[1] <- Xt
  rnorms <- rnorm(n-1,0,1)
  for (i in 2:n){
    #dWt <- rnorm(1,0,1)
    #Xdt <- kappa*(theta-Xt)+sigma2*dWt
    #Xt <- Xt+Xdt
    X[i] <- X[i-1]+kappa*(theta-X[i-1])+sigma2*rnorms[i-1]
  }
  return(X)
}

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
##  Pricing and hedging in the Black-Scholes model
######

# Dette er gjort i rolfs fil.

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

par(mfrow=c(1,1))
qqnorm(X_SPX, pch = 1, frame = FALSE, main="Normal QQ-plot af S&P500")
qqline(X_SPX, col = "steelblue", lwd = 2)
#qqnorm(log(X_SPX), pch = 1, frame = FALSE, main="Normal QQ-plot af ln(S&P500)")
#qqline(log(X_SPX), col = "steelblue", lwd = 2)

MLE <- function(X){
  "
  Tager data, og giver MLE for OU-process tilbage
  
  
  -----
  :param X:
    Type vector
      Vector med floats der er din data
  
  :return:
    Type vector
      Returnere vector med dine tre MLE'er i rækkefølgen
      kappa,theta og sigma^2
  "
  
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
  return(c(hat_kappa,hat_theta,hat_sigma2))
}
MLE_GBM <- function(X){
  "
  Tager data, og giver MLE for OU-process tilbage
  
  
  -----
  :param X:
    Type vector
      Vector med floats der er din data
  
  :return:
    Type vector
      Returnere vector med dine tre MLE'er i rækkefølgen
      beta, sigma^2
  "
  
  data <- X
  X_1 <- head(X,-1)
  X <- tail(X,-1)
  x0 <- sum(X)
  x1 <- sum(X_1)
  x01 <- sum(X*X_1)
  x11 <- sum(X_1*X_1)
  
  hat_b <- 0#sum(X-X_1*(x01/x11) )/sum(1-X_1*(x1/x11) )
  hat_a <- (x01-x1*hat_b)/x11
  hat_v <- mean((X-X_1*hat_a-hat_b)^2)
  
  hat_kappa <- -log(hat_a)
  hat_theta <- hat_b/(1-exp(-hat_kappa))
  hat_sigma2 <- 2*hat_kappa*hat_v/(1-exp(-2*hat_kappa))
  return(c(-hat_kappa,hat_sigma2))
}

MLE_GBM(X_SPX)

window <- 21
lambda <- 0.94
X_VIX <- as.vector(VIX[,6]/100)#head(VIX[,6],100)

for (i in 1:2){
  if (i==1){
    var <- X_SPX
    valname <- "SP500"
    }
  if (i==2){
    var <- X_VIX
    valname <- "VIX"
    }
  
  MLE_sigma_rolling_month <- runner(
    var,
    k=window,
    f= function(X){sqrt(MLE_GBM(X)[2])}
  )

  returns <-  c(NaN,tail(var,-1)/head(var,-1))
  EWMA <- runner(
    MLE_sigma_rolling_month,
    k=300,
    lag=-window,
    f=function(X){
      X <- rev(X)
      out <- 0
      for (i in 1:(window-1+1)){
        out <- lambda*out+(1-lambda)*X[i]^2
        #if (!is.na(out_)){out <- out_}
      }
      return(sqrt(out))
    }
  )
  
  #plot(MLE_sigma_rolling_month)
  print(c("rolling: ",mean(MLE_sigma_rolling_month,na.rm=TRUE)))
  print(c("EMVA",mean(EWMA, na.rm=TRUE)))
  
  #plot(EWMA)
  
  data <- data.frame(cbind(MLE_sigma_rolling_month,EWMA))
  if(i==1){
    data["dag"] = index(X_SPX)
    data <- melt(data ,  id.vars = 'dag', variable.name = 'variabel', value.name = "SP500")
    print(ggplot(data, aes(dag,SP500)) + geom_line(aes(colour = variabel)))
    }
  if(i==2){
    data["dag"] = index(X_VIX)
    data <- melt(data ,  id.vars = 'dag', variable.name = 'variabel', value.name = "VIX")
    print(ggplot(data, aes(dag,VIX)) + geom_line(aes(colour = variabel)))
    }  
}

# Split data into buisness months i.e. 21 days non overlapping

BlackScholesFormula  <- function (spot,timetomat,strike,r, q=0, sigma, opttype=1, greektype=1)
{ 
  "
  Docstring skrevet af Jonas Wolff
  
  -----
  param spot
  ::Type float
    Start kurs
  
  param timetomat
  ::Type float
    løbetiden
    
  param strike
  ::Type float
    Strike
    
  param r
  ::Type float
    risikofri rente
  
  param q
  ::Type float
    tids diskontering på param. spot
    
  param sigma
  ::Type float
    sigma i blackschols (i.e. wienerprocessens koeficient)
    
  param opttype
  ::Type float
    =1 er det call
  
  param greektype
  ::Type float
    =1 er det prisen
    =2 er det tid t delta'et
    
  return
  ::Type float
    Blacksholes formel output
    
  
  "
  #opttype <- 3
  d1<-as.numeric((log(spot/strike)+ ((r-q)+0.5*sigma^2)*timetomat)/(sigma*sqrt(timetomat)))
  d2<-as.numeric(d1-sigma*sqrt(timetomat))
  # spot = Y_o
  # timetomat = T-t
  # q = (beta-alpha)
  if (opttype==1 && greektype==1) result<-spot*exp(-q*timetomat)*pnorm(d1)-strike*exp(-r*timetomat)*pnorm(d2)
  
  if (opttype==2 && greektype==1) result<-spot*exp(-q*timetomat)*pnorm(d1)-strike*exp(-r*timetomat)*pnorm(d2)-spot*exp(-q*timetomat)+strike*exp(-r*timetomat)
  
  if (opttype==1 && greektype==2) result<-exp(-q*timetomat)*pnorm(d1)
  
  if (opttype==2 && greektype==2) result<-exp(-q*timetomat)*(pnorm(d1)-1)
  
  if (greektype==3) result<-exp(-q*timetomat)*dnorm(d1)/(spot*sigma*sqrt(timetomat))
  
  if (greektype==4) result<-exp(-q*timetomat)*spot*dnorm(d1)*sqrt(timetomat)
  
  
  if (opttype==3 && greektype==1) result<- 1*exp(-r*(timetomat))*pnorm(log(strike/spot),(r-1/2*sigma^2)*timetomat,sigma*sqrt(timetomat))
  if (opttype==3 && greektype==2) result<- -1*exp(-r*(timetomat))*dnorm(log(strike/spot),(r-1/2*sigma^2)*timetomat,sigma*sqrt(timetomat))/spot
  return(result)
  BlackScholesFormula<-result
  
}

BlackScholesImpVol  <- function (obsprice,spot,timetomat,strike,r, q=0, opttype=1)
{ difference<- function(sigBS, obsprice,spot,timetomat,strike,r,q,opttype)
{BlackScholesFormula (spot,timetomat,strike,r,q,sigBS, opttype,1)-obsprice
}
uniroot(difference, c(10^-6,10),obsprice=obsprice,spot=spot,timetomat=timetomat,strike=strike,r=r,q=q,opttype=opttype)$root
}


par(mfrow=c(1,1))
hedge_func <- function(data){ #Tyv stjålet fra rolfs kode
  #S0<-100
  S0 <- data[,1]-data[,1]+1000#as.numeric(St[1])
  r <- 0#0.02
  mu <- 0#MLE_GBM(data)[1]#0.02
  
  sigma <- 20.01
  sigma_hedge <- 20.01/sqrt(21)#MLE_GBM(X)[2]
  
  capT<-1
  strike<-1000
  
  Nhedge<-length(data[1,])-1#10#252

  Nrep<-length(data[,1])#1000  
  
  St<-data[,1]-data[,1]+1000#rep(S0, length=Nrep)
  dt<-capT/Nhedge
  initialoutlay<-BlackScholesFormula(S0,capT,strike, r,0,sigma_hedge,1,1)
  #initialoutlay<-BlackScholesFormula(S0,capT,strike, r,0,sigma,1,1)
  
  Vpf<- initialoutlay #rep(initialoutlay,length=Nrep)
  print(Vpf)
  
  a<-BlackScholesFormula(St,capT,strike, r,0,sigma_hedge,1,2)
  b<-Vpf-a*St
  
  for(i in 2:Nhedge){
    #St<-St*exp((mu-0.5*sigma^2)*dt +sigma*sqrt(dt)*rnorm(Nrep))	
    St <- data[,i]-data[,1]+1000
    Vpf<-a*St+b*exp(dt*r)    
    a<-BlackScholesFormula(St,(capT-(i-1)*dt),strike, r,0,sigma_hedge,1,2)
    b<-Vpf-a*St
  }
  
  #St<-St*exp((mu-0.5*sigma^2)*dt +sigma*sqrt(dt)*rnorm(Nrep))
  St <- data[,Nhedge+1]-data[,1]+1000
  Vpf<-a*St+b*exp(dt*r)
  optionpayoff<-pmax(St-strike,0)
  hedgeerror<-Vpf-optionpayoff
  #print(hedgeerror)
  
  
  # SUMMARY STATS & GRAPHS
  # ======================
  
  #print(paste("Initial investment =",round(initialoutlay,4)))
  print(paste("Average discounted option payoff =",round(exp(-r*capT)*mean(optionpayoff),4)))
  print(paste("Average discounted portfolio value =",round(exp(-r*capT)*mean(Vpf),4)))
  print(paste("Terminal Hedge error mean =",mean(abs(hedgeerror))))
  print(paste("Terminal Hedge error SD =",sd(hedgeerror)))
  
  ToFile<-FALSE
  
  if (ToFile) png(file="Scatter.png",res=300,width = 3*480, height = 3*480)
  plot(St,Vpf,col="blue",xlab="S(T)",ylab="Value of hedge portfolio",main="Discrete hedging of a call-option",xlim=c(50,1200),ylim=c(-5,105))
  text(50,100,paste("# hegde points =",Nhedge),adj=0)
  #text(50,95,paste("r-mu =",r-mu),adj=0)
  #text(50,90,paste("sigma-sigma_hedge =",sigma-sigma_hedge),adj=0)
  points(50:1200,pmax(50:1200 - strike,0),type='l',lwd=3)
  
  if (ToFile) dev.off()
  
}

bmonth <- data.frame(split(tail(X_SPX,length(X_SPX)%%21*(-1)),1:21)) #bmonth[1,] would give you the first buisness month

hedge_func(bmonth)




