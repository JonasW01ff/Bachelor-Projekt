
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
  d1<-(log(spot/strike)+ ((r-q)+0.5*sigma^2)*timetomat)/(sigma*sqrt(timetomat))
  d2<-d1-sigma*sqrt(timetomat)
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
  
  BlackScholesFormula<-result
  
}

BlackScholesImpVol  <- function (obsprice,spot,timetomat,strike,r, q=0, opttype=1)
{ difference<- function(sigBS, obsprice,spot,timetomat,strike,r,q,opttype)
{BlackScholesFormula (spot,timetomat,strike,r,q,sigBS, opttype,1)-obsprice
}
uniroot(difference, c(10^-6,10),obsprice=obsprice,spot=spot,timetomat=timetomat,strike=strike,r=r,q=q,opttype=opttype)$root
}

# INITIALIZE
# ==========

S0<-100
r<-0.02
mu<-0.02

sigma<-0.2
sigma_hedge<-0.2

capT<-1
strike<-100

Nhedge<-10#252
Nrep<-1000

# HEDGE
# =====

St<-rep(S0, length=Nrep)
dt<-capT/Nhedge
#initialoutlay<-BlackScholesFormula(S0,capT,strike, r,0,sigma_hedge,1,1)
initialoutlay<-BlackScholesFormula(S0,capT,strike, r,0,sigma,1,1)

Vpf<-rep(initialoutlay,length=Nrep)

a<-BlackScholesFormula(St,capT,strike, r,0,sigma_hedge,1,2)
b<-Vpf-a*St

for(i in 2:Nhedge){
  St<-St*exp((mu-0.5*sigma^2)*dt +sigma*sqrt(dt)*rnorm(Nrep))	
  Vpf<-a*St+b*exp(dt*r)    
  a<-BlackScholesFormula(St,(capT-(i-1)*dt),strike, r,0,sigma_hedge,1,2)
  b<-Vpf-a*St
}

St<-St*exp((mu-0.5*sigma^2)*dt +sigma*sqrt(dt)*rnorm(Nrep))
Vpf<-a*St+b*exp(dt*r)
optionpayoff<-pmax(St-strike,0)
hedgeerror<-Vpf-optionpayoff


# SUMMARY STATS & GRAPHS
# ======================

print(paste("Initial investment =",round(initialoutlay,4)))
print(paste("Average discounted option payoff =",round(exp(-r*capT)*mean(optionpayoff),4)))
print(paste("Average discounted portfolio value =",round(exp(-r*capT)*mean(Vpf),4)))
print(paste("Terminal Hedge error mean =",mean(abs(hedgeerror))))
print(paste("Terminal Hedge error SD =",sd(hedgeerror)))

ToFile<-FALSE

if (ToFile) png(file="Scatter.png",res=300,width = 3*480, height = 3*480)

plot(St,Vpf,col="blue",xlab="S(T)",ylab="Value of hedge portfolio",main="Discrete hedging of a call-option",xlim=c(50,200),ylim=c(-5,105))
text(50,100,paste("# hegde points =",Nhedge),adj=0)
text(50,95,paste("r-mu =",r-mu),adj=0)
text(50,90,paste("sigma-sigma_hedge =",sigma-sigma_hedge),adj=0)
points(50:200,pmax(50:200 - strike,0),type='l',lwd=3)

if (ToFile) dev.off()
