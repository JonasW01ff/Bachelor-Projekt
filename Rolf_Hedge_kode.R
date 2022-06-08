
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
  return(result)
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

Nhedge<-252
Nrep<-1000

reset_var <- function(){
  S0<<-100
  r<<-0.02
  mu<<-0.02
  
  sigma<<-0.2
  sigma_hedge<<-0.2
  
  capT<<-1
  strike<<-100
  
  Nhedge<<-252
  Nrep<<-1000
}
reset_var()
# HEDGE
# =====


hedgeerror_plotter <- function(hedgeerror,St,Vpf,strike,opttype_=1){
  St <- St-strike
  strike <- 0
  lim <- max(abs(hedgeerror))
  histdata <- hist(hedgeerror,plot=FALSE,freq=FALSE,breaks=pretty(c(-lim,lim),n=20))
  histdata <- data.frame(Data=histdata$density,Index=histdata$mids)
  rownames(histdata) <- histdata$Index
  if (opttype_ != 3){xlm <- 1}
  if (opttype_ == 3){xlm <- 5}
  xlm <- max(histdata$Data)*3
  bp <- barplot(height=histdata$Data,width=1,pos=0,xlim=c(xlm,0), axes=FALSE, frame.plot=FALSE,xlab="",ylab="",main="", horiz=TRUE)
  axis(side=4, at = bp, labels=histdata$Index)#,labels=histdata$mids)
  mtext("Hedge error", side = 4, line = 1.9, col = 1)
  axis(side=3, at = pretty(c(xlm*0.9,xlm*0.1),n=10), line=-2, labels = paste(pretty(c(xlm*0.9,xlm*0.1),n=10)*100,"%"))
  par(new=TRUE)
  
  xlim_top <- max(St)*0.8
  ylim_top <- max(abs(Vpf))*0.9
  xlim_bot <- min(St)
  
  plot(St,Vpf,col="blue",xlab="S(T)-K",ylab="",xlim=c(xlim_bot,xlim_top),ylim=c(-ylim_top,ylim_top))
  mtext(paste("# hegde points =",Nhedge), side = 3, line = 2, col = 1)
  mtext(paste("r-mu =",r-mu), side = 3, line = 1, col = 1)
  mtext(paste("sigma-sigma_hedge =",sigma-sigma_hedge), side = 3, line = 0, col = 1)
  mtext("Value of hedge portfolio", side = 2, line = 2, col = 1)
  title("Discrete hedging of a call-option", line = 3.25)
  garbage <- function(X){ X[abs(X)>ylim_top/1.2]=NaN;return(X) }
  if (opttype_ == 1){
    points(xlim_bot:xlim_top,garbage(pmax(xlim_bot:xlim_top - strike,0)),type='l',lwd=3) 
  }
  if (opttype_ == 2){
    points(xlim_bot:xlim_top,pmax(strike-xlim_bot:xlim_top,0),type='l',lwd=3) 
  }
  if (opttype_ == 3){
    points(xlim_bot:xlim_top,pmax(strike-xlim_bot:xlim_top,0)/abs(strike-xlim_bot:xlim_top),type='l',lwd=3) 
  }
}


for (opttype_ in c(1,3)){

  main_code <- function() {
    St<-rep(S0, length=Nrep)
    dt<-capT/Nhedge
    #initialoutlay<-BlackScholesFormula(S0,capT,strike, r,0,sigma_hedge,1,1)
    initialoutlay<-BlackScholesFormula(S0,capT,strike, r,0,sigma,opttype_,1)
    
    Vpf<-rep(initialoutlay,length=Nrep)
    
    a<-BlackScholesFormula(St,capT,strike, r,0,sigma_hedge,opttype_,2)
    b<-Vpf-a*St
    
    for(i in 2:Nhedge){
      St<-St*exp((mu-0.5*sigma^2)*dt +(sigma-i/Nhedge)*sqrt(dt)*rnorm(Nrep))
      Vpf<-a*St+b*exp(dt*r)    
      a<-BlackScholesFormula(St,(capT-(i-1)*dt),strike, r,0,sigma_hedge,opttype_,2)
      b<-Vpf-a*St
    }
    
    St<-St*exp((mu-0.5*sigma^2)*dt +sigma*sqrt(dt)*rnorm(Nrep))
    Vpf<-a*St+b*exp(dt*r)
    if (opttype_ == 1){optionpayoff<-pmax(St-strike,0)}
    if (opttype_ == 2){optionpayoff<-pmax(strike-St,0)}
    if (opttype_ == 3){optionpayoff<-pmax(strike-St,0)/abs(strike-St)}
    hedgeerror<-Vpf-optionpayoff
    
    
    # SUMMARY STATS & GRAPHS
    # ======================
    print("---------------------------------------------------")
    print(c(paste("r-mu =",r-mu),paste("sigma-sigma_hedge =",sigma-sigma_hedge),paste("# hegde points =",Nhedge)))
    
    print(paste("Initial investment =",round(initialoutlay,4)))
    print(paste("Average discounted option payoff =",round(exp(-r*capT)*mean(optionpayoff),4)))
    print(paste("Average discounted portfolio value =",round(exp(-r*capT)*mean(Vpf),4)))
    print(paste("Terminal Hedge error mean =",mean(abs(hedgeerror))))
    print(paste("Terminal Hedge error SD =",sd(hedgeerror)))
    
    ToFile<-FALSE
    #par(mfrow=c(1,2))
    if (ToFile) png(file="Scatter.png",res=300,width = 3*480, height = 3*480)
    #(mfrow=c(1,2))
    hedgeerror_plotter(hedgeerror,St,Vpf,strike,opttype_ = opttype_)
    #text(50,125,paste("# hegde points =",Nhedge),adj=0)
    #text(50,120,paste("r-mu =",r-mu),adj=0)
    #text(50,115,paste("sigma-sigma_hedge =",sigma-sigma_hedge),adj=0)
    if (ToFile) dev.off()
    reset_var() # dårlig kode praksis men meget belejligt her.
  }
  
  par(mfrow=c(1,2),mar=c(5,3,5,3))
  main_code() # standard
  
  sigma_hedge <- 0.1
  main_code()
  
  par(mfrow=c(1,2))
  mu <- 0.12
  main_code()
  
  Nhedge <- 10
  main_code()
}


