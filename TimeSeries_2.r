################ ARCH and GARCH                         

setwd("F:/data/ch2-3")
dat=read.table(file="m-intc7303.txt")              
intc=dat[,2]                                        
intc=ts(log(1+intc),frequency=12,start=c(1973,1))   

par(mfrow=c(2,2))
acf(intc)
pacf(intc^2)

library(FinTS)
help(package="FinTS") 
AutocorTest(intc,12)      
Box.test(intc,12,type="Ljung") 
ArchTest(intc,12)        
### Small p-value in Box-Ljung test and large F-value in LM test imply significance of ARCH effect.

######
exch=ts(read.table(file="exch-perc.txt"))
par(mfrow=c(2,1))
plot(exch) 
plot(exch^2)
acf(exch)
pacf(exch)
pacf(exch^2)

################## ARCH
####Example3.1
   
'''
ARCH(m):
1. Plot data, data^2
2. Box-Ljung test --> small p-value implies no autocorrelation;
   LM test --> large F-value fimplies ARCH effect
3. ACF --> no significant autocorrelation.
       --> order of AR(p) model
   PACF --> order of MA(q) model
   PACF(data^2) --> m order in ARCH model
4. garchOxFit(formula.mean=~(p,q),formula.var=~garch(0,m),series=intc)
   data: ARMA(p,q)
   residules: GARCH(0,m) --> arch(m)=garch(0,m)
5. test: standardlized residules=residules/condvar
   Box-Ljung test of (standardlized residules) -->  small p-value implies no autocorrelation
   Box-Ljung test of (standardlized residules^2) -->  small p-value implies no arch effect
   -->> ARCH model fits data well.
'''


# ACF and PACT(data^2) are to test ARCH effect
# pacf(intc^2) implies ARCH(3) 
acf(intc)
pacf(intc^2)

library(fSeries)       
m3=garchOxFit(formula.mean=~arma(0,0),formula.var=~garch(0,3),series=intc)
m1=garchOxFit(formula.mean=~arma(0,0),formula.var=~garch(0,1),series=intc)
names(m1)
par(mfcol=c(2,1))
plot(intc,type="l")

### test
plot(sqrt(m1$condvars),type="l")
par(mfcol=c(1,1))
sresi=m1$resid/sqrt(m1$condvars) 
acf(sresi)
pacf(sresi^2)
qqnorm(sresi)
qqline(sresi)   

####t-inovatioin
tm1=garchOxFit(formula.mean=~arma(0,0),formula.var=~garch(0,1),series=intc,cond.dist="t") 
sresi=tm1$resid/sqrt(tm1$condvars)
pacf(sresi^2)
qqplot(rt(10000,6.1),sresi)
qqline(sresi)
Box.test(sresi,12,type="Ljung")
Box.test(sresi^2,12,type="Ljung")


####Example 3.2
plot(exch)
acf(exch) 
pacf(exch^2) 
m3=garchOxFit(formula.mean=~arma(0,0),formula.var=~garch(0,3),series=exch)    

#==========================================================================
#####Garch Modle
'''
GARCH:
1. try: GARCH(1,1),GARCH(1，2),GARCH(2，1)
'''
##Example 3.3
sp=as.matrix(read.table(file="sp500.dat"))
par(mfrow=c(3,1))
plot(sp,type="l")
acf(sp)
pacf(sp)
pacf(sp^2)
m1=arima(sp,order=c(0,0,3))
m1
sp=ts(sp)
m2=garchOxFit(formula.mean=~arma(3,0),formula.var=~garch(1,1),series=sp)
# AR order is found to be insignificant

m3=garchOxFit(formula.mean=~arma(0,0),formula.var=~garch(1,1),series=sp)
sresi=m3$resid/sqrt(m3$condvars)
plot(sresi)
acf(sresi)
pacf(sresi^2)
Box.test(sresi,12,type="Ljung")
Box.test(sresi,24,type="Ljung")

#### t innovation
m3=garchOxFit(formula.mean=~arma(0,0),formula.var=~garch(1,1),series=sp,cond.dist="t")


#==================================================================================
#################IGARCH
sp=as.matrix(read.table(file="sp500.dat"))
igarch=garchOxFit(formula.mean=~arma(0,0),formula.var=~igarch(1,1),series=sp)

###GARCH-M
garchm=garchOxFit(formula.mean=~arma(0,0),formula.var=~garch(1,1),series=sp,arch.in.mean=1)###方差 
garchm=garchOxFit(formula.mean=~arma(0,0),formula.var=~garch(1,1),series=sp,arch.in.mean=2)###标准差 

#######E-GARCH
###1926-1997
ibm=read.table(file="m-ibm2697.txt")
ibmln=as.matrix(log(ibm+1))
egarch=garchOxFit(formula.mean=~arma(1,0),formula.var=~egarch(1,1),series=ibmln)
####1926-2003
dat=read.table(file="m-ibm3dx2603.txt")
ibm03=dat[,2]
ibmln03=as.matrix(log(ibm03+1))
egarch=garchOxFit(formula.mean=~arma(0,0),formula.var=~egarch(1,1),series=ibmln03,cond.dist="ged")
########T-GARCH
tgarch=garchOxFit(formula.mean=~arma(0,0),formula.var=~gjr(1,1),series=ibmln03,cond.dist="ged")

############
##EXAMPLE 3.4
dat=read.table(file="m-ibmsplnsu.dat")
ibm=dat[,1]
sp5=dat[,2]
u=dat[,3]
par(mfrow=c(2,1))
plot(ibm,type="l")
plot(sp5,type="l") 
m3=garchOxFit(formula.mean=~arma(1,0),formula.var=~garch(1,1),series=ibm)
library(fGarch)
help(package="fGarch")


####EXAMPLE 3.5
m5=garchOxFit(formula.mean=~arma(0,0),formula.var=~garch(2,1),series=sp5)