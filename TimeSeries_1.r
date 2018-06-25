########## R codes for <Financial Time Series Analysis>, Tsay

########## CHAPTER 2  Linear Time Series Analysis and Application
########## 'FinTS' library is now not officailly available and has to be uploaded by hand.

# R Codes
library('FinTS')
help (package='FinTS')

mibm=read.table(file="data1.txt")
mvw=read.table(file="data2.txt")

par(mfrow=c(2,1))

acf(mibm,main="ACF of IBM simple return")    #plotting ACF of mibm time Series
acf(log(1+mibm),main="ACF of IBM log return")

acf(mvw,main="ACF of Value_weighted index simple return") 
acf(log(1+mvw),main="ACF of Value_weighted index simple return") 

Box.test(mibm,lag=5,type="Ljung")           
#if Box-Ljung test is rejected, time series does not have significant autocorrelation


##### AutoRegression Model AR(p)
pvw=pacf(mvw,lag=10)
pvw
round(pvw$acf,2)
# pacf implies that AR model should be considered

arvw=ar(mvw,order=10)
arvw$aic             # aic here has been adjusted

### estimating parameters in AR(3) model
mvw=ts(mvw)
armvw3=ar.ols(mvw,order=3,demean=F,intercept=T)         # "demean" means adjusting the data by mining means
armvw3
armvw3$asy.se.coef                  
summary(armvw3)                      # return estimated parameters, including intercept and coefficients

### test: whether AR(3) model fits data well or not
# test residules
Box.test(armvw3$resid, lag=10, type='Ljung')



### prediction
length(mvw)
armvw3to858=ar.ols(mvw[1:858],aic=F,order=3,demean=F,intercept=T)
vwp=predict(armvw3to858, n.ahead=6)
vwp

### plotting prediction result
# "U" is upper bound with 0.95 significance level. "L" is lower bound.
U = vwp$pred + 1.96*vwp$se         # upper bound= prediction value + 1.96*standard deviation
L = vwp$pred - 1.96*vwp$se   
par(mfrow=c(1,1))
ts.plot(mvw,vwp$pred,col=1:2, xlim=c(840,864))

lines(U, col="blue", lty="dashed")
lines(L, col="blue", lty="dashed")



#===========================================================================================

##### Moving Average Models

### AR: ar.ols(data,order=3)
### MA: arima(data,order=c(0,0,3),method="ML")        "ML"=maximum likelihood

acf(mvw)
fit = arima(mvw,order=c(0,0,9),method="ML")
Box.test(fit1$resid,10)
1-pchisq(0.1552,1)

### Box-Ljung test can also be applied in following way
library('FinTS')
fit2 = ARIMA(mvw,order=c(0,0,9),method="ML")         # ARIMA in FinTS library
fit2$Box.test                                        # with no need to adjust degree of freedom
tsdiag(fit2)                                         # return plots


### prediction
fit = arima(mvw,order=c(0,0,9),method="ML")
pre=predict(fit,n.ahead=3)




#===========================================================================================

##### ARMA, ARIMA

### simulating ARIMA data
ts.sim = arima.sim(list(order=(1,0,1), ar=0.9, ma=0.8), n=200)
par(mfrow=c(3,1))
acf(ts.sim)
pacf(ts.sim)

### ARMA analysis of 3M corporate
dat=read.table(file="data.txt"）
mm=dat[,2]
ates=as.Data(as.character(dat[,1],format="%Y%m%d))
plot(dates,mm,type="1")
acf(mm)
pacf(mm)

### apply AIC to select order
# model with smallest AIC should be considered
# AIC: take account of negative numbers 
arima(mm, order=c(0,0,0), method="ML")
arima(mm, order=c(1,0,0), method="ML")
arima(mm, order=c(2,0,0), method="ML")
arima(mm, order=c(0,0,1), method="ML")
arima(mm, order=c(0,0,2), method="ML")
arima(mm, order=c(1,0,1), method="ML")
arima(mm, order=c(1,0,2), method="ML")

### test residules to check whether the model fits well
arma1=arima(mm, order=c(1,0,1), method="ML")
Box.test(arma1$resid, 10, type="Ljung")
1-pchisq(8.0674,8)

### prediction
predict(arma1, n.ahead=3)


#================================================================

##### ARIMA 
arima(log(1+mm), order=c(1,0,1), method="ML")
# coefficients in MA and AR are found insignificant. c(0,0,0) is then considered.

arima(log(1+mm), order=c(0,0,0), method="ML")

mu= mean(log(1+mm))
p=cumsum(log(1+mm))           
pstar=cumsum(log(1+mm)-mu)

par(mfrow=c(1,1))
plot(p)



## Example 2.2
dat=read.table(file="data.txt")
gdp=dat[,3]
par(mfrow=c(2,2))
plot(gdp,type="1")               
dgdp=diff(gdp)
plot(dgdp,type="1")                

acf(gdp)
pacf(gdp)

library(urca)
summary(ur.df(gdp,type="trend",selectlags="AIC"))
summary(ur.kpss(gdp))               


library(fUitRoots)
help(package="fUitRoots")
urdfTest(gdp,type="ct")


#### Example 2.2
# data: S&P 500
dat=read.table(file="table.txt")
sp=dat[,2]
plot(sp,type="1")

urdfTest(sp,type="ct")
urdfTest(sp,type="c")
urdfTest(sp,type="nc")

# "value of test" : 44.3803 > 8.43



#=================================================================

##### Season Model

jj=read.table(file="data.txt")
jj=ts(jj)
plot(jj,type="o")                  
plot(log(jj),type="o")                

xt=log(jj)
par(mfrow=c(2,2))
acf(xt,main="x")

dxt=diff(xt)
acf(dxt,main="dx")

dsxt=diff.ts(xt,lag=4)                
acf(dsxt)             


dxds=diff.ts(dxt,lag=4)               
acf(dxds)



#====================================================================

dat=read.table(file="data.txt")
r1=ts(dat[,1])
r3=ts(dat[,2])

plot(r1,col="red",type="1")
lines(r3,col="blue",lty="dashed")


plot(r1,r3)      

lm.r=lm(r3~r1) 
summary(lm.r)

plot(lm.r$resid,type="1")
acf(lm.r$resid)

# r1t=alpha+beta*r2t + et



c1=diff(r1)
c3=diff(r3)

par(mfrow=c(1,2))
plot(r1,r3)
plot(c1,c3)

lm.c=lm(c3~c1)
summary(lm.c)
plot(lm.c$resid, type="1")
acf(lm.c$resid)

fit3=arima(c3,xreg=c1,order=c(0,0,1))
acf(fit3$resid)
