'''
R codes for Chapter4 in <Financial Time Series Analysis>, Tsay
'''

####### Multi Time Series Analysis
####
###Example 8.1
dat=read.table(file="m-ibmsp2699.txt")
ibm=dat[,2]
sp=dat[,3]
ibmln=log(ibm+1)*100 
spln=log(sp+1)*100
ibmln=ts(ibmln)      
spln=ts(spln)
par(mfrow=c(2,1))
plot(ibmln,type="l")
plot(spln,type="l")

ibmln1=lag(ibmln,-1)
spln1=lag(spln,-1)
par(mfrow=c(2,2))
plot(spln,ibmln,main="ibm vs sp")  
plot(spln1,ibmln,main="ibm vs sp(1)")
plot(ibmln1,spln,main="sp vs ibm(1)")
plot(spln1,spln,main="sp vs sp(1)")


library(FinTS)
FinTS.stats(ibmln)
FinTS.stats(spln)


c.cor(ibmln,spln,l=1)#
c.cor(ibmln,spln,l=2)
c.cor(ibmln,spln,l=3)
c.cor(ibmln,spln,l=4)
c.cor(ibmln,spln,l=5)

par(mfrow=c(3,1))
acf(ibmln)
acf(spln)
cc=ccf(ibmln,spln,ylab="cross-correlation") 
c2=ccf(spln,ibmln,ylab="cross-correlation") 

bnd=read.table(file="m-bnd.txt")
b30=bnd[,1]
b20=bnd[,2]
b10=bnd[,3]
b5=bnd[,4]
b1=bnd[,5]

B=cbind(b30,b20,b10,b5,b1)
apply(B,2,mean)
apply(B,2,sd)  
cov(B)
round(cor(B),2)

par(mfrow=c(5,1))
plot(b30,type="l",main="30 years of bond")
plot(b20,type="l",main="20 years of bond")
plot(b10,type="l",main="10 years of bond")
plot(b5,type="l",main="5 years of bond")
plot(b1,type="l",main="1 years of bond")


###########################

x=cbind(ibmln,spln)
mqtest(x,1)
mqtest(x,5)
mqtest(x,10)
###bond of mqtest
Y=cbind(b30,b20,b10,b5,b1)
mqtest(Y,5)


######## VAR
'''
VAR:

ar：
 x.intercept 
      【1】          【2】
【1】IBM-IBM(-3)   IBM-SP
【2】SP-IBM(-3)    SP-SP(-3)

VAR：
VARselect(x) 
summary(fit3)
'''

fit1=ar(x,method="ols")
fit1$aic
mqtest(fit1$resid[-1:-3,],10)
1-pchisq(40.786,28)
fit2=ar(x,method="yule-walker")
##
library(vars)
VARselect(x)   
fit3=VAR(x,p=5) 
summary(fit3)
serial.test(fit3,10)  
predict(fit3,n.ahead=10,ci=0.95)


VARselect(x)   
fit4=VAR(x,p=1)
summary(fit4)
resi=residuals(fit4)
sdd=sd(resi)
mea=apply(resi,2,mean)
U1=mea[1]+2*sdd[1]  
L1=mea[1]-2*sdd[1]   
U2=mea[2]+2*sdd[2]
L2=mea[2]-2*sdd[2] 
par(mfrow=c(2,1))
plot(resi[,1],type="l")
abline(h=U1,col="blue") 
abline(h=L1,col="blue") 
plot(resi[,2],type="l")
abline(h=U2,col="blue")
abline(h=L2,col="blue") 
serial.test(fit4,15)  
pred=predict(fit4,n.ahead=6)
plot(pred)
plot(pred,xlim=c(800,900))

fit4.irf=irf(fit4,n.ahead=6,boot=F)
par(mfrow=c(4,4))
plot(fit4.irf)
fit4.irf1=irf(fit4,impulse="y1",response="y1",n.ahead=6,boot=F)
plot(fit4.irf1)
plot(irf(fit4,impulse="y1",response="y2",n.ahead=6,boot=F))
plot(irf(fit4,impulse="y2",response="y1",n.ahead=6,boot=F))
plot(irf(fit4,impulse="y2",response="y2",n.ahead=6,boot=F))

########VMA

library(fMultivar)
help(package="fMultivar")   