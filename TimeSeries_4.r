'''R Codes for <Financial Time Series Analysis>, Tsay'''

##### VAR Model
dat=read.table(file="m-fac9003.txt") 
colnames(dat)=c("AA","AGE","CAT","F","FDX","GM","HPQ","KMB","MEL","NYT","PG","TRB","TXN","SP5")
round(apply(dat,2,mean),2)
round(apply(dat,2,sd),2)  
G=cbind(rep(1,168),dat[,14])
R=as.matrix(dat[,1:13])
ab.hat=solve(t(G)%*%G)%*%(t(G)%*%R)
alpha.hat=ab.hat[1,]
beta.hat=ab.hat[2,]
E.hat=R-G%*%ab.hat
D.hat=diag(crossprod(E.hat)/(168-2))
r.square=1-(168-2)*D.hat/diag(crossprod(R))
t(rbind(beta.hat,sqrt(D.hat),r.square)) 
barplot(beta.hat,main="Beta")
barplot(r.square,main="R-square")
##
summary(lm(R[,1]~dat[,14]))
summary(lm(R[,2]~dat[,14]))
summary(lm(R[,3]~dat[,14]))
 
###
cov.r=(beta.hat*var(dat[,14]))%*%t(beta.hat)+diag(D.hat) 
sd.r=sqrt(diag(cov.r))
corr.r=cov.r/outer(sd.r,sd.r)
print(corr.r,digits=1,width=2)
print(round(cor(R),1),digits=1,width=2)

###GMVP
w.gmin.model=solve(cov.r)%*%rep(1,nrow(cov.r))
w.gmin.model=w.gmin.model/sum(w.gmin.model)
t(w.gmin.model)

w.gmin.data=solve(var(R))%*%rep(1,nrow(cov.r))
w.gmin.data=w.gmin.data/sum(w.gmin.data)
t(w.gmin.data)

resi.cov=t(E.hat)%*%E.hat/(168-2)
resi.sd=sqrt(diag(resi.cov))
resi.cor=resi.cov/outer(resi.sd,resi.sd)
print(round(resi.cor,2),digits=1,width=2) 

##############
da=read.table(file="m-cpice16-dp7503.txt")
cpi=da[,1]
cen=da[,2]
x1=cbind(cpi,cen)
y1=data.frame(x1)
library(vars)
VARselect(y1)
var3.fit=VAR(y1,p=3)

res=residuals(var3.fit)[166:333,]
G2=cbind(rep(1,168),res)
R=as.matrix(dat[,1:13])
ab.hat=solve(t(G2)%*%G2)%*%(t(G2)%*%R)
alpha.hat=ab.hat[1,]
beta.hat=ab.hat[2:3,]
E.hat=R-G2%*%ab.hat
D.hat=diag(crossprod(E.hat)/(168-2))
r.square=1-(168-2)*D.hat/diag(crossprod(R))
t(rbind(beta.hat,sqrt(D.hat),r.square)) 
par(mfrow=c(1,3))
barplot(beta.hat[1,],main="Beta of CPI",horiz=T)
barplot(beta.hat[2,],main="Beta of CEN",horiz=T)
barplot(r.square,main="R-square",horiz=T)

cov.r=(beta.hat*var(dat[,14]))%*%t(beta.hat)+diag(D.hat) 
sd.r=sqrt(diag(cov.r))
corr.r=cov.r/outer(sd.r,sd.r)
print(corr.r,digits=1,width=2)
print(round(cor(R),1),digits=1,width=2)

resi.cov=t(E.hat)%*%E.hat/(168-2)
resi.sd=sqrt(diag(resi.cov))
resi.cor=resi.cov/outer(resi.sd,resi.sd)
print(round(resi.cor,2),digits=1,width=2) 

################Fundamental model
##BARRA MODEL
da=read.table(file="m-barra-9003.txt")
colnames(da)=c("AGE","C","MWD","MER","DELL","HPQ","IBM","AA","CAT","PG")
round(apply(da,2,mean),2)
round(apply(da,2,sd),2)      
rmean=apply(da,2,mean)
R.rm=da-rmean  
fin=c(rep(1,4),rep(0,6))
tech=c(rep(0,4),rep(1,3),rep(0,3))
oth=c(rep(0,7),rep(1,3))
ind.dum=cbind(fin,tech,oth)  ###beta
ind.dum

cov.R=var(R.rm)
sd.R=sqrt(diag(cov.R))
corr.R=cov.R/outer(sd.R,sd.R)
print(corr.R)

F.hat.o=solve(crossprod(ind.dum))%*%(t(ind.dum)%*%t(R.rm))
E.hat.o=R.rm-ind.dum%*%F.hat.o
diagD.hat.o=diag(apply(E.hat.o,2,var))

Dinv.hat=solve(diagD.hat.o)

H1=t(ind.dum)%*%Dinv.hat%*%ind.dum
Hmtx=solve(H1)%*%t(ind.dum)%*%Dinv.hat
F.hat.g=Hmtx%*%t(R.rm)
par(mfrow=c(3,1))
plot(F.hat.g[1,],type="l")
plot(F.hat.g[2,],type="l")
plot(F.hat.g[3,],type="l")


F.hat.gt=t(F.hat.g)
E.hat.g=R.rm-ind.dum%*%F.hat.g
diagD.hat.g=apply(E.hat.g,2,var)
t(Hmtx)

cov.ind=ind.dum%*%var(F.hat.gt)%*%t(ind.dum)+diag(diagD.hat.g)
sd.ind=sqrt(diag(cov.ind))
corr.ind=cov.ind/outer(sd.ind,sd.ind)
print(corr.ind)

glm(t(R.rm)[,1]~ind.dum[,1]+ind.dum[,2]+ind.dum[,3])

###### PCA
dat=read.table(file="m-5cln.txt")
par(mfrow=c(5,1))
colnames(dat)=c("IBM","HP","INTC","MER","MWD")
plot(dat[,1],type="l",main="IBM")
plot(dat[,2],type="l",main="HPQ")
plot(dat[,3],type="l",main="INTC")
plot(dat[,4],type="l",main="MER")
plot(dat[,5],type="l",main="MWD")

apply(dat,2,mean)
round(cov(dat),2) 
round(cor(dat),2)
pca1=princomp(dat)
summary(pca1,loadings=T)
loadings(pca1)
eigen(cov(dat))
screeplot(pca1,type="l")
biplot(pca1)

pca2=princomp(dat,cor=T)
summary(pca2,loadings=T)
eigen(cor(dat))
screeplot(pca2,type="l")
biplot(pca2)

###########
###EXAMP 9.2 of Tsay's book
dat=read.table(file="m-5cln.txt") 
colnames(dat)=c("IBM","HP","INTC","MER","MWD")
mqtest(dat,8)  

f1=factor.pca(cor(dat),2)    
varimax(f1$loadings,normalize=F) 


factanal(dat,factors=2,rotation="none")
factanal(dat,factors=2,rotation="varimax")


##EXAMPLE 9.3 OF Tsay's book
bnd=read.table(file="m-bnd.txt")
colnames(bnd)=c("30y","20y","10y","5y","1y")
round(cor(bnd),2)
mqtest(bnd,5)

f2=factor.pca(cor(bnd),2)
varimax(f2$loadings,normalize=F) 

factanal(bnd,factors=2,rotation="none")
factanal(bnd,factors=2,rotation="varimax")
f3=factanal(bnd,factors=2,rotation="varimax",scores="regression")
f3$scores
plot(f3$scores,type="n")
text(f3$scores[,1],f3$scores[,2])

####EXAMPLE 9.4
rtn=read.table(file="m-barra-9003.txt")
colnames(rtn)=c("AGE","C","MWD","MER","DELL","HPQ","IBM","AA","CAT","PG")
rtn.fac2=factanal(rtn,factors=2)
rtn.fac2

rtn.fac3=factanal(rtn,factors=3)
rtn.fac3
par(mfrow=c(3,1)) 
barplot(loadings(rtn.fac3)[,1],main="factor1") 
barplot(loadings(rtn.fac3)[,2],main="factor2")
barplot(loadings(rtn.fac3)[,3],main="factor3")

factanal(rtn,factors=3,rotation="varimax")