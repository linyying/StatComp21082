## -----------------------------------------------------------------------------
binorm=function(x,y){
  mu1=0;mu2=0;sigma1=1;sigma2=1;rho=0
  cons=1/(2*pi*sigma1*sigma2*sqrt(1-rho^2))
  kuohao=((x-mu1)/sigma1)^2+
    2*rho*(x-mu1)/sigma1*(y-mu2)/sigma2+
    ((y-mu2)/sigma2)^2
  den=cons*exp(-1/(2*(1-rho^2))*kuohao)
  return(den)
}
x=y=seq(-5,5,by=0.1)
z=outer(x,y,binorm)
persp(x,y,z,col="green")

## -----------------------------------------------------------------------------
library(knitr)
df=data.frame(Name=c("Zhangzhu","Yishi","Wanwu","Zaoliu","Singqi"),
              Gender=c("F","M","F","M","F"),
              Age=c(21,22,23,21,25),
              Height=c(165,168,172,182,167),
              Weight=c(42.0,49.5,61.5,72.0,55.5))
kable(df)

## -----------------------------------------------------------------------------
n=1000
u=runif(n)
sigma=0.5
x=sqrt(-2*sigma^2*log(1-u)) 
hist(x, prob = TRUE, 
     main = expression(f(x)==frac(x,sigma^2)*e^(-frac(x^2,2*sigma^2))))
y=seq(0, 2, .01)
lines(y, y*exp(-y^2/(2*sigma^2))/sigma^2)


## -----------------------------------------------------------------------------
n<-1000
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
Z<-0.75*X1+0.25*X2
hist(Z)
r<-sample(seq(0.05,0.95,length.out=19),4,replace=FALSE)
r[1];r[2];r[3];r[4]
Z1<-r[1]*X1+(1-r[1])*X2
Z2<-r[2]*X1+(1-r[2])*X2
Z3<-r[3]*X1+(1-r[3])*X2
Z4<-r[4]*X1+(1-r[4])*X2
hist(Z1);hist(Z2);hist(Z3);hist(Z4)
D=sample(c(0,1),n,replace=TRUE)
FF=D*X1+(1-D)*X2
hist(FF)
sum(D)/n

## -----------------------------------------------------------------------------
set.seed(1234567)
n=1e4;beta=3;lambda=1;r=4;t=10
N_t=rpois(n, lambda*t)
X_t=length(n)
for (i in 1:n) {
  X_t[i]=sum(rgamma(N_t[i], r, beta))
}
mean(X_t)
var(X_t)


## -----------------------------------------------------------------------------
x<-seq(0.1,0.9,length=9)
mm<-10000
u<-runif(mm)
ff<-numeric(length(x))
for (i in 1:length(x)) {
  g<-30*(x[i])^3*u^2*(1-u*x[i])^2 
  ff[i]<-mean(g)
}

zhi<-pbeta(x,3,3)
print(round(rbind(x,ff,zhi),3))

## -----------------------------------------------------------------------------
sigma=1
M<-function(x,R =10000,anti=TRUE){
  u<-runif(R/2,0,x)
  if (!anti) v<-runif(R/2) else v<-1-u
  u<-c(u,v)
  cdf<-numeric(length(x))
  for (i in 1:length(x)){
    g=u*(x[i]/sigma)^2*exp(-(u*x[i]/sigma)^2/2)
    cdf[i]<-mean(g)
  }
  cdf
}
m=1000
M1<-M2<-numeric(m)
x<-1
for (i in 1:m) {
  M1[i]=M(x,R=1000,anti=F)
  M2[i]=M(x,R=1000,anti=T)
}
v1=var(M1)
v2=var(M2)
rv=round(v2/v1,3)
rv

## -----------------------------------------------------------------------------
k=10000
est<-sd<-numeric(2)
g<-function(x) {
  exp(-1/(2*x^2))/x^4/sqrt(2*pi)*(x>0)*(x<1)
}
x<-runif(k) #f1
pp1<-g(x)/(3*exp(-x^2)/2)
est[1]<-mean(pp1)
sd[1]<-sd(pp1)

x<-runif(k) #f2
pp2<-g(x)/(2*x)
est[2]<-mean(pp2)
sd[2]<-sd(pp2)

res<-rbind(est=round(est,4),sd=round(sd,4))
colnames(res)<-paste0('f',1:2)
res 

## -----------------------------------------------------------------------------
k=10000
g<-function(x) {
  exp(-1/(2*x^2))/x^4/sqrt(2*pi)*(x>0)*(x<1)
}
x<-runif(k) #f1
pf1<-g(x)/(3*exp(-x^2)/2)
gj<-mean(pf1)
gj

## -----------------------------------------------------------------------------
intval=function(x,alpha=0.05){
  n=length(x);df=n-1 
  a=mean(x)-qt(1-alpha/2,df)*sd(x)/sqrt(n)
  b=mean(x)+qt(1-alpha/2,df)*sd(x)/sqrt(n)
  data.frame(a=a,b=b)
}
m=20
left=numeric(m);right=numeric(m)
set.seed(12333)
for(i in 1:m){ 
  x=rchisq(m,df=2)
  estimate=intval(x)
  left[i]=estimate$a
  right[i]=estimate$b
}
sum(left<2&right>2);mean(left<2&right>2)

## -----------------------------------------------------------------------------
#1
n<-20
alpha<-0.05
m<-10000 
kk1<-numeric(m)
for (j in 1:m){
  x<-rchisq(n,df=1)
  t_test<-t.test(x,mu=1)
  kk1[j]<-t_test$p.value
}
pp1<-mean(kk1<alpha)
se1<-sqrt(pp1*(1-pp1)/m)
print(c(pp1,se1))
#2
kk2<-numeric(m)
for (j in 1:m){
  x<-runif(n,0,2)
  t_test<-t.test(x,mu=1)
  kk2[j]<-t_test$p.value
}
pp2<-mean(kk2<alpha)
se2<-sqrt(pp2*(1-pp2)/m)
print(c(pp2,se2))
#3
kk3<-numeric(m)
for (j in 1:m){
  x<-rexp(n,1)
  t_test<-t.test(x,mu=1)
  kk3[j]<-t_test$p.value
}
pp3<-mean(kk3<alpha)
se3<-sqrt(pp3*(1-pp3)/m)
print(c(pp3,se3))

## -----------------------------------------------------------------------------
library(MASS)
SSig<-matrix(c(1,0,0,1),2,2)
n<-c(10, 20, 30, 50)
d=2
alpha=0.1

sk=function(n){
  for(i in 1:n){
    for(j in 1:n){
      x=mvrnorm(n, rep(0,2), SSig)
      x_bar=mean(x)
      cc=cov(x)
      bb[i,j]=(t(x[i,]-x_bar)%*%solve(cc)%*%(x[j,]-x_bar))^3
    }
  }
  bk=sum(bb)/((n)^2)
  return(bk)
}

bijie<-qchisq(1-alpha/2,df= d*(d+1)*(d+2)/6)
bijie
pow_re<-numeric(4)
m <- 100
for (i in 1:length(n)) {
  sket<- numeric(m)
  bb=matrix(0,n[i],n[i])
  for (j in 1:m) {
    sket[j] <- as.integer(n[i]*sk(n[i])/6>=bijie)
  }
  pow_re[i] <- mean(sket)
}

result <- rbind(n, pow_re)
colnames(result)=c("n1", "n2","n3", "n4")
rownames(result)=c("n","p.reject")
result 


## -----------------------------------------------------------------------------
set.seed(0)
library(boot)
library(bootstrap)
n=nrow(scor)
sgm.h=(n-1)/n*cov(scor) 
lab_h=eigen(sgm.h)$values
the_h=max(lab_h)/sum(lab_h)

ad=function(dat,i){
  x=dat[i,]
  n=nrow(scor)
  lab=eigen((n-1)/n*cov(x))$values
  the=max(lab)/sum(lab)
  return(the)
}
bst_res=boot(
  data=cbind(scor$mec,scor$vec,scor$alg,scor$ana,scor$sta),
  statistic=ad,R=2000)
the_b=bst_res$t
bias_b=mean(the_b)-the_h
se_b=sqrt(var(the_b))

the_h
bias_b
se_b

## -----------------------------------------------------------------------------
the_j=rep(0,n)
for(i in 1:n){
  x=scor[-i,]
  lam=eigen((n-1)/n*cov(x))$values
  the_j[i]=max(lam)/sum(lam)
}
bias_j=(n-1)*(mean(the_j)-the_h)
se_j=(n-1)*sqrt(var(the_j)/n)

bias_j
se_j

## -----------------------------------------------------------------------------
boot.ci(bst_res,conf=0.95,type=c('perc','bca'))


## -----------------------------------------------------------------------------
set.seed(0)
library(boot)
n=1000
ksn=function(dat,i){
  m3=mean((dat[i]-mean(dat[i]))^3)
  m2=mean((dat[i]-mean(dat[i]))^2)
  return(m3/m2^1.5)
}

da1=rnorm(n,0,1)
b1=boot(data=da1,statistic=ksn,R=2000)
print(b1)
ci1=boot.ci(b1,type=c("basic","norm","perc"))
print(ci1)
sum(ci1$normal[2]<b1$t&ci1$normal[3]>b1$t)
mean(ci1$normal[2]<b1$t&ci1$normal[3]>b1$t)

da2=rchisq(n,df=5)
b2=boot(data=da2,statistic=ksn,R=2000)
print(b2)
ci2=boot.ci(b2,type = c("basic","norm","perc"))
print(ci2)
sum(ci2$normal[2]<b2$t&ci2$normal[3]>b2$t)
mean(ci2$normal[2]<b2$t&ci2$normal[3]>b2$t)


## -----------------------------------------------------------------------------
set.seed(123456)
R<-999
l<-1:2000
x=rnorm(1000,0,1)
y=runif(1000,3,5)
z<-c(x,y)
n<-length(x)
rs<-numeric(R);t0<-cor.test(x,y,method="spearman")
for (i in 1:R) {
  ll=sample(l,size=n,replace=FALSE)
  x1=z[ll];y1=z[-ll]
  rs[i]=cor(x1,y1,method="spearman")
}
p=mean(abs(rs)>=abs(t0$estimate))
round(c(p,t0$p.value),3)

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(energy)
library(Ball)
alpha<-0.1;k<-3;R<-99;set.seed(12345)
n1<-n2<-5;n<-n1+n2;N=c(n1,n2)
ddg<-function(z,da,sizes,k){
  if(is.vector(z)) z<-data.frame(z,0);
  z<-z[da,];
  NN<-nn2(data=z,k=k+1)
  b1<-NN$nn.idx[1:n1,-1]
  b2<-NN$nn.idx[(n1+1):n,-1]
  i1<-sum(b1<n1+0.5);i2<-sum(b2>n1+0.5)
  (i1+i2)/(k*n)
}
nnq<-function(z,sizes,k){
  bh<-boot(data=z,statistic=ddg,R=R,
           sim="permutation",sizes=sizes,k=3)
  ts<-c(bh$t0,bh$t)
  kkp<-mean(ts>=ts[1])
  list(statistic=ts[1],p.value=kkp)
}
pv1<-matrix(NA,100,3)
for(i in 1:100){
  x=rnorm(5,1,1)
  y=runif(5,0,2)
  z<-c(x,y)
  pv1[i,1]<-nnq(z,N,k)$p.value
  pv1[i,2]<-eqdist.etest(z,sizes=c(n1,n2),R=R)$p.value
  pv1[i,3]<-bd.test(x=x,y=y,num.permutations=99)$p.value
}
pow1<-colMeans(pv1<alpha)
pow1


## -----------------------------------------------------------------------------
n1<-n2<-5;n<-n1+n2;N=c(n1,n2)
pv2<-matrix(NA,100,3)
for(i in 1:100){
  x=rnorm(5,1,1.5)
  y=runif(5,0,3.5)
  z<-c(x,y)
  pv2[i,1]<-nnq(z,N,k)$p.value
  pv2[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
  pv2[i,3]<-bd.test(x=x,y=y,num.permutations=99)$p.value
}
pow2<-colMeans(pv2<alpha)
pow2

## -----------------------------------------------------------------------------
pv3<-matrix(NA,100,3)
for(i in 1:100){
  x=rt(5,df=1)
  y=c(rnorm(2,0,1),rnorm(3,0,1.5))
  z<-c(x,y)
  pv3[i,1]<-nnq(z,N,k)$p.value
  pv3[i,2]<-eqdist.etest(z,sizes=c(n1,n2),R=R)$p.value
  pv3[i,3]<-bd.test(x=x,y=y,num.permutations=99)$p.value
}
pow3<-colMeans(pv3<alpha)
pow3


## -----------------------------------------------------------------------------
n1<-4;n2<-5;n<-n1+n2;N=c(n1,n2)
pv4<-matrix(NA,100,3)
for(i in 1:100){
  x=rnorm(4,0,1)
  y=c(rnorm(2),rnorm(3,0,2))
  z<-c(x,y)
  pv4[i,1]<-nnq(z,N,k)$p.value
  pv4[i,2]<-eqdist.etest(z,sizes=c(n1,n2),R=R)$p.value
  pv4[i,3]<-bd.test(x=x,y=y,num.permutations=99)$p.value
}
pow4<-colMeans(pv4<alpha)
pow4#

## -----------------------------------------------------------------------------
set.seed(12345)
f=function(x,theta,eta){
  if (any(theta<0)) return(0)
  return(1/(theta*pi*(1+((x-eta)/theta)^2)))
}

m=4000
x=numeric(m)
x[1]=rt(1,df=1)
k=0
u=runif(m)

for (i in 2:m){
  xt=x[i-1]
  y=rt(1,df=abs(xt))
  hh=f(y,1,0)*dt(xt,df=abs(y))
  dd=f(xt,1,0)*dt(y,df=abs(xt))
  if (u[i]<=abs(hh/dd)){
    x[i]=y
  }
  else{
    x[i]=xt
    k=k+1    
  }
}
k
y=x[1001:m]
a=ppoints(100)
we=tan(pi*(a-0.5))
ww=quantile(x,a)
qqplot(we,ww, main="",xlab="Cauchy Quantiles", 
       ylab="Sample Quantiles")
h<-seq(0.1,0.9,length=9)
Q=quantile(y,h)
Q
zhen=qt(h,df=1)
summ=round(rbind(Q,zhen),3)
rownames(summ)=c("MH","qt")
knitr::kable(summ)

## -----------------------------------------------------------------------------
N=15000
XY=matrix(0,N,2)
a=2
b=1
n=10
XY[1,]=c(2,0.5) 
for (i in 2:N){
  yy<-XY[i-1,2]
  XY[i,1]<-rbinom(1,n,yy)
  xx<-XY[i,1]
  XY[i,2]<-rbeta(1,xx+a,n-xx+b)
}
x<-XY[1001:N,]

colMeans(x)
cov(x)
cor(x)
plot(x)

## -----------------------------------------------------------------------------
set.seed(12345)
GR=function(psi){
  psi=as.matrix(psi)
  n=ncol(psi)
  k=nrow(psi)
  psi.means=rowMeans(psi) 
  B=n*var(psi.means)
  psi.w=apply(psi,1,"var")
  W=mean(psi.w) 
  v=W*(n-1)/n+(B/n)
  r=v/W 
  return(r)
}
cc=function(N,X1){
  x=numeric(N)
  x[1]=X1
  u=runif(N)
  for (i in 2:N){
    xt=x[i-1]
    y=rt(1,df=abs(xt))
    num=f(y,1,0)*dt(xt,df=abs(y))
    den=f(xt,1,0)*dt(y,df=abs(xt))
    if (u[i]<=abs(num/den)){
      x[i]=y
    } else{
      x[i]=xt
    }
    return(x)
  }
}
n=500
b=100
x0=c(-2,-1,1,2)
X=matrix(0,nrow=4,ncol=n)
for (i in 1:4)
  X[i,]=cc(n,x0[i])
psi=t(apply(X,1,cumsum))
for (i in 1:nrow(psi))
  psi[i,]=psi[i,]/(1:ncol(psi))
print(GR(psi))


## -----------------------------------------------------------------------------
set.seed(12345)
N=5000
XY=matrix(0,N,2)
a=2;b=1;n=100
jc=function(XY,N){
  for (i in 2:N){
    yy<-XY[i-1,2]
    XY[i,1]<-rbinom(1,n,yy)
    xx<-XY[i,1]
    XY[i,2]<-rbeta(1,xx+a,n-xx+b)
  }
  return(XY)
}
psi=t(apply(jc(XY,N), 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,]=psi[i,]/(1:ncol(psi))
print(GR(psi))

## -----------------------------------------------------------------------------
f=function(k,a){
  d=nrow(a)
  l1=(-1/2)^k/factorial(k)
  l2=(norm(a,type="F"))^(2*k+2)/(2*k+1)/(2*k+2)
  l3=exp(lgamma((d+2)/2)+lgamma(k+1.5)-lgamma(k+d/2+1))
  l=l1*l2*l3
  return(l)
}
ss=numeric(length=0)
sumf=function(a){
ss[1]=f(1,a)
i=2
repeat{
  ss[i]=f(i,a)
  if(abs(ss[i]-ss[i-1])<=1e-6) break
  else i=i+1
}
su=sum(ss)
return(su)
}
a=matrix(c(1,2),ncol=1)
round(sumf(a),3) 

## -----------------------------------------------------------------------------
rsq=function(mod) summary(mod)$r.squared
##exercise 3
formulas=list(
  mpg~disp,
  mpg~I(1/disp),
  mpg~disp+wt,
  mpg~I(1/disp)+wt
)
kk3=lapply(1:4,function(i) {lm(formulas[[i]],data=mtcars)})
rt3=matrix(0,4,2)
dimnames(rt3)[[2]]=c("i","r.squared")
for(i in 1:4){
  rt3[i,1]=i
  rt3[i,2]=rsq(kk3[[i]])
}
rt3
##exercise 4
bootstraps=lapply(1:10, function(i) {
  rows=sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
kk4=lapply(1:10,function(i) {lm(mpg~disp,data=bootstraps[[i]])})
rt4=matrix(0,10,2)
dimnames(rt4)[[2]]=c("i","r.squared")
for(i in 1:10){
  rt4[i,1]=i
  rt4[i,2]=rsq(kk4[[i]])
}
rt4

## -----------------------------------------------------------------------------
h1=rnorm(50,1,2)
h2=rt(50,df=3)
hu1=data.frame(h1,h2)
vapply(hu1,sd,numeric(1))
h12=c(h1,h2)
h3=rchisq(100,df=2)
hu2=data.frame(h12,h3)
vapply(hu2,sd,numeric(1))

## -----------------------------------------------------------------------------
mcvapply=function(X,fun,simplify=T){
  answer=lapply(X=as.list(X),FUN=fun)
  if (is.character(X)&&is.null(names(answer))) 
    names(answer)=X
  if (!isFALSE(simplify)) 
    simplify2array(answer,higher=(simplify=="array"))
  else answer
}
mcvapply(mtcars,mean)

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
cppFunction('NumericMatrix fdeC(int a,int b,int N,int x,int y){
  NumericMatrix GX(N,2);
  double X0,X1;
  GX(0,0)=x;GX(0,1)=y;
  for(int i=1;i<N;i++){
  X1=GX(i-1,1);
  GX(i,0)=rbinom(1,25,X1)[0];
  X0=GX(i,0);
  GX(i,1)=rbeta(1,X0+a,25-X0+b)[0];
  }
  return(GX);
}')
fdeR=function(a,b,N,x,y){
  GX=matrix(0,N,2)
  GX[1,]=c(x,y)
  for(i in 2:N){
    X2=GX[i-1,2]
    GX[i,1]=rbinom(1,25,X2)
    X1=GX[i,1]
    GX[i,2]=rbeta(1,X1+a,25-X1+b)
  }
  return(GX)
}
qqplot(fdeR(1,1,1000,0,0.5)[,1],fdeR(1,1,1000,0,0.5)[,2])
qqplot(fdeC(1,1,1000,0,0.5)[,1],fdeC(1,1,1000,0,0.5)[,2])
microbenchmark(fdeR=fdeR(1,1,1000,0,0.5),fdeC=fdeC(1,1,1000,0,0.5))
ts=microbenchmark(fdeR=fdeR(1,1,1000,0,0.5),fdeC=fdeC(1,1,1000,0,0.5))
summary(ts)[,c(1,3,5,6)]
microbenchmark(fdeR=fdeR(1,10,1000,0,0.5),fdeC=fdeC(1,10,1000,0,0.5))
ts=microbenchmark(fdeR=fdeR(1,10,1000,0,0.5),fdeC=fdeC(1,10,1000,0,0.5))
summary(ts)[,c(1,3,5,6)]

