
mrank<-function(x,xdata,a){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  s1<-rep(0,d)
  for(i in 1:n){
    temp<-a[i]*(x-xdata[i,])
    normtemp<-sqrt(sum(temp^2))
    temp<-temp/normtemp
    s1<-s1+temp
  }
  s1/n
}

fsd<-function(x,xdata,l=1){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  KFS<-0
  for(i in 1:n){
    for(j in 1:n){
      kyz<-exp(-sum((xdata[i,]-xdata[j,])^2)/l^2)
      kxy<-exp(-sum((x-xdata[i,])^2)/l^2)
      kxz<-exp(-sum((x-xdata[j,])^2)/l^2)
      num<- 1+kyz-kxy-kxz
      denom<-sqrt(1+1-2*kxy)*sqrt(1+1-2*kxz)
      if(denom != 0.0)KFS<-KFS+num/denom
    }
  }
  KFS<-sqrt(KFS)/n
}

library(MASS)
n <- 100
m <- n
sx <- matrix(c(1,0.3,0.5,0.3,1,0.4,0.5,0.4,1), nrow=3)
mu_1<-matrix(c(1,1,1), nrow=3)
mu_2<-matrix(c(1,0,0), nrow=3)
nsimul<-100
misprob<-NULL
for(isimul in 1:nsimul)
{
  x <- mvrnorm(n, mu_1, sx)
  y <- mvrnorm(n, mu_2, sx)
      r1<-NULL
      for(j in 1:n){ r1<-c(r1, fsd(x[j,], x)) }
      sumr1<-r1/sum(r1)
      r2<-NULL
      for(j in 1:n){ r2<-c(r2, fsd(y[j,], y)) }
      sumr2<-r2/sum(r2)
      z<-mvrnorm(n, mu_1, sx)
      miscount<-0
      for(i in 1:m)
      {
        rx<- 1 - sqrt(sum(mrank(z[i,],x, sumr1)^2))
        ry<- 1 - sqrt(sum(mrank(z[i,],y, sumr2)^2))
	  if(rx<ry)miscount<-miscount+1

	}
  z<-mvrnorm(n, mu_2, sx)
	for(i in 1:m)
  	{
      	rx<- 1 - sqrt(sum(mrank(z[i,],x, sumr1)^2))
	  		ry<- 1 - sqrt(sum(mrank(z[i,],y, sumr2)^2))
	  if(rx>ry)miscount<-miscount+1
	}
	prob<-miscount/(2.0*m)
	misprob<-c(misprob,prob)
}
quantile(misprob)
mean(misprob)
sd(misprob)