mrank<-function(x,xdata,l=9){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  KFS<-0
  for(i in 1:n){
    for(j in 1:n){
      kyz<-exp(-sum((xdata[i,]-xdata[j,])^2)/l)
      kxy<-exp(-sum((x-xdata[i,])^2)/l)
      kxz<-exp(-sum((x-xdata[j,])^2)/l)
      num<- 1+kyz-kxy-kxz
      denom<-sqrt(1+1-2*kxy)*sqrt(1+1-2*kxz)
      KFS<-KFS+num/denom
    }
  }
  KFS<-sqrt(KFS)/n
}

library(MASS)
n <- 50
m <- n
sx <- matrix(c(1,0.3,0.5,0.3,1,0.4,0.5,0.4,1), nrow=3)
mu_1<-matrix(c(1,1,1), nrow=3)
mu_2<-matrix(c(1,0,0), nrow=3)
misprob<-NULL
for(isimul in 1:100)
{
  
      x<-mvrnorm(n, mu_1, sx)
      y<-mvrnorm(n, mu_2, sx)
      z<-mvrnorm(n, mu_1, sx)
      miscount<-0
      for(i in 1:m)
      {
        rx<- 1 - mrank(z[i,],x)
        ry<- 1 - mrank(z[i,],y)
	  if(rx<ry)miscount<-miscount+1

	}
	z<-mvrnorm(n, mu_2, sx)
	for(i in 1:m)
  	{
      		rx<- 1 - mrank(z[i,],x)
	  		ry<- 1 - mrank(z[i,],y)
	  if(rx>ry)miscount<-miscount+1
	}
	prob<-miscount/(2.0*m)
	misprob<-c(misprob,prob)
}
quantile(misprob)
mean(misprob)
sd(misprob)
