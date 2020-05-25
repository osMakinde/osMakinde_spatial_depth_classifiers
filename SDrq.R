
mrank<-function(x,xdata,c=1.5){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  KFS<-0
  for(i in 1:n){
    for(j in 1:n){
      kyz<-1 - sum((xdata[i,]-xdata[j,])^2)/(c+sum((xdata[i,]-xdata[j,])^2))
      kxy<-1 - sum((x-xdata[i,])^2)/(c+sum((x-xdata[i,])^2))
      kxz<-1 - sum((x-xdata[j,])^2)/(c+sum((x-xdata[j,])^2))
      kxx<-1 - sum((x-x)^2)/(c+sum((x-x)^2))
      kyy<-1 - sum((xdata[i,]-xdata[i,])^2)/(c+sum((xdata[i,]-xdata[i,])^2))
      kzz<-1 - sum((xdata[j,]-xdata[j,])^2)/(c+sum((xdata[j,]-xdata[j,])^2))
      num<- kxx+kyz-kxy-kxz
      if(num > 0 && kxx+kyy-2*kxy > 0 && kxx+kzz-2*kxz > 0)KFS<-KFS + num/(sqrt(kxx+kyy-2*kxy)*sqrt(kxx+kxz-2*kxz))
      
    }
  }
  KFS<-1 - sqrt(KFS)/n
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
      m<-n
      z<-mvrnorm(n, mu_1, sx)
      miscount<-0
      for(i in 1:m)
      {
      	rx<- mrank(z[i,],x)
	  		ry<- mrank(z[i,],y)
	      if(rx<ry)miscount<-miscount+1
	    }
      z<-mvrnorm(n, mu_2, sx)
	    for(i in 1:n)
  	  {
      	rx<- mrank(z[i,],x)
	  		ry<- mrank(z[i,],y)
	      if(rx>ry)miscount<-miscount+1
	    }
	prob<-miscount/(2.0*n)
	misprob<-c(misprob,prob)
}
quantile(misprob)
mean(misprob)
sd(misprob)
