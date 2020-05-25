mrank<-function(x,xdata){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  s1<-rep(0,d)
  for(i in 1:n){
    temp<-x-xdata[i,]
    normtemp<-sqrt(sum(temp^2))
    temp<-temp/normtemp
    s1<-s1+temp
  }
  s1/n
}

library(MASS)
n <- 100
m <- n
sx <- matrix(c(1,0.3,0.5,0.3,1,0.4,0.5,0.4,1), nrow=3)
mu_1<-matrix(c(1,1,1), nrow=3)
mu_2<-matrix(c(1,0,0), nrow=3)
misprob<-NULL
for(isimul in 1:100)
{

  x<-mvrnorm(n, mu_1, sx)
  y<-mvrnorm(n, mu_2, sx)
 	xs<-var(x)
	xs1<-chol(xs)
	xsigma<-solve(xs1)
	ys<-var(y)
	ys1<-chol(ys)
	ysigma<-solve(ys1); x<-x%*%xsigma; y<-y%*%ysigma
  z<-mvrnorm(n, mu_1, sx); z1<-z%*%xsigma; z2<-z%*%ysigma
      miscount<-0
      for(i in 1:m)
      {
	  		rx<- 1 - sqrt(sum(mrank(z1[i,],x)^2))
	  		ry<- 1 - sqrt(sum(mrank(z2[i,],y)^2))
	  if(rx<ry)miscount<-miscount+1

	}
  z<-mvrnorm(n, mu_2, sx); z1<-z%*%xsigma; z2<-z%*%ysigma
	for(i in 1:m)
  	{
      		rx<- 1 - sqrt(sum(mrank(z1[i,],x)^2))
      		ry<- 1 - sqrt(sum(mrank(z2[i,],y)^2))
	  if(rx>ry)miscount<-miscount+1
	}
	prob<-miscount/(2.0*m)
	misprob<-c(misprob,prob)
}
quantile(misprob)
mean(misprob)
sd(misprob)
