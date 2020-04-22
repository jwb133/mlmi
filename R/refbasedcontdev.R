library(MASS)
library(nlme)

expit <- function(x) {
  exp(x)/(1+exp(x))
}

n <- 50000
#we will make correlation with baseline the same to visit 1 and visit 2
corr <- matrix(1, nrow=4, ncol=4) + diag(0.5, nrow=4)
corr
data <- mvrnorm(n, mu=c(0,0,0,0), Sigma=corr)

trt <- 1*(runif(n)<0.5)

y0 <- data[,1]
y1 <- data[,2]
y2 <- data[,3]
y3 <- data[,4]

#add in effect of treatment
y1 <- y1+trt*0.5
y2 <- y2+trt*1
y3 <- y3+trt*1.5

#now make some patients dropout before visit 1
#r1=1 indicates visit 1 observed
r1 <- 1*(runif(n)<expit(1-y0))

#dropout before visit 2, based on change from y0 t0 y1
r2 <- 1*(runif(n)<expit(2-(y1-y0)))
r2[r1==0] <- 0

r3 <- 1*(runif(n)<expit(2-(y2-y0)))
r3[r2==0] <- 0

mean(r1)
mean(r2)
mean(r3)

y1[r1==0] <- NA
y2[r2==0] <- NA
y3[r3==0] <- NA

wideData <- data.frame(id=1:n, trt=trt, y0=y0, y1=y1, y2=y2, y3=y3)
