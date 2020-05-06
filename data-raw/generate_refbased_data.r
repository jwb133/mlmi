#this code generates a simulated dataset to illustrate the reference based continuous imputation methods

set.seed(1234)
n <- 500
corr <- matrix(1, nrow=5, ncol=5) + diag(0.5, nrow=5)
corr[,1] <- c(1.5, 1, 0.75, 0.5, 0.25)
corr[1,] <- corr[,1]
corr
data <- mvrnorm(n, mu=c(0,0,0,0,0), Sigma=corr)

trt <- 1*(runif(n)<0.5)
v <- data[,1]
y0 <- data[,2]
y1 <- data[,3]
y2 <- data[,4]
y3 <- data[,5]

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

y1[r1==0] <- NA
y2[r2==0] <- NA
y3[r3==0] <- NA

ctsTrialWide <- data.frame(id=1:n, trt=trt, v=v, y0=y0, y1=y1, y2=y2, y3=y3)

usethis::use_data(ctsTrialWide)
