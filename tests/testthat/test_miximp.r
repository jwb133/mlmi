context("Testing mixed variables imputation functions")
library(mlmi)

test_that("Imputation of one continous variable no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- mixImp(temp, nCat=1, M=10, pd=FALSE, rseed=4423)
  }, NA)
})

test_that("Imputation of one binary variable no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- x+rnorm(n)
    x[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- mixImp(temp, nCat=1, M=10, pd=FALSE, rseed=4423)
  }, NA)
})

test_that("Imputation of continous variable with PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- mixImp(temp, nCat=1, M=10, pd=TRUE, rseed=4423)
  }, NA)
})

test_that("Imputation of binary variable with PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- x+rnorm(n)
    x[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- mixImp(temp, nCat=1, M=10, pd=TRUE, rseed=4423)
  }, NA)
})

test_that("Imputation of both variables no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 500
    x <- 1+(runif(n)<0.5)
    y <- x+rnorm(n)
    x[runif(n)<0.25] <- NA
    y[runif(n)<0.25] <- NA
    temp <- data.frame(x,y)
    imps <- mixImp(temp, nCat=1, M=10, pd=FALSE, rseed=4423)
  }, NA)
})

test_that("Restricted imputation with more variables no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 500
    x1 <- 1+(runif(n)<0.2)
    x2 <- 1+(runif(n)<0.5)
    x3 <- 1+(runif(n)<0.7)
    y <- x1+x2+x3+rnorm(n)
    x1[runif(n)<0.25] <- NA
    y[runif(n)<0.25] <- NA
    temp <- data.frame(x1,x2,x3,y)
    #specify margins for 2-way associations between categorical variables
    mymargins <- c(1,2,0,1,3,0,2,3)
    #specify design matrix for main effects only of x1 to x3
    mydesign <- matrix(c(1, 0, 0, 0,
                1, 1, 0, 0,
                1, 0, 1, 0,
                1, 1, 1, 0,
                1, 0, 0, 1,
                1, 1, 0, 1,
                1, 0, 1, 1,
                1, 1, 1, 1), byrow=TRUE, nrow=8)
    imps <- mixImp(temp, nCat=3, margins=mymargins,
                   design=mydesign, M=10, pd=FALSE, rseed=4423)
  }, NA)
})

test_that("Restricted imputation with more variables with PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 500
    x1 <- 1+(runif(n)<0.2)
    x2 <- 1+(runif(n)<0.5)
    x3 <- 1+(runif(n)<0.7)
    y <- x1+x2+x3+rnorm(n)
    x1[runif(n)<0.25] <- NA
    y[runif(n)<0.25] <- NA
    temp <- data.frame(x1,x2,x3,y)
    #specify margins for 2-way associations between categorical variables
    mymargins <- c(1,2,0,1,3,0,2,3)
    #specify design matrix for main effects only of x1 to x3
    mydesign <- matrix(c(1, 0, 0, 0,
                         1, 1, 0, 0,
                         1, 0, 1, 0,
                         1, 1, 1, 0,
                         1, 0, 0, 1,
                         1, 1, 0, 1,
                         1, 0, 1, 1,
                         1, 1, 1, 1), byrow=TRUE, nrow=8)
    imps <- mixImp(temp, nCat=3, margins=mymargins,
                   design=mydesign, M=10, pd=TRUE, rseed=4423)
  }, NA)
})


test_that("Restricted imputation using marginsType and designType no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 500
    x1 <- 1+(runif(n)<0.2)
    x2 <- 1+(runif(n)<0.5)
    x3 <- 1+(runif(n)<0.7)
    y <- x1+x2+x3+rnorm(n)
    x1[runif(n)<0.25] <- NA
    y[runif(n)<0.25] <- NA
    temp <- data.frame(x1,x2,x3,y)
    imps <- mixImp(temp, nCat=3, marginsType=1,
                   designType=1, M=10, pd=FALSE, rseed=4423)
  }, NA)
})

test_that("Restricted imputation testing marginsType default", {
  expect_error({
    set.seed(1234)
    n <- 500
    x1 <- 1+(runif(n)<0.2)
    x2 <- 1+(runif(n)<0.5)
    x3 <- 1+(runif(n)<0.7)
    y <- x1+x2+x3+rnorm(n)
    x1[runif(n)<0.25] <- NA
    y[runif(n)<0.25] <- NA
    temp <- data.frame(x1,x2,x3,y)
    imps <- mixImp(temp, nCat=3,
                   designType=1, M=10, pd=FALSE, rseed=4423)
  }, NA)
})

test_that("Restricted imputation testing designType and marginsType defaults", {
  expect_error({
    set.seed(1234)
    n <- 500
    x1 <- 1+(runif(n)<0.2)
    x2 <- 1+(runif(n)<0.5)
    x3 <- 1+(runif(n)<0.7)
    y <- x1+x2+x3+rnorm(n)
    x1[runif(n)<0.25] <- NA
    y[runif(n)<0.25] <- NA
    temp <- data.frame(x1,x2,x3,y)
    imps <- mixImp(temp, nCat=3, M=10, pd=FALSE, rseed=4423)
  }, NA)
})

test_that("Restricted imputation gives unbiased estimates when it should", {
  expect_equal({
    set.seed(1234)
    n <- 500000
    x1 <- 1+(runif(n)<0.5)
    x2 <- 1+1*(runif(n)<plogis(-x1))
    x3 <- 1+1*(runif(n)<plogis(0.5*x1-0.5*x2))
    y <- x1+x2+x3+rnorm(n)
    x2[runif(n)<plogis(x1)] <- NA
    x3[runif(n)<plogis(-0.25*x1)] <- NA
    y[runif(n)<plogis(0.4*x1)] <- NA
    temp <- data.frame(x1,x2,x3,y)
    imps <- mixImp(temp, nCat=3, M=1, pd=FALSE, rseed=4423)
    mod <- lm(y~x1+x2+x3, data=imps)
    (abs(sum((coefficients(mod)-c(0,1,1,1))^2))<0.01)
  }, TRUE)
})

test_that("Unrestricted imputation gives unbiased estimates when it should", {
  expect_equal({
    set.seed(1234)
    n <- 500000
    x1 <- 1+(runif(n)<0.5)
    x2 <- 1+1*(runif(n)<plogis(-x1))
    x3 <- 1+1*(runif(n)<plogis(0.5*x1-0.5*x2))
    y <- x1+x2+x3+rnorm(n)
    x2[runif(n)<plogis(x1)] <- NA
    x3[runif(n)<plogis(-0.25*x1)] <- NA
    y[runif(n)<plogis(0.4*x1)] <- NA
    temp <- data.frame(x1,x2,x3,y)
    imps <- mixImp(temp, nCat=3, M=1, pd=FALSE, marginsType=3, designType=2, rseed=4423)
    mod <- lm(y~x1+x2+x3, data=imps)
    (abs(sum((coefficients(mod)-c(0,1,1,1))^2))<0.01)
  }, TRUE)
})


test_that("Categorical variables must be numerics", {
  expect_error({
    set.seed(1234)
    n <- 100
    x1 <- 1+(runif(n)<0.5)
    x2 <- 1+1*(runif(n)<plogis(-x1))
    x3 <- 1+1*(runif(n)<plogis(0.5*x1-0.5*x2))
    y <- x1+x2+x3+rnorm(n)
    x2[runif(n)<plogis(x1)] <- NA
    x3[runif(n)<plogis(-0.25*x1)] <- NA
    y[runif(n)<plogis(0.4*x1)] <- NA
    temp <- data.frame(as.factor(x1),x2,x3,y)
    imps <- mixImp(temp, nCat=3, M=1, pd=FALSE, marginsType=3, designType=2, rseed=4423)
  }, NULL)
})

test_that("Categorical variables must be consecutive integers starting at 1", {
  expect_error({
    set.seed(1234)
    n <- 100
    x1 <- 0+(runif(n)<0.5)
    x2 <- 1+1*(runif(n)<plogis(-x1))
    x3 <- 1+1*(runif(n)<plogis(0.5*x1-0.5*x2))
    y <- x1+x2+x3+rnorm(n)
    x2[runif(n)<plogis(x1)] <- NA
    x3[runif(n)<plogis(-0.25*x1)] <- NA
    y[runif(n)<plogis(0.4*x1)] <- NA
    temp <- data.frame(x1,x2,x3,y)
    imps <- mixImp(temp, nCat=3, M=1, pd=FALSE, marginsType=3, designType=2, rseed=4423)
  }, NULL)
})

test_that("Variables must be consecutive integers starting at 1", {
  expect_error({
    set.seed(1234)
    n <- 100
    x1 <- 0.5+(runif(n)<0.5)
    x2 <- 1+1*(runif(n)<plogis(-x1))
    x3 <- 1+1*(runif(n)<plogis(0.5*x1-0.5*x2))
    y <- x1+x2+x3+rnorm(n)
    x2[runif(n)<plogis(x1)] <- NA
    x3[runif(n)<plogis(-0.25*x1)] <- NA
    y[runif(n)<plogis(0.4*x1)] <- NA
    temp <- data.frame(x1,x2,x3,y)
    imps <- mixImp(temp, nCat=3, M=1, pd=FALSE, marginsType=3, designType=2, rseed=4423)
  }, NULL)
})

test_that("Variables must be consecutive integers starting at 1", {
  expect_error({
    set.seed(1234)
    n <- 100
    x1 <- 2+(runif(n)<0.5)
    x2 <- 1+1*(runif(n)<plogis(-x1))
    x3 <- 1+1*(runif(n)<plogis(0.5*x1-0.5*x2))
    y <- x1+x2+x3+rnorm(n)
    x2[runif(n)<plogis(x1)] <- NA
    x3[runif(n)<plogis(-0.25*x1)] <- NA
    y[runif(n)<plogis(0.4*x1)] <- NA
    temp <- data.frame(x1,x2,x3,y)
    imps <- mixImp(temp, nCat=3, M=1, pd=FALSE, marginsType=3, designType=2, rseed=4423)
  }, NULL)
})
