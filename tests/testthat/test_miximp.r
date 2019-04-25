context("Testing mixed variables imputation functions")
library(mlmi)

test_that("Unrestricted imputation of continous variable no PD draw runs", {
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

test_that("Unrestricted imputation of binary variable no PD draw runs", {
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

test_that("Unrestricted imputation of continous variable with PD draw runs", {
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

test_that("Unrestricted imputation of binary variable with PD draw runs", {
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

test_that("Unrestricted imputation of both variables no PD draw runs", {
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

test_that("Restricted imputation no PD draw runs", {
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

test_that("Restricted imputation with PD draw runs", {
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

