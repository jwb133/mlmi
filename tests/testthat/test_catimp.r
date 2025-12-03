context("Testing categorical imputation functions")
library(mlmi)

test_that("Saturated imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- catImp(temp, M=10, pd=FALSE, type=3, rseed=4423)
  }, NA)
})

test_that("Saturated imputation with PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- catImp(temp, M=10, pd=TRUE, type=3, rseed=4423)
  }, NA)
})

test_that("Two-way association imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z)
    imps <- catImp(temp, M=10, pd=FALSE, type=1, rseed=4423)
  }, NA)
})

test_that("Two-way association imputation with PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 1000
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z)
    imps <- catImp(temp, M=10, pd=TRUE, type=1, rseed=4423)
  }, NA)
})

test_that("Three-way association imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    w <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z,w)
    imps <- catImp(temp, M=10, pd=FALSE, type=2, rseed=4423)
  }, NA)
})


test_that("Three-way association imputation PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    w <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z,w)
    imps <- catImp(temp, M=10, pd=TRUE, type=2, rseed=4423)
  }, NA)
})

test_that("Variables must be numerics", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    w <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z,as.factor(w))
    imps <- catImp(temp, M=10, pd=TRUE, type=2, rseed=4423)
  }, NULL)
})

test_that("Variables must be consecutive integers starting at 1", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 0+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    w <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z,w)
    imps <- catImp(temp, M=10, pd=TRUE, type=2, rseed=4423)
  }, NULL)
})

test_that("Variables must be consecutive integers starting at 1", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 0.5+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    w <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z,w)
    imps <- catImp(temp, M=10, pd=TRUE, type=2, rseed=4423)
  }, NULL)
})

test_that("Variables must be consecutive integers starting at 1", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- 2+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    w <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z,w)
    imps <- catImp(temp, M=10, pd=TRUE, type=2, rseed=4423)
  }, NULL)
})

