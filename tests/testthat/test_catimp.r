context("Testing categorical imputation functions")
library(mlmi)

test_that("Saturated imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    norm::rngseed(123126)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- catImp(temp, M=10, pd=FALSE, type=3)
  }, NA)
})

test_that("Saturated imputation with PD draw runs", {
  expect_error({
    set.seed(1234)
    norm::rngseed(123126)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- catImp(temp, M=10, pd=TRUE, type=3)
  }, NA)
})

test_that("Two-way association imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    norm::rngseed(123126)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z)
    imps <- catImp(temp, M=10, pd=FALSE, type=1)
  }, NA)
})

test_that("Two-way association imputation with PD draw runs", {
  expect_error({
    set.seed(1234)
    norm::rngseed(123126)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z)
    imps <- catImp(temp, M=10, pd=TRUE, type=1)
  }, NA)
})

test_that("Three-way association imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    norm::rngseed(123126)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    w <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z,w)
    imps <- catImp(temp, M=10, pd=FALSE, type=2)
  }, NA)
})


test_that("Three-way association imputation PD draw runs", {
  expect_error({
    set.seed(1234)
    norm::rngseed(123126)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    z <- 1+(runif(n)<0.5)
    w <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y,z,w)
    imps <- catImp(temp, M=10, pd=TRUE, type=2)
  }, NA)
})
