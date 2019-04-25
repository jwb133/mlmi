context("Testing categorical imputation functions")
library(mlmi)

test_that("Saturated imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    rngseed(123126)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- catSatImp(temp, M=10, pd=FALSE)
  }, NA)
})

test_that("Saturated imputation with PD draw runs", {
  expect_error({
    set.seed(1234)
    rngseed(123126)
    n <- 100
    x <- 1+(runif(n)<0.5)
    y <- 1+(runif(n)<0.5)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- catSatImp(temp, M=10, pd=TRUE)
  }, NA)
})
