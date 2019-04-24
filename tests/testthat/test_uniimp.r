context("Testing univariate imputation functions")
library(mlmi)

test_that("Univariate normal imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 10
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:5] <- NA
    temp <- data.frame(x,y)
    imps <- normUniImp(temp, y~x, M=10, pd=FALSE)
  }, NA)
})

test_that("Univariate normal imputation with PD runs", {
  expect_error({
    set.seed(1234)
    n <- 10
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:5] <- NA
    temp <- data.frame(x,y)
    imps <- normUniImp(temp, y~x, M=10, pd=TRUE)
  }, NA)
})

test_that("Univariate normal imputation with M=1 returns data frame not a list", {
  expect_equal({
    set.seed(1234)
    n <- 10
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:5] <- NA
    temp <- data.frame(x,y)
    imps <- normUniImp(temp, y~x, M=1, pd=TRUE)
    is.data.frame(imps)
  }, TRUE)
})

test_that("Within between function runs no PD, scalar model", {
  expect_error({
    set.seed(1234)
    n <- 10
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:5] <- NA
    temp <- data.frame(x,y)
    imps <- normUniImp(temp, y~x, M=10, pd=FALSE)

    myanalysis <- function(data) {
      n <- dim(data)[1]
      list(est=mean(data$y), var=var(data$y)/n)
    }

    result <- withinBetween(imps, myanalysis)
  }, NA)
})
