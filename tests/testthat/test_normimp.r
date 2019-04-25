context("Testing multivariate normal imputation functions")
library(mlmi)

test_that("MVN imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    rngseed(123126)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=10, pd=FALSE)
  }, NA)
})

test_that("MVN imputation with PD draw runs", {
  expect_error({
    set.seed(1234)
    rngseed(123126)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=10, pd=TRUE)
  }, NA)
})

test_that("MVN imputation returns a data frame when M=1", {test_that("MVN imputation returns a data frame when M=1", {
  expect_equal({
    set.seed(1234)
    rngseed(123126)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=1, pd=FALSE)
    is.data.frame(imps)
  }, TRUE)
})
  expect_equal({
    set.seed(1234)
    rngseed(123126)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=1, pd=FALSE)
    is.data.frame(imps)
  }, TRUE)
})

test_that("MVN imputation returns a list with correct pd value attribute", {
  expect_equal({
    set.seed(1234)
    rngseed(123126)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=10, pd=TRUE)
    attr(imps, "pd")
  }, TRUE)
})
