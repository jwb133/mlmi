context("Testing multivariate normal imputation functions")

test_that("MVN imputation no PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=10, pd=FALSE, rseed=4423)
  }, NA)
})

test_that("If you don't pass a seed, you get an error", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=10, pd=FALSE)
  })
})

test_that("Manually setting seed prior to calling normImp", {
  expect_error({
    set.seed(1234)
    norm::rngseed(51312)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=10, pd=FALSE, rseed=NULL)
  }, NA)
})

test_that("MVN imputation with PD draw runs", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=10, pd=TRUE, rseed=4423)
  }, NA)
})

test_that("MVN imputation returns a data frame when M=1", {test_that("MVN imputation returns a data frame when M=1", {
  expect_equal({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=1, pd=FALSE, rseed=4423)
    is.data.frame(imps)
  }, TRUE)
})
  expect_equal({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=1, pd=FALSE, rseed=4423)
    is.data.frame(imps)
  }, TRUE)
})

test_that("MVN imputation returns a list with correct pd value attribute", {
  expect_equal({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normImp(temp, M=10, pd=TRUE, rseed=4423)
    attr(imps, "pd")
  }, TRUE)
})
