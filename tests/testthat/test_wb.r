context("Testing within between function")

test_that("Runs a simple example without errors", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normUniImp(temp, y~x, M=10, pd=FALSE)

    analysisFun <- function(inputData) {
      mod <- lm(y~x, data=inputData)
      list(est=coef(mod), var=vcov(mod))
    }
    withinBetween(imps,analysisFun)
  }, NA)
})

test_that("If the imps don't have a pd attribute and user doesn't specify pd argument, you get an error", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normUniImp(temp, y~x, M=10, pd=FALSE)
    attributes(imps) <- NULL

    analysisFun <- function(inputData) {
      mod <- lm(y~x, data=inputData)
      list(est=coef(mod), var=vcov(mod))
    }
    withinBetween(imps,analysisFun)
  })
})

