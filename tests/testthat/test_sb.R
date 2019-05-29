context("Testing score based approach")

test_that("Runs a simple example without errors", {
  expect_error({
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    y <- x+rnorm(n)
    y[1:50] <- NA
    temp <- data.frame(x,y)
    imps <- normUniImp(temp, y~x, M=10, pd=FALSE)

    yonx <- function(inputData) {
      fitmod <- lm(y~x, data=inputData)
      list(est=c(fitmod$coef,sigma(fitmod)^2))
    }

    myScore <- function(inputData, parm) {
      beta0 <- parm[1]
      beta1 <- parm[2]
      sigmasq <- parm[3]
      res <- inputData$y - beta0 - beta1*inputData$x
      cbind(res/sigmasq, (res*inputData$x)/sigmasq, res^2/(2*sigmasq^2)-1/(2*sigmasq))
    }

    scoreBased(imps, analysisFun=yonx, scoreFun=myScore)
  }, NA)
})
