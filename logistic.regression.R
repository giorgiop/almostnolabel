library(doMC)
registerDoMC()

#------------------------------------------------------------------------------
# Fully supervised logistic regression in Mean Map's fashion
#------------------------------------------------------------------------------

logistic.regression <- function(data, lambda) {
  
  f <- function(w) {
    xw    <- X %*% w
    ai    <- log(exp(xw) + exp(-xw))
    ai    <- ifelse(is.finite(ai), ai, xw) # Handle numerical overflow for log(exp(..))
    lterm <- sum(ai)
    rterm <- w %*% true.mean.op
    loss  <- as.numeric(lterm - rterm + (0.5*lambda*(w %*% w)))
    #cat(sprintf("loss: %f\n", loss))
    loss
  }
  
  g <- function(w) {
    xw    <- as.vector(X %*% w)
    lterm <- colSums(X*tanh(xw))
    lterm - true.mean.op + lambda*w
  }
  
  #number of sample
  M <- nrow(data)
  #number of features
  X <- as.matrix(data[,-(1:2)])
  K <- ncol(X)
  
  true.mean.op <- colSums(data[,1] * data[,-c(1,2)]) #Do not divide by M
  
  w0 <- rep(0.001,K)
  result <- optim(w0, fn=f, gr=g, method="L-BFGS-B")#, control=list(trace=3, maxit=100))
  result$par
}