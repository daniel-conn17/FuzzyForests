#' Function to generate multi-class data from a multinomial logistic
#' regression.  Assumes there are 5 classes.  Only supports two modules for now.
#' Currently this function is used for testing.
#'
#' @title Multinomial Logistic Regression
#' @name multi_class_lr
#' @param n                 Sample size.
#' @param mod1_size         Size of first module.
#' @param mod2_size         Size of second module.
#' @param rho               Correlation of covariates.
#' @param beta              A matrix of parameters.
#' @return list with design matrix X, outcome y, and beta.
#' @note
#' This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
multi_class_lr <- function(n, mod1_size=10, mod2_size=10, rho=.8, beta=NULL){
  #Say there are 5 significant features.
  #We will assume a multinomial logistic regression.
  J <- 5 #J=number of outcomes
  #signalp refers to the dimensionality of module which
  #the significant covariates are a part of.
  signalp <- mod1_size
  if(is.null(beta)){
    beta <- diag(rep(5, 4))
    beta <- rbind(beta, matrix(0, signalp - 4, 4))
  }
  sigma <- matrix(rho, signalp, signalp)
  diag(sigma) <- 1
  signalX <- rmvnorm(n, mean=rep(0, signalp), sigma=sigma)

  y <- rep(NA, n)
  pimat <- matrix(NA, n, J)
  for(i in 1:n){
    lps <- signalX[i, ]%*%beta
    den <- 1 + sum(exp(lps))
    #first category is the reference category
    #pi1 = 1/(1 + sum(exp(x_{j}'b_{j}))
    pi <- c(1, exp(lps))/den
    pimat[i, ] <- pi
    outcm <- rmultinom(1, size=1, prob=pi)
    y[i] <- which(outcm == 1)
  }

  #additional noise covariates
  noisep <- mod2_size
  sigma <- matrix(rho, noisep, noisep)
  diag(sigma) <- 1
  noise <- rmvnorm(n, mean=rep(0, noisep), sigma=sigma)
  X <- cbind(signalX, noise)
  X <- as.data.frame(X)
  y <- as.factor(y)
  out <- list(X=X, y=y, beta=beta)
}

