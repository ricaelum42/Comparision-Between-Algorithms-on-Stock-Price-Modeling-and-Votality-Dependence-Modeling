source("group88.R")
library(xts)
library(optimCheck)
load("stockprice.Rda")

time <- c("2012-01-03", "2017-12-31")
S <- stockprice[paste0(time, collapse = "/"),] # data
S <- na.omit(S)

Y <- getreturns(S, "log")
T <- nrow(Y)
D <- ncol(Y)

garch_basic <- ugarchspec(variance.model =list(model = "sGARCH", garchOrder=c(1,1)),
                          mean.model = list(armaOrder=c(0,0)),
                          distribution.model = "std")

fit.garch <- lapply(Y, function(x) ugarchfit(garch_basic, x))

eps <- as.matrix(do.call(merge, lapply(fit.garch, residuals)))
sigma2 <- as.matrix(do.call(merge, lapply(fit.garch, sigma)))^2


# MLE checck rugarch package-----------------------------------------------------------
objfun <- function(x) {
  garch.loglik(x[1], x[2], x[3], x[4], eps, sigma2)
}

for(i in 1:D) {
  Xfit <- fit.garch[[i]]@fit$coef[c("omega", "alpha1", "beta1", "shape")]
  oproj <- optim_proj(fun = objfun,              # objective function
                      xsol = Xfit,           # potential solution
                      maximize = FALSE,          # indicates that a local minimum is sought
                      xrng = 0.5)                 # range of projection plot: x_i +/- .5*|x_i|
  
}


# test residual distribution------------------------------------------------------------
X <- as.matrix(do.call(merge, lapply(fit.garch, residuals, standardize = TRUE)))
par(mfrow=c(4, 5))

for(i in 1:D) {
  hist(X[, i], probability = TRUE, main = sprintf("Residual Dist of stock %d", i))
  lines(seq(-5, 5, length = 100), dnorm(seq(-5, 5, length = 100)))
}


U <- sapply(1:D, function(i) U <- pt(X[,i], df = fit.garch[[i]]@fit$coef["shape"]))

for(i in 1:D) {
  hist(U[, i], probability = TRUE, main = sprintf("Residual Dist of stock %d", i))
  lines(seq(-1, 1, length = 100), dunif(seq(-1, 1, length = 100)))
}



Z <- qnorm(U)

for(i in 1:D) {
  hist(Z[, i], probability = TRUE, main = sprintf("Trans Residual Dist %d", i))
  lines(seq(-5, 5, length = 100), dnorm(seq(-5, 5, length = 100)))
}

par(mfrow=c(1, 1))
