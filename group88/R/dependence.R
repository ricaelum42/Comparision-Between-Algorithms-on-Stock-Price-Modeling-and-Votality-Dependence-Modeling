#' Compute returns from given data
#'
#' @param S A nxd matrix, vector or \eqn{xts} object input into function to calculate returns from.
#' @param type Type of calculated return, one of "absolute", "relative" and "log".
#' @return A matrix, vector or \eqn{xts} object of returns.
#' @details Absolute returns: \eqn{y){ti} = s_{ti} - s_{t-1,i}}
#' Relative returns: \eqn{y_s_{ti} - (s_{ti} - s_{t-1,i})/s_{t-1,i}}
#' Log returns: \eqn{y_{ti} = log(s_{ti}) - log(s_{t-1,i})}
#' @export
getreturns <- function(S, type) {
  if (missing(type)) {
    type = "absolute"
  }

  stopifnot(type %in% c("absolute", "relative", "log"))

  if (type == "absolute") {
    diff(S)[-1]
  } else if(type == "relative") {
    diff(S)[-1] / (S[-nrow(S)])
  } else if(type == "log") {
    diff(log(S))[-1]
  }
}


#' Generate random Gaussian Copula variates
#'
#' @param n Number of observations generated
#' @param sigma Sigma value for the multivariate gaussian model
#' @return Generated random gaussian-copula variate
rgCopula <- function(n, sigma) {
  pnorm(rmvnorm(n, sigma = sigma))
}


#' Conputes distribution function of Gaussian Copula
#' @param U Input Variables. Distribution is U(0,1)
#' @param sigma Covariance of gaussian copula
#' @return Distribution function of input
#' @export
pgCopula <- function(U, sigma) {
  mg.c <- qnorm(U)
  apply(mg.c, 1, function(z) {pmvnorm(upper = z, sigma = sigma)})
}

#' Compute dendity of Gaussian Copula
#' @param Input Variable. Distribution is U(0,1)
#' @param sigma Covariate of gaussian copula
#' @param log Logical variable. If true, return log density
#' @export
dgCopula <- function(U, sigma, log = FALSE) {
  Z <- qnorm(U)

  d.log <- dmvnorm(Z, sigma=sigma, log=TRUE) - rowSums(dnorm(X, log=TRUE))
  if (log) {
    d.log
  }else {
      exp(d.log)
  }
}


#' Gaussian-Copula Log-Likelihood
#'
#' @param sigma Upper triangular of sigma for multi-variate gaussian model
#' @param U nxd matrix of observations with distribution U~(0,1)
#' @return Log-likelihood of Gaussian Copula
#' @export
gCopula.LogLik <- function(sigma, U) {

  sum(dgCopula(U, sigma, log = TRUE))

}



#' Generate random t-Copula variates
#'
#' @param n Number of observations generated
#' @param sigma Sigma value for the multivariate t model
#' @param Degrees of freedom for the t-copula
#' @return Generated random t-copula variate
rtCopula <- function(n, sigma, df) {
  pt(rmvt(n, sigma = sigma, df = df), df = df)
}

#' Conputes distribution function of t Copula
#' @param U Input Variables. Distribution is U(0,1)
#' @param sigma Covariance of t copula
#' @param df Degrees of freedom of t copula
#' @return Distribution function of input
#' @export
ptCopula <- function(U, sigma, df) {
  mg.c <- qt(U, df=df)

  apply(mg.c, 1, function(z) {pmvt(upper = z, df=df, sigma = sigma)})
}

#' Compute dendity of t Copula
#' @param Input Variable. Distribution is U(0,1)
#' @param sigma Covariate of t copula
#' @param df Degrees of freedom of t copula
#' @param log Logical variable. If true, return log density
#' @export
dtCopula <- function(u, sigma, df, log = FALSE) {
  Z <- qt(U, df = df)

  d.log <- dmvt(X, sigma = sigma, df= df, log = TRUE) - rowSums(dt(X, df = df, log=TRUE))

  if (log) {
    d.log
  }else {
    exp(d.log)
  }
}


#' t-Copula Log-Likelihood
#'
#' @param sigma Upper triangular of sigma for multi-variate t model
#' @param df Degrees of freedom for multi-variate t model
#' @param U nxd matrix of observations with distribution U~(0,1)
#' @return Log-likelihood of t-Copula
#' @export
tCopula.Loglik <- function(sigma, df, U) {
  sum(dtCopula(U, sigma, df, log = TRUE))
}


#' Fit a gaussian or t copula and return copula parameters using maximum likelihood method
#'
#' @param U nxd matrix of observations with distribution U~(0,1)
#' @param type Type of copula to fit. One of "gaussian" and "t"
#' @return Returns fitted copula variables. Sigma for both gaussian and t copula and df for t copula only
#' @details Using \eqn{optim} to optimize loglikelihood function. The method used in optimize is "BFGS".
#' @export
fitCopula.GT <- function(U, type) {
  D = ncol(U)
  param <- cor(U)
  param <- param[upper.tri(param)]

  if (type == "t") {param <- c(param, 4)}

  gc.fit <- optim(param, fn = function(para) {
    if (type == "gaussian") {
      sig <- diag(0, nrow=D)
      sig[upper.tri(sig)] <- para
      sig <- sig + t(sig)
      diag(sig) <- rep.int(1, D)

      -gCopula.LogLik(sig, U)
      }else if(type == "t") {
      sig <- diag(0, nrow = D)
      sig[upper.tri(sig)] <- para[-length(para)]
      sig <- sig+ t(sig)
      diag(sig) <- rep.int(1, D)

      df <- para[length(para)]

      -tCopula.Loglik(sig, df, U)
    }
  }, method = "BFGS")
  print(sprintf("Convergence:%d, LogLik:%.2f, Counts:%d", gc.fit$convergence, gc.fit$value, gc.fit$counts[1]))

  copula.mle <- gc.fit$par

  if (type == "gaussian") {
    sig.mle <- copula.mle
    }else if (type == "t") {
    sig.mle <- copula.mle[-length(copula.mle)]
    df.mle <- copula.mle[length(copula.mle)]
  }

  P <- diag(0, nrow=D)
  P[upper.tri(P)] <- sig.mle
  P <- P + t(P)
  diag(P) <- rep.int(1,D)

  if (type == "gaussian") {
    list(Sigma = P)
  }else if (type == "t") {
      list(Sigma = P, df = df.mle)
    }
}

#' GARCH loglikelihood function.
#'
#' @param omega,alpha,beta GARCH parameters (scalars; see \code{\link{garch.sim}}).
#' @param eps Vector of GARCH observations \code{eps1, ..., epsN}.
#' @return The GARCH log-likelihood.
#' @details Validates parameters by checking that \code{all(sig2 > 0) == TRUE}.
#' @export
garch.loglik <- function(omega, alpha, beta, df, eps, sig2) {
  N <- length(eps)
  if(any(sig2 < 0)) return(-Inf) # prevent parameters giving negative variances
  N * log(gamma((df+1)/2) / (sqrt(pi*(df-2))*gamma(df/2))) -
    1/2 * sum(log(sig2)) - (df+1)/2*sum(log(1+ (eps^2) / (sig2 * (df - 2))))
}


#' Log-likelihood function of Copula-Garch(1,1), with Gaussian Copula or t Copula.
#'
#' @param alpha A vector contains alpha values for Garch model for each stock.
#' @param beta A vector contains beta valus for Garch model for each stock.
#' @param omega A vector contains omega values for Garch model for each stock.
#' @param mu A vector contains mu values for each Garch model for each stock.
#' @param shape A vector contains the degrees of freedom for Garch residual distribution for each stock.
#' @param rho A vector contains the upper triangular values of correlation matrix for the copula.
#' @param df A vector contains the degrees of freedom for the copula. Only needed for t-Copula.
#' @param type Type of copula. "gaussian" or "t".
#' @param Y A matrix, vector or xts object of returns.
#' @return Loglikelihood of Copula-Garch
#' @export
Garchcop.LogLik <- function(alpha, beta, omega, mu, shape, rho, df, type, Y) {
  # Only works for gaussian and t copula
  stopifnot(type %in% c("gaussian", "t"))
  # Parameter Checking
  stopifnot(!missing(df) | (type != "t"))

  D <- ncol(Y)
  T <- nrow(Y)

  # fit parameters to a garch model
  garch.fits <- lapply(1:D, function(i){
    spec <- ugarchspec(variance.model =list(model = "sGARCH", garchOrder=c(1,1)),
                       mean.model = list(armaOrder=c(0,0)),
                       distribution.model = "std")
    setfixed(spec) <- list(alpha1 = alpha[i], beta1 = beta[i], omega = omega[i], mu = mu[i], shape = shape[i])
    ugarchfilter(spec, Y[ ,i])
  })

  # get sigma for each model
  sig <- lapply(garch.fits, sigma)

  sig <- matrix(unlist(sig), ncol = D)
  colnames(sig) <- colnames(Y)

  # get residuals for each model
  X <- lapply(garch.fits, function(fit) {
    residuals(fit, standardize = TRUE)
  })

  X <- matrix(unlist(X), ncol=D)
  colnames(X) <- colnames(Y)

  # transform residuals to norm scale
  Z <- sapply(1:D, function(i){
    U <- pt(X[,i], df = shape[i])
    if (type == "gaussian") {
      qnorm(U)
    } else if (type == "t") {
      qt(U, df = df)
    }
  })

#  Z <- matrix(unlist(Z), ncol = D)
  colnames(Z) <- colnames(Y)

  L <- sapply(1:D, function(i){
    if (type == "gaussian") {
      log(dt(X[, i], df = shape[i])) - log(sig[,i]) - log(dnorm(Z[,i]))
    } else if (type == "t") {
      log(dt(X[, i], df = shape[i])) - log(sig[,i]) - log(dt(Z[,i], df = df))
    }
  })

#  L <- sum(matrix(unlist(L), ncol = D))

  if (type == "gaussian") {
    L <- sum(dmvnorm(Z, sigma = rho)) + sum(L)
  }else if (type == "t") {
      L <- sum(dmvt(Z, sigma = rho, df = df)) + sum(L)}

  L
}




#' Fit a Corpula-Garch(1,1) using maximum likelihood method
#'
#' @param Y a univariate data object.
#' @param type Type of Copula-Garch model, one of "gaussian" and "t"
#' @return A list of Copula-Garch(1,1) parameters
#' @details The garch model used in the Copula-Garch is \eqn{"sGarch"} with a standard t residual distribution.
#' The maximum-likeihood is optimized using \eqn{optim} function, with parameters set to default.
#' @export
fitgarchCopula <- function(Y, type){
  garch_basic <- ugarchspec(variance.model =list(model = "sGARCH", garchOrder=c(1,1)),
                            mean.model = list(armaOrder=c(0,0)),
                            distribution.model = "std")
  fit.garch <- lapply(Y, function(x) ugarchfit(garch_basic, x))
  X <- as.matrix(do.call(merge, lapply(fit.garch, residuals, standardize = TRUE)))
  T <- nrow(X)
  D <- ncol(X)

  coef <- lapply(1:D, function(i) {
    fit.garch[[i]]@fit$coef
  })
  coef <- matrix(unlist(coef), nrow=5)
  rownames(coef) <- c("mu", "omega", "alpha1", "beta1", "shape")
  mu <- coef["mu",]
  omega <- coef["omega",]
  alpha <- coef["alpha1",]
  beta <- coef["beta1",]
  shape <- coef["shape", ]
  rho <- cor(Y)
  rho <- rho[upper.tri(rho)]
  df <- 4

  para <- c(alpha, beta, omega, mu, shape, rho, df)

  ggc.fit <- optim(par = para,
                   fn = function(ltheta) {
                     alphas <- ltheta[1:D]
                     betas <- ltheta[(D+1) : (D*2)]
                     omegas <- ltheta[(D*2+1) : (D*3)]
                     mus <- ltheta[(D*3+1): (D*4)]
                     shapes <- ltheta[(D*4+1) : (D*5)]
                     rhos <- diag(0, nrow = D)
                     rhos[upper.tri(rhos)] <- ltheta[(D*5+1) : (length(ltheta)-1)]
                     rhos <- rhos + t(rhos)
                     diag(rhos) <- rep(1, D)
                     df <- ltheta[length(ltheta)]

                     -Garchcop.LogLik(alphas, betas, omegas, mus, shapes, rhos, df, type, Y)
                   })

  print(sprintf("Convergence:%d, LogLik:%.2f, Counts:%d", ggc.fit$convergence, ggc.fit$value, ggc.fit$counts[1]))

  ggc.mle <- ggc.fit$par


  alpha.mle <- ggc.mle[1:D]
  beta.mle <- ggc.mle[(D+1) : (D*2)]
  omega.mle <- ggc.mle[(D*2+1) : (D*3)]
  mu.mle <- ggc.mle[(D*3+1): (D*4)]
  shape.mle <- ggc.mle[(D*4+1) : (D*5)]
  rho.mle <- diag(0, nrow=D)
  rho.mle[upper.tri(rho.mle)] <- ggc.mle[(D*5+1) : (length(ggc.mle) - 1)]
  rho.mle <- rho.mle + t(rho.mle)
  diag(rho.mle) <- rep(1, D)
  df.mle <- ggc.mle[length(ggc.mle)]

  if (type == "gaussian") {
    list(alpha = alpha.mle, beta=beta.mle, omega=omega.mle, mu=mu.mle, shape = shape.mle, rho = rho.mle)
  } else if (type == "t") {
    list(alpha = alpha.mle, beta=beta.mle, omega=omega.mle, mu=mu.mle, shape = shape.mle, rho = rho.mle, df = df.mle)
  }
}


#' Simulate from copula-garch models given start values or continue from given data
#'
#'  @param n Number of simulated data points
#'  @param type type of copula. One of "gaussian" or "t"
#'  @param param List of copula-garch model parameters.
#'  @param sig0 Vector of start values of sigma for each stock. If not provided, will continue from Y.
#'  @param eps0 Vector of start values of residual for each stock. If not provided, will continue from Y.
#'  @param returns0 Vector of start values of returns for each stock. If not provided, will continue from Y.
#'  @param Y Data to continue with. Either Y or all of sig0, eps0 and returns0 need to provided
#'  @return Simulated copula-garch data.
#'  @export
copulaGarch.sim <- function(n, type, param, sig0, eps0, returns0, Y){
  D <- length(param$alpha)


  stopifnot(type %in% c("gaussian", "t"))
  stopifnot(!missing(Y) | (!missing(sig0) & !missing(eps0) & !missing(returns0)))

  if (type == "gaussian") {
    U <- rgCopula(n, param$rho)
  }else if (type == "t") {
      U <- rtCopula(n, param$rho, param$df)
  }
#  U <- matrix(runif(D * n), ncol = D)

  Z <- sapply(1:D, function(i) {
    qt(U[,i], df = param$shape[i]) / (sqrt(param$shape[i] / (param$shape[i] - 2)))
  })

  colnames(Z) <- colnames(Y)

  if (!missing(Y)) {
    T <- nrow(Y)
    garch.fits <- lapply(1:D, function(i) {
      spec <- ugarchspec(variance.model =list(model = "sGARCH", garchOrder=c(1,1)),
                         mean.model = list(armaOrder=c(0,0)),
                         distribution.model = "std")
      setfixed(spec) <- list(alpha1 = param$alpha[i], beta1 = param$beta[i], omega = param$omega[i],
                             mu = param$mu[i], shape = param$shape[i])
      filter <- ugarchfilter(spec, Y[ ,i])
      sigmas <- sigma(filter)

      epss <- residuals(filter)

      list(spec = spec, sigma0 = sigmas[T], eps0 = epss[T])
  })}

  if (missing(sig0)) { sig0 <-  sapply(1:D, function(i) {garch.fits[[i]]$sigma0}) }
  if (missing(eps0)) { eps0 <- sapply(1:D, function(i) {garch.fits[[i]]$eps0}) }
  if (missing(returns0)) { returns0 <- sapply(1:D, function(i) {Y[T, i]}) }


  garch.sim <- lapply(1:D, function(i) {
    ugarchpath(garch.fits[[i]]$spec, n.sim = n, presigma = sig0[i],
               prereturns = returns0[i], preresiduals = eps0[i],
               custom.dist = list(name = "sample", distfit = Z[ ,i ,drop = FALSE]))
  })

  Y.sim <- sapply(garch.sim, function(x) fitted(x))
  colnames(Y.sim) <- colnames(Y)

  Y.sim
}


#' A Garch GC simulate function that takes fitted garch models
#' @param n Number of simulated data points
#' @param model List of fitted garch models
#' @param List contains the sectir-wised covariance matrix
#' @param D number of stocks in total
#' @return Simulated data
#' @export
Garchgc.sim <- function(n, model, var, D) {
  df.lst <- unlist(lapply(1:D, function(i) {fit.garch[[i]]@fit$coef["shape"]}))

  U <- lapply(var, function(v) {
    v <- (v + t(v)) / 2
    rgCopula(m, v)
  })

  U <- matrix(unlist(U), ncol = D)
  Z <- sapply(1:D, function(i) {
    qt(U[,i], df = df.lst[i]) / (sqrt(df.lst[i] / (df.lst[i] - 2)))
  })

  sim <- lapply(1:D, function(i) {
    ugarchsim(model[[i]], n.sim = n, custom.dist = list(name = "sample", distfit = Z[,i, drop = FALSE]))
  })

  sapply(sim, function(x) fitted(x))
}



#--- Functions for the hierarchical multivariate normal model ------------------

# Model is:
# S = (S1, ..., S_K) ~ N(0, Sigma)
# V_ik | S ~ind N(S_k, Psi_k)

#' Intergroup and intragroup variance matrix.
#'
#' @param Sigma \code{K x K} intergroup variance matrix.
#' @param Psi List of \code{K} intragroup variance matrices.
#' @param group Length-\code{N} vector of integers between 1 and \code{K} of class labels.
#' @return Variance matrix of size \code{(N+K) x (N+K)}.
#' @details Returns the variance of \code{(S, V)}, where
#' \preformatted{
#'   S = (S_1, ..., S_K) ~ N(0, Sigma)
#'   V_k = (V(k)_1, ..., V(K)_NK) | S ~ind N(S_k, Psi_k),
#' }
#' and \code{N = N1 + ... + NK}.
hier.var <- function(Sigma, Psi, group) {
  K <- ncol(Sigma)
  if(!all(sort(unique(group)) == 1:K)) {
    stop("group must be integer vector with entries in 1:ncol(Sigma), and at least one of each.")
  }
  N <- length(group)
  Vfull <- matrix(NA, N+K, N+K) # full variance matrix
  Vfull[N + 1:K,N + 1:K] <- Sigma # var(S)
  # cov(S, V) and cov(V, S)
  ind <- as.matrix(expand.grid(1:K, group))
  Vfull[N + 1:K,1:N] <- Sigma[ind]
  Vfull[1:N,N + 1:K] <- t(Vfull[N + 1:K,1:N])
  # var(V)
  ind <- as.matrix(expand.grid(group, group))
  Vfull[1:N,1:N] <- Sigma[ind]
  ind <- rep(FALSE, N+K)
  for(ii in 1:K) {
    ind[1:N] <- group == ii # which(group == ii) resizes ind each time
    Vfull[ind,ind] <- Vfull[ind,ind] + Psi[[ii]]
  }
  Vfull
}

#' Compute the Log Likelihood of hierarchical multivariate normal distribution
#' @param S Matrix denote the mean of sector K for each day
#' @param Z Matrix denote the normalized residual of of each stock for each day
#' @param Sigma \code{K x K} intergroup variance matrix.
#' @param Psi List of \code{K} intragroup variance matrices.
#' @param group Length-\code{N} vector of integers between 1 and \code{K} of class labels.
#' @return Loglikelihood of sector-level hierarchical modelling
#' @export
hier.logLik <- function(S, Z, Sigma, Psi, group) {
  S.mu <- S$mu
  S.sigma <- S$sigma
  T <- nrow(S.mu)
  K <- length(unique(group))
  L1 <- -1/2 * (T * log(det(Sigma)) +
                sum(sapply(1:T, function(t) sum(diag(solve(Sigma) %*% (S.sigma + S.mu[t,] %*% t(S.mu[t,])))))))

  L2 <- sapply(1:K, function(k) {
    Psik <- Psi[[k]]
    groupk <- group == k
    Zk <- Z[,groupk]
    Sk <- S.mu[,k]
    Tk <- Reduce("+", lapply(1:T, function(t) Zk[t,] %*% t(Zk[t,])))
    Jn <- rep(1, ncol(Zk))
    v1 <- T * log(det(Psik))
    v2 <- sum(diag(solve(Psik) %*% Tk))

    Tk2 <- Reduce("+", lapply(1:T, function(t) Sk[t]*(Zk[t, ] %*% t(Jn) + Jn %*% t(Zk[t, ]))))
    v3 <- -sum(diag(solve(Psik) %*% Tk2))



    Tk3 <- Reduce("+", lapply(1:T, function(t) {
      SS <- S.sigma + S.mu[t,] %*% t(S.mu[t,])
      M3 <- SS[k,k] * Jn
      M3 %*% t(Jn)
    }))
    v4 <- sum(diag(solve(Psik) %*% Tk3))




    v1 + v2 + v3 + v4
  })
  L2 <- -1/2 * sum(L2)
  L1 + L2
}

#' Computes the conditional expectation of stock price given normalized residual
#' @param Sigma \code{K x K} intergroup variance matrix.
#' @param Psi List of \code{K} intragroup variance matrices.
#' @param Z Matrix denote the normalized residual of of each stock for each day
#' @param group Length-\code{N} vector of integers between 1 and \code{K} of class labels.
#' @export
S.expect <- function(Sigma, Psi, Z, group) {
  K <- ncol(Sigma)
  N <- length(group)
  Var <- hier.var(Sigma, Psi, group)
  mu.hat <- Var[N+1:K, 1:N] %*% solve(Var[1:N, 1:N]) %*% t(Z)
  Sigma.hat <- Var[N+1:K, N+1:K] - Var[N+1:K, 1:N] %*% solve(Var[1:N, 1:N]) %*% Var[1:N, N+1:K]
  list(mu = t(mu.hat), sigma = Sigma.hat)
}

#' Use in the maximization step of hierarchical modeling for updating the value of Sigma
#' @param S list denotes the distribution parameter of S
#' @return A
hier.A <- function(S) {
  mu.hat <- S$mu
  T <- nrow(mu.hat)
  sigma.hat <- S$sigma
  Reduce("+", lapply(1:T, function(t) sigma.hat + mu.hat[t,] %*% t(mu.hat[t,])))
}

#' Used in the maximization step of hierarchical modelling for updating the value of Psi
#' @param S List denotes the distribution of S
#' @param Z Matrix denote the normalized residual of of each stock for each day
#' @param group Length-\code{N} vector of integers between 1 and \code{K} of class labels.
#' @return B for each sector
hier.B <- function(S, Z, group) {
  S.mu <- S$mu
  S.sigma <- S$sigma
  K <- ncol(S.mu)
  T <- nrow(S.mu)
  lapply(1:K, function(k) {
    groupk <- group == k
    Zk <- Z[,groupk]
    Sk <- S.mu[,k]
    Tk <- Reduce("+", lapply(1:T, function(t) Zk[t,] %*% t(Zk[t,])))
    Jn <- rep(1, ncol(Zk))

    Tk2 <- Reduce("+", lapply(1:T, function(t) Sk[t]*(Zk[t, ] %*% t(Jn) + Jn %*% t(Zk[t, ]))))



    Tk3 <- Reduce("+", lapply(1:T, function(t) {
      SS <- S.sigma + S.mu[t,] %*% t(S.mu[t,])
      M3 <- SS[k,k] * Jn
      M3 %*% t(Jn)
    }))

    Tk - Tk2 + Tk3
  })
}
