source("group88.R")
library(xts)
load("stockprice.Rda")


stocks <- c("BAM.A", "BMO", "BNS", "CM", "KFS", "MFC", "RY", "SLF", "TD") # electricity
stocks <- c("BB", "CAE", "GIBA")
time <- c("2012-01-03", "2017-12-31")
S <- stockprice[paste0(time, collapse = "/"),] # data
S <- na.omit(S)

S.elec <- stockprice[paste0(time, collapse = "/"), stocks] # data
S.elec <- na.omit(S.elec)


Y <- getreturns(S, "log")

param <- fitgarchCopula(Y, type="gaussian")


B <- 200
Y.sim <- lapply(1:B, function(x) {
  copulaGarch.sim(60, "gaussian", param, Y=Y)
})


Xs <- rowSums(Y)
Xs. <- sapply(Y.sim, rowSums)
Xs.mean <- rowMeans(Xs.)
Xs.CI <- apply(Xs., 1, function(x) quantile(x, probs = c(0.025, 0.975))) # CIs; (2, m)-matrix
alpha <- 0.99 # confidence level
VaR <- apply(Xs., 1, function(x) quantile(x, probs = alpha)) # VaR_alpha; m-vector

n <- nrow(Y)
tm <- index(stockprice)
start <- match(time[1], as.character(tm))
past <- tm[start:(start+n-1)]
future <- tm[(start+n):(start+n+60-1)]
plot(past, Xs, type = "l", xlim = range(c(past, future)), xlab = "", ylab = "") # actual (past) losses
polygon(c(future, rev(future)), c(Xs.CI[1,], rev(Xs.CI[2,])),
        border = NA, col = "grey80") # CI region
lines(future, Xs.mean, col = "royalblue3") # predicted aggregated loss
lines(future, Xs.CI[1,], col = "grey50") # lower CI
lines(future, Xs.CI[2,], col = "grey50") # upper CI
lines(future, VaR, col = "maroon3") # VaR_alpha
legend("bottomright", bty = "n", lty = rep(1, 4),
       col = c("black", "royalblue3", "grey50", "maroon3"),
       legend = c("(Aggregated) loss", "(Simulated) predicted loss",
                  "95% CIs", as.expression(substitute("Simulated"~VaR[a], list(a = alpha)))))


#-----------------------Sector Level hierarchical modelling-------------------------------------
group <- c(rep(1,2), rep(2, 17), rep(3, 9), rep(4, 23), rep(5, 2), rep(6, 3), rep(7, 4), rep(8, 3))


T <- nrow(S)
D <- ncol(S)


garch_basic <- ugarchspec(variance.model =list(model = "sGARCH", garchOrder=c(1,1)),
                          mean.model = list(armaOrder=c(0,0)),
                          distribution.model = "std")

fit.garch <- lapply(Y, function(x) ugarchfit(garch_basic, x))
X <- as.matrix(do.call(merge, lapply(fit.garch, residuals, standardize = TRUE)))
Z <- sapply(1:D, function(i) {
  U <- pt(X[,i], df = fit.garch[[i]]@fit$coef["shape"])
  qnorm(U)
})

S.sector <- sapply(1: length(unique(group)), function(i) {
  groupi <- group == i
  Si <- Z[, groupi]
  rowMeans(Si)
})

row.names(S.sector) <- as.character(row.names(X))


sigma0 <- var(S.sector)
Psi0 <- lapply(1:8, function(x) {
  Zk <- Z[,group == x]
  cor(Zk)
})

m <- 500
L.cur <- -Inf
L.prev <-0
while(abs(L.cur - L.prev) > 0.1) {
  L.prev <- L.cur
  
  S.cur <- S.expect(sigma0, Psi0, Z, group)
  A <- hier.A(S.cur)
  B <- hier.B(S.cur, Z, group)
  
  sigma0 <- A / T
  Psi0 <- lapply(B, function(x) x/T)
  
  L.cur <- hier.logLik(S.cur, Z, sigma0, Psi0, group)
  print(L.cur)
}

Vfull <- hier.var(sigma0, Psi0, group)

Var.sector <- lapply(1:length(unique(group)), function(i) {
  ind <- rep(FALSE, ncol(Vfull))
  ind[1:length(group)] <- group == i
  Vfull[ind, ind]
  
})


require(optimCheck)
objfun <- function(ltheta) {
    D <- ncol(Y)
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
  
    -Garchcop.LogLik(alphas, betas, omegas, mus, shapes, rhos, df, "gaussian", Y)
}

xfit <- c(param$alpha, param$beta, param$omega, param$mu, param$shape, param$rho[upper.tri(param$rho)], 4)




oproj <- optim_proj(fun = objfun,              # objective function
                    xsol = xfit,           # potential solution
                    maximize = FALSE,          # indicates that a local minimum is sought
                    xrng = .5)                 # range of projection plot: x_i +/- .5*|x_i|


garch.fits <- lapply(1:D, function(i){
  spec <- ugarchspec(variance.model =list(model = "sGARCH", garchOrder=c(1,1)),
                     mean.model = list(armaOrder=c(0,0)),
                     distribution.model = "std")
  setfixed(spec) <- list(alpha1 = alpha[i], beta1 = beta[i], omega = omega[i], mu = mu[i], shape = shape[i])
  ugarchfilter(spec, Y[ ,i])
})

X <- lapply(garch.fits, function(fit) {
  residuals(fit, standardize = TRUE)
})

X <- matrix(unlist(X), ncol=D)

U <- sapply(1:D, function(i){
  pt(X[,i], df = shape[i])
})


copu <- fitCopula.GT(U, "gaussian")

Xfit <- copu$Sigma[upper.tri(copu$Sigma)]

objfun <- function(x) {
  P <- diag(0, nrow=D)
  P[upper.tri(P)] <- x
  P <- P + t(P)
  diag(P) <- rep.int(1,D)
  -gCopula.LogLik(P, U)
}

oproj <- optim_proj(fun = objfun,              # objective function
                    xsol = Xfit,           # potential solution
                    maximize = FALSE,          # indicates that a local minimum is sought
                    xrng = .5)                 # range of projection plot: x_i +/- .5*|x_i|

