source("group88.R")
library(xts)
load("pricecopy.Rda")

#-----------------------------sector level hierarchical modelling ----------------------
pricecopy <- as.xts(pricecopy)
time <- c("2012-01-03", "2017-12-31")
S <- pricecopy[paste0(time, collapse = "/"),] # data
S <- na.omit(S)

Y <- getreturns(S, "log")
group <- c(rep(1,2), rep(2, 17), rep(3, 9), rep(4, 23), rep(5, 2), rep(6, 3), rep(7, 4), rep(8, 3))

T <- nrow(Y)
D <- ncol(Y)


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


#Test mdodel performance---------------------------------------------------------------------------------------------
B <- 200
m <- 40

stock.name <- colnames(S)

Y.sim <- lapply(1:B, function(b) {
  sim <- Garchgc.sim(m, fit.garch, Var.sector, D)
  colnames(sim) <- colnames(Y)
  sim
})

par(mfrow = c(2,2))

for(i in 1:D) {
  Ys <- Y[,i]
  Ys. <- sapply(Y.sim, function(y) y[,i])
  Ys.mean <- rowMeans(Ys.)
  Ys.CI <- apply(Ys., 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  tm <- index(pricecopy)
  start <- match(time[1], as.character(tm))
  past <- tm[start:(start+T-1)]
  future <- tm[(start+T):(start+T+m-1)]
  plot(past, Ys, type = "l", xlim = range(c(past, future)), xlab = "", ylab = "", main = stock.name[i]) 
  polygon(c(future, rev(future)), c(Ys.CI[1,], rev(Ys.CI[2,])),
        border = NA, col = "grey80") 
  lines(future, Ys.mean, col = "royalblue3")
  lines(future, Ys.CI[1,], col = "grey50") 
  lines(future, Ys.CI[2,], col = "grey50")
  legend("bottomright", bty = "n", lty = rep(1, 3),
       col = c("black", "royalblue3", "grey50"),
       legend = c("(Aggregated) loss", "(Simulated) predicted loss",
                  "95% CIs"))
}




tm <- index(pricecopy)
start <- match(time[1], as.character(tm))
past <- tm[start:(start+T-1)]
future <- tm[(start+T+1):(start+T+m+1)]
S.actual <- as.matrix(pricecopy[future, ])

S.start <- log(S)
S.start <- as.matrix(S.start[nrow(S), ])

S.sim <- lapply(Y.sim, function(y){
  log.sim <- log(S.actual[-41,]) + y
  exp(log.sim)
})

pt <- lapply(S.sim, function(s) {
  s <= S.actual[-1,]
})


pt <- Reduce("+", pt)/ 200
zt <- qnorm(pt)

par(mfrow = c(3,3))

stock.name <- colnames(S)

stock.names <- numeric(0)
p.val <- numeric(0)


for(i in 1:D) {
  zz <- zt[,i]
  zz <- zz[!is.infinite(zz)]
  if(length(zz) == 0) next
  if(is.na(zz)) next
  k <- shapiro.test(zz) 
  stock.names <- c(stock.names,  stock.name[i])
  p.val <- c(p.val, k$p.value)
  qqnorm(zz, main = stock.name[i])
}

library(knitr)
shapiro.result <- matrix(c(stock.names, p.val), ncol=2)
colnames(shapiro.result) <- c("stock", "p-val")
kable(shapiro.result)
