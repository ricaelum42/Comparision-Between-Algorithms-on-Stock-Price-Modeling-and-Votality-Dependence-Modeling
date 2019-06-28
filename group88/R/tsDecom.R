#' Fit time series model via decomposition
#'
#' @import stats splines FNN tseries optimCheck
#'
#' @description a different process of time series prediction: decompose the time series into three pieces,
#' seasonality, trend and random residuals; build a model to each of them and forecast; unit all three pieces.
#'
#' @param s a stock dataset, should be a vector
#' @param start the time of the first observation. See \link{ts}
#' @param end the time of the last observation. See \link{ts}
#' @param frequency the number of observations per unit of time. See \link{ts}
#' @param alpha critical value of \link{runs.test}
#' @param trendSmooth smoothing method to trend. Currently the method should be one of these three \code{natural spline}, \code{KNN} and \code{lm}
#' @param seed the random seed
#' @param resHandle method to residuals. Currently the method should be one of these three \code{AR(1)}, \code{AR(1)} with \code{t} distribution and \code{Regime switching model AR(1)}
#' @param bound Give the bound of the residuals generation. If \code{NULL}, the residuals would be generated without restriction.
#' @param optimCheck If \code{TRUE}, the optim value would be checked. See \link{optimCheck}
#' @param HardMargin If \code{TRUE}, the prediction of residuals would be ignored.
#'
#' @return An object of the class \code{tsDecom} with following components:
#' \describe{
#'   \item{stress}{measurement of fit. The smaller the better}
#'   \item{pred}{prediction of the model. pred = pred_trend + pred_season + pred_res}
#'   \item{pred_season}{prediction of season}
#'   \item{pred_trend}{prediction of trend}
#'   \item{pred_res}{prediction of residuals}
#'   \item{s}{the original time series}
#'   \item{times}{the time of the time series}
#' }
#' @examples
#' tsD <- tsDecom(co2, start = c(1959, 1), frequency = 12)
#' plot(tsD$s, col = "red", type = "l", ylim = extendrange(c(tsD$s, tsD$pred)))
#' Xorder <- order(tsD$times)
#' lines(tsD$times[Xorder], tsD$pred[Xorder])
#'
#' @export


tsDecom <- function(s, start  = NULL, end = NULL, frequency = 365, alpha = 0.05,
                    trendSmooth = c("Natural", "KNN", "lm"), seed = 2018440,
                    resHandle = c("ARGaussian", "ARt", "RS"), bound = c(-3, 3),
                    optimCheck = FALSE, HardMargin = FALSE){
  set.seed(seed)
  s <- na.omit(s)
  if(is.null(start) && is.null(end) ){stop("start or end time should not be NULL at the same time")}
  else if(is.null(start)){
    s <- ts(s, end = end, frequency = frequency)
  }else if(is.null(end) ){
    s <- ts(s, start = start, frequency = frequency)
  }
  DecomSTL <- stl(s, "periodic")
  season <- DecomSTL$time.series[,1]
  trend <- DecomSTL$time.series[,2]
  res <- DecomSTL$time.series[,3]
  times <- time(s)
  # season
  fit_season <- loess(season ~ times, span = 0.01, control=loess.control(surface="direct"))
  pred_season <- predict(fit_season, data.frame(times))
  # trend
  trendSmooth <- match.arg(trendSmooth)
  if(trendSmooth == "Natural"){
    fitSmooth <- lm(trend ~ ns(times , df = 10))
    pred_trend <- predict(fitSmooth, newdata=data.frame(x = times))
  }else if(trendSmooth == "KNN"){
    fitKNN <- knn.reg(times , y= trend, k=5)
    pred_trend <- fitKNN$pred
  }else if(trendSmooth == "lm"){
    fitLinear <- lm(trend~times)
    pred_trend <- predict(fitLinear, newdata = data.frame(times))
  }
  # residuals
  run <- runs.test(factor(sign(res)))
  isDiff <- FALSE
  if(run$p.value < alpha){
    res <- diff(res)
    isDiff <- TRUE
  }
  # AR(1) model for res

  resHandle <- match.arg(resHandle)
  if (resHandle == "ARGaussian"){
    model <- arima(res, c(1,0,0), include.mean = F)
    pred_res <- predAR(Rt1 = 0, Rt2 = 0, phi = model$coef, len = length(times),
                       sd = sqrt(model$sigma2), bound = bound, isDiff)
  }else if(resHandle == "ARt"){
    optL <- optim(c(0.5, 0.5,3), Lfun,Y= res)
    pred_res <- predLogAR(Rt1 = 0, Rt2 = 0, phi = optL$par[2], len = length(times),
                          df = optL$par[3], sigma = optL$par[1], bound = bound, isDiff)
    if(optimCheck){
      optim_proj(optL$par,fun = function(par) Lfun(par = par, Y = res) )
    }

  }else if(resHandle == "RS"){
    optR <- optim(c(1, 1,1), Rfun, Y= res)
    pred_res <- predRSAR(Rt1 = 0, Rt2 = 0, phiE = optR$par[1], phiR = optR$par[2],
                         len = length(times), sigma = optR$par[3], bound = bound, isDiff)
    if(optimCheck){
      optim_proj(optR$par,fun = function(par) Rfun(par = par, Y = res) )
    }
  }
  if(HardMargin){
    pred <- pred_season + pred_trend
  }else{
    pred <- pred_season + pred_trend + pred_res
  }
  stress <- sum( (s - pred)^2)/ sum( (s - mean(s))^2 )
  list(stress = stress, pred = pred, pred_trend = pred_trend, pred_season = pred_season,
       pred_res = pred_res, s = s, times = times)
}




############ helper function

predAR <- function(Rt1 = 0, Rt2 = 0, phi = model$coef, len = length(times),
                   sd = sqrt(model$sigma2),  bound = NULL, isDiff){
  R <- rep(0, len-2)
  R <- c(Rt1, Rt2, R)
  if(is.null(bound)){
    for(i in 3 : len){
      if(isDiff){
        R[i] <- R[i-1] + phi*(R[i - 1] - R[i-2]) + rnorm(1, 0, sd)
      }else{
        R[i] <- phi*(R[i - 1]) + rnorm(1, 0, sd)
      }
    }
  }else{
    for(i in 3 : len){
      if(isDiff){
        R[i] <- R[i-1] + phi*(R[i - 1] - R[i-2]) + rnorm(1, 0, sd)
      }else{
        R[i] <- phi*(R[i - 1]) + rnorm(1, 0, sd)
      }
      if(R[i] > bound[2] || R[i] < bound[1]){
        while(R[i] > bound[2] || R[i] < bound[1]){
          if(isDiff){
            R[i] <- R[i-1] + phi*(R[i - 1] - R[i-2]) + rnorm(1, 0, sd)
          }else{
            R[i] <- phi*(R[i - 1]) + rnorm(1, 0, sd)
          }
        }
      }
    }
  }
  R
}

Lfun <- function(par = par, Y = Y){
  sig <- par[1]
  phi <- par[2]
  v <- par[3]
  Yn <- Y[-1]
  Yn_1 <- Y[-length(Y)]
  n <- length(Y)
  L <- n*(log(gamma((v+1)/2)) - log(gamma(v/2)) - 0.5*log(v*pi)) - n*log(sig) - (v+1)/2*sum(log(1+(Yn - phi*Yn_1)^2/(v*sig^2))) - (v+1)/2*log(1+(1-phi^2)*Y[1]^2/(sig^2*v))
  L <- -L
  return(L)
}


predLogAR <- function(Rt1 = 0, Rt2 = 0, phi = optL$par[2], len = length(times),
                      df = optL$par[3], sigma = optL$par[1], bound = NULL, isDiff ){
  R <- rep(0, len-2)
  R <- c(Rt1, Rt2, R)
  if(is.null(bound)){
    for(i in 3 : len){
      if(isDiff){
        R[i] <- R[i-1] + phi*(R[i - 1] - R[i-2]) + sigma * rt(1, df)
      }else{
        R[i] <- phi*(R[i - 1]) + sigma * rt(1, df)
      }
    }
  }else{
    for(i in 3 : len){
      if(isDiff){
        R[i] <- R[i-1] + phi*(R[i - 1] - R[i-2]) + sigma * rt(1, df)
      }else{
        R[i] <- phi*(R[i - 1]) + sigma * rt(1, df)
      }
      if(R[i] > bound[2] || R[i] < bound[1]){
        while(R[i] > bound[2] || R[i] < bound[1]){
          if(isDiff){
            R[i] <- R[i-1] + phi*(R[i - 1] - R[i-2]) + sigma * rt(1, df)
          }else{
            R[i] <- phi*(R[i - 1]) + sigma * rt(1, df)
          }
        }
      }
    }
  }
  R
}

Rfun <- function(par = par, Y = Y){
  phiE <- par[1]
  phiR <- par[2]
  sig <- par[3]
  L <- 0
  for(i in 2:length(Y)){
    if(Y[i-1]>=0){
      delta <- 1
    }else{delta <- 0}
    L <- L- log(sig) -(Y[i] - (delta*phiE*Y[i-1] + (1 - delta)*phiR*Y[i-1]) )^2/(2*sig^2)
  }
  L <- -L
  return(L)
}

predRSAR <- function(Rt1 = 0, Rt2 = 0, phiE = optR$par[1], phiR = optR$par[2],
                     len = length(times), sigma = optR$par[3], bound = NULL, isDiff ){
  R <- rep(0, len-2)
  R <- c(Rt1, Rt2, R)
  if(is.null(bound)){
    for(i in 3 : len){
      if(isDiff){
        if(R[i-1] - R[i-2] >= 0 ){delta <- 1}
        else{delta <- 0}
        R[i] <- R[i-1] + delta*phiE * (R[i - 1] - R[i-2]) +
          (1-delta) * phiR * (R[i - 1] - R[i-2]) + rnorm(1, 0, sigma)
      }else{
        if(R[i-1] >= 0 ){delta <- 1}
        else{delta <- 0}
        R[i] <-  delta*phiE * (R[i - 1] ) +
          (1-delta) * phiR * (R[i - 1] ) + rnorm(1, 0, sigma)
      }
    }
  }else{
    for(i in 3 : len){
      if(isDiff){
        if(R[i-1] - R[i-2] >= 0 ){delta <- 1}
        else{delta <- 0}
        R[i] <- R[i-1] + delta*phiE * (R[i - 1] - R[i-2]) +
          (1-delta) * phiR * (R[i - 1] - R[i-2]) + rnorm(1, 0, sigma)
      }else{
        if(R[i-1] >= 0 ){delta <- 1}
        else{delta <- 0}
        R[i] <-  delta*phiE * (R[i - 1] ) +
          (1-delta) * phiR * (R[i - 1] ) + rnorm(1, 0, sigma)
      }
      if(R[i] > bound[2] || R[i] < bound[1]){
        while(R[i] > bound[2] || R[i] < bound[1]){
          if(isDiff){
            if(R[i-1] - R[i-2] >= 0 ){delta <- 1}
            else{delta <- 0}
            R[i] <- R[i-1] + delta*phiE * (R[i - 1] - R[i-2]) +
              (1-delta) * phiR * (R[i - 1] - R[i-2]) + rnorm(1, 0, sigma)
          }else{
            if(R[i-1] >= 0 ){delta <- 1}
            else{delta <- 0}
            R[i] <-  delta*phiE * (R[i - 1] ) +
              (1-delta) * phiR * (R[i - 1] ) + rnorm(1, 0, sigma)
          }
        }
      }
    }
  }
  R
}
