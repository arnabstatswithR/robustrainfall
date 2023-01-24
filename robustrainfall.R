# first install the package fitdistrplus

#--------------------

# mdpde.exp, mdpde.gamma, mdpde.lnorm, mdpde.weibull

# Inputs for these functions are data and a value of the tuning parameter. 
# They return the MDPDE estimates of the model parameters for exponential, gamma, lognormal and Weibull distributions respectively.

# WD.mdpde.exp, WD.mdpde.gamma, WD.mdpde.lnorm, WD.mdpde.weibull

# Inputs for these functions are data and a value of the tuning parameter. 
# They return the corresponding Wasserstein distances (WDs) for exponential, gamma, lognormal and Weibull distributions respectively.

# optim.alpha.exp, optim.alpha.gamma, optim.alpha.lnorm, optim.alpha.weibull

# Inputs for these functions are data only. 
# They return the optimal values of tuning parameter $\alpha$ by minimizing WDs for 
# exponential, gamma, lognormal and Weibull distributions respectively.

# sd.mdpde.exp, sd.mdpde.gamma, sd.mdpde.lnorm, sd.mdpde.weibull

# Inputs for these functions are data and a value of the tuning parameter. 
# They return the standard errors of the MDPDE estimates for exponential, gamma, lognormal and 
# Weibull distributions respectively. We compute the standard error using $B=1000$ bootstrap samples. 
# The user can set a different value of $B$ if required.

# RIC.exp, RIC.gamma, RIC.lnorm, RIC.weibull

# Inputs for these functions are data and a value of the tuning parameter. 
# They return the corresponding RIC for exponential, gamma, lognormal and Weibull distributions respectively.

#--------------------

mdpde.exp <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  
  lambda.mle <- 1 / mean(X)
  if(alpha == 0){
    return(lambda.mle)
  }else{
    log.lambda.mle <- log(lambda.mle)
    robust_est <- function(log.lambda){
      (exp(log.lambda)^alpha) * 
        (1 / (1 + alpha) - (1 + 1 / alpha) * mean(exp(-alpha * exp(log.lambda) * X))) + 1 / alpha}
    optim.out <- suppressWarnings(optim(par = log.lambda.mle, robust_est))
    return(exp(optim.out$par))}
}

#--------------------

mdpde.gamma <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  
  library(fitdistrplus)
  fit.gamma <- fitdist(X, distr = "gamma", method = "mle")
  shape_rate_hat <- fit.gamma$estimate
  if(alpha == 0){
    return(shape_rate_hat)
  }else{
    log.shape.rate <- log(shape_rate_hat)
    robust_est <- function(par){
      ((exp(par[2])^alpha) * gamma(1 + (exp(par[1]) - 1) * (alpha + 1)) / 
         gamma(exp(par[1]))^(1 + alpha) /
         (1 + alpha)^(1 + (exp(par[1]) - 1) * (alpha + 1)) - 
         (1 + 1 / alpha) * mean(dgamma(X, shape = exp(par[1]), 
                                       rate = exp(par[2]))^alpha)) + 1 / alpha}
    optim.out <- suppressWarnings(optim(par = log.shape.rate, robust_est))
    out <- optim.out$par
    names(out) <- c("shape", "rate")
    return(exp(out))}
}

#--------------------

mdpde.lnorm <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  
  library(fitdistrplus)
  fit.lnorm <- fitdist(X, distr = "lnorm", method = "mle")
  mu_sigma_hat <- fit.lnorm$estimate
  if(alpha == 0){
    return(mu_sigma_hat)
  }else{
    mu.log.sigma <- c(mu_sigma_hat[1], log(mu_sigma_hat[2]))
    robust_est <- function(par){
      (exp(-alpha * par[1] + 0.5 * (alpha^2) * exp(2 * par[2])/(1 + alpha)) / 
          ((exp(alpha * par[2]) * (sqrt(2 * pi)^alpha) * sqrt(alpha + 1))) -
          (1 + 1 / alpha) * mean(dlnorm(X, meanlog = par[1], sdlog = exp(par[2]))^alpha)) + 1 / alpha}
    optim.out <- suppressWarnings(optim(par = mu.log.sigma, robust_est))
    out <- optim.out$par
    out <- c(out[1], exp(out[2]))
    names(out) <- c("meanlog", "sdlog")
    return(out)}
}

#--------------------

mdpde.weibull <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  
  library(fitdistrplus)
  fit.weibull <- fitdist(X, distr = "weibull", method = "mle")
  shape_scale_hat <- fit.weibull$estimate
  shape_rate_hat <- c(shape_scale_hat[1], 1 / shape_scale_hat[2])
  names(shape_rate_hat) <- c("shape", "rate")
  if(alpha == 0){
    return(shape_rate_hat)
  }else{
    log.shape.rate <- log(shape_rate_hat)
    robust_est <- function(par){
      exp(par[1])^alpha * exp(par[2])^alpha * 
        gamma(1 + (exp(par[1]) - 1) * alpha / exp(par[1])) / 
        ((alpha + 1)^(1 + (exp(par[1]) - 1) * alpha / exp(par[1]))) - 
          (1 + 1 / alpha) * mean(dweibull(X, shape = exp(par[1]), 
                                          scale = exp(-par[2]))^alpha) + 1 / alpha}
    optim.out <- suppressWarnings(optim(par = log.shape.rate, robust_est))
    out <- optim.out$par
    names(out) <- c("shape", "rate")
    return(exp(out))}
}

#--------------------

WD.mdpde.exp <- function(X, alpha){
    X <- X[is.na(X) == 0]
    X <- X[X > 0]
    X <- sort(X)
    n <- length(X)
    
    mdpde.ests <- sapply(1:n, function(i){mdpde.exp(X[-i], alpha = alpha)})
    cdf.emp <- (seq_len(n) - 0.5) / n
    cdf.th <- sapply(1:n, function(i){pexp(X[i], rate = mdpde.ests[i])})
    mean(abs(cdf.emp - cdf.th))}

#--------------------

WD.mdpde.gamma <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  X <- sort(X)
  n <- length(X)
  
  mdpde.ests <- t(sapply(1:n, function(i){mdpde.gamma(X[-i], alpha = alpha)}))
  cdf.emp <- (seq_len(n) - 0.5) / n
  cdf.th <- sapply(1:n, function(i){
    pgamma(X[i], shape = mdpde.ests[i, 1], rate = mdpde.ests[i, 2])})
  mean(abs(cdf.emp - cdf.th))}

#--------------------

WD.mdpde.lnorm <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  X <- sort(X)
  n <- length(X)
  
  mdpde.ests <- t(sapply(1:n, function(i){mdpde.lnorm(X[-i], alpha = alpha)}))
  cdf.emp <- (seq_len(n) - 0.5) / n
  cdf.th <- sapply(1:n, function(i){
    plnorm(X[i], meanlog = mdpde.ests[i, 1], sdlog = mdpde.ests[i, 2])})
  mean(abs(cdf.emp - cdf.th))}

#--------------------

WD.mdpde.weibull <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  X <- sort(X)
  n <- length(X)
  
  mdpde.ests <- t(sapply(1:n, function(i){mdpde.weibull(X[-i], alpha = alpha)}))
  cdf.emp <- (seq_len(n) - 0.5) / n
  cdf.th <- sapply(1:n, function(i){
    pweibull(X[i], shape = mdpde.ests[i, 1], scale = 1 / mdpde.ests[i, 2])})
  mean(abs(cdf.emp - cdf.th))}

#--------------------

optim.alpha.exp <- function(X){
    out <- optimize(function(alpha){WD.mdpde.exp(X, alpha)}, lower = 0, upper = 1)
    out <- list(optim.alpha = out$minimum, min.WD = out$objective)
    out}

optim.alpha.gamma <- function(X){
  out <- optimize(function(alpha){WD.mdpde.gamma(X, alpha)}, lower = 0, upper = 1)
  out <- list(optim.alpha = out$minimum, min.WD = out$objective)
  out}

optim.alpha.lnorm <- function(X){
  out <- optimize(function(alpha){WD.mdpde.lnorm(X, alpha)}, lower = 0, upper = 1)
  out <- list(optim.alpha = out$minimum, min.WD = out$objective)
  out}

optim.alpha.weibull <- function(X){
  out <- optimize(function(alpha){WD.mdpde.weibull(X, alpha)}, lower = 0, upper = 1)
  out <- list(optim.alpha = out$minimum, min.WD = out$objective)
  out}

#--------------------

J.K.mdpde.exp <- function(rate, alpha){
  J <- (alpha^2 + 1) / (alpha + 1)^3 * rate^(alpha - 2)
  K <- ((4 * alpha^2 + 1) / (2 * alpha + 1)^3 - alpha^2 / (alpha + 1)^4) * rate^(2 * alpha - 2)
list(J = J, K = K)}

J.K.mdpde.gamma <- function(shape, rate, alpha, nrand = 1000, seed = 1){
  set.seed(seed)
  rsamp <- rgamma(nrand, shape = shape, rate = rate)
  u1 <- log(rate) - digamma(shape) + log(rsamp)
  u2 <- shape / rate - rsamp
  
  J <- matrix(NA, 2, 2)
  density.vals <- dgamma(rsamp, shape = shape, rate = rate)
  J[1, 1] <- mean(u1^2 * density.vals^alpha)
  J[1, 2] <- J[2, 1] <- mean(u1 * u2 * density.vals^alpha)
  J[2, 2] <- mean(u2^2 * density.vals^alpha)
  
  xi <- c(mean(u1 * density.vals^alpha),
          mean(u2 * density.vals^alpha))
  
  K <- matrix(NA, 2, 2)
  K[1, 1] <- mean(u1^2 * density.vals^(2 * alpha))
  K[1, 2] <- K[2, 1] <- mean(u1 * u2 * density.vals^(2 * alpha))
  K[2, 2] <- mean(u2^2 * density.vals^(2 * alpha))
  K <- K - tcrossprod(xi)
  list(J = J, K = K)}

J.K.mdpde.lnorm <- function(meanlog, sdlog, alpha, nrand = 1000, seed = 1){
  set.seed(seed)
  rsamp <- rlnorm(1e6, meanlog = meanlog, sdlog = sdlog)
  
  u1 <- (log(rsamp) - meanlog) / sdlog^2
  u2 <- -1/sdlog + (log(rsamp) - meanlog)^2 / sdlog^3
  
  J <- matrix(NA, 2, 2)
  J[1, 1] <- mean(u1^2 * dlnorm(rsamp, meanlog = meanlog, sdlog = sdlog)^alpha)
  J[1, 2] <- J[2, 1] <- mean(u1 * u2 * dlnorm(rsamp, meanlog = meanlog, sdlog = sdlog)^alpha)
  J[2, 2] <- mean(u2^2 * dlnorm(rsamp, meanlog = meanlog, sdlog = sdlog)^alpha)
  
  xi <- c(mean(u1 * dlnorm(rsamp, meanlog = meanlog, sdlog = sdlog)^alpha),
          mean(u2 * dlnorm(rsamp, meanlog = meanlog, sdlog = sdlog)^alpha))
  
  K <- matrix(NA, 2, 2)
  K[1, 1] <- mean(u1^2 * dlnorm(rsamp, meanlog = meanlog, sdlog = sdlog)^(2 * alpha))
  K[1, 2] <- K[2, 1] <- mean(u1 * u2 * dlnorm(rsamp, meanlog = meanlog, sdlog = sdlog)^(2 * alpha))
  K[2, 2] <- mean(u2^2 * dlnorm(rsamp, meanlog = meanlog, sdlog = sdlog)^(2 * alpha))
  K <- K - tcrossprod(xi)
  list(J = J, K = K)}

J.K.mdpde.weibull <- function(shape, rate, alpha, nrand = 1000, seed = 1){
  set.seed(seed)
  rsamp <- rweibull(nrand, shape = shape, scale = 1 / rate)
  
  u1 <- (1 / shape) + log(rate) + log(rsamp) - (rate * rsamp)^shape * (log(rate) + log(rsamp))
  u2 <- shape / rate - shape * (rate^(shape - 1)) * rsamp^shape
  
  J <- matrix(NA, 2, 2)
  J[1, 1] <- mean(u1^2 * dweibull(rsamp, shape = shape, scale = 1 / rate)^alpha)
  J[1, 2] <- J[2, 1] <- mean(u1 * u2 * dweibull(rsamp, shape = shape, scale = 1 / rate)^alpha)
  J[2, 2] <- mean(u2^2 * dweibull(rsamp, shape = shape, scale = 1 / rate)^alpha)
  
  xi <- c(mean(u1 * dweibull(rsamp, shape = shape, scale = 1 / rate)^alpha),
          mean(u2 * dweibull(rsamp, shape = shape, scale = 1 / rate)^alpha))
  
  K <- matrix(NA, 2, 2)
  K[1, 1] <- mean(u1^2 * dweibull(rsamp, shape = shape, scale = 1 / rate)^(2 * alpha))
  K[1, 2] <- K[2, 1] <- mean(u1 * u2 * dweibull(rsamp, shape = shape, scale = 1 / rate)^(2 * alpha))
  K[2, 2] <- mean(u2^2 * dweibull(rsamp, shape = shape, scale = 1 / rate)^(2 * alpha))
  K <- K - tcrossprod(xi)
  list(J = J, K = K)}

#--------------------

sd.mdpde.exp <- function(X, alpha, B = 1000){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  n <- length(X)
  set.seed(B)
  sd(apply(t(replicate(B, sample(X, n, replace = T))), 1, mdpde.exp, alpha = alpha))}

sd.mdpde.gamma <- function(X, alpha, B = 1000){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  n <- length(X)
  set.seed(B)
  apply(apply(t(replicate(B, sample(X, n, replace = T))), 1, mdpde.gamma, alpha = alpha), 1, sd)}

sd.mdpde.lnorm <- function(X, alpha, B = 1000){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  n <- length(X)
  set.seed(B)
  apply(apply(t(replicate(B, sample(X, n, replace = T))), 1, mdpde.lnorm, alpha = alpha), 1, sd)}

sd.mdpde.weibull <- function(X, alpha, B = 1000){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  n <- length(X)
  set.seed(B)
  apply(apply(t(replicate(B, sample(X, n, replace = T))), 1, mdpde.weibull, alpha = alpha), 1, sd)}

#--------------------

RIC.exp <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  n <- length(X)
  
  lambda.mdpde <- mdpde.exp(X, alpha = alpha)
  
  if(alpha == 0){
    print("AIC is returned as the tuning parameter is zero")
    H.mle <- lambda.mdpde * mean(X) - log(lambda.mdpde)
    out <- H.mle + 1 / n
  }else{
    H.mdpde <- (lambda.mdpde^alpha) * (1 / (1 + alpha) - (1 + 1 / alpha) * 
                                         mean(exp(-alpha * lambda.mdpde * X)))
    J.K <- J.K.mdpde.exp(lambda.mdpde, alpha = alpha)
    out <- H.mdpde + (J.K$K / J.K$J) / (1 + alpha) / n
  }
  out}

RIC.gamma <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  n <- length(X)
  
  shape.rate.mdpde <- mdpde.gamma(X, alpha = alpha)
  
  if(alpha == 0){
    print("AIC is returned as the tuning parameter is zero")
    H.mle <- -mean(dgamma(X, shape = shape.rate.mdpde[1], rate = shape.rate.mdpde[2], log = T))
    out <- H.mle + 2 / n
  }else{
    H.mdpde <- (shape.rate.mdpde[2]^alpha) * gamma(1 + (shape.rate.mdpde[1] - 1) * (alpha + 1)) / 
      gamma(shape.rate.mdpde[1])^(1 + alpha) /
      (1 + alpha)^(1 + (shape.rate.mdpde[1] - 1) * (alpha + 1)) - 
      (1 + 1 / alpha) * mean(dgamma(X, shape = shape.rate.mdpde[1], 
                                    rate = shape.rate.mdpde[2])^alpha)
    J.K <- J.K.mdpde.gamma(shape = shape.rate.mdpde[1], rate = shape.rate.mdpde[2], alpha = alpha)
    out <- as.vector(H.mdpde + sum(diag(solve(J.K$J) %*% J.K$K)) / (1 + alpha) / n)
  }
  out}

RIC.lnorm <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  n <- length(X)
  
  mu.sigma.mdpde <- mdpde.lnorm(X, alpha = alpha)
  
  if(alpha == 0){
    print("AIC is returned as the tuning parameter is zero")
    H.mle <- -mean(dlnorm(X, meanlog = mu.sigma.mdpde[1], sdlog = mu.sigma.mdpde[2], log = T))
    out <- H.mle + 2 / n
  }else{
    H.mdpde <- exp(-alpha * mu.sigma.mdpde[1] + 0.5 * (alpha^2) * mu.sigma.mdpde[2]^2 /(1 + alpha)) / 
       ((mu.sigma.mdpde[2]^alpha * (sqrt(2 * pi)^alpha) * sqrt(alpha + 1))) -
       (1 + 1 / alpha) * mean(dlnorm(X, meanlog = mu.sigma.mdpde[1], sdlog = mu.sigma.mdpde[2])^alpha)
    J.K <- J.K.mdpde.lnorm(meanlog = mu.sigma.mdpde[1], sdlog = mu.sigma.mdpde[2], alpha = alpha)
    out <- as.vector(H.mdpde + sum(diag(solve(J.K$J) %*% J.K$K)) / (1 + alpha) / n)
  }
  out}

RIC.weibull <- function(X, alpha){
  X <- X[is.na(X) == 0]
  X <- X[X > 0]
  n <- length(X)
  
  shape.rate.mdpde <- mdpde.weibull(X, alpha = alpha)
  
  if(alpha == 0){
    print("AIC is returned as the tuning parameter is zero")
    H.mle <- -mean(dweibull(X, shape = shape.rate.mdpde[1], scale = 1 / shape.rate.mdpde[2], log = T))
    out <- H.mle + 2 / n
  }else{
    H.mdpde <- shape.rate.mdpde[1]^alpha * shape.rate.mdpde[2]^alpha * 
      gamma(1 + (shape.rate.mdpde[1] - 1) * alpha / shape.rate.mdpde[1]) / 
      ((alpha + 1)^(1 + (shape.rate.mdpde[1] - 1) * alpha / shape.rate.mdpde[1])) - 
      (1 + 1 / alpha) * mean(dweibull(X, shape = shape.rate.mdpde[1], scale = 1 / shape.rate.mdpde[2])^alpha)
    J.K <- J.K.mdpde.weibull(shape = shape.rate.mdpde[1], rate = shape.rate.mdpde[2], alpha = alpha)
    out <- as.vector(H.mdpde + sum(diag(solve(J.K$J) %*% J.K$K)) / (1 + alpha) / n)
  }
  out}



