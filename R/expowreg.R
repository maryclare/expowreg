# library(tmvtnorm)
library(truncdist) # Use to sample from truncated gamma directly
library(actuar) # Use to access inverse gamma probability function from truncdist functons
library(BayesBridge) # Use to sample from truncated univ. normal

sample.sigma.sq <- function(y, X, beta, q, a.pr) {
  n <- length(y)
  p <- ncol(X)
  bb <- crossprod(y - crossprod(t(X), beta))/2 + b.pr
  aa <- (n + 1)/2 + a.pr
  return(rinvgamma(1, aa, 1/bb))
}

sample.tau.sq <- function(beta, q) {
  p <- length(beta)
  aa <- (p + 1)/q + 1
  bb <- sum(abs(beta*(gamma(1/q)/gamma(3/q))^(-1/2))^q)
  return((rinvgamma(1, aa, 1/bb))^(2/q))
}

log.lik <- function(y, X, beta, tau.sq, sig.sq, q, samp.sig.sq, samp.tau.sq) {
  p <- ncol(X)
  n <- nrow(X)
  r <- y - crossprod(t(X), beta)
  dnorm(crossprod(r),
        sd = sqrt(sig.sq), log = TRUE) + lpbeta(beta, tau.sq, q) +
    samp.tau.sq*dgamma(tau.sq^(-q/2), (p + 1)/q + 1,
           (gamma(1/q)/gamma(3/q))^(-q/2)*sum(abs(beta)^q),
           log = TRUE) +
    samp.sig.sq*dgamma(1/sig.sq, (n + 1)/2 + 1,
           crossprod(r)/2, log = TRUE)
}

# Found method for MVNORM simulation in Polson, Scott, Windle (2014)
samp.singtmvnorm <- function(beta.start, DUty, delta, d, Vt, sig.sq, W) {

  p <- length(delta)
  z <- crossprod(t(Vt), beta.start)
  for (i in 1:p) {

    right <- rep(delta, 2) - W[, -i]%*%z[-i]
    left <- W[, i]
    upp <- (right/left)[left > 0]
    low <- (right/left)[left < 0]

    upp.lim <- min(upp)
    low.lim <- max(low)

    mm <- DUty[i]/d[i]^2
    ss <- sig.sq/d[i]^2

    if (d[i] > 10^(-16)) {
      z[i] <- BayesBridge::rtnorm(1, left = low.lim, right = upp.lim,
                                  mu = mm, sig = sqrt(ss))
    } else {
      z[i] <- runif(1, low.lim, upp.lim)
    }

  }
  beta <- crossprod(Vt, z)
  return(as.numeric(beta))
}

lpbeta <- function(beta, tau.sq, q) {
  p <- length(beta)
  p*log(sqrt(gamma(3/q))*q) -
    sum((gamma(3/q)/gamma(1/q))^(q/2)*abs((beta)/sqrt(tau.sq))^q) -
    p*log(2*sqrt(tau.sq)*sqrt(gamma(1/q)^3))
}

sample.beta <- function(A, b, gamma, sig.sq, tau.sq, q, beta.old, sing.A,
                        d, Vt, DUty, W, diag.A) {

  delta <- sqrt(gamma(1/q)/gamma(3/q))*sqrt(tau.sq/2)*gamma^(1/q)

  # Check that design matrix is diagonal, if not use samp.signtmvnorm
  # otherwise can just sample from truncated univariate normals
  if (!diag.A) {
    beta <- samp.singtmvnorm(beta.start = beta.old, DUty = DUty,
                             delta = delta, d = d,
                             Vt = Vt,
                             sig.sq = sig.sq, W = W)
  } else {
    beta <- numeric(length(gamma))
    too.small <- sqrt(diag(A)/sig.sq) <= 10^(-16)
    beta[!too.small] <- BayesBridge::rtnorm(sum(!too.small),
                                            left = -delta[!too.small], right = delta[!too.small],
                                            mu = b[!too.small], sig = sqrt(sig.sq/diag(A)[!too.small]))

    beta[too.small] <- runif(sum(too.small), -delta[too.small], delta[too.small])

  }

  return(beta)
}

sample.gamma <- function(beta, tau.sq, q) {
  eta <- (sqrt(gamma(3/q)/gamma(1/q))*sqrt(2/tau.sq)*abs(beta))^q
  if (TRUE) { # if (q <= 4) {
    gamma.inv <- numeric(length(beta))
    for (i in 1:length(beta)) {
      lim.a <- 0
      lim.b <- 1/eta[i]
      p.a <- pinvgamma(lim.a, shape = 1 + 1/q, scale = 2^(q/2))
      p.b <- pinvgamma(lim.b, shape = 1 + 1/q, scale = 2^(q/2))
      if (p.a != p.b) {
        gamma.inv[i] <- rtrunc(1, spec = "invgamma", a = lim.a, b = lim.b,
                               shape = 1 + 1/q, scale = 2^(q/2))
      } else {
        gamma.inv[i] <- runif(1, lim.a, lim.b)
      }
    }
    gamma <- 1/gamma.inv
  } else {
    gamma <- numeric(length(beta))
    for (i in 1:length(beta)) {

      lim.a <- eta[i]
      lim.b <- Inf
      p.a <- pgamma(lim.a, shape = 1 + 1/q, rate = 2^(-q/2))
      p.b <- pgamma(lim.b, shape = 1 + 1/q, rate = 2^(-q/2))
      if (p.a != p.b) {
      gamma[i] <- rtrunc(1, spec = "gamma", a = lim.a, b = lim.b,
                         shape = 1 + 1/q, rate = 2^(-q/2))
      } else {
        gamma[i] <- runif(1, lim.a, 10^(100))
      }
    }
  }
  return(gamma)
}
#' @export
epr.sampler <- function(X, y,
                       num.samp = 1000, print.iter = FALSE,
                       num.chains = 1,
                       q,
                       sig.sq = NULL, tau.sq = NULL,
                       burn.in = 500, comp.lik = TRUE,
                       b.start = NULL, a.pr = a.pr, b.pr = b.pr) {
  p <- ncol(X)
  svd.X <- svd(X, nv = ncol(X))
  U <- svd.X$u
  Vt <- t(svd.X$v)
  d <- c(svd.X$d, rep(0, nrow(Vt) - length(svd.X$d)))

  DUty <- crossprod(crossprod(t(U), diag(svd.X$d)), y)
  W <- crossprod(t(rbind(diag(rep(1, p)), diag(rep(-1, p)))), t(Vt))
  A <- tcrossprod(tcrossprod(Vt, diag(d)))
  sing.A <- length(svd.X$d) < ncol(X)
  diag.A <- sum(abs(A[lower.tri(A, diag = FALSE)]) <= 10^(-14)) == p*(p - 1)/2

  if (sing.A) {
    del <- (1 - min(eigen(A)$values))
    b <- as.numeric(crossprod(solve(A + del*diag(ncol(X))), crossprod(X, y)))
  } else {
    b <- as.numeric(crossprod(solve(A), crossprod(X, y)))
  }
  if (!is.null(b.start)) {
    b <- b.start
  }

  betas <- gammas <- matrix(nrow = num.samp + 1, ncol = ncol(X))
  sig.sqs <- tau.sqs <- lls <- numeric(num.samp + 1)
  if (!is.null(tau.sq)) {
    sig.sqs <- rep(sig.sq, num.samp + 1)
  } else {
    sig.sqs[1] <- 1
  }
  if (!is.null(tau.sq)) {
    tau.sqs <- rep(tau.sq, num.samp + 1)
  } else {
    tau.sqs[1] <- 1
  }
  betas[1, ] <- b
  if (FALSE) { # Need to be careful about starting values when q is large, as q gets large range of beta gets restricted by prior
    betas[1, ] <- 0.001*rnorm(p)
  } else {
    betas[1, ] <- b
  }

  for (i in 2:(num.samp + 1)) {
    if (print.iter) {cat("i = ", i, "\n")}

    if (is.null(tau.sq)) {
      tau.sqs[i] <- sample.tau.sq(beta = betas[i - 1, ], q = q)
    }
    if (is.null(sig.sq)) {
      sig.sqs[i] <- sample.sigma.sq(y = y, X = X, betas[i - 1, ], q = q,
                                    a.pr = a.pr, b.pr = b.pr)
    }
    gammas[i, ] <- sample.gamma(beta = betas[i - 1, ], tau.sq = tau.sqs[i], q = q)
    # betas[i, ] <- sample.beta(A = A, b = b, gamma = gammas[i, ],
    #                           sig.sq = sig.sqs[i],
    #                           tau.sq =  tau.sqs[i], q = q,
    #                           beta.old = betas[i - 1, ], sing.A = sing.A, d = d,
    #                           Vt = Vt, DUty = DUty, W = W, diag.A = diag.A)
    delta <- sqrt(gamma(1/q)/gamma(3/q))*sqrt(tau.sqs[i]/2)*gammas[i, ]^(1/q)
    betas[i, ] <- powreg::sampleBeta(start = betas[i - 1, ],
                                     DUty = DUty, delta = delta, d = d, Vt = Vt,
                                     sigsq = sig.sq, W = W)
    if (comp.lik) {
      lls[i] <- log.lik(y = y, X = X, beta = betas[i, ], tau.sq = tau.sqs[i],
                        sig.sq = sig.sqs[i], q = q, samp.sig.sq = is.null(sig.sq),
                        samp.tau.sq = is.null(tau.sq))
    }
  }
  return(list("beta" = betas[!1:(num.samp + 1) %in% 1:burn.in, ],
              "gamma" = gammas[!1:(num.samp + 1) %in% 1:burn.in, ],
              "sig.sq" = sig.sqs[!1:(num.samp + 1) %in% 1:burn.in],
              "tau.sq" = tau.sqs[!1:(num.samp + 1) %in% 1:burn.in],
              "ll" = lls[!1:(num.samp + 1) %in% 1:burn.in]))
}
