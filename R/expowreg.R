library(powreg)

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
      z[i] <- powreg::rtnormrej(l = low.lim, r = upp.lim,
                                  mu = mm, sd = sqrt(ss))
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
    beta[!too.small] <- powreg::rtnormrej(l = -delta[!too.small], r = delta[!too.small],
                                            mu = b[!too.small], sd = sqrt(sig.sq/diag(A)[!too.small]))

    beta[too.small] <- runif(sum(too.small), -delta[too.small], delta[too.small])

  }

  return(beta)
}

sample.gamma <- function(beta, tau.sq, q) {
  etaq <- (sqrt(gamma(3/q)/gamma(1/q))*sqrt(2/tau.sq)*abs(beta))^q
  d <- 2^(-q/2)
  gamma <- rexp(length(beta), d) + etaq
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
    betas[i, ] <- sample.beta(A = A, b = b, gamma = gammas[i, ],
                              sig.sq = sig.sqs[i],
                              tau.sq =  tau.sqs[i], q = q,
                              beta.old = betas[i - 1, ], sing.A = sing.A, d = d,
                              Vt = Vt, DUty = DUty, W = W, diag.A = diag.A)
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

# Function for estimating variance parameters from likelihood
rrmmle<-function(y,X,emu=FALSE,s20=1,t20=1)
{
  sX<-svd(X, nu = nrow(X), nv = ncol(X))
  lX<-sX$d^2
  tUX<-t(sX$u)
  xs<-apply(X,1,sum)

  if(nrow(X)>ncol(X))
  {
    lX<-c(lX,rep(0,nrow(X)-ncol(X)))
  }

  objective<-function(ls2t2,emu, lX, xs, tUX, y)
  {
    s2<-exp(ls2t2[1]) ; t2<-exp(ls2t2[2])
    mu<-emu*sum((tUX%*%xs)*(tUX%*%y)/(lX*t2+s2))/sum((tUX%*%xs)^2/(lX*t2+s2))
    ev <- lX*t2 + s2
    lev <- log(ifelse(ev > 10^(-300), ev, 1))
    squares <- ifelse(is.infinite((tUX%*%(y-mu*xs))^2/ev), 10^(300),
                      (tUX%*%(y-mu*xs))^2/ev)
    # Useful cat statements for debugging
    # cat("s2=", s2, " t2=", s2, "\n")
    # cat("slev=", sum(lev), " ssquares=", sum(squares), "\n")
    sum(lev) + sum(squares) # + 1/s2^(1/2) # Could keep s2 = 0 from being a mode but is a pretty artificial fix

  }


  # I include an upper bound to keep this from misbehaving, sometimes when
  # s2 and t2 get really small, the gradient goes crazy and a huge
  # value of s2 and t2 can get suggested
  fit<-optim(log(c(s20,t20)),objective,emu=emu, method = "L-BFGS-B", lX = lX,
             xs = xs, tUX = tUX, y = y, upper = rep(log(10^(30)), 2))
  s2<-exp(fit$par[1]) ; t2<-exp(fit$par[2])
  mu<-emu*sum((tUX%*%xs)*(tUX%*%y)/(lX*t2+s2))/sum((tUX%*%xs)^2/(lX*t2+s2))
  if (fit$convergence == 0) {
    return(c(mu,t2,s2))
  } else {
    return(rep(NA, 3))
  }
}

fq <- function(q, kurt) {
  gamma(5/q)*gamma(1/q)/(gamma(3/q)^2) - kurt
}

fpq <- function(q) {
  (gamma(5/q)*gamma(1/q)/(q^2*gamma(3/q)^2))*(6*digamma(3/q) - digamma(1/q) - 5*digamma(5/q))
}

# Use Newton's method: https://en.wikipedia.org/wiki/Newton%27s_method
#' @export
nrq <- function(kurt, sval = 0.032, tol = 10^(-12)) { # This starting value is the lowest possible
  # Kurtosis is bounded below by 1.8, so round if needed
  kurt <- ifelse(kurt <= 1.8, 1.81, kurt)
  # Kurtosis greater than 1.8 gives a q value of 1086.091
  # Value of fpq at q = 1086.091 is about -10^(-8), so the curve *is* pretty flat at this point
  if (kurt < 6) {
    sval <- 1
  } else if (kurt < 3) {
    sval <- 2
  }
  x.old <- Inf
  x.new <- sval
  while (abs(x.new - x.old) > tol) {
    x.old <- x.new
    x.new <- x.old - fq(x.old, kurt)/fpq(x.old)
  }
  return(x.new)
}
#' Function for estimating tuning parameters
#'
#' \code{estRegPars}
#'
#' @param \code{y} regression response
#' @param \code{X} regression design matrix
#' @param \code{delta} ridge regression parameter for when X is not full rank
#' @return Estimates \code{sigma.beta.sq.hat}, \code{sigma.epsi.sq.hat} and \code{kappa.hat}
#' @export
estRegPars <-function(y, X, delta.sq = NULL, precomp = NULL, comp.q = FALSE, mom = TRUE) {


  p <- ncol(X)
  n <- nrow(X)

  XtX <- crossprod(X)
  C <- cov2cor(XtX)
  V <- diag(sqrt(diag(XtX/C)))
  ceval <- eigen(C)$values
  if (!is.null(delta.sq)) {
    delta.sq <- delta.sq
  } else {
    delta.sq <- max(1 - min(ceval), 0)
  }
  D.inv <- tcrossprod(crossprod(V, (C + delta.sq*diag(p))), V)
  D <- solve(D.inv)
  DXtX <- crossprod(D, XtX)
  XD <- crossprod(t(X), D)
  DXtXD <- crossprod(XD)

  b <- crossprod(D, crossprod(X, y))

  if (mom) {

    XXt <- tcrossprod(X)
    XXt.eig <- eigen(XXt)
    XXt.val <- XXt.eig$values
    XXt.vec <- XXt.eig$vectors
    A <- diag(n)
    # B <- tcrossprod(tcrossprod(XXt.vec, diag(ifelse(XXt.val > 1, 1/XXt.val, 0))), XXt.vec)
    B <- tcrossprod(tcrossprod(XXt.vec, diag(1/(XXt.val + 1))), XXt.vec)

    E <- rbind(c(sum(diag(crossprod(t(XXt), A))), sum(diag(A))),
               c(sum(diag(crossprod(t(XXt), B))), sum(diag(B))))

    sig.2.ests <- solve(E)%*%matrix(c(crossprod(y, crossprod(t(A), y)),
                                      crossprod(y, crossprod(t(B), y))), nrow = 2, ncol = 1)

    sigma.beta.sq.hat <- sig.2.ests[1]
    sigma.epsi.sq.hat <- sig.2.ests[2]

  } else {
    vpars <- rrmmle(y = y, X = X)
    sigma.beta.sq.hat <- vpars[2]
    sigma.epsi.sq.hat <- vpars[3]
  }


  alpha.beta <- sum(diag(crossprod(DXtX)))/p
  gamma.beta <- sum(diag(DXtX^4))/p
  omega.beta <- 3*(sum(diag(crossprod(DXtX)^2)) - sum(diag(DXtX^4)))/p

  test.stat <- (mean(b^4))/(mean(b^2)^2)

  kappa.hat <- (alpha.beta^2/gamma.beta)*(test.stat - omega.beta/alpha.beta^2)
  q.hat <- ifelse(comp.q, nrq(kappa.hat), NA)

  return(list("sigma.beta.sq.hat" = sigma.beta.sq.hat,
              "sigma.epsi.sq.hat" = sigma.epsi.sq.hat,
              "kappa.hat" = kappa.hat,
              "q.hat" = q.hat,
              "test.stat" = test.stat,
              "DXtX" = DXtX,
              "DXtXD" = DXtXD,
              "delta.sq" = delta.sq))

}



