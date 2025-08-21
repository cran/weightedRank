wgtRank<-function (y, phi = "u868", phifunc = NULL, gamma = 1)
{
  stopifnot(is.matrix(y) | is.data.frame(y))
  stopifnot(0 == sum(is.na(as.vector(y))))
  stopifnot(min(dim(y)) >= 2)
  stopifnot(is.vector(gamma) & (length(gamma) == 1) & (gamma >=
                                                         1))
  if (is.null(phifunc)) {
    stopifnot(is.element(phi, c("u868", "u878", "u888", "u858",
                                "quade", "wilc", "mixed")))
  }
  multrnksU <- function(pk, m1 = 2, m2 = 2, m = 2) {
    n <- length(pk)
    q <- rep(0, n)
    q <- rep(0, n)
    for (l in m1:m2) {
      q <- q + (l * choose(m, l) * (pk^(l - 1)) * ((1 -
                                                      pk)^(m - l)))
    }
    q
  }
  u868 <- function(pk) {
    q<-multrnksU(pk, m1 = 6, m2 = 8, m = 8)
    q/max(q)
  }
  u878 <- function(pk) {
    q<-multrnksU(pk, m1 = 7, m2 = 8, m = 8)
    q/max(q)
  }
  u888 <- function(pk) {
    q<-multrnksU(pk, m1 = 8, m2 = 8, m = 8)
    q/max(q)
  }
  u858 <- function(pk) {
    q<-multrnksU(pk, m1 = 5, m2 = 8, m = 8)
    q/max(q)
  }
  quade <- function(pk) {
    q<-pk
    q/max(q)
  }
  wilc <- function(pk) {
    q<-rep(1, length(pk))
    q/max(q)
  }
  mixed <- function(pk) {
    q<-multrnksU(pk, m1 = 19, m2 = 20, m = 20)+
      multrnksU(pk, m1 = 19, m2 = 19, m = 20)
    q/max(q)
  }
  if (is.null(phifunc)) {
    if (phi == "u868")
      phifunc <- u868
    else if (phi == "u878")
      phifunc <- u878
    else if (phi == "quade")
      phifunc <- quade
    else if (phi == "wilc")
      phifunc <- wilc
    else if (phi == "u888")
      phifunc <- u888
    else if (phi == "u858")
      phifunc <- u858
    else if (phi == "mixed")
      phifunc <- mixed
  }
  separable1kA <- function(ymat, gamma = 1) {
    n <- dim(ymat)[1]
    m <- dim(ymat)[2]
    o <- t(apply(ymat, 1, sort))
    allmu <- matrix(NA, n, m - 1)
    allsigma2 <- matrix(NA, n, m - 1)
    maxmu <- rep(-Inf, n)
    maxsig2 <- rep(-Inf, n)
    for (j in 1:(m - 1)) {
      pr <- c(rep(1, j), rep(gamma, m - j))/(j + ((m -
                                                     j) * gamma))
      mu <- as.vector(o %*% pr)
      sigma2 <- as.vector((o * o) %*% pr) - (mu * mu)
      chgmu <- (mu > maxmu)
      samemu <- (mu == maxmu)
      if (sum(chgmu) > 0) {
        maxmu[chgmu] <- mu[chgmu]
        maxsig2[chgmu] <- sigma2[chgmu]
      }
      if (sum(samemu) > 0) {
        maxsig2[samemu] <- pmax(sigma2[samemu], maxsig2[samemu])
      }
    }
    list(maxmu = maxmu, maxsig2 = maxsig2)
  }
  J <- dim(y)[2]
  nset <- dim(y)[1]
  rg <- apply(y, 1, max) - apply(y, 1, min)
  rkrg <- phifunc(rank(rg)/nset)
  rky <- t(apply(y, 1, rank))
  mv <- separable1kA(rky, gamma = gamma)
  ts <- mean(rky[, 1] * rkrg)
  ex <- sum(mv$maxmu * rkrg)/nset
  va <- sum(mv$maxsig2 * rkrg * rkrg)/(nset * nset)
  dev <- (ts - ex)/sqrt(va)
  pval <- 1 - stats::pnorm(dev)
  detail <- c(dev, ts, ex, va, gamma)
  names(detail) <- c("Deviate", "Statistic", "Expectation",
                     "Variance", "Gamma")
  list(pval = pval, detail = detail)
}
