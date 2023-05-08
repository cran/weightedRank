dwgtRank <-
function (y, gamma=1, m=2, m1=2, m2=2, phifunc = NULL,
                    alternative="greater", scores=NULL, range=TRUE)
{
  # Adjusted from wgtRank in weightedRank package, version 0.1.6
  # The four specified phi functions are removed and replaced by m's
  # Within sets scores are now used if scores is specified
  # A lower tailed test is now an explicit option
  # Within-rank scores are now an option
  # The gap may now optionally replace with within-block range

  # Check input
  stopifnot(is.matrix(y) | is.data.frame(y))
  if (!is.null(scores)) stopifnot(length(scores)==(dim(y)[2]))
  stopifnot(0 == sum(is.na(as.vector(y))))
  stopifnot(min(dim(y)) >= 2)
  stopifnot(is.vector(gamma) & (length(gamma) == 1)
            & (gamma >=  1))
  stopifnot((alternative=="greater")|(alternative=="less"))

  # Define several functions that will be used in computations

  # This is the function used to score the ranks of the blocks.
  # It is given by the right side of expression (9) in
  # Rosenbaum (2011) Biometrics, 67, 1017–1027
  multrnksU <- function(pk, m1 = 2, m2 = 2, m = 2) {
    n <- length(pk)
    q <- rep(0, n)
    q <- rep(0, n)
    for (l in m1:m2) {
      q <- q + (l * choose(m, l) * (pk^(l - 1)) * ((1 -
                                                      pk)^(m - l)))
    }
    q/max(q)
  }

  # If the user did not specify phifunc, then use multrnksU instead
  if (is.null(phifunc)) {
    stopifnot((m1>=1)&(m2>=m1)&(m>=m2))
    phifunc<-function(pk){multrnksU(pk,m1=m1,m2=m2,m=m)}
  }

  # This function computes the asymptotic separable approximation
  # to the upper bound on the P-value, as discussed in
  # Gastwirth, Krieger and Rosenbaum (2000) JRSS-B
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

  # Begin computations
  J <- dim(y)[2]
  nset <- dim(y)[1]

  # Use the range or the gap?
  if (range) rg <- apply(y, 1, max) - apply(y, 1, min)
  else {
    mx <- apply(y, 1, max)
    sm <- apply(y, 1, sum) - mx
    rg <- mx - (sm/(J-1))
  }

  # Score the block ranks
  rkrg <- phifunc(rank(rg)/nset)

  # Determine the within-block ranks
  if (is.null(scores)){
    rky <- t(apply(y, 1, rank)) # uses average ranks
  }
  if (!is.null(scores)){
    rky<-t(apply(y, 1, rank, ties.method = "min")) # uses min ranks
    orky<-rky
    for (i in 1:J)  rky[orky==i] <- scores[i]
  }

  if (alternative=="less") rky<-(-rky)

  mv <- separable1kA(rky, gamma = gamma)
  ts <- mean(rky[, 1] * rkrg)
  ex <- sum(mv$maxmu * rkrg)/nset
  va <- (sum(mv$maxsig2 * rkrg * rkrg)/nset)/nset
  dev <- (ts - ex)/sqrt(va)
  pval <- 1 - stats::pnorm(dev)
  if (alternative=="less"){
    ts<-(-ts)
    ex<-(-ex)
    dev<-(-dev)
  }

  detail <- c(dev, ts, ex, va, gamma)
  names(detail) <- c("Deviate", "Statistic", "Expectation",
                     "Variance", "Gamma")
  list(pval = pval, detail = detail)
}
