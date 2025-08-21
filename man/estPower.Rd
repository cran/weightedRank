\name{estPower}
\alias{estPower}

\title{
   Estimate the Power of a Sensitivity Analysis
}
\description{
  Uses the jackknife to estimate the power of a sensitivity analysis based on pilot data for an I x J observational block design in which each of I blocks contains one treated individual and J-1 controls.
}
\usage{
estPower(y,gammas,ssratio=1,phi="u868",alpha=0.05)
}
\arguments{
  \item{y}{
 A matrix or data frame with I rows and J columns containing the pilot data.  Column 1 contains the response of the treated individuals and columns 2 throught J contain the responses of controls in the same block.  An error will result if y contains NAs.
}
\item{gammas}{A vector of values of the sensitivity parameter gamma at which power will be estimated.  All values in gammas must be greater than or equal to 1.}
\item{ssratio}{The sample size ratio.  One positive value indicating the ratio of the sample size at which power is to be computed divided by the sample size of the pilot data.  If the pilot data contain 100 blocks, and power is to be computed for a study with 1000 blocks, then ssratio is 10 = 1000/100.}
\item{phi}{
     The weight function to be applied to the ranks of the within block ranges. The options are the same as the corresponding options in the wgtRank() function in this package. The options are: (i) "wilc" for the stratified Wilcoxon test, which gives every block the same weight, (ii) "quade" which ranks the within block ranges from 1 to I, and is closely related to Quade's (1979) statistic; see also Tardif (1987), (iii) "u858", "u868", "u878", "u888" or "mixed" based on Rosenbaum (2011, 2025a,b).
}
\item{alpha}{The level of the test, conventionally alpha=0.05.  Power is computed for a one-sided level-alpha test that rejects the hypothesis of no effect when the weighted rank statistic is large.  Commonly, the power of a two-sided 0.05-level test is computed by setting alpha=0.025 and viewing a two-sided test as two one-sided tests with a Bonferroni correction for conducting two tests; see Cox (1977).}
}
\details{
For a weighted rank statistic of the type computed by wgtRank() in this package, the estPower() function uses the jackknife to estimate its expecation and variance under the alternative hypothesis that the blocks are an iid sample from a population of blocks in which there is a treatment effect and there is no bias from unmeasured covariates.  Then, ssratio is used to adjust the estimated variance to a possibly different sample size, and power is estimated using a Normal approximation.

The method is discussed in Rosenbaum (2025b) where it is observed that the weighted rank statistics in Rosenbaum (2024) are U-statistics that are asymptotically Normal under both the null and alternative hypotheses, thereby justifying the use of a Normal approximation in the estimate of power.

One should not estimate the power for several test statistics using y, and then analyze y using the most powerful test, because that uses the data twice.  Instead, use adaptive inference in the wgtRanktt() function in this package; see Rosenbaum (2024, Section 5), Rosenbaum (2025a, Section 9.5) or
Rosenbaum (2025b, Section 4.4).  Alternatively, estimate the power using estPower for pilot data provided by split samples; see Heller et al. (2009).

For general discussion of the jackknife, see Chapter 11 of Efron and Tibshirani (1993).  For discussion of the jackknife applied to a U-statistic, see Arvesen (1969) and Chapter 5 of Lee (1990).
}
\value{
 \item{power }{A vector with the same length as gammas, where power[i] is the estimated power in a sensitivity analysis with gamma set to gammas[i].}
 \item{jackm }{The jackknife estimate of the expectation of the test statistic in the pilot data.}
 \item{jackv }{The jackknife estimate of the variance of the test statistic in the pilot data.}
}
\references{

Arvesen, J. N. (1969) Jackknifing U-statistics.  Annals of Mathematical Statistics, 40, 2076-2100.

Cox, D. R. (1977). The role of significance tests [with discussion and reply]. Scandinavian Journal of Statistics, 4, 49-70.

Efron, B. and Tibshirani, R. J. (1993) An Introduction to the Bootstrap.  Chapman and Hall/CRC, Boca Raton, FL.

Gastwirth, J. L., Krieger, A. M., and Rosenbaum, P. R. (2000).
<doi:10.1111/1467-9868.00249> Asymptotic separability in sensitivity analysis.  Journal of the Royal Statistical Society B 2000, 62, 545-556.

Heller, R., Small, D. S. and Rosenbaum, P. R. (2009)
<doi:10.1198/jasa.2009.tm08338> Split samples and design sensitivity in observational studies.  Journal of the American Statistical Association,
104, 1090-1101.

Lee, A. J. (1990) U-Statistics: Theory and Practice.  Chapman and Hall/CRC, Boca Raton, FL.

Lehmann, E. L. (1975). Nonparametrics: Statistical Methods Based on Ranks. San Francisco: Holden-Day.

Quade, D. (1979). <doi:10.2307/2286991> Using weighted rankings in the analysis of complete blocks with additive block effects. Journal of the American Statistical Association, 74, 680-683.

Rosenbaum, P. R. (1987). <doi:10.2307/2336017> Sensitivity analysis for certain permutation inferences in matched observational studies. Biometrika, 74(1), 13-26.

Rosenbaum, P. R. (2011). <doi:10.1111/j.1541-0420.2010.01535.x> A new U‐Statistic with superior design sensitivity in matched observational studies. Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2013). <doi:10.1111/j.1541-0420.2012.01821.x> Impact of multiple matched controls on design sensitivity in observational studies.  Biometrics, 2013, 69, 118-127.

Rosenbaum, P. R. (2014) <doi:10.1080/01621459.2013.879261> Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls.  Journal of the American Statistical Association, 109(507), 1145-1158.

Rosenbaum, P. R. (2015). <doi:10.1080/01621459.2014.960968> Bahadur efficiency of sensitivity analyses in observational studies. Journal of the American Statistical Association, 110(509), 205-217.

Rosenbaum, P. (2017). <doi:10.4159/9780674982697> Observation and Experiment: An Introduction to Causal Inference.  Cambridge, MA: Harvard University Press.

Rosenbaum, P. R. (2018). <doi:10.1214/18-AOAS1153> Sensitivity analysis for stratified comparisons in an observational study of the effect of smoking on homocysteine
levels. The Annals of Applied Statistics, 12(4), 2312-2334.

Rosenbaum, P. R. (2024) <doi:10.1080/01621459.2023.2221402> Bahadur efficiency of observational block designs. Journal of the American Statistical Association,
119, 1871-1881.

Rosenbaum, P. R. (2025a) <doi:10.1007/978-3-031-90494-3> An Introduction to the Theory of Observational Studies.  Switzerland: Springer.

Rosenbaum, P. R. (2025b) A jackknife estimate of the power of a sensitivity analysis in an observational study.  Manuscript.

Tardif, S. (1987). <doi:10.2307/2289476> Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82(398), 637-644.
}
\author{
Paul R. Rosenbaum
}

\examples{
data(aHDL)
hdl0<-t(matrix(aHDL$hdl,4,406))

\donttest{


# The code that follows reproduces the computations and figures in
# Rosenbaum (2025b), including the plots of power and the plots of
# phi-functions.

gammas<-(10:110)/10

plot(gammas,c(0,rep(1,length(gammas)-1)),type="n",
     cex.axis=.7,xlim=c(min(gammas)-.5,max(gammas)+.2),
     ylab="Power",xlab=expression(Gamma),las=1,
     cex.lab=.8,main="I=406",cex.main=.8)
legend(0,.85,c("Wilcoxon","Quade","U868","U878","U888","Mixed"),
       col=c("gray","gray","black","black","black","black"),
       lty=c(2,1,3,2,4,1),lwd=c(2,2,2,1,1,1),cex=.6,bty="n")
lines(gammas,estPower(hdl0,gammas,ssratio=1,phi="wilc")$power,
      col="gray",lty=2,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=1,phi="quade")$power,
      col="gray",lty=1,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=1,phi="u868")$power,
      col="black",lty=3,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=1,phi="u878")$power,
      col="black",lty=2,lwd=1)
lines(gammas,estPower(hdl0,gammas,ssratio=1,phi="u888")$power,
      col="black",lty=4,lwd=1)
lines(gammas,estPower(hdl0,gammas,ssratio=1,phi="mixed")$power,
      col="black",lty=1,lwd=1)


gammas<-(10:110)/10
ssratio<-1000/406


plot(gammas,c(0,rep(1,length(gammas)-1)),type="n",
     cex.axis=.7,xlim=c(min(gammas)-.5,max(gammas)+.2),
     ylab="Power",xlab=expression(Gamma),las=1,
     cex.lab=.8,main="I=1000",cex.main=.8)
legend(0,.85,c("Wilcoxon","Quade","U868","U878","U888","Mixed"),
       col=c("gray","gray","black","black","black","black"),
       lty=c(2,1,3,2,4,1),lwd=c(2,2,2,1,1,1),cex=.6,bty="n")
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="wilc")$power,
      col="gray",lty=2,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="quade")$power,
      col="gray",lty=1,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="u868")$power,
      col="black",lty=3,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="u878")$power,
      col="black",lty=2,lwd=1)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="u888")$power,
      col="black",lty=4,lwd=1)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="mixed")$power,
      col="black",lty=1,lwd=1)


plotPhi<-
  function (phi = "u868")
  {
    m<-NULL
    if (is.null(m)) {
      stopifnot(is.element(phi, c("u868", "u878", "u888", "u858",
                                  "quade", "wilc","mixed")))
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
    mixed<-function(pk){
      q<-multrnksU(pk, m1 = 19, m2 = 20, m = 20)+
        multrnksU(pk, m1 = 19, m2 = 19, m = 20)
      q/max(q)
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
      pk/max(pk)
    }
    wilc <- function(pk) {
      rep(1, length(pk))
    }

    if (is.null(m)) {
      if (phi == "mixed")
        phifunc <- mixed
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
    }
    u<-(0:300)/300
    phifunc(u)
  }




u<-(0:300)/300
plot(u,u,type="n",ylab=expression(varphi(v[i])),las=1,
     cex.axis=.6,xlab=expression(v[i]),cex.lab=.7,
     ylim=c(0,1.1),xlim=c(0,1),cex.main=.7)
legend(0,.90,c("Wilcoxon","Quade","U868","U878"),
       col=c("gray","gray","black","black"),
       lty=c(2,1,3,2),lwd=c(2,2,2,1),cex=.5,bty="n")
lines(u,plotPhi(phi="wilc"),
      col="gray",lty=2,lwd=2)
lines(u,plotPhi(phi="quade"),
      col="gray",lty=1,lwd=2)
lines(u,plotPhi(phi="u868"),
      col="black",lty=3,lwd=2)
lines(u,plotPhi(phi="u878"),
      col="black",lty=2,lwd=1)

u<-(0:300)/300
plot(u,u,type="n",ylab=expression(varphi(v[i])),las=1,
     cex.axis=.6,xlab=expression(v[i]),cex.lab=.7,
     ylim=c(0,1.1),xlim=c(.5,1),cex.main=.7)
legend(.5,.95,c("U888","Mixed"),
       col=c("black","black"),
       lty=c(4,1),lwd=c(1,1),cex=.5,bty="n")
lines(u,plotPhi(phi="u888"),
      col="black",lty=4,lwd=1)
lines(u,plotPhi(phi="mixed"),
      col="black",lty=1,lwd=1)

# Larger sample sizes for the appendix

gammas<-(10:110)/10
ssratio<-2000/406


plot(gammas,c(0,rep(1,length(gammas)-1)),type="n",
     cex.axis=.7,xlim=c(min(gammas)-.5,max(gammas)+.2),
     ylab="Power",xlab=expression(Gamma),las=1,
     cex.lab=.8,main="I=2000",cex.main=.8)
legend(0,.85,c("Wilcoxon","Quade","U868","U878","U888","Mixed"),
       col=c("gray","gray","black","black","black","black"),
       lty=c(2,1,3,2,4,1),lwd=c(2,2,2,1,1,1),cex=.6,bty="n")
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="wilc")$power,
      col="gray",lty=2,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="quade")$power,
      col="gray",lty=1,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="u868")$power,
      col="black",lty=3,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="u878")$power,
      col="black",lty=2,lwd=1)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="u888")$power,
      col="black",lty=4,lwd=1)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="mixed")$power,
      col="black",lty=1,lwd=1)

gammas<-(10:110)/10
ssratio<-5000/406


plot(gammas,c(0,rep(1,length(gammas)-1)),type="n",
     cex.axis=.7,xlim=c(min(gammas)-.5,max(gammas)+.2),
     ylab="Power",xlab=expression(Gamma),las=1,
     cex.lab=.8,main="I=5000",cex.main=.8)
legend(0,.85,c("Wilcoxon","Quade","U868","U878","U888","Mixed"),
       col=c("gray","gray","black","black","black","black"),
       lty=c(2,1,3,2,4,1),lwd=c(2,2,2,1,1,1),cex=.6,bty="n")
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="wilc")$power,
      col="gray",lty=2,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="quade")$power,
      col="gray",lty=1,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="u868")$power,
      col="black",lty=3,lwd=2)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="u878")$power,
      col="black",lty=2,lwd=1)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="u888")$power,
      col="black",lty=4,lwd=1)
lines(gammas,estPower(hdl0,gammas,ssratio=ssratio,phi="mixed")$power,
      col="black",lty=1,lwd=1)

rm(u,gammas,ssratio,plotPhi)
}
}

\keyword{ htest }
\concept{Sensitivity analysis}
\concept{Causal inference}
\concept{Observational study}
\concept{Quade's test}
\concept{Jackknife}
\concept{rbounds}
