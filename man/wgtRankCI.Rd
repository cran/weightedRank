\name{wgtRankCI}
\alias{wgtRankCI}

\title{
   Sensitivity Analysis for Confidence Intervals and Point Estimates from Weighted Rank Statistics in Block Designs
}
\description{
  Uses a weighted rank statistic to perform a sensitivity analysis for an I x J observational block design in which each of I blocks contains one treated individual and J-1 controls.  Inverts the test to obtain an estimate and a confidence interval for an additive treatment effect, tau; see Rosenbaum (1993, 2002, 2007).
}
\usage{
wgtRankCI(y, phi = "u868", phifunc = NULL, gamma = 1,
     alternative="greater", alpha=0.05, eps = 0.00001)
}
\arguments{
  \item{y}{
 A matrix or data frame with I rows and J columns.  Column 1 contains the response of the treated individuals and columns 2 throught J contain the responses of controls in the same block.  An error will result if y contains NAs.
}
  \item{phi}{
     The weight function to be applied to the ranks of the within block ranges.  The options are: (i) "wilc" for the stratified Wilcoxon test, which gives every block the same weight, (ii) "quade" which ranks the within block ranges from 1 to I, and is closely related to Quade's (1979) statistic; see also Tardif (1987), (iii) "u868" based on Rosenbaum (2011), (iv) u878 based on Rosenbaum (2011).  Note that phi is ignored if phifunc is not NULL.
}
  \item{phifunc}{ If not NULL, a user specified weight function for the ranks of the within block ranges.  The function should map [0,1] into [0,1].  The function is applied to the ranks divided by the sample size.
}
  \item{gamma}{
   A single number greater than or equal to 1.  gamma is the sensitivity parameter.  Two individuals with the same observed covariates may differ in their odds of treatment by at most a factor of gamma; see Rosenbaum (2002, Chapter 4; 2017, Chapter 9).
}
  \item{alternative}{Must equal "greater", "less", or "twosided".  This determines
  whether the confidence interval is one-sided or two-sided.}
  \item{alpha}{Coverage of the confidence interval is 1-alpha.}
  \item{eps}{
  A small number.  In searching for a root, a difference of eps will be regarded as negligible.  A larger eps will produce an answer more quickly, while a smaller eps will produce a more accurate answer.
}
}
\details{
This test is developed and evaluated in Rosenbaum (2024), and it is inverted for
point estimates and confidence intervals in the usual way.  Understand the test first.

The two-sided interval is the intersection of two one-sided intervals, each with coverage 1-alpha/2; see Cox (1977, Section 4.2).
}
\value{
 \item{estimate }{The interval of point estimates, reducing to a single number when gamma=1.}
 \item{confidence }{The one-sided or two-sided confidence interval.}
}
\references{
Brown, B. M. (1981). <doi:10.1093/biomet/68.1.235> Symmetric quantile averages and related estimators. Biometrika, 68(1), 235-242.

Cox, D. R. (1977). The role of significance tests [with discussion and reply]. Scandinavian Journal of Statistics, 4, 49-70.

Gastwirth, J. L., Krieger, A. M., and Rosenbaum, P. R. (2000).
<doi:10.1111/1467-9868.00249> Asymptotic separability in sensitivity analysis.  Journal of the Royal Statistical Society B 2000, 62, 545-556.

Hodges, J. L. and Lehmann, E. L. (1963). <doi:10.1214/aoms/1177704172> Estimates of Location Based on Rank Tests. The Annals of Mathematical Statistics, 34(2), 598-611.

Lehmann, E. L. (1975). Nonparametrics: Statistical Methods Based on Ranks. San Francisco: Holden-Day.

Quade, D. (1979). <doi:10.2307/2286991> Using weighted rankings in the analysis of complete blocks with additive block effects. Journal of the American Statistical Association, 74, 680-683.

Rosenbaum, P. R. (1993). <doi:10.2307/2291264> Hodges-Lehmann point estimates of treatment effect in observational studies. Journal of the American Statistical Association, 88(424), 1250-1253.

Rosenbaum, P. R. (2002). <doi:10.1007/978-1-4757-3692-2> Observational Studies (2nd Edition).  New York: Springer.

Rosenbaum, P. R. (2007). <doi:10.1111/j.1541-0420.2006.00717.x> Sensitivity analysis for M‐estimates, tests, and confidence intervals in matched observational studies. Biometrics, 63(2), 456-464.

Rosenbaum, P. R. (2011). <doi:10.1111/j.1541-0420.2010.01535.x> A new U‐Statistic with superior design sensitivity in matched observational studies. Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2014) <doi:10.1080/01621459.2013.879261> Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls.  Journal of the American Statistical Association, 109(507), 1145-1158.

Rosenbaum, P. R. (2015). <doi:10.1080/01621459.2014.960968> Bahadur efficiency of sensitivity analyses in observational studies. Journal of the American Statistical Association, 110(509), 205-217.

Rosenbaum, P. R. (2017). <doi:10.4159/9780674982697> Observation and Experiment: An Introduction to Causal Inference.  Cambridge, MA: Harvard University Press.

Rosenbaum, P. R. (2018). <doi:10.1214/18-AOAS1153> Sensitivity analysis for stratified comparisons in an observational study of the effect of smoking on homocysteine levels. The Annals of Applied Statistics, 12(4), 2312-2334.

Rosenbaum, P. R. (2024) <doi:10.1080/01621459.2023.2221402> Bahadur efficiency of observational block designs. Journal of the American Statistical Association.

Tardif, S. (1987). <doi:10.2307/2289476> Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82(398), 637-644.
}
\author{
Paul R. Rosenbaum
}
\note{
The computations use the separable approximation discussed in Gastwirth et al. (2000) and Rosenbaum (2018).  Compare with the method in Rosenbaum (2014) and the R package sensitivitymw.
}

\seealso{
An alternative approach avoids rank tests and uses weighted M-statistics instead, as in the sensitivitymw package and Rosenbaum (2014).  However, Bahadur efficiency calculations are available for weighted rank statistics; see Rosenbaum (2024).
}
\examples{
\dontrun{
data(aHDL)
y<-t(matrix(aHDL$hdl,4,406))

# Without sensitivity analysis
wgtRankCI(y,phi="wilc",gamma=1) # Stratified Wilcoxon rank sum test

# With sensitivity analysis

# The CI excludes 0 if the test rejects 0.
wgtRank(y,phi="wilc",gamma=3.7) # Stratified Wilcoxon rank sum test
wgtRankCI(y,phi="wilc",gamma=3.7)

wgtRank(y,phi="u878",gamma=3.7) # U-statistic with weights (8,7,8)
wgtRankCI(y,phi="u878",gamma=3.7)

# A user defined weight function, brown, analogous to Brown (1981).
brown<-function(v){((v>=.333)+(v>=.667))/2}
wgtRankCI(y,phifunc=brown,gamma=3.7)


#
# COMPARISONS WITH STANDARD METHODS
#
y2<-c(-17,  -3,  -1,   2,   8,   9,  10,  11,  12,  16,  20)
y2a<-cbind(y2,-y2)/2
#
# Compare with the usual Hodges-Lehmann estimate,
# namely the median of the Walsh averages (yi+yj)/2
#
hl<-function(y){
  w<-NULL
  I<-length(y)
  for (i in 1:I) for (j in i:I) w<-c(w,(y[i]+y[j])/2)
  median(w)
}
hl(y2)
wgtRankCI(y2a,phi="quade")$estimate
#
# Compare with wilcox.test confidence interval in the stats package
#
stats::wilcox.test(y2,conf.int=TRUE,correct=FALSE,exact=FALSE,
        alternative = "greater")$conf.int
wgtRankCI(y2a,phi="quade",alternative="greater")$confidence
stats::wilcox.test(y2,conf.int=TRUE,correct=FALSE,exact=FALSE,
        alternative = "less")$conf.int
wgtRankCI(y2a,phi="quade",alternative="less")$confidence
stats::wilcox.test(y2,conf.int=TRUE,correct=FALSE,exact=FALSE,
        alternative = "two.sided")$conf.int
wgtRankCI(y2a,phi="quade",alternative="twosided")$confidence
#
# Compare with senWilcox in the DOS2 package
#
# Note: senWilcox reports only one-sided estimates with 1-sided tests
#
DOS2::senWilcox(y2,conf.int=TRUE,alternative="greater",gamma=1.25)
wgtRankCI(y2a,phi="quade",gamma=1.25,alternative = "greater")
DOS2::senWilcox(y2,conf.int=TRUE,alternative="less",gamma=1.25)
wgtRankCI(y2a,phi="quade",gamma=1.25,alternative = "less")
DOS2::senWilcox(y2,conf.int=TRUE,alternative="twosided",gamma=1.25)
wgtRankCI(y2a,phi="quade",gamma=1.25,alternative = "twosided")
#
# Compare with senU in the DOS2 package
#
# Note: senU reports only one-sided estimates with 1-sided tests
#
DOS2::senU(y2,gamma=1.25,m=8,m1=6,m2=8,alternative = "greater",
          exact=FALSE,conf.int = TRUE)
wgtRankCI(y2a,phi="u868",gamma=1.25,alternative = "greater")
DOS2::senU(y2,gamma=1.25,m=8,m1=6,m2=8,alternative = "less",
          exact=FALSE,conf.int = TRUE)
wgtRankCI(y2a,phi="u868",gamma=1.25,alternative = "less")
DOS2::senU(y2,gamma=1.25,m=8,m1=6,m2=8,alternative = "twosided",
          exact=FALSE,conf.int = TRUE)
wgtRankCI(y2a,phi="u868",gamma=1.25,alternative = "twosided")
}
}

\keyword{ htest }
\concept{Sensitivity analysis for confidence intervals}
\concept{Sensitivity analysis for point estimates}
\concept{Causal inference}
\concept{Observational study}
\concept{Quade's test}
\concept{rbounds}
