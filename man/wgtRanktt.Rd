\name{wgtRanktt}
\alias{wgtRanktt}
\title{
Adaptive Inference Using Two Test Statistics in a Block Design
}
\description{
   Tests twice, using the better of two test statistics; see Rosenbaum (2012, 2022).
}
\usage{
wgtRanktt(y, phi1 = "u868", phi2 = "u878", phifunc1 = NULL, phifunc2 = NULL, gamma = 1)
}
\arguments{
  \item{y}{
 A matrix or data frame with I rows and J columns.  Column 1 contains the response of the treated individuals and columns 2 throught J contain the responses of controls in the same block.  A error will result if y contains NAs.
}
  \item{phi1}{
The weight function to be applied to the ranks of the within block ranges.  The options are: (i) "wilc" for the stratified Wilcoxon test, which gives every block the same weight, (ii) "quade" which ranks the within block ranges from 1 to I, and is closely related to Quade's (1979) statistic; see also Tardif (1987), (iii) "u868" based on Rosenbaum (2011), (iv) u878 based on Rosenbaum (2011).  Note that phi is ignored if phifunc is not NULL.
}
  \item{phi2}{
  See phi1.
}
  \item{phifunc1}{
If not NULL, a user specified weight function for the ranks of the within block rates.  The function should map [0,1] into [0,1].  The function is applied to the ranks divided by the sample size.  See the example.
}
  \item{phifunc2}{
 See phifunc1.
}
  \item{gamma}{
A single number greater than or equal to 1.  gamma is the sensitivity parameter.  Two individuals with the same observed covariates may differ in their odds of treatment by at most a factor of gamma; see Rosenbaum (1987; 2017, Chapter 9).
}
}

\value{
\item{jointP }{Upper bound on the one-sided joint P-value obtained
from two test statistics in the presence of a bias of at most gamma.}
\item{cor12 }{Correlation of the two test statistics at the treatment assignment distribution that provides the joint upper bound.  Often, this correlation is high, so the joint distribution that is used here is much less conservative than use of the Bonferroni inequality when testing twice.}
\item{detail }{Details about the two statistics separately.  Equivalent to the result from wgtRank() run twice with different test statistics.}
}

\references{
Berk, R. H. and Jones, D. H. (1978)
<https://www.jstor.org/stable/4615706> Relatively optimal combinations of test statistics. Scandinavian Journal of Statistics, 5, 158-162.

Brown, B. M. (1981) <doi:10.1093/biomet/68.1.235> Symmetric quantile averages and related estimators. Biometrika, 68(1), 235-242.

Gastwirth, J. L., Krieger, A. M., and Rosenbaum, P. R. (2000)
<doi:10.1111/1467-9868.00249> Asymptotic separability in sensitivity analysis.  Journal of the Royal Statistical Society B 2000, 62, 545-556.

Lehmann, E. L. (1975)  Nonparametrics: Statistical Methods Based on Ranks. San Francisco: Holden-Day.

Quade, D. (1979) <doi:10.2307/2286991> Using weighted rankings in the analysis of complete blocks with additive block effects. Journal of the American Statistical Association, 74, 680-683.

Rosenbaum, P. R. (1987) <doi:10.2307/2336017> Sensitivity analysis for certain permutation inferences in matched observational studies. Biometrika, 74(1), 13-26.

Rosenbaum, P. R. (2011) <doi:10.1111/j.1541-0420.2010.01535.x> A new U‐Statistic with superior design sensitivity in matched observational studies. Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2012) <doi:10.1093/biomet/ass032> Testing one hypothesis twice in observational studies. Biometrika, 99(4), 763-774.

Rosenbaum, P. R. (2014) <doi:10.1080/01621459.2013.879261> Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls.  Journal of the American Statistical Association, 109(507), 1145-1158.

Rosenbaum, P. R. (2015) <doi:10.1080/01621459.2014.960968> Bahadur efficiency of sensitivity analyses in observational studies. Journal of the American Statistical Association, 110(509), 205-217.

Rosenbaum, P. (2017) <doi:10.4159/9780674982697> Observation and Experiment: An Introduction to Causal Inference.  Cambridge, MA: Harvard University Press.

Rosenbaum, P. R. (2018) <doi:10.1214/18-AOAS1153> Sensitivity analysis for stratified comparisons in an observational study of the effect of smoking on homocysteine
levels. The Annals of Applied Statistics, 12(4), 2312-2334.

Rosenbaum, P. R. (2022)  Bahadur efficiency of observational block designs.  Manuscript.

Tardif, S. (1987) <doi:10.2307/2289476> Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82(398), 637-644.
}
\author{
Paul R. Rosenbaum
}
\note{
For discussion of testing twice in matched pairs, see Rosenbaum (2012).

Testing twice is also possible in block designs using weighted rank statistics because the same value of the unobserved covariate provides the upper bound for both statistics when using the separable approximation in Gastwirth et al. (2000) and Rosenbaum (2018, Remarks 4 and 5).  See also Rosenbaum (2022) where the Bahadur efficiency of such tests is computed.

Other packages that use testing twice in a different way are "sensitivity2x2xk"" and "testtwice".  The "testtwice" package is restricted to matched pairs, and "sensitivity2x2xk"" is for binary outcomes.  With some attention to detail (e.g., the handling of zero pair differences), in the case of matched pairs, the "testtwice" package and the wgtRanktt() function will yield identical results.  In that sense, wgtRanktt() extends the method to blocks designs.

Testing twice achieves the larger Bahadur efficiency of the two component statistics; see Berk and Jones (1978).
}


\examples{
data(aHDL)
y<-t(matrix(aHDL$hdl,4,406))

# This is the simplest example of a general property.  The
# example simply illustrates, but does not fully exploit
# the property.  In this case, use of the stratified
# Wilcoxon statistic is a mistake, because Quade's
# statistic correctly reports insensitivity to a bias
# of gamma=4.5, but the stratified Wilcoxon statistic
# is sensitive at gamma=3.5.  The adaptive procedure
# that does both tests and corrects for multiple testing
# is insensitive to gamma=4.4; so, it is almost as good
# as knowing what you cannot know, namely that Quade's
# statistic is the better choice in this one example.
# The price paid for testing twice is very small;
# see Berk and Jones (1978) and Rosenbaum (2012, 2022).
wgtRank(y,phi="wilc",gamma=3.5)
wgtRank(y,phi="quade",gamma=3.5)
wgtRank(y,phi="wilc",gamma=4.5)
wgtRank(y,phi="quade",gamma=4.5)
wgtRanktt(y,phi1="wilc",phi2="quade",gamma=4.4)

# Sensitivity to gamma=3.5 is very different from
# sensitivity to gamma=4.4; see documentation for amplify.
amplify(3.5,8)
amplify(4.4,8)


# In this example, u878 exhibits greater insensitivity to bias
# than u868.  However, adaptive inference using both is almost
# as good as the better statistic, yet it strongly controls the
# family-wise error rate despite testing twice;
# see Rosenbaum (2012,2022).
wgtRank(y,phi="u868",gamma=6) # New U-statistic weights (8,6,8)
wgtRank(y,phi="u878",gamma=6) # New U-statistic weights (8,7,8)
wgtRanktt(y,phi1="u868",phi2="u878",gamma=5.9)

# A user defined weight function, brown, analogous to Brown (1981).
brown<-function(v){((v>=.333)+(v>=.667))/2}
# In this example, the joint test rejects based on u878
wgtRanktt(y,phi1="u878",phifunc2=brown,gamma=5.8)
}
\keyword{ htest }
\concept{Sensitivity analysis}
\concept{Causal inference}
\concept{Observational study}
\concept{Adaptive inference}
\concept{Testing twice}
