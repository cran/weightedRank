\name{gwgtRank}
\alias{gwgtRank}

\title{
   Generalized Sensitivity Analysis for Weighted Rank Statistics in Block Designs
}
\description{
  Uses a weighted rank statistic to perform a sensitivity analysis for an I x J observational block design in which each of I blocks contain J individuals, some of whom are treated individuals and others are controls.
}
\usage{
gwgtRank(y, z, phi = "u868", phifunc = NULL, gamma = 1, detail=FALSE)
}
\arguments{
  \item{y}{
 A matrix or data frame of responses with I rows and J columns.  An error will result if y contains NAs.
}
\item{z}{
 A matrix or data frame of treatment indicators with I rows and J columns.  In z, the ith row and jth column is a 1 if this individual is treated, or a 0 if this individual is a control.  If a row of z contains J ones or J zeros, that block will be removed from further computations as there is no one to compare, and a warning will appear.  An error will result if z contains NAs.
}
  \item{phi}{
     The weight function to be applied to the ranks of the within block ranges.  The options are: (i) "wilc" for the stratified Wilcoxon test, which gives every block the same weight, (ii) "quade" which ranks the within block ranges from 1 to I, and is closely related to Quade's (1979) statistic; see also Tardif (1987), (iii) "u868" based on Rosenbaum (2011), (iv) u878 based on Rosenbaum (2011).  Note that phi is ignored if phifunc is not NULL.
}
  \item{phifunc}{ If not NULL, a user specified weight function for the ranks of the within block ranges.  The function should map [0,1] into [0,1].  The function is applied to the ranks divided by the sample size.  See the example.
}
  \item{gamma}{
   A single number greater than or equal to 1.  gamma is the sensitivity parameter.  Two individuals with the same observed covariates may differ in their odds of treatment by at most a factor of gamma; see Rosenbaum (1987; 2017, Chapter 9).
}
  \item{detail}{
   If detail=FALSE, then brief output is given, similar to that provided by wgtRank().  In both cases, the asymptotically separable approximation of Gastwirth et al. (2000) is used. If detail=TRUE, then the separable approximation is compared with a bound in Rosenbaum (2018), and the output from the senstrat package is given.  For moderate J and large I, the bound is likely to agree exactly with the separable approximation.
}
}
\details{
The function gwgtRank()
differs from wgtRank() in that wgtRank requires one treated individual and J-1 controls,
whereas gwgtRank permits 1 to J-1 treated indivuals with the rest as controls.

This function uses the senstrat() function in the senstrat package to perform the calculations.

To test in the lower tail -- to test against the alternative that treated responses are lower than control responses, apply the function to -y.  For a two-sided test, do both one-sided tests and apply the Bonferroni inequality, doubling the smaller of the two one-sided P-value bounds; see Cox (1977, Section 4.2).
}
\value{
 \item{pval }{Upper bound on the one-sided P-value when testing
 the null hypothesis of no treatment effect against the alternative hypothesis that treated responses are higher than control responses.}
 \item{detail }{Details of the computation of pval: the standardized deviate, the test statistic, its null expectation, its null variance and the value of gamma.}
}
\references{
Brown, B. M. (1981). <doi:10.1093/biomet/68.1.235> Symmetric quantile averages and related estimators. Biometrika, 68(1), 235-242.

Cox, D. R. (1977). The role of significance tests [with discussion and reply]. Scandinavian Journal of Statistics, 4, 49-70.

Gastwirth, J. L., Krieger, A. M., and Rosenbaum, P. R. (2000).
<doi:10.1111/1467-9868.00249> Asymptotic separability in sensitivity analysis.  Journal of the Royal Statistical Society B 2000, 62, 545-556.

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

Rosenbaum, P. R. (2024) <doi:10.1080/01621459.2023.2221402> Bahadur efficiency of observational block designs. Journal of the American Statistical Association.

Tardif, S. (1987). <doi:10.2307/2289476> Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82(398), 637-644.
}
\author{
Paul R. Rosenbaum
}

\seealso{
An alternative approach avoids rank tests and uses weighted M-statistics instead, as in the sensitivitymw package and Rosenbaum (2014).  However, Bahadur efficiency calculations are available for weighted rank statistics; see Rosenbaum (2024).
}
\examples{
data(aHDL)
y<-t(matrix(aHDL$hdl,4,406))
z<-matrix(0,dim(y)[1],dim(y)[2])
z[,1]<-1
gwgtRank(y,z,phi="quade",gamma=3.5) # Quade's test
gwgtRank(y,z,phi="quade",gamma=3.5,detail=TRUE) # Alternative output

# A user defined weight function, brown, analogous to Brown (1981).
brown<-function(v){((v>=.333)+(v>=.667))/2}
gwgtRank(y,z,phifunc=brown,gamma=4.7)
}

\keyword{ htest }
\concept{Sensitivity analysis}
\concept{Causal inference}
\concept{Observational study}
\concept{Quade's test}
\concept{rbounds}
