\name{wgtRank}
\alias{wgtRank}

\title{
   Sensitivity Analysis for Weighted Rank Statistics in Block Designs
}
\description{
  Uses a weighted rank statistic to perform a sensitivity analysis for an I x J observational block design in which each of I blocks contains one treated individual and J-1 controls.
}
\usage{
wgtRank(y, phi = "u868", phifunc = NULL, gamma = 1)
}
\arguments{
  \item{y}{
 A matrix or data frame with I rows and J columns.  Column 1 contains the response of the treated individuals and columns 2 throught J contain the responses of controls in the same block.  A error will result if y contains NAs.
}
  \item{phi}{
     The weight function to be applied to the ranks of the within block ranges.  The options are: (i) "wilc" for the stratified Wilcoxon test, which gives every block the same weight, (ii) "quade" which ranks the within block ranges from 1 to I, and is closely related to Quade's (1979) statistic; see also Tardif (1987), (iii) "u868" based on Rosenbaum (2011), (iv) u878 based on Rosenbaum (2011).  Note that phi is ignored if phifunc is not NULL.
}
  \item{phifunc}{ If not NULL, a user specified weight function for the ranks of the within block rates.  The function should map [0,1] into [0,1].  The function is applied to the ranks divided by the sample size.  See the example.
}
  \item{gamma}{
   A single number greater than or equal to 1.  gamma is the sensitivity parameter.  Two individuals with the same observed covariates may differ in their odds of treatment by at most a factor of gamma; see Rosenbaum (1987; 2017, Chapter 9).
}
}
\details{
This method is developed and evaluated in Rosenbaum (2022).

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

Rosenbaum, P. R. (2022).  Bahadur efficiency of observational block designs.  Manuscript.

Tardif, S. (1987). <doi:10.2307/2289476> Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82(398), 637-644.
}
\author{
Paul R. Rosenbaum
}
\note{
The computations use the separable approximation discussed in Gastwirth et al. (2000) and Rosenbaum (2018).  Compare with the method in Rosenbaum (2014) and the R package sensitivitymw.
}

\seealso{
An alternative approach avoids rank tests and uses weighted M-statistics instead, as in the sensitivitymw package and Rosenbaum (2014).  However, Bahadur efficiency calculations are available for weighted rank statistics; see Rosenbaum (2022).
}
\examples{
data(aHDL)
y<-t(matrix(aHDL$hdl,4,406))
wgtRank(y,phi="wilc",gamma=3.5) # Stratified Wilcoxon rank sum test
wgtRank(y,phi="quade",gamma=3.5) # Quade's test
wgtRank(y,phi="quade",gamma=4.5) # Quade's test, larger gamma
wgtRank(y,phi="quade",gamma=4.6) # Quade's test, larger gamma
wgtRank(y,phi="u868",gamma=5.4) # New U-statistic weights (8,6,8)
wgtRank(y,phi="u878",gamma=6) # New U-statistic weights (8,7,8)

# As an aid to interpreting gamma, see the amplify function.
amplify(3.5,8)
amplify(4.6,8)
amplify(5.4,8)
amplify(6,8)

# A user defined weight function, brown, analogous to Brown (1981).
brown<-function(v){((v>=.333)+(v>=.667))/2}
wgtRank(y,phifunc=brown,gamma=4.7)
}

\keyword{ htest }
\concept{Sensitivity analysis}
\concept{Causal inference}
\concept{Observational study}
\concept{Quade's test}
\concept{rbounds}
