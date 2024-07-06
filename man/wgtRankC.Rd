\name{wgtRankC}
\alias{wgtRankC}
\title{
Sensitivity Analysis for a Conditional Weighted Rank Test
}
\description{
  Uses a conditional weighted rank statistic to perform a sensitivity analysis for an I x J observational block design in which each of I blocks contains one treated individual and J-1 controls.  The size J of a block must be at least 3.
}
\usage{
wgtRankC(y, phi = "u878", phifunc = NULL, gamma = 1, alternative = "greater")
}
\arguments{
  \item{y}{
 A matrix or data frame with I rows and J columns, where J is at least 3.  Column 1 contains the response of the treated individuals and columns 2 throught J contain the responses of controls in the same block.  An error will result if y contains NAs.
 An error will result if y has fewer than 3 columns.
}
  \item{phi}{
     The weight function to be applied to the ranks of the within block ranges.  The options are: (i) "wilc" for the stratified Wilcoxon test, which gives every block the same weight, (ii) "quade" which ranks the within block ranges from 1 to I, and is closely related to Quade's (1979) statistic; see also Tardif (1987), (iii) "u868", "u878", or "u888" based on Rosenbaum (2011).  Note that phi is ignored if phifunc is not NULL.
}
  \item{phifunc}{ If not NULL, a user specified weight function for the ranks of the within block ranges.  The function should map [0,1] into [0,1].  The function is applied to the ranks divided by the sample size.  See the example.
}
  \item{gamma}{
   A single number greater than or equal to 1.  gamma is the sensitivity parameter.  Two individuals with the same observed covariates may differ in their odds of treatment by at most a factor of gamma; see Rosenbaum (1987; 2017, Chapter 9).
}
  \item{alternative}{
  The null hypothesis asserts that the treatment has no effect.  If
alternative is "greater", then the null hypothesis is tested against
the alternative that the treatment increases the response of treated
individuals.  If alternative is "less", then the null hypothesis is tested against
the alternative that the treatment decreases the response of treated
individuals.
}
}
\details{
The conditional test restricts attention to blocks in which the treated
individual has either the highest or lowest response in a block.  This
tactic may be shown to increase design sensitivity; see Rosenbaum (2024).
}
\value{
\item{pval }{The upper bound on the one-sided P-value}
\item{detail }{Details of the computation of the P-value}
\item{block.information }{The conditional test uses a block only if
the treated individual has either the highest or lowest
respones in the block.  This vector indicates: (i) the
number of blocks that were used, (ii) in those blocks, the number
in which the treated individual had the highest response, and (iii)
in those block, the number of blocks in which there was a tie
for the maximum or minimum response.}
}
\references{
Brown, B. M. (1981) <doi:10.1093/biomet/68.1.235> Symmetric quantile averages and related estimators. Biometrika, 68, 235-242.

Cox, D. R. (1977). The role of significance tests (with discussion and reply). Scandinavian Journal of Statistics, 4, 49-70.

Noether, G. E. (1963) <doi:10.2307/2283321> Efficiency of the Wilcoxon two-sample statistic for randomized blocks. Journal of the American Statistical Association, 58, 894-898.

Noether, G. E. (1973) <doi:10.2307/2284805> Some simple distribution-free confidence intervals for the center of a symmetric distribution. Journal of the American Statistical Association, 68, 716-719.

Quade, D. (1979) <doi:10.2307/2286991> Using weighted rankings in the analysis of complete blocks with additive block effects. Journal of the American Statistical Association, 74, 680-683.

Rosenbaum, P. R. (1987). <doi:10.2307/2336017> Sensitivity analysis for certain permutation inferences in matched observational studies. Biometrika, 74(1), 13-26.

Rosenbaum, P. R. (2004) <doi:10.1093/biomet/91.1.153> Design sensitivity in observational studies. Biometrika, 91, 153-164.

Rosenbaum, P. R. (2011). <doi:10.1111/j.1541-0420.2010.01535.x> A new U‐Statistic with superior design sensitivity in matched observational studies. Biometrics, 67(3), 1017-1027.

Rosenbaum, P. (2017). <doi:10.4159/9780674982697> Observation and Experiment: An Introduction to Causal Inference.  Cambridge, MA: Harvard University Press.

Rosenbaum, P. R. (2024) <doi:10.1080/01621459.2023.2221402> Bahadur efficiency of observational block designs. Journal of the American Statistical Association.

Rosenbaum, P. R. (2024) A conditioning tactic that increases design
sensitivity in observational block designs.  Manuscript.

Tardif, S. (1987) <doi:10.2307/2289476> Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82, 637-644.
}
\author{
Paul R. Rosenbaum
}
\note{
Do a two-sided test by testing in both tails
and rejecting at level alpha if the smaller P-value is less than alpha/2; see
Cox (1977, section 4.2).
Aside from labeling the output, setting alternative to "less" is the same
as applying the test to -y with alternative set to "greater".
}

\examples{
data(aHDLe)
y<-t(matrix(aHDLe$hdl,4,722))
colnames(y)<-c("D","R","N","B")
y<-y[,c(1,3,2,4)]
boxplot(y,ylab="HDL Cholesterol",las=1,xlab="Group")
wgtRankC(y,gamma=6.65,phi="u878")
wgtRankC(y,gamma=6.24,phi="u868")
wgtRankC(y,phi="quade",gamma=5.28)
wgtRankC(y,phi="wilc",gamma=4.9)
#
#  Examples with a user-defined phi-function
#  Rank scores of Brown (1981) and Noether (1973)
brown<-function(v){(v>=(1/3))+(v>=(2/3))}
wgtRankC(y,phifunc=brown,gamma=5.5)
noether<-function(v){v>=(2/3)}
wgtRankC(y,phifunc=noether,gamma=6.5)
}
\keyword{ htest }
\concept{Sensitivity analysis}
\concept{Causal inference}
\concept{Observational study}
\concept{Quade's test}
\concept{rbounds}
\concept{Conditional weighted rank test}
