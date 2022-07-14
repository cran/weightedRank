\name{weightedRank-package}
\alias{weightedRank-package}
\alias{weightedRank}
\docType{package}
\title{
\packageTitle{weightedRank}
}
\description{
\packageDescription{weightedRank}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{weightedRank}
\packageIndices{weightedRank}
The package conducts either fixed or adaptive sensitivity analyses for observational studies with I blocks and J individuals in each block, one treated and J-1 controls.  The two main functions are wgtRank() for a fixed test statistic, and wgtRanktt() for an adaptive choice of one of two test statistics.
}
\author{
\packageAuthor{weightedRank}

Maintainer: \packageMaintainer{weightedRank}
}
\references{
Berk, R. H. and Jones, D. H. (1978) <https://www.jstor.org/stable/4615706> Relatively optimal combinations of test statistics. Scandinavian Journal of Statistics, 5, 158-162.

Quade, D. (1979) <doi:10.2307/2286991> Using weighted rankings in the analysis of complete blocks with additive block effects. Journal of the American Statistical Association, 74, 680-683.

Rosenbaum, P. R. (2011) <doi:10.1111/j.1541-0420.2010.01535.x> A new U‐Statistic with superior design sensitivity in matched observational studies. Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2012) <doi:10.1093/biomet/ass032> Testing one hypothesis twice in observational studies. Biometrika, 99(4), 763-774.

Rosenbaum, P. R. (2022)  Bahadur efficiency of observational block designs.  Manuscript.

Tardif, S. (1987) <doi:10.2307/2289476> Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82(398), 637-644.
}
\keyword{ package }

\examples{
data(aHDL)
y<-t(matrix(aHDL$hdl,4,406))
wgtRank(y,phi="u878",gamma=6) # New U-statistic weights (8,7,8)
wgtRanktt(y,phi1="u868",phi2="u878",gamma=5.9)
}
\concept{Sensitivity analysis}
\concept{Causal inference}
\concept{Observational study}
\concept{Adaptive inference}
\concept{Testing twice}
\concept{Block design}
