\name{stepSolve}
\alias{stepSolve}
\title{
Root of a Monotone Decreasing Step Function
}
\description{
   Of limited interest to most users, stepSolve finds a root of a monotone decreasing step function by bisection search.  It is called by other functions in this package,
including wgtRankCI.  If the function f passes but never equals zero, then the argument where f passes zero is returned.  If the function equals zero for an interval, the average of the endpoints of the interval are returned.  The function is used to obtain confidence intervals and point estimates by inverting a rank test, in the spirit of the Hodges-Lehmman estimate.
}
\usage{
stepSolve(f, int, eps = 1e-05)
}
\arguments{
  \item{f}{
 A nonincreasing function of a single argument.
}
  \item{int}{
A vector of length 2 with int[1]<int[2].  A search for a root of f will be done on the interval int.
}
  \item{eps}{
  A small number.  In searching for a root, a difference of eps will be regarded as negligible.  A larger eps will produce an answer more quickly, while a smaller eps will produce a more accurate answer.
}
}

\value{
\item{average }{The average of low and high, defined below.}
\item{low}{See high, defined below.}
\item{high}{If the monotone decreasing step function f equals zero for an interval,
then that interval is approximately [low, high].  If f crosses zero but never
equals zero, then the crossing point is between low and high, which should be
close.  Reducing eps will search longer to improve the accuracy of [low, high].}
}


\references{
Hodges Jr, J. L. and Lehmann, E. L. (1963) <doi:10.1214/aoms/1177704172> Estimates of location Bbased on rank tests. The Annals of Mathematical Statistics, 34(2), 598-611.

Lehmann, E. L. (1963). Nonparametric confidence intervals for a shift parameter. The Annals of Mathematical Statistics, 34(4), 1507-1512.

Rosenbaum, P. R. (1993) <doi:10.1080/01621459.1993.10476405> Hodges-Lehmann point estimates of treatment effect in observational studies. Journal of the American Statistical Association, 88 (424), 1250-1253.
}

\author{
Paul R. Rosenbaum
}



\examples{
#  Uses stepSolve to find the one-sample Hodges-Lehmann estimate
\dontrun{
set.seed(1)
y<-rnorm(4)+7
wilcox.test(y,conf.int=TRUE,conf.level=.5)$estimate
f<-function(hl){wilcox.test(y-hl)$statistic-(length(y)*(1+length(y)))/4}
stepSolve(f,c(-100,100))
}
}
\keyword{ htest }
\concept{Optimization}
\concept{Binary search}

