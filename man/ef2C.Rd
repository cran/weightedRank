\name{ef2C}
\alias{ef2C}

\title{
Evidence Factors For Matched Triples With Two Control Groups
}
\description{
In an observational complete block design, with bocks of size three, each containing
a treated individual and one control from each of two control groups, ef2C performs
the evidence factor analysis suggested in Rosenbaum (2023).  One factor compares
the treated group to the first control group in a matched pairs analysis.  The
other factor pools the treated group and the first control group and compares it
to the second control group.
}
\usage{
ef2C(y, gamma = 1, upsilon = 1, alternative = "greater",
     trunc = 0.2)
}

\arguments{
  \item{y}{
With I blocks and 3 individuals in each block, y is and I x 3 matrix or dataframe containing the outcomes.  The first column is the response of the treated individual.  The second response is the response of the control from the first control group.  The third response is the response of the control from the second control group.
}
  \item{gamma}{
A real number >=1 giving the value of the sensitivity parameter for the comparison of the treated group and the first control group.  gamma=1 yields a randomization test.
}
  \item{upsilon}{
A real number >=1 giving the value of the sensitivity parameter for the comparison of the second control group and the combination of the treated group plus the first control group.  upsilon=1 yields a randomization test.
}
  \item{alternative}{
Use alternative=greater if the treatment is expect to cause an increase in the response in y.
Use alternative=less if the treatment is expect to cause an decrease in the response in y.  In this context, a two-sided test is best viewed as two one-sided tests with a Bonferroni correction, e.g., testing in both tails at level 0.025 to ensure overall level of 0.05; see Cox (1977).  For more information, see the notes.
}
  \item{trunc}{
The two P-values from the two factors are combined using the trucated product of P-values due to Zaykin et al. (2002): it is the P-value derived from the product of those P-values that are less
than trunc.  For more information, see the notes.
}
}
\details{
This is the recommended analysis in Rosenbaum (2023).  The example below reproduces
some results from the example in that paper.  That paper considered 40 test statistics in terms of
the Bahadur efficiency of a sensitivity -- all of these analyses can be
reproduced by the more flexible but more complicated dwgtRank function.

The comparison of the treated group and the first control group is equivalent to
dwgtRank(y[,1:2],gamma=gamma,m=8,m1=7,m2=8,range=TRUE,alternative="greater"), and
these settings are motivated by results in Rosenbaum (2011, 2015).  Notice that
y[,1:2] uses the first two columns of y.

The comparison of the second control group and the merger of the treated group
with the first control group is equivalent to
dwgtRank(y[,3:1], gamma=upsilon, m=8,m1=8,m2=8, range=FALSE, alternative="less", scores=c(1,2,5)), and  these settings are motivated by results in Rosenbaum (2023).  Notice that
y[,3:1] compares the third column to the pooled group consisting of columns 1 and 2.
}
\value{
\item{pvals }{Upper bounds on the one-sided P-values for the two factors and their combination.}
\item{detail }{A matrix with some details of the computations that produced the
P-values.}
}

\note{
The two P-values from the two factors are combined using the trucated product of P-values due to Zaykin et al. (2002): it is the P-value derived from the product of those P-values that are less
than trunc.  Taking trunc=1 yields Fisher's method for combining independent P-values.  Fisher's method is not ideal when combining P-value bounds produced by sensitivity analyses; see Hsu et al. (2013).  Reasonable values are trunc=.1, truc=.15 and trunc=.2.  As illustrated in the example
below, lower truncation values produce smaller combined P-values when the P-values are
below the truncation point, but a P-value that barely exceeds the truncation point
is effectively discarded.  Hsu et al. (2013) compare truncation values when used in a sensitivity analysis.  For discussion of combining sensitivity analyses as
independent, see the required conditions in Rosenbaum (2011b, 2021).  These conditions
hold for the comparison performed by ef2C.
}

\note{
The setting alternative = "less" simply replaces y by -y before testing in the
upper tail.
}

\note{
For a deeper understanding, see the documentation of dwgtRank.  That function
is more general, but it requires more attention to detail by the user.  The
documentation for dwgtRank also produces additional analyses from Rosenbaum (2023).
}

\references{

Cox, D. R. (1977). The role of significance tests [with discussion and reply]. Scandinavian Journal of Statistics, 4, 49-70.

Hsu, J. Y., Small, D. S., & Rosenbaum, P. R. (2013). <doi:10.1080/01621459.2012.742018> Effect modification and design sensitivity in observational studies. Journal of the American Statistical Association, 108(501), 135-148.

Rosenbaum, P. R. (1987). <doi:10.1214/ss/1177013232> The role of a second control group in an observational study. Statistical Science, 2, 292-306.

Rosenbaum, P. R. (2011a). <doi:10.1111/j.1541-0420.2010.01535.x> A new U‐Statistic with superior design sensitivity in matched observational studies. Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2011). <doi:10.1198/jasa.2011.tm10422> Some approximate evidence factors in observational studies. Journal of the American Statistical Association, 106(493), 285-295.

Rosenbaum, P. R. (2015). <doi:10.1080/01621459.2014.960968> Bahadur efficiency of sensitivity analyses in observational studies. Journal of the American Statistical Association, 110(509), 205-217.

Rosenbaum, P. R. (2021) <doi:10.1201/9781003039648> Replication and Evidence Factors in Observational Studies.  Chapman and Hall/CRC.

Rosenbaum, P. R. (2023) A second evidence factor for a second control group.  Manuscript.  This paper develops a new statistic that performs well when there are two control groups and the second control group forms a second evidence factor.  See also Rosenbaum (2011b, 2021).

Zaykin, D. V., Zhivotovsky, L. A., Westfall, P. H. and Weir, B. S. (2002) <doi:10.1002/gepi.0042> Truncated product method of combining P-values. Genetic Epidemiology, 22, 170-185.

}
\author{
Paul R. Rosenbaum
}

\examples{
# The calculation below reproduce analyses from Rosenbaum (2023).

data(aBP)
attach(aBP)

yD<-t(matrix(bpDiastolic,3,207))
yS<-t(matrix(bpSystolic,3,207))
vS<-c(yS[,1]-yS[,2],yS[,1]-yS[,3],yS[,2]-yS[,3])
vD<-c(yD[,1]-yD[,2],yD[,1]-yD[,3],yD[,2]-yD[,3])
y<-(yD/median(abs(vD)))+(yS/median(abs(vS)))

# EVIDENCE FACTOR ANALYSIS, COMBINING TWO FACTORS
# The evidence factor analysis compares treated to the first control group,
# then compares the second control group to the pooled group consisting of
# treated and first control, then combines the two analyses using meta-analysis.
# Treated/first-control matched pairs are compared using the method in
# Rosenbaum (2011).

ef2C(y,gamma=2.3,upsilon=1.45)
amplify(2.3,4)
amplify(1.45,2.5)

# THE COMBINED ANALYSIS IS INSENSITIVE TO LARGER BIASES
# The combined analysis is insensitive to larger biases
# than are its components
ef2C(y,gamma=2.6,upsilon=1.7)
amplify(2.6, 5)
amplify(1.7,c(2.7,3))

# The calculations above are also produced in the
# example for dwgtRank, where alternative
# analyses from Rosenbaum (2023) are compared.

####################################################

# Comparing trucation points to understand trunc:
ef2C(y,gamma=2.6,upsilon=1.7,trunc=.2) # Default
ef2C(y,gamma=2.6,upsilon=1.7,trunc=1) # Fisher's method
ef2C(y,gamma=2.6,upsilon=1.7,trunc=.1)
ef2C(y,gamma=2.5,upsilon=1.6,trunc=.2)
ef2C(y,gamma=2.5,upsilon=1.6,trunc=.1)
# See Hsu et al. (2013) for discussion of the
# truncation point for a sensitivity analysis.
}

\keyword{ nonparametric }
\keyword{ htest }
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Evidence factors}
\concept{Multiple control groups}
\concept{rbounds}
