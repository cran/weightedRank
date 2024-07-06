\name{dwgtRank}
\alias{dwgtRank}

\title{
Weighted Rank Statistics for Evidence Factors with Two Control Groups
}
\description{
In an observational complete block design, dwgtRank computes a sensitivity analysis for a weighted rank statistic designed to perform well when the comparison of a treated group and two control
groups is conducted as two nearly independent evidence factors; see Rosenbaum (2023).  For this task, a suggested setting of m, m1, m2, scores and range is given in a note in the documentation below.  A simpler way to use the suggested settings is to use the function ef2C instead.
}
\usage{
dwgtRank(y, gamma = 1, m = 2, m1 = 2, m2 = 2, phifunc = NULL,
   alternative = "greater", scores = NULL, range = TRUE)
}

\arguments{
  \item{y}{
With I blocks and J individuals in each block, y is and I x J matrix or dataframe containing the outcomes.  The first column of y is compared to columns 2, ..., J. J must be at least 2.
}
  \item{gamma}{
A real number >=1 giving the value of the sensitivity parameter.  gamma=1 yields a randomization test.
}
  \item{m}{
One of three parameters that define the weights that attach to blocks.  The three parameters are integers with 1 <= m1 <= m2 <= m.  See Details.
}
  \item{m1}{
See m.
}
  \item{m2}{
See m.
}
  \item{phifunc}{
An optional function that can be used to substitute your own weights for the weights defined by (m, m1, m2).  The function must map [0,1] into [0,1].  If phifunc is NULL, then the weight function is defined by (m, m1, m2).  If phifunc is not NULL, then it defines the weights and (m, m1, m2) are ignored.
}
  \item{alternative}{
For an upper-tailed test, use the default, alternative="greater".  For an lower-tailed test, use alternative="less".  An error will result if alternative is something besides "greater" or "less".  In this context, a two-sided test is best viewed as two one-sided tests with a Bonferroni correction, e.g., testing in both tails at level 0.025 to ensure overall level of 0.05; see Cox (1977).  For more information, see the notes.
}
  \item{scores}{
If scores is NULL, the scores are 1, 2, ..., J.  Otherwise, scores should specify the J scores for the J within-block ranks.  If scores are specified, there must be J scores, but the J scores need not be distinct.
}
  \item{range}{
If range=TRUE, then the within-block ranges are calculated, ranked from 1 to I, and scored (m, m1, m2) or phifunc.  If range=FALSE, then the within-block gap between the largest response and the average of the remaining J-1 responses is used instead.
}
}
\details{
The method uses a weighted rank statistic to compare the first column of y to the rest; see Rosenbaum (2023, 2024).  Weighted rank statistics generalize the methods of Quade (1979) and Tardif (1987).  Quade (1979) applied unscored ranks to the I within block ranges, and used unscored ranks within-blocks.  In contrast, here, the scores of ranks of ranges or gaps are based on expression (9) in Rosenbaum (2011a); see also Rosenbaum (2014) where weighted M-statistics are used instead of weighted rank statistics.  If J=2, the method agrees exactly with the method for pairs in Rosenbaum (2011a).

Using m=1, m1=1, m2=1, is the same as the stratified Wilcoxon rank sum with I strata, ignoring the ranges or gaps; see Lehmann 1975, Chapter 3).  Using m=2, m1=2, m2=2 applies unscored ranks to the I ranges or gaps.  Using m=5, m1=5, m2=5 is the suggestion of Conover and Salsburg (1988), and m=8, m1=8, m2=8 is a more extreme version of the same theme.  In pairs, J=2, m=8, m1=7, m2=8 performs well in Rosenbaum (2011a), as does m=8, m1=6, m2=8.  Detailed evaluations in terms of design sensitivity and Bahadur efficiency are in Rosenbaum (2023, 2024).
}
\value{
\item{pval }{Upper bound on the one-sided P-value.}
\item{detail }{A vector with the standarized deviate, the statistic, its null expectation and variance and the value of gamma.}
}

\note{
SUGGESTED SETTINGS FOR m, m1, m2, range AND scores WHEN USED WITH TWO CONTROL GROUPS.  These suggested settings are more conveniently implemented in the function ef2C.  Rosenbaum (2023) considered a matched block design with I blocks of size 3,
containing one treated individual and one control from each of two control
groups.  The two evidence factors are: (1) compare treated to the first control
group, and (2) compare the second control group to the pooled group that does
not distinguish the treated individual and the control from the first control
group.  For the matched pair comparison (1) with one control group,
the suggested settings are (m=8, m1=7, m2=8), together with the defaults of range=TRUE and scores=NULL; see Rosenbaum (2011a).  For comparison (2), Rosenbaum (2023) evaluated 40 statistics, judging best the
statistic with (m=8, m1=8, m2=8), range=FALSE, scores=c(1,2,5), as illustrated
below.  This statistic had good Bahadur efficiency of a sensitivity analysis
against several simple alternative hypotheses involving a treatment effect
and no unmeasured bias.

If we expect the treated group to have higher responses than controls,
then comparison (1) sets alternative to greater and comparison (2) sets
alternative to less.  If we expect the treated group to have lower responses than controls,
then comparison (1) sets alternative to less and comparison (2) sets
alternative to greater.  See also the note about alternatives.

Suppose that the data are initially in an Ix3 matrix with outcomes
for treated in the first column, control group 1 in the second column, and
control group 2 in the third column.  The dwgtRank function always compares
the first column to the remaining columns.  So, the first factor applies the
function to y[,1:2] and the second factor applies the function to y[,3:1].  Note
carefully here that y[,3:1] has reversed the order of the columns, so column
3 is compared with the other two columns.
If the treatment is expected to cause an increase in the response, then
comparison (1) applies the function to y[,1:2] with alternative =
greater, and comparison (2) applies the function to y[,3:1] with
alternative = less.  This is illustrated in the example below which
reproduces analyses from Rosenbaum (2023).  If the treatment is
expected to cause a decrease in the response, then repeat these steps
with y replaced by -y, so the treatment is expected to increase -y.
}

\note{
ALTERNATIVE.  Setting alternative to less is the same as changing the sign of
the within-block scored rank, that is, changing phi(a_ij) to -phi(a_ij) in Rosenbaum (2023); see especially equation (1) in the on-line supplement to that paper.  Note carefully that setting alternative to less does not change the between block ranks, even when range=FALSE.  In general, in dwgtRank, changing the alternative will give a different answer from changing y to -y
if range=FALSE because, unlike the range, the gap is not invariant to the sign
change.  I suggest using function ef2C -- the recommended analysis -- before trying out variations on that analysis using dwgtRank.

The dwgtRank function was designed for the evidence factor analysis
with two different control groups in blocks of size 3, and this is reflected in the way alternative is defined when range=FALSE.  For a simple analysis in the suggested form, use ef2C instead of dwgtRank; it calls dwgtRank with appropriate settings.  If you wish to explore alternative settings for this problem, use dwgtRank.  For several controls from a single control group, use wgtRank instead of dwgtRank.
}

\note{

  TIES WITHIN BLOCKS.  If there are ties within blocks, then these are resolved as follows.
If scores are not specified, so the within block ranks are intended to be 1, 2, ..., J, then average ranks are used for ties.  If scores are specified, then ties are resolved by the ties.method="min" in the rank function in base R.  This means that tied observations are all given the same rank, hence the same score, and that score corresponds with the smallest rank to which a tied group is entitled.  Suppose the scores are scores=c(1,2,5) for J=3.  If all three observations are different, then the smallest observation gets score 1, the middle gets 2, and the largest gets 5.  If the three values are, say, 16, 14, 14 in a block, they get ranks 3, 1, 1, with scores 5, 1, 1.  If the three values are 16, 16, 14 in a block, they get ranks 2, 2, 1, with scores 2, 2, 1.  This is in keeping with the idea that we want to emphasize those blocks in which one observation stands well above the rest.
  In the example, there are no within-block ties, so the issue does not arise.

  TIES BETWEEN BLOCKS.  If there are ties among the I blocks in the within-block ranges or gaps, then average ranks are used for ties.
}
\references{
Conover, W. J. and Salsburg, D. S. (1988). <doi:10.2307/2531906> Locally most powerful tests for detecting treatment effects when only a subset of patients can be expected to "respond" to treatment. Biometrics, 44, 189-196.

Cox, D. R. (1977). The role of significance tests [with discussion and reply]. Scandinavian Journal of Statistics, 4, 49-70.

Gastwirth, J. L., Krieger, A. M., and Rosenbaum, P. R. (2000).
<doi:10.1111/1467-9868.00249> Asymptotic separability in sensitivity analysis.  Journal of the Royal Statistical Society B 2000, 62, 545-556.

Lehmann, E. L. (1975). Nonparametrics: Statistical Methods Based on Ranks. San Francisco: Holden-Day.

Quade, D. (1979). <doi:10.2307/2286991> Using weighted rankings in the analysis of complete blocks with additive block effects. Journal of the American Statistical Association, 74, 680-683.

Rosenbaum, P. R. (1987). <doi:10.2307/2336017> Sensitivity analysis for certain permutation inferences in matched observational studies. Biometrika, 74(1), 13-26.

Rosenbaum, P. R. (2010). <doi:10.1093/biomet/asq019> Evidence factors in observational studies. Biometrika, 97(2), 333-345.

Rosenbaum, P. R. (2011a). <doi:10.1111/j.1541-0420.2010.01535.x> A new U‐Statistic with superior design sensitivity in matched observational studies. Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2011b). <doi:10.1198/jasa.2011.tm10422> Some approximate evidence factors in observational studies. Journal of the American Statistical Association, 106(493), 285-295.

Rosenbaum, P. R. (2013). <doi:10.1111/j.1541-0420.2012.01821.x> Impact of multiple matched controls on design sensitivity in observational studies.  Biometrics, 2013, 69, 118-127.

Rosenbaum, P. R. (2014) <doi:10.1080/01621459.2013.879261> Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls.  Journal of the American Statistical Association, 109(507), 1145-1158.

Rosenbaum, P. R. (2015). <doi:10.1080/01621459.2014.960968> Bahadur efficiency of sensitivity analyses in observational studies. Journal of the American Statistical Association, 110(509), 205-217.

Rosenbaum, P. (2017). <doi:10.4159/9780674982697> Observation and Experiment: An Introduction to Causal Inference.  Cambridge, MA: Harvard University Press.

Rosenbaum, P. R. (2018). <doi:10.1214/18-AOAS1153> Sensitivity analysis for stratified comparisons in an observational study of the effect of smoking on homocysteine
levels. The Annals of Applied Statistics, 12(4), 2312-2334.

Rosenbaum, P. R. (2021) <doi:10.1201/9781003039648> Replication and Evidence Factors in Observational Studies.  Chapman and Hall/CRC.

Rosenbaum, P. R. (2024) <doi:10.1080/01621459.2023.2221402> Bahadur efficiency of observational block designs. Journal of the American Statistical Association.

Rosenbaum, P. R. (2023) <doi:10.1111/biom.13921> A second evidence factor for a second control group.  Biometrics 79, 3968-3980.

Tardif, S. (1987). <doi:10.2307/2289476> Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82(398), 637-644.

}
\author{
Paul R. Rosenbaum
}

\note{
  This function compares the first column of y to the other columns.  To implement the second evidence factor analysis in Rosenbaum (2023), the second control group must be placed in the first column.  See the examples, where y is y[,1:2] for the first evidence factor, but y becomes y[,3:1] for the second evidence factor.  All of this is automated in the function ef2C.
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

# The following analysis contrasts the second P control
# with the pooled group consisting of the treated B
# binge drinker and the N control.  That contrast creates
# a second evidence factor, not redundant with the comparison
# of B and N, yet avoids dilution of the effect by using
# a statistic like that of Conover and Salsburg (1988)
# that allows for nonresponders.

# SECOND EVIDENCE FACTOR: CONTROL2 VS TREATED+CONTROL1
dwgtRank(y[,3:1],gamma=1.45,alternative="less",
     scores=c(1,2,5),range=FALSE,m=8,m1=8,m2=8)

amplify(1.45, 2.5)

# This is a much less sensitive result than is obtained
# from the stratified Wilcoxon rank sum statistic wiht
# I strata,

dwgtRank(y[,3:1],gamma=1.45,alternative="less",
          scores=c(1,2,3),m=1,m1=1,m2=1)

# and theory leads us to expect this difference
# in performance of the two statistics; see Rosenbaum (2023)

# EVIDENCE FACTOR ANALYSIS, COMBINING TWO FACTORS
# The evidence factor analysis compares treated to the first control group,
# then compares the second control group to the pooled group consisting of
# treated and first control, then combines the two analyses using meta-analysis.
# Treated/first-control matched pairs are compared using the method in
# Rosenbaum (2011).

p1<-dwgtRank(y[,1:2],gamma=2.3,alternative="greater",
          m=8,m1=7,m2=8)$pval
p2<-dwgtRank(y[,3:1],gamma=1.45,alternative="less",
          scores=c(1,2,5),range=FALSE,m=8,m1=8,m2=8)$pval
c(p1,p2)
sensitivitymv::truncatedP(c(p1,p2))
amplify(2.3,4)
amplify(1.45,2.5)

# THE COMBINED ANALYSIS IS INSENSITIVE TO LARGER BIASES
# The combined analysis is insensitive to larger biases
# than are its components
p1<-dwgtRank(y[,1:2],gamma=2.6,alternative="greater",
          m=8,m1=7,m2=8)$pval
p2<-dwgtRank(y[,3:1],gamma=1.7,alternative="less",
          scores=c(1,2,5),range=FALSE,m=8,m1=8,m2=8)$pval
c(p1,p2)
sensitivitymv::truncatedP(c(p1,p2))

amplify(2.6,5)
amplify(1.7,3)

# CONNECTION WITH OTHER PACKAGES
# Although dwgtRank() computes the matched pair P-value bound,
dwgtRank(y[,1:2],gamma=2.3,alternative="greater",
         m=8,m1=7,m2=8)$pval
# a simpler way to do it uses senU() in the DOS2 package
DOS2::senU(y[,1]-y[,2],m=8,m1=7,m2=8,gamma=2.3)
# where senU also provides bounds on point estimates and confidence
# intervals for each gamma
DOS2::senU(y[,1]-y[,2],m=8,m1=7,m2=8,gamma=1.5,conf.int=TRUE)

detach(aBP)
rm(p1,p2)
rm(y)

# USING SIMULATION TO GET THE GENERAL IDEAS OF DESIGN SENSITIVITY
# AND THE BAHADUR EFFICIENCY OF A SENSITIVITY ANALYSIS

# IN THIS LARGE SAMPLE SIZE, THE DESIGN SENSITIVITY PREDICTS
# U888 WILL HAVE MORE POWER THAN U555, AND IT DOES.
# SEE TABLE 2 OF ROSENBAUM (2023), NORMAL tau=1/2
# FOR U888/125/GAP AND U555/125/GAP
set.seed(1)
ss<-10000
ysim<-matrix(rnorm(3*ss),ss,3)
ysim[,1]<-ysim[,1]+sqrt(2)/2 # This is tau=1/2 for Normal errors
# Compare U888/125/gap and U555/125/gap
dwgtRank(ysim[,3:1],gamma=3,alternative="less",scores=c(1,2,5),
  range=FALSE,m=8,m1=8,m2=8)$pval
dwgtRank(ysim[,3:1],gamma=3,alternative="less",scores=c(1,2,5),
  range=FALSE,m=5,m1=5,m2=5)$pval
# IF YOU INCREASED ss FROM 10000, AS ABOVE, TO INFINITY, THE
# POWER FUNCTION WOULD TEND TO A STEP FUNCTION WITH A SINGLE
# STEP DOWN FROM POWER 1 TO POWER 0 AT THE DESIGN SENSITIVITY.

# IN THIS SMALLER SAMPLE SIZE, THE BAHADUR EFFICIENCY PREDICTS
# U555 WILL HAVE MORE POWER THAN U888, AND IT DOES.
# SEE TABLE 3 OF ROSENBAUM (2023), NORMAL tau=1/2
# FOR U888/125/GAP AND U555/125/GAP AT UPSILON = 1.5
set.seed(1)
ss<-100
ysim<-matrix(rnorm(3*ss),ss,3)
ysim[,1]<-ysim[,1]+sqrt(2)/2 # This is tau=1/2 for Normal errors
# Compare U888/125/gap and U555/125/gap
dwgtRank(ysim[,3:1],gamma=1.5,alternative="less",scores=c(1,2,5),
  range=FALSE,m=8,m1=8,m2=8)$pval
dwgtRank(ysim[,3:1],gamma=1.5,alternative="less",scores=c(1,2,5),
  range=FALSE,m=5,m1=5,m2=5)$pval

}

\keyword{ nonparametric }
\keyword{ htest }
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Evidence factors}
\concept{Multiple control groups}
\concept{rbounds}
