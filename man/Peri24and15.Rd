\name{Peri24and15}
\alias{Peri24and15}
\docType{data}
\title{
Periodontal Disease and Smoking in a Mixed Block Design
}
\description{
An observational study of smoking and periodontal disease in which there are
819 block, where 606 blocks contain 2 smokers and 2 controls, and 213 blocks contain 1 smoker and 4 controls.
}
\usage{data("Peri24and15")}
\format{
  A data frame with 3489 observations on the following 21 variables.
  \describe{
    \item{\code{SEQN}}{NHANES ID number}
    \item{\code{female}}{1=female, 0=male}
    \item{\code{age}}{Age in years, capped at 80 for confidentiality}
    \item{\code{ageFloor}}{Age decade = floor(age/10)}
    \item{\code{educ}}{Education as 1 to 5.  1 is less than 9th grade, 2 at least 9th grade with no high school degree, 3 is a high school degree, 4 is some college, such as a 2-year associates degree, 5 is at least a 4-year college degree.}
    \item{\code{noHS}}{No high school degree.  1 if educ is 1 or 2, 0 if educ is 3 or more}
    \item{\code{income}}{Ratio of family income to the poverty level, capped at 5 for confidenditality}
    \item{\code{nh}}{The specific NHANES survey.  A factor \code{nh0910} < \code{nh1112} < \code{nh1314}}
    \item{\code{cigsperday}}{Number of cigarettes smoked per day.  0 for nonsmokers.}
    \item{\code{z}}{Daily smoker.  1 indicates someone who smokes everyday.  0 indicates a never-smoker who smoked fewer than 100 cigarettes in their life.}
    \item{\code{pd}}{A percent indicating periodontal disease.  See details. }
    \item{\code{prop}}{A propensity score created in the example for PeriUnmatched.  This propensity score decided which smokers would have 1 control and which would have 5 controls.}
    \item{\code{pr}}{A second propensity score used to create matched pairs or matched 1-to-4 sets, after the split based on prop}
    \item{\code{mset}}{Indicator of the matched set, 1, 2, ..., 1425}
    \item{\code{treated}}{The SEQN for the smoker in this matched set.  Contains the same information as mset, but in a different form.}
    \item{\code{pair}}{1 for a matched pair, 0 for a 1-to-4 matched set}
    \item{\code{grp2}}{An ordered factor with the same information as z: S=daily smoker, N=never smoker. \code{S} < \code{N}}
    \item{\code{grp3}}{A factor with the joint information in pair and grp2.  \code{1-1:S} \code{1-1:N} \code{1-4:S} \code{1-4:N}}
    \item{\code{ageStart}}{Age at which a smoker began smoking.  See details.}
    \item{\code{packY}}{Pack-years of smoking.  Missing for nonsmokers.}
    \item{\code{block}}{Block identifiers, 1, 2, ..., 819.}
  }
}
\details{
The data in Peri24and15 rearranges the data in PeriMatched in the aamatch
package.  Specifically, the 1212 matched pairs in PeriMatched are paired
to form 606 blocks of size J=4 containing two smokers and two nonsmoking
controls.  The variable pair indicates whether an individual was in the
1212 pairs in PeriMatched or in one of the 213 1-4 matched sets.  The original matched sets in Peri24and15 are indicated by mset, but the 606+213 blocks
in Peri24and15 are indicated by block.

The construction of Peri24and15 from PeriUnmatched is in a donttest section of the example below.  The pairing of pairs is done using the nbpMatching package in R; see Greevy and Lu (2023), Lu et al. (2011) and Derigs (1988).

Additionally, Peri24and15 adds two new variables, ageStart and packY.  Seven
of 1425 smokers were missing ageStart, and the missing value was replaced by
the median starting age of 17.

Measurements were made for up to 28 teeth, 14 upper, 14 lower, excluding 4 wisdom teeth. Pocket depth and loss of attachment are two complementary measures of the degree to which the gums have separated from the teeth; see Wei, Barker and Eke (2013). Pocket depth and loss of attachment are measured at six locations on each tooth, providing the tooth is present. A measurement at a location was taken to exhibit disease if it had either a loss of attachement >=4mm or a pocked depth >=4mm, so each tooth contributes six binary scores, up to 6x28=168 binary scores.  The variable pd is the percent of these binary scores indicating periodontal disease, 0 to 100 percent.

The data from three NHANES surveys (specifically 2009-2010, 2011-2012, and 2013-2014) contain periodontal data and are used as an example in Rosenbaum (2026).  The data from one survey, 2011-2012, were used in Rosenbaum (2016).
The example replicates analyses from Rosenbaum (2026); see also the documentation for gwgtRankC in this package where further analyses are replicated.
}
\note{Analyses should distinguish blocks of different sizes, but the information they can be combined in various ways: see the documentation for gwgtRankC in this package.
In contrast, some care is required in plots and descriptive statistics.  One can straightforwardly but separately plot the blocks of the same size, but one cannot ignore block size by merging blocks of different sizes.  Suppose, however, that one merges the two treated groups from 2-to-2 blocks and 1-to-4 blocks, and merges the two control groups also; then marginal distributions of outcomes from the pooled treated and control groups are no longer comparable.  See Pimentel, Yoon and Keele (2015).  For instance, in the example, there is exact matching for sex; however, most 2-to-2 blocks are men and most 1-to-4 blocks are women.  Pool them and the pooled control group has proportionately more women than the pooled treated group.
The simple, often enlightening, solution is to plot 2-to-2 and 1-to-4 blocks in parallel but separately, and to do the same with descriptive statistics.}

\source{
US National Health and Nutrition Examination Survey (NHANES).
https://www.cdc.gov/nchs/nhanes/
}
\references{
Derigs, U. (1988) <doi:10.1007/BF02288324> Solving non-bipartite matching problems via shortest path techniques. Annals of Operations Research, 13(1), 225-261.

Greevy, R. A., & Lu, B. (2023) <doi:10.1201/9781003102670> Optimal nonbipartite matching. Handbook of Matching and Weighting Adjustments for Causal Inference, Chapman and Hall/CRC, pp. 227-238.

Lu, B., Greevy, R., Xu, X., & Beck, C. (2011) <doi:10.1198/tast.2011.08294> Optimal nonbipartite matching and its statistical applications. The American Statistician, 65(1), 21-30.

Pimentel, S. D., Yoon, F., & Keele, L. (2015) <doi:10.1002/sim.6593> Variable‐ratio matching with fine balance in a study of the Peer Health
Exchange. Statistics in Medicine, 34(30), 4070-4082.

Rosenbaum, P. R. (2015) <doi:10.1353/obs.2015.0000> Two R packages for sensitivity analysis in observational studies. Observational Studies, 1(2),
1-17.  Available on-line at: muse.jhu.edu/article/793399/summary

Rosenbaum, P. R. (2016) <doi:10.1214/16-AOAS942> Using Scheffe projections for multiple outcomes in an observational study of smoking and periondontal disease. Annals of Applied Statistics, 10, 1447-1471.

Rosenbaum, Paul R. (2025) A design for observational studies in which some people avoid treatment.  Manuscript.

Rosenbaum, Paul R. (2026)  Simple, widely applicable observational designs that improve upon the matched pairs design.  Manuscript.

Tomar, S. L. and Asma, S. (2000). Smoking attributable periodontitis in the United States: Findings from NHANES III. J. Periodont. 71, 743-751.

Wei, L., Barker, L. and Eke, P. (2013). Array applications in determining periodontal disease measurement. SouthEast SAS User's Group. (SESUG2013) Paper CC-15, analytics.ncsu.edu/ sesug/2013/CC-15.pdf.
}
\examples{
data(Peri24and15)
# There are 606 blocks of size J=4 and 213 blocks of size J=5
table(table(Peri24and15$block))
# Blocks of size J=4 contain k=2 smokers,
# and blocks of size J=5 contain k=1 smoker.
table(tapply(Peri24and15$z,Peri24and15$block,sum),table(Peri24and15$block))
# The data are analyzed in the example for the gwgtRankC() function in
# this package.

\donttest{
# This donttest portion documents the creation of Peri24and15 from PeriMatched.
# There is no need to run this portion if you only wish to use Peri24and15.
# This portion is slow but instructive if you wish to convert a paired
# design into a design with blocks of size 4 with 2 treated
# individuals and 2 controls in each block.
# The steps that follow convert the pairs in PeriMatched into
# pairs of pairs in Peri24and15 using optimal nonbipartite matching.
# The 1-to-4 matched sets are the same in PeriMatched and Peri24and15.

library(aamatch)
data("PeriMatched")
pairs<-PeriMatched[PeriMatched$pair==1,]

pairsMean<-cbind(
  tapply(pairs$female,pairs$mset,mean),
  tapply(pairs$ageFloor,pairs$mset,mean),
  tapply(pairs$age,pairs$mset,mean),
  tapply(pairs$educ,pairs$mset,mean),
  tapply(pairs$income,pairs$mset,mean))
colnames(pairsMean)<-c("female","ageFloor","age","educ","income")
npairs<-dim(pairsMean)[1]
dist<-matrix(NA,npairs,npairs)
icov<-MASS::ginv(stats::cov(pairsMean))
icov2<-MASS::ginv(stats::cov(pairsMean[,1:2]))
for (i in 1:npairs){
  mh<-stats::mahalanobis(pairsMean,pairsMean[i,],icov,inverted=TRUE)
  mh2<-stats::mahalanobis(pairsMean[,1:2],pairsMean[i,1:2],icov2,inverted=TRUE)
  dist[i,]<-mh+30*mh2
}
mset<-pairs$mset[pairs$z==1]
rownames(pairsMean)<-mset
dist<-cbind(mset,dist)
dist2<-nbpMatching::distancematrix(dist)
pp<-nbpMatching::nonbimatch(dist2)
pppairs<-NULL
halves<-pp$halves
halves1<-as.numeric(halves$Group1.ID)
halves2<-as.numeric(halves$Group2.ID)
for (i in 1:(dim(halves)[1])){
  pppairs<-rbind(pppairs,pairs[pairs$mset==halves1[i],])
  pppairs<-rbind(pppairs,pairs[pairs$mset==halves2[i],])
}
block<-as.vector(t(matrix(rep(1:(dim(halves)[1]),4),(dim(halves)[1]),4)))
pppairs<-cbind(pppairs,block)

# The pairs of pairs in pppairs and Peri24and15 are the same.
table(pppairs$SEQN==Peri24and15$SEQN[1:2424])

# The blocks are fairly homogeneous in the matched covariates.
# Of course, covariate balance is unchanged by pairing pairs,
# because the treated and control groups have not changeds.
range2<-function(v){max(v)-min(v)}
summary(tapply(pppairs$ageFloor,pppairs$block,range2))
summary(tapply(pppairs$female,pppairs$block,range2))
summary(tapply(pppairs$age,pppairs$block,range2))
summary(tapply(pppairs$educ,pppairs$block,range2))
summary(tapply(pppairs$income,pppairs$block,range2))

rm(npairs,halves1,halves2,dist2,i,block,pp,mh,mh2,mset,
   range2,icov,icov2,dist,halves)
}
}
\keyword{datasets}
