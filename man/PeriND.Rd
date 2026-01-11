\name{PeriND}
\alias{PeriND}
\docType{data}
\title{
New Designs for an Observational Study of Periodonal Disease and Smoking
}
\description{
A new design for data from NHANES 2009-2010, 2011-2012, 2013-2014 concerning smoking and periodontal disease.  The blocked data were built from the matched
data in PeriMatched in package aamatch by pairing its 1212 matched pairs to form 606 blocks of size four, and adding the 213 1-to-4 sets from PeriMatched.
}
\usage{data("PeriND")}
\format{
  A data frame with 3489 observations on the following 19 variables.
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
    \item{\code{prop}}{A propensity score created in the example for PeriUnmatched.  This propensity score decided which smokers would have 1 control and which would have 4 controls.}
    \item{\code{pr}}{A second propensity score used to create matched pairs or matched 1-to-4 sets, after the split based on prop}
    \item{\code{mset}}{Indicator of the matched set, 1, 2, ..., 1425}
    \item{\code{treated}}{The SEQN for the smoker in this matched set.  Contains the same information as mset, but in a different form.}
    \item{\code{pair}}{1 for a matched pair, 0 for a 1-to-4 matched set}
    \item{\code{grp2}}{An ordered factor with the same information as z: S=daily smoker, N=never smoker. \code{S} < \code{N}}
    \item{\code{grp3}}{A factor with the joint information in pair and grp2.  \code{1-1:S} \code{1-1:N} \code{1-4:S} \code{1-4:N}}
    \item{\code{block}}{Block indicators, 1 to 816, for 606 blocks of two treated individuals and two controls, and 213 blocks of one treated individual and five controls, where 816 = 606 +213.}
  }
}
\details{
  The PeriMatched data in package aamatch contains 1212 matched treated-control pairs and 213 1-to-4 matched sets.  The PeriND data uses optimal nonbipartite matching from the nbpMatching package to pair the 1212 pairs into 606 pairs-of-pairs or blocks of size 4 with 2 treated and 2 controls in each block.  For optimal nonbipartite matching, see Lu et al. (2011).

Briefly, PeriND rearranges the rows of PeriMatched and adds one variable, blocks, that identifies the 606 2-to-2 blocks and the 213 1-to-4 blocks.

The steps of the construction of PeriND from PeriMatched by optimal nonbipartite matching is given below in the example.

Measurements were made for up to 28 teeth, 14 upper, 14 lower, excluding 4 wisdom teeth. Pocket depth and loss of attachment are two complementary measures of the degree to which the gums have separated from the teeth; see Wei, Barker and Eke (2013). Pocket depth and loss of attachment are measured at six locations on each tooth, providing the tooth is present. A measurement at a location was taken to exhibit disease if it had either a loss of attachement >=4mm or a pocked depth >=4mm, so each tooth contributes six binary scores, up to 6x28=168 binary scores.  The variable pd is the percent of these binary scores indicating periodontal disease, 0 to 100 percent.

The data from three NHANES surveys (specifically 2009-2010, 2011-2012, and 2013-2014) contain periodontal data and are used as an example in Rosenbaum (2025).  The data from one survey, 2011-2012, were used in Rosenbaum (2016).
The example replicates analyses from Rosenbaum (2025).
}
\source{
US National Health and Nutrition Examination Survey (NHANES).
https://www.cdc.gov/nchs/nhanes/
}
\references{
Beck C, Lu B, Greevy R (2024). _nbpMatching: Functions for Optimal Non-Bipartite Matching_. R package version 1.5.6.
<https://CRAN.R-project.org/package=nbpMatching>.

Lu, Bo, Robert Greevy, Xinyi Xu, and Cole Beck. (2011)
<doi:10.1198/tast.2011.08294> Optimal nonbipartite matching and its statistical applications. American Statistician 65(1), 21-30.

Rosenbaum, P. R. (2007)  Sensitivity analysis for m-estimates, tests, and confidence intervals in matched observational studies. Biometrics, 63(2), 456-464. <doi:10.1111/j.1541-0420.2006.00717.x>

Rosenbaum, P. R. (2015) <doi:10.1353/obs.2015.0000> Two R packages for sensitivity analysis in observational studies. Observational Studies, 1(2), 1-17.  Available on-line at: muse.jhu.edu/article/793399/summary

Rosenbaum, P. R. (2016) <doi:10.1214/16-AOAS942> Using Scheffe projections for multiple outcomes in an observational study of smoking and periondontal disease. Annals of Applied Statistics, 10, 1447-1471.

Rosenbaum, P. R., & Small, D. S. (2017) <doi:10.1111/biom.12591> An adaptive Mantel–Haenszel test for sensitivity analysis in observational studies. Biometrics, 73(2), 422-430.

Rosenbaum, Paul R. (2025) A Design for Observational Studies in Which Some People Avoid Treatment.  Manuscript.

Tomar, S. L. and Asma, S. (2000). Smoking attributable periodontitis in the United States: Findings from NHANES III. J. Periodont. 71, 743-751.

Wei, L., Barker, L. and Eke, P. (2013). Array applications in determining periodontal disease measurement. SouthEast SAS User's Group. (SESUG2013) Paper CC-15, analytics.ncsu.edu/ sesug/2013/CC-15.pdf.
}
\examples{
data(PeriND)
# The calculations that follow show that the block, while not perfectly
# homogeneous in the covariates, are typically quite homogeneous.
# the 606+213 within-block ranges are summarized for each covariate.
range2<-function(v){max(v)-min(v)}
summary(tapply(PeriND$ageFloor,PeriND$block,range2)) # age decade
summary(tapply(PeriND$female,PeriND$block,range2)) # female indicator
summary(tapply(PeriND$age,PeriND$block,range2)) # age in years
sum(tapply(PeriND$age,PeriND$block,range2)>10)
summary(tapply(PeriND$educ,PeriND$block,range2)) # 5 categories of education
sum(tapply(PeriND$educ,PeriND$block,range2)>1)
summary(tapply(PeriND$income,PeriND$block,range2)) # income
sum(tapply(PeriND$income,PeriND$block,range2)>2)


\donttest{
rm(PeriND)
# The following code creates PeriND from PeriMatched
# using optimal nonbipartite matching in package nbpMatching.

data("PeriMatched",package="aamatch")
pairs<-PeriMatched[PeriMatched$pair==1,] # 1212 matched pairs
sets<-PeriMatched[PeriMatched$pair==0,] # 213 1-to-4 matched sets
#
# For each pair, compute the covariate mean of the two individuals
# in that pair.
pairsMean<-cbind(
  tapply(pairs$female,pairs$mset,mean),
  tapply(pairs$ageFloor,pairs$mset,mean),
  tapply(pairs$age,pairs$mset,mean),
  tapply(pairs$educ,pairs$mset,mean),
  tapply(pairs$income,pairs$mset,mean))
colnames(pairsMean)<-c("female","ageFloor","age","educ","income")
npairs<-dim(pairsMean)[1]
#
# Construct a 1212 x 1212 distance matrix between the 1212
# covariate mean vectors for the 1212 matched pairs
dist<-matrix(NA,npairs,npairs)
icov<-MASS::ginv(stats::cov(pairsMean))
icov2<-MASS::ginv(stats::cov(pairsMean[,1:2]))
for (i in 1:npairs){
  mh<-stats::mahalanobis(pairsMean,pairsMean[i,],icov,inverted=TRUE)
  mh2<-stats::mahalanobis(pairsMean[,1:2],pairsMean[i,1:2],icov2,inverted=TRUE)
  dist[i,]<-mh+30*mh2
}
# Note that this distance matrix is the sum of two distance matrices, where
# mh2 emphasizes age-decade and sex with weight 30, and mh uses all values
# of all covariates with weight 1.
#
# Set up and call nbpMatching to do the match
mset<-pairs$mset[pairs$z==1]
rownames(pairsMean)<-mset
dist<-cbind(mset,dist)
dist2<-nbpMatching::distancematrix(dist)
pp<-nbpMatching::nonbimatch(dist2)
#
# Reorganize the pairs into pairs-of-pairs using match results
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
idx<-(max(block)+1):(max(block)+(dim(sets)[1]/5))
rm(block)
#
#  Now add the 1-to-4 matched sets
block<-as.vector(t(matrix(idx,length(idx),5)))
sets<-cbind(sets,block)
PeriND<-rbind(pppairs,sets)
rm(npairs,halves1,halves2,dist2,i,block,pp,mh,mh2,mset,
   icov,icov2,dist,halves,pairs,pairsMean,pppairs,sets,idx)
dim(PeriMatched)
dim(PeriND)
length(unique(PeriND$SEQN))
sum(is.element(PeriND$SEQN,PeriMatched$SEQN))
sum(is.element(PeriMatched$SEQN,PeriND$SEQN))
}
}
\keyword{datasets}
