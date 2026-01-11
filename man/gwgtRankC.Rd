\name{gwgtRankC}
\alias{gwgtRankC}
\title{
Generalized Conditional Weighted Rank Test
}
\description{
By focusing on decisive blocks in an observational block design with I blocks of size J, the gwgtRankC() function seeks to increase design sensitivity using a conditional test.  Unlike wgtRankC() in this package, gwgtRankC() does not require each block to contain 1 treated individual and J-1 controls.
}
\usage{
gwgtRankC(y, z, gamma = 1, m = 8, m1 = 7, m2 = 8,
      m2plus = NULL, alternative="greater",warn0=FALSE)
}
\arguments{
  \item{y}{
An I x J matrix or data.frame of outcomes in a block design with I blocks of size J.  An error will result if y contains NAs.
}
  \item{z}{
z is a matrix, data.frame or a single integer.  If z is
I x J matrix or data.frame of binary treatment indicators in a block design, then z[i,j]=1 if the jth person in the ith block is treated, and
z[i,j]=0 if this person is a control.  An error will result unless all z[i,j] are 1 or 0.  If z is a single integer greater than or equal to 1 and strictly less than the number of columns of y, then an I x J binary treatment/control matrix is contructed in which the first z columns are 1 and the last J-z columns are zero, meaning the first z individuals in each block are treated and the last J-z individuals are controls.
}
  \item{gamma}{
A number greater than or equal to 1.  gamma is the sensitivity parameter.
}
  \item{m}{
See m2 below.
}
  \item{m1}{
See m2 below.
}
  \item{m2}{
The triple (m,m1,m2) defines a test statistic.  Specifically, (m,m1,m2) defines the weights that are attached to blocks.  Each coordinate of (m,m1,m2) is a positive integer.  An error will result unless m is greater than or equal to m2 which in turn is greater than or equal to m1.  See the Details. The default weights are often reasonable.
}
  \item{m2plus}{
If m2plus is not NULL, then the weights for (m,m1,m2) are added to the weights for (m,m1,m2plus).  See the Details.
}
  \item{alternative}{
alternative must equal "greater" or "less".  For alternative="greater", the null hypothesis of no effect is tested against the one-sided alternative of a treatment effect that increases responses.  To test against a treatment effect that reduces responses, set alternative="less".  See the Note.
}
  \item{warn0}{If warn0=TRUE, then a warning is printed if some within block
ranges are zero.  The number of blocks with zero ranges is also printed.  (In warn0, 0 is zero, no oh.)
  }
}
\details{
The conditional test focuses on so-called ``decisive'' blocks in which one of the following two conditions holds: (i) the highest response is treated and the lowest response is control, or (ii) the highest response is control and the lowest response is treated.  This generalizes the conditioning tactic in Rosenbaum (2025a), which is restricted to matching one treated individual to J-1 controls and is implemented in wgtRankC() in this package.

If a block has a zero within-block range, then by definition, it is not
decisive.  However, like all blocks, such a block is included when ranking the ranges.

The weight function (m,m1,m2) is derived from a U-statistic for matched pairs in Rosenbaum (2011), and was applied to block designs in Rosenbaum (2024, 2025a).  A possible alternative to the default weights is m=20, m1=19, m2=20.  Setting m=20, m1=19, m2=20 and m2plus=19 produces continuous weights that resemble with winning tailored step-function in Rosenbaum (2015): they are not monotone increasing in the within block ranges, declining slightly for the largest ranges.

The weights are calculated from the within-block ranges from the I rows of y, as in Quade (1979) and Tardiff (1987).  Quade's (1979) ranks -- ordinary ranks of the ranges -- correspond with m=2,
m1=2, m2=2.  The ranges are ranked from 1 to I, with average ranks for ties, and these ranks are scored (or transformed) using the weight function.  The default weight function is monotone increasing but pays little attention to blocks with small ranges.

Ties are handled by creating for each block a 2x2 contingency table for observations with the highest or lowest responses in the block.  In the untied case, a decisive block has row and column totals of 1 in this 2x2 table, but ties may increase the number of observations with the maximum or minimum response.  The sensitivity bound is obtained as a weighted combination of extended hypergeometric random variables with parameter gamma; see Proposition 2 in Rosenbaum (1995).
}
\value{
\item{pval }{The upper bound on the one-sided P-value testing the null hypothesis of no treatment effect against the alternative hypothesis that the treatment causes responses to increase.  The upper bound is the largest P-value that can be produced by a bias of at most gamma when the null hypothesis is
true.}
\item{detail }{A vector of five numbers: the test statistic T, its maximum expectation under the null hypothesis with a bias of at most gamma, its associated variance, the standardized deviate that is compared to the Normal distribution to calculate pval, and the value of gamma.}
\item{comp2 }{A vector of four counts: the number of blocks, the number of blocks with relevant ties, the number of decisive blocks, and the number of decisive blocks in which the highest response is a treated response.}
}
\references{
Cox, D. R. (1977) The role of significance tests (with discussion and reply). Scandinavian Journal of Statistics, 4, 49-70.

Quade, D. (1979) <doi:10.2307/2286991> Using weighted rankings in the analysis of complete blocks with additive block effects. Journal of the American Statistical Association, 74, 680-683.

Rosenbaum, P. R. (1995) <doi:10.1080/01621459.1995.10476648> Quantiles in nonrandom samples and observational
studies. Journal of the American Statistical Association, 90(432), 1424-1431.

Rosenbaum, P. R. (2011) <doi:10.1111/j.1541-0420.2010.01535.x> A new U‐Statistic with superior design sensitivity in matched observational studies. Biometrics, 67(3), 1017-1027.

Rosenbaum, P. R. (2015). <doi:10.1080/01621459.2014.960968> Bahadur efficiency of sensitivity analyses in observational studies. Journal of the American Statistical Association, 110(509), 205-217.

Rosenbaum, P. R. (2024) <doi:10.1080/01621459.2023.2221402> Bahadur efficiency of observational block designs. Journal of the American Statistical Association,
119, 1871-1881.

Rosenbaum, P. R. (2025a) <doi:10.1093/jrsssb/qkaf007> A conditioning tactic that increases design sensitivity in observational block designs.  Journal of the Royal Statistical Society Series B: Statistical Methodology, 85, 1085-1099.

Rosenbaum, P. R. (2025b) <doi:10.1007/978-3-031-90494-3> An Introduction to the Theory of Observational Studies.  Switzerland: Springer.

Rosenbaum, P. R. (2026) Two Simple, Widely Applicable Observational Designs that Improve Upon the Matched Pairs Design.  Manuscript.

Tardif, S. (1987) <doi:10.2307/2289476> Efficiency and optimality results for tests based on weighted rankings. Journal of the American Statistical Association, 82, 637-644.

Zaykin, D. V., Zhivotovsky, L. A., Westfall, P. H. and Weir, B. S. (2002)
<doi:10.1002/gepi.0042> Truncated product method of combining P-values. Genetic Epidemiology, 22, 170-185.
}
\author{
Paul R. Rosenbaum
}
\note{
Properties of the method are developed in Rosenbaum (2025c).
}

\note{
Two sided tests:  Do a two-sided test by testing in both tails
and rejecting at level alpha if the smaller P-value is less than alpha/2; see
Cox (1977, section 4.2). Setting alternative to "less" is the same as applying gwgtRankC() to -y rather than to y.
}

\note{
There is a minor detail about the relationship between gwgtRankC and wgtRankC  in this package:  They handle ties differently.  Both functions can be applied to a block design with 1 treated individual and J-1 controls in each block, but gwgtRankC() is more general: it permits more than one treated individual, or varied numbers of treated individuals, in blocks of size J.  The two functions give the same answer with untied data when there is 1 treated individual and J-1 controls, but they give ever so slightly different answers with tied data.
With 1 treated individual and J-1 controls, there is no reason to prefer one way of scoring ties to the other; so, either function can be used.  The two functions were not harmonized because wgtRankC() replicates certain published
results.
}

\examples{
# The example replicates results from Rosenbaum (2026).

data(Peri24and15)
pppairs<-Peri24and15[Peri24and15$pair==1,] # Just the pairs of pairs
perD<-t(matrix(pppairs$pd,4,606)) # Periodontal disease outcomes
z<-matrix(0,606,4)
z[,c(1,3)]<-1 # Treatment/control indicator matrix
weightedRank::gwgtRankC(perD,z,gamma=11,m1=19,m2=20,m=20)
#___________________________________________________________________
# Selected results from Table 5 of Rosenbaum (2026):

# This uses just the 606 blocks of size J=4 with k=2 smokers
# Using decisive pairs, Tc
weightedRank::gwgtRankC(perD,z,gamma=11,m1=19,m2=20,m=20)$pval
weightedRank::gwgtRankC(perD,z,gamma=7,m1=7,m2=8,m=8)$pval
weightedRank::gwgtRankC(perD,z,gamma=4.5,m1=2,m2=2,m=2)$pval

# Using weighted ranks without decisive pairs, T
weightedRank::wgtRank(perD,phi="u878",gamma=7)$pval
weightedRank::wgtRank(perD,phi="u878",gamma=5)$pval
weightedRank::wgtRank(perD,phi="u878",gamma=3)$pval

# This analyzes the same observations as pairs, not blocks of size 4.
uperD<-rbind(perD[,1:2],perD[,3:4])
weightedRank::gwgtRankC(uperD,1,gamma=7,m1=19,m2=20,m=20)$pval
weightedRank::gwgtRankC(uperD,1,gamma=5,m1=19,m2=20,m=20)$pval
weightedRank::gwgtRankC(uperD,1,gamma=3,m1=19,m2=20,m=20)$pval

#__________________________________________________________________

# Replication of Rosenbaum (2026, section 8.2)

# Combining blocks of size J=4 with k=2 smokers
# and blocks (i.e., matched set) of size J=5 with k=1 smoker

sets<-Peri24and15[Peri24and15$pair==0,] # Just the 1-to-4 matched sets
sperD<-t(matrix(sets$pd,5,213)) # Periodontal outcomes for the sets
zs<-matrix(0,213,5)
zs[,1]<-1  # Treatment indicators for the sets

# Results for sets alone.
weightedRank::gwgtRankC(sperD,zs,gamma=8.9,m1=7,m2=8,m=8)

# Combining (J,k)=(4,2) and (5,1), or all 819 blocks:

pool<-function(y1,z1,y2,z2,gamma=1,m=8,m1=7,m2=8){
  res1<-weightedRank::gwgtRankC(y1,z1,m=m,m1=m1,m2=m2,gamma=gamma)
  res2<-weightedRank::gwgtRankC(y2,z2,m=m,m1=m1,m2=m2,gamma=gamma)
  fisher<-sensitivitymv::truncatedP(c(res1$pval,res2$pval),trunc=1)
  trunc<-sensitivitymv::truncatedP(c(res1$pval,res2$pval))
  TS<-res1$detail[1]+res2$detail[1]
  EX<-res1$detail[2]+res2$detail[2]
  VA<-res1$detail[3]+res2$detail[3]
  dev<-(TS-EX)/sqrt(VA)
  pval<-1-pnorm(dev)
  detail<-c(TS,EX,VA,dev,gamma)
  names(detail)<-c("T","E(T)","var(T)","Deviate","Gamma")
  pvals<-c(pval,fisher,trunc)
  names(pvals)<-c("Pooled","Fisher","Truncated")
  list(pval=pvals,detail=detail,result1=res1,result2=res2)
}


# Pooling (J,k)=(4,2) and (J,k)=(5,1), or blocks and sets
pool(perD,z,sperD,zs,m=20,m1=19,m2=20,gamma=11)
pool(perD,z,sperD,zs,m=20,m1=19,m2=20,gamma=12)
# Pooling (J,k)=(2,1) and (J,k)=(5,1), or pairs and sets
pool(rbind(perD[,1:2],perD[,3:4]),rbind(z[,1:2],z[,3:4]),sperD,zs,
     m=20,m1=19,m2=20,gamma=6.6)

#________________________________________________________________________
#  Make Figure 1 in Rosenbaum (2026)
decisive<-function(y,z){
  stopifnot(is.matrix(y)|is.data.frame(y))
  stopifnot(is.matrix(z)|is.data.frame(z))
  stopifnot(all(dim(y)==dim(z)))
  stopifnot(all((as.vector(z)==1)|(as.vector(z)==0)))
  I<-dim(y)[1]
  J<-dim(y)[2]
  ys<-matrix(NA,I,J)
  zs<-matrix(NA,I,J)
  os<-t(apply(y,1,order))
  for (i in 1:I){
    ys[i,]<-y[i,os[i,]]
    zs[i,]<-z[i,os[i,]]
  }
  relevantTie<-(ys[,1]==ys[,2])|(ys[,J-1]==ys[,J])
  dif<-(ys[,J]-ys[,1])*(zs[,J]-zs[,1])
  zero<-abs(dif)==0
  decDif<-dif
  decDif[zero|relevantTie]<-NA
  o<-data.frame(decDif,zero,dif,relevantTie)
  if (!is.null(rownames(y))) rownames(o)<-rownames(y)
  o
}

temp<-decisive(perD,z)
pdif<-c(perD[,1]-perD[,2],perD[,3]-perD[,4])

old.par <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
barplot(table(factor(round(pdif/10),levels=(-10):10,ordered=TRUE)),
        cex.names=.8,horiz = TRUE,las=1,cex.axis=.8,xlab="Count",cex.lab=.8,
        ylab="Rounded Difference/10",xlim=c(0,400),
        main="Pair Difference",cex.main=.8)
barplot(table(factor(round(temp$decDif/10),levels=(-10):10,ordered=TRUE)),
        cex.names=.8,horiz = TRUE,las=1,cex.axis=.8,cex.lab=.8,
        xlab="Count",ylab="Rounded Difference/10",
        main="Decisive Pair Difference",cex.main=.8)
par(old.par)
rm(temp,pdif)
# Note that the plot excludes "relevant" ties for decisive pairs,
# but the analyses include "relevant" ties.
#
#______________________________________________________________________
#
# This is the evidence factor analysis in Section 8.3 of Rosenbaum (2026).
#
# First, a data frame is constructed for pack-years with one
# observation per block.
#
s<-pppairs$z==1
packyears<-tapply((pppairs$age[s]-pppairs$ageStart[s])*
          pppairs$cigsperday[s]/20,pppairs$block[s],mean)
pdDif<-tapply(pppairs$pd[s],pppairs$block[s],mean)-
       tapply(pppairs$pd[!s],pppairs$block[!s],mean)
block<-tapply(pppairs$block[s],pppairs$block[s],median)
dEF<-data.frame(block,packyears,pdDif)
rm(block,packyears,pdDif,s)
library(DOS2)
crosscutplot(dEF$packyears,dEF$pdDif)
crosscut(dEF$packyears,dEF$pdDif,gamma=3.75)

# This function will create Figure 3 from Rosenbaum (2026)
# The function slightly edits the crosscutplot function
# in the DOS2 package.
crosscutplot2<-function(x, y, ct = 0.25, xlab = "",cex=.5,rnd=0,
                        ylab = "", main = "", ylim = NULL)
{
  par(mar=c(5,4,4,4))
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(length(x) == length(y))
  stopifnot((ct > 0) & (ct <= 0.5))
  qx1 <- stats::quantile(x, ct)
  qx2 <- stats::quantile(x, 1 - ct)
  qy1 <- stats::quantile(y, ct)
  qy2 <- stats::quantile(y, 1 - ct)
  use <- ((x <= qx1) | (x >= qx2)) & ((y <= qy1) | (y >= qy2))
  if (is.null(ylim))
    graphics::plot(x, y, xlab = xlab, ylab = ylab, main = main,
                   type = "n",cex.lab=.9,cex.axis=.9,las=1)
  else graphics::plot(x, y, xlab = xlab, ylab = ylab, ylim = ylim,
                      main = main, type = "n",las=1)
  graphics::points(x[use], y[use], pch = 16, cex=cex)
  graphics::points(x[!use], y[!use], col = "gray", pch = 16,cex=cex)
  graphics::abline(h = c(qy1, qy2))
  graphics::abline(v = c(qx1, qx2))
  lines(lowess(x,y),lwd=2)
  cout<-DOS2::crosscut(x, y, ct = 0.25)
  axis(3,at=cout$quantiles[1:2],rnd,labels=round(cout$quantiles[1:2],rnd),
       cex.axis=.9,tick=TRUE)
  axis(4,at=cout$quantiles[3:4],rnd,labels=round(cout$quantiles[3:4],rnd),
       cex.axis=.9,las=1)
}



s<-sets$z==1
pdDif<-tapply(sets$pd[s],sets$mset[s],mean)-
    tapply(sets$pd[!s],sets$mset[!s],mean)
block<-(max(dEF$block)+1):(max(dEF$block)+length(pdDif))
packyears<-sets$packY[s]
Jk<-c(rep("J4k2",dim(dEF)[1]),rep("J5k1",length(pdDif)))
Jk<-as.factor(Jk)
dEF<-rbind(dEF,data.frame(block,packyears,pdDif))
dEF<-cbind(dEF,Jk)
rm(pdDif,block,packyears,Jk,s)
#
# The crosscut analysis.
DOS2::crosscut(dEF$packyears,dEF$pdDif,gamma=4)
#
# Create the plot
crosscutplot2(dEF$packyears,dEF$pdDif,ylab="Periodontal Disease",
              xlab="Pack-Years")
text(2,-50,"n=62",cex=.75)
text(2,90,"n=22",cex=.75)
text(90,90,"n=94",cex=.75)
text(90,-50,"n=37",cex=.75)

# Combine the two evidence factors using the truncated
# product of P-values (or Fisher's method if trunc=1)
joint<-function(gamma1,gamma2,trunc=.2){
  p1<-pool(perD,z,sperD,zs,m=20,m1=19,m2=20,gamma=gamma1)$pval[1]
  p2<-crosscut(dEF$packyears,dEF$pdDif,gamma=gamma2)$output$pval
  pf<-sensitivitymv::truncatedP(c(p1,p2),trunc=trunc)
  o<-c(p1,p2,pf)
  names(o)<-c("P1","P2","Pooled")
  list(pval=o)
}

joint(13,5,trunc=.2)
joint(99,3.4,trunc=.2)
joint(9.4,10,trunc=.2)

rm(dEF,perD,pppairs,sets,sperD,uperD,z,zs)
rm(crosscutplot2,decisive,joint,pool)
par(old.par)
rm(old.par)
}
\keyword{ htest }
\concept{Sensitivity analysis}
\concept{Causal inference}
\concept{Observational study}
\concept{Quade's test}
\concept{rbounds}
\concept{Conditional weighted rank test}
\concept{Conditioning tactic}
