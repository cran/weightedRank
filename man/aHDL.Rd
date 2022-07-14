\name{aHDL}
\alias{aHDL}
\docType{data}
\title{
Alcohol and HDL Cholesterol
}
\description{
A small observational study of light daily alcohol consumption and HDL cholesterol -- so-called good cholesterol -- derived from NHANES 2013-2014
and 2015-2016.  There are 406 matched sets of four individuals, making 1624 individuals in total.  Sets were matched for age, female and education in five ordered categories.
}
\usage{data("aHDL")}
\format{
  A data frame with 1624 observations on the following 11 variables.
  \describe{
    \item{\code{nh}}{NHANES 2013-2014 is 1314, and NHANES 2015-2016 is 1516}
    \item{\code{SEQN}}{NHANES ID number}
    \item{\code{age}}{Age in years}
    \item{\code{female}}{1=female, 0=male}
    \item{\code{education}}{1 is <9th grade, 3 is high school, 5 is a BA degree}
    \item{\code{z}}{1=light almost daily alcohol, 0=little or no alcohol last year.}
    \item{\code{grp}}{Treated group and control groups.  Daily=light almost daily alcohol, Never=fewer than 12 drinks during entire life, Rarely=more than 12 drinks in life, but fewer than 12 in the past year, and never had a period of daily binge drinking, PastBinge = a past history of binge drinking on most days, but currently drinks once a week or less.  For details, see Rosenbaum (2022a, Appendix).}
    \item{\code{grpL}}{Short labels for plotting formed as the first letters of grp. \code{D} < \code{N} < \code{R} < \code{B}}
    \item{\code{hdl}}{HDL cholesterol level mg/dL}
    \item{\code{mmercury}}{Methylmercury level ug/L}
    \item{\code{mset}}{Matched set indicator, 1, 2, ..., 406.  The 1624 observations are in 406 matched sets, each of size 4.}
  }
}
\details{
There is a debate about whether light daily alcohol consumption -- a single glass of red wine -- shortens or lengthens life.  LoConte et al. (2018) emphasize that alcohol is a carcinogen.  Suh et al. (1992) claim reduced cardiovascular mortality brought about by an increase in high density high-density lipoprotein (HDL) cholesterol, the so-called good cholesterol.  There is on-going debate about whether there are cardiovascular benefits, and if they exist, whether they are large enough to offset an increased risk of cancer.  This example looks at a small corner of the larger debate, namely the effect on HDL cholesterol.

The example contains several attempts to detect unmeasured confounding bias, if present.  There is a secondary outcome thought to be unaffected by alcohol consumption, namely methylmercury levels in the blood, likely an indicator of the consumption of fish, not of alcohol; see Pedersen et al. (1994) and WHO (2021).  There are also three control groups, all with little present alcohol consumption, but with different uses of alcohol in the past; see the definition of variable grp above.

The appendix to Rosenbaum (2022a) describes the data and matching in detail.  It is used as an example in Rosenbaum (2022b).

The help file for boxplotTT() applies the tail transformation to this example, reproducing a plot from Rosenbaum (2022b).

This data set is also included in the tailTransform package.  See also the informedSen package which contains a part of this data set.
}
\source{
US National Health and Nutrition Examination Survey (NHANES), 2013-2014 and 2015-2016.
}
\references{
LoConte, N. K., Brewster, A. M., Kaur, J. S., Merrill, J. K., and Alberg, A. J. (2018). Alcohol and cancer: a statement of the American Society of Clinical Oncology. Journal of Clinical Oncology 36, 83-93.
<doi:10.1200/JCO.2017.76.1155>

Pedersen, G. A., Mortensen, G. K. and Larsen, E. H. (1994) Beverages
as a source of toxic trace element intake. Food Additives and Contaminants, 11, 351–363. <doi:10.1080/02652039409374234>

Rosenbaum, P. R. (1987). The role of a second control group in an observational study. Statistical Science, 2, 292-306.
<doi:10.1214/ss/1177013232>

Rosenbaum, P. R. (1989). The role of known effects in observational studies. Biometrics, 45, 557-569. <doi:10.2307/2531497>

Rosenbaum, P. R. (1989). On permutation tests for hidden biases in observational studies. The Annals of Statistics, 17, 643-653.
<doi:10.1214/aos/1176347131>

Rosenbaum, P. R. (2014) Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls. Journal of the American Statistical Association, 109(507), 1145-1158
<doi:10.1080/01621459.2013.879261>

Rosenbaum, P. R. (2022a). Sensitivity analyses informed by tests for bias in observational studies. Biometrics. <doi:10.1111/biom.13558>

Rosenbaum, P. R. (2022b). A new transformation of treated-control matched-pair differences for graphical display.  American Statistician, to appear. <doi:10.1080/00031305.2022.2063944>

Suh, I., Shaten, B. J., Cutler, J. A., and Kuller, L. H. (1992). Alcohol use and mortality from coronary heart disease: the role of high-density lipoprotein cholesterol. Annals of Internal Medicine 116, 881-887.
<doi:10.7326/0003-4819-116-11-881>

World Health Organization (2021). Mercury and Health,
<https://www.who.int/news-room/fact-sheets/detail/mercury-and-health>, (Accessed 30 August 2021).
}
\examples{
data(aHDL)
table(aHDL$grp,aHDL$grpL) # Short labels for plotting
boxplot(aHDL$age~aHDL$grp,xlab="Group",ylab="Age")
boxplot(aHDL$education~aHDL$grp,xlab="Group",ylab="Education")
table(aHDL$female,aHDL$grpL)
table(aHDL$z,aHDL$grpL)

# The sets were also matched for is.na(aHDL$mmercury), for use
# in Rosenbaum (2022a).  About half of the matched sets
# have values for mmercury.
table(is.na(aHDL$mmercury),aHDL$grp)

# See also the informedSen package for additional analysis
}
\keyword{datasets}
\concept{Observational study}
\concept{Causal inference}
\concept{Multiple control groups}
\concept{Unaffected outcome}
\concept{Placebo outcome}
\concept{Negative control}
\concept{Control outcome}
