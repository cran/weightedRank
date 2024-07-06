\name{aHDLe}
\alias{aHDLe}
\docType{data}
\title{
HDL Cholesterol and Light Daily Alcohol 2013-2020
}
\description{
A blocked observational study of light daily alcohol consumption and HDL
cholesterol from NHANES 2013-2020.  This is an enlarged version of the aHDL data
in this package, adding data from NHANES 2017-2020.
}
\usage{data("aHDLe")}
\format{
  A data frame with 2888 observations on the following 10 variables.
  \describe{
    \item{\code{nh}}{NHANES 2013-2014 is 1314, NHANES 2015-2016 is 1516, and
    NHANES 2017-2020 is 1720}
    \item{\code{SEQN}}{NHANES ID number}
    \item{\code{age}}{Age in years}
    \item{\code{female}}{1=female, 0=male}
    \item{\code{education}}{1 is <9th grade, 3 is high school, 5 is a BA degree}
    \item{\code{z}}{1=light almost daily alcohol, 0=little or no alcohol last year.}
    \item{\code{grpL}}{Treated group and control groups.  D=Daily=light almost daily alcohol, Never=never drinks -- see note, Rarely drinks -- see note, B = PastBinge = a past history of binge drinking on most days, but currently drinks once a week or less.  For details, see Rosenbaum (2023, Appendix).}
    \item{\code{hdl}}{HDL cholesterol level mg/dL}
    \item{\code{mmercury}}{Methylmercury level ug/L}
    \item{\code{mset}}{Matched set indicator, 1, 2, ..., 722.  The 2888 observations are in 722 matched sets, each of size 4.}
  }
}
\details{
This data set enlarges the aHDL data in this package to include also data from
NHANES 2017-2020.  The aHDL data included 406 blocks of 4 people from NHANES
2013-2016, and they are included in aHDLe.  From NHANES 2017-2020,
an additional 316 blocks of 4 people were added, making 722 blocks of size 4
in total.

The alcohol questions changed slightly in NHANES 2017-2020, forcing small
changes in the definitions of two of the control groups.  See the Note for specifics.

There is a debate about whether light daily alcohol consumption -- a single glass of red wine -- shortens or lengthens life.  LoConte et al. (2018) emphasize that alcohol is a carcinogen.  Suh et al. (1992) claim reduced cardiovascular mortality brought about by an increase in high density high-density lipoprotein (HDL) cholesterol, the so-called good cholesterol.  There is on-going debate about whether there are cardiovascular benefits, and if they exist, whether they are large enough to offset an increased risk of cancer.  This example looks at a small corner of the larger debate, namely the effect on HDL cholesterol.

The example contains several attempts to detect unmeasured confounding bias, if present.  There is a secondary outcome thought to be unaffected by alcohol consumption, namely methylmercury levels in the blood, likely an indicator of the consumption of fish, not of alcohol; see Pedersen et al. (1994) and WHO (2021).  There are also three control groups, all with little present alcohol consumption, but with different uses of alcohol in the past; see the definition of variable grp above.

The appendix to Rosenbaum (2023) describes the 2013-2016 data and matching in detail.  See also the documentation in this package for the aHDL data.

}
\source{
US National Health and Nutrition Examination Survey, 2013-2014, 2015-2016
and 2017-2020.
}
\note{
There is a treated group and three control groups, with one member of each group in each of 722 blocks.  Of these, 406 blocks came from NHANES 2013-2016 and 316 blocks
came from NHANES 2017-2020, making 722=406+316 blocks in total.

In all three NHANES surveys, 2013-2014, 2015-2016 and 2017-2020, groups D and B
had the same definition, as described in the data appendix to Rosenbaum (2023).
Specifically:

D=Daily alcohol.  Drank alcohol on most days, meaning at least 260 = 5 x 52 days in the past year, drinking between 1 and 3 alcoholic drinks on drinking days.

B=Past binge drinker.  An individual in this control group had a period in
their life when they drank at least 4 or 5 drinks on most days, but quit,
and now drinks on at most one day a week, or at most 52 days in the past year.

NHANES slightly changed its alcohol questions in 2017-2020, so that small changes
were required in the definitions of the other two control groups.  NHANES 2013-2016 had asked whether a person had "12 drinks in their life," but in 2017-2020 NHANES
asked whether a person had "ever had a drink".  These were screen questions, so little was asked of a person about alcohol consumption if the person said at the outset that they had never had much alcohol.  As it turned out, a yes to the old question, "fewer than 12 drinks in their life," was fairly common, but a yes to the new question, "never had one drink," was much less common. This affects the definitions of groups N=never and R=rare drinking.

As described in greater detail in the data appendix to Rosenbaum (2023), for NHANES 2013-2016, group N=never
was comprised of people who said they had fewer than 12 drinks in their life.  In
the expansion to include data from NHANES 2017-2020, group N=never was comprised
of people who either (i) never had a drink of alcohol, or (ii) had 0 alcoholic
drinks in the past year and never in their life had a period of regular binge
drinking on most days, defined to be at least 4 or 5 drinks per day.

As described in greater detail in the data appendix to Rosenbaum (2023), for NHANES 2013-2016, group R="rare drinking" was comprised of people who did have 12 or more drinks in their life,
but had fewer than 12 drinks in the past year, and never in their life had a period of regular binge drinking on most days, defined to be at least 4 or 5 drinks per day.
In the expansion to include data from NHANES 2017-2020, group R="rare drinking"
drank between 1 and 11 alcoholic drinks in the past year and never in their life had a period of regular binge drinking on most days, defined to be at least 4 or 5 drinks per day.

In brief, everyone in control groups N and R reports drinking little or no alcohol
in the past year.  No one in groups N and R reports a period in their lives
when they engaged in binge drinking on most days.  As nearly as possible given
the questions asked, group N appears to be abstaining from alcohol, at least in the
past year.  As nearly as possible given the questions asked, group R appears
to be open to having a rare drink at a wedding or a wake, but drank a negligible
amount of alcohol in the past year.  Yet, the exact definitions of control groups N and R are slightly different for NHANES 2017-2020. In each block, everyone was given the same structured questions by NHANES, so differences between questions occur between
blocks, not within blocks.
}
\references{
LoConte, N. K., Brewster, A. M., Kaur, J. S., Merrill, J. K., and Alberg, A. J. (2018). Alcohol and cancer: a statement of the American Society of Clinical Oncology. Journal of Clinical Oncology 36, 83-93.
<doi:10.1200/JCO.2017.76.1155>

Pedersen, G. A., Mortensen, G. K. and Larsen, E. H. (1994) Beverages
as a source of toxic trace element intake. Food Additives and Contaminants, 11, 351–363. <doi:10.1080/02652039409374234>

Rosenbaum, P. R. (1987) <doi:10.1214/ss/1177013232> The role of a second control group in an observational study. Statistical Science, 2, 292-306.  Discusses multiple
control groups, as in this example.

Rosenbaum, P. R. (1989) <doi:10.2307/2531497> The role of known effects in observational studies. Biometrics, 45, 557-569.  Discusses a known effect,
such as the effect of alcohol on methylmercury.

Rosenbaum, P. R. (1989) <doi:10.1214/aos/1176347131> On permutation tests for hidden biases in observational studies. The Annals of Statistics, 17, 643-653.  Abstractly discusses multiple control groups and known effects.


Rosenbaum, P. R. (2014) <doi:10.1080/01621459.2013.879261> Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls. Journal of the American Statistical Association, 109(507), 1145-1158.  An alternative to weighted rank statistics for weighted block analyses.


Rosenbaum, P. R. (2023) <doi:10.1111/biom.13558> Sensitivity analyses informed by tests for bias in observational studies. Biometrics 79, 475–487. Uses the known effect of alchol on methylmercury.  Data appendix contains detail information about the
NHANES data.

Rosenbaum, P. R. (2024) <doi:10.1080/01621459.2023.2221402> Bahadur efficiency of observational block designs. Journal of the American Statistical Association.  Discusses properties of weighted rank statistics.

Suh, I., Shaten, B. J., Cutler, J. A., and Kuller, L. H. (1992)
<doi:10.7326/0003-4819-116-11-881> Alcohol use and mortality from coronary heart disease: the role of high-density lipoprotein cholesterol. Annals of Internal Medicine 116, 881-887.

World Health Organization (2021). Mercury and Health,
<https://www.who.int/news-room/fact-sheets/detail/mercury-and-health>, (Accessed 30 August 2021).
}
\examples{
data(aHDLe)
nh20172020<-aHDLe$nh==1720
table(nh20172020)
table(nh20172020,aHDLe$grpL)
par(mfrow=c(1,2))
boxplot(aHDLe$hdl[!nh20172020]~aHDLe$grpL[!nh20172020],main="NHANES 2013-2016",
    ylim=c(15,230),ylab="HDL Cholesterol",xlab="Group",las=1)
boxplot(aHDLe$hdl[nh20172020]~aHDLe$grpL[nh20172020],main="NHANES 2017-2020",
    ylim=c(15,230),ylab="HDL Cholesterol",xlab="Group",las=1)
par(mfrow=c(1,1))
}
\keyword{datasets}
