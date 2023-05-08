\name{aBP}
\alias{aBP}
\docType{data}
\concept{Observational study}
\concept{Causal inference}
\concept{Multiple control groups}
\concept{Evidence factors}
\concept{Second control group}
\title{
Binge Drinking and Blood Pressure
}
\description{
A matched observational study from NHANES with two control groups, examining the possible effects of on blood pressure of frequent binge drinking of alcohol.
}
\usage{data("aBP")}
\format{
  A data frame with 621 observations on the following 13 variables.
  \describe{
    \item{\code{SEQN}}{NHANES identification number}
    \item{\code{age}}{Age in years}
    \item{\code{female}}{1=female, 0=male}
    \item{\code{education}}{Education, with levels: "<9th" = less than 9th grade, "9-11" = grades 9 to 11, "HS" = high school, "SomeCol" = Some College, ">=BA" = BA degree or more}
    \item{\code{bmi}}{BMI or body-mass index}
    \item{\code{waisthip}}{Waist-to-hip ratio}
    \item{\code{vigorR}}{Engages in vigorous recreational exercise, 1=yes, 0=no}
    \item{\code{smokenow}}{Do you smoke now? Answers: Everyday, Some days, No}
    \item{\code{bpRX}}{Reports currently taking medication for high blood pressure}
    \item{\code{bpSystolic}}{Systolic blood pressure, mm Hg.  Average of up to three readings.}
    \item{\code{bpDiastolic}}{Diastolic blood pressure, mm Hg.  Average of up to three readings.}
    \item{\code{group}}{Drinking group, B=currently engages in frequent binge drinking, N=never binged regularly, and drank at most one drink per week in the last year, P=binged on most days for some period in the past but stopped, nevery binged in the last year, and drank at most one drink per week in the last year.  See Details.}
    \item{\code{mset}}{Matched set indicator, 1, 2, ..., 207.  There are 207 blocks of size 3, each containing one B, one N and one P.}
  }
}

\note{
The data are used as an example in Rosenbaum (2023).
}
\details{
    The data are from data from the 2017-2020 National Health and Nutrition Examination Survey (which was interrupted by COVID-19, so it is not a survey).  There were 5624 people who were at least 20 years of age, with an alcohol use survey, blood pressure measurements and covariates used here.  Blood pressure measurements are the average of up to three measurements.  One question asked about binge drinking in the past, defined as 4 drinks for women or 5 drinks for men.  Question ALQ151 asks (essentially): "Was there ever a time or times in your life when you drank 4/5 or more drinks of any kind of alcoholic beverage almost every day?"  Another question ALQ142 asked about binge drinking last year: "During the past 12 months, about how often did you have 4/5 or more drinks of any alcoholic beverage?"  Question ALQ121 about the overall frequency of alcohol consumption in the past 12 months. Proper use of ALQ121 and ALQ142 accounts for certain screening questions. By definition, group "binge" responded by saying that they engaged in binge drinking on 3 or more days each week in the past 12 months.  By definition, group "never" responded to ALQ142 saying they never binged in the past 12 months, responded to ALQ151 saying they had no past time when they binged almost every day, and drank any alcohol on at most one day a week in the past 12 months.  By definition, group "past" said yes to question ALQ151, so there was a period in their life of binge drinking almost every day, but they never binged in the past 12 months, and drank alcohol on at most one day a week in the past 12 months.
    There were 9232 people aged 20 or more.  Of these, 9187 had covariate information, aside from BMI and waist/hip ratio.  Of these, 7876 had an alcohol survey.  Of these, 7281 had at least one measurement of diastolic and systolic blood pressure.  Of these, 7076 had body measurements, namely BMI and waist/hip ratio.  The three treatment groups --- binge, never and past --- are mutually exclusive but not exhaustive, and 5624 people fell in one of the groups.  All 207 members of the binge group were matched to one control from each control group.  Before matching, the never group had 3995 people and the past group had 505 people.  The never group was large enough to closely match two or three controls to each member of the binge group, but that was not done in this illustrative example.
    Up to 3 repeated measures of blood pressure were often present, and the analysis uses their average.
}
\source{
US National Health and Nutrition Examination Survey
(https://www.cdc.gov/nchs/nhanes/index.htm)
}
\references{
Roerecke, M., Kaczorowski, J., Tobe, S. W., Gmel, G., Hasan, O. S. and Rehm, J. (2017). <doi:10.1016/S2468-2667(17)30003-8> The effect of a reduction in alcohol consumption on blood pressure: a systematic review and meta-analysis. Lancet Public Health, 2, e108-e120.

Rosenbaum, P. R. (2023) A second evidence factor associated with a second control group.  Manuscript.
}
\examples{

# The following code creates Figure 2 in Rosenbaum (2023)

data(aBP)
attach(aBP)

yD<-t(matrix(bpDiastolic,3,207))
yS<-t(matrix(bpSystolic,3,207))
vS<-c(yS[,1]-yS[,2],yS[,1]-yS[,3],yS[,2]-yS[,3])
vD<-c(yD[,1]-yD[,2],yD[,1]-yD[,3],yD[,2]-yD[,3])
y<-(yD/median(abs(vD)))+(yS/median(abs(vS)))


par(mfrow=c(1,3))
graphics::boxplot(yD[,1]-yD[,2],yD[,1]-yD[,3],yD[,2]-yD[,3],las=1,
        main="",ylab="Difference mm Hg",
        names=c("B-N","B-P","N-P"),cex.main=.9,
        cex.axis=.8,cex.lab=.9,xlab="Diastolic Difference")
graphics::abline(h=0)
wx<-round(stats::wilcox.test(yD[,1]-yD[,2],conf.int=TRUE)$conf.int,1)
graphics::segments(1,wx[1],1,wx[2],col="black",lwd=2)
wx<-round(stats::wilcox.test(yD[,1]-yD[,3],conf.int=TRUE)$conf.int,1)
graphics::segments(2,wx[1],2,wx[2],col="black",lwd=2)
wx<-round(stats::wilcox.test(yD[,2]-yD[,3],conf.int=TRUE)$conf.int,1)
graphics::segments(3,wx[1],3,wx[2],col="black",lwd=2)


graphics::boxplot(yS[,1]-yS[,2],yS[,1]-yS[,3],yS[,2]-yS[,3],las=1,
        main="",ylab="Difference mm Hg",
        names=c("B-N","B-P","N-P"),cex.main=.9,
        cex.axis=.8,cex.lab=.9,xlab="Systolic Difference")
graphics::abline(h=0)
wx<-round(stats::wilcox.test(yS[,1]-yS[,2],conf.int=TRUE)$conf.int,1)
graphics::segments(1,wx[1],1,wx[2],col="black",lwd=2)
wx<-round(stats::wilcox.test(yS[,1]-yS[,3],conf.int=TRUE)$conf.int,1)
graphics::segments(2,wx[1],2,wx[2],col="black",lwd=2)
wx<-round(stats::wilcox.test(yS[,2]-yS[,3],conf.int=TRUE)$conf.int,1)
graphics::segments(3,wx[1],3,wx[2],col="black",lwd=2)

graphics::boxplot(y[,1]-y[,2],y[,1]-y[,3],y[,2]-y[,3],las=1,
        main="",ylab="(Diastolic/10.7)+(Systolic/14.7)",
        names=c("B-N","B-P","N-P"),cex.main=.9,
        cex.axis=.8,cex.lab=.9,xlab="Combined Difference")
graphics::abline(h=0)
wx<-round(stats::wilcox.test(y[,1]-y[,2],conf.int=TRUE)$conf.int,1)
graphics::segments(1,wx[1],1,wx[2],col="black",lwd=2)
wx<-round(stats::wilcox.test(y[,1]-y[,3],conf.int=TRUE)$conf.int,1)
graphics::segments(2,wx[1],2,wx[2],col="black",lwd=2)
wx<-round(stats::wilcox.test(y[,2]-y[,3],conf.int=TRUE)$conf.int,1)
graphics::segments(3,wx[1],3,wx[2],col="black",lwd=2)
graphics::abline(h=0)
par(mfrow=c(1,1))
detach(aBP)
}
\keyword{datasets}
