wgtRankCI<-function(y, phi = "u868", phifunc = NULL, gamma = 1,
                    alternative="greater", alpha=0.05, eps = 0.00001){

  ############## #
  # check input  #
  ################


  stopifnot(is.matrix(y)|is.data.frame(y))
  stopifnot(is.vector(eps)&(length(eps)==1)&(eps>0))
  stopifnot(0 == sum(is.na(as.vector(y))))
  stopifnot(min(dim(y))>=2)
  stopifnot(is.vector(gamma)&(length(gamma)==1)
            &(gamma>=1))
  stopifnot(is.vector(alpha)&(length(alpha)==1)&(alpha>0)&(alpha<1))
  stopifnot(is.element(alternative,c("greater","less","twosided")))
  if (is.null(phifunc))
    {stopifnot(is.element(phi,c("u868","u878",
                              "quade","wilc")))}

  if (alternative!="twosided") crit<-stats::qnorm(1-alpha)
  else crit<-stats::qnorm(1-(alpha/2))

  ############## #
  # subfunctions #
  ################

  deviateU<-function(taus){
    ntaus<-length(taus)
    o<-rep(NA,ntaus)
    for (i in 1:ntaus){
      yt<-y
      yt[,1]<-yt[,1]-taus[i]
      o[i]<-wgtRank(yt,phi=phi,phifunc=phifunc,gamma=gamma)$detail[1]
    }
    round(o-crit,5) #edit
  }

  deviateL<-function(taus){
    ntaus<-length(taus)
    o<-rep(NA,ntaus)
    for (i in 1:ntaus){
      yt<-(-y)
      yt[,1]<-yt[,1]-(-taus[i])
      o[i]<-wgtRank(yt,phi=phi,phifunc=phifunc,gamma=gamma)$detail[1]
    }
    -round(o-crit,5)  # edit
  }

  deviateU0<-function(taus){
    ntaus<-length(taus)
    o<-rep(NA,ntaus)
    for (i in 1:ntaus){
      yt<-y
      yt[,1]<-yt[,1]-taus[i]
      o[i]<-wgtRank(yt,phi=phi,phifunc=phifunc,gamma=gamma)$detail[1]
    }
    round(o,5)
  }

  deviateL0<-function(taus){
    ntaus<-length(taus)
    o<-rep(NA,ntaus)
    for (i in 1:ntaus){
      yt<-(-y)
      yt[,1]<-yt[,1]-(-taus[i])
      o[i]<-wgtRank(yt,phi=phi,phifunc=phifunc,gamma=gamma)$detail[1]
    }
    round(-o,5)  # edit
  }

  ################
  # Begin ########
  ################

  rg<-c(min(y[,1])-max(y[,-1]),max(y[,1])-min(y[,-1])) # range of taus

  CIl<-(-Inf)
  CIh<-Inf


  if (alternative!="less") CIl<-stepSolve(deviateU,rg,eps=eps)$low

  if (alternative!="greater") CIh<-stepSolve(deviateL,rg,eps=eps)$high

  HLl<-stepSolve(deviateU0,rg,eps=eps)$average

  HLh<-stepSolve(deviateL0,rg,eps=eps)$average

  CI<-c(CIl,CIh)
  names(CI)<-c("Low","High")
  HL<-c(HLl,HLh)
  names(HL)<-c("Low","High")

  list(estimate=HL,confidence=CI)
}
