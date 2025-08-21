
estPower<-function(y,gammas,ssratio=1,phi="u868",alpha=0.05){
  stopifnot(is.matrix(y)|is.data.frame(y))
  stopifnot(is.vector(gammas)&(all(gammas>=1)))
  stopifnot(is.vector(ssratio)&(length(ssratio)==1)&(ssratio>0))
  stopifnot(is.element(phi, c("u868", "u878", "u888", "u858",
                              "quade", "wilc", "mixed")))
  stopifnot(is.vector(alpha)&(length(alpha)==1))
  stopifnot((alpha>0)&(alpha<1))
  jack<-function(y, phi = "u868"){
    I<-dim(y)[1]
    del<-rep(NA,I)
    for (i in 1:I){
      del[i]<-wgtRank(y[-i,], phi = phi)$detail[2]
    }
    delm<-mean(del)
    delv<-((I-1)/I)*sum((del-delm)^2)
    list(delm=delm,delv=delv,del=del)
  }
  o<-rep(NA,length(gammas))
  names(o)<-gammas
  crit<-qnorm(1-alpha)
  jackout<-jack(y, phi = phi)
  for (j in 1:length(gammas)){
    res<-wgtRank(y,gamma=gammas[j],phi=phi)
    detail<-res$detail
    ncp<-((detail[3]-jackout$delm)+
            crit*sqrt(detail[4])/ssratio)/
      sqrt(jackout$delv/ssratio)
    o[j]<-1-pnorm(ncp)
  }
  list(power=o,jackm=jackout$delm, jackv=jackout$delv)
}
