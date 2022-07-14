wgtRanktt <-
function(y,phi1="u868",phi2="u878",phifunc1=NULL,phifunc2=NULL,gamma=1){
  
  ############## #
  # check input  #
  ################
  
  stopifnot(is.matrix(y)|is.data.frame(y))
  stopifnot(0 == sum(is.na(as.vector(y))))
  stopifnot(min(dim(y))>=2)
  stopifnot(is.vector(gamma)&(length(gamma)==1)
            &(gamma>=1))
  if (is.null(phifunc1)) {
    stopifnot(is.element(phi1,c("u868","u878",
                                "quade","wilc")))}
  if (is.null(phifunc2)) {
    stopifnot(is.element(phi2,c("u868","u878",
                                "quade","wilc")))}
  
  
  ############## #
  # subfunctions #
  ################
  
  multrnksU <- function(pk, m1 = 2, m2 = 2, m = 2) {
    # This is the right side of expression (9) 
    # in Rosenbaum (2011) Biometrics page 1022
    # Note that 0 <= pk <= 1
    n <- length(pk)
    q <- rep(0, n)
    q <- rep(0, n)
    for (l in m1:m2) {
      q <- q + (l * choose(m, l) * (pk^(l - 1)) * ((1 - pk)^(m - l)))
    }
    q/max(q)
  }
  
  u868<-function(pk){multrnksU(pk,m1=6,m2=8,m=8)}
  u878<-function(pk){multrnksU(pk,m1=7,m2=8,m=8)}
  quade<-function(pk){pk}
  wilc<-function(pk){rep(1,length(pk))}
  if (is.null(phifunc1)){
    if (phi1=="u868") phifunc1<-u868
    else if (phi1=="u878") phifunc1<-u878
    else if (phi1=="quade") phifunc1<-quade
    else if (phi1=="wilc") phifunc1<-wilc
  }
  if (is.null(phifunc2)){
    if (phi2=="u868") phifunc2<-u868
    else if (phi2=="u878") phifunc2<-u878
    else if (phi2=="quade") phifunc2<-quade
    else if (phi2=="wilc") phifunc2<-wilc
  }
  
  separable1kA <- function (ymat, gamma = 1) 
  {
    # Modified from separable1k in sensitivitymw package
    # Instead of producing the final inference,
    # this version computes the max expectation and var 
    n <- dim(ymat)[1]
    m <- dim(ymat)[2]
    o <- t(apply(ymat, 1, sort))
    allmu <- matrix(NA, n, m - 1)
    allsigma2 <- matrix(NA, n, m - 1)
    maxmu <- rep(-Inf, n)
    maxsig2 <- rep(-Inf, n)
    for (j in 1:(m - 1)) {
      pr <- c(rep(1, j), rep(gamma, m - j))/(j + ((m - j) * 
                                                    gamma))
      mu <- as.vector(o %*% pr)
      sigma2 <- as.vector((o * o) %*% pr) - (mu * mu)
      chgmu <- (mu > maxmu)
      samemu <- (mu == maxmu)
      if (sum(chgmu) > 0) {
        maxmu[chgmu] <- mu[chgmu]
        maxsig2[chgmu] <- sigma2[chgmu]
      }
      if (sum(samemu) > 0) {
        maxsig2[samemu] <- pmax(sigma2[samemu], maxsig2[samemu])
      }
    }
    list(maxmu=maxmu,maxsig2=maxsig2)
  }
  
  #######################
  # Begin main function #
  #######################
  
  
  J<-dim(y)[2]
  nset<-dim(y)[1]
  rg<-apply(y,1,max)-apply(y,1,min) # ranges within blocks
  rkrg1<-phifunc1(rank(rg)/nset) # scored ranks of ranges
  rkrg2<-phifunc2(rank(rg)/nset) # scored ranks of ranges
  rky<-t(apply(y,1,rank)) # ranks within blocks
  mv<-separable1kA(rky,gamma=gamma) # separable calculation
  ts1<-mean(rky[,1]*rkrg1) # test statistic
  ts2<-mean(rky[,1]*rkrg2) # test statistic
  ex1<-sum(mv$maxmu*rkrg1)/nset
  va1<-sum(mv$maxsig2*rkrg1*rkrg1)/(nset*nset)
  ex2<-sum(mv$maxmu*rkrg2)/nset
  va2<-sum(mv$maxsig2*rkrg2*rkrg2)/(nset*nset)
  cov12<-sum(rkrg1*rkrg2*mv$maxsig2)/(nset*nset)
  cor12<-cov12/sqrt(va1*va2)
  dev1<-(ts1-ex1)/sqrt(va1)
  dev2<-(ts2-ex2)/sqrt(va2)
  devmx<-max(dev1,dev2)
  cmat<-matrix(c(1,cor12,cor12,1),2,2)
  jointP<-1-mvtnorm::pmvnorm(lower=c(-Inf,-Inf),mean=c(0,0),
                             keepAttr=FALSE,sigma=cmat,upper=c(devmx,devmx))
  pval1<-1-pnorm(dev1)
  pval2<-1-pnorm(dev2)
  phi1<-c(pval1,dev1,ts1,ex1,va1,gamma)
  phi2<-c(pval2,dev2,ts2,ex2,va2,gamma)
  detail<-rbind(phi1,phi2)
  colnames(detail)<-c("Pval","Deviate","Statistic","Expectation","Variance","Gamma")
  list(jointP=jointP,cor12=cor12,detail=detail)
}
