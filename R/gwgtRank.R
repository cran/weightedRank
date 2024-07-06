
gwgtRank<-function(y,z,phi="u868",phifunc=NULL,gamma=1,detail=FALSE){

  # Check input
  stopifnot(is.logical(detail))
  stopifnot(is.matrix(y)|is.data.frame(y))
  stopifnot(is.matrix(z)|is.data.frame(z))
  stopifnot(0 == sum(is.na(as.vector(y))))
  stopifnot(0 == sum(is.na(as.vector(z))))
  stopifnot(all((1==as.vector(z))|(0==as.vector(z))))
  stopifnot(all(dim(y)==dim(z)))
  stopifnot(is.vector(gamma)&(length(gamma)==1)&(gamma>=1))
  stopifnot(min(dim(y)) >= 2)
  if (is.null(phifunc)) {
    stopifnot(is.element(phi, c("u868", "u878", "quade",
                                "wilc")))
  }

  # Delete blocks with no treated or no controls
  J<-dim(z)[2]
  tot<-apply(z,1,sum)
  who<-(tot>0)&(tot<J)
  if (sum(!who)>0){
    warning(paste(sum(!who)," blocks have no 1's or no 0s and were removed"))
    y<-y[who,]
    z<-z[who,]
  }
  I<-dim(z)[1]
  mset<-matrix(rep(1:I,J),I,J)

  # Define the score function from Biometrics 2011;67:1017–1027 formula (9)
  multrnksU <- function(pk, m1 = 2, m2 = 2, m = 2) {
    n <- length(pk)
    q <- rep(0, n)
    q <- rep(0, n)
    for (l in m1:m2) {
      q <- q + (l * choose(m, l) * (pk^(l - 1)) * ((1 -
                                                      pk)^(m - l)))
    }
    q/max(q)
  }
  u868 <- function(pk) {
    multrnksU(pk, m1 = 6, m2 = 8, m = 8)
  }
  u878 <- function(pk) {
    multrnksU(pk, m1 = 7, m2 = 8, m = 8)
  }
  quade <- function(pk) {
    pk
  }
  wilc <- function(pk) {
    rep(1, length(pk))
  }
  if (is.null(phifunc)) {
    if (phi == "u868")
      phifunc <- u868
    else if (phi == "u878")
      phifunc <- u878
    else if (phi == "quade")
      phifunc <- quade
    else if (phi == "wilc")
      phifunc <- wilc
  }

  # Begin calculations
  rng<-rank(apply(y,1,max)-apply(y,1,min))/(dim(y)[1]) # rank within-block ranges
  rng<-phifunc(rng) # score the range ranks
  rk<-t(apply(y,1,rank)) # within-block ranks
  o<-rk
  for (j in 1:J) o[,j]<-rk[,j]*rng # combine within block ranks and range ranks

  # Call senstrat package to obtain the sensitivity analysis
  o<-as.vector(t(o))
  z<-as.vector(t(z))
  mset<-as.vector(t(mset))
  if (detail) senstrat::senstrat(o,z,mset,detail=TRUE,gamma=gamma)
  else {
    ot<-senstrat::senstrat(o,z,mset,detail=TRUE,gamma=gamma)
    pval<-as.vector(ot$Separable[1])
    list(pval=pval,detail=ot$Separable[-1],description=ot$Description)
  }
}
