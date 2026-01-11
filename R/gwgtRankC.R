
gwgtRankC<-function(y,z,gamma=1,m=8,m1=7,m2=8,m2plus=NULL,
                    alternative="greater",warn0=FALSE){
  # Define functions
  multrnksU <- function(pk, m1 = 2, m2 = 2, m = 2) {
    n <- length(pk)
    q <- rep(0, n)
    q <- rep(0, n)
    for (l in m1:m2) {
      q <- q + (l * choose(m, l) * (pk^(l - 1)) * ((1 -
                                                      pk)^(m - l)))
    }
    q
  }
  # End define functions
  # Check input
  stopifnot(is.logical(warn0))
  stopifnot((alternative=="greater")|(alternative=="less"))
  if (alternative=="less") y<-(-y)
  if (is.null(m)){
    m<-20
    m1<-19
    m2<-20
    m2plus<-19
  }
  else stopifnot((!is.null(m1))&(!is.null(m2)))
  stopifnot((m>=m2)&(m2>=m1)&(m1>=1))
  if (!is.null(m2plus)) stopifnot((m>=m2plus)&(m2plus>=m1))
  stopifnot(is.matrix(y)|is.data.frame(y))
  stopifnot(all(!is.na(as.vector(y))))
  if (length(as.vector(z))==1){
    stopifnot((z>=1)&(z<(dim(y)[2])))
    oldz<-z
    z<-matrix(0,dim(y)[1],dim(y)[2])
    z[,1:oldz]<-1
  }
  stopifnot(is.matrix(z)|is.data.frame(z))
  stopifnot(all(!is.na(as.vector(z))))
  stopifnot(all(dim(y)==dim(z)))
  I<-dim(y)[1]
  oldI<-I
  J<-dim(y)[2]
  stopifnot(all((as.vector(z)==1)|(as.vector(z)==0)))
  # End check input
  #
  # Begin computations
  yo<-t(apply(y,1,order))
  ys<-t(apply(y,1,sort))
  zs<-matrix(NA,I,J)
  for (i in 1:I) zs[i,]<-z[i,yo[i,]]
  # Compute ranks of block ranges
  rg<-ys[,J]-ys[,1]
  rkrg<-rank(rg)/length(rg)
  sc<-multrnksU(rkrg,m1=m1,m2=m2,m=m)
  if (!is.null(m2plus)) sc<-sc+multrnksU(rkrg,m1=m1,m2=m2plus,m=m)
  sc<-sc/max(sc)
  # Remove any blocks that are completely tied
  who<-(rg>0)
  I<-sum(who)
  if (I<oldI){
    yo<-yo[who,]
    ys<-ys[who,]
    zs<-zs[who,]
    rg<-rg[who]
    rkrg<-rkrg[who]
    sc<-sc[who]
    if (warn0) warning(paste(sum(1-who)," blocks had 0 range and were removed."))
  }
  # Find all the maximum and minumum responses, even if tied
  mx<-matrix(0,I,J)
  for (j in 1:J) mx[ys[,j]==ys[,J],j]<-1
  mn<-matrix(0,I,J)
  for (j in 1:J) mn[ys[,j]==ys[,1],j]<-1
  bkTies<-sum(apply(mn+mx,1,sum)>2)

  # Determine decisive pairs
  dec<-(rg>0)&(apply((mn+mx)*zs,1,sum)>=1)&(apply((mn+mx)*(1-zs),1,sum)>=1)
  ndec<-sum(dec)
  dsc<-sc[dec]
  mn<-mn[dec,]
  mx<-mx[dec,]
  dz<-zs[dec,]
  w<-((mn+mx)==1)
  TS<-0
  EX<-0
  VR<-0
  nn<-apply(w*dz,1,sum)
  nm1<-apply(mx*w,1,sum)
  nm2<-apply(mn*w,1,sum)
  cnr<-apply(mx*dz,1,sum)
  dsc<-dsc/nn
  for (i in 1:ndec){
    TS<-TS+cnr[i]*dsc[i]
    EX<-EX+BiasedUrn::meanFNCHypergeo(nm1[i], nm2[i], nn[i], gamma)*dsc[i]
    VR<-VR+BiasedUrn::varFNCHypergeo(nm1[i], nm2[i], nn[i], gamma)*dsc[i]*dsc[i]
  }
  dev<-(TS-EX)/sqrt(VR)
  detail<-c(TS,EX,VR,dev,gamma)
  names(detail)<-c("T","E(T)","var(T)","Deviate","Gamma")
  pval<-1-pnorm(dev)
  treatedEx<-sum(apply(mx*dz,1,sum)>=1)
  BlockCounts<-c(oldI,bkTies,ndec,treatedEx)
  if (alternative=="greater") {
    names(BlockCounts)<-c("Blocks","Relevant Ties","Decisive",
                          "  A Treated is Max")
  }
  else {
    names(BlockCounts)<-c("Blocks","Relevant Ties","Decisive",
                          "  A Treated is Min")
  }
  list(pval = pval, detail = detail, block.counts = BlockCounts)
}
