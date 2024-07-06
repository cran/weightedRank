ef2C <-
function (y, gamma=1, upsilon=1, alternative="greater", trunc=0.2,
          m=c(8,8),m1=c(7,8),m2=c(8,8),scores=c(1,2,5),range=FALSE)
{
  stopifnot((alternative=="greater")|(alternative=="less"))
  stopifnot(is.matrix(y) | is.data.frame(y))
  stopifnot(0 == sum(is.na(as.vector(y))))
  stopifnot((dim(y)[2]) == 3)
  stopifnot(is.vector(gamma) & (length(gamma) == 1)
            & (gamma >=  1))
  stopifnot(is.vector(upsilon) & (length(upsilon) == 1)
            & (upsilon >=  1))
  stopifnot((length(m)==2)&(length(m1)==2)&(length(m2)==2))
  stopifnot((1<=m1[1])&(m1[1]<=m2[1])&(m2[1]<=m[1]))
  stopifnot((1<=m1[2])&(m1[2]<=m2[2])&(m2[2]<=m[2]))
  stopifnot(length(scores)==3)
  stopifnot((range==TRUE)|(range==FALSE))

  if (alternative=="less") y<- (-y)

  f1<-dwgtRank(y[,1:2],gamma=gamma,m=m[1],m1=m1[1],m2=m2[1], phifunc = NULL,
               alternative="greater", scores=NULL, range=TRUE)
  f2<-dwgtRank(y[,3:1],gamma=upsilon,m=m[2],m1=m1[2],m2=m2[2], phifunc = NULL,
               alternative="less", scores=scores, range=range)
  TreatedVSControl1<-f1$detail
  Control2vsOthers<-f2$detail
  pvals<-c(f1$pval,f2$pval,sensitivitymv::truncatedP(c(f1$pval,f2$pval),trunc=trunc))
  names(pvals)<-c("TreatedVSControl1","Control2vsOthers","Combined")
  detail<-rbind(TreatedVSControl1,Control2vsOthers)
  colnames(detail)[5]<-"Gamma/Upsilon"
  names(dimnames(detail))<-c("Evidence Factor","Test Details")
  list(pvals=pvals,detail=detail)
}
