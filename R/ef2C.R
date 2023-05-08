ef2C <-
function (y, gamma=1, upsilon=1, alternative="greater", trunc=0.2)
{
  stopifnot((alternative=="greater")|(alternative=="less"))
  stopifnot(is.matrix(y) | is.data.frame(y))
  stopifnot(0 == sum(is.na(as.vector(y))))
  stopifnot((dim(y)[2]) == 3)
  stopifnot(is.vector(gamma) & (length(gamma) == 1)
            & (gamma >=  1))
  stopifnot(is.vector(upsilon) & (length(upsilon) == 1)
            & (upsilon >=  1))

  if (alternative=="less") y<- (-y)

  f1<-dwgtRank(y[,1:2],gamma=gamma,m=8,m1=7,m2=8, phifunc = NULL,
               alternative="greater", scores=NULL, range=TRUE)
  f2<-dwgtRank(y[,3:1],gamma=upsilon,m=8,m1=8,m2=8, phifunc = NULL,
               alternative="less", scores=c(1,2,5), range=FALSE)
  TreatedVSControl1<-f1$detail
  Control2vsOthers<-f2$detail
  pvals<-c(f1$pval,f2$pval,sensitivitymv::truncatedP(c(f1$pval,f2$pval),trunc=trunc))
  names(pvals)<-c("TreatedVSControl1","Control2vsOthers","Combined")
  detail<-rbind(TreatedVSControl1,Control2vsOthers)
  colnames(detail)[5]<-"Gamma/Upsilon"
  names(dimnames(detail))<-c("Evidence Factor","Test Details")
  list(pvals=pvals,detail=detail)
}
