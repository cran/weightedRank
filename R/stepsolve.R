
stepSolve<-function(f,int,eps=0.00001){
  stopifnot(is.vector(int)&(length(int)==2))
  stopifnot(int[1]<int[2])
  xL1<-int[1]
  xL2<-int[2]
  xH1<-int[1]
  xH2<-int[2]
  y1<-f(xL1)
  y2<-f(xL2)
  if (y1<=y2) stop("Error in stepSolve: f is not decreasing.")
  while ((xL2-xL1)>eps) {
    xm<-(xL2+xL1)/2
    fm<-f(xm)
    if (fm>0) xL1<-xm
    else xL2<-xm
  }
  while ((xH2-xH1)>eps) {
    xm<-(xH2+xH1)/2
    fm<-f(xm)
    if (fm<0) xH2<-xm
    else xH1<-xm
  }
  x<-c(xL1,xL2,xH1,xH2)
  list(average=mean(x),low=(xL1+xL2)/2,high=(xH1+xH2)/2)
}
