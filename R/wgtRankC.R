wgtRankC <-
function (y, phi = "u878", phifunc = NULL, gamma = 1,
                         alternative="greater")
{
  # Check input
  stopifnot((alternative=="greater")|(alternative=="less"))
  stopifnot(is.matrix(y) | is.data.frame(y))
  stopifnot(0 == sum(is.na(as.vector(y))))
  stopifnot(3<=dim(y)[2])
  stopifnot(is.vector(gamma) & (length(gamma) == 1))
  stopifnot(gamma >= 1)
  if (is.null(phifunc)) {
    stopifnot(is.element(phi, c("u868", "u878", "u888",
                                "quade", "wilc")))
  }

  # Define functions

  multrnksU <- function(pk, m1 = 2, m2 = 2, m = 2) {
    # Expressin (9) in Rosenbaum Biometrics 2011, 67, 1017–1027
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
  u888 <- function(pk) {
    multrnksU(pk, m1 = 8, m2 = 8, m = 8)
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
    else if (phi == "u888")
      phifunc <- u888
    else if (phi == "quade")
      phifunc <- quade
    else if (phi == "wilc")
      phifunc <- wilc
  }

  # Begin
  J <- dim(y)[2]
  nset <- dim(y)[1]
  #
  # Determine the scores that attach to each block
  # Also, remove blocks with range=0
  #
  mx<-apply(y, 1, max)
  mn<-apply(y, 1, min)
  rg <- mx-mn # within block range
  rkrg <- phifunc(rank(rg)/nset) # scores of ranks of ranges
  #
  #  Remove any blocks that are completely tied,
  #  with max=min and range=0
  #
  if (sum(mx==mn)>0){
    who<-mx!=mn
    y<-y[who,]
    mn<-mn[who]
    mx<-mx[who]
    rkrg<-rkrg[who]
    warning(paste(sum(1-who)," blocks had 0 range and were removed."))
    nset<-dim(y)[1]
  }
  #
  #  Determine rank scores for individuals
  #
  sc<-t(apply(y,1,rank)) # Uses average ranks for ties
  for (j in 1:J) sc[,j]<-sc[,j]*rkrg
  #
  #  Which individuals are max or min?  Ei=1
  #  Due to ties, there may be more than 2 in a block.
  #
  I<-nset
  E<-matrix(0,I,J)
  for (j in 1:J) E[(y[,j]==mx)|(y[,j]==mn),j]<-1
  #
  # Save some summary information for output later
  if (alternative=="greater") tmax<-sum(y[,1]==mx)
  else if (alternative=="less") tmax<-sum(y[,1]==mn)
  #
  #  Remove blocks in which treated is neither max nor min
  #
  sc<-sc[E[,1]==1,]
  E1<-E[E[,1]==1,]
  I1<-dim(sc)[1]
  #
  #  Prepare input for senstrat
  #
  #  For senstrat, z identifies treated
  z<-matrix(0,I1,J)
  z[,1]<-1
  #  For senstrat, block identifiers
  blocks<-matrix(rep(1:I,J),I,J)
  blocks<-blocks[E[,1]==1,]
  blocks<-t(blocks)[t(E1)==1]   # block ID numbers
  z1<-t(z)[t(E1)==1]
  sc1<-t(sc)[t(E1)==1]
  o<-senstrat::senstrat(sc1,z1,blocks,gamma=gamma,alternative=alternative,detail=TRUE)
  o<-o$Separable
  iblocks<-dim(E1)[1]
  tied<-sum(apply(E1,1,sum)>2)
  if (!is.na(tmax)){
    bkinfo<-c(iblocks,tmax,tied)
    if (alternative=="greater")
      names(bkinfo)<-c("Blocks Used","Treated Largest","Used Blocks with Ties")
    else names(bkinfo)<-c("Blocks Used","Treated Smallest","Used Blocks with Ties")
  }
  else{
    bkinfo<-c(iblocks,tied)
    names(bkinfo)<-c("Blocks Used","Used Blocks with Ties")
  }
  list(pval=as.vector(o[1]),detail=o,block.information=bkinfo)
}
