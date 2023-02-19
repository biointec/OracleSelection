## Get Covariance matrix for UC
#-------------------------------------------------------------------------
#Copyright (c) 2020 Antoine Allier
#-------------------------------------------------------------------------

GenCovProgeny <- function(ObjectGenoP1,ObjectGenoP2,ObjectMap){
  # Computing expected frequency of recombinants c1
  myDist <- sapply(1:nrow(ObjectMap),
                   function(x) abs(ObjectMap$POSITION[x] - ObjectMap$POSITION))
  myCHR1 <- do.call(rbind,lapply(1:nrow(ObjectMap),function(x)
    rep(ObjectMap$CHROMOSOME[x],nrow(ObjectMap))))
  myCHR2 <- t(myCHR1)
  c1 <- 0.5*(1-exp(-2*(myDist/100)))
  c1[myCHR1!=myCHR2] <- 0.5
  Recomb <- (1-2*c1)
  # Computing diseqilibrium
  tmp <- rbind(ObjectGenoP1,ObjectGenoP2)/2
  D <- crossprod(scale(tmp,scale=F))/2
  # Sigma matrix
  MyVarCov <- 4*D*Recomb
  return(MyVarCov)
}