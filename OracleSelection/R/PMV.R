PMV <- function(Bmat,ObjectGenoP1,ObjectGenoP2,ObjectMap){
  # Variance-covariance matrix of genotypes in progeny
  Gcov <- GenCovProgeny(ObjectGenoP1,ObjectGenoP2,ObjectMap)
  # Posterior variance-covariance matrix of marker effects
  vB <- 1/nrow(Bmat)*crossprod(scale(Bmat,TRUE,FALSE))
  # Posterior mean of marker effects
  Bhat <- as.matrix(apply(Bmat,2,mean))
  # Posterior mean variance in progeny
  PMV <- sum(diag(Gcov %*% vB)) + t(Bhat) %*% Gcov %*% Bhat
  return(PMV)
}