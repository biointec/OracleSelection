UC <- function(Bmat,ObjectGenoP1,ObjectGenoP2,ObjectMap,Psel=0.05,h=1){
  Bhat <- colMeans(Bmat)
  i <- dnorm(qnorm(1-Psel))/Psel
  # Posterior mean variance
  VarG <- PMV(Bmat,ObjectGenoP1,ObjectGenoP2,ObjectMap)
  # Mean GEBVs of parents
  mu <- mean(crossprod(ObjectGenoP1,Bhat),crossprod(ObjectGenoP2,Bhat))
  return(mu + i*h*sqrt(VarG))
}