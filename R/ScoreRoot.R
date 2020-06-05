#' @importFrom Rdpack reprompt
#' @import parallel
#' @import tictoc
#' @import rootSolve
#' @export
ScoreRoot=function(t, parms){
  var = parms$var
  fitted = parms$fitted
  family = fitted$family$family
  coefB = coef(fitted)
  nonA <- !is.na(coefB)
  Pnames <- names(coefB)
  pv0 <- t(as.matrix(coefB))
  mf <- model.frame(fitted)
  ci = parms$ci
  Y <- model.response(mf)
  n <- NROW(Y)
  O <- model.offset(mf)
  if (!length(O)) O <- rep(0, n)

  W <- model.weights(mf)
  if (length(W) == 0L) W <- rep(1, n)

  X <- model.matrix(fitted)
  fam <- family(fitted)
  B <- coefB[var]
  LP <- X[, nonA, drop = FALSE] %*% coefB[nonA] + O
  a <- nonA
  a[which(names(a) == var)] <- FALSE
  Xi = X[, a, drop = FALSE]
  pi = Pnames[which(Pnames == var)]
  bi <- t
  o <- O + X[, var] * bi
  fm <- glm.fit(x = Xi, y = Y, weights = W, etastart = LP,
                offset = o, family = fam, control = fitted$control)
  LP <- Xi %*% fm$coefficients + o
  ri <- pv0
  ri[, names(coef(fm))] <- coef(fm)
  ri[, pi] <- bi
  d = length(ri)
  u = In = s = vector()
  IF = matrix(rep(0,d^2),d,d)
  r = as.vector(ri);

if(family == "binomial"){
    for (k in 1:d){
      u[k] = sum((Y - (exp(X%*%r)/(1 + exp(X%*%r))))*X[,k]);
      for (l in 1:d){
        In = (sum(X[,l]*X[,k]*exp(X%*%r)/(1 + exp(X%*%r))^2))
        IF[k,l] = In
      }}
  }else{
    stop('Unknown family.');
  }
  S = u%*%solve(IF)%*%u
  S = max(S, 0)
  z = S - qchisq(ci, df=1)
  return(z)
}
