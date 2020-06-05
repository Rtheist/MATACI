#' @importFrom Rdpack reprompt
#' @import parallel
#' @import tictoc
#' @import rootSolve
#' @export
startpoints = function(fitted, modnames, alpha){
  coefB = lapply(fitted, coef)
  aicB = lapply(fitted, function(x) x$aic)
  AICList = list()
  mn1 = length(fitted);
  pm = length(modnames)
  for (g in 1:mn1){
    AICList[[g]] = rep(aicB[[g]], pm)
    names(AICList[[g]]) = names(coefB[[g]])
  }
  AICe = Mrank1 = SMrank = w = wu = ww = wuu = indi = list()
  for (g in 1:pm){
    AICe[[g]] = lapply(AICList, function(x) unlist(x[which(names(x) == modnames[g])]))
  }
  mACIe = lapply(AICe, function(x) rep(min(unlist(x)), mn1))

  for (g in 1:pm){
    Mrank1[[g]] = Map('-', AICe[[g]], mACIe[[g]])
  }
  SMrank = lapply(Mrank1, function(x) sum(exp(-0.5*unlist(x))))
  for (g in 1:pm){
    wu[[g]] = lapply(Mrank1[[g]], function(x){
      if(length(x) == 0) {
        x = NA
      }else{
        x = x
      }})}
  for (g in 1:pm){
    wuu[[g]] = lapply(wu[[g]], function(x) exp(-0.5*((unlist(x)))))
  }
  for (g in 1:pm){
    w[[g]] = lapply(wuu[[g]], function(x) unlist(x)/SMrank[[g]])
  }
  m = matrix(unlist(w), ncol = pm)
  m[is.na(m)] = 0
  colnames(m) = modnames
  summ = lapply(fitted, summary)
  std.err = lapply(summ, function(x) as.numeric(t(x$coefficients[, "Std. Error", drop = FALSE])))
  MRX = as.data.frame(matrix(c(rep(0, (pm)*mn1)), ncol=pm))
  colnames(MRX) = c(modnames)
  MRX.se = MRX
  for (i in 1:(mn1)){
    coef = coefB[[i]]
    l = length(coef)
    se1 = std.err[[i]]
    nam = names(se1) = names(coef)
    coef = coef[modnames]
    names(coef) = modnames
    coef[is.na(coef)] = 0
    coefo = coef[which(coef != 0)]
    if (length(coefo) == 0){
      coefo = rep(0,l)
      names(coefo) = nam
    }
    MRX[i,] = coef
    se1 = se1[modnames]
    names(se1) = modnames
    se1[is.na(se1)] = 0
    MRX.se[i,] = se1
  }
  Lower = as.matrix(MRX/MRX.se - qnorm(1 - alpha, 0, 1))
  Lower[which(!is.finite(Lower))] = 0
  Lower = m*Lower
  Upper = as.matrix(MRX/MRX.se - qnorm(alpha, 0, 1))
  Upper[which(!is.finite(Upper))] = 0
  Upper = m*Upper
  W.se = as.matrix(m/MRX.se);
  W.se[which(!is.finite(W.se))] = 0
  MATAF.lower = apply(Lower, 2, sum)/apply(W.se, 2, sum)
  MATAF.lower[which(!is.finite(MATAF.lower))] = 0
  MATAF.lower = MATAF.lower[modnames]
  names(MATAF.lower) = modnames
  MATAF.lower[is.na(MATAF.lower)] = 0
  MATAF.upper = apply(Upper, 2, sum)/apply(W.se, 2, sum)
  MATAF.upper[which(!is.finite(MATAF.upper))] = 0
  MATAF.upper = MATAF.upper[modnames]
  names(MATAF.upper) = modnames
  MATAF.upper[is.na(MATAF.upper)] = 0
  MATAF.Est = apply(m*MRX, 2, sum)
  MATAF.Est = MATAF.Est[modnames]
  names(MATAF.Est) = modnames
  MATAF.Est[is.na(MATAF.Est)] = 0
  results = data.frame(Estimates = MATAF.Est, CI.lower = MATAF.lower, CI.upper = MATAF.upper);
  colnames(results) = c("Estimates", paste(round(100*alpha, 2), "%",
                                           sep=""),  paste(round(100*(1-alpha), 2), "%", sep=""))
  return(results)
}
