#' @importFrom Rdpack reprompt
#' @import parallel
#' @import tictoc
#' @import rootSolve
#' @export
WaldCI = function (fitted, fm, alpha = 0.05, startL, startU) {
  coefB = lapply(fitted, coef)
  aicB = lapply(fitted, function(x) x$aic)
  alpha = alpha/2;
  AICList = list()
  for (g in 1:length(coefB)){
    len = length(coefB[[g]])
    re = rep(aicB[[g]], len)
    AICList[[g]] = re
    names(AICList[[g]]) = names(coefB[[g]])
  }
  nonA = lapply(coefB, function(x) !is.na(x))
  Pnames = lapply(coefB, names)
  pv0 = lapply(coefB, function(x) t(as.matrix(x)))
  p = lapply(Pnames, length)
  which = lapply(p, function(x) 1:x)
  summ = lapply(fitted, summary)
  std.err = lapply(summ, function(x) as.numeric(t(x$coefficients[, "Std. Error", drop = FALSE])))
  for (i in 1:length(nonA)){
    names(std.err[[i]]) = Pnames[[i]]
  }
  mf = lapply(fitted, model.frame)
  mfy = model.frame(fm)
  Fnames = names(B0 <- coef(fm))
  Fwhich = 1:length(Fnames)
  Y = model.response(mfy)
  n = NROW(Y)
  O = lapply(mf, model.offset)
  O = lapply(O, function(x){
    if (!length(x)) x = rep(0, n)})
  W = lapply(mf, model.weights)
  W = lapply(W, function(x){
    if (length(x) == 0L) x = rep(1, n)})
  X = lapply(fitted, model.matrix)
  fam = lapply(fitted, function(x) x$family[1])
  scor = vector("list", length = length(Fwhich))
  names(scor) = Fnames[Fwhich]
  for (i in Fwhich) {
    a = nonA
    var = Fnames[i]
    an = length(a)
    AICe = lapply(AICList, function(x) x[which(names(x) == var)])
    pe = lapply(coefB, function(x) x[which(names(x) == var)])
    spe = lapply(std.err, function(x) x[which(names(x) == var)])
    tl = startL[i]
    tu = startU[i]
    parametersu = list(an=an,AICe=AICe,var=var,nonA=nonA,n=n,coefB=coefB,O=O,a=a,X=X,p=p,Pnames=Pnames,W=W,fitted=fitted,fam=fam,pv0=pv0,alpha=alpha,errb="u",pointest=pe,spe=spe)
    parametersl = list(an=an,AICe=AICe,var=var,nonA=nonA,n=n,coefB=coefB,O=O,a=a,X=X,p=p,Pnames=Pnames,W=W,fitted=fitted,fam=fam,pv0=pv0,alpha=alpha,errb="l",pointest=pe,spe=spe)
    topot = function(t, parms){
      an = parms$an
      var = parms$var
      nonA = parms$nonA
      coefB = parms$coefB
      O = parms$O
      a = parms$a
      X = parms$X
      n = parms$n
      errb = parms$errb
      Pnames = parms$Pnames
      AICe = parms$AICe
      W = parms$W
      fitted = parms$fitted
      fam = parms$fam
      pv0 = parms$pv0
      alpha = parms$alpha
      pointest = parms$pointest
      spe = parms$spe
      p = parms$p
      Xi=pi=bi=d=u=IF=LP=LPm=mark=B=LPj=S=Sm=SE=r=ri=o=fm=indi=list()
      maicB = min(unlist(AICe))
      Mrank1 = Map('-', AICe, maicB)
      SMrank = sum(exp(-0.5*unlist(Mrank1)))
      w = lapply(Mrank1, function(x) exp(-0.5*(x))/SMrank)
      for (g in 1:an){
        if(length(pointest[[g]]) == 0) w[[g]] = 0
      }
      for (j in 1:an){
        if (is.element(var,names(a[[j]])) == T){
          B[[j]] = coefB[[j]][var]
          SE[[j]] = spe[[j]]
          bi[[j]] = t
          S[[j]] = Sm[[j]] = (B[[j]] - bi[[j]])/SE[[j]]
        }else{
          S[[j]] = 0
          Sm[[j]] = 0
        }
        Sm[[j]] = pt(Sm[[j]], df = n - length(coefB[[j]]))
      }
      if (errb == "u"){
        z = sum(unlist(Sm)*unlist(w)) - alpha
      }else{
        z = sum((1 - unlist(Sm))*unlist(w)) - alpha
      }
      return(z)
    }
    options(warn = 2)
    test = try(multiroot(topot, start = tl, maxiter = 1000, useFortran = T,
                         parms = parametersl))
    options(warn = 1)
    war = inherits(test, "try-error")
    if (war == T){
      l = uniroot(topot, lower = tl, upper = tu, extendInt = "yes",
                  maxiter = 1000, tol = 1e-5, parms = parametersl)$root
    }else{
      l = multiroot(topot, start = tl, maxiter = 1000, useFortran = T,
                    ctol = 1e-5, parms = parametersl)$root
    }
    options(warn = 2)
    test = try(multiroot(topot, start = tu, maxiter = 1000, useFortran = T,
                         parms = parametersu))
    options(warn = 1)
    war = inherits(test, "try-error")
    if (war == T){
      u = uniroot(topot, lower = tl, upper = tu, maxiter = 1000, extendInt = "yes",
                  tol = 1e-5, parms = parametersu)$root
    }else{
      u = multiroot(topot, start = tu, maxiter = 1000, useFortran = T,
                    ctol = 1e-5, parms = parametersu)$root
    }
    med = u
    if (l > u){
      u = l;
      l = med;
    }
    scor[[i]] = data.frame(rbind(l, u))
    colnames(scor[[i]]) = NULL
  }
  CI = matrix(rep(0, 2*length(Fwhich)), ncol = 2)
  for(i in Fwhich){
    CI[i,] = t(scor[[i]])
  }
  rownames(CI) = Fnames
  colnames(CI) = c(paste(round(100*alpha, 2), "%", sep=""),
                   paste(round(100*(1-alpha), 2), "%", sep=""))
  return(CI)
}
