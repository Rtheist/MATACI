#' @importFrom Rdpack reprompt
#' @import parallel
#' @import tictoc
#' @import rootSolve
#' @export
ScoreCI = function (fitted, fm, alpha = 0.05, startL, startU){
  coefB = lapply(fitted, coef)
  nonA = lapply(coefB, function(x) !is.na(x))
  aicB = lapply(fitted, function(x) x$aic)
  AICList = list()
  alpha = alpha/2;
  for (g in 1:length(coefB)){
    len = length(coefB[[g]])
    re = rep(aicB[[g]], len)
    AICList[[g]] = re
    names(AICList[[g]]) = names(coefB[[g]])
  }
  Pnames = lapply(coefB, names)
  pv0 = lapply(coefB, function(x) t(as.matrix(x)))
  p = lapply(Pnames, length)
  which = lapply(p, function(x) 1:x)
  summ = lapply(fitted, summary)
  mf = lapply(fitted, model.frame)
  mfy = model.frame(fm)
  Fnames=names(B0 <- coef(fm))
  Fwhich=1:length(Fnames)
  Y = model.response(mfy)
  n = NROW(Y)
  O = lapply(mf, model.offset)
  O = lapply(O, function(x){
    if (!length(x)) x = rep(0, n)})
  W = lapply(mf, model.weights)
  W = lapply(W, function(x){
    if (length(x) == 0L) x = rep(1, n)})
  X = lapply(fitted, model.matrix)
  fam = lapply(fitted, function(x) x$family)
  for (i in 1:length(fam)){
    switch(fam[[i]]$family, binomial = , poisson = , 'Negative Binomial' = {
      zmax = sqrt(qchisq(1 - alpha, 1))
      profName = 'z'
    }, gaussian = , quasi = , inverse.gaussian = , quasibinomial = ,
    quasipoisson = , {
      zmax = sqrt(qf(1 - alpha, 1, n - p[[i]]))
      profName = 'tau'
    })}
  scor = vector("list", length = length(Fwhich))
  names(scor) = Fnames[Fwhich]
  for (i in Fwhich) {
    a = nonA
    var = Fnames[i]
    an = length(a)
    AICe = lapply(AICList, function(x) x[which(names(x) == var)])
    pe = lapply(coefB, function(x) x[which(names(x) == var)])
    tl = startL[i]
    tu = startU[i]
    parametersu = list(an=an,var=var,AICe=AICe,nonA=nonA,coefB=coefB,O=O,a=a,X=X,Pnames=Pnames,W=W,fitted=fitted,fam=fam,pv0=pv0,alpha=alpha,errb="u",pointest=pe)
    parametersl = list(an=an,var=var,AICe=AICe,nonA=nonA,coefB=coefB,O=O,a=a,X=X,Pnames=Pnames,W=W,fitted=fitted,fam=fam,pv0=pv0,alpha=alpha,errb="l",pointest=pe)
    topot = function(t, parms){
      an = parms$an
      var = parms$var
      nonA = parms$nonA
      errb = parms$errb
      coefB = parms$coefB
      O = parms$O
      a = parms$a
      X = parms$X
      Pnames = parms$Pnames
      W = parms$W
      fitted = parms$fitted
      fam = parms$fam
      pv0 = parms$pv0
      alpha = parms$alpha
      pointest = parms$pointest
      Xi=pi=bi=d=u=IF=IFlm=LP=mark=B=S=Sm=r=ri=o=fm=indi=list()
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
          LP[[j]] = X[[j]][, nonA[[j]], drop = FALSE] %*% coefB[[j]][nonA[[j]]] + O[[j]]
          a[[j]][which(names(a[[j]]) == var)] = FALSE
          Xi[[j]] = X[[j]][, a[[j]], drop = FALSE]
          pi[[j]] = Pnames[[j]][which(Pnames[[j]] == var)]
          bi[[j]] = t
          o[[j]] = O[[j]] + X[[j]][, var] * bi[[j]]
          fm[[j]] = glm.fit(x = Xi[[j]], y = Y, weights = W[[j]], etastart = LP[[j]],
                             offset = o[[j]], family = fam[[j]], control = fitted[[j]]$control)
          LP[[j]] = Xi[[j]] %*% fm[[j]]$coefficients + o[[j]]
          ri[[j]] = pv0[[j]]
          ri[[j]][, names(coef(fm[[j]]))] = coef(fm[[j]])
          ri[[j]][, pi[[j]]] = bi[[j]]
          d[[j]] = length(ri[[j]])
          u[[j]] = IFlm[[j]] = vector()
          IF[[j]] = matrix(rep(0, d[[j]]^2), d[[j]], d[[j]])
          r[[j]] = as.vector(ri[[j]]);
          if(fam[[j]]$family == "binomial"){
            for (k in 1:d[[j]]){
              u[[j]][k] = sum((Y - (exp(X[[j]]%*%r[[j]])/(1 + exp(X[[j]]%*%r[[j]]))))*X[[j]][,k]);
              for (l in 1:d[[j]]){
                In = (sum(X[[j]][,l]*X[[j]][,k]*exp(X[[j]]%*%r[[j]])/(1+exp(X[[j]]%*%r[[j]]))^2))
                IF[[j]][k,l] = In
              }}}else{
                break;
                print('Unknown family.');
              }
          S[[j]] = u[[j]]%*%solve(IF[[j]])%*%u[[j]]
          S[[j]] = max(S[[j]], 0)
          Sm[[j]] = sign(B[[j]] - bi[[j]]) * sqrt(S[[j]])
          mark[[j]] = 1
        }else{
          S[[j]] = 0
          Sm[[j]] = 0
          mark[[j]] = 0
        }
        Sm[[j]] = pnorm(Sm[[j]])
      }
      if (errb == "u"){
        z = sum(unlist(Sm)*unlist(w)*unlist(mark)) - alpha
      }else{
        z = sum((1 - unlist(Sm))*unlist(w)*unlist(mark)) - alpha
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
      l = multiroot(topot, start=tl, maxiter = 1000, useFortran = T,
                    ctol = 1e-5, parms = parametersl)$root
    }
    options(warn = 2)
    test = try(multiroot(topot, start = tu, maxiter = 1000, useFortran = T,
                         parms = parametersu))
    options(warn = 1)
    war = inherits(test, "try-error")
    if (war == T){
      u = uniroot(topot, lower = tl, upper = tu, maxiter = 1000, extendInt="yes",
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
