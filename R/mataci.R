#' @importFrom Rdpack reprompt
#' @import parallel
#' @import tictoc
#' @import rootSolve
#' @importFrom stats as.formula coef confint deviance
#' @importFrom stats glm glm.fit logLik model.frame model.matrix
#' @importFrom stats model.offset model.response model.weights pnorm pt
#' @importFrom stats qchisq qf qnorm step terms uniroot
#' @export
mataci = function(formula, data, nboot = 1000, family = "binomial", weights = NULL, inc=0.5, cim="Wald", ci=0.95, par=F,...){
  '%!in%' = Negate('%in%');
  if (cim %!in% c("approx", "Wald", "Score", "PL", "Wald-S", "Wald-PL")){
    stop("Unknown method.");
  }
  trms = terms(formula, data = data);
  variables = attr(trms, "term.labels");
  out = all.vars(formula)[1];
  formula = as.formula(paste(paste(out, "~"), paste(variables, collapse = "+")))
  n = dim(data)[1];
  if(is.null(weights)){weights = as.vector((rep(1,n)))}
  data$weights = weights
  Full.Model = glm(formula = formula, family = family, data = data, weights = weights,...)
  rownames(data) = NULL;
  data = data[,c(variables, out)];
  alpha = (1-ci)/2;
  forname = c("(Intercept)", variables);
  tic("Total")
  tic("Bootstrapping")
  if(length(inc) != 1 & length(inc) != length(variables)){
    stop("Inc has to be a single value or vector of length equal to number of considered variables.")}
  if(length(inc) == 1) inc = rep(inc, length(variables))
  cl = makeCluster(detectCores());
  clusterExport(cl, c("data", "n", "formula", "variables", "weights", "family"), envir = environment());
  repl1 = parLapply(cl = cl, 1:nboot, function(i, dataA = data, smpl = n, ...){
    #Resampling
    dataA$weights = weights
    dataB = data.frame(dataA[sample(nrow(dataA), size = smpl, replace = TRUE), ]);
    mus = glm(formula, family = family, weights = weights, data = dataB,...);
    #Stepwise selection
    Sl.us = step(mus, direction = "backward", trace = F, k = 2);
    options(warn = 2);
    test = try(glm(Sl.us$formula, family = family, weights = weights,data = data.frame(dataB),...),...);
    options(warn = 1);
    war = inherits(test, "try-error");
    Slstep.us = summary(Sl.us);
    selected = attr(terms(Slstep.us), "term.labels");
    estimates = names(Sl.us$coefficients)
    Est = rep(1, length(estimates))
    names(Est) = estimates
    estimates.all = names(mus$coefficients)
    Sel = rep(1, length(selected))
    names(Sel) = selected
    ZC0us = setdiff(attr(terms(summary(mus)), "term.labels"), names(Sel));
    NonEst = setdiff(names(mus$coef), estimates);
    ZCzero.us = rep(0, length(ZC0us));
    NonEst.us = rep(0, length(NonEst));
    names(ZCzero.us) = ZC0us;
    names(NonEst.us) = NonEst;
    indication = c(ZCzero.us, Sel);
    Prob = c(NonEst.us, Est);
    indication = indication[attr(terms(summary(mus)), "term.labels")]
    Prob = Prob[names(mus$coef)]
    return(list(indication = indication, war = war, Prob = Prob))})
  stopCluster(cl);
  p = length(repl1[[1]]$indication);
  pProb = length(repl1[[1]]$Prob);
  beta.Step.zero = t(matrix(unlist(lapply(repl1, "[[", "indication"), use.names = T),nrow = p));
  pProb = t(matrix(unlist(lapply(repl1, "[[", "Prob"), use.names = T), nrow = pProb));
  pProb = apply(pProb,2,mean)*100
  war = as.numeric(t(matrix(unlist(lapply(repl1, "[[", "war")), nrow = 1)));
  toc()
  tic("Model averaging")
  message("The proportion of misconvergence is ", round(mean(war),2), "%");
  b = data.frame(cbind(beta.Step.zero, war));
  beta.Step.zero = b[which(b$war == 0), ];
  beta.Step.zero = beta.Step.zero[, 1:(p)];
  freq = apply(beta.Step.zero, 2, mean);
  names(freq)=variables
  Prob = round(freq * 100, 2);
  inc.frac = names(freq)[which(freq >= inc)];
  if (length(inc.frac) == length(forname)){
    minInc = names(freq)[which(freq == min(freq))];
    inc.frac = setdiff(forname, minInc)}
  exc.frac = names(freq)[which(freq < inc)];

  lst = rep(list(c(T, F)), length(exc.frac))
  regMat = expand.grid(lst);
  names(regMat) = exc.frac
  MT = data.frame(matrix(rep(TRUE, dim(regMat)[1]*length(inc.frac)), nrow = dim(regMat)[1], ncol = length(inc.frac)))
  names(MT) = inc.frac
  MGrid = data.frame(regMat,MT)
  MGrid = MGrid[variables]
  allModelsList = apply(MGrid, 1, function(x) as.formula(paste(paste(c(out, "~"), collapse=""), paste(variables[x], collapse="+"), collapse="")))
  if (par == F){
    atr.inc = lapply(allModelsList, function(x, data) glm(x, data = data.frame(data), family = family, weights = weights,...), data = data);
  }else{
    MAdata = data.frame(data)
    cl = makeCluster(detectCores());
    clusterExport(cl,c("MAdata", "allModelsList", "family", "weights"), envir = environment());

    atr.inc = parLapply(cl = cl, allModelsList, function(x, data) glm(x, data = data.frame(MAdata), family = family, weights = weights,...), data = MAdata);
    stopCluster(cl)
  }
  message("The model averaging is done over ", length(allModelsList), " models");

  ######################
  #MATA after INCLUSION#
  ######################
  modnames = names(Full.Model$coefficients)
  start.Inc = data.frame(startpoints(fitted = atr.inc, modnames = modnames, alpha = alpha));
  if (cim != "approx"){
    if(cim == "Wald"){
      confi = WaldCI(fitted = atr.inc, fm = Full.Model, alpha = alpha,
                    startL = start.Inc[, 2], startU = start.Inc[, 3]);
    }else if(cim == "Score"){
      confi = ScoreCI(fitted = atr.inc, fm = Full.Model, alpha = alpha,
                     startL = start.Inc[, 2], startU = start.Inc[, 3]);
    }else if(cim == "PL"){
      confi = ProfLCI(fitted = atr.inc, fm = Full.Model, alpha = alpha,
                     startL = start.Inc[, 2], startU = start.Inc[, 3]);
    }else if(cim == "Wald-S"){
      confi = WaldCorCI(fitted = atr.inc, fm = Full.Model, alpha = alpha,
                       startL = start.Inc[, 2], startU = start.Inc[, 3], CImethod = "Score");
    }else if(cim == "Wald-PL"){
      confi = WaldCorCI(fitted = atr.inc, fm = Full.Model, alpha = alpha,
                       startL = start.Inc[, 2], startU = start.Inc[, 3], CImethod = "PL");
    }
    results = data.frame(Estimates = start.Inc[, 1], CI.lower = confi[, 1],CI.upper = confi[, 2], Prop = pProb);
  }else{
    results = data.frame(Estimates = start.Inc[, 1], CI.lower = start.Inc[, 2],CI.upper = start.Inc[, 3], Prop = pProb);
    rownames(results) = rownames(start.Inc)
    }
  toc()
  toc()
  colnames(results) = c("Estimates", paste(round(100*alpha, 2), "%",
                                           sep=""),  paste(round(100*(1-alpha), 2), "%", sep=""),"Prop")
  return(results)
}
