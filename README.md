# MATACI

Frequentist Model-Averaged Tail Area Confidence Intervals Estimation

## Description

This package allows to estimate valid confidence intervals for logistic regression coefficients estimated from model averaging procedure. The model averaging can be done over prespecified set of models or the set can be automatically constructed by bootstrap procedure with specified threshold.

## Getting Started

### Dependencies

The package needs **parallel** library for parallel computing of set of models before averaging process. The **tictoc** package is used to report the time that each part of algorithm in the main function **mataci** took. The **rootSolve** package is used for confidence intervals estimation.

### Installing

* Install as you would normally install a package from GitHub repository.

```{r, eval = FALSE}
	  # install.packages("devtools")
	  devtools::install_github("Rtheist/MATACI")
```

* Check that needed packages are installed and activated.

* If R Documentation is not loading after installation of the package, try to restart R.

### Executing program

* The main function **mataci**, allows to select the set of models by bootstrap and estimate confidence intervals with desired method.

```{r, eval = FALSE}
	  set.seed(616)
	  X1 = rnorm(100)
	  X2 = rnorm(100)
	  X3 = rnorm(100)
	  z = -1 - 0.5*X1 + 2*X2
	  pr = 1/(1+exp(-z))
	  Y = rbinom(100, 1, pr)
	  data = data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)
	  mataci(Y~., data = data, nboot = 1000, cim = "approx", ci = 0.95)
```

* If the models are already selected then the confidence intervals can be obtained by implementing one of the following commands. First, simulate the data and set of models

```{r, eval = FALSE}
	# Create data set
	  set.seed(616)
	  X1 = rnorm(100)
	  X2 = rnorm(100)
	  X3 = rnorm(100)
	  z = -1 - 0.5*X1 + 2*X2
	  pr = 1/(1+exp(-z))
	  Y = rbinom(100, 1, pr)
	  data = data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)

	# Define list of fitted models that were selected
	  fitted = list()
	  fitted[[1]] = fm = glm(Y ~ X1 + X2 + X3, data = data, family = "binomial")
	  fitted[[2]] = glm(Y ~ X1 + X2, data = data, family = "binomial")
	  fitted[[3]] = glm(Y ~ X1 + X3, data = data, family = "binomial")
	  fitted[[4]] = glm(Y ~ X2 + X3, data = data, family = "binomial")
	# Define vector of names - modnames = names(fm$coef)
	# Note that if you define it manually you must include "(Intercept)"
	  modnames = c("(Intercept)", "X1", "X2", "X3")
```

* For asymptotic version of the transformation-based model-averaged tail area (ATMATA) confidence intervals ([Yu et al. 2014](https://link.springer.com/article/10.1007/s00180-014-0514-1)):

```{r, eval = FALSE}
	# ATMATA confidence intervals
	   sp = startpoints(fitted, modnames, alpha = 0.05)
```

* For Wald-type MATA confidence intervals ([Turek and Fletcher 2012](https://www.sciencedirect.com/science/article/abs/pii/S0167947312001144)):

```{r, eval = FALSE}
	# Find initial points
	   sp = startpoints(fitted, modnames, alpha = 0.05)
	# Estimate Wald MATA intervals
	   WaldCI(fitted, fm, alpha = 0.05, startL = sp[,2] , startU = sp[,3])
```
* For profile-likelihood based MATA confidence intervals ([Fletcher and Turek 2012](https://link.springer.com/article/10.1007/s13253-011-0064-8)):

```{r, eval = FALSE}
	# Find initial points
	   sp = startpoints(fitted, modnames, alpha = 0.05)
	# Estimate profile-likelihood MATA intervals
	   ProfLCI(fitted, fm, alpha = 0.05, startL = sp[,2] , startU = sp[,3])
```


* For score function based MATA confidence intervals ([Uvarov 2019](https://ir.lib.uwo.ca/cgi/viewcontent.cgi?article=8759&context=etd)):

```{r, eval = FALSE}
	# Find initial points
	   sp = startpoints(fitted, modnames, alpha = 0.05)
	# Estimate score function based MATA intervals
	   ScoreCI(fitted, fm, alpha = 0.05, startL = sp[,2] , startU = sp[,3])
```

* For modifyed Wald-type MATA confidence intervals ([Uvarov 2019](https://ir.lib.uwo.ca/cgi/viewcontent.cgi?article=8759&context=etd)):

```{r, eval = FALSE}
	# Find initial points
	   sp = startpoints(fitted, modnames, alpha = 0.05)
	# Estimate modifyed Wald MATA intervals
	   WaldCorCI(fitted, fm, alpha = 0.05, startL = sp[,2], startU = sp[,3], CImethod = "Score")
```

## Authors

[Artem Uvarov](mailto:auvarov@uwo.ca)
