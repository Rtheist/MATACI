\name{ScoreCI}
\alias{ScoreCI}

\title{Score based MATA Confidence Intervals}

\description{Estimates score based MATA over a set of selected models}

\usage{ScoreCI(fitted, fm, alpha = 0.05, startL, startU)}

\arguments{

  \item{fitted}{List of fitted regression models.}
  \item{fm}{Fitted full model.}
  \item{alpha}{Type I error.}
  \item{startL}{A vector containing the initial points of the lower limit to be searched.}
  \item{startU}{A vector containing the initial points of the upper limit to be searched.}
}

\value{A data.frame with confidence intervals.}
\references{
\insertRef{turek2012model}{MATACI}

\insertRef{fletcher2012model}{MATACI}

\insertRef{yu2014transformation}{MATACI}

\insertRef{fletcher2018}{MATACI}

\insertRef{uvarov2019towards}{MATACI}
}

\author{Artem Uvarov \email{auvarov@uwo.ca}}

\examples{
## Create data set
   set.seed(616)
   X1 = rnorm(1000)
   X2 = rnorm(1000)
   X3 = rnorm(1000)
   z = -1 -0.5*X1 + 2*X2
   pr = 1/(1+exp(-z))
   Y = rbinom(1000, 1, pr)
   data = data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)

## Define list of fitted models that were selected
   fitted = list()
   fitted[[1]] = fm = glm(Y ~ X1 + X2 + X3, data = data, family = "binomial")
   fitted[[2]] = glm(Y ~ X1 + X2, data = data, family = "binomial")
   fitted[[3]] = glm(Y ~ X1 + X3, data = data, family = "binomial")
   fitted[[4]] = glm(Y ~ X2 + X3, data = data, family = "binomial")
## Define vector of names - modnames = names(fm$coef)
## Note that if you define it manually you must include "(Intercept)"
   modnames = c("(Intercept)", "X1", "X2", "X3")

## Find initial points
   sp = startpoints(fitted,modnames,alpha = 0.05)
## Estimate profile-likelihood MATA interval
   ScoreCI(fitted, fm, alpha = 0.05, startL = sp[,2] , startU = sp[,3])
}
