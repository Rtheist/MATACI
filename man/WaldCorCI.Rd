\name{WaldCorCI}
\alias{WaldCorCI}

\title{Corrected Wald type MATA Confidence Intervals}

\description{Estimates corrected Wald type MATA over a set of selected models.It replaces Wald standard errors by ones obtained from single-model profile-likelihood and score confidence intervals, respectively; the estimated intervals are approximation to profile-likelihood or score based MATA intervals.}

\usage{WaldCorCI(fitted, fm, alpha = 0.05, startL, startU, CImethod = "Score")}

\arguments{
  \item{fitted}{List of fitted regression models.}

  \item{fm}{Fitted full model.}

  \item{alpha}{Type I error.}

  \item{startL}{A vector containing the initial points of the lower limit to be searched.}

  \item{startU}{A vector containing the initial points of the upper limit to be searched.}

  \item{CImethod}{A specification of which method are to be used to modify Wald type MATA confidence intervals.}
}
\details{If \code{CImethod = "Score"} then Wald standard errors will be substituted by the standard errors obtained from score confidence intervals. Otherwise, if \code{CImethod = "PL"} the profile-likelihood standard errors will be used.
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
   X1 = rnorm(100)
   X2 = rnorm(100)
   X3 = rnorm(100)
   z = -1 - 0.5*X1 + 2*X2
   pr = 1/(1+exp(-z))
   Y = rbinom(100, 1, pr)
   data = data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)

## Define list of fitted models that were selected
   fitted = list()
   fitted[[1]] = fm = glm(Y ~ X1 + X2 + X3, data=data, family = "binomial")
   fitted[[2]] = glm(Y ~ X1 + X2, data=data, family = "binomial")
   fitted[[3]] = glm(Y ~ X1 + X3, data=data, family = "binomial")
   fitted[[4]] = glm(Y ~ X2 + X3, data=data, family = "binomial")
## Define vector of names - modnames = names(fm$coef)
## Note that if you define it manually you must include "(Intercept)"
   modnames = c("(Intercept)", "X1", "X2", "X3")

## Find initial points
   sp = startpoints(fitted, modnames, alpha = 0.05)
## Estimate modifyed Wald MATA intervals
   WaldCorCI(fitted, fm, alpha = 0.05, startL = sp[,2],
             startU = sp[,3], CImethod = "Score")
}
