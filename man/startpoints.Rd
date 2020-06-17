\name{startpoints}
\alias{startpoints}

\title{Initial points for faster optimization}
\description{
This function estimates confidence intervals with asymptotic version of the transformation-based model-averaged tail area (ATMATA) confidence intervals \insertCite{yu2014transformation}{MATACI}.
}
\usage{startpoints(fitted,modnames,alpha)}

\arguments{
  \item{fitted}{List of fitted regression models}
  \item{modnames}{Vector of names of fitted variables. Categorical variables must be well defined with prespecified reference category. Vector should not contain reference category.}
  \item{alpha}{\code{(1-ci)/2}}
}
\details{
modnames is extended vector of variable names with specified groups for categorical variables. If one want to fit \code{y~X1+X2+X3}, where \code{X2} is categorical variable with three categories - low, middle and high, then use \code{modnames = c("(Intercept)","X1","X2high","X2middle","X3")}
}

\value{A data.frame with point estimates and confidence intervals.}

\references{
\insertRef{turek2012model}{MATACI}

\insertRef{yu2014transformation}{MATACI}

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
   fitted[[1]] = fm = glm(Y ~ X1 + X2 + X3, data = data, family = "binomial")
   fitted[[2]] = glm(Y ~ X1 + X2, data = data, family = "binomial")
   fitted[[3]] = glm(Y ~ X1 + X3, data = data, family = "binomial")
   fitted[[4]] = glm(Y ~ X2 + X3, data = data, family = "binomial")
## Define vector of names - modnames = names(fm$coef)
## Note that if you define it manually you must include "(Intercept)"
   modnames = c("(Intercept)", "X1", "X2", "X3")

## ATMATA
   sp = startpoints(fitted, modnames, alpha = 0.05)
}



