\name{mataci}
\alias{mataci}

\title{MATA Confidence Intervals}

\description{Selects models for model averaging and computes confidence intervals for parameters.}

\usage{
mataci(formula, data, nboot = 1000, family = "binomial", weights = NULL, inc = 0.5,
       cim = "Wald", ci = 0.95, par = F, ...)
}

\arguments{
  \item{formula}{An object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{A data frame (or object coercible by \code{\link[base]{as.data.frame}} to a data frame) containing the variables in the model.}
  \item{nboot}{The number of bootstrap replicates.}
  \item{family}{A description of the error distribution and link function to be used in the model. Currently, \code{family = "binomial"} is only available.}
  \item{weights}{An optional vector of ‘prior weights’ to be used in the fitting process of \code{\link[stats]{glm}}. Should be NULL or a numeric vector.}
  \item{inc}{An optional treshhold number (0\eqn{\le}inc\eqn{\le}1) or vector that splits variables into a group of fixed variables that will be included in each model that is used in model averaging, and a group of variable that will be permuted across the models.}
  \item{cim}{A specification of which method are to be used for confidence intervals estimation. The details of available methods are given under ‘Details’}
  \item{ci}{The confidence level required.}
  \item{par}{If \code{TRUE}, models will be fitted using  \code{\link[parallel]{parLapply}} function for parallel computation.}
  \item{\dots}{additional arguments to pass to a fitted \code{\link[stats]{glm}} models.}
}

\details{By default \code{mataci} estimates "\code{Wald}" type MATA confidence intervals of \insertCite{turek2012model}{MATACI}.

Method "\code{PL}" is that of \insertCite{fletcher2012model}{MATACI} which uses profile-likelihood procedure to estimate confidence intervals.

Method "\code{Score}" allows one to estimate score function based confidence intervals \insertCite{uvarov2019towards}{MATACI}.

Methods "\code{Wald-PL}" and "\code{Wald-S}" are modifyed versions of original Wald type MATA intervals that omits the restriction of the identity of the lower and upper standard errors. "\code{Wald-PL}" suggests replacing Wald standard errors by ones obtained from single-model profile-likelihood and score confidence intervals, respectively; these options can be used as an approximation to profile-likelihood and score based MATA intervals.

Method "\code{approx}" activates only \code{startpoint} function that provides asymptotic transformation-based MATA confidence intervals \insertCite{yu2014transformation}{MATACI}.
}

\value{A data.frame with point estimates, confidence intervals and the percentages of time that the variables were selected in the bootstrap process.
}
\references{
\insertRef{turek2012model}{MATACI}

\insertRef{fletcher2012model}{MATACI}

\insertRef{yu2014transformation}{MATACI}

\insertRef{fletcher2018}{MATACI}

\insertRef{uvarov2019towards}{MATACI}
}

\author{Artem Uvarov \email{auvarov@uwo.ca}}

\note{The "\code{Wald}", "\code{PL}" and "\code{approx}" options can also provide confidence intervals for liner regression model.}


\examples{
\dontrun{
## Create data set
   set.seed(616)
   X1 = rnorm(100)
   X2 = rnorm(100)
   X3 = rnorm(100)
   z = -1 - 0.5*X1 + 2*X2
   pr = 1/(1+exp(-z))
   Y = rbinom(100, 1, pr)
   data = data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)

   mataci(Y~., data = data, nboot = 1000, cim = "approx", ci = 0.95)

## In the next example 'inc' for the second variable X2 equals to 0, which means
## that this variable will always be selected in the model despite the bootstrap
## results.

   mataci(Y~X1+X2+X3, data = data, nboot = 1000, cim = "Wald", inc = c(0.4, 0, 0.5))
}}
