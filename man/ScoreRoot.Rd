\name{ScoreRoot}
\alias{ScoreRoot}

\title{Score confidence intervals for single model}
\description{
This function prepares the score function for further optimization. The optimization is then done by \code{\link[rootSolve]{multiroot}}.
}
\usage{ScoreRoot(t, parms)}

\arguments{
  \item{t}{optimization parameter that presents confidence interval estimate}
  \item{parms}{list of parameters used in optimization. See Details.}
}
\details{
\code{parms} should include fitted model as '\code{fitted}', variable name '\code{var}' for which confidence limits are estimated, '\code{bound}' that indicates limit that is estimated - "l" for lower limit, "u" for upper, and '\code{ci}' that presents confidence level.
}
\value{After \code{multiroot} optimization function returns lower or upper limit of confidence interval for selected variable.}
\references{

\insertRef{rao1948large}{MATACI}

\insertRef{uvarov2019towards}{MATACI}

}
\author{Artem Uvarov \email{auvarov@uwo.ca}}

\examples{
\dontrun{
## Suppose we fit logistic model with three variables,
## where X2 has three categories - low, middle and high.

## First, fit the regression model
   fitted = glm(y ~ X1 + X2 + X3, data, family = "binomial");
## Let's estimate score confidence intervals for 'middle' category of X2.
## Define list of parameters
   paraml=list(fitted = fitted, var = "X2middle", bound = "l", ci = 0.95);
   paramu=list(fitted = fitted, var = "X2middle", bound = "u", ci = 0.95);

## Define tl and tu

## Use multiroot to get CI limits
   LL = multiroot(ScoreRoot, start = tl, maxiter = 1000, useFortran = T,
               parms = paraml)$root;
   UL = multiroot(ScoreRoot, start = tu, maxiter = 1000, useFortran = T,
               parms = paramu)$root;
   CI = c(LL, UL)
}}
