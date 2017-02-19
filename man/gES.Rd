\name{ges.mcmc}
\alias{ges.mcmc}
\title{MCMC algorithm for graphical exponential screening estimator}
\description{
Function of MCMC algorithm to compute graphical exponential screening estimator.
}
\usage{
ges.mcmc(x1, x2, hard.thres = 0.05, num.burn.in = 4000, num.estimation = 4000, plot.mcmc = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
  An \code{n} by \code{p} data matrix used for constructing individual graph estimators, where \code{n} is the sample size and \code{p} is the dimension.
  }
  \item{x2}{
  An \code{m} by \code{p} data matrix used for aggregating individual estimators, where \code{m} is the sample size and \code{p} is the dimension.
  }
  \item{hard.thres}{
  Hard thresholding parameter for covariance matrix estimation.
  }
  \item{num.burn.in}{
  Number of burn-in iterations in MCMC algorithm.
  }
  \item{num.estimation}{
  Number of MCMC iterations for precision matrix estimation.
  }
  \item{plot.mcmc}{
  If \code{TRUE}, plot the number of selected edges as a function of MCMC steps.
  }
}
\value{
Estimated precision matrix via graphical exponential screening.
}
\examples{
library(huge)

n <- 400
p <- 50
model <- huge.generator(n = n, d = p, graph = "hub", g = 5, v = 0.3)
x <- model$data
est.ges <- ges.mcmc(x[1:(n / 2), ], x[(n / 2 + 1):n, ])
}
