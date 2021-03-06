% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{tauStarTest}
\alias{tauStarTest}
\title{Test of Independence Using the Tau* Measure}
\usage{
tauStarTest(x, y, mode = "auto", resamples = 1000)
}
\arguments{
\item{x}{a vector of sampled values.}

\item{y}{a vector of sampled values corresponding to x, y must be the same
length as x.}

\item{mode}{should be one of five possible values: "auto", "continuous",
"discrete", "mixed", or "permutation". If "auto" is selected then the
function will attempt to automatically determine whether x,y are
discrete or continuous and then perform the appropriate asymptotic
test. In cases "continuous", "discrete", and "mixed" we perform the
associated asymptotic test making the given assumption. Finally
if "permutation" is selected then the function runs a Monte-Carlo
permutation test for some given number of resamplings.}

\item{resamples}{the number of resamplings to do if mode = "permutation".
Otherwise this value is ignored.}
}
\value{
a list with class "tstest" recording the outcome of the test.
}
\description{
Performs a (consistent) test of independence between two input vectors using
the asymptotic (or permutation based) distribution of the test statistic t*.
The asymptotic results hold in the case that x is generated from either a
discrete or continous distribution and similarly for y (in particular it is
allowed for one to be continuous while the other is discrete). The asymptotic
distributions were computed in Nandy, Weihs, and Drton (2016)
<http://arxiv.org/abs/1602.04387>.
}
\examples{
set.seed(123)
x = rnorm(100)
y = rnorm(100)
testResults = tauStarTest(x,y)
print(testResults$pVal) # big p-value

y = y + x # make x and y correlated
testResults = tauStarTest(x,y)
print(testResults$pVal) # small p-value
}
\references{
Preetam Nandy, Luca Weihs, and Mathias Drton. Large-Sample Theory for the
Bergsma-Dassios Sign Covariance. arXiv preprint arXiv:1602.04387. 2016.
}
