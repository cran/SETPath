\name{unbias.eigens}
\alias{unbias.eigens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Unbiased estimation of leading eigenvalues
}
\description{
A function called by setpath(), used to attain unbiased estimates of leading eigenvalues and to estimate the leading eigenvalues under the null hypothesis.}
\usage{
unbias.eigens(L, g, w, minalpha = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
A vector of length two containing a leading eigenvalue from each of the two datasets.}
  \item{g}{
A vector of length two containing the gamma (p/n) values from the two datasets.}
  \item{w}{
The weights to assign to the two classes when estimating common leading eigenvalues under the null hypothesis.
}
  \item{minalpha}{
Can be used to tweak the estimation of the leading eigenvalues under the null hypothesis.  If NULL, eigenvalues are truncated below at 1 + 2*sqrt(gamma).  Otherwise, eigenvalues are truncated below at 1 + sqrt(gamma) + minalpha.
}
}
\details{
Called by the setpath() function, not useful on its own.
}
\value{
\item{QLcorrection}{The correction factor to remove the bias in the difference between the two eigenvalues.}
\item{a0}{The common eigenvalue estimated under the null hypothesis}
\item{a}{The unbiased eigenvalue estimates from each class}
%%    return(list(QLcorrection = QLcorrection, a0 = a0, a = a))
}
\references{
Patrick Danaher, Debashis Paul, and Pei Wang.  "Covariance-based analyses of biological pathways."  Biometrika (2015)
}
\author{
Patrick Danaher
}
\examples{

## The function is currently defined as
function (L, g, w, minalpha = NULL) 
{
    if (length(minalpha) == 0) {
        minalpha = sqrt(max(g))
    }
    a = c()
    for (k in 1:2) {
        if (L[k] < (1 + sqrt(g[k]))^2) {
            a[k] = 1 + sqrt(max(g)) + minalpha
        }
        if (L[k] >= (1 + sqrt(g[k]))^2) {
            a[k] = ((1 + L[k] - g[k]) + sqrt((1 + L[k] - g[k])^2 - 
                4 * L[k]))/2
        }
    }
    a0 = sum(a * w)
    a0 = max(c(a0, 1 + sqrt(g)))
    QLcorrection = (g[1] - g[2]) * a0/(a0 - 1)
    return(list(QLcorrection = QLcorrection, a0 = a0, a = a))
  }
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
