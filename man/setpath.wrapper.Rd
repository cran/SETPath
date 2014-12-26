\name{setpath.wrapper}
\alias{setpath.wrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs the Spiked Eigenvalue Test for Pathway data (SETPath) on multiple pathways in a dataset
}
\description{
For pathway in a list, compares the pathway's gene expression data from two classes, testing the null hypothesis that the first eigenvalue and the sum of the eigenvalues are equal between classes.  Returns a matrix of results.
}
\usage{
setpath.wrapper(d1, d2, pathwaygenes, pathwaynames, M = 1, transform = NULL,
	minalpha = NULL, normalize = TRUE, pvalue = "chisq", npermutations = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d1}{
A n1*p matrix from the pathway's data in the first class}
  \item{d2}{
A n2*p matrix from the pathway's data in the second class, with the same genes (column names) as d1}
  \item{pathwaygenes}{
A list of vectors containing the names of the genes in each of the pathways}
  \item{pathwaynames}{
A vector of the names of the pathways
}
  \item{M}{
The null hypothesis specifies that the first M eigenvalues are equal between the two classes.  SETPath was conceived as a test of just the first eigenvalue and the sum of the eigenvalues (M=1), but cases where multiple leading eigenvalues are of biological interest may indicate M>1.
}
  \item{transform}{
Default NULL.  Otherwise, a matrix or vector specifying a linear transformation of the null hypothesis.  The use of the transform argument modifies the standard SETPath test of equality of the first eigenvalue and the sum of eigenvalues between classes.  Specifically, if the null hypothesis is that (alpha0.1, ..., alpha0.M,T0) = (alpha1.1, ..., alpha1.M,T1), the transform argument changes the null hypothesis to t(transform)\%*\%(alpha0.1, ..., alpha0.M,T0) = t(transform)\%*\%(alpha1.1, ..., alpha1.M,T1).
}
  \item{minalpha}{
Can be used to tweak the estimation of the leading eigenvalues under the null hypothesis.  If NULL, eigenvalues are truncated below at 1 + 2*sqrt(gamma).  Otherwise, eigenvalues are truncated below at 1 + sqrt(gamma) + minalpha.}
  \item{normalize}{
SETPath assumes that the unspiked eigenvalues of each class's dataset are equal to 1.  If the normalize argument is set to TRUE, the data will be rescaled by the average of the median of the non-zero eigenvalues of the two datasets.  A further adjustment is performed when n<p.   
}
  \item{pvalue}{
If pvalue=="chisq", the theoretical distribution of the test statistic will be used to compute a p-value.  Alternatively, use pvalue=="permutation".
}
  \item{npermutations}{
If pvalue == "permutation", the number of permutations.
}
}
\value{
A matrix of results, with a row for each pathway and (2*(M+2)) columns.  For each pathway, the following is reported: the number of genes in the pathway, the leading M eigenvalues in each class, the sum of the eigenvalues in each class, and the p-value for the test.}
\references{
Patrick Danaher, Debashis Paul, and Pei Wang.  "Covariance-based analyses of biological pathways."  Biometrika (2015)
}
\author{
Patrick Danaher
}

\examples{
# use the function setpath.wrapper to analyze several pathways simultaneously
data(setpath.data)
setpath.wrapper(d1,d2,pathwaygenes,pathwaynames,M=1,transform=NULL,minalpha=NULL,normalize=TRUE,
	pvalue="chisq",npermutations=10000)

## The function is currently defined as
function (d1, d2, pathwaygenes, pathwaynames, M = 1, transform = NULL, 
    minalpha = NULL, normalize = TRUE, pvalue = "chisq", npermutations = 10000) 
{
    K = length(pathwaynames)
    results = matrix(NA, K, 2 * (M + 1) + 2)
    dimnames(results)[[1]] = pathwaynames
    dimnames(results)[[2]] = c("n.genes", paste("alpha.0", 1:M, 
        sep = "."), "T.0", paste("alpha.0", 1:M, sep = "."), 
        "T.0", "pval")
    if (!identical(dimnames(d1)[[2]], dimnames(d2)[[2]])) {
        stop("d1 and d2 have different feature (column) names.")
    }
    for (k in 1:K) {
        missinggenes = setdiff(pathwaygenes[[k]], dimnames(d1)[[2]])
        if (length(missinggenes) > 0) {
            warning(c("The following pathway genes are missing from the dataset:", 
                missinggenes))
            pathwaygenes[[k]] = intersect(pathwaygenes[[k]], 
                dimnames(d1)[[2]])
        }
        temp = setpath(d1[, pathwaygenes[[k]]], d2[, pathwaygenes[[k]]], 
            M = M, transform = transform, verbose = TRUE, minalpha = minalpha, 
            normalize = normalize, pvalue = pvalue, npermutations = npermutations)
        results[k, ] = c(length(pathwaygenes[[k]]), temp$stats[, 
            1], temp$stats[, 2], temp$pval)
    }
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
