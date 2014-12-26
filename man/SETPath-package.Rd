\name{SETPath-package}
\alias{SETPath-package}
\alias{SETPath}
\docType{package}
\title{
Spiked Eigenvalue Test for Pathway data
}
\description{
Tests gene expression data from a biological pathway for biologically meaningful differences in the eigenstructure between two classes.  Specifically, it tests the null hypothesis that the two classes' leading eigenvalues and sums of eigenvalues are equal.  A pathway's leading eigenvalue arguably represents the total variability due to variability in pathway activity, while the sum of all its eigenvalues represents the variability due to pathway activity and to other, unregulated causes. 
}
\details{
\tabular{ll}{
Package: \tab SETPath\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-26\cr
License: \tab GPL-2\cr
}
The SETPath provides two functions for the implementation of the SETPath method.  The function setpath() runs the test on data from a single pathway.  The function setpath.wrapper() runs the test on a list of pathways within a dataset.  
}
\author{
Patrick Danaher

Maintainer: Patrick Danaher <patrickjdanaher at gmail.com>
}
\references{
Patrick Danaher, Debashis Paul, and Pei Wang.  "Covariance-based analyses of biological pathways."  Biometrika (2015)
}
\keyword{ package }

\examples{
#load data:
data(setpath.data)
# identify desired gene list:
genes.in.pathway = pathwaygenes[[1]]
# run test using theoretical quantiles to derive a p-value:
setpath(d1[,genes.in.pathway],d2[,genes.in.pathway],M=1,transform=NULL,verbose=TRUE,minalpha=NULL,
	normalize=TRUE,pvalue="chisq")
# now using a permutation test:
setpath(d1[,genes.in.pathway],d2[,genes.in.pathway],M=1,transform=NULL,verbose=TRUE,minalpha=NULL,
	normalize=TRUE,pvalue="permutation",npermutations=1000)
# now using the "transform" argument to test the null hypothesis that variability unrelated to the
#  first principal component (i.e. the sum of the second through final eigenvalues) is the same 
#  between classes:
setpath(d1[,genes.in.pathway],d2[,genes.in.pathway],M=1,transform=c(-1,1),verbose=TRUE,
	minalpha=NULL,normalize=TRUE,pvalue="chisq",npermutations=1000)
# now using the "transform" argument to test the compound null hypothesis that the second and third
#  eigenvalues are the same between classes:
linear.transformation = matrix(c(0,1,0,0,0,0,1,0),4)
print(linear.transformation)
setpath(d1[,genes.in.pathway],d2[,genes.in.pathway],M=3,transform=linear.transformation,
	verbose=TRUE,minalpha=NULL,normalize=TRUE,pvalue="chisq",npermutations=1000)

# use the function setpath.wrapper to analyze several pathways simultaneously
setpath.wrapper(d1,d2,pathwaygenes,pathwaynames,M=1,transform=NULL,minalpha=NULL,normalize=TRUE,
	pvalue="chisq",npermutations=10000)
}
