\name{setpath}
\alias{setpath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs the Spiked Eigenvalue Test for Pathway data (SETPath) on data from a single pathway
}
\description{
Compares a pathway's gene expression data from two classes, testing the null hypothesis that the first eigenvalue and the sum of the eigenvalues are equal between classes.  
}
\usage{
setpath(d1, d2, M = 1, transform = NULL, verbose = FALSE, minalpha = NULL, 
	normalize = TRUE, pvalue = "chisq", npermutations = 10000)
}
\arguments{
  \item{d1}{
A n1*p matrix from the pathway's data in the first class}
  \item{d2}{
A n2*p matrix from the pathway's data in the second class, with the same genes (column names) as d1}
  \item{M}{
The null hypothesis specifies that the first M eigenvalues are equal between the two classes.  SETPath was conceived as a test of just the first eigenvalue and the sum of the eigenvalues (M=1), but cases where multiple leading eigenvalues are of biological interest may indicate M>1.
}
  \item{transform}{
Default NULL.  Otherwise, a matrix or vector specifying a linear transformation of the null hypothesis.  The use of the transform argument modifies the standard SETPath test of equality of the first eigenvalue and the sum of eigenvalues between classes.  Specifically, if the null hypothesis is that (alpha0.1, ..., alpha0.M,T0) = (alpha1.1, ..., alpha1.M,T1), the transform argument changes the null hypothesis to t(transform)\%*\%(alpha0.1, ..., alpha0.M,T0) = t(transform)\%*\%(alpha1.1, ..., alpha1.M,T1).
}
  \item{verbose}{
Indicates whether to return more than the p-value.
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
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{pval}{The p-value from SETPath}
If verbose==TRUE, the additional values are output:
\item{a0}{Estimates of the leading M eigenvalues, assuming the null hypothesis is true.}
\item{correction}{The correction factors applied to the between-class differences between the M leading eigenvalues.  (Correction factors are needed to account for the sample size-dependent bias in empirical eigenvalues.)}
\item{covQ}{An (M+1)*(M+1) matrix giving the estimated covariance of the between-class differences in the leading M eigenvalues and the sum of the eigenvalues.}
\item{m}{M, the number of leading eigenvalues included in the null hypothesis.}
}
\references{
Patrick Danaher, Debashis Paul, and Pei Wang.  "Covariance-based analyses of biological pathways."  Biometrika (2015)
}
\author{
Patrick Danaher}


%% ~Make other sections like Warning with \section{Warning }{....} ~


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



## The function is currently defined as
function (d1, d2, M = 1, transform = NULL, verbose = FALSE, minalpha = NULL, 
    normalize = TRUE, pvalue = "chisq", npermutations = 10000) 
{
    p = dim(d1)[2]
    n1 = dim(d1)[1]
    n2 = dim(d2)[1]
    if (normalize) {
      e1 = eigen(cov(d1),symmetric=TRUE,only.values=TRUE)$values
      e2 = eigen(cov(d2),symmetric=TRUE,only.values=TRUE)$values
      if(n1>p){medianeigen1 = median(e1)}
      if(n2>p){medianeigen2 = median(e2)}
      if(n1<=p){medianeigen1 = median(e1[e1>1e-12])*n1/p}
      if(n2<=p){medianeigen2 = median(e2[e2>1e-12])*n2/p}
      scaling.factor = mean(medianeigen1,medianeigen2)
      d1 = d1/sqrt(scaling.factor)
      d2 = d2/sqrt(scaling.factor)
    }
    p = dim(d1)[2]
    n1 = dim(d1)[1]
    n2 = dim(d2)[1]
    n = c(n1, n2)
    w = n/sum(n)
    d = list(d1, d2)
    g1 = p/n1
    g2 = p/n2
    g = c(g1, g2)
    covs = list()
    covs[[1]] = cov(d1)
    covs[[2]] = cov(d2)
    sighat = list(cov(d1), cov(d2))
    e1 = eigen(sighat[[1]], symmetric = TRUE, only.values = TRUE)$values
    e2 = eigen(sighat[[2]], symmetric = TRUE, only.values = TRUE)$values
    L = cbind(e1[1:M], e2[1:M])
    T1 = sum(e1)
    T2 = sum(e2)
    T = c(T1, T2)
    alphabar = matrix(NA, M, 2)
    a0 = QLcorrection = c()
    for (m in 1:M) {
        eigencorrect = unbias.eigens(L[m, ], g, w, minalpha)
        alphabar[m, ] = eigencorrect$a
        a0[m] = eigencorrect$a0
        QLcorrection[m] = eigencorrect$QLcorrection
    }
    thresh = (1 + sqrt(g))^2 + sqrt(2 * log(n)/n)
    mhat = c()
    mhat[1] = max(sum(e1 > thresh[1]), M)
    mhat[2] = max(sum(e2 > thresh[2]), M)
    spikes = list()
    for (k in 1:2) {
        spikes[[k]] = rep(NA, mhat[k])
    }
    for (k in 1:k) {
        for (m in 1:mhat[k]) {
            tempeigencorrect = unbias.eigens(c(e1[m], e2[m]), 
                g, w, minalpha)
            spikes[[k]][m] = tempeigencorrect$a[k]
        }
    }
    covQ = matrix(0, M + 1, M + 1)
    varT = c()
    for (k in 1:2) {
        varT[k] = 2 * (sum(a0^2)/n[k] + (p - M)/n[k])
        if (mhat[k] > M) {
            varT[k] = varT[k] + 2/n[k] * (sum((spikes[[k]][(M + 
                1):mhat[k]])^2) - (mhat[k] - M))
        }
    }
    covQ[M + 1, M + 1] = sum(varT)
    for (m in 1:M) {
        rho = theta = varLs = c()
        c0 = (1/g[1] + 1/g[2])^2 * (a0[m] - 1)^2
        for (k in 1:2) {
            rho[k] = a0[m] * (1 + g[k]/(a0[m] - 1))
            deriv.f.k = 0.5 * (1 + (rho[k] - 1 - g[k])/sqrt((rho[k] - 
                1 - g[k])^2 - 4 * g[k]))
            theta[k] = 1 + (g[1] - g[2])/c0 * deriv.f.k/g[k]
            varLs[k] = 2 * a0[m]/n[k] * theta[k]^2 * rho[k]/(1 + 
                a0[m] * g[k]/((a0[m] - 1)^2 - g[k]))
        }
        covQ[m, m] = sum(varLs)
    }
    for (m in 1:M) {
        rho = theta = covLTs = c()
        c0 = (1/g[1] + 1/g[2])^2 * (a0[m] - 1)^2
        for (k in 1:2) {
            rho[k] = a0[m] * (1 + g[k]/(a0[m] - 1))
            deriv.f.k = 0.5 * (1 + (rho[k] - 1 - g[k])/sqrt((rho[k] - 
                1 - g[k])^2 - 4 * g[k]))
            theta[k] = 1 + (g[1] - g[2])/c0 * deriv.f.k/g[k]
            covLTs[k] = 2 * a0[m]/n[k] * theta[k] * rho[k]/(1 + 
                a0[m] * g[k]/((a0[m] - 1)^2 - g[k]))
        }
        covLT = sum(covLTs)
        covQ[m, M + 1] = covQ[M + 1, m] = covLT
    }
    Q = c(L[, 1] - L[, 2] - QLcorrection, T1 - T2)
    A = transform
    if (length(transform) == 0) {
        A = diag(M + 1)
    }
    A=as.matrix(A)
    if(dim(A)[1]!=M+1){stop("Dimension of linear transformation does not match
					   the dimension of the null hypothesis.")}
    stat = t(Q) \%*\% A \%*\% solve(t(A) \%*\% covQ \%*\% A) \%*\% t(A) \%*\% 
        Q
    out = list()
    out$stat = stat
    if (pvalue == "chisq") {
        out$pval = 1 - pchisq(stat, dim(A)[2])
    }
    if (pvalue == "permutation") {
        d.combined = rbind(d1, d2)
        permstats = c()
        for (i in 1:npermutations) {
            prows1 = sample(1:dim(d.combined)[1], dim(d1)[1], 
                replace = FALSE)
            prows2 = setdiff(1:dim(d.combined)[1], prows1)
            permstats[i] = setpath(d1 = d.combined[prows1, ], d2 = d.combined[prows2, 
                ], M = M, transform = transform, verbose = FALSE, 
                minalpha = minalpha, normalize = FALSE, pvalue = "chisq")$stat
        }
        out$pval = mean(as.vector(stat) < permstats)
    }
    if (verbose) {
        out$stats = rbind(L, c(T1, T2))
        out$a0 = a0
        out$correction = QLcorrection
        out$covQ = covQ
        out$m = m
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ setpath }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
