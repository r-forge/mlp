\name{pperm}
\alias{pperm}
\title{The function computes a matrix of gene set p values.}
\usage{pperm(x, p0, p1, ind.exact=TRUE, np=0)}
\description{The function computes a matrix of gene set p values. Each column in the matrix is a set of p values (for the rows of
the matrix) corresponding to a permutation of the columns of the original expression matrix.
The original permutation is identified and removed from the simulated data.}
\value{matrix of p values corresponding to a permutation of the columns of the expression matrix}
\arguments{\item{x}{is the expression data for which the original and column-permuted p-values are calculated.}
\item{p0}{vector of column indices (of x) indicating the control columns.}
\item{p1}{vector of column indices (of x) indicating the treatment columns.}
\item{ind.exact}{logical indicating whether exact permutations should be used (TRUE; default) or
whether random permutations should be used (FALSE)}
\item{np}{is the number of random permutations, if ind.exact != TRUE.
TODO turn into real (executable) example
Example: ypp <- pperm(dma.norm, 1:4, 5:8, ind.exact = TRUE, 1))
Example: ypp <- pperm(dma.norm, 1:4, 5:8, ind.exact = FALSE, 100))}}
