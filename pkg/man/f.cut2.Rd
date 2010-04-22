\name{f.cut2}
\alias{f.cut2}
\title{function to define the smoothed monotonic cutoffs leveraging across gene sets of similar size...}
\usage{f.cut2(w0, w, q.cutoff)}
\description{function to define the smoothed monotonic cutoffs leveraging across gene sets of similar size}
\value{vector of p values for each gene set}
\arguments{\item{w0}{vector of observed gene set statistics
TODO rename to observedGenesetStats}
\item{w}{matrix of gene set statistics based on the permutation procedure
TODO rename to permutationGenesetStats}
\item{q.cutoff}{vector of quantiles at which p values for each gene set are desired (to be
specified by the users)
TODO rename to criticalValuesGeneset}}
