\name{getSmoothedPValues}
\alias{getSmoothedPValues}
\title{Function to define the smoothed monotonic cutoffs leveraging across gene sets of similar size;...}
\usage{getSmoothedPValues(observedGenesetStats, permutationGenesetStats,
    q.cutoff, df=9)}
\description{Function to define the smoothed monotonic cutoffs leveraging across gene sets of similar size;
it calls ee1 and ctpval1 to create smoothed quantile curves and calculate the corresponding
p values for the geneset.}
\value{vector of p values for each gene set}
\arguments{\item{observedGenesetStats}{dataframe as returned by mlpStatistic}
\item{permutationGenesetStats}{output of MLP except for the geneset size column; 
matrix of gene set statistics based on the permutation procedure}
\item{q.cutoff}{vector of quantiles at which p values for each gene set are desired (to be
specified by the users)}
\item{df}{degrees of freedom used for the smoothing splines used to calculate the smoothed quantile curves; defaults to 9.}}
