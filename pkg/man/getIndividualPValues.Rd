\name{getIndividualPValues}
\alias{getIndividualPValues}
\title{Function to define the cutoffs to be used for each individual gene set statistic;...}
\usage{getIndividualPValues(observedGeneSetStats, permutationGeneSetStats)}
\description{Function to define the cutoffs to be used for each individual gene set statistic;}
\value{vector of p values for each gene set}
\arguments{\item{observedGeneSetStats}{dataframe as returned by mlpStatistic}
\item{permutationGeneSetStats}{output of MLP except for the gene set size column; 
matrix of gene set statistics based on the permutation procedure}}
