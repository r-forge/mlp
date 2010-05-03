\name{getIndividualPValues}
\alias{getIndividualPValues}
\title{Function to define the cutoffs to be used for each individual gene set statistic;...}
\usage{getIndividualPValues(observedGenesetStats, permutationGenesetStats)}
\description{Function to define the cutoffs to be used for each individual gene set statistic;}
\value{vector of p values for each gene set}
\arguments{\item{observedGenesetStats}{dataframe as returned by mlpStatistic}
\item{permutationGenesetStats}{output of MLP except for the geneset size column; 
matrix of gene set statistics based on the permutation procedure}}
