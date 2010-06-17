\name{mlpStatistic}
\alias{mlpStatistic}
\title{This function calculates the mean of the gene-statistic y[, 2] (the MLP statistic)...}
\usage{mlpStatistic(geneSetPValues)}
\description{This function calculates the mean of the gene-statistic y[, 2] (the MLP statistic)
for each unique gene set in geneSetPValues.}
\value{data frame with two columns; first column contains the number of genes
in the gene set; the second column contains the MLP statistic (u) for the gene set.}
\arguments{\item{geneSetPValues}{is a list of gene sets; each component of the list (i.e. each gene
set) contains a vector of p values for the genes in that gene set
respectively}}
