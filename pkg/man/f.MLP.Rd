\name{f.MLP}
\alias{f.MLP}
\title{This function calculates p-values for each gene set based on row permutations...}
\usage{f.MLP(geneSet, genePValue, minGenes=5, maxGenes=100,
    rowPermutations=TRUE, nPermutations=10, smoothPValues=TRUE)}
\description{This function calculates p-values for each gene set based on row permutations
of the gene p values or column permutations of the expression matrix; the p values
can be obtained either as individual gene set p values or p values based on smoothing
across gene sets of similar size.}
\value{data frame with three columns: genesetSize, genesetStatistic and genesetPValue.}
\references{TODO}
\arguments{\item{geneSet}{is the input table (two-column numeric matrix) of gene-sets and gene-ids. (GO category and Gene ID)
TODO matrix --> should become a list: numeric -> character matrix (-> list later on)}
\item{genePValue}{is the input table (two-column numeric matrix) of gene-ids and pvalues (or other summary statistic).
TODO gene-ids consistent with gene ID's in x (i.e. character vector)}
\item{minGenes}{minimum number of genes in a gene set for it to be considered (lower threshold
for gene set size)}
\item{maxGenes}{maximum number of genes in a gene set for it to be considered (upper threshold
for gene set size)}
\item{rowPermutations}{logical indicating whether to use row permutations (TRUE; default) or column 
permutations (FALSE)}
\item{nPermutations}{is the number of simulations.}
\item{smoothPValues}{logical indicating whether one wants to calculate smoothed cut-off thresholds (TRUE; default)
or not (FALSE).}}
