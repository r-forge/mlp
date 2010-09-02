\name{plotGeneSetSignificance}
\alias{plotGeneSetSignificance}
\title{Plot the Significance for the Genes of a Given Gene Set...}
\usage{plotGeneSetSignificance(geneSet, geneSetIdentifier, geneStatistic,
    annotationPackage, barColors)}
\description{Plot the Significance for the Genes of a Given Gene Set}
\value{no return value}
\arguments{\item{geneSet}{object of class 'geneSetMLP' as produced by function getGeneSets}
\item{geneSetIdentifier}{identifier of the gene set for which a significance plot should be produced;
character of length one}
\item{geneStatistic}{vector of gene statistics (e.g. p values)}
\item{annotationPackage}{name of the annotation package to be used (without .db extension);
character of length one}
\item{barColors}{vector of colors to use for the bars of the barplot; defaults to NULL
in which case 'grey50' is used}}
