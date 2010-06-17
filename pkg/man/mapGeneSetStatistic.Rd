\name{mapGeneSetStatistic}
\alias{mapGeneSetStatistic}
\title{Function that will replace all gene identifiers in a gene set with the p values associated to the gene identifiers...}
\usage{mapGeneSetStatistic(geneSet, geneStatistic)}
\description{Function that will replace all gene identifiers in a gene set with the p values associated to the gene identifiers}
\value{list of gene sets; each component of the list (gene set) is a numeric vector of p values; the names
of the numeric vector are the gene identifiers corresponding to the p values}
\arguments{\item{geneSet}{list of gene sets; each gene set (component of the list) is a character vector of gene identifiers}
\item{geneStatistic}{named vector; each name is a gene identifier and the corresponding element contains the p value
corresponding to its gene identifier}}
