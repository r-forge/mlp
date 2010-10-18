\name{plotGeneSetSignificance}
\alias{plotGeneSetSignificance}
\title{Plot the Significance for the Genes of a Given Gene Set...}
\usage{plotGeneSetSignificance(geneSet, geneSetIdentifier, geneStatistic, annotationPackage, barColors)
}
\description{Plot the Significance for the Genes of a Given Gene Set}
\value{no return value}
\arguments{\item{geneSet}{object of class 'geneSetMLP' as produced by function getGeneSets}
\item{geneSetIdentifier}{identifier of the gene set for which a significance plot should be produced;
character of length one}
\item{geneStatistic}{vector of gene statistics (e.g. p values)}
\item{annotationPackage}{name of the annotation package to be used (without .db extension);
character of length one}
\item{barColors}{vector of colors to use for the bars of the barplot; defaults to NULL
in which case 'grey50' is used}
}
\examples{pathExamplePValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
pathExampleGeneSet <- system.file("exampleFiles", "exampleGeneSet.rda", package = "MLP")
pathExampleMLPResult <- system.file("exampleFiles", "exampleMLPResult.rda", package = "MLP")
load(pathExampleGeneSet)
load(pathExamplePValues)
load(pathExampleMLPResult) 
annotationPackage <- if (require(mouse4302mmentrezg)) "mouse4302mmentrezg" else "mouse4302"
geneSetID <- rownames(exampleMLPResult)[1]
dev.new(width = 10, height = 10)
op <- par(mar = c(25, 10, 6, 2))
plotGeneSetSignificance(
geneSet = exampleGeneSet, 
geneSetIdentifier = geneSetID, 
geneStatistic = examplePValues, 
annotationPackage = annotationPackage
)
par(op)}

