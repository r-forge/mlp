\name{MLP}
\alias{MLP}
\title{This function calculates p-values for each gene set based on row permutations...}
\usage{MLP(geneSet, geneStatistic, minGenes=5, maxGenes=100, rowPermutations=TRUE,
    nPermutations=100, smoothPValues=TRUE, probabilityVector=c(0.5,
    0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999), df=9,
    addGeneSetDescription=TRUE)
}
\description{This function calculates p-values for each gene set based on row permutations
of the gene p values or column permutations of the expression matrix; the p values
can be obtained either as individual gene set p values or p values based on smoothing
across gene sets of similar size.}
\value{data frame with four (or five) columns: totalGeneSetSize, testedGeneSetSize, geneSetStatistic and geneSetPValue
and (if addDescription is set to TRUE) geneSetDescription; the rows of the data frame are ordered by 
ascending geneSetPValue.}
\references{Raghavan, Nandini et al. (2007). The high-level similarity of some disparate gene expression measures,
Bioinformatics, 23, 22, 3032-3038.}
\arguments{\item{geneSet}{is the input list of gene sets (components) and gene IDs (character vectors). 
A gene set can, for example, be a GO category with for each category Entrez gene identifiers; 
The \link{getGeneSets} function can be used to construct the geneSet argument for different
pathway sources.}
\item{geneStatistic}{is either a named numeric vector (if rowPermutations is TRUE)
or a numeric matrix of pvalues (if rowPermutations is FALSE). The names of the numeric vector
or row names of the matrix should represent the gene IDs.}
\item{minGenes}{minimum number of genes in a gene set for it to be considered (lower threshold
for gene set size)}
\item{maxGenes}{maximum number of genes in a gene set for it to be considered (upper threshold
for gene set size)}
\item{rowPermutations}{logical indicating whether to use row permutations (TRUE; default) or column 
permutations (FALSE)}
\item{nPermutations}{is the number of simulations. By default 100 permutations are conducted.}
\item{smoothPValues}{logical indicating whether one wants to calculate smoothed cut-off thresholds (TRUE; default)
or not (FALSE).}
\item{probabilityVector}{vector of quantiles at which p values for each gene set are desired}
\item{df}{degrees of freedom for the smooth.spline function used in getSmoothedPValues}
\item{addGeneSetDescription}{logical indicating whether a column with the gene set description be added to
the output data frame; defaults to TRUE.}
}
\examples{if (require(GO.db)){
pathExampleGeneSet <- system.file("exampleFiles", "exampleGeneSet.rda", package = "MLP")
pathExamplePValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
load(pathExampleGeneSet)
load(pathExamplePValues)
head(examplePValues)
head(exampleGeneSet)
mlpResult <- MLP(geneSet = exampleGeneSet, geneStatistic = examplePValues)
head(mlpResult)
}}

