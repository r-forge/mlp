\name{getGeneSets}
\alias{getGeneSets}
\title{Prepare Pathway Data for the MLP Function...}
\usage{getGeneSets(species="Mouse", geneSetSource, entrezIdentifiers)
}
\description{Prepare Pathway Data for the MLP Function}
\details{The return value of the getGeneSets function has as primary use
to serve as geneSet argument for the MLP function}
\value{object of class geneSetMLP which is essentially a named list of pathway categories. 
Each list component contains a vector of Entrez IDs related to that particular pathway}
\arguments{\item{species}{character vector of length one indicating the species, one of
'Mouse', 'Human' or 'Rat'; defaults to 'Mouse'.}
\item{geneSetSource}{source to be used to construct the list of pathway categories; 
for public data sources, the user can specify a string (one of 'GOBP', 'GOMF', 'GOCC' or 'KEGG')
and BioC packages will be used to construct the list of pathway categories; 
for non-public data sources, the user can pass the pathway data as a dataframe with (at least) 
the following four columns: PATHWAYID, TAXID, PATHWAYNAME and GENEID. It is assumed all columns
are of type character.}
\item{entrezIdentifiers}{Entrez ID identifiers used to subset the relevant gene set}
}
\examples{if (require(GO.db) && require(org.Mm.eg.db)){
pathExampleData <- system.file("exampleFiles", "expressionSetGcrma.rda", package = "MLP")
pathExamplePValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
load(pathExampleData)
load(pathExamplePValues)
geneSet <- getGeneSets(species = "Mouse", geneSetSource = "GOBP", entrezIdentifiers = names(examplePValues)[1:2000])
head(geneSet)
}}

