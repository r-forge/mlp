\name{addGeneSetDescription}
\alias{addGeneSetDescription}
\title{Utility function which adds the biological description of the gene sets as...}
\usage{addGeneSetDescription(object, geneSetSource)}
\description{Utility function which adds the biological description of the gene sets as
a column to the return value of the MLP function (data frame)}
\value{the data frame as returned by MLP enriched with an additional column geneSetDescription, providing
a concise description of the gene set}
\seealso{\link{MLP}}
\arguments{\item{object}{object of class 'MLP' as produced by the 'MLP' function}
\item{geneSetSource}{source to be used to construct the list of pathway categories; 
for public data sources, the user can specify a string (one of 'GOBP', 'GOMF', 'GOCC' or 'KEGG')
and BioC packages will be used to construct the list of pathway categories; 
for non-public data sources, the user can pass the pathway data as a dataframe with (at least) 
the following four columns: PATHWAYID, TAXID, PATHWAYNAME and GENEID. It is assumed all columns
are of type character. The 'geneSetSource' argument should be the same as the argument
provided to the getGeneSets function; defaults to NULL}}
