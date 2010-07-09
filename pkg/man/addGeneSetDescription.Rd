\name{addGeneSetDescription}
\alias{addGeneSetDescription}
\title{Utility function which adds the biological description of the gene sets as...}
\usage{addGeneSetDescription(object, pathwaySource, ...)}
\description{Utility function which adds the biological description of the gene sets as
a column to the return value of the MLP function (data frame) and orders
the MLP results by increasing p value for the gene set}
\value{TODO}
\arguments{\item{object}{object of class 'MLP' as produced by the 'MLP' function}
\item{pathwaySource}{source to be used to construct the list of pathway categories; 
for public data sources, the user can specify a string (one of 'GOBP', 'GOMF', 'GOCC' or 'KEGG')
and BioC packages will be used to construct the list of pathway categories; 
for non-public data sources, the user can pass the pathway data as a dataframe with (at least) 
the following four columns: PATHWAYID, TAXID, PATHWAYNAME and GENEID. It is assumed all columns
are of type character. The 'pathwaySource' argument should be the same as the argument
provided to the getGeneSets function; defaults to NULL}
\item{...}{further arguments; currently none are used}}
