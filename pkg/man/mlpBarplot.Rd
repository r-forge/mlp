\name{mlpBarplot}
\alias{mlpBarplot}
\title{Draw a Barplot for MLP Results...}
\usage{mlpBarplot(object, pathwaySource, nRow=20, main)}
\description{Draw a Barplot for MLP Results}
\value{the midpoints of all the bars are returned invisibly (using the conventions of barplot); 
an MLP-specific barplot is drawn to the current device;}
\seealso{barplot}
\arguments{\item{object}{object of class MLP}
\item{pathwaySource}{source to be used to construct the list of pathway categories; 
for public data sources, the user can specify a string (one of 'GOBP', 'GOMF', 'GOCC' or 'KEGG')
and BioC packages will be used to construct the list of pathway categories; 
for non-public data sources, the user can pass the pathway data as a dataframe with (at least) 
the following four columns: PATHWAYID, TAXID, PATHWAYNAME and GENEID. It is assumed all columns
are of type character.}
\item{main}{main title; if NULL (default) "Effect of the treatment on <pathwaySource> gene sets"
will be used}}
