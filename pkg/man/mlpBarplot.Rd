\name{mlpBarplot}
\alias{mlpBarplot}
\title{Draw a Barplot for MLP Results...}
\usage{mlpBarplot(object, nRow=20, main)}
\description{Draw a Barplot for MLP Results}
\value{the midpoints of all the bars are returned invisibly (using the conventions of barplot); 
an MLP-specific barplot is drawn to the current device;}
\seealso{barplot}
\arguments{\item{object}{object of class MLP}
\item{nRow}{number of rows of the MLP data frame to depict in the barplot; defaults to 20.}
\item{main}{main title; if NULL (default) "Effect of the treatment on <pathwaySource> gene sets"
will be used}}
