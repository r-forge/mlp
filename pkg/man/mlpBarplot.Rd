\name{mlpBarplot}
\alias{mlpBarplot}
\title{Draw a Barplot for MLP Results...}
\usage{mlpBarplot(object, nRow=20, barColors, main)}
\description{Draw a Barplot for MLP Results}
\value{the midpoints of all the bars are returned invisibly (using the conventions of barplot); 
an MLP-specific barplot is drawn to the current device;}
\seealso{barplot}
\arguments{\item{object}{object of class MLP}
\item{nRow}{number of rows of the MLP data frame to depict in the barplot; defaults to 20.}
\item{barColors}{vector of colors to use for the bars of the barplot; defaults to NULL; 
if NULL, three gray shades are used reflecting the proportion of tested genes of a gene set
versus the total number of genes in a geneset. If the proportion exceeds 75\%, the darkest
shade is used; between 50 and 75\% a moderately dark shade is used; below 50\% a lighter gray
shade is used.}
\item{main}{main title; if NULL (default) "Effect of the treatment on <geneSetSource> gene sets"
will be used}}
