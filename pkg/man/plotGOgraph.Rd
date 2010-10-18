\name{plotGOgraph}
\alias{plotGOgraph}
\title{Graphical Representation of GO Based MLP Results...}
\usage{plotGOgraph(object, nRow=5, main)
}
\description{Graphical Representation of GO Based MLP Results}
\value{GO graph is plotted to the current device}
\arguments{\item{object}{object of class MLP (as produced by the MLP function)}
\item{nRow}{number of GO IDs for which to produce the plot}
\item{main}{main title of the graph; if NULL (default) the main title is set to 'GO graph'}
}
\examples{if (require(GO.db) && require(Rgraphviz)){
pathExampleMLPResult <- system.file("exampleFiles", "exampleMLPResult.rda", package = "MLP")
load(pathExampleMLPResult)
plotGOgraph(exampleMLPResult, main = "GO Graph")
}}

