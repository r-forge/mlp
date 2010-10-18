\name{plot.MLP}
\alias{plot.MLP}
\title{Plot the Results of an MLP Run...}
\usage{plot.MLP(x, y, type=c("barplot", "GOgraph", "quantileCurves"), ...)
}
\description{Plot the Results of an MLP Run}
\value{for type = "barplot", the midpoints of the barplot}
\arguments{\item{x}{object of class 'MLP'}
\item{y}{argument added to comply with generic; not used, defaults to NULL}
\item{type}{character of length one; one of 'barplot', 'GOgraph' or 'quantileCurves'}
\item{...}{further arguments for the plot functions for each type}
}
\examples{pathExampleMLPResult <- system.file("exampleFiles", "exampleMLPResult.rda", package = "MLP")
load(pathExampleMLPResult)
dev.new(width = 10, height = 10)
op <- par(mar = c(30, 10, 6, 2))
plot(exampleMLPResult, type = "barplot")
par(op)
plot(exampleMLPResult, type = "quantileCurves")
if (require(GO.db) && require(Rgraphviz)){
plot(exampleMLPResult, type = "GOgraph")
}}

