#' Graphical Representation of GO Based MLP Results
#' @param object object of class MLP (as produced by the MLP function) 
#' @param nRow number of GO IDs for which to produce the plot
#' @param mainTitle main title of the graph; defaults to 'GO graph' 
#' @return GO graph is plotted to the current device
#' @export
plotGOgraph <- function (object, nRow = 5, mainTitle = "GO graph") {
  if (!inherits(object, "MLP")) 
    stop("The 'object' argument should be an object of class 'MLP' as produced by the MLP function")
  
  require(GO.db)
  require(Rgraphviz)
  require(gplots)
  require(gmodels)
  require(gdata)
  require(gtools)
  require(GOstats)
  require(annotate)
  
  goids <- rownames(object)[1:nRow]
  ontology <- sub("GO", "", attributes(object)$pathwaySource)
  basicGraph <- GOGraph(goids, get(paste("GO", ontology, "PARENTS", 
              sep = "")))
  basicGraph <- removeNode("all", basicGraph)
  basicGraph <- removeNode(setdiff(nodes(basicGraph), rownames(object)), 
      basicGraph)
  basicGraph <- layoutGraph(basicGraph)
  pvalues <- object[nodes(basicGraph), "geneSetPValue"]
  names(pvalues) <- nodes(basicGraph)
  pvalues <- pvalues[!is.na(pvalues)]
  pvalues[pvalues == 0] <- min(pvalues[pvalues != 0])/10
  scores <- -log(pvalues)
  scores[scores <= 0.1] <- 0.1
  nColors <- round(max(scores) * 10)
  gocolors <- colorpanel(nColors, low = "lightyellow", high = "olivedrab")
  nodeFillColor <- rep("white", length(nodes(basicGraph)))
  names(nodeFillColor) <- nodes(basicGraph)
  nodeFillColor[names(scores)] <- gocolors[trunc(scores * 10)]
  nodeRenderInfo(basicGraph) <- list(fill = nodeFillColor)
  inbin <- object$totalGeneSetSize
  names(inbin) <- rownames(object)
  onchip <- object$testedGeneSetSize
  names(onchip) <- rownames(object)
  allcounts <- matrix(ncol = 2, nrow = length(inbin), dimnames = list(c(names(inbin)), 
          c("onchip", "inbin")))
  allcounts[names(onchip), "onchip"] <- onchip
  allcounts[names(inbin), "inbin"] <- inbin
  counts <- apply(allcounts, 1, function(x) {
        paste(x[1], x[2], sep = " - ")
      })
  terms <- getGOTerm(nodes(basicGraph))
  goTerm <- terms[[1]]
  goTerm <- sapply(goTerm, function(x) {
        paste(substr(x, start = 1, stop = 30), substr(x, start = 31, 
                stop = 60), sep = "\\\n")
      })
  nodeLabel <- paste(nodes(basicGraph), goTerm[nodes(basicGraph)], 
      counts[nodes(basicGraph)], sep = "\\\n")
  names(nodeLabel) <- nodes(basicGraph)
  nodeRenderInfo(basicGraph) <- list(label = nodeLabel)
  edgeRenderInfo(basicGraph)$arrowhead <- "none"
  nodeRenderInfo(basicGraph) <- list(shape = "ellipse")
  nodeRenderInfo(basicGraph) <- list(cex = 0.5)
  nodeRenderInfo(basicGraph) <- list(lWidth = 60)
  nodeRenderInfo(basicGraph) <- list(labelJust = "c")
  graphRenderInfo(basicGraph)$main <- mainTitle
  renderGraph(basicGraph)
  smartlegend("right", "bottom", legend = paste(c(" least", 
              " medium", " most"), " (scores ", round(max(scores)) * 
              c(2, 5, 8)/10, ")", sep = ""), fill = gocolors[round(max(scores)) * 
              c(2, 5, 8)], cex = 0.7)
}


