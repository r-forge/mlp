#' Graphical Representation of GO Based MLP Results
#' @param object object of class MLP (as produced by the MLP function) 
#' @param ontology ontology 
#' @param annotation annotation
#' @param nRow nRow
#' @param mainTitle main title of the graph; defaults to 'GO graph' 
#' @return GO graph is plotted to the current device
#' @export
plotGOgraph <- function(object, ontology, annotation, nRow = 5, mainTitle = "GO graph"){
  
  if (!inherits(object, "MLP"))
    stop("The 'object' argument should be an object of class 'MLP' as produced by the MLP function")
  
  require(GO.db)
  require(paste(annotation, ".db", sep = ""), character.only = TRUE)
  require(Rgraphviz)
  require(gplots)
  require(gmodels)
  require(gdata)
  require(gtools)
  require(GOstats)
  require(annotate)
  
  goToEntrezArray <- as.list(get(paste(annotation, "GO2ALLPROBES", sep = "")))
  
  species <- attr(object, "species")
  switch(species, 
      Mouse = {
        require(org.Mm.eg.db)
        goToEntrez <- as.list(org.Mm.egGO2ALLEGS)
      }, Human = {
        require(org.Hs.eg.db)
        goToEntrez <- as.list(org.Hs.egGO2ALLEGS)
      }, Rat = {
        require(org.Rn.eg.db)
        goToEntrez <- as.list(org.Rn.egGO2ALLEGS)
      })
  allGOTerms <- as.list(Term(GOTERM))
  goids <- rownames(object)[1:nRow]
  basicGraph <- GOGraph(goids, get(paste("GO", ontology, "PARENTS", sep = ""))) ## Biological process
  basicGraph <- removeNode("all", basicGraph)
  basicGraph <- removeNode(setdiff(nodes(basicGraph), rownames(object)), basicGraph)
  basicGraph <- layoutGraph(basicGraph)
  
  # Node colors
  pvalues <- object[nodes(basicGraph), "geneSetPValue"]
  names(pvalues) <- nodes(basicGraph)
  pvalues <- pvalues[!is.na(pvalues)]
  pvalues[pvalues == 0] <- min(pvalues[pvalues != 0])/10
  scores <- -log(pvalues)
  # min score set at 0.1 because the coloring will be based on a color vector using round(score*10) as the index
  scores[scores <= 0.1] <- 0.1
  nColors <- round(max(scores)*10)
  gocolors <- colorpanel(nColors,low = "lightyellow", high = "olivedrab")
  nodeFillColor <-  rep("white", length(nodes(basicGraph)))
  names(nodeFillColor) <- nodes(basicGraph)
  nodeFillColor[names(scores)] <- gocolors[trunc(scores*10)]
  nodeRenderInfo(basicGraph) <- list(fill = nodeFillColor)
  
  # Node labels
  allGOontol <- as.list(Ontology(GOTERM))
  GOs <- names(which(allGOontol == ontology))
  
  # Gene Counts based on inhouse data sources
  # #inbin per GOID
  go <- goToEntrez[names(goToEntrez) %in% GOs]
  inbin <- sapply(go, function(x){length(unique(x))})
  # #on chip per GOID
  go <- goToEntrezArray[names(goToEntrezArray) %in% GOs]
  onchip <- sapply(go, function(x){length(unique(x))})
  
  # merge counts in one table
  allcounts <- matrix(ncol = 2, nrow = length(inbin), dimnames = list(c(names(inbin)),c("onchip", "inbin")))
  allcounts[names(onchip),"onchip"] <- onchip
  allcounts[names(inbin),"inbin"] <- inbin
  
  # make a text vector to add as annotation in the GO graph
  counts <- apply(allcounts, 1, function(x) {
        paste(x[1], x[2], sep = " - ")
      })
  
  # GOIDs and GOterms and merge with the counts
  terms <- getGOTerm(nodes(basicGraph))
  goTerm <- terms[[1]]
  goTerm <- sapply(goTerm, function(x) {
        paste(substr(x, start = 1, stop = 30), substr(x, start = 31, stop = 60), sep = "\\\n")
      })
  
  nodeLabel <- paste(nodes(basicGraph), goTerm[nodes(basicGraph)], counts[nodes(basicGraph)], sep = "\\\n")
  names(nodeLabel) <- nodes(basicGraph)
  nodeRenderInfo(basicGraph) <- list(label = nodeLabel)
  
  # other settings
  edgeRenderInfo(basicGraph)$arrowhead <- "none"
  nodeRenderInfo(basicGraph) <- list(shape = "ellipse")
  nodeRenderInfo(basicGraph) <- list(cex = 0.5)
  nodeRenderInfo(basicGraph) <- list(lWidth = 60)
  nodeRenderInfo(basicGraph) <- list(labelJust = "c")
  graphRenderInfo(basicGraph)$main <- mainTitle
  renderGraph(basicGraph)
  smartlegend("right", "bottom", 
      legend = paste(c(" least", " medium", " most"), " (scores ", round(max(scores))*c(2,5,8)/10, ")", sep = ""), 
      fill = gocolors[round(max(scores))*c(2, 5, 8)], cex = 0.7)
}


