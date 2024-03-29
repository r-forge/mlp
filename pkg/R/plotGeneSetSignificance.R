#' Plot the Significance for the Genes of a Given Gene Set 
#' @param geneSet object of class 'geneSetMLP' as produced by function getGeneSets
#' @param geneSetIdentifier identifier of the gene set for which a significance plot should be produced;
#'   character of length one  
#' @param geneStatistic vector of gene statistics (e.g. p values) 
#' @param annotationPackage name of the annotation package to be used (without .db extension);
#'   character of length one
#' @param barColors vector of colors to use for the bars of the barplot; defaults to NULL
#'  in which case 'grey50' is used 
#' @return no return value 
#' @examples
#' pathExamplePValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
#' pathExampleGeneSet <- system.file("exampleFiles", "exampleGeneSet.rda", package = "MLP")
#' pathExampleMLPResult <- system.file("exampleFiles", "exampleMLPResult.rda", package = "MLP")
#' load(pathExampleGeneSet)
#' load(pathExamplePValues)
#' load(pathExampleMLPResult) 
#' annotationPackage <- if (require(mouse4302mmentrezg)) "mouse4302mmentrezg" else "mouse4302"
#' geneSetID <- rownames(exampleMLPResult)[1]
#' dev.new(width = 10, height = 10)
#' op <- par(mar = c(25, 10, 6, 2))
#' plotGeneSetSignificance(
#'     geneSet = exampleGeneSet, 
#'     geneSetIdentifier = geneSetID, 
#'     geneStatistic = examplePValues, 
#'     annotationPackage = annotationPackage
#' )
#' par(op)
#' @export
plotGeneSetSignificance <- function(geneSet, geneSetIdentifier, geneStatistic, annotationPackage, barColors = NULL){
  
  if (!inherits(geneSet, "geneSetMLP"))
    stop("geneSet should be an object of class 'geneSetMLP' as produced by 'getGeneSets'")
  
  if (!geneSetIdentifier %in% names(geneSet))
    stop("Please provide as 'geneSetName' a gene set name belonging to 'geneSets', i.e. the group of gene sets specified.")
  
  require(annotate)
  require(paste(annotationPackage, ".db", sep = ""), character.only = TRUE)
  
  
  entrezids <- geneSet[[geneSetIdentifier]]
  entrezids <- entrezids[entrezids %in% names(geneStatistic)]
  genePvalues <- geneStatistic[entrezids]
  genePvalues <- sort(genePvalues)
  psids <- if (length(grep("_at", names(genePvalues))) == 0){
    paste(names(genePvalues), "_at", sep = "")
  } else {
    names(genePvalues)
  }
  names(genePvalues) <- paste(
      unlist(lookUp(psids, annotationPackage, "SYMBOL")),
      unlist(lookUp(psids, annotationPackage, "GENENAME")),
      sep = ":")
  names(genePvalues) <- substr(names(genePvalues), 1, 60)
  
  barColors <- if (is.null(barColors)) "grey50" else barColors
  
  barplot(-log10(genePvalues), xlab = "", 
      main = paste("Significance of tested genes involved in gene set", geneSetIdentifier), 
      border = "white", col = barColors,
      las = 3, ylab = "Significance")
}
