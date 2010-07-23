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
  if (length(grep("_at", names(pvalues))) == 0){
    psids <- paste(names(genePvalues), "_at", sep = "")
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
