#' This function calculates p-values for each gene set based on row permutations
#' of the gene p values or column permutations of the expression matrix; the p values
#' can be obtained either as individual gene set p values or p values based on smoothing
#' across gene sets of similar size.
#' @param geneSet is the input list of gene sets (components) and gene IDs (character vectors). 
#' A gene set can, for example, be a GO category with for each category Entrez gene identifiers; 
#' The \link{getGeneSets} function can be used to construct the geneSet argument for different
#' pathway sources.
#' @param geneStatistic is either a named numeric vector (if rowPermutations is TRUE)
#'   or a numeric matrix of pvalues (if rowPermutations is FALSE). The names of the numeric vector
#'   or row names of the matrix should represent the gene IDs.
#' @param minGenes minimum number of genes in a gene set for it to be considered (lower threshold
#'    for gene set size)
#' @param maxGenes maximum number of genes in a gene set for it to be considered (upper threshold
#'    for gene set size)
#' @param rowPermutations logical indicating whether to use row permutations (TRUE; default) or column 
#'    permutations (FALSE) 
#' @param nPermutations is the number of simulations. By default 100 permutations are conducted.
#' @param smoothPValues logical indicating whether one wants to calculate smoothed cut-off thresholds (TRUE; default)
#'    or not (FALSE).
#' @param criticalValues vector of quantiles at which p values for each gene set are desired
#' @param df degrees of freedom for the smooth.spline function used in getSmoothedPValues 
#' @param addGeneSetDescription logical indicating whether a column with the gene set description be added to
#' the output data frame; defaults to TRUE.
#' @return data frame with four (or five) columns: totalGeneSetSize, testedGeneSetSize, geneSetStatistic and geneSetPValue
#' and (if addDescription is set to TRUE) geneSetDescription; the rows of the data frame are ordered by 
#' ascending geneSetPValue.
#'
#' @references Raghavan, Nandini et al. (2007). The high-level similarity of some disparate gene expression measures,
#' Bioinformatics, 23, 22, 3032-3038.
#' @export
MLP <- function (geneSet, geneStatistic, minGenes = 5, maxGenes = 100, 
    rowPermutations = TRUE, nPermutations = 100, smoothPValues = TRUE, 
    criticalValues = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999), 
    df = 9, addGeneSetDescription = TRUE) 
{
  if (!inherits(geneSet, "geneSetMLP")) {
    stop("The 'geneSet' should be an object of class 'geneSetMLP' as produced by getGeneSets")
  }
  species <- attr(geneSet, "species")
  if (is.vector(geneStatistic)) 
    if (!rowPermutations) 
      stop("If rowPermutations is set to FALSE, the 'geneStatistic' needs to be a matrix")
  geneStatistic <- as.matrix(geneStatistic)
  if (!is.numeric(geneStatistic) & !is.matrix(geneStatistic)) {
    warning("Argument 'geneStatistic' should be a numeric vector or matrix")
  }
  if (!rowPermutations) {
    nPermutations <- ncol(geneStatistic) - 2
  }
  
  totalGeneSetSize <- sapply(geneSet, length)
  geneSets <- sapply(geneSet, function(x) x[x %in% rownames(geneStatistic)])
  
  geneSetSize <- sapply(geneSets, function(x) length(x))
  geneSetIndices <- which(geneSetSize >= minGenes & geneSetSize <= 
          maxGenes)
  geneSets <- geneSets[geneSetIndices]
  geneStatistic <- na.omit(geneStatistic)
  geneStatistic <- geneStatistic[grep("^[[:digit:]]+$", rownames(geneStatistic)), 
      , drop = FALSE]
  mapResult <- mapGeneSetStatistic(geneSets, geneStatistic)
  pValueNAlist <- lapply(mapResult, function(x) which(!is.na(x)))
  filterFunction <- function(x, y) {
    x[y]
  }
  reducedMapResult <- mapply(filterFunction, mapResult)
  reducedGeneSet <- mapply(filterFunction, geneSets)
  n2 <- nrow(geneStatistic)
  w0 <- t(mlpStatistic(reducedMapResult))
  n11 <- length(geneSets)
  if (!rowPermutations) {
    w <- mlpStatistic(reducedMapResult)
  }
  else {
    p1 <- apply(matrix(1:n2, n2, nPermutations), 2, sample)
    y1 <- matrix(geneStatistic[p1, 1], n2, nPermutations)
    row.names(y1) <- row.names(geneStatistic)
    rm(p1)
    mapResultPermuted <- mapGeneSetStatistic(reducedGeneSet, 
        y1)
    w <- t(mlpStatistic(mapResultPermuted))
  }
  if (smoothPValues) {
    if (!rowPermutations) {
      warning("Cannot smooth p-values if ind.sim = FALSE")
    }
    pw0 <- getSmoothedPValues(w0, w[, -1], q.cutoff = criticalValues, 
        df = df)
  }
  else {
    pw0 <- getIndividualPValues(w0, w)
  }
  
  res <- data.frame(testedGeneSetSize = w0[, 1], geneSetStatistic = w0[, 
          2], geneSetPValue = pw0)
  res <- res[order(res$geneSetPValue), ]
  res <- data.frame(totalGeneSetSize = totalGeneSetSize[rownames(res)], res)
  class(res) <- c("MLP", class(res))
  if (addGeneSetDescription){
    pathwaySource <- attr(geneSet, "pathwaySource")
    res <- addGeneSetDescription(object = res, pathwaySource = pathwaySource)
  }
  if (is.null(attr(res, "species"))) 
    attr(res, "species") <- species 
  attr(res, "pathwaySource") <- attr(geneSet, "pathwaySource")
  if (!is.null(attr(pw0, "quantileCurveInformation")))
    attr(res, "quantileCurveInformation") <- attr(pw0, "quantileCurveInformation")
  return(res)
}
