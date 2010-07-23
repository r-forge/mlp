#' Draw a Barplot for MLP Results 
#' @param object object of class MLP 
#' @param nRow number of rows of the MLP data frame to depict in the barplot; defaults to 20.
#' @param barColors vector of colors to use for the bars of the barplot; defaults to NULL; 
#'   if NULL, three gray shades are used reflecting the proportion of tested genes of a gene set
#'   versus the total number of genes in a geneset. If the proportion exceeds 75\%, the darkest
#'   shade is used; between 50 and 75\% a moderately dark shade is used; below 50\% a lighter gray
#'   shade is used.
#' @param main main title; if NULL (default) "Effect of the treatment on <geneSetSource> gene sets"
#' will be used
#' @return the midpoints of all the bars are returned invisibly (using the conventions of barplot); 
#'   an MLP-specific barplot is drawn to the current device;  
#' @seealso barplot
#' @export
mlpBarplot <- function (object, nRow = 20, barColors = NULL, main = NULL) {
  
  if (!inherits(object, "MLP")) 
    stop("'object' should be an object of class 'MLP' as produced by the MLP function")
  
  geneSetSource <- attr(object, "geneSetSource")
  mlpResults <- head(object, nRow)
  dat <- -log10(mlpResults$geneSetPValue)
  names(dat) <- rownames(mlpResults)
  
  if (is.null(barColors)){
    barColors <- rep("grey", length(dat))
    percentTested <- 100*(object$testedGeneSetSize/object$totalGeneSetSize)
    percentTested <- percentTested[1:length(dat)]
    
    barColors[percentTested >= 75] <- "grey50"
    barColors[percentTested < 75 & percentTested > 50] <- "grey60"
    barColors[percentTested <= 50] <- "grey70"
  } else {
    barColors <- barColors
  }
  
  if (is.null(object$geneSetDescription)) {
    if (geneSetSource %in% c("GOBP", "GOMF", "GOCC")) {
      allGOTerms <- as.list(Term(GOTERM))
      if (!all(names(dat) %in% names(allGOTerms))) 
        stop("Check the geneSetSource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
      descr <- allGOTerms[names(dat)]
    }
    else {
      if (geneSetSource == "KEGG") {
        allKEGGterms <- as.list(KEGGPATHID2NAME)
        geneSetNames <- gsub("^[[:alpha:]]{3}", "", names(dat))
        if (!all(geneSetNames %in% names(allKEGGterms))) 
          stop("Check the geneSetSource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
        descr <- allKEGGterms[geneSetNames]
      }
      else {
        if (!(geneSetSource %in% c("GOBP", "GOMF", "GOCC", 
                  "KEGG"))) {
          if (!all(rownames(object) %in% geneSetSource$PATHWAYID)) 
            stop("Check the geneSetSource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
          idx <- match(names(dat), geneSetSource$PATHWAYID)
          descr <- geneSetSource$PATHWAYNAME[idx]
        }
      }
    }
  }
  else {
    descr <- mlpResults$geneSetDescription
  }
  descriptionLength <- 60
  descr <- substr(descr, 1, descriptionLength)
  names(dat) <- paste(descr, " (", mlpResults$testedGeneSetSize, 
      "-", mlpResults$totalGeneSetSize, ")", sep = "")
  bottomMar <- 30
  op <- par(mar = c(bottomMar, 10, 6, 2))
  mp <- barplot(dat, xlab = "", main = "", border = "white", 
      las = 3, ylab = "", col = barColors)
  if (is.null(main)) {
    if (!(geneSetSource %in% c("GOBP", "GOMF", "GOCC", "KEGG"))){
      mainTitle <- "Effect of the treatment on the gene sets"
    } else {
      mainTitle <- paste("Effect of the treatment on", geneSetSource, "gene sets")
    }
  } else {
    mainTitle <- main
  }
  mtext(mainTitle, side = 2, line = 5)
  par(op)
  invisible(mp)
}
