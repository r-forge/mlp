#' Draw a Barplot for MLP Results 
#' @param object object of class MLP 
#' @param nRow number of rows of the MLP data frame to depict in the barplot; defaults to 20.
#' @param main main title; if NULL (default) "Effect of the treatment on <pathwaySource> gene sets"
#' will be used
#' @return the midpoints of all the bars are returned invisibly (using the conventions of barplot); 
#'   an MLP-specific barplot is drawn to the current device;  
#' @seealso barplot
#' @export
mlpBarplot <- function (object, nRow = 20, main = NULL) 
{
  if (!inherits(object, "MLP"))
    stop("'object' should be an object of class 'MLP' as produced by the MLP function")
  
  pathwaySource <- attr(object, "pathwaySource")
  
  mlpResults <- head(object, nRow)
  dat <- -log(mlpResults$geneSetPValue)
  names(dat) <- rownames(mlpResults)
  barColors <- rep("grey", length(dat))
  barColors[1:5] <- c("grey10", "grey20", "grey30", "grey40", 
      "grey50")
  
  if (is.null(object$geneSetDescription)){  
    if (pathwaySource %in% c("GOBP", "GOMF", "GOCC")) {
      allGOTerms <- as.list(Term(GOTERM))
      if (!all(names(dat) %in% names(allGOTerms))) 
        stop("Check the pathwaySource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
      descr <- allGOTerms[names(dat)]
    }
    else {
      if (pathwaySource == "KEGG") {
        allKEGGterms <- as.list(KEGGPATHID2NAME)
        geneSetNames <- gsub("^[[:alpha:]]{3}", "", names(dat))
        if (!all(geneSetNames %in% names(allKEGGterms))) 
          stop("Check the pathwaySource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
        descr <- allKEGGterms[geneSetNames]
      }
      else {
        if (!(pathwaySource %in% c("GOBP", "GOMF", "GOCC", 
                  "KEGG"))) {
          if (!all(rownames(object) %in% pathwaySource$PATHWAYID)) 
            stop("Check the pathwaySource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
          idx <- match(names(dat), pathwaySource$PATHWAYID)
          descr <- pathwaySource$PATHWAYNAME[idx]
        }
      }
    }
  } else {
    descr <- mlpResults$geneSetDescription
  }    
  
  descriptionLength <- 60 # TODO reimplement in grid
  descr <- substr(descr, 1, descriptionLength)
  names(dat) <- paste(descr, " (", 
      mlpResults$testedGeneSetSize, "-",
      mlpResults$totalGeneSetSize, 
      ")", sep = "")
  
  # bottomMar <- 6 + min(descriptionLength/2, nchar(descr))
  bottomMar <- 30
  op <- par(mar = c(bottomMar, 10, 6, 2))
  mp <- barplot(dat, xlab = "", main = "", border = "white", 
      las = 3, ylab = "", col = barColors)
  if (is.null(main)) {
    mainTitle <- paste("Effect of the treatment on", pathwaySource, 
        "gene sets")
  }
  else {
    mainTitle <- main
  }
  mtext(mainTitle, side = 2, line = 5)
  par(op)
  invisible(mp)
}
