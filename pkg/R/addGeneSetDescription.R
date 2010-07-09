#' Utility function which adds the biological description of the gene sets as
#' a column to the return value of the MLP function (data frame) and orders
#' the MLP results by increasing p value for the gene set 
#' @param object object of class 'MLP' as produced by the 'MLP' function 
#' @param pathwaySource source to be used to construct the list of pathway categories; 
#' for public data sources, the user can specify a string (one of 'GOBP', 'GOMF', 'GOCC' or 'KEGG')
#' and BioC packages will be used to construct the list of pathway categories; 
#' for non-public data sources, the user can pass the pathway data as a dataframe with (at least) 
#' the following four columns: PATHWAYID, TAXID, PATHWAYNAME and GENEID. It is assumed all columns
#' are of type character. The 'pathwaySource' argument should be the same as the argument
#' provided to the getGeneSets function; defaults to NULL 
#' @param ... further arguments; currently none are used
#' @return the data frame as returned by MLP enriched with an additional column geneSetDescription, providing
#' a concise description of the gene set
#' @seealso \link{MLP}
#' @export
addGeneSetDescription <- 
    function(object, 
        pathwaySource = NULL, 
        ...)
{
  species <- attr(object, "species")
  if (is.null(pathwaySource))
    stop("Please provide the same source of gene sets as provided to the getGeneSets function. More info, see help.")
  if (any(!is.data.frame(pathwaySource) & !(pathwaySource %in% c("GOBP", "GOMF", "GOCC", "KEGG")))) 
    stop("Please provide the same source of gene sets as provided to the getGeneSets function. More info, see help.")
  if (any(is.data.frame(pathwaySource) & !(c("PATHWAYID", "TAXID", "PATHWAYNAME", "GENEID") %in% colnames(pathwaySource))))
    stop("Please provide the same source of gene sets as provided to the getGeneSets function. More info, see help.")
  if (any(pathwaySource %in% c("GOBP", "GOMF", "GOCC"))) {
    allGOTerms <- as.list(Term(GOTERM))
    geneSetNames <- rownames(object)
    if (!all(geneSetNames %in% names(allGOTerms))) stop("Check the pathwaySource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
    returnValue <- data.frame(object, 
        geneSetDescription = unlist(allGOTerms[geneSetNames]))
  }
  if (any(pathwaySource %in% c("KEGG"))) {
    allKEGGterms <- as.list(KEGGPATHID2NAME)
    geneSetNames <- gsub("^[[:alpha:]]{3}", "", rownames(object))
    if (!all(geneSetNames %in% names(allKEGGterms))) stop("Check the pathwaySource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
    returnValue <- data.frame(object, 
        geneSetDescription = unlist(allKEGGterms[geneSetNames]))
  }
  if (all(!(pathwaySource %in% c("GOBP", "GOMF", "GOCC", "KEGG")))) {
    if (!all(rownames(object) %in% pathwaySource$PATHWAYID)) stop("Check the pathwaySource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
    idx <- match(rownames(object), pathwaySource$PATHWAYID)
    returnValue <- data.frame(object, 
        geneSetDescription = pathwaySource$PATHWAYNAME[idx])
  }
  attr(returnValue, "species") <- species
  class(returnValue) <- c("MLP", class(returnValue))
  return(returnValue)
}
