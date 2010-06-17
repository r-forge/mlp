#' This function extracts the GO terms from one of the following BioC packages
#' and maps them to the subset of Entrez IDs that are present in the list of feature
#' names provided (org.Mm.eg.db, org.Hs.eg.db, org.Rn.eg.db) It also thresholds the GO terms to those
#' that have between minGenes and maxGenes number of Entrez IDs.
#' @param organism character vector of length one indicating the organism, one of
#' Mouse, Human, Rat or Dog; defaults to Mouse. 
#' @param ontology one of three ontologies from the GO database 
#' @param featureNames character vector of the feature names (Entrez IDs) to be included
#' @param minGenes minimum number of genes in a gene set; defaults to 5
#' @param maxGenes maximum number of genes in a gene set; defaults to 100
#' @return list of GO terms. Each list component contains a vector of Entrez IDs related
#'   to the pa
#' @import AnnotationDbi
#' @export
goAnnotation <- function(organism = "Mouse", ontology = "BP", 
		featureNames, minGenes = 5, maxGenes = 100){

  if (!organism %in% c("Mouse", "Human", "Rat", "Dog"))
	  stop("The 'organism' argument should be one of 'Mouse', 'Human', 'Rat', 'Dog'.")
  if (!ontology %in% c("MF", "BP", "CC"))
	  stop("The 'ontology' argument should be one of 'MF', 'BP' or 'CC'.")
  
  switch(organism,
      Mouse = {
        require(org.Mm.eg.db)
        goToEntrez <- as.list(org.Mm.egGO2ALLEGS)
      },
      Human = {
        require(org.Hs.eg.db)
        goToEntrez <- as.list(org.Hs.egGO2ALLEGS)
      },
      Rat = {
        require(org.Rn.eg.db)
        goToEntrez <- as.list(org.Rn.egGO2ALLEGS)
      },
      Dog = {
        require(org.Cf.eg.db)
        goToEntrez <- as.list(org.Cf.egGO2ALLEGS)
      }
  )
  
  # create first input object with GO info
  allGOontol <- as.list(Ontology(GOTERM))  
  
  # filter for GO terms related to specific ontology
  GOs <- names(which(allGOontol == ontology))
  go <- goToEntrez[names(goToEntrez) %in% GOs] 
  
  # Filter for GO terms related to specific Ontology and 
  # specified list of feature names (potentially a subset of a chip)
  goInFeatureNames <- vector("list", length(go))
  for (j in 1:length(go)){
	  goInFeatureNames[[j]] <- unique(go[[j]][go[[j]] %in% featureNames])
  }
  names(goInFeatureNames) <- names(go) 
  anyGenesInGeneSet <- ifelse(unlist(lapply(goInFeatureNames, length)) > 0, TRUE, FALSE)
  goInFeatureNames <- goInFeatureNames[anyGenesInGeneSet]
  passMinAndMaxThresholds <- ifelse(lapply(goInFeatureNames,length) >= minGenes & lapply(goInFeatureNames,length) <= maxGenes, TRUE, FALSE)
  goInFeatureNames <- goInFeatureNames[passMinAndMaxThresholds]
  
  return(goInFeatureNames)
}

#' TODO (low priority currently) add print method for MLP objects


#' Summary function for MLP objects; this function combines
#'   the MLP results with the biological description of the
#'   gene sets
#' @param object object of class 'MLP' as produced by the 'MLP' function 
#' @param geneSets list of gene sets; each gene set (component) contains
#'   the member genes as a character vector. Typically, the output of goAnnotation
#'   is used.  
#' @param ... further arguments; currently none are used
#' @return TODO
#' @export
addGeneSetDescription <- function(object, geneSets, ...){
	allGOTerms  <- as.list(Term(GOTERM))
	geneSetNames <- names(geneSets)
	nGenesInGeneSet <- unlist(lapply(geneSets, length))
	
	returnValue <- data.frame(geneSetSize = nGenesInGeneSet, object[, c("geneSetStatistic", "geneSetPValue")],
			geneSetDescription = unlist(allGOTerms[geneSetNames]))
	orderByPValue <- order(returnValue[,"geneSetPValue"])
	returnValue <- returnValue[orderByPValue,]
	row.names(returnValue) <- geneSetNames[orderByPValue]
	
	return(returnValue)
}
