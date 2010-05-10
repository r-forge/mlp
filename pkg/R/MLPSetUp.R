# TODO integrate goAnnotation inside MLP and change MLP to have genes and
# GO terms as character rather than numeric

#' This function extracts the GO terms from one of the following BioC packages
#' and maps them to the subset of Entrez IDs that are present in the list of feature
#' names provided (org.Mm.eg.db, org.Hs.eg.db, org.Rn.eg.db) It also thresholds the GO terms to those
#' that have between minGenes and maxGenes number of Entrez IDs.
#' @param organism character vector of length one indicating the organism, one of
#' Mouse, Human or Rat; defaults to Mouse. 
#' @param ontology one of three ontologies from the GO database 
#' @param featureNames character vector of the feature names (Entrez IDs) to be included
#' @param minGenes minimum number of genes in a gene set; defaults to 5
#' @param maxGenes maximum number of genes in a gene set; defaults to 100
#' @return data structure to be used as geneSet argument for goInputMLP; list
#' of GO terms each of which has a list of Entrez IDs as a vector
#' @import AnnotationDbi
#' @export
goAnnotation <- function(organism = "Mouse", ontology = "BP", 
		featureNames, minGenes = 5, maxGenes = 100){

  if (!organism %in% c("Mouse", "Human", "Rat"))
	  stop("The 'organism' argument should be one of 'Mouse', 'Human' or 'Rat'.")
  if (!ontology %in% c("MF", "BP", "CC"))
	  stop("The 'ontology' argument should be one of 'MF', 'BP' or 'CC'.")
  
#  switch(organism,
#      Mouse = {
#        require(org.Mm.eg.db)
#        goToEntrez <- as.list(org.Mm.egGO2ALLEGS) # TODO
#      },
#      Human = {
#        require(org.Hs.eg.db)
#        goToEntrez <- as.list(org.Hs.egGO2ALLEGS)
#      },
#      Rat = {
#        require(org.Rn.eg.db)
#        goToEntrez <- as.list(org.Rn.egGO2ALLEGS)
#      }
#  )
#  
  if (organism == "Mouse"){
    require(org.Mm.eg.db)
    goToEntrez <- as.list(org.Mm.egGO2ALLEGS) # TODO
  } else if (organism == "Human"){
    require(org.Hs.eg.db)
    goToEntrez <- as.list(org.Hs.egGO2ALLEGS)
  } else { # Rat
    require(org.Rn.eg.db)
    goToEntrez <- as.list(org.Rn.egGO2ALLEGS)
  }
  # create first input object with GO info
  allGOontol <- eapply(GOTERM, Ontology)  
  allGOTerm  <- eapply(GOTERM, Term)
  
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

#' Prepare the input geneSet for the MLP function 
#' @param goInFeatureNames output from goAnnotation
#' @return data structure ready to be used as geneSet argument for the MLP
#'   function
#' @export
goInputMLP <- function(goInFeatureNames){

	out <- lapply(names(goInFeatureNames), function(goid) {
				i.genes <- unique(goInFeatureNames[[goid]])
				matrix(c(rep(x = goid, times = length(i.genes)), i.genes), ncol = 2, byrow = FALSE)
			})
	out <- do.call(rbind, out)
	out <- as.data.frame(out, stringsAsFactors = FALSE)
	out <- out[!is.na(out[,2]), ]
	colnames(out) <- c('GO', 'Gene.ID')
	out[,1] <- sub('GO:', '', out[,1])
	out[,1] <- as.numeric(out[,1])
	out[,2] <- as.numeric(out[,2])
	out <- unique(out)
	out <- as.matrix(out) ### listed as a matrix as needed for ML
	return(out)	
}

#' TODO fix S3 export
#' TODO 
#' TODO (low priority currently) add print method for MLP objects


#' Summary function for MLP objects; this function combines
#'   the MLP results with the biological description of the
#'   gene sets
#' @param object object of class 'MLP' as produced by the 'MLP' function 
#' @param goInFeatureNames list of gene sets; each gene set (component) contains
#'   the member genes as a character vector. Typically, the output of goAnnotation
#'   is used.  
#' @param ... further arguments; currently none are used
#' @return TODO
#' @method summary MLP
#' @S3method summary MLP
#' @export
summary.MLP <- function(object, goInFeatureNames, ...){
	
	allGOTerm  <- eapply(GOTERM, Term)
	geneSetNames <- names(goInFeatureNames)
	nGenesInGeneset <- unlist(lapply(goInFeatureNames, length))
	
	returnValue <- data.frame(genesetSize = nGenesInGeneset, object[, c("genesetStatistic", "genesetPValue")],
			genesetDescription = unlist(allGOTerm[geneSetNames]))
	
	returnValue <- returnValue[order(returnValue[,"genesetPValue"]),]
	row.names(returnValue) <- geneSetNames
	
	return(returnValue)
}
