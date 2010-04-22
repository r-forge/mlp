



###==============================================ANNOTATION FUNCTIONS===========================================


###Required Libraries:

###Organism library:
# library(affy)
# library(org.Hs.eg.db)   ### library(org.Mm.eg.db)  ### library(org.Rn.eg.db)
# library(GO.db)
# library(affy)
# library(AnnotationDbi)

#' TODO
#' @param Org 
#' @param Onto 
#' @param fNames 
#' @param nMin 
#' @param nMax 
#' @return 
#' @export
f.GOAnnotation <- function(Org="Mouse", Onto="BP", fNames=fNames, nMin=5,nMax=100){
  
  ###Authors: An De Bondt, Tine Casneuf, Nandini Raghavan.
  ###This code creates the necessary inputs for MLP analysis using libraries available from Bioconductor.
  
  ###Inputs: 
  ###Org=c("Mouse",Human","Rat")
  ###Onto=c("MF",BP","CC")
  ###fNames = feature names in the dataset (Entrez IDs).
  ###nMin = minimum number of genes in geneset.
  ###nMax = maximum number of genes in geneset.
  switch(Org,
      Mouse = {require(org.Mm.eg.db)
        go2entrez <- as.list(org.Mm.egGO2ALLEGS)},
      Human = {require(org.Mm.eg.db)
        go2entrez <- as.list(org.Mm.egGO2ALLEGS)},
      Rat = {require(org.Rn.eg.db)
        go2entrez <- as.list(org.Rn.egGO2ALLEGS)}
  )
  
  ## create first input object with GO info
  allGOontol <- eapply(GOTERM, Ontology)  
  allGOTerm  <- eapply(GOTERM, Term)
  
  ###Filter for GO terms related to specific Ontology:
  GOs <- names(which(allGOontol == Onto))
  go <- go2entrez[names(go2entrez) %in% GOs] 
  
  ###Filter for GO terms related to specific Ontology and specific Chip:
  go.eSet <- vector("list", length(go))
  for (j in 1:length(go)){go.eSet[[j]] <- unique(go[[j]][go[[j]] %in% fNames]  )}
  names(go.eSet) <- names(go) 
  ind0 <- ifelse(unlist(lapply(go.eSet,length)) >0,T,F)
  go.eSet <- go.eSet[ind0==T]
  ind1 <- ifelse(lapply(go.eSet,length) >= nMin & lapply(go.eSet,length) <= nMax, TRUE, FALSE)
  go.eSet <- go.eSet[ind1==T]
  
  return(go.eSet)
  
}

#' TODO 
#' @param x 
#' @return 
#' @export
f.GOInputMLP <- function(x = go.eSet){
  
  ###Authors: An De Bondt, Tine Casneuf, Nandini Raghavan.
  ###This code converts the output of f.GOAnnotation to the necessary format for input for MLP analysis.
  
  ###Input: go.eSet output from f.GOAnnotation.
  out <- lapply(names(go.eSet), function(goid) {
        i.genes <- unique(go.eSet[[goid]])
        matrix(c(rep(x = goid, times = length(i.genes)), i.genes), ncol = 2, byrow = FALSE)
      })
  out <- do.call(rbind, out)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out <- out[!is.na(out[,2]), ]
  colnames(out) <- c('GO', 'Gene.ID')
  out[, 1] <- sub('GO:', '', out[, 1])
  out[, 1] <- as.numeric(out[, 1])
  out[, 2] <- as.numeric(out[, 2])
  out <- unique(out)
  out <- as.matrix(out) ### listed as a matrix as needed for ML
  return(out)
  
}
