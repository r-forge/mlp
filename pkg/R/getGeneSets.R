#' Prepare Pathway Data for the MLP Function
#' 
#' The return value of the getGeneSets function has as primary use
#' to serve as geneSet argument for the MLP function
#' @param species character vector of length one indicating the species, one of
#' 'Mouse', 'Human' or 'Rat'; defaults to 'Mouse'. 
#' @param pathwaySource source to be used to construct the list of pathway categories; 
#' for public data sources, the user can specify a string (one of 'GOBP', 'GOMF', 'GOCC' or 'KEGG')
#' and BioC packages will be used to construct the list of pathway categories; 
#' for non-public data sources, the user can pass the pathway data as a dataframe with (at least) 
#' the following four columns: PATHWAYID, TAXID, PATHWAYNAME and GENEID. It is assumed all columns
#' are of type character.
#' @param eset ExpressionSet object which will be used to subset the relevant Entrez ID identifiers
#' in the relevant pathway category list components; the subset taken corresponds to the feature
#' names of this ExpressionSet object
#' @return object of class geneSetMLP which is essentially a named list of pathway categories. 
#' Each list component contains a vector of Entrez IDs related to that particular pathway
#' @import AnnotationDbi
#' @export
getGeneSets <- 
    function (species = "Mouse", pathwaySource = NULL, eset) 
{
  if (!species %in% c("Mouse", "Human", "Rat")) 
    stop("The 'species' argument should be one of 'Mouse', 'Human' or 'Rat)")
  if (is.null(pathwaySource)) 
    stop("Please provide a source of gene sets. For more info, see the help page.")
  if (!is.data.frame(pathwaySource) & !(pathwaySource %in% 
        c("GOBP", "GOMF", "GOCC", "KEGG"))) 
    stop("The 'pathwaySource' argument should be one of 'GOBP', 'GOMF', 'GOCC', 'KEGG' or a data.frame.  More info, see help.")
# if (is.data.frame(pathwaySource) & !(c("PATHWAYID", "TAXID", 
#           "PATHWAYNAME", "GENEID") %in% colnames(pathwaySource))) 
  # because otherwise there is a warning which is not needed
  if (is.data.frame(pathwaySource)) if (!(c("PATHWAYID", "TAXID", 
              "PATHWAYNAME", "GENEID") %in% colnames(pathwaySource))) 
      stop("The pathwaySource as data.frame should have at least the 4 columns 'PATHWAYID', 'TAXID', 'PATHWAYNAME' and 'GENEID'. More info on their content, see help.")
  if (pathwaySource %in% c("GOBP", "GOMF", "GOCC")) {
    require(GO.db)
    ontology <- sub("GO", "", pathwaySource)
    switch(species, Mouse = {
          require(org.Mm.eg.db)
          geneSetToEntrez <- as.list(org.Mm.egGO2ALLEGS)
        }, Human = {
          require(org.Hs.eg.db)
          geneSetToEntrez <- as.list(org.Hs.egGO2ALLEGS)
        }, Rat = {
          require(org.Rn.eg.db)
          geneSetToEntrez <- as.list(org.Rn.egGO2ALLEGS)
        })
    switch(ontology, BP = {
          GOs <- names(as.list(GOBPANCESTOR))
        }, MF = {
          GOs <- names(as.list(GOMFANCESTOR))
        }, CC = {
          GOs <- names(as.list(GOCCANCESTOR))
        })
    go <- geneSetToEntrez[names(geneSetToEntrez) %in% GOs]
    geneSets <- lapply(go, unique)
    anyGenesInGeneSet <- ifelse(unlist(lapply(geneSets, length)) > 
            0, TRUE, FALSE)
    geneSets <- geneSets[anyGenesInGeneSet]
  }
  else {
    if (pathwaySource == "KEGG") {
      require(KEGG.db)
      geneSetToEntrez <- as.list(KEGGPATHID2EXTID)
      switch(species, Mouse = {
            prefix <- "mmu"
          }, Human = {
            prefix <- "hsa"
          }, Rat = {
            prefix <- "rno"
          })
      geneSets <- geneSetToEntrez[grep(prefix, names(geneSetToEntrez))]
    }
    else {
      switch(species, Mouse = {
            TAXID <- "10090"
          }, Human = {
            TAXID <- "9606"
          }, Rat = {
            TAXID <- "10116"
          })
      if (!(TAXID %in% pathwaySource$TAXID)) 
        stop(paste("Please check the available TAXIDs in your pathwaySource. The pathwaySource object does not contain gene sets for ", 
                species, sep = ""))
      pathwaySource <- pathwaySource[pathwaySource$TAXID == 
              TAXID, ]
      geneSets <- by(pathwaySource$GENEID, INDICES = pathwaySource$PATHWAYID, 
          FUN = list)
      geneSets <- lapply(geneSets, as.character)
    }
  }
  entrezIds <- sub("_at", "", featureNames(eset))
  # correct way (= not deleting genes for which there are no expression data) 
  # which currently still gives an error when running MLP
  tfidx <- sapply(geneSets, function(geneSet){sum(geneSet %in% entrezIds) > 0})
  geneSets <- geneSets[tfidx]
  # wrong way
  #geneSets <- sapply(geneSets, function(x) x[x %in% entrezIds])
  attr(geneSets, "species") <- species
  class(geneSets) <- c("geneSetMLP", class(geneSets))
  return(geneSets)
}
