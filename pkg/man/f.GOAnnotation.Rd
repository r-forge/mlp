\name{f.GOAnnotation}
\alias{f.GOAnnotation}
\title{This function extracts the GO terms from one of the following BioC packages...}
\usage{f.GOAnnotation(organism="Mouse", ontology="BP", featureNames,
    minGenes=5, maxGenes=100)}
\description{This function extracts the GO terms from one of the following BioC packages
and maps them to the subset of Entrez IDs that are present in the list of feature
names provided (org.Mm.eg.db, org.Hs.eg.db, org.Rn.eg.db) It also thresholds the GO terms to those
that have between minGenes and maxGenes number of Entrez IDs.}
\value{data structure to be used as geneSet argument for f.GOInputMLP; list
of GO terms each of which has a list of Entrez IDs as a vector}
\arguments{\item{organism}{character vector of length one indicating the organism, one of
Mouse, Human or Rat; defaults to Mouse.}
\item{ontology}{one of three ontologies from the GO database}
\item{featureNames}{character vector of the feature names (Entrez IDs) to be included}
\item{minGenes}{minimum number of genes in a gene set}
\item{maxGenes}{maximum number of genes in a gene set}}
