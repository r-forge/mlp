
require(MLP)	
set.seed(479)

# This is just the expressionset for this experiment.

pathExampleData <- system.file("exampleFiles", "expressionSetGcrma.rda", package = "MLP")
load(pathExampleData)

#Libraries needed
library(limma)
library(org.Mm.eg.db) # for mouse

exprs(expressionSetGcrma)[1:2,]
#              2760     2763     2765     2766     2768     2769     2761     2762    2764     2767
#100009600 2.371111 2.170060 2.233383 2.180717 2.325886 2.239441 2.297301 2.409001 2.49458 2.115814
#100012    2.176163 2.318876 2.419263 2.223307 2.585125 2.346060 2.292061 2.336415 2.47979 2.361981
#              2770     2771
#100009600 2.371262 2.267459
#100012    2.330418 2.520918

pData(expressionSetGcrma)
#     sample subGroup sampleColor
#2760      1        1     #FF0000
#2763      4        1     #FF0000
#2765      6        1     #FF0000
#2766      7        1     #FF0000
#2768      9        1     #FF0000
#2769     10        1     #FF0000
#2761      2        2     #0000FF
#2762      3        2     #0000FF
#2764      5        2     #0000FF
#2767      8        2     #0000FF
#2770     11        2     #0000FF
#2771     12        2     #0000FF

pData(expressionSetGcrma)$subGroup1 <- ifelse(pData(expressionSetGcrma)$subGroup==1,"WT","KO")

###==============================================GENERATING LIMMA p-VALUES=================================

# boxplot(data.frame(exprs(expressionSetGcrma))
normDat  <- normalizeQuantiles(exprs(expressionSetGcrma), ties=TRUE)
subGroup <- pData(expressionSetGcrma)$subGroup
design <- model.matrix(~ -1 +factor(subGroup ))

colnames(design) <- c("group1", "group2")
contrast.matrix <- makeContrasts(group1-group2, levels=design)
fit <- lmFit(normDat,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
normDat.p <- fit2$p.value

normDat.p[1:5]
#[1] 0.4328583 0.7448996 0.6088859 0.1845008 0.2312761

##=================================================INPUTS FOR MLP==========================================

go.eSet <- goAnnotation(organism="Mouse", ontology="BP", featureNames = featureNames(expressionSetGcrma),
    minGenes = 5, maxGenes = 100)
go.eSet[1:3]
#$`GO:0000002`
#[1] "18975"  "27393"  "27395"  "27397"  "83945"  "382985"
#
#$`GO:0000018`
# [1] "12053"  "12144"  "12487"  "13990"  "15978"  "16183"  "16189"  "17350"  "17685"  "17688" 
#[11] "19264"  "20371"  "20852"  "21803"  "21939"  "50931"  "57765"  "65113"  "104271"
#
#$`GO:0000038`
#[1] "15488"  "19305"  "54326"  "94180"  "171281" "217698"

# x <- goInputMLP(go.eSet) ###First Input for MLP (both columns are numeric).
# x[1:10,]
#   GO Gene.ID
#1   2   18975
#2   2   27393
#3   2   27395
#4   2   27397
#5   2   83945
#6   2  382985
#7  18   12053
#8  18   12144
#9  18   12487
#10 18   13990

y <- normDat.p[,1]
names(y) <- featureNames(expressionSetGcrma)

y[1:10]
# 100009600    100012    100017    100019 100034251 100036521 100037258 100037278 
# 0.4328583 0.7448996 0.6088859 0.1845008 0.2312761 0.7865153 0.7772888 0.1037431 
# 100038570 100038635 
# 0.1368744 0.3272610 

out.MLP <- MLP(geneSet = go.eSet, geneStatistic = y, minGenes = 5, maxGenes = 100, rowPermutations = TRUE, 
    nPermutations = 6, smoothPValues = TRUE)

tmp <- summary(object = out.MLP, goInFeatureNames = go.eSet)

tmp[1:10,]
#            genesetSize genesetStatistic genesetPValue
# GO:0000002          47        0.8484050  1.812210e-08
# GO:0000018          17        1.0826887  2.794698e-06
# GO:0000038          10        1.2423499  1.654631e-05
# GO:0000041          17        1.0291455  1.685600e-05
# GO:0000045          15        1.0367976  4.477027e-05
# GO:0000050          15        1.0367976  4.477027e-05
# GO:0000052          17        0.9905982  5.507421e-05
# GO:0000060          43        0.7476719  6.437538e-05
# GO:0000070          25        0.8604682  1.135816e-04
# GO:0000075          25        0.8570210  1.288780e-04
#                                                                                          genesetName
# GO:0000002                                                                          female pregnancy
# GO:0000018                                                             monocarboxylic acid transport
# GO:0000038                                                                         tRNA modification
# GO:0000041                                                                       embryo implantation
# GO:0000045                   antigen processing and presentation of peptide antigen via MHC class II
# GO:0000050         antigen processing and presentation of exogenous peptide antigen via MHC class II
# GO:0000052 antigen processing and presentation of peptide or polysaccharide antigen via MHC class II
# GO:0000060                                                       antigen processing and presentation
# GO:0000070                                             establishment or maintenance of cell polarity
# GO:0000075                                                                      vacuole organization

