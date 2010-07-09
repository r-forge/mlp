
require(MLP)	
set.seed(479)

# This is just the expressionset for this experiment.

pathExampleData <- system.file("exampleFiles", "expressionSetGcrma.rda", package = "MLP")
load(pathExampleData)

# Libraries needed
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

system.time(goGeneSet <- getGeneSets(species = "Mouse", pathwaySource = "GOBP", eset = expressionSetGcrma))
goGeneSet[1:3]
# $`GO:0000002`
# [1] "18975"  "19819"  "27393"  "27395"  "27397"  "57813"  "83945"  "226153"
# [9] "382985"
# 
# $`GO:0000018`
#  [1] "12053"  "12144"  "12487"  "13990"  "15978"  "16183"  "16189"  "17350" 
#  [9] "17685"  "17688"  "19264"  "20371"  "20852"  "21803"  "21939"  "22210" 
# [17] "50931"  "57765"  "58186"  "65113"  "104271"
# 
# $`GO:0000038`
# [1] "15488"  "19305"  "26458"  "26569"  "54326"  "94180"  "171281" "217698"


y <- normDat.p[,1]
names(y) <- featureNames(expressionSetGcrma)

y[1:10]
# 100009600    100012    100017    100019 100034251 100036521 100037258 100037278 
# 0.4328583 0.7448996 0.6088859 0.1845008 0.2312761 0.7865153 0.7772888 0.1037431 
# 100038570 100038635 
# 0.1368744 0.3272610 

outMLP <- MLP(geneSet = goGeneSet, geneStatistic = y, minGenes = 5, maxGenes = 100, rowPermutations = TRUE, 
    nPermutations = 6, smoothPValues = TRUE)

tmp <- addGeneSetDescription(object = outMLP, pathwaySource = "GOBP",
    eset = expressionSetGcrma)

tmp[1:10,]
#            geneSetSize geneSetStatistic geneSetPValue
# GO:0007565          47        0.8484050  8.438006e-05
# GO:0006400          10        1.2423499  2.534188e-04
# GO:0015718          17        1.0826887  3.335550e-04
# GO:0007566          17        1.0291455  7.618846e-04
# GO:0002495          15        1.0367976  1.023345e-03
# GO:0019886          15        1.0367976  1.023345e-03
# GO:0002504          17        0.9905982  1.340495e-03
# GO:0035088          12        1.0522212  1.511658e-03
# GO:0019882          43        0.7476719  2.152273e-03
# GO:0002478          22        0.8934001  2.468006e-03
#                                                                                   geneSetDescription
# GO:0007565                                                                          female pregnancy
# GO:0006400                                                                         tRNA modification
# GO:0015718                                                             monocarboxylic acid transport
# GO:0007566                                                                       embryo implantation
# GO:0002495                   antigen processing and presentation of peptide antigen via MHC class II
# GO:0019886         antigen processing and presentation of exogenous peptide antigen via MHC class II
# GO:0002504 antigen processing and presentation of peptide or polysaccharide antigen via MHC class II
# GO:0035088                                establishment or maintenance of apical/basal cell polarity
# GO:0019882                                                       antigen processing and presentation
# GO:0002478                          antigen processing and presentation of exogenous peptide antigen


plotGOgraph(object = outMLP, ontology = "BP", annotation = "mouse4302mmentrezg")

