
require(MLP)	

# This is just the expressionset for this experiment.

pathExampleData <- system.file("exampleFiles", "expressionSetGcrma.rda", package = "MLP")
load(pathExampleData)

#Libraries needed
library(limma)
library(affy)
library(GO.db)
library(AnnotationDbi)
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

colnames(design) <- c("group1","group2")
contrast.matrix <- makeContrasts(group1-group2,levels=design)
fit <- lmFit(normDat,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
normDat.p <- fit2$p.value

normDat.p[1:5]
#[1] 0.4328583 0.7448996 0.6088859 0.1845008 0.2312761

##=================================================INPUTS FOR MLP==========================================

go.eSet <- f.GOAnnotation(Org="Mouse", Onto="BP", fNames=featureNames(expressionSetGcrma),nMin=5,nMax=100)
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

x <- f.GOInputMLP(go.eSet) ###First Input for MLP (both columns are numeric).
x[1:10,]
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

y <- data.frame("Entrez.ID"= as.numeric(featureNames(expressionSetGcrma)),"p.Value"=normDat.p)
#Warning message:
#In data.frame(Entrez.ID = as.numeric(featureNames(expressionSetGcrma)),  :
#  NAs introduced by coercion

dim(y)
#[1] 16395     2

y <- as.matrix(na.omit(y))  ###Deletes all rows with non-numeric featureNames.
dim(y)
#[1] 16331     2

y[1:10,]
#          Entrez.ID group1...group2
#100009600 100009600       0.4328583
#100012       100012       0.7448996
#100017       100017       0.6088859
#100019       100019       0.1845008
#100034251 100034251       0.2312761
#100036521 100036521       0.7865153
#100037258 100037258       0.7772888
#100037278 100037278       0.1037431
#100038570 100038570       0.1368744
#100038635 100038635       0.3272610

###==============================================RUN MLP===========================================

out.MLP <- f.GOFiSH(x = x, y = y, nmin = 5, nmax = 100, ind.sim = 1, nsim = 200, ind.smooth = 1)
tmp <- f.GOOutputMLP(x1=go.eSet,x2 = out.MLP)

tmp[1:10,]
#        Geneset
#943  GO:0015718
#660  GO:0007565
#661  GO:0007566
#202  GO:0002504
#201  GO:0002495
#1074 GO:0019886
#590  GO:0007163
#1499 GO:0035088
#200  GO:0002478
#1072 GO:0019882
#                                                                                  Geneset.Name
#943                                                              monocarboxylic acid transport
#660                                                                           female pregnancy
#661                                                                        embryo implantation
#202  antigen processing and presentation of peptide or polysaccharide antigen via MHC class II
#201                    antigen processing and presentation of peptide antigen via MHC class II
#1074         antigen processing and presentation of exogenous peptide antigen via MHC class II
#590                                              establishment or maintenance of cell polarity
#1499                                establishment or maintenance of apical/basal cell polarity
#200                           antigen processing and presentation of exogenous peptide antigen
#1072                                                       antigen processing and presentation
#     Geneset.Size MLP.Observed      P.Value
#943             6    2.4976051 4.357394e-07
#660            38    0.8221092 1.724650e-03
#661            17    1.0222076 1.785178e-03
#202            16    1.0209136 2.217215e-03
#201            15    1.0367976 2.272937e-03
#1074           15    1.0367976 2.272937e-03
#590            23    0.8920877 3.231891e-03
#1499           12    1.0522212 3.629602e-03
#200            22    0.8934001 3.684589e-03
#1072           43    0.7543953 4.324502e-03



###==============================================================================================