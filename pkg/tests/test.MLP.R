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

system.time(goGeneSet <- getGeneSets(species = "Mouse", pathwaySource = "GOBP", eset = expressionSetGcrma))
goGeneSet[1:3]
# output changes with annotation version !

y <- normDat.p[,1]
names(y) <- featureNames(expressionSetGcrma)

y[1:10]
# 100009600    100012    100017    100019 100034251 100036521 100037258 100037278 
# 0.4328583 0.7448996 0.6088859 0.1845008 0.2312761 0.7865153 0.7772888 0.1037431 
# 100038570 100038635 
# 0.1368744 0.3272610 

mlpObject <- MLP(geneSet = goGeneSet, geneStatistic = y, minGenes = 5, maxGenes = 100, rowPermutations = TRUE, 
    nPermutations = 6, smoothPValues = TRUE)



mlpObject[1:10, ]
# output changes with annotation version !

plotGOgraph(object = mlpObject, ontology = "BP", annotation = "mouse4302mmentrezg")

pdf(file = "test10.pdf", width = 10, height = 10)
plot(mlpObject, nRow = 10) # by default:  type = "barplot"
dev.off()

unlink("test10.pdf")

if (FALSE){
  pdf(file = "test5.pdf", width =10, height = 10)
  mlpBarplot(object = mlpObject, pathwaySource = "GOBP", nRow = 10, descriptionLength = 5)
  dev.off()
  
  unlink("test5.pdf")
  
  pdf(file = "test100.pdf", width =10, height = 20)
  mlpBarplot(object = mlpObject, pathwaySource = "GOBP", nRow = 10, descriptionLength = 100)
  dev.off()
  
  unlink("test100.pdf")
}

plot(mlpObject, type = "quantileCurves")
