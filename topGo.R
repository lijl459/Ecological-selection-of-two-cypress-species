#install.packages("BiocManager")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

BiocManager::install("topGO", force = TRUE)

library("topGO")

setwd("\\PSGs\\go")


gene2go <- readMappings("LJL17009_target.emapper.go.txt")

DEGs_path <- "Ccas.csv"
prefix <- "Ccash.PSG"


DEGs <- read.table(DEGs_path)


genenames <- names(gene2go)

genelist <- factor(as.integer(genenames %in% DEGs[, 1])) 

names(genelist) <- genenames
GOdata <- new("topGOdata", ontology = "BP", allGenes = genelist, 
              annot = annFUN.gene2GO, gene2GO = gene2go)




GOdataMF <- new("topGOdata", ontology = "MF", allGenes=genelist, annot = annFUN.gene2GO, gene2GO = gene2go, nodeSize=10)
GOdataBP <- new("topGOdata", ontology = "BP", allGenes=genelist, annot = annFUN.gene2GO, gene2GO = gene2go, nodeSize=10)
GOdataCC <- new("topGOdata", ontology = "CC", allGenes=genelist, annot = annFUN.gene2GO, gene2GO = gene2go, nodeSize=10)
resultMF <- runTest(GOdataMF, algorithm = "weight01", statistic = "fisher")
resultBP <- runTest(GOdataBP, algorithm = "weight01", statistic = "fisher")
resultCC <- runTest(GOdataCC, algorithm = "weight01", statistic = "fisher")
#resultMF <- runTest(GOdataMF, algorithm = "classic", statistic = "fisher")
#resultBP <- runTest(GOdataBP, algorithm = "classic", statistic = "fisher")
#resultCC <- runTest(GOdataCC, algorithm = "classic", statistic = "fisher")
allGOMF = usedGO(object = GOdataMF)
allGOBP = usedGO(object = GOdataBP)
allGOCC = usedGO(object = GOdataCC)
tableMF <- GenTable(GOdataMF, classicFisher = resultMF, ranksOf = "classicFisher",topNodes = length(allGOMF), numChar=1000)
tableBP <- GenTable(GOdataBP, classicFisher = resultBP, ranksOf = "classicFisher",topNodes = length(allGOBP), numChar=1000)
tableCC <- GenTable(GOdataCC, classicFisher = resultCC, ranksOf = "classicFisher",topNodes = length(allGOCC), numChar=1000)

write.table(tableMF,file=paste0(prefix, "_MF.csv"),sep=",", row.name = FALSE, col.names=TRUE)
write.table(tableBP,file=paste0(prefix, "_BP.csv"),sep=",", row.name = FALSE, col.names=TRUE)
write.table(tableCC,file=paste0(prefix, "_CC.csv"),sep=",", row.name = FALSE, col.names=TRUE)



printGraph(GOdataMF, resultMF, firstSigNodes = 10, fn.prefix = "union_BayeScenv_RDA.MF", useInfo = "all", pdfSW = TRUE)
printGraph(GOdataBP, resultBP, firstSigNodes = 10, fn.prefix = "union_BayeScenv_RDA.BP", useInfo = "all", pdfSW = TRUE)
printGraph(GOdataCC, resultCC, firstSigNodes = 10, fn.prefix = "union_BayeScenv_RDA.CC", useInfo = "all", pdfSW = TRUE)




