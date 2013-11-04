#
# Copyright John Reid 2008
#


library(topGO)
library("Rgraphviz")


buildGoData <- function(
	geneNames,
	genes2Go,
	ontology = "BP"
) {
	#
	# Builds the data structure for GO enrichment testing
	#
    selGenes <- sample(geneNames, 1)
	inGenes <- factor(as.integer(geneNames %in% selGenes))
    names(inGenes) <- geneNames
	new(
		"topGOdata",
		ontology = ontology,
		allGenes = inGenes,
	    annot = annFUN.gene2GO,  ## the new annotation function
      	gene2GO = genes2Go     ## the gene ID to GOs dataset
	)
}


updateGoData <- function(
	goData,
	selGenes
) {
	#
	# Update a topGOdata object with a new list of interesting genes
	#
	geneNames <- genes(goData)
	inGenes <- factor(as.integer(geneNames %in% selGenes))
	names(inGenes) <- geneNames
	updateGenes(goData, inGenes)
}


goEnrichment <- function(
	goData, 
	pValueThreshold,
	prefix,			## prefix for filenames
	pdfSW = FALSE
) {
	#
	# Analyse
	#
	test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
	resultFis <- getSigGroups(goData, test.stat)

	test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test", cutOff = 0.01)
	resultElim <- getSigGroups(goData, test.stat)

	test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", cutOff = 0.01, sigRatio = "ratio")
	resultWeight <- getSigGroups(goData, test.stat)



	#
	# Look at results
	#
	allRes <- GenTable(goData, classic=resultFis, elim=resultElim, weight=resultWeight, 
		ranksOf='classic', orderBy="weight")
	allRes <- subset(allRes, weight <= pValueThreshold)
	showSigOfNodes(goData, score(resultFis), firstTerms = 5, useInfo = "all")
	printGraph(goData, resultWeight, firstSigNodes = 5, fn.prefix = prefix,
		pdfSW = pdfSW)

	print(allRes)
	allRes
}

subsetResultsByClassic <- function(results, threshold) {
	subset(results, classic <= threshold)
}

subsetResultsByWeight <- function(results, threshold) {
	subset(results, weight <= threshold)
}

getPValues <- function(...) {
	resList <- list(...)
	resList <- lapply(resList, score)
}



# goEnrichment(goData, selGenes, 0.5, 'test')

test <- function() {
	library(ALL)
	data(ALL)
	affyLib <- annotation(ALL)
	library(package = affyLib, character.only = TRUE)
	library(genefilter)
	f1 <- pOverA(0.25, log2(100))
	f2 <- function(x) (IQR(x) > 0.5)
	ff <- filterfun(f1, f2)
	eset <- ALL[genefilter(ALL, ff), ]
	geneNames <- featureNames(eset)
	geneNames <- geneNames[1:30]
	length(geneNames)
	myInterestedGenes <- sample(geneNames, 10)
	geneList <- factor(as.integer(geneNames %in% myInterestedGenes))
	names(geneList) <- geneNames
	str(geneList)
	goData <- new("topGOdata", ontology = "MF", allGenes = geneList,
		annot = annFUN.hgu, affyLib = affyLib)

}

