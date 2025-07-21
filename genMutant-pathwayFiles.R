library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(tidyverse)


#'Filters significant genes using p-value and Spearman correlation cutoff,
#'returns a vector of gene symbols. 
getSigGenes = function(geneSet, 
                       pvalueCutoff = .05, 
                       spearmanCutOff = .75){
  filteredGenes=geneSet %>%
    filter(padjust < pvalueCutoff & 
             abs(Spearman_Correlation) > spearmanCutOff) %>%
    select(Gene)%>%
    distinct()%>%
    pull(Gene)
    
}
#'runs enrichment on msigdb using kegg, reactome and pid genesets
runMsigDb = function(geneSet, 
                     pvalueCutoff = .1, 
                     minGSSize = 25, 
                     maxGSSize = 300, 
                     backgroundGenes,
                     countCutOff = 4, 
                     spearmanRes,
                     tissue){
  subCollections = list('PID',
                        'REACTOME',
                        'KEGG_LEGACY',
                        'KEGG_MEDICUS'
  )
  for(s in subCollections){
    pathways = msigdbr( species = 'mouse', collection = "C2", subcollection = paste0('CP:',
                                                                                     s)) %>%
      select(gs_name, gene_symbol)
    
    pathwaysGeneSets = enricher(
      gene = geneSet,
      TERM2GENE = pathways,
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      universe = backgroundGenes
    ) %>%
      as.data.frame()
    getMutants(spearmanRes = spearmanRes, 
               goOutput = pathwaysGeneSets, 
               tissue = tissue, 
               enrichment = s)
  }
  
}

#'Runs Gene Ontology enrichment analysis on a set of genes, uses the 
#'pvalue, minGSS and maxGSS cutoffs defined in Sambahwas paper
runGeneOntology = function(geneSet, 
                           pvalueCutoff = .1, 
                           minGSSize = 25, 
                           maxGSSize = 300, 
                           backgroundGenes,
                           countCutOff = 4){
  
  ego <- enrichGO(gene          = geneSet,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "SYMBOL",
                  ont           = "BP",
                  universe      = backgroundGenes,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = pvalueCutoff,
                  minGSSize     = minGSSize,
                  maxGSSize     = maxGSSize)%>%
    clusterProfiler::simplify(cutoff = .6)%>%
    as.data.frame()%>%
    filter(Count > countCutOff)
  
  return(ego)
}

#'Function to get mutants associated with each pathway and write the gene-mutants associations to .tsv
#'calculates the absolute spearman coef 
getMutants = function(spearmanRes, 
                      goOutput, 
                      pvalueCutoff = .05, 
                      spearmanCutOff = .75, 
                      tissue, 
                      enrichment = ''){
  directory = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/',tissue,'ClusterIntegration/spearmanRes/pathways/', enrichment)
  
  if(!dir.exists(directory)) 
    #Create the directory if it does not exist
  dir.create(directory, 
             recursive = TRUE)

  for (i in 1:10){
   genes = goOutput$geneID[i] %>%
     strsplit(split = "/") %>%
     unlist()
   pathwayAssocations=spearmanRes%>%
     filter(Gene %in% genes, 
            padjust < pvalueCutoff, 
            abs(Spearman_Correlation) >spearmanCutOff)%>%
     mutate(absoluteSpearman = abs(Spearman_Correlation))
    
  pathway = goOutput$Description[i]
  path = paste0(directory, '/', pathway, '.tsv')
  write_tsv(pathwayAssocations, 
            path)
  }
}


#Wrapper script for running the entire process of generating clusters
genClusters = function(backgroundSet, 
                       inputSet,
                       tissue){
  
  filteredInput = getSigGenes(inputSet)
  goOutput=runGeneOntology(filteredInput, 
                 backgroundGenes = backgroundSet$SYMBOL)
  getMutants(spearmanRes = inputSet, 
             goOutput = goOutput, 
             tissue = tissue, 
             enrichment = 'geneOntology')
  
  runMsigDb(geneSet = filteredInput,
            pvalueCutoff = .1, 
            minGSSize = 25, 
            maxGSSize = 300, 
            backgroundGenes = backgroundSet$SYMBOL,
            spearmanRes = inputSet,
            tissue = tissue)
  
}

spearmanRes = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/djClusterIntegration/spearmanRes/sigDjBarseqGenesStrongMutantSpearmanCorrelationEstimatedPvalues.tsv')
backgroundGenes = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqJejunumOutputs/lrtResults/allTestedGenes.tsv')

genClusters(backgroundSet = backgroundGenes, 
             inputSet = spearmanRes, 
             tissue = 'jejunum')

spearmanRes = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/colonClusterIntegration/spearmanRes/sigRes.tsv')
backgroundGenes = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/allGenesTab.tsv')

genClusters(backgroundSet = backgroundGenes, 
            inputSet = spearmanRes, 
            tissue = 'colon')

