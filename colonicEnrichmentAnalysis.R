library(clusterProfiler)
library(msigdbr)
library(tidyverse)

getMsigDb = function(inputGenes, 
                     backgroundGenes, 
                     category, 
                     subCategory, 
                     padjCutoff = 0.1,
                     minGSSCutoff = 25,
                     maxGSSCutoff = 500, 
                     ...) {
  
  geneSet = msigdbr( species = 'mouse', category = category, subcollection = subCategory) %>%
    select(gs_name, gene_symbol)
  
  pathwaysGeneSets = enricher(
    gene = inputGenes,
    TERM2GENE = geneSet,
    pvalueCutoff = padjCutoff,
    minGSSize = minGSSCutoff,
    maxGSSize = maxGSSCutoff,
    universe = backgroundGenes,
    ...
  ) %>%
    as.data.frame()
  
  return(pathwaysGeneSets)
}

plotTop10Terms = function(enricherOutput, 
                          title, 
                          countCutOff = 4, 
                          ...) {
  p = enricherOutput %>%
    filter(Count > countCutOff) %>%
    arrange(p.adjust, desc(Count)) %>%
    head(n = 10) %>%
    ggplot(aes(y = reorder(Description, -p.adjust), 
               x = -log10(p.adjust), 
               fill = Count)) +
    geom_col() +
    scale_fill_gradient(
      low = "darkblue", 
      high = "red"
    ) +
    labs(title = title, 
         y = 'Pathway', 
         x = '-log10(p.adjust)',
         fill = 'Gene Count')
  
  plot(p)
}

filterSignificantGenes = function(deseqDf, 
                                  foldChangeCutoff = .75,
                                  padjCutoff = 0.1) {
  
  sigGenes=deseqDf %>%
    filter(padj < padjCutoff) %>%
    filter(abs(log2FoldChange) > foldChangeCutoff) 
  return(sigGenes)
}

runEnrichment = function(deseqResults, 
                         universeGenes,
                         compairsonName){
  
  
  subCollections = list('PID',
                        'REACTOME',
                        'KEGG_LEGACY',
                        'KEGG_MEDICUS'
                  )
  sigGenes = filterSignificantGenes(deseqResults)
  
  
   'Just set up this way bc any later time point was set to the denominator in DeSeq2. '
  decreasing = sigGenes%>%
    filter(log2FoldChange > 0)%>%
    pull(SYMBOL)
  
  increasing = sigGenes%>%
    filter(log2FoldChange < 0)%>%
    pull(SYMBOL)
  
  for(s in subCollections){
    pathwaysGeneSets = getMsigDb(
      inputGenes = increasing,
      backgroundGenes = universeGenes,
      category = 'C2',
      subCategory = s,
      padjCutoff = 0.1,
      minGSSCutoff = 25,
      maxGSSCutoff = 500
    )
    
    if(nrow(pathwaysGeneSets) > 0){
      p = plotTop10Terms(
        enricherOutput = pathwaysGeneSets, 
        title = paste0(compairsonName, ' - ', s)
      )
      
      ggsave(
        filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', compairsonName,'/temporalIncrease/',s, '_', s, '.png'),
        plot = p,
        width = 8, height = 6
      )
      write_tsv(
        pathwaysGeneSets,
        file = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', compairsonName,'/temporalIncrease/',s, '_', s, '.tsv')
      )
    }
  }
  
  for(s in subCollections){
    pathwaysGeneSets = getMsigDb(
      inputGenes = decreasing,
      backgroundGenes = universeGenes,
      category = 'C2',
      subCategory = s,
      padjCutoff = 0.1,
      minGSSCutoff = 25,
      maxGSSCutoff = 500
    )
    
    if(nrow(pathwaysGeneSets) > 0){
      p = plotTop10Terms(
        enricherOutput = pathwaysGeneSets, 
        title = paste0(compairsonName, ' - ', s)
      )
      ggsave(
        filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', compairsonName,'/temporalDecrease/',s, '_', s, '.png'),
        plot = p,
        width = 8, height = 6
      )
      write_tsv(
        pathwaysGeneSets,
        file = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', compairsonName,'/temporalDecrease/',s, '_', s, '.tsv')
      )
    }
  }
  
}


# Define background genes
allGenes = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/allGenesTab.tsv')
allGenes = allGenes$SYMBOL

day1VDay3 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day1Vday3.tsv')

runEnrichment(
  deseqResults = day1VDay3,
  universeGenes = allGenes,
  compairsonName = 'day1Vday3'
)

day3VDay7 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day3Vday7.tsv')

runEnrichment(
  deseqResults = day3VDay7,
  universeGenes = allGenes,
  compairsonName = 'day3Vday7'
)


day7VDay14 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day7Vday14.tsv')

runEnrichment(
  deseqResults = day7VDay14,
  universeGenes = allGenes,
  compairsonName = 'day7Vday14'
)

#' cluster enrichment is not returning functions for a lot of the pathways when using all genes 
#' tested as the background set. It might make more sense to try and the significant genes 
#' instead. I think logically, the question becomes instead of the what genes are overall enriched 
#' at each time point, of those significant genes, what pathways are enriched? 
#' JK that only really improved cluster 1 i think, the rest did not return any results at all. 

backgroundSet=rbind(day1VDay3, 
      day3VDay7, 
      day7VDay14)%>%
  filterSignificantGenes()%>%
  .$SYMBOL%>%
  unique()

#'Setting min cluster size for consideration. 
clusterSize = 5
clustersDf = read_tsv('DeSeq2/barseqColonicOutputs/manualGeneClusteringClusters.tsv')

names(clustersDf)[1] <- "SYMBOL"

  
clusters=clustersDf%>%
  distinct()%>%
  group_by(clusterId)%>%
  summarise(nGenes = n())%>%
  filter(nGenes > clusterSize)%>%
  pull(clusterId)
  
for (c in clusters){
  clusterGenes = clustersDf %>%
    filter(clusterId == c) %>%
    distinct()%>%
    pull(SYMBOL)
  
  subCollections = list('PID',
                        'BIOCARTA',
                        'REACTOME',
                        'KEGG_LEGACY',
                        'KEGG_MEDICUS'
  )
  
  for(s in subCollections){
    pathwaysGeneSets = getMsigDb(
      inputGenes = clusterGenes,
      backgroundGenes = allGenes,
      category = 'C2',
      subCategory = s,
      padjCutoff = 0.1,
      minGSSCutoff = 25,
      maxGSSCutoff =300
    )
    
    if(nrow(pathwaysGeneSets) > 0){
      p = plotTop10Terms(
        enricherOutput = pathwaysGeneSets, 
        title = paste0('Cluster ', c, ' - ', s),
        countCutOff = 4
      )
      if(dir.exists(paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/cluster', c, '/')) == FALSE){
        dir.create(paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/cluster', c, '/'))
      }
      ggsave(
        filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/cluster', c, '/', s, '_', s, '.png'),
        plot = p,
        width = 8, height = 6
      )
      
      write_tsv(
        pathwaysGeneSets,
        file = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/cluster', c, '/', s, '_', s, '.tsv')
      )
    }
  }
}
  
  

