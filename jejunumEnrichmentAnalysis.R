library(clusterProfiler)
library(msigdbr)
library(tidyverse)
library(stringr)

getMsigDb = function(inputGenes, 
                     backgroundGenes, 
                     category, 
                     subCategory, 
                     padjCutoff = 0.1,
                     minGSSCutoff = 10,
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
    mutate(wrappedDescription = 
             str_wrap(str_replace_all(Description, "_", " "), width = 30))%>%
    ggplot(aes(
      y = fct_reorder(wrappedDescription, -p.adjust),
      x = -log10(p.adjust),
      fill = Count
    )) +
    geom_col()+
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
      maxGSSCutoff = 300
    )
    
    if(nrow(pathwaysGeneSets) > 0){
      p = plotTop10Terms(
        enricherOutput = pathwaysGeneSets, 
        title = paste0(compairsonName, ' - ', s)
      )
      
      ggsave(
        filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqJejunumOutputs/lrtResults/jejunumEnrichmentAnalysis/', compairsonName,'/temporalIncrease/',s, '_', s, '.png'),
        plot = p,
        width = 12, height = 8
      )
      write_tsv(
        pathwaysGeneSets,
        file = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqJejunumOutputs/lrtResults/jejunumEnrichmentAnalysis/', compairsonName,'/temporalIncrease/',s, '_', s, '.tsv')
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
      minGSSCutoff = 10,
      maxGSSCutoff = 500
    )
    
    if(nrow(pathwaysGeneSets) > 0){
      p = plotTop10Terms(
        enricherOutput = pathwaysGeneSets, 
        title = paste0(compairsonName, ' - ', s)
      )
      ggsave(
        filename = paste0('DeSeq2/barseqJejunumOutputs/lrtResults/jejunumEnrichmentAnalysis/', compairsonName,'/temporalDecrease/',s, '_', s, '.png'),
        plot = p,
        width = 12, height = 8
      )
      write_tsv(
        pathwaysGeneSets,
        file = paste0('DeSeq2/barseqJejunumOutputs/lrtResults/jejunumEnrichmentAnalysis/', compairsonName,'/temporalDecrease/',s, '_', s, '.tsv')
      )
    }
  }
  
}
# Define background genes
allGenes = read_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/allTestedGenes.tsv')
allGenes = allGenes$SYMBOL

day1VDay3 = read_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day1Day3.tsv')

runEnrichment(
  deseqResults = day1VDay3,
  universeGenes = allGenes,
  compairsonName = 'day1Vday3'
)

day3VDay7 = read_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day3Day7.tsv')

runEnrichment(
  deseqResults = day3VDay7,
  universeGenes = allGenes,
  compairsonName = 'day3Vday7'
)


day7VDay14 = read_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day7Day14.tsv')

runEnrichment(
  deseqResults = day7VDay14,
  universeGenes = allGenes,
  compairsonName = 'day7Vday14'
)


