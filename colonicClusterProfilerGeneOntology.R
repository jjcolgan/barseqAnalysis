library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)

plotTop10Terms = function(enricherOutput, 
                          countCutOff = 4, 
                          title = "Top 10 GO Terms",
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

runClusterProfiler = function(geneList, 
                              backgroundGenes,
                              pvalueCutOff = .1, 
                              minGSSize = 10, 
                              maxGSSize = 500, 
                            ...){
    goRes=enrichGO(geneList, 
             OrgDb = org.Mm.eg.db, 
             ont= 'BP',
             keyType = 'SYMBOL', 
             pAdjustMethod = 'fdr', 
             pvalueCutoff = pvalueCutOff, 
             qvalueCutoff = pvalueCutOff,
             universe = backgroundGenes)
    
    if(is.null(goRes) || nrow(goRes) == 0){
      return(NULL)
    }else{
      goRes=goRes %>%
         clusterProfiler::simplify()
      return(as.data.frame(goRes))
    }
}

runGO = function(inputDf, 
                 backgroundGenes,
                 pvalueCutOff = .1, 
                 foldChangeCutOff = .75,
                 comparison,
                 ...){
  sigGenes=inputDf%>%
    filter(is.na(padj) | padj < pvalueCutOff,
           abs(log2FoldChange) > foldChangeCutOff)
  
  increasingGenes = sigGenes%>%
    filter(log2FoldChange < 0)%>%
    pull(SYMBOL) 
  
  decreasingGenes = sigGenes%>%
    filter(log2FoldChange > 0)%>%
    pull(SYMBOL) 
  
  increasingGo=runClusterProfiler(increasingGenes, 
                   backgroundGenes = backgroundGenes,
                   pvalueCutOff = pvalueCutOff, 
                   ...)
  if(!is.null(increasingGo)){
    p=plotTop10Terms(increasingGo, 
                     title = paste0('GO Enrichment for ', comparison, ' Increasing Genes'))
    ggsave(
      filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', comparison,'/temporalIncrease/', '_goBp.png'),
      plot = p,
      width = 8, height = 6
    )
    write_tsv(increasingGo, 
              paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', comparison, '/temporalIncrease/', '_goBp.tsv'))
  }
  
  decreasingGo=runClusterProfiler(decreasingGenes, 
                   backgroundGenes = backgroundGenes,
                   pvalueCutOff = pvalueCutOff, 
                   ...)
  if(!is.null(decreasingGo)){
    p=plotTop10Terms(decreasingGo, 
                     title = paste0('GO Enrichment for ', comparison, ' Decreasing Genes'))
    ggsave(
      filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', comparison, '/temporalDecrease/', '_goBp.png'),
      plot = p,
      width = 8, height = 6
    )
    write_tsv(decreasingGo, 
              paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', comparison, '/temporalDecrease/', '_goBp.tsv'))
  }
}

allGenes = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/allGenesTab.tsv')
allGenes = allGenes$SYMBOL
day1VDay3 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day1Vday3.tsv')

runGO(day1VDay3, 
       backgroundGenes = allGenes, 
       pvalueCutOff = .1, 
       foldChangeCutOff = .75, 
       comparison = 'day1Vday3')
day3VDay7 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day3Vday7.tsv')
runGO(day3VDay7, 
       backgroundGenes = allGenes, 
       pvalueCutOff = .1, 
       foldChangeCutOff = .75, 
       comparison = 'day3Vday7')

day7VDay14 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day7Vday14.tsv')
runGO(day7VDay14, 
       backgroundGenes = allGenes, 
       pvalueCutOff = .1, 
       foldChangeCutOff = .75, 
       comparison = 'day7Vday14')