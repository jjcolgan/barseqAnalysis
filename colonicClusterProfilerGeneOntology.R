# colonicClusterProfilerGeneOntology.R
# This script performs Gene Ontology (GO) analysis on colonic gene expression data
# using clusterProfiler with a specific background gene set for colonic samples.
# NOTE: This is distinct from other gene ontology analyses in this repository, as it uses a separate
# background database (colonic gene set), not the default or previously used database.

library(clusterProfiler) # For GO enrichment analysis
library(org.Mm.eg.db)    # Mouse gene annotation database
library(tidyverse)       # For data manipulation and plotting

# Function to plot the top 10 GO terms from clusterProfiler output
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

# Function to run clusterProfiler's enrichGO on a given gene list and background
runClusterProfiler = function(geneList, 
                              backgroundGenes,
                              pvalueCutOff = .1, 
                              minGSSize = 10, 
                              maxGSSize = 500, 
                              ...){
    # Note: Using colonic-specific gene universe as backgroundGenes argument!
    goRes = enrichGO(geneList, 
             OrgDb = org.Mm.eg.db, 
             ont= 'BP',
             keyType = 'SYMBOL', 
             pAdjustMethod = 'fdr', 
             pvalueCutoff = pvalueCutOff, 
             qvalueCutoff = pvalueCutOff,
             universe = backgroundGenes)
    
    if(is.null(goRes) || nrow(goRes) == 0){
      return(NULL)
    } else {
      goRes = goRes %>%
         clusterProfiler::simplify()
      return(as.data.frame(goRes))
    }
}

# Main wrapper to run GO analysis for up- and down-regulated genes in a comparison
runGO = function(inputDf, 
                 backgroundGenes,
                 pvalueCutOff = .1, 
                 foldChangeCutOff = .75,
                 comparison,
                 ...){
  # Identify significantly changed genes using p-value and fold change cutoffs
  sigGenes = inputDf %>%
    filter(is.na(padj) | padj < pvalueCutOff,
           abs(log2FoldChange) > foldChangeCutOff)
  
  # Genes increased in expression
  increasingGenes = sigGenes %>%
    filter(log2FoldChange < 0) %>%
    pull(SYMBOL) 
  
  # Genes decreased in expression
  decreasingGenes = sigGenes %>%
    filter(log2FoldChange > 0) %>%
    pull(SYMBOL) 
  
  # Run GO analysis for increased genes
  increasingGo = runClusterProfiler(increasingGenes, 
                   backgroundGenes = backgroundGenes,
                   pvalueCutOff = pvalueCutOff, 
                   ...)
  if(!is.null(increasingGo)){
    p = plotTop10Terms(increasingGo, 
                      title = paste0('GO Enrichment for ', comparison, ' Increasing Genes'))
    ggsave(
      filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', comparison, '_increasing.pdf'),
      plot = p,
      width = 8, height = 6
    )
    write_tsv(increasingGo, 
              paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', comparison, '_increasing.tsv'))
  }
  
  # Run GO analysis for decreased genes
  decreasingGo = runClusterProfiler(decreasingGenes, 
                   backgroundGenes = backgroundGenes,
                   pvalueCutOff = pvalueCutOff, 
                   ...)
  if(!is.null(decreasingGo)){
    p = plotTop10Terms(decreasingGo, 
                      title = paste0('GO Enrichment for ', comparison, ' Decreasing Genes'))
    ggsave(
      filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', comparison, '_decreasing.pdf'),
      plot = p,
      width = 8, height = 6
    )
    write_tsv(decreasingGo, 
              paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/', comparison, '_decreasing.tsv'))
  }
}

# --- Main Script Execution ---

# Read in the colonic background gene universe (distinct from other scripts!)
allGenes = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/allGenesTab.tsv')
allGenes = allGenes$SYMBOL

# Run analysis on day1 vs day3 comparison
day1VDay3 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day1Vday3.tsv')
runGO(day1VDay3, 
       backgroundGenes = allGenes, 
       pvalueCutOff = .1, 
       foldChangeCutOff = .75, 
       comparison = 'day1Vday3')

# Run analysis on day3 vs day7 comparison
day3VDay7 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day3Vday7.tsv')
runGO(day3VDay7, 
       backgroundGenes = allGenes, 
       pvalueCutOff = .1, 
       foldChangeCutOff = .75, 
       comparison = 'day3Vday7')

# Run analysis on day7 vs day14 comparison
day7VDay14 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day7Vday14.tsv')
runGO(day7VDay14, 
       backgroundGenes = allGenes, 
       pvalueCutOff = .1, 
       foldChangeCutOff = .75, 
       comparison = 'day7Vday14')

runGO(day7VDay14, 
       backgroundGenes = allGenes, 
       pvalueCutOff = .1, 
       foldChangeCutOff = .75, 
       comparison = 'day7Vday14')
