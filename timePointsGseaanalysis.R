library(clusterProfiler)
library(msigdbr)
library(tidyverse)

# Function to get GSEA results using revised ranking score
getMsigDbGSEA = function(rankedGenes, 
                         category, 
                         subCategory, 
                         minGSSCutoff = 0,
                         maxGSSCutoff = 1000, 
                         pvalueCutoff = .1,
                         ...) {
  
  geneSet = msigdbr(
    species = 'mouse', 
    collection = category,
    subcollection = subCategory
  ) %>%
    select(gs_name, gene_symbol)
  
  gseaResult = GSEA(
    geneList = rankedGenes,
    TERM2GENE = geneSet,
    minGSSize = minGSSCutoff,
    maxGSSize = maxGSSCutoff,
    pvalueCutoff = pvalueCutoff,
    verbose = FALSE,
    ...
  ) %>%
    as.data.frame()
  
  return(gseaResult)
}


# Plotting function
plotTopGSEA = function(gseaOutput, 
                       title, 
                       topN = 10) {
  p = gseaOutput %>%
    arrange(p.adjust) %>%
    slice_head(n = topN) %>%
    ggplot(aes(y = reorder(Description, NES),
               x = NES,
               fill = p.adjust)) +
    geom_col() +
    scale_fill_gradient(low = "red", high = "darkblue") +
    labs(title = title, x = 'Normalized Enrichment Score (NES)', y = 'Pathway', fill = 'padj')
  
  plot(p)
}

# Wrapper to run GSEA with log2FC * -log10(padj) score
runGSEA = function(deseqResults, 
                   compairsonName,
                   subCollections = c('PID','REACTOME','KEGG_LEGACY','KEGG_MEDICUS'),
                   pvalueCutoff = 0.1) {
  
  # Build ranked gene list
  rankedGeneList = deseqResults %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    arrange(desc(log2FoldChange)) %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    #'score set like this rather than using just pvalue or just foldchange to avoid
    #'issues where there are genes with good pvalues and small foldchanges
    mutate(score = -log10(padj + 1e-10) * log2FoldChange) %>%
    select(SYMBOL, score) %>%
    deframe() %>%
    sort(decreasing = TRUE)
  
  # Loop through each subcollection
  for (s in subCollections) {
    gseaRes = getMsigDbGSEA(
      rankedGenes = rankedGeneList,
      category = 'C2',
      subCategory = s,
      pvalueCutoff = pvalueCutoff
    )
    
    if (nrow(gseaRes) > 0) {
      p = plotTopGSEA(
        gseaOutput = gseaRes,
        title = paste0(compairsonName, ' - ', s)
      )
      
      outdir = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicGSEA/', compairsonName, '/', s, '/')
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      ggsave(
        filename = paste0(outdir, s, '.png'),
        plot = p,
        width = 8, height = 6
      )
      write_tsv(
        gseaRes,
        file = paste0(outdir, s, '.tsv')
      )
    }
  }
}

allGenes = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/allGenesTab.tsv')$SYMBOL

day1VDay3 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day1Vday3.tsv')
runGSEA(day1VDay3, 'day1Vday3')

day3VDay7 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day3Vday7.tsv')
runGSEA(day3VDay7, 'day3Vday7')

day7VDay14 = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/day7Vday14.tsv')
runGSEA(day7VDay14, 'day7Vday14')

