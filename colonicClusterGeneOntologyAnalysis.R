library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)

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

#'Setting min cluster size for consideration. 
clusterSize = 5
clustersDf = read_tsv('DeSeq2/barseqColonicOutputs/manualGeneClusteringClusters.tsv')

allGenes = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/allGenesTab.tsv')
allGenes = allGenes$SYMBOL

names(clustersDf)[1] <- "SYMBOL"


clusters=clustersDf%>%
  distinct()%>%
  group_by(clusterId)%>%
  summarise(nGenes = n())%>%
  filter(nGenes > clusterSize)%>%
  pull(clusterId)

for(c in clusters){
  genes = clustersDf %>%
    distinct()%>%
    filter(clusterId == c) %>%
    pull(SYMBOL)
    
    goRes= enrichGO(gene = genes, 
             OrgDb = org.Mm.eg.db, 
             universe = allGenes,
             keyType = "SYMBOL", 
             ont = "BP", 
             pAdjustMethod = "BH", 
             pvalueCutoff = 0.1, 
             minGSSize = 25, 
             maxGSSize = 300)%>%
      clusterProfiler::simplify()
    
    p = goRes%>% plotTop10Terms(title = paste0("GO Enrichment for Cluster ", c), 
                                 enricherOutput = goRes, 
                                 countCutOff = 4)
  goRes%>%
    as.data.frame()%>%
    write_tsv(paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/cluster', c, '/goEnrichment.tsv'))
  
  ggsave(filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/cluster', c, '/geneOntology.png'), 
         plot = p, 
         width = 8, 
         height = 6)
  
}


