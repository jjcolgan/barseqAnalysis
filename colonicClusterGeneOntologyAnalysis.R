# Load required libraries for enrichment analysis, gene annotation, and data wrangling
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)

# Function to plot the top 10 enriched GO terms
plotTop10Terms = function(enricherOutput, 
                          title, 
                          countCutOff = 4, 
                          ...) {
  # Filter GO terms with gene count above the cutoff, sort by adjusted p-value and count
  p = enricherOutput %>%
    filter(Count > countCutOff) %>%
    arrange(p.adjust, desc(Count)) %>%
    head(n = 10) %>% # Take top 10 terms
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
  
  # Display the plot
  plot(p)
}

# Set the minimum cluster size for inclusion in the analysis
clusterSize = 5

# Read the cluster assignments for each gene
clustersDf = read_tsv('DeSeq2/barseqColonicOutputs/manualGeneClusteringClusters.tsv')

# Read the full list of gene symbols for the universe
allGenes = read_tsv('DeSeq2/barseqColonicOutputs/lrtResults/allGenesTab.tsv')
allGenes = allGenes$SYMBOL

# Ensure the first column is named "SYMBOL" for consistency
names(clustersDf)[1] <- "SYMBOL"

# Identify clusters with more than the minimum number of genes
clusters = clustersDf %>%
  distinct() %>%
  group_by(clusterId) %>%
  summarise(nGenes = n()) %>%
  filter(nGenes > clusterSize) %>%
  pull(clusterId)

# Loop through each qualifying cluster to perform GO enrichment and plot results
for(c in clusters){
  # Get the unique gene symbols for this cluster
  genes = clustersDf %>%
    distinct() %>%
    filter(clusterId == c) %>%
    pull(SYMBOL)
    
  # Perform GO enrichment analysis for Biological Process (BP)
  goRes = enrichGO(gene = genes, 
                   OrgDb = org.Mm.eg.db, 
                   universe = allGenes,
                   keyType = "SYMBOL", 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.1, 
                   minGSSize = 25, 
                   maxGSSize = 300) %>%
    clusterProfiler::simplify()
    
  # Plot the top 10 enriched terms
  p = goRes %>% plotTop10Terms(title = paste0("GO Enrichment for Cluster ", c), 
                               enricherOutput = goRes, 
                               countCutOff = 4)
  # Save the full enrichment results table as a TSV file
  goRes %>%
    as.data.frame() %>%
    write_tsv(paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/cluster', c, '/GO_enrichment.tsv'))
  
  # Save the plot to a file
  ggsave(filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/colonicEnrichmentAnalysis/cluster', c, '/GO_enrichment_plot.png'), 
         plot = p, 
         width = 8, 
         height = 6)
}

