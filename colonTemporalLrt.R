# ' I feel pretty good about this model.'

# Load required libraries for analysis and plotting
library(org.Mm.eg.db)           # Mouse genome annotation
library(clusterProfiler)         # For gene annotation and enrichment analysis
library(tidyverse)               # Data manipulation and visualization
library(DESeq2)                  # RNA-seq differential expression
library(ashr)                    # Adaptive shrinkage (not directly used here)
library(ComplexHeatmap)          # Complex heatmaps

# 'Not entirely sure that this will work'

# Function to plot and save results for significant genes from DESeq2 results
plotRes=function(deSeqResults, comparison, dds, dayList){
  sigRes = deSeqResults %>%
    as.data.frame()%>%
    filter(padj < .1,
           abs(log2FoldChange) > .75)%>%         # Filter significant genes
    arrange(desc(-log10(padj)))
  if (nrow(sigRes) < 100) {
    i = nrow(sigRes)
  }
  else{
    i = 100
  }
  path = paste0('DeSeq2/barseqColonicOutputs/lrtResults/', comparison)
  dir.create(path)                                # Create output directory
  for (j in 1:i){
    temp = sigRes[j,]
    temp = temp%>%
      rownames_to_column('gene')
    gene = temp$gene
    fc = temp$log2FoldChange
    padj = temp$padj

    p = plotCounts(dds, intgroup = 'Day', gene = gene, returnData = T, normalized = T, transform = T)%>%
      as.data.frame()%>%
      filter(Day %in% dayList)%>%
      ggplot(aes(y = count,
                 x = Day))+
      geom_boxplot(outliers = F)+
      geom_jitter()+
      labs(title = paste0(gene , ' ' , comparison),
           caption = paste0('Padj = ', padj, '\nFC = ', fc))
    ggsave(path = path, filename= paste0(gene,'.pdf') ,p, width = 5, height = 5, units = 'in')

  }

}

# Load data files
geneExRaw<- read_csv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/gene_matrix_count.csv')
metadata <- read.csv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/metadata.csv')
biomart = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/mart_export (1).txt')

# Clean up column names and biomart for later filtering
colnames(geneExRaw) <-sub("_S\\d+.*$", "", colnames(geneExRaw))
colnames(biomart) = c('Geneid', 'transcriptType')

# Keep only protein coding genes
proteinCoding = biomart %>%
  filter(transcriptType == 'protein_coding')

# Generate a new column for groupings in metadata
metadata=metadata%>%
  mutate(tissueDay = paste0(Tissue, Day))

# Reshape and join count and metadata; filter for relevant samples (Colon, Bar-seq)
fullDataTab=geneExRaw %>%
  pivot_longer(cols = 2:ncol(geneExRaw),
               names_to = 'library',
               values_to = 'unNormalizedCounts')%>%
  filter(library != 'EC_CG_29')%>%
  left_join(metadata, by = "library")%>%
  filter(Tissue == 'Colon'& Treatment =='Bar-seq')

# Prepare metadata for DESeq2
metaDeSeq=fullDataTab %>%
  select(-c( Geneid, unNormalizedCounts, library))%>%
  distinct()%>%
  #filter(sample != '2308Co')%>%
  column_to_rownames('sample')

# Keep the sample with more reads if duplicated
desired_order <- rownames(metaDeSeq)

# Prepare expression data matrix for DESeq2; keep only protein coding genes and order columns
expressionDeseq=fullDataTab %>%
  filter(Geneid %in% proteinCoding$Geneid)%>%
  pivot_wider(id_cols = 'Geneid',
              values_from = 'unNormalizedCounts',
              names_from = 'sample')%>%
  relocate(Geneid) %>%
  select(Geneid, all_of(desired_order))%>%
  column_to_rownames('Geneid')

# Map ENSEMBL gene IDs to gene symbols using clusterProfiler
geneIds = clusterProfiler::bitr(geneID = rownames(expressionDeseq), fromType = 'ENSEMBL', toType = 'SYMBOL',org.Mm.eg.db, drop = T)

# Check for genes with multiple mappings to symbols
expressionDeseq%>%
  rownames_to_column('ENSEMBL')%>%
  left_join(geneIds , by = "ENSEMBL")%>%
  filter(is.na(SYMBOL)==F)%>%
  dplyr::select(-ENSEMBL)%>%
  group_by(SYMBOL)%>%
  summarise(observations = n())%>%
  filter(observations > 1)

# Collapse multiple ENSEMBL IDs per symbol, summing counts if needed
expressionDeseq=expressionDeseq %>%
  rownames_to_column('ENSEMBL') %>%
  left_join(geneIds, by = "ENSEMBL") %>%
  filter(!is.na(SYMBOL)) %>%
  dplyr::select(-ENSEMBL) %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE),
            across(where(~!is.numeric(.)), ~ first(.)),
            .groups = "drop")%>%
  column_to_rownames('SYMBOL')

# Check group size in metadata
metaDeSeq%>%
  group_by(tissueDay)%>%
  summarise(n())

# Bin age columns into 3 quantile groups
metaDeSeq$Age.at.sac = cut(metaDeSeq$Age.at.sac, 3)
metaDeSeq$Age.at.transplant = cut(metaDeSeq$Age.at.transplant, 3)

# 'Adding cage results in a model with a dist of pvalues which is much worse than the model without pvalues.'

# Create DESeq2 dataset, using only Day as design
dds = DESeqDataSetFromMatrix(countData = expressionDeseq,
                             colData = metaDeSeq,
                             design = ~Day)

# 'Filtering by group'

# Get normalized counts and reshape for filtering
normalizedCounts=counts(dds)%>%
  as.data.frame()%>%
  rownames_to_column('gene')

# Filter genes: must have at least 10 counts in 4+ samples per group
genesPassingFilter=normalizedCounts%>%
  pivot_longer(cols = 2:ncol(normalizedCounts),
               names_to = 'sample',
               values_to = 'normalizedCount')%>%
  left_join(metadata, by = 'sample')%>%
  filter(normalizedCount >= 10)%>%
  group_by(gene,
           Day)%>%
  summarise(nPass = n())%>%
  filter(nPass >= 4)%>%
  select(gene)%>%
  distinct()

# Calculate variance for each gene and sort
geneVariance=normalizedCounts%>%
  pivot_longer(cols = 2:ncol(normalizedCounts),
               names_to = 'sample',
               values_to = 'normalizedCount')%>%
  left_join(metadata, by = 'sample')%>%
  filter(gene %in% genesPassingFilter$gene)%>%
  group_by(gene)%>%
  summarise(variance = var(normalizedCount))%>%
  arrange(variance)

# Keep top 95% variable genes
lowerBound=nrow(geneVariance)*.05
upperBound = nrow(geneVariance)
genesPassingFilter=geneVariance$gene[lowerBound:upperBound]

# Filter low count genes in dds
dds <- dds[rownames(dds)%in%genesPassingFilter,]
dds$Day =  relevel(dds$Day, ref = "1")   # Set reference level for Day

# Run DESeq2 using likelihood ratio test (LRT)
dds <- DESeq(dds,test="LRT", reduced = ~1)

# Plot dispersion estimates
plotDispEsts(dds)

# Perform transformations for PCA and visualization
vsd = vst(dds, blind = T)
plotPCA(vsd, intgroup = 'Day')
plotPCA(vsd, intgroup = 'Day',returnData = T)

rlogD=rlog(dds)
plotPCA(rlogD, intgroup = 'Day')

resultsNames(dds) # Show available contrasts

# Conduct pairwise comparisons between days and plot results
day1VDay3 = results(dds, contrast = c('Day','1', '3'), alpha = .1)
plotMA(day1VDay3)
p=day1VDay3%>%
  as.data.frame()%>%
  na.omit()%>%
  ggplot(aes(x = log2FoldChange,
             col = padj < .1 & abs(log2FoldChange)> .75,
             y = -log10(padj)))+
  geom_point()+
  labs(title = 'day 1 v day3', 
       col = 'Significant')+
  xlim(-5,5)
ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/volcanoPlot1v3.pdf')

# Plot histogram of p-values for visual QC
day1VDay3%>%
  as.data.frame()%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Day 1 versus day 3')
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/histogramOfPvalues1v3.pdf')

# Show top significant results
day1VDay3%>%
  as.data.frame()%>%
  filter(padj < .1,
         abs(log2FoldChange) > 1.5)%>%
  arrange(desc(-log10(padj)))
# 'Fold change close to cut off, pvalue is fine.'

# Repeat for other day comparisons
day1VDay7 = results(dds, contrast = c('Day','1', '7'), alpha = .1)
plotMA(day1VDay7)
day1VDay7%>%
  as.data.frame()%>%
  ggplot(aes(x = log2FoldChange,
             col = padj < .1,
             y = -log10(padj),
  ))+
  geom_point()+
  labs(title = 'day 1 v day7')+
  xlim(-5,5)
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/volcanoPlot1v7.pdf')

day1VDay7%>%
  as.data.frame()%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Day 1 versus day 7')
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/histogramOfPvalues1v7.pdf')

day1VDay7%>%
as.data.frame()%>%
filter(padj < .1,
abs(log2FoldChange) > 1.5)%>%
arrange(desc(-log10(padj)))

day1VDay14 = results(dds, contrast = c('Day','1', '14'), alpha = .1)
plotMA(day1VDay14)

day1VDay14%>%
as.data.frame()%>%
ggplot(aes(x = pvalue))+
geom_histogram(bins = 100)+
  labs(title = 'day 1 versus day 14')
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/histogramOfPvalues1v14.pdf')

day1VDay14%>%
  as.data.frame()%>%
  ggplot(aes(x = log2FoldChange,
             col = padj < .1,
             y = -log10(padj)))+
  geom_point()+
  labs(title = 'day 1 v day14')+
  xlim(-5,5)
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/volcanoPlot1v14.pdf')

day1VDay14 %>%
  as.data.frame()%>%
  filter(padj < .1,
         abs(log2FoldChange) > 1.5)%>%
  arrange(desc(-log10(padj)))

day3VDay7 = results(dds, contrast = c('Day','3', '7'))
plotMA(day3VDay7)

day3VDay7%>%
  as.data.frame()%>%
  na.omit()%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'day 3 v day 7')
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/histogramOfPvalues3v7.pdf')

 p=day3VDay7%>%
   na.omit()%>%
  as.data.frame()%>%
  ggplot(aes(x = log2FoldChange,
             col = padj < .1 & abs(log2FoldChange) > .75,
             y = -log10(padj)))+
  geom_point()+
  labs(title = 'day 3 v day 7', 
       col = "Significant")
ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/volcanoPlot3v7.pdf')

day3VDay7%>%
  as.data.frame()%>%
  filter(padj < .1,
         abs(log2FoldChange) > 1.5)%>%
  arrange(desc(-log10(padj)))

day7VDay14 = results(dds, contrast = c('Day','7', '14'))
plotMA(day7VDay14)
p=day7VDay14%>%
  as.data.frame()%>%
  na.omit()%>%
  ggplot(aes(x = log2FoldChange,
             col = padj < .1 & abs(log2FoldChange) > .75,
             y = -log10(padj)))+
  geom_point()+
  labs(title = 'day 7 v day 14', 
       col = 'Significant')
ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/volcanoPlot7v14.pdf')

day7VDay14%>%
as.data.frame()%>%
ggplot(aes(x = pvalue))+
geom_histogram(bins = 100)
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/histogramOfPvalues7v14.pdf')

day7VDay14%>%
  as.data.frame()%>%
  filter(padj < .1,
         abs(log2FoldChange) > 1.5)%>%
  arrange(desc(-log10(padj)))

as.data.frame(day7VDay14)[1,]

day3VDay14  = results(dds, contrast = c('Day','3', '14'))

p = day3VDay14%>%
as.data.frame()%>%
ggplot(aes(x = pvalue))+
geom_histogram(bins = 100)
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/histogramOfPvalues3v14.pdf')

day3VDay14%>%
  as.data.frame()%>%
  ggplot(aes(x = log2FoldChange,
             col = padj < .1,
             y = -log10(padj)))+
  geom_point()+
  labs(title = 'day 3 v day 14')
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/volcanoPlot7v14.pdf')

# 'I might be using the wrong table to perform my correlation analysis.
# I think either the vst dataframe or rlog dataframe might be better '

# Preview transformed data for correlation
assay(vsd)%>%
  as.data.frame()%>%
  head()

expressionRlogNorm=assay(rlogD)%>%
  as.data.frame()

expressionVsNorm = assay(vsd)%>%
  as.data.frame()

# Plot and save top results for each comparison using the earlier function
plotRes(deSeqResults = day7VDay14,
        comparison = 'Day7VDay14',
        dayList = c('7','14'),
        dds = dds)

plotRes(deSeqResults = day1VDay3,
        comparison = 'Day1VDay3',
        dayList = c('1', '3'),
        dds = dds)

plotRes(deSeqResults = day3VDay7,
        comparison = 'Day3VDay7',
        dayList = c('3', '7'),
        dds = dds)

plotRes(deSeqResults = day1VDay7,
        comparison = 'Day1VDay7',
        dayList = c('1', '7'),
        dds = dds)

plotRes(deSeqResults = day1VDay14,
        comparison = 'Day1VDay14',
        dayList = c('1', '14'),
        dds = dds)

# Save all results tables to TSV files for further analysis
day1VDay3%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqColonicOutputs/lrtResults/Day1VDay3/deseqResults.tsv')

day1VDay7%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqColonicOutputs/lrtResults/Day1VDay7/deseqResults.tsv')

day1VDay14%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqColonicOutputs/lrtResults/Day1VDay14/deseqResults.tsv')

day3VDay7%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqColonicOutputs/lrtResults/Day3VDay7/deseqResults.tsv')

day7VDay14%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqColonicOutputs/lrtResults/Day7VDay14/deseqResults.tsv')

# Extract significant genes for each comparison (abs(log2FC) > .75, padj < .1)
sigDay1Day3=day1VDay3%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange) >.75 & padj < .1)

sigDay3Day7=day3VDay7%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange) >.75 & padj < .1)

sigDay7Day14=day7VDay14%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange) >.75 & padj < .1)

sigDay1Day7 = day1VDay7%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange) >.75 & padj < .1)

sigDay1Day14 = day1VDay14%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange) >.75 & padj < .1)

sigDay3Day14 = day3VDay14%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange) >.75 & padj < .1)

# Combine all significant genes from a subset of comparisons
sigAll=rbind(sigDay1Day3, sigDay3Day7, sigDay7Day14)

# Extract normalized counts for significant genes for downstream analysis
sigGeneTab=assay(vsd) %>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  filter(SYMBOL %in% rownames(sigAll))

write_tsv(sigGeneTab,
          file = 'DeSeq2/barseqColonicOutputs/lrtResults/sigGeneTab.tsv')

# Save additional result tables for each comparison
day1VDay3%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv(file = 'DeSeq2/barseqColonicOutputs/lrtResults/day1Vday3.tsv')

day3VDay7%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv(file = 'DeSeq2/barseqColonicOutputs/lrtResults/day3Vday7.tsv')

day7VDay14%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv(file = 'DeSeq2/barseqColonicOutputs/lrtResults/day7Vday14.tsv')

assay(vsd) %>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv(file = 'DeSeq2/barseqColonicOutputs/lrtResults/allGenesTab.tsv')

# Plot histogram of log2 fold changes for all significant genes
rbind(day1VDay3,
      day3VDay7,
      day7VDay14
      )%>%
  as.data.frame()%>%
  filter(padj <.1,
         abs(log2FoldChange) >.75)%>%
  ggplot(aes(x = abs(log2FoldChange)))+
  geom_histogram(bins = 25)

# Get integration genes with strong changes across any comparison
integrationGenes=rbind(day1VDay3,
      day3VDay7,
      day7VDay14
)%>%
  as.data.frame()%>%
  filter(padj <.1,
         abs(log2FoldChange) >1)%>%
  rownames_to_column('gene')%>%
  select(gene)%>%
  distinct()%>%
  .$gene

assay(vsd) %>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  filter(SYMBOL %in% integrationGenes)%>%
  write_tsv(file = 'DeSeq2/barseqColonicOutputs/lrtResults/integrationGenes.tsv')

