' I feel pretty good about this model.'
library(org.Mm.eg.db)
library(clusterProfiler)
library(tidyverse)
library(DESeq2)
library(ashr)
library(ComplexHeatmap)

'Not entirely sure that this will work'

plotRes=function(deSeqResults, comparison, dds, dayList){
  sigRes = deSeqResults %>%
    as.data.frame()%>%
    filter(padj < .1,
           abs(log2FoldChange) > .75)%>%
    arrange(desc(-log10(padj)))
  if (nrow(sigRes) < 100) {
    i = nrow(sigRes)
  }
  else{
    i = 100
  }
  path = paste0('DeSeq2/barseqColonicOutputs/lrtResults/', comparison)
  dir.create(path)
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


geneExRaw<- read_csv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/gene_matrix_count.csv')
metadata <- read.csv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/metadata.csv')
biomart = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/mart_export (1).txt')

colnames(geneExRaw) <-sub("_S\\d+.*$", "", colnames(geneExRaw))
colnames(biomart) = c('Geneid', 'transcriptType')

proteinCoding = biomart %>%
  filter(transcriptType == 'protein_coding')

metadata=metadata%>%
  mutate(tissueDay = paste0(Tissue, Day))

fullDataTab=geneExRaw %>%
  pivot_longer(cols = 2:ncol(geneExRaw),
               names_to = 'library',
               values_to = 'unNormalizedCounts')%>%
  filter(library != 'EC_CG_29')%>%
  left_join(metadata, by = "library")%>%
  filter(Tissue == 'Colon'& Treatment =='Bar-seq')

metaDeSeq=fullDataTab %>%
  select(-c( Geneid, unNormalizedCounts, library))%>%
  distinct()%>%
  #filter(sample != '2308Co')%>%
  column_to_rownames('sample')

'Keeping the sample which has more reads in the duplicated set'
desired_order <- rownames(metaDeSeq)

expressionDeseq=fullDataTab %>%
  filter(Geneid %in% proteinCoding$Geneid)%>%
  pivot_wider(id_cols = 'Geneid',
              values_from = 'unNormalizedCounts',
              names_from = 'sample')%>%
  relocate(Geneid) %>%
  select(Geneid, all_of(desired_order))%>%
  column_to_rownames('Geneid')

geneIds = clusterProfiler::bitr(geneID = rownames(expressionDeseq), fromType = 'ENSEMBL', toType = 'SYMBOL',org.Mm.eg.db, drop = T)

expressionDeseq%>%
  rownames_to_column('ENSEMBL')%>%
  left_join(geneIds , by = "ENSEMBL")%>%
  filter(is.na(SYMBOL)==F)%>%
  dplyr::select(-ENSEMBL)%>%
  group_by(SYMBOL)%>%
  summarise(observations = n())%>%
  filter(observations > 1)


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

metaDeSeq%>%
  group_by(tissueDay)%>%
  summarise(n())

metaDeSeq$Age.at.sac = cut(metaDeSeq$Age.at.sac, 3)
metaDeSeq$Age.at.transplant = cut(metaDeSeq$Age.at.transplant, 3)




'Adding cage results in a model with a dist of pvalues which is much worse than the model
without pvalues.'
dds = DESeqDataSetFromMatrix(countData = expressionDeseq,
                             colData = metaDeSeq,
                             design = ~Day)

'Filtering by group'
normalizedCounts=counts(dds)%>%
  as.data.frame()%>%
  rownames_to_column('gene')
'This is just a more sophisticated version of the deseq2 filter, where the gene needs to be pesent in at least 4 samples for 20 counts. '
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
geneVariance=normalizedCounts%>%
  pivot_longer(cols = 2:ncol(normalizedCounts),
               names_to = 'sample',
               values_to = 'normalizedCount')%>%
  left_join(metadata, by = 'sample')%>%
  filter(gene %in% genesPassingFilter$gene)%>%
  group_by(gene)%>%
  summarise(variance = var(normalizedCount))%>%
  arrange(variance)

lowerBound=nrow(geneVariance)*.05

upperBound = nrow(geneVariance)

genesPassingFilter=geneVariance$gene[lowerBound:upperBound]

'filter low count'
#smallestGroupSize <- 4
#keep <- rowSums(counts(dds) >= 35) >= smallestGroupSize
dds <- dds[rownames(dds)%in%genesPassingFilter,]
dds$Day =  relevel(dds$Day, ref = "1")


dds <- DESeq(dds,test="LRT", reduced = ~1)

plotDispEsts(dds)

vsd = vst(dds, blind = T)
plotPCA(vsd, intgroup = 'Day')

plotPCA(vsd, intgroup = 'Day',returnData = T)

rlogD=rlog(dds)
plotPCA(rlogD, intgroup = 'Day')

resultsNames(dds)

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



day1VDay3%>%
  as.data.frame()%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Day 1 versus day 3')
#ggsave(p,filename = 'DeSeq2/barseqColonicOutputs/lrtResults/histogramOfPvalues1v3.pdf')

day1VDay3%>%
  as.data.frame()%>%
  filter(padj < .1,
         abs(log2FoldChange) > 1.5)%>%
  arrange(desc(-log10(padj)))
'Fold change close to cut off, pvalue is fine.'


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


'I might be using the wrong table to perform my correlation analysis.
I think either the vst dataframe or rlog dataframe might be better '
assay(vsd)%>%
  as.data.frame()%>%
  head()

expressionRlogNorm=assay(rlogD)%>%
  as.data.frame()

expressionVsNorm = assay(vsd)%>%
  as.data.frame()

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


sigAll=rbind(sigDay1Day3, sigDay3Day7, sigDay7Day14)

sigGeneTab=assay(vsd) %>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  filter(SYMBOL %in% rownames(sigAll))

write_tsv(sigGeneTab,
          file = 'DeSeq2/barseqColonicOutputs/lrtResults/sigGeneTab.tsv')

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


rbind(day1VDay3,
      day3VDay7,
      day7VDay14
      )%>%
  as.data.frame()%>%
  filter(padj <.1,
         abs(log2FoldChange) >.75)%>%
  ggplot(aes(x = abs(log2FoldChange)))+
  geom_histogram(bins = 25)

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


