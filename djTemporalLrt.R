library(org.Mm.eg.db)
library(clusterProfiler)
library(tidyverse)
library(DESeq2)
library(ashr)

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
  path = paste0('DeSeq2/barseqJejunumOutputs/lrtResults/', comparison)
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
  filter(Tissue == 'Jejunum'& Treatment =='Bar-seq')

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

geneIds = clusterProfiler::bitr(geneID = rownames(expressionDeseq),
                                fromType = 'ENSEMBL',
                                toType = 'SYMBOL',
                                org.Mm.eg.db,
                                drop = T)

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


expressionDeseq%>%
  as.matrix()%>%
  MatrixGenerics::rowVars()%>%
  as.data.frame()%>%
  rename(geneVar = c(1))%>%
  ggplot(aes(x = geneVar))+
  geom_density()

expressionDeseq%>%
  as.matrix()%>%
  MatrixGenerics::rowVars()%>%
  as.data.frame()%>%
  rename(geneVar = c(1))%>%
  summary(.$geneVar)



'Adding cage results in a model with a dist of pvalues which is much worse than the model
without pvalues.'
dds = DESeqDataSetFromMatrix(countData = expressionDeseq,
                             colData = metaDeSeq,
                             design = ~Day)

metaDeSeq %>%
  group_by(Day)%>%
  summarise(n())

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

upperBound=nrow(geneVariance)

genesPassingFilter=geneVariance$gene[lowerBound:upperBound]

'filter low count'
#smallestGroupSize <- 4
#keep <- rowSums(counts(dds) >= 35) >= smallestGroupSize
dds <- dds[rownames(dds)%in%genesPassingFilter,]
dds$Day =  relevel(dds$Day, ref = "1")


dds <- DESeq(dds,test="LRT", reduced = ~1)
plotDispEsts(dds)


vsd = vst(dds, blind = F)
plotPCA(vsd, intgroup = 'Day', ntop = 1000)

plotPCA(vsd, intgroup = 'Day', return =T)

day1VDay3 = results(dds, contrast = c('Day','1', '3'), alpha = .1)
plotMA(day1VDay3)
day1VDay3%>%
  as.data.frame()%>%
  ggplot(aes(x = log2FoldChange,
             col = padj < .1 & abs(log2FoldChange) > .75,
             y = -log10(padj)))+
  geom_point()+
  labs(title = 'day 1 v day3', 
       col = 'Significant')+
  xlim(-5,5)

day1VDay3%>%
  as.data.frame()%>%
  nrow()


day1VDay3%>%
  as.data.frame()%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Day 1 versus day 3')

day3VDay7 = results(dds, contrast = c('Day','3', '7'), alpha = .1)

day3VDay7%>%
  as.data.frame()%>%
  na.omit()%>%
  ggplot(aes(x = log2FoldChange,
             col = padj < .1 & abs(log2FoldChange) > .75,
             y = -log10(padj)))+
  geom_point()+
  labs(title = 'day 3 v day 7', 
       col = "Significant")+
  xlim(-5,5)

day3VDay7%>%
  as.data.frame()%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Day 3 versus day 7')

day7VDay14 = results(dds, contrast = c('Day','7', '14'), alpha = .1)

day7VDay14%>%
  as.data.frame()%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Day 7 versus day 14')

day7VDay14%>%
  as.data.frame()%>%
  na.omit()%>%
  ggplot(aes(x = log2FoldChange,
             col = padj < .1 & abs(log2FoldChange) > .75, 
             y = -log10(padj)))+
  geom_point()+
  labs(title = 'day 7 v day 14', 
       col = 'Significant')+
  xlim(-5,5)

plotRes(deSeqResults = day1VDay3,
        comparison = 'day1VDay3',
        dds,
        dayList = c('1', '3'))

plotRes(deSeqResults = day3VDay7,
        comparison = 'day3VDay7',
        dds,
        dayList = c('3', '7'))

plotRes(deSeqResults = day7VDay14,
        comparison = 'day7VDay14',
        dds,
        dayList = c('7', '14'))

day1VDay7 = results(dds, contrast = c('Day','1', '7'), alpha = .1)

day1VDay14 = results(dds, contrast = c('Day','1', '14'), alpha = .1)

day1VDay7%>%
  as.data.frame()%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Day 1 versus day 7')

day3VDay14 = results(dds, contrast = c('Day','3', '14'), alpha = .1)


plotRes(deSeqResults = day1VDay7,
        comparison = 'day1VDay7',
        dds,
        dayList = c('1', '7'))

day1VDay14%>%
  as.data.frame()%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Day 1 versus day 14')

plotRes(deSeqResults = day1VDay14,
        comparison = 'day1VDay14',
        dds,
        dayList = c('1', '14'))

day1VDay3%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day1VDay3/deseqResults.tsv')

day3VDay7%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day3VDay7/deseqResults.tsv')

day7VDay14%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day7VDay14/deseqResults.tsv')

sigDay1Day7=day1VDay7%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange)> .75 & padj < .1 )%>%
  rownames_to_column('SYMBOl')

sigDay1Day14=day1VDay7%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange)> .75 & padj < .1 )%>%
  rownames_to_column('SYMBOl')

sigDay3Day7=day3VDay7%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange)> .75 & padj < .1 )%>%
  rownames_to_column('SYMBOl')

sigDay3Day14 = day3VDay14%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange)> .75 & padj < .1 )%>%
  rownames_to_column('SYMBOl')

sigDay7Day14=day7VDay14%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange)> .75 & padj < .1 )%>%
  rownames_to_column('SYMBOl')

sigDay1Day3=day1VDay3%>%
  as.data.frame()%>%
  filter(abs(log2FoldChange)> .75 & padj < .1 )%>%
  rownames_to_column('SYMBOl')

sigAll=rbind(sigDay1Day3,
              sigDay3Day7,
             sigDay7Day14)

vsd %>%
  assay()%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  filter(SYMBOL %in% sigAll$SYMBOl)%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/significantGenes.tsv')

sigDay1Day3%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/sigDay1Day3.tsv')

sigDay1Day7%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/sigDay1Day7.tsv')

sigDay1Day14%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/sigDay1Day14.tsv')

sigDay3Day7%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/sigDay3Day7.tsv')

sigDay3Day14%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/sigDay3Day14.tsv')

sigDay7Day14%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/sigDay7Day14.tsv')

day1VDay3%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day1Day3.tsv')

day1VDay7%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day1Day7.tsv')

day1VDay14%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day1Day14.tsv')

day3VDay7%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day3Day7.tsv')

day3VDay14%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day3Day14.tsv')

day7VDay14%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/day7Day14.tsv')


vsd %>%
  assay()%>%
  as.data.frame()%>%
  rownames_to_column('SYMBOL')%>%
  write_tsv('DeSeq2/barseqJejunumOutputs/lrtResults/allTestedGenes.tsv')