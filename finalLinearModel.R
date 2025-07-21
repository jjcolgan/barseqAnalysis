library(tidyverse)
library(lmerTest)

# Simple QQ plot function
simpleQQPlot = function (observedPValues) {
  plot(-log10(1:length(observedPValues)/length(observedPValues)),
       -log10(sort(observedPValues)))
  abline(0, 1, col = "red")
}

calculateMeanDeltaFitness = function(inputDf){
  groupMeans = inputDf %>%
    group_by(tissue) %>%
    summarise(groupMean = mean(fitness))
  absoluteDelta = abs(groupMeans$groupMean[1]-groupMeans$groupMean[2])
  return(absoluteDelta)
}

# Load data
fitnessScores = read_tsv('barseqAdjustedParams/fit_logratios.tab')
metadata = read_tsv('fullbarseqMeta.txt')
annotations = read_tsv('genesWithAnvioAnnotations.tsv')

# Format metadata
metadata=metadata%>%
  filter(tissue != 'T0')
metadata$day = factor(metadata$day, levels = c('t0', 'day1', 'day3', 'day7', 'day14'))
metadata$phase = ifelse(metadata$dayNumeric < 7, 'early', 'late')
metadata$cage = as.factor(metadata$cage)

# KEGG annotations
keggs = annotations %>%
  select(locusId, kofamAccession, kofamFunction) %>%
  distinct()

# Clean up fitness score sample names
colnames(fitnessScores) <- sub("setA", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("_.*", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("CO$", "Co", colnames(fitnessScores))
colnames(fitnessScores) <- sub("DJ$", "Dj", colnames(fitnessScores))

# Subset metadata by tissue
colonMeta = metadata %>%
  filter(tissue == 'colon') %>%
  column_to_rownames('sample')

djMeta = metadata %>%
  filter(tissue == 'dj') %>%
  column_to_rownames('sample')

# Prepare input for linear modeling
lmIn = fitnessScores %>%
  filter(desc != 'hypothetical protein',
         grepl("^tRNA-[A-Za-z]{3}$", desc)!=T)%>%
  select(-c(sysName, desc)) %>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  filter(sample %in% metadata$sample)%>%
  column_to_rownames('sample')%>%
  t()

loci = rownames(lmIn)

#plot PCA to tune filters
var <-.5
lmIn_filtered=lmIn[apply(lmIn, 1, var) >= var, ]
pcaOut=prcomp(t(lmIn_filtered))
pcaSummary=summary(pcaOut)
pcaSummary$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = "sample")%>%
  ggplot(aes(x = PC1,
             y = PC2,
             shape = day,
             col = tissue,
             ))+
  geom_point()+
  stat_ellipse()

loci = rownames(lmIn_filtered)

# Fit GLM for each gene
output = data.frame(gene = character(), 
                    pvalue = double(), 
                    coef = double(), 
                    absoluteDelta = double())

for (l in seq_along(loci)) {
  input = lmIn[rownames(lmIn) == loci[l], ] %>%
    as.data.frame()%>%
    rownames_to_column('sample')%>%
    rename(fitness = c(2))%>%
    merge(metadata, by = 'sample')
  input$tissue = as.factor(input$tissue)

  #glm = glm(fitness ~ tissue + percentPerfectBarcode + cage+millionBases, data = input)
  
  lmer_fit <- lmer(
    fitness ~ tissue + millionBases + (1|lane:cage:mouse),
    data = input,
    REML = FALSE
  )
  # glmRes = summary(lmer_fit)
  # coef = glmRes$coefficients[2, 1]
  # p = glmRes$coefficients[2, 4]
  # absoluteDelta = calculateMeanDeltaFitness(input)
  
  coefs <- summary(lmer_fit)$coefficients
  p <- coefs[2, "Pr(>|t|)"]
  coef = coefs[2,1]
  absoluteDelta = calculateMeanDeltaFitness(input)

  output = rbind(output, data.frame(gene = loci[l], 
                                    pvalue = p, coef = coef,
                                    absoluteDelta = absoluteDelta))
}

# Multiple testing correction
output$padj = p.adjust(output$pvalue, method = 'fdr')
summary(output$padj)

# Significant hits
sig = output %>%
  filter(padj < 0.1) %>%
  arrange(padj)

# Plotting
output %>%
  ggplot(aes(x = pvalue)) +
  geom_histogram(bins = 100) +
  labs(title = 'Histogram of unadjusted p-values') +
  theme_bw()

output %>%
  ggplot(aes(x = pvalue)) +
  geom_density()

simpleQQPlot(output$pvalue)

# Mean tissue fitness per gene
tissueScatterIn = fitnessScores %>%
  pivot_longer(cols = 4:ncol(.), names_to = 'sample', values_to = 'fitnessScore') %>%
  left_join(metadata, by = 'sample') %>%
  filter(tissue != 'T0') %>%
  group_by(locusId, tissue) %>%
  summarise(meanTissueFitness = mean(fitnessScore), .groups = 'drop') %>%
  pivot_wider(id_cols = locusId, names_from = tissue, values_from = meanTissueFitness)

# Color significant genes
tissueScatterIn$col = 'P-adjust > .1'
tissueScatterIn$col[tissueScatterIn$locusId %in% sig$gene] = 'P-adjust < .1'

# Scatter plot of mean fitness by tissue
tissueScatterIn %>%
  ggplot(aes(x = dj, y = colon, col = col)) +
  geom_point(alpha = 0.15) +
  geom_abline() +
  theme_bw() +
  labs(col = 'Significance')

# Annotate significant genes with KEGG
sig = sig %>%
  rename(locusId = gene) %>%
  left_join(keggs, by = 'locusId')

# Save results
write_tsv(output, 'linear models/lmTissueComparisionsAllTimePoints/linearModelRes/linearModelRes.tsv')
