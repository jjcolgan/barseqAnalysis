library(tidyverse)
library(lme4)
'this might be broken'
run_lm <- function(inputMatrix) {
  loci <- rownames(inputMatrix)

  modelStatsCategoricalDay <- data.frame(
    gene = character(0),
    pvalue = numeric(0),
    coef = numeric(0),
    absoluteDelta = numeric(0)
  )

  for (l in seq_along(loci)) {
    input <- inputMatrix[rownames(inputMatrix) == loci[l], ] %>%
      pivot_longer(
        cols = everything(),
        names_to = "sample",
        values_to = "fitness"
      ) %>%
      merge(djMeta, by = "sample") %>%
      filter(tissue != "T0")

    input$millionBases <- scale(input$millionBases)

    lm_fit <- lm(
      fitness ~ day + percentPerfectBarcode + millionBases,
      data = input
    )

    coefs <- summary(lm_fit)$coefficients
    if (nrow(coefs) >= 2) {
      pvalue <- coefs[2, "Pr(>|t|)"]
      coef   <- coefs[2, "Estimate"]
    } else {
      pvalue <- NA
      coef   <- NA
    }
    absoluteDelta <- calculateMeanDeltaFitness(input)

    modelStatsCategoricalDay <- rbind(
      modelStatsCategoricalDay,
      data.frame(
        gene = loci[l],
        pvalue = pvalue,
        coef = coef,
        absoluteDelta = absoluteDelta
      )
    )
  }

  return(modelStatsCategoricalDay)
}

calculateMeanDeltaFitness = function(inputDf){
  groupMeans = inputDf %>%
    group_by(day) %>%
    summarise(groupMean = mean(fitness))
  absoluteDelta = abs(groupMeans$groupMean[1]-groupMeans$groupMean[2])
  return(absoluteDelta)
}

calculateAbsoluteFoldChange = function(inputDf){
  groupMeans = inputDf %>%
    group_by(day) %>%
    summarise(groupMean = mean(fitness))
  absoluteFoldChange = abs(groupMeans$groupMean[1]/groupMeans$groupMean[2])
  log2AbsoluteFoldChange = log2(absoluteFoldChange)
  return(log2AbsoluteFoldChange)
}

fitnessScores = read_tsv('barseqAdjustedParams/fit_logratios.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$day = factor(metadata$day, levels = c('t0', 'day1', 'day3', 'day7', 'day14'))
annotations = read_tsv('genesWithAnvioAnnotations.tsv')

metadata$phase = NA
metadata$phase[metadata$dayNumeric < 7] ='early'
metadata$phase[metadata$dayNumeric >= 7] ='late'

metadata$cage = as.factor(metadata$cage)

keggs = annotations %>%
  select(locusId, kofamAccession, kofamFunction)%>%
  distinct()
colnames(fitnessScores) <- sub("setA", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("_.*", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("CO$", "Co", colnames(fitnessScores))
colnames(fitnessScores) <- sub("DJ$", "Dj", colnames(fitnessScores))

djMeta = metadata%>%
  filter(tissue == 'dj')

fitnessScores%>%
  select(-c(sysName,
            desc))%>%
  pivot_longer(cols = 2:ncol(.),
               values_to = 'fitnessScores',
               names_to = 'sample')%>%
  filter(sample %in% djMeta$sample)%>%
  group_by(locusId)%>%
  summarise(locusVariance = var(fitnessScores))%>%
  ggplot(aes(x= log10(locusVariance)))+
  geom_histogram(bins = 100)+
  labs(title = 'Dist of dj variance')

passVarainceFilter=fitnessScores%>%
  filter(desc != 'hypothetical protein',
         grepl("^tRNA-[A-Za-z]{3}$", desc)!=T)%>%
  select(-c(sysName,
            desc))%>%
  pivot_longer(cols = 2:ncol(.),
               values_to = 'fitnessScores',
               names_to = 'sample')%>%
  filter(sample %in% djMeta$sample)%>%
  group_by(locusId)%>%
  summarise(locusVariance = var(fitnessScores))%>%
  filter(locusVariance > 0.5)


loci = passVarainceFilter$locusId

pcaOut=fitnessScores%>%
  filter(locusId %in% loci)%>%
  select(-c(sysName,
            desc))%>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  filter(sample %in% djMeta$sample)%>%
  column_to_rownames('sample')%>%
  prcomp(center = T, scale = T)

pcaSummary=summary(pcaOut)



pcaOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(djMeta,
            by = 'sample')%>%
  ggplot(aes(x = PC1,
             col = day,
             y = PC2))+
  geom_point()+
  labs(x = paste('PC1 - ', pcaSummary$importance[2,1],'%'),
       y = paste('PC2 - ', pcaSummary$importance[2,2],'%'),
       title = 'DJ fitness pca'
  )+
  theme_bw()

lmIn = fitnessScores%>%
  filter(locusId %in% loci)%>%
  select(-c(sysName,
            desc))%>%
  column_to_rownames('locusId')

modelStatsNumericDay <- data.frame(
  gene = character(0),
  pvalue = numeric(0),
  coef = numeric(0),
  absoluteFoldChange = numeric(0),
  absoluteDelta = numeric()
)

for (l in seq_along(loci)) {
  input <- lmIn[rownames(lmIn) == loci[l], ] %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "fitness"
    ) %>%
    merge(djMeta, by = "sample") %>%
    filter(tissue != "T0")

  input$tissue <- as.factor(input$tissue)
  input$millionBases <- scale(input$millionBases)

  lm_fit <- lm(
    fitness ~ dayNumeric + percentPerfectBarcode + millionBases,
    data = input
  )

  coefs <- summary(lm_fit)$coefficients

  pvalue <- coefs["dayNumeric", "Pr(>|t|)"]
  coef <- coefs["dayNumeric", "Estimate"]
  absoluteFoldChange = input%>%
    filter(dayNumeric %in% c(1, 14))%>%
    calculateAbsoluteFoldChange()

  absoluteDelta = input%>%
    filter(dayNumeric %in% c(1, 14))%>%
    calculateMeanDeltaFitness()

  modelStatsNumericDay <- rbind(
    modelStatsNumericDay,
    data.frame(
      gene = loci[l],
      pvalue = pvalue,
      coef = coef,
      absoluteFoldChange = absoluteFoldChange,
      absoluteDelta = absoluteDelta
    )
  )
}

modelStatsNumericDay$padj = p.adjust(modelStatsNumericDay$pvalue, method = 'fdr')

sigNumericDay=modelStatsNumericDay%>%
  filter(padj < .1,
         absoluteDelta> 1)

sigNumericDay%>%
  ggplot(aes(x = absoluteDelta))+
  geom_histogram(bins=100)

'Model is working well'
modelStatsNumericDay%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 100)+
  labs(title = 'Dist of raw p-values numeric day DJ linear model')

sigNumericDay = sigNumericDay %>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in unique(sigNumericDay$locusId)){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(djMeta, by = 'sample')%>%
    filter(tissue != 'T0',
           locusId == s)%>%
    ggplot(aes(x = day, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sigNumericDay$coef[sigNumericDay$locusId==s],
                          '\np-adj :', sigNumericDay$padj[sigNumericDay$locusId==s]))
  plot(p)
}

write_tsv(modelStatsNumericDay, 'linear models/linearRegresssionDjOnly/dayNumericaResDj.tsv')

day1Day3Samples=djMeta%>%
  filter(day %in% c('day1', 'day3'))%>%
  .$sample

day1Day3Fitness=lmIn[colnames(lmIn)%in% day1Day3Samples]

day1Day3LmRes=run_lm(day1Day3Fitness)

day1Day3LmRes$padj = p.adjust(day1Day3LmRes$pvalue, method = 'fdr')

day1Day3LmRes%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 50)+
  labs(title = 'Linear model')

day1Day3LmRes$comparison = 'day1Vday3'

day1Day3LmRes$significant = 'No'
day1Day3LmRes$significant[day1Day3LmRes$padj < .1] = 'Yes'

day1Day3LmRes%>%
  ggplot(aes(x = absoluteDelta,
             y = -log10(padj),
             col = significant,
  ))+
  geom_point()+
  labs(title = 'Linear model results')

day1Day3LmRes%>%
  ggplot(aes(x = coef,
             y = -log10(padj),
             col = significant,
  ))+
  geom_point()+
  labs(title = 'Linear model results')

day1Day3LmRes%>%
  filter(padj < .1,
         absoluteDelta> 1)

day1Day3LmRes%>%
  filter(padj < .1)

sigDay1Day3= day1Day3LmRes %>%
  filter(padj < .1,
         absoluteDelta> 1)%>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

sigDay1Day3$comparison = 'day1Vday3'
sigDay1Day3  %>%
  write_tsv('linear models/linearRegresssionDjOnly/day1Vday3Res.tsv')

for (s in sigDay1Day3$locusId){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(djMeta, by = 'sample')%>%
    filter(tissue != 'T0',
           day %in% c('day1', 'day3'),
           locusId == s)%>%
    ggplot(aes(x = day, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    ylim(-10,10)+
    labs(title = s,
         caption = paste0('Coef: ', sigDay1Day3$coef[sigDay1Day3$locusId==s],
                          '\np-adj :', sigDay1Day3$padj[sigDay1Day3$locusId==s]))
  plot(p)
}

day3Day7Samples=djMeta%>%
  filter(day %in% c('day3', 'day7'))%>%
  .$sample

day3Day7Fitness=lmIn[colnames(lmIn)%in% day3Day7Samples]

day3Day7LmRes=run_lm(day3Day7Fitness)

day3Day7LmRes$padj = p.adjust(day3Day7LmRes$pvalue, method = 'fdr')

day3Day7LmRes%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 50)

day3Day7LmRes%>%
  ggplot(aes(x = absoluteDelta,
             y = -log10(padj)))+
  geom_point()

day3Day7LmRes%>%
  filter(padj < .1,
         absoluteDelta> 1)

day3Day7LmRes$comparison = 'day3Vday7'
day3Day7LmRes  %>%
  write_tsv('linear models/linearRegresssionDjOnly/day3Vday7Res.tsv')

sigDay3Day7= day3Day7LmRes %>%
  filter(padj < .1,
         absoluteDelta> 1)%>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in sigDay3Day7$locusId[1:10]){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(djMeta, by = 'sample')%>%
    filter(tissue != 'T0',
           day %in% c('day3', 'day7'),
           locusId == s)%>%
    ggplot(aes(x = day, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    ylim(-10, 10)+
    labs(title = s,
         caption = paste0('Coef: ', sigDay3Day7$coef[sigDay3Day7$locusId==s],
                          '\np-adj :', sigDay3Day7$padj[sigDay3Day7$locusId==s]))
  plot(p)
}


day7Day14Samples=djMeta%>%
  filter(day %in% c('day7', 'day14'))%>%
  .$sample

day7Day14Fitness=lmIn[colnames(lmIn)%in% day7Day14Samples]

day7Day14LmRes=run_lm(day7Day14Fitness)

day7Day14LmRes$padj = p.adjust(day7Day14LmRes$pvalue, method = 'fdr')

day7Day14LmRes$comparison = 'day7Vday14'
day7Day14LmRes  %>%
  write_tsv('linear models/linearRegresssionDjOnly/day7Vday14Res.tsv')


day7Day14LmRes%>%
  ggplot(aes(x = pvalue))+
  geom_histogram(bins = 50)

day7Day14LmRes%>%
  filter(padj < .1,
         absoluteDelta> 1)

sigDay7Day14= day7Day14LmRes %>%
  filter(padj < .1,
         absoluteDelta> 1)%>%
  rename(locusId = gene)%>%
  left_join(keggs, by = 'locusId')%>%
  arrange(desc(-log10(padj)))

for (s in sigDay7Day14$locusId){
  p=fitnessScores%>%
    pivot_longer(cols = c(4:ncol(.)), names_to = 'sample', values_to = 'fitnessScore')%>%
    left_join(djMeta, by = 'sample')%>%
    filter(tissue != 'T0',
           day %in% c('day7', 'day14'),
           locusId == s)%>%
    ggplot(aes(x = day, y = fitnessScore))+
    geom_boxplot(outliers = F)+
    geom_point(aes(col = day))+
    labs(title = s,
         caption = paste0('Coef: ', sigDay7Day14$coef[sigDay7Day14$locusId==s],
                          '\np-adj :', sigDay7Day14$padj[sigDay7Day14$locusId==s]))
  plot(p)
}

sigCategoricalDay=sigDay1Day3%>%
  select(-significant)%>%
  rbind(sigDay3Day7)%>%
  rbind(sigDay7Day14)

colnames(sigCategoricalDay)
colnames(sigNumericDay)

sigNumericDay=sigNumericDay %>%
  select(-absoluteFoldChange)

sigCategoricalDay %>%
  select(-comparison)%>%
  rbind(sigNumericDay)%>%
  write_tsv('linear models/linearRegresssionDjOnly/sigRes.tsv')
