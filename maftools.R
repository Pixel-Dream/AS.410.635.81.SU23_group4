setwd('C:/Users/ashwi/ToolsforGenomeAnalysis')
library(maftools)
#path to TCGA LAML MAF file
my_maf = "C:/Users/ashwi/ToolsforGenomeAnalysis/merged_100_filtered.maf"
laml = read.maf(maf = my_maf)
getGeneSummary(laml)
oncoplot(maf = laml, top = 10)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
lollipopPlot(
  maf = laml,
  gene = 'BRCA2',
  showMutationRate = TRUE,
  
)