#GET THE TOOLS
library(RnBeads)
#GET THE TOOLS


#SETUP WORKING ENVIROMENT
#setwd("D:/Zycie_zawodowe/Fede_seq/rnbeads")
setwd(getwd())
data_dir <- paste0(getwd(), "/data")
dataset_dir <- data_dir
sample_annotation <- paste0(data_dir, "/annot.csv")
analysis_dir <- paste0(getwd(), "/analysis_dir")
report_dir <- paste0(analysis_dir, "/reports")
#SETUP WORKING ENVIROMENT


#SETUP PARALLELISM
parallel.setup(4)
parallel.isEnabled()
#SETUP PARALLELISM


#SETUP RUN OPTIONS
rnb.options(
  analysis.name = "Fede_Brain1",
  #import.bed.test.only = T,
  assembly = "mm10",
  region.aggregation = "coverage.weighted",
  identifiers.column="sampleID", ###Can behave wierdly. doesnt line "_"? integers?
  #IMPORT
  #import.default.data.type = "bs.bed.dir",
  import.bed.style = "bismarkCov",
  import.gender.prediction = F,
  #import.bed.frame.shift
  #QC
  qc.coverage.plots = TRUE,
  qc.coverage.histograms = TRUE,
  qc.coverage.violins = FALSE,
  qc.coverage.threshold.plot = c(1, 2, 3, 5, 10),
  #FILTERING
  ###filtering.snp = "no",
  ###filtering.missing.value.quantile = 0.4,
  #filtering.coverage.threshold = 5,
  ###filtering.low.coverage.masking = T,
  #filtering.high.coverage.outliers = F, 
  #filtering.deviation.threshold = 0,
  #INFERENCE
  inference = TRUE,
  inference.targets.sva = c("prep_batch", "seq_batch"),
  #inference.reference.methylome.column = 
  #inference.max.cell.type.markers = 
  #inference.top.cell.type.markers = 
  #EXPLORATORY
  #exploratory.columns
  #exploratory.gene.symbols  = 
  #exploratory.deviation.plots = 
  #exploratory.custom.loci.bed = 
  exploratory.clustering.top.sites = 200,
  #DIFFERENTIAL
  #min.group.size = 1,
  #differential.comparison.columns.all.pairwise = 
  #differential.variability = 
  #differential.variability.method = 
  #differential.permutations = 
  differential.enrichment.go = F,
  differential.adjustment.sva = T,
  covariate.adjustment.columns = c("prep_batch", "seq_batch"),
  differential.comparison.columns = "treatment"
  #OTHER
  #export.to.ewasher
  )
#SETUP RUN OPTIONS


#STEP-BY-STEP ANALYSIS
## DATA_IMPORT -> (QUALITY_CONTROL, PREPROCESSING), PREPROCESSING -> INFERENCE -> (EXPLORATION, DIFFERENTIAL)
result <- rnb.run.import(data.source=c(dataset_dir, sample_annotation),
                         data.type="bs.bed.dir", 
                         dir.reports=report_dir) ## results in a list with two elements: the dataset (rnb.set) and a report

rnb.run.qc(result$rnb.set, report_dir) ## Outputs just report i think. I think it needs to write in /dev
PP_result <- rnb.run.preprocessing(result$rnb.set, report_dir)
I_result <- rnb.run.inference(PP_result$rnb.set, report_dir)
E_result <- rnb.run.exploratory(I_result$rnb.set, report_dir)
D_result <- rnb.run.differential(I_result$rnb.set, report_dir)
rnb.run.tnt(RNBset, report_dir)




2018-10-08 12:10:25    12.9   ERROR                                 Could not create file. cannot popen '/usr/bin/which 'gs' 2>/dev/null', probable reason 'Cannot allocate memory'
Error in logger.error(em) : 
  Could not create file. cannot popen '/usr/bin/which 'gs' 2>/dev/null', probable reason 'Cannot allocate memory'
