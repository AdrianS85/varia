#GET THE TOOLS
source("https://bioconductor.org/biocLite.R")
biocLite("RnBeads.mm9")
biocLite("RnBeads.mm10")
biocLite("RnBeads")
biocLite("doParallel")
biocLite("ggbio")
biocLite("rtracklayer")
biocLite("impute")
biocLite("sva")
biocLite("GOstats")
biocLite("org.Mm.eg.db")
library(RnBeads)
library(GOstats)
#GET THE TOOLS


#SETUP WORKING ENVIROMENT
setwd("D:/Zycie_zawodowe/Fede_seq/rnbeads")
data_dir <- paste0(getwd(), "/data_dir")
dataset_dir <- paste0(data_dir, "/dataset_dir")
sample_annotation <- paste0(data_dir, "/sample_annotationL.csv")
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
  identifiers.column="sampleID",
  
  #import.default.data.type = "bs.bed.dir",
  import.bed.style = "bismarkCov",
  import.gender.prediction = F,
  #import.bed.frame.shift
  
  qc.coverage.plots = TRUE,
  qc.coverage.histograms = TRUE,
  qc.coverage.violins = FALSE,
  qc.coverage.threshold.plot = c(1, 2, 3, 5, 10),
  
  filtering.snp = "no",
  filtering.missing.value.quantile = 0.4,
  #filtering.coverage.threshold = 5,
  filtering.low.coverage.masking = T,
  #filtering.high.coverage.outliers = F, 
  #filtering.deviation.threshold = 0,
  
  inference = TRUE,
  inference.targets.sva = c("prep_batch", "seq_batch"),
  #inference.reference.methylome.column = 
  #inference.max.cell.type.markers = 
  #inference.top.cell.type.markers = 
  
  #exploratory.columns
  #exploratory.gene.symbols  = 
  #exploratory.deviation.plots = 
  #exploratory.custom.loci.bed = 
  exploratory.clustering.top.sites = 200,
  
  #min.group.size = 1,
  #differential.comparison.columns.all.pairwise = 
  #differential.variability = 
  #differential.variability.method = 
  #differential.permutations = 
  differential.enrichment.go = F,
  differential.adjustment.sva = T,
  covariate.adjustment.columns = c("prep_batch", "seq_batch"),
  differential.comparison.columns = "treatment"
  
  #export.to.ewasher
  )
#SETUP RUN OPTIONS


#SETUP RUN ENVIROMENT
Sys.getenv()
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.23/bin/gswin64c.exe")
Sys.setenv(R_ZIPCMD="C:/Rtools/bin/zip.exe")
Sys.setenv(R_UNZIPCMD="C:/Rtools/bin/unzip.exe")
#SETUP RUN ENVIROMENT


#RUN
rnb.run.analysis(
  data.type="bs.bed.dir",
  dir.reports=report_dir,
  sample.sheet=sample_annotation,
  data.dir=dataset_dir
)
#RUN


#STEP-BY-STEP DEBUGGING
result <- rnb.run.import(data.source=c(dataset_dir, sample_annotation),
                         data.type="bs.bed.dir", 
                         dir.reports=report_dir) ## results in a list with two elements: the dataset (rnb.set) and a report
RNBset <- result$rnb.set
QC_RNBset <- rnb.run.qc(RNBset, report_dir) ## Outputs just report i think
PP_RNBset <- rnb.run.preprocessing(QC_RNBset, report_dir)
I_RNBset <- rnb.run.inference(PP_RNBset$rnb.set, report_dir)
E_RNBset <- rnb.run.exploratory(I_RNBset$rnb.set, report_dir)
rnb.run.differential(RNBset, report_dir)
rnb.run.tnt(RNBset, report_dir)

Sys.getenv()
Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.23/bin/gswin64.exe")
#https://groups.google.com/forum/#!topic/epigenomicsforum/TY7q2togFX4
rnb.run.qc(rnb.set, report.dir)
#RUN
debugger()
rnb.run.example(index = 2L, dir.output = "example2")
#STEP-BY-STEP DEBUGGING
