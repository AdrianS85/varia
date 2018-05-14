print: 
	#The part before ":" part can be either name of the file to recruit? or just like a name of the function to run. Also, these are called "rules" for some reason. Tutaj też można dodać prerequesites (files that need to exist) for running this function. Jeśli zrobimy z tego file, to zostanie on stworzony jeśli nie istnieje z początku	
	@echo "ass"

#for analyzing all the files in given folder; make fastqcMake; put all needed scripts into "Programs" folder; questions: czy napewno dobre kodowanie?,
#Struktura folderów: Programs, Genomes, #foldery_z_rawseqfiles#
#!! ADD nohup TO PREVENT SHUTTING DOWN PROCESSES AFTER EXITING TERMINAL
PrepareMake: 
#PREPARE FOLDERS
	export CORES=32 #This is global variable
	mkdir Raw_Data Fastqc_Raw Trim_Galore Diversity_Cut Fastq_Screen_Div_Cut Bismark Bismark_Strip_and_Dedup Bedgraph && 
	mkdir ./Trim_Galore/Trim_Galore_Raports ./Diversity_Cut/Fastqc_Trimmed ./Bismark/Bismark_Report ./Bismark/Bismark_Raw && 
	mkdir ./Bismark/Bismark_Raw/Bismark_Raw_Bamqc ./Bismark_Strip_and_Dedup/Bismark_Bamqc ./Bismark_Strip_and_Dedup/Strip_and_Dedup_Report && 
	cp ../Programs/FastQC_aggregate.sh ./Fastqc_Raw/ && 
	cp ../Programs/FastQC_aggregate.sh ./Diversity_Cut/Fastqc_Trimmed && 
	cp ../Programs/trimRRBSdiversityAdaptCustomers.py ./Trim_Galore/ && 
	cp ../Programs/strip_bismark_sam.sh ./Bismark && 
	cp ../Programs/nugentechnologies*/nudup.py ./Bismark && 
	mv *_R2_* ./Bismark && 
	parallel gzip -d ::: ./Bismark/*.gz && 
	cd Bismark && 
		rename 's/\.fastq$/.fq/' *.fastq && 
	cd ..
#PREPARE FOLDERS

# Programs needed: FastQC, FastQ Screen

GetGenomesMake:
#GET GENOMES
	#A contig sequence that is released outside of the full assembly release cycle. These 
	#sequences are meant to add information to the assembly without disrupting the stable coordinate system. There are 
	#two types of patches, FIX and NOVEL. FIX patches are released to correct an error 
	#in the assembly and will be removed when the new full assembly is released. 
	#NOVEL sequences are sequences that were not in the last full assembly release 
	#and will be retained with the next full assembly release.
	#http://hgdownload.soe.ucsc.edu/downloads.html
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.fna.gz
	parallel gzip -d ::: ./*.fna.gz ### From here to checkup
	rename 's/\.fna$/.fa/' ./*.fna #Change names so that they are compatible with Bismark
	mkdir ../Genomes/Human_38_12 ../Genomes/Ecoli ../Genomes/Rat_6 ../Genomes/Mouse_38_6 
	mv *GRCh38* ../Genomes/Human_38_12; mv *GRCm38* ../Genomes/Mouse_38_6; mv *ASM584v2* ../Genomes/Ecoli mv *Rnor* ../Genomes/Rat_6;
#GET GENOMES

BismarkGenomeMake:
#MAKE CONVERTED GENOMES
	#cp -R ../Programs/*ismark-* ../Genomes
	#
	parallel ../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ::: ../Genomes/Rat_6
	parallel ../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ::: ../Genomes/Mouse*
	parallel ../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ::: ../Genomes/Human_38_12
	parallel ../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ::: ../Genomes/Ecoli
#MAKE CONVERTED GENOMES

TrimMake:
#FASTQC FROM RAW DATA
	export CORES=32
	fastqc --threads ${CORES} --outdir ./Fastqc_Raw *fastq*
	cd ./Fastqc_Raw 
	chmod 755 FastQC_aggregate.sh
	./multiqc .
	cd ../
#FASTQC FROM RAW DATA
#TRIM_GALORE
	ls *R1* | sort >> r1; ls *R3* | sort >> r2; paste r1 r2 >> read_pairs; rm r1 r2 #Outputs list of paired reads files
	parallel --colsep '\t' "trim_galore --paired  --retain_unpaired --output_dir ./Trim_Galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC {1} {2}" :::: read_pairs #Using paired, columned list we pair the names
	mv ./Trim_Galore/*trimming_report.txt ./Trim_Galore/Trim_Galore_Raports && rm read_pairs 
	#ADD multiqc HERE
#TRIM_GALORE
#DIVERSITY CUTTING
	cd Trim_Galore
	ls *R1*val* | sort >> r1; ls *R3*val* | sort >> r2; paste r1 r2 >> read_pairs; rm r1 r2
	parallel --colsep '\t' "python ./trimRRBSdiversityAdaptCustomers.py -1 {1} -2 {2}" :::: read_pairs
	mv *trimmed* ../Diversity_Cut
	rm read_pairs; cd ../Diversity_Cut
	fastqc --threads $CORES --outdir ./Fastqc_Trimmed *trimmed*
	#ADD multiqc HERE
	cd ..
	mv *fastq.gz ./Raw_Data/ #Moving raw data to appropriate folder
#DIVERSITY CUTTING

#LATER##DOESNT WORK ON SERWER THOUGH IT WORKS ON DESKTOP 
FastqcScreenMake:
#FASTQ SCREEN SETUP #Also checks if we have the best reference genome?
	cp ../../Programs/fastq_screen*/fastq_screen.conf.example ../../Programs/fastq_screen*/fastq_screen.conf
	vim ../../Programs/fastq_screen*/fastq_screen.conf 
	#THREADS 	$CORES
	#BOWTIE2	/usr/bin/bowtie2 	
	#BISMARK 	~/data/Programs/*ismark-*/bismark
	#DATABASE	~/data/Genomes/Mouse*/Bisulfite_Genome
	
#FASTQ SCREEN SETUP
#FASTQ SCREEN 
#needs bowtie2 and bismark; needs to configure fastq_screen.config; download contaminating genomes (downloaded from ncbi Assembly?); download contaminating adapters (see config file to see from where); for bismark option you have to set config databases to "Bisulfite_Genome" folders created by bismark_genome_preparation tool;
	cd Diversity_Cut # cd ~/data/FE/Diversity_Cut
	ls *trimmed* >> reads
	parallel ../../Programs/fastq_screen*/fastq_screen --aligner bowtie2 --bisulfite -outdir ../Fastq_Screen_Div_Cut/ :::: reads
	cd ..
#FASTQ SCREEN
#LATER##DOESNT WORK ON SERWER THOUGH IT WORKS ON DESKTOP


BismarkMake:
#MAKE BAM
	cd Diversity_Cut
	ls *val_1* | sort >> r1; ls *val_2* | sort >> r2; paste r1 r2 >> read_pairs; rm r1 r2
	parallel --colsep '\t' ../../Programs/*ismark-*/bismark --bowtie2  --genome_folder ../../Genomes/Mouse*/ -1 {1} -2 {2} :::: read_pairs
	#Library is assumed to be strand-specific (directional), alignments to strands complementary to the original top or bottom strands will be ignored (i.e. not performed!)
	mv ./*.bam ../Bismark; 
	mv ./*bt2_PE_report.txt ../Bismark;
	rm read_pairs
#MAKE BAM
#LATER##GENERATE BISMARK REPORTS - ADD MULTIQC STEP!
	cd ../Bismark/
	#mv bismark_summary_report* ./Bismark_Summary
	parallel ../../Programs/*ismark-*/bam2nuc --genome_folder ../../Genomes/Mouse_mm10/ ::: *_pe.bam
	parallel ../../Programs/*ismark-*/bismark2summary ::: *_pe.bam
	parallel --colsep '\t' ../../Programs/*ismark-*/bismark2report --dir ./Bismark_Report ::: *bam #This is final report after everything
	mv *PE_report.* ./Bismark_Report
	mv *nucleotide_stats* ./Bismark_Report
	#ADD multiqc HERE
#LATER##GENERATE BISMARK REPORTS
#STRIP OVATION-SPECIFIC
	ls *.bam >> r1; cp r1 r2; sed 's/\.bam$/.sam/' r2 >> r3; paste r3 r1 >> read_pairs; rm r1 r2 r3 # Getting nice names for sam files
	parallel --colsep '\t' "samtools view -h -o {1} {2}"  :::: read_pairs
	#mv ./*.bam ./Bismark/Bismark_Raw; rm read_pairs #We dont need raw bismark files anymore
	parallel ./strip_bismark_sam.sh {} ::: *sam
	rm *pe.sam read_pairs
#STRIP OVATION-SPECIFIC
#DEDUPULICATION OVATION-SPECIFIC
	ls *_R2_* >> r1; ls *_stripped.sam >> r3; sed 's/_stripped.sam$/_stripped_dedup.sam/' r3 >> r2; paste r1 r2 r3  >> read_pairs; rm r1 r2 r3
	parallel --colsep '\t' "python ./nudup.py --paired-end -f {1} -o {2} {3}" :::: read_pairs #Can also use .bam
	mv *sam_dup_log.txt ../Bismark_Strip_and_Dedup/Strip_and_Dedup_Report
#DEDUPULICATION OVATION-SPECIFIC
#BAMQC
	parallel bamqc ::: *.pe.bam
	mv *bamqc.html ../Bismark_Raw/Bismark_Raw_Bamqc; mv *bamqc.zip ./Bismark_Raw/Bismark_Raw_Bamqc;
	parallel bamqc ::: *sorted.dedup*
	mv *bamqc.html ../Bismark_Strip_and_Dedup/Bismark_Bamqc; 
	mv *bamqc.zip ../Bismark_Strip_and_Dedup/Bismark_Bamqc;
	
#BAMQC
#METHYLATION CALLING
	ls *sorted.dedup.bam >> r2 && cp r2 r1_1 && sed 's/sorted.dedup.bam$/sorted.dedup.final.bam/' r1_1 >> r1 && paste r1 r2 >> read_pairs && rm r1 r1_1 r2
	parallel --colsep '\t' samtools sort -n -o {1} {2} :::: read_pairs ## Needs to sort the files cause they are desorted by nugen deduplication script
	parallel "../../Programs/*ismark-*/bismark_methylation_extractor --bedGraph --output ../Bedgraph --paired-end --comprehensive --merge_non_CpG" ::: *final.bam
#The IDs of Read 1 (NB500931:147:HC5VCBGX5:3:21604:26239:14679) and Read 2 (NB500931:147:HC5VCBGX5:2:11204:15206:10902) are not the same.                                      
#This might be the result of sorting the paired-end SAM/BAM files by chromosomal position which is not compatible with correct methylation                                     
#extraction. Please use an unsorted file instead or sort the file using 'samtools sort -n' (by read name). This may also occur using                                      
#samtools merge as it does not guarantee the read order. To properly merge files please use 'samtools merge -n' or 'samtools cat'.
#METHYLATION CALLING
#CLEANUP
	mv *sam.sorted.dedup* ../Bismark_Strip_and_Dedup
	mv *sam.sorted.markdup* ../Bismark_Strip_and_Dedup
	mv bismark_summary_report* ./Bismark_Report
	rm ./*sam_stripped.sam
#CLEANUP



#RNBEADS SCRIPT (analysis.dir) -> report.dir; (data.dir) -> sample.annotation, dataset.dir; Patched genomes dont work!
#GET THE TOOLS
source("https://bioconductor.org/biocLite.R")
biocLite("RnBeads.mm9")
biocLite("RnBeads.mm10")
biocLite("RnBeads")
biocLite("doParallel")
biocLite("ggbio")
biocLite("rtracklayer")
biocLite("impute")
library(RnBeads)
#GET THE TOOLS


#SETUP WORKING ENVIROMENT
setwd("D:/Zycie_zawodowe/Fede_seq/rnbeads")
data_dir <- paste0(getwd(), "/data_dir")
dataset_dir <- paste0(data_dir, "/dataset_dir")
sample_annotation <- paste0(data_dir, "/sample_annotation.csv")
analysis_dir <- paste0(getwd(), "/analysis_dir")
report_dir <- paste0(analysis_dir, "/reports")
#SETUP WORKING ENVIROMENT


#SETUP PARALLELISM
parallel.setup(4)
parallel.isEnabled()
#SETUP PARALLELISM


#SETUP RUN OPTIONS
rnb.options(
  assembly = "mm10",
  import.bed.style = "bismarkCov",
  identifiers.column="sampleID",
  qc.coverage.plots = TRUE,
  qc.coverage.histograms = TRUE,
  inference = TRUE,
  #min.group.size = 1,
  differential.comparison.columns = "treatment")
#SETUP RUN OPTIONS


Sys.getenv()
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.23/bin/gswin64c.exe")
Sys.setenv(R_ZIPCMD="C:/Rtools/bin/zip.exe")
Sys.setenv(R_UNZIPCMD="C:/Rtools/bin/unzip.exe")


#try creating .Rprofile and adding it there

#C:\Rtools\bin
#RUN
rnb.run.analysis(
  dir.reports=report_dir,
  sample.sheet=sample_annotation,
  data.dir=dataset_dir,
  data.type="bs.bed.dir"
)

#STEP-BY-STEP DEBUGGING
result <- rnb.run.import(data.source=c(dataset_dir, sample_annotation),
                         data.type="bs.bed.dir", 
                         dir.reports=report_dir)


ADD!
rnb.execute.low.coverage.masking
rnb.execute.na.removal
rnb.execute.high.coverage.removal
rnb.execute.snp.removal
rnb.execute.variability.removal

DOWNLOADING ANNOTATIONS:
https://genome.ucsc.edu/cgi-bin/hgTables


methyAnalysis - doesnt have specific functions for seq?
M3D - perhaps
methylPipe - perhaps






setwd("D:/Zycie_zawodowe/Fede_seq/rnbeads")

data.dir <- paste0(getwd(), "/data.dir")
dataset.dir <- file.path(data.dir, "idat")
sample.annotation <- file.path(data.dir, "sample_annotation.csv")
analysis.dir <- "~/RnBeads/analysis" # Directory where the output should be written to
report.dir <- file.path(analysis.dir, "reports") # Directory where the report files should be written to
	#FOLDER STRUCTURE
parallel.setup(num.cores)
parallel.isEnabled()

rnb.options(
	assembly = "mm10", #Is this the same genome which was used?
	import.bed.style = "bismarkCov",
	identifiers.column="Sample_ID",
	qc.coverage.plots = TRUE,
	qc.coverage.histograms = TRUE,
	inference = TRUE)
rnb.run.analysis(
	data.source = c(dataset.dir, sample.annotation, "the index of the column of the sample annotation sheet that contains the names or full
paths to the bed files")
	dir.reports = report.dir, 
	sample.sheet = sample.annotation,
	data.dir=dataset.dir, 
	data.type="bs.bed.dir"),
	rnb.run.qc(rnb.set, report.dir)
#RNBEADS SCRIPT












### RNBEADS EXAMPLE
################################################################################
# DNA methylation analysis with RnBeads
# Epigenomics 2016 Workshop
# ------------------------------------------------------------------------------
#  Vanilla analysis of the Ziller2011 450K dataset
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# created: 2016-01-22
# author:  Fabian Mueller <rnbeads@mpi-inf.mpg.de>
# http://rnbeads.mpi-inf.mpg.de/
################################################################################

################################################################################
# (0) Preliminaries
################################################################################
# load the package
library(RnBeads)

# define the directory structure
# setwd(".")
dataDir <- file.path(getwd(), "data")
resultDir <- file.path(getwd(), "results")

# dataset and file locations
datasetDir <- file.path(dataDir, "Bock2012_MolCell_RRBS")
bedDir <- file.path(datasetDir, "dataset", "bed")
sampleSheet <- file.path(datasetDir, "dataset", "sample_annotation.csv")
reportDir <- file.path(resultDir, "report_Bock2012_vanilla")
################################################################################
# (1) Set analysis options
################################################################################
rnb.options(
	identifiers.column                = "sampleID",
	import.bed.style                  = "EPP",
	assembly                          = "mm9",
	filtering.low.coverage.masking    = TRUE,
	filtering.greedycut               = FALSE,
	filtering.missing.value.quantile  = 0.5,
	filtering.high.coverage.outliers  = TRUE
)
# optionally disable some parts of the analysis to reduce runtime
rnb.options(
	exploratory.intersample           = FALSE,
	exploratory.region.profiles       = character(0),
	differential.comparison.columns   = c("differentiation_level",
		"blood_lineage", "cmp_blood_CLP_CMP"),
	differential.report.sites         = FALSE
)
################################################################################
# (2) Run the analysis
################################################################################
rnb.run.analysis(
	dir.reports=reportDir,
	sample.sheet=sampleSheet,
	data.dir=bedDir,
	data.type="bs.bed.dir"
)

################################################################################
# Link to finished analysis
################################################################################
# see the results at:
# http://rnbeads.mpi-inf.mpg.de/reports/tutorial/epigenomics2016/results/report_Bock2012_vanilla/index.html



