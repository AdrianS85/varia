print: 
	#The part before ":" part can be either name of the file to recruit? or just like a name of the function to run. Also, these are called "rules" for some reason. Tutaj też można dodać prerequesites (files that need to exist) for running this function. Jeśli zrobimy z tego file, to zostanie on stworzony jeśli nie istnieje z początku	
	@echo "ass"

#for analyzing all the files in given folder; make fastqcMake; put all needed scripts into "Programs" folder; questions: czy napewno dobre kodowanie?,
#Struktura folderów: Programs, Genomes, #foldery_z_rawseqfiles#

PrepareMake: 
#PREPARE FOLDERS
	export CORES=32 #This is global variable
	mkdir Raw_Data Fastqc_Raw Trim_Galore Diversity_Cut Fastq_Screen_Div_Cut Bismark Bismark_Strip_and_Div Bismark_Extracted Bedgraph 
	mkdir ./Trim_Galore/Trim_Galore_Raports ./Diversity_Cut/Fastqc_Trimmed ./Bismark/Bismark_Report ./Bismark/Bismark_Summary ./Bismark/Bismark_Raw
	cp ../Programs/FastQC_aggregate.sh ./Fastqc_Raw/; #cd Fastqc_Raw; FastQC_aggregate.sh; cd ..;
	cp ../Programs/FastQC_aggregate.sh ./Diversity_Cut/Fastqc_Trimmed;
	cp ../Programs/trimRRBSdiversityAdaptCustomers.py ./Trim_Galore/;
	cp ../Programs/strip_bismark_sam.sh ./Bismark
#PREPARE FOLDERS

# Programs needed: FastQC, FastQ Screen

GetGenomesMake:
#GET GENOMES
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.fna.gz
	mv *fna* ../Genomes
	gzip -d ../Genomes/*.gz
	rename 's/\.fna$/.fa/' * #Change names so that they are compatible with Bismark
#GET GENOMES

BismarkGenomeMake:
#MAKE CONVERTED GENOMES
	#cp -R ../Programs/*ismark-* ../Genomes
	#
	../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ../Genomes/Rat_6
	../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ../Genomes/Mouse*
	../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ../Genomes/Human_38_12
	../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ../Genomes/Ecoli
#MAKE CONVERTED GENOMES

TrimMake:
#FASTQC FROM RAW DATA
	fastqc --threads $(CORES) --outdir ./Fastqc_Raw *fastq*
	# Aggregate the results https://gist.github.com/danielecook/8e9afb2d2df7752efd8a#file-fastqc_aggregate-sh
#FASTQC FROM RAW DATA
#TRIM_GALORE
	ls *R1* | sort >> r1; ls *R3* | sort >> r2; paste r1 r2 >> read_pairs; rm r1 r2 #Outputs list of paired reads files
	parallel --bar --colsep '\t' trim_galore --paired  --retain_unpaired --output_dir ./Trim_Galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC {1} {2} >> trim_galore_report  :::: read_pairs #Using paired, columned list we pair the names
	mv ./Trim_Galore/*trimming_report.txt ./Trim_Galore/Trim_Galore_Raports; rm read_pairs
#TRIM_GALORE
#DIVERSITY CUTTING
	cd Trim_Galore
	ls *R1*val* | sort >> r1; ls *R3*val* | sort >> r2; paste r1 r2 >> read_pairs; rm r1 r2
	parallel --bar --colsep '\t' python trimRRBSdiversityAdaptCustomers.py -1 {1} -2 {2} >> div_cut_report :::: read_pairs
	mv *trimmed* ../Diversity_Cut
	rm read_pairs; cd ../Diversity_Cut
	fastqc --threads $CORES --outdir ./Fastqc_Trimmed *trimmed* #On server I dont need the paranthesis in CORES but on my desktop I do? wtf?
	cd ..
	mv *fastq.gz ./Raw_Data/ #Moving raw data to appropriate folder
#DIVERSITY CUTTING

#DOESNT WORK ON SERWER THOUGH IT WORKS ON DESKTOP
FastqcScreenMake:
#FASTQ SCREEN SETUP
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
#DOESNT WORK ON SERWER THOUGH IT WORKS ON DESKTOP


BismarkMake:
#MAKE BAM
	cd Diversity_Cut
	ls *val_1* | sort >> r1; ls *val_2* | sort >> r2; paste r1 r2 >> read_pairs; rm r1 r2
	parallel --bar --colsep '\t' ../../Programs/*ismark-*/bismark --bowtie2  --genome_folder ../../Genomes/Mouse*/ -1 {1} -2 {2} >> bismark_raport :::: read_pairs
	#Library is assumed to be strand-specific (directional), alignments to strands complementary to the original top or bottom strands will be ignored (i.e. not performed!)
	#?rm *G_to_A* *C_to_T*
	mv ./*.bam ../Bismark; rm read_pairs
	mv ./*bt2_PE_report.txt ../Bismark/Bismark_Report
#MAKE BAM
#STRIP OVATION-SPECIFIC
	cd ../Bismark
	ls *.bam >> r1; cp r1 r2; sed 's/\.bam$/.sam/' r2 >> r3; paste r1 r3 >> read_pairs; rm r1 r2 r3 # Getting nice names for sam files
	parallel samtools view -h -o {2} {1} :::: read_pairs
	#ls *bam | parallel strip_bismark_sam.sh
	#mv ./*.bam ./Bismark/Bismark_Raw; rm read_pairs #We dont need raw bismark files anymore
#STRIP OVATION-SPECIFIC
#DEDUPULICATION OVATION-SPECIFIC
	#ls *sam_stripped | python nudup.py –f index.fq
#DEDUPULICATION OVATION-SPECIFIC
#BAMQC
	#ls *pe.bam | parallel bamqc
#BAMQC
#GENERATE BISMARK REPORTS
	#cd ../Bismark
	#ls | parallel --bar --colsep '\t' ../../Programs/*ismark-*/bismark2summary *bam
	#mv bismark_summary_report* ./Bismark_Summary
	#ls *bam | parallel --bar --colsep '\t' ../../Programs/*ismark-*/bismark2report --dir ./Bismark_Report *bam #This is final report after everything
#GENERATE BISMARK REPORTS
	#ls | bam2nuc --genome_folder ../../Genome/Mouse*/ --dir ../../Programs/*ismark-*/bismark2report
#METHYLATION CALLING
	#ls *sam >> files
	#parallel --colsep '\t' ../../Programs/*ismark-*/bismark_methylation_extractor --output ../Bedgraph --paired-end --comprehensive --merge_non_CpG >> bismark2bedGraph_report :::: files
#METHYLATION CALLING

#RNBEADS SCRIPT
	#FOLDER STRUCTURE
(analysis.dir) -> report.dir; (data.dir) -> sample.annotation, dataset.dir
data.dir <- "~/RnBeads/data/Ziller2011_PLoSGen_450K"
dataset.dir <- file.path(data.dir, "idat")
sample.annotation <- file.path(data.dir, "sample_annotation.csv")
analysis.dir <- "~/RnBeads/analysis" # Directory where the output should be written to
report.dir <- file.path(analysis.dir, "reports") # Directory where the report files should be written to
	#FOLDER STRUCTURE
parallel.setup(num.cores)
parallel.isEnabled()
"hg19" %in% rnb.get.assemblies()
	RnBeads.mm10
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


