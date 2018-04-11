print: 
	#The part before ":" part can be either name of the file to recruit? or just like a name of the function to run. Also, these are called "rules" for some reason. Tutaj też można dodać prerequesites (files that need to exist) for running this function. Jeśli zrobimy z tego file, to zostanie on stworzony jeśli nie istnieje z początku	
	@echo "ass"

#for analyzing all the files in given folder; make fastqcMake; put all needed scripts into "Programs" folder; questions: czy napewno dobre kodowanie?,
#Struktura folderów: Programs, Genomes, #foldery_z_rawseqfiles#

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

PrepareMake: 
#PREPARE FOLDERS
	export CORES=32 #This is global variable
	mkdir Fastqc_Raw Trim_Galore Diversity_Cut Fastq_Screen_Div_Cut
	mkdir ./Trim_Galore/Trim_Galore_Raports ./Diversity_Cut/Fastqc_Trimmed
	cp ../Programs/FastQC_aggregate.sh ./Fastqc_Raw/; #cd Fastqc_Raw; FastQC_aggregate.sh; cd ..;
	cp ../Programs/FastQC_aggregate.sh ./Diversity_Cut/Fastqc_Trimmed;
	cp ../Programs/trimRRBSdiversityAdaptCustomers.py ./Trim_Galore/;
#PREPARE FOLDERS
FastqcMake:
#FASTQC
	fastqc --threads $(CORES) --outdir ./Fastqc_Raw *fastq*
	# Aggregate the results https://gist.github.com/danielecook/8e9afb2d2df7752efd8a#file-fastqc_aggregate-sh
#FASTQC

TrimMake:
#TRIM_GALORE
	ls *R1* | sort >> r1; ls *R3* | sort >> r2; paste r1 r2 >> read_pairs; rm r1 r2 #Outputs list of paired reads
	parallel --bar --colsep '\t' trim_galore --paired  --retain_unpaired --output_dir ./Trim_Galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC {1} {2}  :::: read_pairs #Using paired, columned list we pair the names
	mv ./Trim_Galore/*trimming_report.txt ./Trim_Galore/Trim_Galore_Raports; rm read_pairs
#TRIM_GALORE
#DIVERSITY CUTTING
	cd Trim_Galore
	ls *R1*val* | sort >> r1; ls *R3*val* | sort >> r2; paste r1 r2 >> read_pairs; rm r1 r2 #Outputs list of paired reads
	parallel --bar --colsep '\t' python trimRRBSdiversityAdaptCustomers.py -1 {1} -2 {2} 2>error :::: read_pairs
	mv *trimmed* ../Diversity_Cut; rm read_pairs; cd ../Diversity_Cut
	fastqc --threads 4 --outdir ./Fastqc_Trimmed *trimmed*; cd ..;
#DIVERSITY CUTTING

FastqcScreenMake:
#FASTQ SCREEN 
#needs bowtie2 and bismark; needs to configure fastq_screen.config; download contaminating genomes (downloaded from ncbi Assembly?); download contaminating adapters (see config file to see from where); for bismark option you have to set config databases to "Bisulfite_Genome" folders created by bismark_genome_preparation tool;
	cd Diversity_Cut
	../Programs/fastq_screen_v0.11.1/fastq_screen --aligner bowtie2 --bisulfite -outdir ../Fastq_Screen_Div_Cut/ *trimmed*
	
#FASTQ SCREEN

BismarkGenomeMake:
#MAKE CONVERTRD GENOMES
	cp -R ../Programs/*ismark-* ../Genomes
	../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ../Genomes/Rat_6
	../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ../Genomes/Mouse_38_6
	../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ../Genomes/Human_38_12
	../Programs/*ismark-*/bismark_genome_preparation --bowtie2 --verbose --yes ../Genomes/Ecoli

#MAKE CONVERTED GENOMES









