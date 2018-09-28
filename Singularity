From: ubuntu:18.04
MirrorURL: http://us.archive.ubuntu.com/ubuntu/

#sudo singularity build --writable RRBS_Singularity.simg ./Singularity
#SANDBOX: sudo singularity build -s sandbox1.simg ./Singularity
#sudo singularity shell --cleanenv --bind ./Analysis:/Analysis --workdir /Analysis --writable RRBS_Singularity.simg
#File folder needs to be structured as such: Folder Analysis to be Binded with container. Inside this folder: 1)Files to be analyzed, .gz packed, 2)Folder "Genome" with Bismark index "Bisulfite_Genome", Reference genome e.g. "Mus_musculus_genomic_refseq_GRCm38.p5.fa"... also there is "doc_Mus_musculus_db_refseq.txt", 3) Snakefile?

%environment
        export PATH="/samtools-1.9:$PATH"



%post
        apt update
        apt -y install vim wget perl unzip default-jdk bowtie2 python-pip libcurl3
        
        cd /
        
        wget https://github.com/nugentechnologies/NuMetRRBS/blob/master/strip_bismark_sam.sh
        mv strip_bismark_sam.sh /usr/bin
        
        wget https://github.com/nugentechnologies/NuMetRRBS/blob/master/trimRRBSdiversityAdaptCustomers.py
        
        
        wget -qO- https://get.nextflow.io | bash
        mv nextflow /usr/bin
        
        wget https://github.com/FelixKrueger/Bismark/archive/master.zip
        unzip master.zip
        export PATH="/Bismark-master:$PATH"
        rm master.zip
        
        wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
        unzip fastqc_v0.11.7.zip
        chmod 755 FastQC/fastqc
        export PATH="/FastQC:$PATH"
        rm fastqc_v0.11.7.zip
        
        wget https://github.com/FelixKrueger/TrimGalore/archive/master.zip
        unzip master.zip
        export PATH="/TrimGalore-master:$PATH"
        rm master.zip
        
        wget https://github.com/s-andrews/BamQC/archive/master.zip
        unzip master.zip
        export PATH="/BamQC-master:$PATH"
        chmod 755 BamQC-master/bamqc
        rm master.zip
        
        wget https://github.com/nugentechnologies/nudup/archive/master.zip
        unzip master.zip
        export PATH="/nudup-master:$PATH"
        rm master.zip
        
        pip install multiqc cutadapt

