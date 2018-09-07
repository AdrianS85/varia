From: ubuntu:18.04
MirrorURL: http://us.archive.ubuntu.com/ubuntu/

#sudo singularity build --writable RRBS_Singularity.simg ./Singularity
#SANDBOX: sudo singularity build -s sandbox1.simg ./Singularity
#sudo singularity shell --cleanenv --bind ./Analysis:/Analysis --workdir /Analysis --writable RRBS_Singularity.simg
#File folder needs to be structured as such: Folder Analysis to be Binded with container. Inside this folder: 1)Files to be analyzed, .gz packed, 2)Folder "Genome" with Bismark index "Bisulfite_Genome", Reference genome e.g. "Mus_musculus_genomic_refseq_GRCm38.p5.fa"... also there is "doc_Mus_musculus_db_refseq.txt", 3) Snakefile?

%files
        ProgramsCopy/* /.
        usr/lib/locale/locale-archive /usr/lib/locale/locale-archive

%environment
        export PATH="/samtools-1.9:$PATH"

%post
        apt update
        apt -y install perl default-jdk bowtie2 python-pip python3-pip libcurl3
        pip install multiqc cutadapt
        pip3 install snakemake
        mkdir Analysis
        cd ../Analysis
