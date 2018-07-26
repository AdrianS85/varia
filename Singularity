BootStrap: docker
From: ubuntu:16.04
MirrorURL: http://us.archive.ubuntu.com/ubuntu/

#sudo singularity build --writable test2.simg ./Singularity
#sudo singularity shell --bind ./Analysis:/Analysis --workdir /Analysis --writable test2.simg

%files
	ProgramsCopy/* /.
	/usr/lib/locale/locale-archive /usr/lib/locale/locale-archive
%post
	apt update
	apt -y install default-jdk perl bowtie2 snakemake python-pip
	pip install multiqc cutadapt
	mkdir Analysis
