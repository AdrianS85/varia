#!/usr/bin/env nextflow

/* 
wget -qO- https://get.nextflow.io | bash
-with-report -with-trace -with-timeline -with-dag
ps, date, sed, grep, egrep, awk
sendMail

https://groups.google.com/forum/#!topic/nextflow/a1fTBd1bPYw
output:
         file '*_fastqc.zip' into fastqc_results
nextflow run RRBS.nf
//TG_out1 = set file(d1), file(d2) from TG_out  
//TG_out.subscribe { println "TG_out : $it" } 
//RP.subscribe { println "RP : $it" }
*/



params.all = "*.fastq.gz"
ALL = file(params.all).flatten()
R1 = Channel.fromPath("*_R1_001.fastq.gz").toSortedList().flatten()
R3 = Channel.fromPath("*_R3_001.fastq.gz").toSortedList().flatten()
RP = Channel.fromFilePairs("*_R{1,2,3}_001.fastq.gz", size: -1, flat: true) //https://groups.google.com/forum/#!topic/nextflow/X4YyYmLTbTo

params.outputFR = "Fastqc_Raw"
params.outputTG = "Trim_Galore"
params.outputD = "Diversity"
params.outputFD = "Fastqc_Div"
params.outputB = "Bismark"

/*
process Fastqc_Raw {
   publishDir path: params.outputFR, mode: 'copy' //INFO: for some reason path dir can only be established by referencing parameter  

   input:
        file ALL

    output:
         file '*_fastqc.{zip,html}' //INFO: {} brackets include OR operator
 
    """
    fastqc $ALL
    """
}*/



process Trim_Galore {
   publishDir path: params.outputTG, mode: 'copy'  

    tag { ID }

   input:
        set val(ID), file(tg1), file(tg2), file(tg3) from RP
        
    output:
        set ID, file("${ID}_R1_001_val_1.fq.gz"), file("${ID}_R3_001_val_2.fq.gz") into TG_out //INFO: we can transfer values from input to output by directly referencing value: "ID"
        file "${ID}_R{1,3}_001.fastq.gz_trimming_report.txt"
    """
    trim_galore --paired -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC ${tg1} ${tg3}
    """
}



process Diversity {
   publishDir path: params.outputD, mode: 'copy'  

   input:
        set val(ID), file(d1), file(d2) from TG_out	

    output:
        set ID, file("${ID}_R1_001_val_1.fq_trimmed.fq.gz"), file("${ID}_R3_001_val_2.fq_trimmed.fq.gz") into D_out
 
    """
    python /home/adrian/bin/trimRRBSdiversityAdaptCustomers.py -1 $d1 -2 $d2
    """
}

//D_out.subscribe { println "D_out : $it" } 




process Fastqc_Div {
   publishDir path: params.outputFD, mode: 'copy' 

   input:
        set val(ID), file(fd1), file(fd2) from (D_out, D_outB)

    output:
         file '*_fastqc.{zip,html}'
 
    """
    fastqc $fd1; fastqc $fd2
    """
}



process Bismark {
   publishDir path: params.outputB, mode: 'copy' 

   input:
        set val(ID), file(b1), file(b2) from D_outB	

    output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt") into B_out
 
    """
    /home/adrian/bismark_v0.18.1/bismark --bowtie2 --genome_folder /media/adrian/One1/Nextflow/Analysis/Genome/ -1 $b1 -2 $b2
    """
}


