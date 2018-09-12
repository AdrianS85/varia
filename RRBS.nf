#!/usr/bin/env nextflow

/* 
wget -qO- https://get.nextflow.io | bash
-with-report -with-trace -with-timeline -with-dag
ps, date, sed, grep, egrep, awk
sendMail
-resume
https://groups.google.com/forum/#!topic/nextflow/a1fTBd1bPYw

nextflow run RRBS.nf

//RP.subscribe { println "RP : $it" }


https://github.com/FelixKrueger/Bismark/issues/111

#export LANG=en_US.UTF-8
	#PWD=/root???
#locale-gen en_US.UTF-8

10.121.4.107
*/




params.all = "*.fastq.gz"
ALL = file(params.all).flatten()
R1 = Channel.fromPath("*_R1_001.fastq.gz").toSortedList().flatten()
R2 = Channel.fromPath("*_R2_001.fastq.gz").toSortedList().flatten()
R3 = Channel.fromPath("*_R3_001.fastq.gz").toSortedList().flatten()
RP = Channel.fromFilePairs("*_R{1,2,3}_001.fastq.gz", size: -1, flat: true) //https://groups.google.com/forum/#!topic/nextflow/X4YyYmLTbTo
RP2 = Channel.fromFilePairs("*_R2_001.fastq.gz", size: -1, flat: true) 

params.gen_fol = "/Analysis/Genome/"
BM_Switch = false
params.outputFR = "Fastqc_Raw"
params.outputFRM = "Fastqc_Raw/Multiqc"
params.outputTG = "Trim_Galore"
params.outputTGM = "Trim_Galore/Multiqc"
params.outputD = "Diversity"
params.outputFD = "Fastqc_Div"
params.outputFDM = "Fastqc_Div/Multiqc"
params.outputB = "Bismark"
params.outputBM = "Bismark/Multiqc"
params.outputDD = "Deduplication"
params.outputC = "Calling"




process Fastqc_Raw {
         publishDir path: params.outputFR, mode: 'copy' //INFO: for some reason path dir can only be established by referencing parameter  

         input:
         file ALL

         output:
         file '*_fastqc.{zip,html}' into (FR_out1, FR_out2) //INFO: {} brackets include OR operator

         """
         fastqc $ALL
         """
}

FR_out2.subscribe onComplete: {
         FR_out3 = Channel.create()
         FR_out3 = FR_out1.collect()

         process Fastqc_Raw_Multi {
                  publishDir path: params.outputFRM, mode: 'copy'

                  input:
                  file FR_out3

                  output:
                  set file("Raw_Fastqc_multiqc.html"), file("Raw_Fastqc_multiqc_data.tar.gz")

                  """
                  multiqc -n Raw_Fastqc_multiqc.html ${FR_out3}; tar czf Raw_Fastqc_multiqc_data.tar.gz Raw_Fastqc_multiqc_data
                  """
}}




process Trim_Galore {
         publishDir path: params.outputTG, mode: 'copy'  

         input:
         set val(ID), file(tg1), file(tg2), file(tg3) from RP

         output:
         set ID, file("${ID}_R1_001_val_1.fq.gz"), file("${ID}_R3_001_val_2.fq.gz") into (TG_out) //INFO: we can transfer values from input to output by directly referencing value: "ID"
         file "${ID}_R{1,3}_001.fastq.gz_trimming_report.txt" into (TG_out1, TG_out2)
         """
         trim_galore --paired -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC ${tg1} ${tg3}
         """
}

////WARNING! This will work only if Calling will be performed during given run
TG_out2.subscribe onComplete: {
         TG_out3 = Channel.create()
         TG_out3 = TG_out1.collect()

         process Trim_Multi {
         publishDir path: params.outputTGM, mode: 'copy'

         input:
         file TG_out3

         output:
         set file("Trim_multiqc.html"), file("Trim_multiqc_data.tar.gz")

         """
         multiqc -n Trim_multiqc.html ${TG_out3}; tar czf Trim_multiqc_data.tar.gz Trim_multiqc_data
         """
}}




process Diversity {
         publishDir path: params.outputD, mode: 'copy'  

         input:
         set val(ID), file(d1), file(d2) from TG_out	

         output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed.fq.gz"), file("${ID}_R3_001_val_2.fq_trimmed.fq.gz") into (D_out, D_outB)

         """
         python /trimRRBSdiversityAdaptCustomers.py -1 $d1 -2 $d2
         """
} // pythons scripts are not binaries and cannot be accesed via PATH




process Fastqc_Div {
         publishDir path: params.outputFD, mode: 'copy' 

         input:
         set val(ID), file(fd1), file(fd2) from D_out

         output:
         file '*_fastqc.{zip,html}' into (FD_out1, FD_out2)

         """
         fastqc $fd1; fastqc $fd2
         """
}

////WARNING! This will work only if Calling will be performed during given run
FD_out2.subscribe onComplete: {
         FD_out3 = Channel.create()
         FD_out3 = FD_out1.collect()

         process Fastqc_Div_Multi {
                  publishDir path: params.outputFDM, mode: 'copy'

                  input:
                  file FD_out3

                  output:
                  set file("Fastqc_Div_multiqc.html"), file("Fastqc_Div_multiqc_data.tar.gz")

                  """
                  multiqc -n Fastqc_Div_multiqc.html ${FD_out3}; tar czf Fastqc_Div_multiqc_data.tar.gz Fastqc_Div_multiqc_data
                  """
}}




process Bismark {
         publishDir path: params.outputB, mode: 'copy' 

         input:
         set val(ID), file(b1), file(b2) from D_outB	

         output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt") into (B_out, B_outDD, Rep_out1)

         """
         bismark --bowtie2 --genome_folder ${params.gen_fol}  -1 $b1 -2 $b2
         """ //gen_fol should be always aviable based on how the singularity is made and run
}




process Bismark_Report {
         publishDir path: params.outputB, mode: 'copy' 

         input:
         set val(ID), file(br1), file(br2) from B_out	

         output:
         set val(ID), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.nucleotide_stats.txt"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.html"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.zip") into Rep_out2

         """
         bam2nuc --genome_folder ${params.gen_fol} $br1; bamqc $br1
         """
}




BRP2_outDD = Channel.create()
BRP2_outDD = B_outDD.join(RP2)

process DeDuplication {
         publishDir path: params.outputDD, mode: 'copy' 

         input:
         set val(ID), file(dd1), file(dd2), file(dd3) from BRP2_outDD

         output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bam") into DD_out
         set file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted_bamqc.zip"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted_bamqc.html"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped_dup_log.txt") 

         """
         samtools view -h -o ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam ${dd1}; 
         strip_bismark_sam.sh ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam; 
         python /nugentechnologies-nudup-7a126eb/nudup.py -T /Analysis/tmp/ --paired-end -f ${dd3} -o ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sam; 
         samtools sort -n -o ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bam ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.bam;
         bamqc ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bam 
         """
}




process Calling {
         publishDir path: params.outputC, mode: 'copy' 
         publishDir path: params.outputB, mode: 'copy', pattern: "*txt" 

         input:
         set val(ID), file(c1) from DD_out	

         output:
         set file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bismark.cov.gz"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bedGraph.gz") 
         set val(ID), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.M-bias.txt"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted_splitting_report.txt") into (Rep_out3, Rep_out4)

         """
         bismark_methylation_extractor --bedGraph --paired-end --comprehensive --merge_non_CpG ${c1}
         """
}




////WARNING! This will work only if Calling will be performed during given run
Rep_out4.subscribe onComplete: {
         Rep_out5 = Channel.create()
         Rep_out12 = Channel.create()
         Rep_out125 = Channel.create()
         Rep_out5 = Rep_out3
         Rep_out12 = Rep_out1.join(Rep_out2)
         Rep_out125 = Rep_out12.join(Rep_out5)

         process Bismark_Report2 {
                  publishDir path: params.outputB, mode: 'copy'

                  input:
                  set val(ID), file(bam), file(rep), file(nucl), file(bqch), file(bqcz), file(mbias), file(split) from Rep_out125

                  output:
                  set val(ID), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.html") into BM_out

                  """
                  bismark2report --alignment_report ${rep} --splitting_report ${split} --mbias_report  ${mbias} --nucleotide_report ${nucl}  
                  """
}

         BM_out.subscribe onComplete:{
                  process Bismark_Multi {
                           echo true	

                           """
                           mkdir /Analysis/Bismark/Multiqc;
                           bismark2summary /Analysis/Bismark/*.bam; 
                           mv bismark_summary_report.html bismark_summary_report.txt /Analysis/Bismark/Multiqc;
                           multiqc -n Bismark_multiqc.html /Analysis/Bismark/.; 
                           mv -R /Analysis/Bismark/Bismark_multiqc.html /Analysis/Bismark/Bismark_multiqc_data /Analysis/Bismark/Multiqc/ 
                           """
}}}

workflow.onComplete {
         def msg = """\
                  Pipeline execution summary
                  ---------------------------
                  Completed at: ${workflow.complete}
                  Duration    : ${workflow.duration}
                  Success     : ${workflow.success}
                  workDir     : ${workflow.workDir}
                  exit status : ${workflow.exitStatus}
                  """
                  .stripIndent()

         sendMail(to: 'adrianstankiewicz85@gmail.com', subject: 'RRBS pipeline finished', body: msg)
}
