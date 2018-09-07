ingularity RRBS_Singularity.simg:/Analysis> snakemake -j 30 --printshellcmds --resources mem_mb=100000 --verbose
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 30
Rules claiming more threads will be scaled down.
Provided resources: mem_mb=100000
Job counts:
        count   jobs
        5       Bismark
        5       Calling
        5       Deduplication
        5       Div_Fastqc
        5       Diversity
        1       Multiqc
        15      Raw_Fastqc
        5       Trim
        1       all
        47
Resources before job selection: {'mem_mb': 100000, '_cores': 30, '_nodes': 9223372036854775807}
Ready jobs (20):
        Trim
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Trim
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Trim
        Raw_Fastqc
        Trim
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Trim
Selected jobs (3):
        Trim
        Trim
        Trim
Resources after job selection: {'mem_mb': 4000, '_cores': 0, '_nodes': 9223372036854775804}

[Thu Sep  6 19:15:18 2018]
rule Trim:
    input: B_Control-NM_86_S1_R1_001.fastq.gz, B_Control-NM_86_S1_R3_001.fastq.gz
    output: Trim_Galore/B_Control-NM_86_S1_R1_001_val_1.fq.gz, Trim_Galore/B_Control-NM_86_S1_R3_001_val_2.fq.gz, Reports/Trim_Galore/B_Control-NM_86_S1_R1_001.fastq.gz_trimming_report.txt, Reports/Trim_Galore/B_Control-NM_86_S1_R3_001.fastq.gz_trimming_report.txt
    jobid: 8
    wildcards: b=B_Control-NM_86_S1
    threads: 10
    resources: mem_mb=32000

/TrimGalore*/trim_galore --paired --output_dir ./Trim_Galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC B_Control-NM_86_S1_R1_001.fastq.gz B_Control-NM_86_S1_R3_001.fastq.gz &>> Reports/Trim_Galore/Trim_Log.txt || error_exit && mv Trim_Galore/B_Control-NM_86_S1_R1_001.fastq.gz_trimming_report.txt Trim_Galore/B_Control-NM_86_S1_R3_001.fastq.gz_trimming_report.txt Reports/Trim_Galore/ || error_exit

[Thu Sep  6 19:15:18 2018]
rule Trim:
    input: B_Control-NM_34_S9_R1_001.fastq.gz, B_Control-NM_34_S9_R3_001.fastq.gz
    output: Trim_Galore/B_Control-NM_34_S9_R1_001_val_1.fq.gz, Trim_Galore/B_Control-NM_34_S9_R3_001_val_2.fq.gz, Reports/Trim_Galore/B_Control-NM_34_S9_R1_001.fastq.gz_trimming_report.txt, Reports/Trim_Galore/B_Control-NM_34_S9_R3_001.fastq.gz_trimming_report.txt
    jobid: 29
    wildcards: b=B_Control-NM_34_S9
    threads: 10
    resources: mem_mb=32000

/TrimGalore*/trim_galore --paired --output_dir ./Trim_Galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC B_Control-NM_34_S9_R1_001.fastq.gz B_Control-NM_34_S9_R3_001.fastq.gz &>> Reports/Trim_Galore/Trim_Log.txt || error_exit && mv Trim_Galore/B_Control-NM_34_S9_R1_001.fastq.gz_trimming_report.txt Trim_Galore/B_Control-NM_34_S9_R3_001.fastq.gz_trimming_report.txt Reports/Trim_Galore/ || error_exit

[Thu Sep  6 19:15:18 2018]
rule Trim:
    input: B_Control-NM_42_S13_R1_001.fastq.gz, B_Control-NM_42_S13_R3_001.fastq.gz
    output: Trim_Galore/B_Control-NM_42_S13_R1_001_val_1.fq.gz, Trim_Galore/B_Control-NM_42_S13_R3_001_val_2.fq.gz, Reports/Trim_Galore/B_Control-NM_42_S13_R1_001.fastq.gz_trimming_report.txt, Reports/Trim_Galore/B_Control-NM_42_S13_R3_001.fastq.gz_trimming_report.txt
    jobid: 46
    wildcards: b=B_Control-NM_42_S13
    threads: 10
    resources: mem_mb=32000

/TrimGalore*/trim_galore --paired --output_dir ./Trim_Galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC B_Control-NM_42_S13_R1_001.fastq.gz B_Control-NM_42_S13_R3_001.fastq.gz &>> Reports/Trim_Galore/Trim_Log.txt || error_exit && mv Trim_Galore/B_Control-NM_42_S13_R1_001.fastq.gz_trimming_report.txt Trim_Galore/B_Control-NM_42_S13_R3_001.fastq.gz_trimming_report.txt Reports/Trim_Galore/ || error_exit
[Thu Sep  6 19:45:10 2018]
Finished job 29.
1 of 47 steps (2%) done
Resources before job selection: {'mem_mb': 36000, '_cores': 10, '_nodes': 9223372036854775805}
Ready jobs (18):
        Diversity
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Trim
        Raw_Fastqc
        Trim
        Raw_Fastqc
Selected jobs (1):
        Trim
Resources after job selection: {'mem_mb': 4000, '_cores': 0, '_nodes': 9223372036854775804}

[Thu Sep  6 19:45:10 2018]
rule Trim:
    input: B_Control-NM_18_S5_R1_001.fastq.gz, B_Control-NM_18_S5_R3_001.fastq.gz
    output: Trim_Galore/B_Control-NM_18_S5_R1_001_val_1.fq.gz, Trim_Galore/B_Control-NM_18_S5_R3_001_val_2.fq.gz, Reports/Trim_Galore/B_Control-NM_18_S5_R1_001.fastq.gz_trimming_report.txt, Reports/Trim_Galore/B_Control-NM_18_S5_R3_001.fastq.gz_trimming_report.txt
    jobid: 36
    wildcards: b=B_Control-NM_18_S5
    threads: 10
    resources: mem_mb=32000

/TrimGalore*/trim_galore --paired --output_dir ./Trim_Galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC B_Control-NM_18_S5_R1_001.fastq.gz B_Control-NM_18_S5_R3_001.fastq.gz &>> Reports/Trim_Galore/Trim_Log.txt || error_exit && mv Trim_Galore/B_Control-NM_18_S5_R1_001.fastq.gz_trimming_report.txt Trim_Galore/B_Control-NM_18_S5_R3_001.fastq.gz_trimming_report.txt Reports/Trim_Galore/ || error_exit
[Thu Sep  6 19:48:07 2018]
Finished job 46.
2 of 47 steps (4%) done
Resources before job selection: {'mem_mb': 36000, '_cores': 10, '_nodes': 9223372036854775805}
Ready jobs (18):
        Diversity
        Diversity
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Trim
        Raw_Fastqc
Selected jobs (1):
        Diversity
Resources after job selection: {'mem_mb': 4000, '_cores': 0, '_nodes': 9223372036854775804}

[Thu Sep  6 19:48:07 2018]
rule Diversity:
    input: Trim_Galore/B_Control-NM_42_S13_R1_001_val_1.fq.gz, Trim_Galore/B_Control-NM_42_S13_R3_001_val_2.fq.gz
    output: Diversity/B_Control-NM_42_S13_R1_001_val_1.fq_trimmed.fq.gz, Diversity/B_Control-NM_42_S13_R3_001_val_2.fq_trimmed.fq.gz
    jobid: 41
    wildcards: b=B_Control-NM_42_S13
    threads: 10
    resources: mem_mb=32000

python /trimRRBSdiversityAdaptCustomers.py -1 Trim_Galore/B_Control-NM_42_S13_R1_001_val_1.fq.gz -2 Trim_Galore/B_Control-NM_42_S13_R3_001_val_2.fq.gz &>> Reports/Diversity_Log.txt || error_exit && mv Trim_Galore/B_Control-NM_42_S13_R1_001_val_1.fq_trimmed.fq.gz Trim_Galore/B_Control-NM_42_S13_R3_001_val_2.fq_trimmed.fq.gz ./Diversity/ || error_exit
[Thu Sep  6 19:56:27 2018]
Finished job 8.
3 of 47 steps (6%) done
Resources before job selection: {'mem_mb': 36000, '_cores': 10, '_nodes': 9223372036854775805}
Ready jobs (18):
        Diversity
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Diversity
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Trim
        Raw_Fastqc
Selected jobs (1):
        Diversity
Resources after job selection: {'mem_mb': 4000, '_cores': 0, '_nodes': 9223372036854775804}

[Thu Sep  6 19:56:27 2018]
rule Diversity:
    input: Trim_Galore/B_Control-NM_86_S1_R1_001_val_1.fq.gz, Trim_Galore/B_Control-NM_86_S1_R3_001_val_2.fq.gz
    output: Diversity/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed.fq.gz, Diversity/B_Control-NM_86_S1_R3_001_val_2.fq_trimmed.fq.gz
    jobid: 25
    wildcards: b=B_Control-NM_86_S1
    threads: 10
    resources: mem_mb=32000

python /trimRRBSdiversityAdaptCustomers.py -1 Trim_Galore/B_Control-NM_86_S1_R1_001_val_1.fq.gz -2 Trim_Galore/B_Control-NM_86_S1_R3_001_val_2.fq.gz &>> Reports/Diversity_Log.txt || error_exit && mv Trim_Galore/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed.fq.gz Trim_Galore/B_Control-NM_86_S1_R3_001_val_2.fq_trimmed.fq.gz ./Diversity/ || error_exit
[Thu Sep  6 20:14:51 2018]
Finished job 36.
4 of 47 steps (9%) done
Resources before job selection: {'mem_mb': 36000, '_cores': 10, '_nodes': 9223372036854775805}
Ready jobs (18):
        Diversity
        Diversity
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Trim
        Raw_Fastqc
Selected jobs (1):
        Trim
Resources after job selection: {'mem_mb': 4000, '_cores': 0, '_nodes': 9223372036854775804}

[Thu Sep  6 20:14:51 2018]
rule Trim:
    input: B_Control-NM_9_S1_R1_001.fastq.gz, B_Control-NM_9_S1_R3_001.fastq.gz
    output: Trim_Galore/B_Control-NM_9_S1_R1_001_val_1.fq.gz, Trim_Galore/B_Control-NM_9_S1_R3_001_val_2.fq.gz, Reports/Trim_Galore/B_Control-NM_9_S1_R1_001.fastq.gz_trimming_report.txt, Reports/Trim_Galore/B_Control-NM_9_S1_R3_001.fastq.gz_trimming_report.txt
    jobid: 9
    wildcards: b=B_Control-NM_9_S1
    threads: 10
    resources: mem_mb=32000

/TrimGalore*/trim_galore --paired --output_dir ./Trim_Galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC B_Control-NM_9_S1_R1_001.fastq.gz B_Control-NM_9_S1_R3_001.fastq.gz &>> Reports/Trim_Galore/Trim_Log.txt || error_exit && mv Trim_Galore/B_Control-NM_9_S1_R1_001.fastq.gz_trimming_report.txt Trim_Galore/B_Control-NM_9_S1_R3_001.fastq.gz_trimming_report.txt Reports/Trim_Galore/ || error_exit
[Thu Sep  6 20:39:13 2018]
Finished job 9.
5 of 47 steps (11%) done
Resources before job selection: {'mem_mb': 36000, '_cores': 10, '_nodes': 9223372036854775805}
Ready jobs (18):
        Diversity
        Diversity
        Diversity
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
Selected jobs (1):
        Diversity
Resources after job selection: {'mem_mb': 4000, '_cores': 0, '_nodes': 9223372036854775804}

[Thu Sep  6 20:39:13 2018]
rule Diversity:
    input: Trim_Galore/B_Control-NM_34_S9_R1_001_val_1.fq.gz, Trim_Galore/B_Control-NM_34_S9_R3_001_val_2.fq.gz
    output: Diversity/B_Control-NM_34_S9_R1_001_val_1.fq_trimmed.fq.gz, Diversity/B_Control-NM_34_S9_R3_001_val_2.fq_trimmed.fq.gz
    jobid: 6
    wildcards: b=B_Control-NM_34_S9
    threads: 10
    resources: mem_mb=32000

python /trimRRBSdiversityAdaptCustomers.py -1 Trim_Galore/B_Control-NM_34_S9_R1_001_val_1.fq.gz -2 Trim_Galore/B_Control-NM_34_S9_R3_001_val_2.fq.gz &>> Reports/Diversity_Log.txt || error_exit && mv Trim_Galore/B_Control-NM_34_S9_R1_001_val_1.fq_trimmed.fq.gz Trim_Galore/B_Control-NM_34_S9_R3_001_val_2.fq_trimmed.fq.gz ./Diversity/ || error_exit
[Thu Sep  6 21:30:20 2018]
Finished job 41.
6 of 47 steps (13%) done
Resources before job selection: {'mem_mb': 36000, '_cores': 10, '_nodes': 9223372036854775805}
Ready jobs (19):
        Diversity
        Diversity
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Div_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Bismark
        Raw_Fastqc
        Raw_Fastqc
Selected jobs (1):
        Diversity
Resources after job selection: {'mem_mb': 4000, '_cores': 0, '_nodes': 9223372036854775804}

[Thu Sep  6 21:30:20 2018]
rule Diversity:
    input: Trim_Galore/B_Control-NM_18_S5_R1_001_val_1.fq.gz, Trim_Galore/B_Control-NM_18_S5_R3_001_val_2.fq.gz
    output: Diversity/B_Control-NM_18_S5_R1_001_val_1.fq_trimmed.fq.gz, Diversity/B_Control-NM_18_S5_R3_001_val_2.fq_trimmed.fq.gz
    jobid: 42
    wildcards: b=B_Control-NM_18_S5
    threads: 10
    resources: mem_mb=32000

python /trimRRBSdiversityAdaptCustomers.py -1 Trim_Galore/B_Control-NM_18_S5_R1_001_val_1.fq.gz -2 Trim_Galore/B_Control-NM_18_S5_R3_001_val_2.fq.gz &>> Reports/Diversity_Log.txt || error_exit && mv Trim_Galore/B_Control-NM_18_S5_R1_001_val_1.fq_trimmed.fq.gz Trim_Galore/B_Control-NM_18_S5_R3_001_val_2.fq_trimmed.fq.gz ./Diversity/ || error_exit
[Thu Sep  6 22:05:40 2018]
Finished job 25.
7 of 47 steps (15%) done
Resources before job selection: {'mem_mb': 36000, '_cores': 10, '_nodes': 9223372036854775805}
Ready jobs (20):
        Diversity
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Div_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Bismark
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Raw_Fastqc
        Bismark
        Div_Fastqc
        Raw_Fastqc
        Raw_Fastqc
Selected jobs (1):
        Bismark
Resources after job selection: {'mem_mb': 4000, '_cores': 0, '_nodes': 9223372036854775804}

[Thu Sep  6 22:05:40 2018]
rule Bismark:
    input: Diversity/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed.fq.gz, Diversity/B_Control-NM_86_S1_R3_001_val_2.fq_trimmed.fq.gz
    output: Bismark_Raw/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam, Reports/Bismark_Raw/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt
    jobid: 5
    wildcards: b=B_Control-NM_86_S1
    threads: 10
    resources: mem_mb=32000

/Bismark*/bismark --bowtie2 --genome_folder ./Genome/ -1 Diversity/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed.fq.gz -2 Diversity/B_Control-NM_86_S1_R3_001_val_2.fq_trimmed.fq.gz &>> Reports/Bismark_Raw/Bismark_Log.txt || error_exit && /Bismark*/bam2nuc --genome_folder ./Genome/ B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam &>> Reports/Bismark_Raw/Bam2nuc_Log.txt || error_exit && /Bismark*/bismark2summary B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2sum_Log.txt && /Bismark*/bismark2report B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2rep_Log.txt && /Bamqc/bamqc B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bamqc_Raw_Log.txt && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam ./Bismark_Raw/ && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.nucleotide_stats.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.zip ./Reports/Bismark_Raw/
/bin/bash: error_exit: command not found
/bin/bash: error_exit: command not found
Full Traceback (most recent call last):
  File "/usr/local/lib/python3.6/dist-packages/snakemake/executors.py", line 1313, in run_wrapper
    singularity_args, use_singularity, None, jobid, is_shell)
  File "/Analysis/Snakefile", line 218, in __rule_Bismark
  File "/usr/local/lib/python3.6/dist-packages/snakemake/shell.py", line 133, in __new__
    raise sp.CalledProcessError(retcode, cmd)
subprocess.CalledProcessError: Command ' set -euo pipefail;  /Bismark*/bismark --bowtie2 --genome_folder ./Genome/ -1 Diversity/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed.fq.gz -2 Diversity/B_Control-NM_86_S1_R3_001_val_2.fq_trimmed.fq.gz &>> Reports/Bismark_Raw/Bismark_Log.txt || error_exit && /Bismark*/bam2nuc --genome_folder ./Genome/ B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam &>> Reports/Bismark_Raw/Bam2nuc_Log.txt || error_exit && /Bismark*/bismark2summary B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2sum_Log.txt && /Bismark*/bismark2report B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2rep_Log.txt && /Bamqc/bamqc B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bamqc_Raw_Log.txt && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam ./Bismark_Raw/ && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.nucleotide_stats.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.zip ./Reports/Bismark_Raw/ ' returned non-zero exit status 127.

    [Thu Sep  6 22:05:41 2018]
    Error in rule Bismark:
        jobid: 5
        output: Bismark_Raw/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam, Reports/Bismark_Raw/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt

Full Traceback (most recent call last):
  File "/usr/local/lib/python3.6/dist-packages/snakemake/executors.py", line 1313, in run_wrapper
    singularity_args, use_singularity, None, jobid, is_shell)
  File "/Analysis/Snakefile", line 218, in __rule_Bismark
  File "/usr/local/lib/python3.6/dist-packages/snakemake/shell.py", line 133, in __new__
    raise sp.CalledProcessError(retcode, cmd)
subprocess.CalledProcessError: Command ' set -euo pipefail;  /Bismark*/bismark --bowtie2 --genome_folder ./Genome/ -1 Diversity/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed.fq.gz -2 Diversity/B_Control-NM_86_S1_R3_001_val_2.fq_trimmed.fq.gz &>> Reports/Bismark_Raw/Bismark_Log.txt || error_exit && /Bismark*/bam2nuc --genome_folder ./Genome/ B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam &>> Reports/Bismark_Raw/Bam2nuc_Log.txt || error_exit && /Bismark*/bismark2summary B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2sum_Log.txt && /Bismark*/bismark2report B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2rep_Log.txt && /Bamqc/bamqc B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bamqc_Raw_Log.txt && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam ./Bismark_Raw/ && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.nucleotide_stats.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.zip ./Reports/Bismark_Raw/ ' returned non-zero exit status 127.

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/usr/local/lib/python3.6/dist-packages/snakemake/executors.py", line 346, in _callback
    raise ex
  File "/usr/lib/python3.6/concurrent/futures/thread.py", line 56, in run
    result = self.fn(*self.args, **self.kwargs)
  File "/usr/local/lib/python3.6/dist-packages/snakemake/executors.py", line 1325, in run_wrapper
    show_traceback=True))
snakemake.exceptions.RuleException: CalledProcessError in line 119 of /Analysis/Snakefile:
Command ' set -euo pipefail;  /Bismark*/bismark --bowtie2 --genome_folder ./Genome/ -1 Diversity/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed.fq.gz -2 Diversity/B_Control-NM_86_S1_R3_001_val_2.fq_trimmed.fq.gz &>> Reports/Bismark_Raw/Bismark_Log.txt || error_exit && /Bismark*/bam2nuc --genome_folder ./Genome/ B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam &>> Reports/Bismark_Raw/Bam2nuc_Log.txt || error_exit && /Bismark*/bismark2summary B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2sum_Log.txt && /Bismark*/bismark2report B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2rep_Log.txt && /Bamqc/bamqc B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bamqc_Raw_Log.txt && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam ./Bismark_Raw/ && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.nucleotide_stats.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.zip ./Reports/Bismark_Raw/ ' returned non-zero exit status 127.
  File "/Analysis/Snakefile", line 119, in __rule_Bismark

RuleException:
CalledProcessError in line 119 of /Analysis/Snakefile:
Command ' set -euo pipefail;  /Bismark*/bismark --bowtie2 --genome_folder ./Genome/ -1 Diversity/B_Control-NM_86_S1_R1_001_val_1.fq_trimmed.fq.gz -2 Diversity/B_Control-NM_86_S1_R3_001_val_2.fq_trimmed.fq.gz &>> Reports/Bismark_Raw/Bismark_Log.txt || error_exit && /Bismark*/bam2nuc --genome_folder ./Genome/ B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam &>> Reports/Bismark_Raw/Bam2nuc_Log.txt || error_exit && /Bismark*/bismark2summary B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2sum_Log.txt && /Bismark*/bismark2report B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bismark2rep_Log.txt && /Bamqc/bamqc B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam  &>> Reports/Bismark_Raw/Bamqc_Raw_Log.txt && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam ./Bismark_Raw/ && mv B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe.nucleotide_stats.txt B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.html B_Control-NM_86_S1_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.zip ./Reports/Bismark_Raw/ ' returned non-zero exit status 127.
  File "/Analysis/Snakefile", line 119, in __rule_Bismark
  File "/usr/lib/python3.6/concurrent/futures/thread.py", line 56, in run
[Thu Sep  6 22:08:09 2018]
Finished job 6.
8 of 47 steps (17%) done
[Thu Sep  6 22:58:49 2018]
Finished job 42.
9 of 47 steps (19%) done
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: /Analysis/.snakemake/log/2018-09-06T191518.240813.snakemake.log
unlocking
removing lock
removing lock
removed all locks
