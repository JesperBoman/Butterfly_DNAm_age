#!/bin/bash -l

#I run this pipeline in a screen sessions

ml bioinfo-tools Nextflow/latest
export NXF_HOME="/crex/proj/uppstore2017185/b2014034_nobackup/Jesper/EpiClock/methylseq_pipe"

NXF_OPTS='-Xms1g -Xmx4g'

genome="EpiClock/reference/GCA_905220365.1_ilVanCard2.1_genomic_chroms_SH.fasta"

#Lambda phage genome to check bisulfite conversion error rates
#genome="lambda_phage.fa"


nextflow run nf-core/methylseq -name "mtDNA" --input NF_samp_sheet_upd.csv -profile uppmax --project "your project ID" --max_cpus 20 --max_memory 128GB --fasta $genome --outdir ./Results --save_reference --comprehensive --cytosine_report --clip_r1 8 --clip_r2 8 --skip_preseq
