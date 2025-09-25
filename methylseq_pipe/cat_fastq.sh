#!/bin/bash -l
#SBATCH -J cat_fastq
#SBATCH -o cat_fastq.out
#SBATCH -e cat_fastq.error
#SBATCH --mail-user "your email"
#SBATCH --mail-type=ALL
#SBATCH -t 00-10:00:00
#SBATCH -A "your project ID"
#SBATCH -p core
#SBATCH -n 2


while IFS= read -r sample
do

cat $(grep -w "$sample"  NF_sample_sheet.csv | cut -f2 -d ",") > fastq_tmp/$sample.R1.fastq.gz &
cat $(grep -w "$sample"  NF_sample_sheet.csv | cut -f3 -d ",") > fastq_tmp/$sample.R2.fastq.gz &

wait
done < "samples.list"
