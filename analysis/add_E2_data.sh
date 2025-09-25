#!/bin/bash -l
#SBATCH -J add_E2_data
#SBATCH -o add_E2_data.output
#SBATCH -e add_E2_data.error
#SBATCH --mail-user "your email"
#SBATCH --mail-type=ALL
#SBATCH -t 00-01:00:00
#SBATCH -A "your project ID"
#SBATCH -p core

#Also known as add_older_data.sh
#This script combines the data from (comb_data.sh) with data from the E2 experiment (see Boman et al 2023: https://onlinelibrary.wiley.com/doi/10.1111/mec.16957)

# HDAL, LDAL and HDLI: Age = 0 --> added sample size: 12 (9 females, 3 males)

# AL, LI: Age = 5 --> Added sample size: 6



refdir="EpiClock/reference"
inputdir="VanessaDNAm/stranded_CpG_reports"


while IFS= read -r sample_double
do

sample=$(cut -f1 <(echo $sample_double) -d " ")
age=$(cut -f2 <(echo $sample_double) -d " ")

date
echo $sample

report=$(grep $sample"_" <(ls $inputdir))
echo $report


awk -v sample=$sample -v age=$age 'NR==FNR{a[$2]=$1} NR!=FNR && (FNR % 2 != 0) {meth=$4; unmeth=$5} NR!=FNR && (FNR % 2 == 0){cov=(meth+$4+unmeth+$5); if(cov > 9){print a[$1] "\t" $2-1 "\t" (meth+$4)/(meth+$4+unmeth+$5) "\t" sample "\t" age}}' $refdir/chrom.match.list <(zcat $inputdir/$report) >> bigdata.comb.10

done < "earlyStudySamples_and_age_list.txt"


grep -v "Chr_W" bigdata.comb.10 > bigdata.comb.10.noW
