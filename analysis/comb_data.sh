#!/bin/bash -l
#SBATCH -J comb_data
#SBATCH -o comb_data.output
#SBATCH -e comb_data.error
#SBATCH --mail-user ""
#SBATCH --mail-type=ALL
#SBATCH -t 00-05:00:00
#SBATCH -A "your project ID"
#SBATCH -p core
#SBATCH -n 1


#Basically combines data from different individuals into one dataset
inputdir="EpiClock/methylseq_pipe/Results/bismark/coverage2cytosine/reports"


while IFS= read -r sample_double
do

sample=$(cut -f1 <(echo $sample_double) -d " " )
age=$(cut -f2 <(echo $sample_double) -d " "  )

date
echo $sample

awk -v sample=$sample -v age=$age 'NR % 2 != 0 {meth=$4; unmeth=$5} NR % 2 == 0{cov=(meth+$4+unmeth+$5); if(cov > 9){print $1 "\t" $2-1 "\t" (meth+$4)/(meth+$4+unmeth+$5) "\t" sample "\t" age}}' <(zcat $inputdir/$sample.CpG_report.txt.gz) >> data.comb.10

done < "samples_info.txt"
