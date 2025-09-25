#!/bin/bash -l
#SBATCH -J DMR_prep
#SBATCH -o DMR_prep.output
#SBATCH -e DMR_prep.error
#SBATCH --mail-user jesper.boman@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 00-01:00:00
#SBATCH -A "project ID"
#SBATCH -p core
#SBATCH -n 5

#THis script is used to make bsseq formatted data from the coverage2cytosine reports retrieved from the Bismark pipeline

grep "GTcoll19\|GTcoll18" samples_info.txt | cut -f1 > wild.samples.list

dir="EpiClock/methylseq_pipe/Results/bismark/coverage2cytosine/reports"

mkdir bsseq_formatted_data


while IFS= read -r sample
do


zcat $dir/$sample.CpG_report.txt.gz | awk '{if($6 == "CG" && prevPos == $2-1){CpG++} else{CpG=0}; if($6 == "CG"  && prevPos == $2-1 && $1 == prevChr && CpG == 2  ){CpG=0; print $1 "\t" $2 "\t" $4+prevM "\t" $4+$5+prevM+prevC}; if($6 == "CG"){prevM=$4; prevC=$5; prevChr=$1; prevPos=$2; CpG++}}' > bsseq_formatted_data/${sample}_data_formatted &

p=$(($p+1))
remainder=$(( p  % 5 ))
if [[ $remainder -eq 0 ]] ; then
wait
fi


done <"wild.samples.18.list"

wait
