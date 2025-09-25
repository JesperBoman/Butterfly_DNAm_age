#!/bin/bash -l
#SBATCH -J DMR_analysis_top_epiClock
#SBATCH -o DMR_analysis_top_epiClock.output
#SBATCH -e DMR_analysis_top_epiClock.error
#SBATCH --mail-user "your email"
#SBATCH --mail-type=ALL
#SBATCH -t 00-10:00:00
#SBATCH -A "project ID"
#SBATCH -p core
#SBATCH -n 20

ml bioinfo-tools R_packages/4.1.1




dir="/crex/proj/uppstore2017185/b2014034_nobackup/Jesper/EpiClock/DMRs"


Group1="Meconium"
Group2="Long_distance"

argsG1=$(ls bsseq_formatted_data | grep "M" | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')
argsG2=$(ls bsseq_formatted_data | grep "B\|H" | awk -v dir="$dir" '{printf dir "/" "bsseq_formatted_data" "/" $1 "\t"}')

#Step 1: Mostly just reading data into bsseq R format
Rscript DMR_analysis.R $argsG1 $argsG2

mv BS_data.rda BS_data_${Group1}v${Group2}.rda
mv BS_data.fit.rda BS_data_${Group1}v${Group2}.fit.rda

mv BS_data_${Group1}v${Group2}.rda bsseq_read_rda
mv BS_data_${Group1}v${Group2}.fit.rda bsseq_fit_rda

mkdir $dir/results/${Group1}v${Group2}

cd $dir/results/${Group1}v${Group2}

#Step 2: calling DMRs
Rscript ../../DMR_analysis_findDMRs.R $dir/bsseq_fit_rda/BS_data_${Group1}v${Group2}.fit.rda $dir/bsseq_read_rda/BS_data_${Group1}v${Group2}.rda




#NOTES
#10 cores not enough when smoothing 14 samples. Use 20 cores instead.
