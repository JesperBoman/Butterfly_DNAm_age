#!/bin/bash -l

fro=$1
to=$2

module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

outdir=$3
inputdata=$4
threshold_age=11

for i in $(seq $fro $to);
do

#For running the E.clock
#seed outdir inputdata threshold_age
Rscript fancier_varying_seed_clock5_F.R $i $outdir $inputdata $threshold_age

#For running the W.clock
#Rscript fancier_varying_seed_clock5_F_WCLOCK.R $i $outdir $inputdata $adult_age


done
