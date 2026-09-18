#!/bin/bash -l

fro=$1
to=$2

module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

outdir=$3
inputdata=$4
samplesize=$5
threshold_age=11

for i in $(seq $fro $to);
do

#seed outdir inputdata samplesize threshold_age

#When subsampling the E.clock
Rscript fancier_varying_seed_clock5_F_SUBSAMPLE.R $i $outdir $inputdata $samplesize $threshold_age

#When subsampling the W.clock
Rscript fancier_varying_seed_clock5_F_WCLOCK_SUBSAMPLE.R $i $outdir $inputdata $samplesize $threshold_age


done
