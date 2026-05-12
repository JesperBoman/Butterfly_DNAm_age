#!/bin/bash -l

fro=$1
to=$2

module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

outdir=$3
inputdata=$4
samplesize=$5

for i in $(seq $fro $to);
do

#seed outdir inputdata

Rscript fancier_varying_seed_clock_SUBSAMPLE.R $i $outdir $inputdata $samplesize

wait


done
