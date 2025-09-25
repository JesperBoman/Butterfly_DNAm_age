#!/bin/bash -l

fro=$1
to=$2

ml R_packages/4.3.1

outdir=$3
inputdata=$4


for i in $(seq $fro $to);
do

#seed outdir inputdata

Rscript fancier_varying_seed_clock.R $i $outdir $inputdata

wait


done
