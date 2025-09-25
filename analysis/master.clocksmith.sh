#!/bin/bash -l

outdir="varying_seed_CDS_dinuc"
mkdir $outdir
inputdata="bigdata.comb.10.noW.CDSonly"
id="bigdata.comb.10.noW.CDSonly"


unset to #If starting from anything but 1, then "to" needs to be set to a relevant variable
unset fro

for i in $(seq 1 100);
do


fro=$(( to+1 ))
to=$(( 10*i ))


sbatch -J smither.$id.$fro.$to -o slurm/smither.$id.$fro.$to.output -e slurm/smither.$id.$fro.$to.error -A "project_ID" -t 24:00:00 -p core -n 5 smither.sh $fro $to $outdir $inputdata

done
