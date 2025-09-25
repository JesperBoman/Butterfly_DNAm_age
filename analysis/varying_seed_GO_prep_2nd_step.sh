#!/bin/bash -l


for i in $(seq 1 1000);
do

paste -d "\t"  varying_seed_CDS_dinuc.$i.CDS.geneIDs.list varying_seed_CDS_dinuc.$i.CDS.loci.list >> vs_geneID_and_loci.list

done

sort vs_geneID_and_loci.list | uniq -c > vs_geneID_and_loci_counts.list


awk '{if($1 >= 100) print $2}' vs_geneID_and_loci_counts.list > vs_geneID_atLeast100seeds.list
awk '{if($1 >= 500) print $2}' vs_geneID_and_loci_counts.list > vs_geneID_atLeast500seeds.list
