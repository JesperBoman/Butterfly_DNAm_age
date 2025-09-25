#!/bin/bash -l
#SBATCH -J varying_seed_GO_prep
#SBATCH -o varying_seed_GO_prep.output
#SBATCH -e varying_seed_GO_prep.error
#SBATCH --mail-user "your_email"
#SBATCH --mail-type=ALL
#SBATCH -t 00-10:00:00
#SBATCH -A "project_ID"
#SBATCH -p core
#SBATCH -n 1


#Note: one possible interesting change that could be made is to count each gene only once per seed (to control for gene length)
#Or implement a calculation so that they get weighted by the total CDS length per gene

module load bioinfo-tools BEDTools/2.29.2

dir="EpiClock/analysis/varying_seed_CDS_dinuc"
mRNA_list="EpiClock/reference/mRNA.gff"
prom_list="EpiClock/reference/promoters.list"

name="varying_seed_CDS_dinuc"
fro=1
to=1000

for i in $(seq $fro $to);
do

cut -f1-2 $dir/meth.df.elastic_nonZero.df.$i | sort -u | awk '{print $1 "\t" $2-1 "\t" $2}' > tmp.sites



annot_double="CDS,EpiClock/reference/cds.bed"

shortAnnot=$(cut -f1 <(echo $annot_double) -d ",")
annot=$(cut -f2 <(echo $annot_double) -d ",")

#Bedtools merge here is a convenient heuristic since annotation features can overlap but are bound to a specific gene such as is the case for UTRs
#Genes can also overlap and it is probably fine to consider them as one long gene region instead of counting them as a double overlap.
awk '{print $1 "\t" $2 "\t" $3}' $annot | sort -k1,1 -k2,2n | bedtools merge > tmp_annot_sp
annot=tmp_annot_sp




bedtools intersect -a tmp.sites -b $annot -wa > tmp.sites.bed
bedtools intersect -a tmp.sites.bed -b $mRNA_list -wb > tmp_sites_mrna.bed

cut -f12 tmp_sites_mrna.bed |  sed -E 's/ID.*;Parent=(.*);Name.*;/\1/g' > GO_analysis/$name.$i.$shortAnnot.geneIDs.list
cut -f1,3 tmp_sites_mrna.bed > GO_analysis/$name.$i.$shortAnnot.loci.list



done

rm tmp_sites_mrna.bed
rm tmp_annot_sp
rm tmp.sites
rm tmp.sites.bed
