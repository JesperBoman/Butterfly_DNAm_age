This is a pipeline to call differentially methylated regions (DMRs).

It uses the R package [bsseq](https://www.bioconductor.org/packages/devel/bioc/vignettes/bsseq/inst/doc/bsseq.html) and the [BSmooth method](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r83).

1. First we reformat coverage2cytosine reports from Bismark to bsseq format using DMR_prep.sh
2. Then we read the data into R and perform smoothing of methylation values using DMR_analysis_step1.R.
3. Then we call DMRs using: DMR_analysis_findDMRs.R.

Step 2 and 3 can be run on the cluster using the DMR_analysis_top_epiClock.sh run script.
