| File              | Description                                                      |
| ----------------- | ---------------------------------------------------------------- |
| `run_pipe.sh`     | Main methylation-processing pipeline.                            |
| `cat_fastq.sh`    | FASTQ file preparation and concatenation.                        |
| `samples.list`    | Sample definition file for pipeline execution.                   |
| `lambda_phage.fa` | Lambda phage reference used for bisulfite conversion assessment. |

This is a directory for the [Nextflow methylseq pipeline](https://github.com/nf-core/methylseq), which is the pipeline I use to go from sequence reads (fastq-files) to methylation calls.

Step 1: Use cat_fastq.sh to combine technical replicate fastq files (sequenced on different lanes, for example) of the same biological sample.

Step 2: Run pipeline with run_pipe.sh
