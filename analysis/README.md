**Varying seed clock pipeline**

Main contents
| File                                               | Description                                     |
| -------------------------------------------------- | ----------------------------------------------- |
| `master.clocksmith.sh`                             | Master workflow for epigenetic clock analyses.  |
| `fancier_varying_seed_clock5_F.R`                  | Main DNAm age model generation script.          |
| `fancier_varying_seed_clock5_F_SUBSAMPLE.R`        | Subsampling version of clock construction.      |
| `fancier_varying_seed_clock5_F_WCLOCK.R`           | Wild-clock implementation.                  |
| `fancier_varying_seed_clock5_F_WCLOCK_SUBSAMPLE.R` | Wild-clock analysis with subsampling.       |
| `smither.sh`                                       | Helper script for running clock analyses.       |
| `smither_SUBSAMPLE.sh`                             | Subsampling workflow helper.                    |
| `topGO.R`                                          | GO enrichment analyses of age-associated loci.  |
| `varying_seed_GO_prep_*.sh`                        | Preparation scripts for GO enrichment analyses. |
| `add_E2_data.sh`                                   | Data integration utility.                       |
| `comb_data.sh`                                     | Combines analysis datasets.                     |

Supporting files
| File                                 | Purpose                                    |
| ------------------------------------ | ------------------------------------------ |
| `samples_info.txt`                   | Sample metadata.                           |
| `earlyStudySamples_and_age_list.txt` | Age information for study samples.         |
| `example_input_data.txt`             | Example analysis input.                    |
| `R_sessionInfo_seed_clock.txt`       | R session information for reproducibility. |



**Data preparation**

1. Comb_data.sh
2. add_E2_data.sh

See example data structure in: example_input_data.txt



**Varying seed clock**

This pipeline randomly splits training+test data based on a particular seed and constructs a clock for each seed. We do this to evaluate the effect of a particular training+test split on predictiveness of the model. In a sense, we then get sets of clocks. 

master.clocksmith.sh runs smither.sh which in turn runs fancier_varying_seed_clock.R

It is fancier_varying_seed_clock.R that is the script that is training the penalized regression model. The master.clocksmith and smither are used to parallelize the process without sending all jobs at once, to avoid sending too many concurrent jobs to a cluster.


**GO prep**

Also included here are scripts that create files that can be used downstream for GO entology enrichment analysis.
