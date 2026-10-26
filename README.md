
[![DOI](https://zenodo.org/badge/1048927647.svg)](https://doi.org/10.5281/zenodo.23561348)

# Butterfly_DNAm_age

This repository contains code related to the publication "An epigenetic aging clock in butterflies".

DOI: 10.1098/rspb.2026.1679

A typical workflow is:

1. Process methylation sequencing data using scripts in **methylseq_pipe/**.
2. Construct DNA methylation clocks using scripts in **analysis/**.
3. Perform DMR analyses using scripts in **DMRs/** if required.
4. Generate figures and statistical summaries using **plots_and_stats/**.
5. Run isotope-based geolocation analyses from **isotope_geolocation/** if required.

![alt text](https://github.com/JesperBoman/Butterfly_DNAm_age/blob/main/Butterfly_DNAm_age_github2.jpg)



## Repository structure and analytical workflow

```mermaid
flowchart LR

    A[Raw bisulfite sequencing data]

    B[methylseq_pipe/
    Read processing and
    methylation calling]

    C[analysis/
    Epigenetic clock
    construction]

    D[DMRs/
    Differentially methylated
    region analyses]

    F[plots_and_stats/
    Manuscript figures and
    statistical summaries]

    G[isotope_geolocation/
    Geographic origin analyses]

    H[Final manuscript results]

    A --> B
    B --> C
    B --> D
    D --> C
    C --> F
    G --> H
    F --> H
```


**License**

This project is distributed under the terms of the license provided in the LICENSE file.
