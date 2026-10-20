# Butterfly_DNAm_age

This repository contains code related to the publication "An epigenetic aging clock in butterflies".

10.1098/rspb.xxxx.xxxx (to be updated)

A typical workflow is:

1. Process methylation sequencing data using scripts in **methylseq_pipe/**.
2. Construct DNA methylation clocks using scripts in **analysis/**.
3. Perform DMR analyses using scripts in **DMRs/**.
4. Generate figures and statistical summaries using **plots_and_stats/**.
5. Run isotope-based geolocation analyses from **isotope_geolocation/** if required.

## Repository workflow

```mermaid
flowchart TD
    A[methylseq_pipe<br/>Methylation sequencing processing]

    B[analysis<br/>Epigenetic clock construction<br/>and age prediction]

    C[DMRs<br/>Differentially methylated<br/>region analyses]

    D[plots_and_stats<br/>Statistical analyses<br/>and manuscript figures]

    E[isotope_geolocation<br/>Stable isotope-based<br/>geolocation analyses]

    F[Publication results<br/>Figures • Tables • Statistics]

    A --> B
    B --> C
    B --> D
    C --> D
    D --> F
    E --> F

```


![alt text](https://github.com/JesperBoman/Butterfly_DNAm_age/blob/main/Butterfly_DNAm_age_github.jpg)

**License**

This project is distributed under the terms of the license provided in the LICENSE file.
