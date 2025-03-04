# Nuclei Isolation Initial Submission

Code used to analyze data and generate figures for Kersey, Acri, Dabin et al., _Comparative analysis of nuclei isolation protocols for transcriptomic profiling of brain tissue_ (DOI: XXXXXXXX)

We performed {brief summary of methods and overview}

Raw data (fastq and cellranger output) can be found GEO (GSE:XXXXXXXX). Private until published, please use reviewer key to access before publication.

The code in this repository is everything which can analyze the data as well as generate all the figures for the paper, both main and supplementary. A full list of figures and where you can find the code for them is below: 

Folder | Script | Figures
--- | --- | ---
scripts/wet-lab | NuleiRecoveryPlots | 1E, 1F
scripts/preprocess | Step1_IsolationPreProcess_DFandSX | NA
scripts/preprocess | Step2_job_IntegrationPostSXandDFv2 | NA
scripts/preprocess | Step3_ Clustering2DetermineResolution | NA
scripts/analysis | HOLLY_ANNOTATION_AND_PROPS | 2A, 2B, 2C, 2D
scripts/analysis | HOLLY_QC | 3A-#
scripts/analysis | SAHA | 4A
scripts/analysis | CellIdentity | 4B-N