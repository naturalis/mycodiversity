# PROFUNGIS
Processing of Fungal ITS Sequences
This pipeline has been created to download SRA studies of fungal ITS markers that will be used to construct a unique set of Z(ero radiance) OTUS to analyse fungal diversity.
## Installing PROFUNGIS
Install before usage:
  - Snakemake
  - BLAST+
  - usearch v11
    - rename the executable file 'usearch11' and place in /deps/ folder
To check if all dependencies are placed correctly, execute the checkDependencies.sh script.

## The PROFUNGIS workflow
PROFUNGIS is divided into three main steps: preprocessing, denoising and postprocessing. In the preprocessing step, the reads corresponding to the provided SRA run id(s) will be downloaded into the ./sample/ directory. Next, the presence of the provided primers in the primer datafile is checked. If one or both of the primers are missing, the user is asked to provide the sequence and the part of ITS it amplifies.
Based on the collection of the user-provided parameters, a configuration file is created (config.yml). This file is required by the denoising part of the pipeline
