#!/bin/bash
 
#PBS -P rs55
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=190gb
#PBS -l jobfs=15GB
#PBS -l walltime=20:00:00
#PBS -M clara.savary@newcastle.edu.au
#PBS -m ae
#PBS -j oe
#PBS -N PPK_COMBO_beta9
#PBS -o /scratch/rs55/cs4309/PPK_COMBO/scMultiome_pip/logs
#PBS -e /scratch/rs55/cs4309/PPK_COMBO/scMultiome_pip/stderr


#__________________________________________________________________________________
# LOAD CONDA ENVIRONMENT ----
#__________________________________________________________________________________


#PATH="$HOME/tools/miniconda3/bin:$PATH"

source activate env_scMultiome

#module load R


#__________________________________________________________________________________
# CUSTOMIZABLE PARAMETERS ----
#__________________________________________________________________________________

# Script file
RMD="/scratch/rs55/cs4309/PPK_COMBO/scMultiome_pip/code/scripts/scMultiome_preprocessing_beta_9.Rmd"

# Configuration file
CONFIG="/scratch/rs55/cs4309/PPK_COMBO/scMultiome_pip/config/scMultiome_config.tsv"

# Format date
current_date=$(date "+%Y_%m_%d")

# Output file
OUT="/scratch/rs55/cs4309/PPK_COMBO/output/scMultiome_preprocessing_beta_9_$current_date.html"


#__________________________________________________________________________________
# COMMAND LINES ----
#__________________________________________________________________________________


## LAUNCH SCMULTIOME PREPROCESSING ----
Rscript -e "rmarkdown::render(input = '$RMD', params = list(config = '$CONFIG'), output_file = '$OUT')"