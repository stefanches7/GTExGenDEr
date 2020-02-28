#!/bin/bash
# 
# author: Stefan Dvoretskii
#
# run Snakemake on SLURM 
#
####
# SLURM OPTION
####
# number of CPU/Threads
#SBATCH -n 10
#SBATCH -N 1
#
# memory limit in Mb
#SBATCH --mem 120000
#
# time need to run the script [HH:MM:SS]
#SBATCH --time 10:00:00
#
# STDERR/STDOUT redirection %J == jobid from slurm
# SBATCH --error ./output/job.%J.err
# SBATCH --output ./output/job.%J.out
#
# email report
#SBATCH --mail-type ALL
#SBATCH --mail-user stefan.dvoretskii@campus.lmu.de
#
# set the working directory
#SBATCH -D /data/ouga/home/ag_gagneur/dvoretsk/projects/GTExGenDEr
#
########
# END OF SLURM CONFIG
########
#
# run program
#

echo "Slurm job started"
snakemake --unlock
snakemake --cores 6 --reason
echo "Slurm job ended"
