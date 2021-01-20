#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:00:00
#SBATCH --job-name=RemoveFiles
#SBATCH --output=./slurm_out/RemoveFiles.%A-%a.out
#SBATCH --array=1-22

rm OUTPUT_BED.bed
rm OUTPUT_BIM.bim