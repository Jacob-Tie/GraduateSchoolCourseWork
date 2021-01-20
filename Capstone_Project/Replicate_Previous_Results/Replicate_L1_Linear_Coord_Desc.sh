#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --time=20:00:00
#SBATCH --mem=180gb
#SBATCH --nodes=2
#SBATCH --ntasks=24
#SBATCH --job-name=Replicate_L1_Linear_Coord_Desc
#SBATCH --output=Replicate_L1_Linear_Coord_Desc.%j.out

source PATH_TO_ANACONDA
conda PATH_TO_YOUR_ENVIRONMENT
python Replicate_L1_Linear_Coord_Desc.py
