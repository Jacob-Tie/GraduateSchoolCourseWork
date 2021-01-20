#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --time=04:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --job-name=plinkToZarr
#SBATCH --output=plinkToZarr.%j.out



source PATH_TO_ANACONDA
conda PATH_TO_YOUR_ENVIRONMENT

python3 plinkToZarr.py