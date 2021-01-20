#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=plinkToZarrQC
#SBATCH --output=plinkToZarrQC.%j.out



source PATH_TO_ANACONDA
conda PATH_TO_YOUR_ENVIRONMENT

python3 plinkToZarrQC.py
