#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=4:00:00
#SBATCH --job-name=magmaCommands
#SBATCH --output=magmaCommands.%j.out

for chr in {1..23}; do \
PATH_TO_MAGMA --annotate --snp-loc PATH_TO_BIM.bim --gene-loc PATH_TO_MAGMA_GENE_LOC.gene.loc --out MAGMA_OUTPUT_PATH_${chr}
done
