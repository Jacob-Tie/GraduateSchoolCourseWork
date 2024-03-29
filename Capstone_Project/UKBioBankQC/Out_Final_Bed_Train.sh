#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=9:00:00
#SBATCH --job-name=OutputBed
#SBATCH --output=./slurm_out/OutputBed.%A-%a.out
#SBATCH --array=1-22

gunzip -c PATH_TO_ZIPPED_UKBIOBANK_DATA_SPLIT_BY_CHROM_${SLURM_ARRAY_TASK_ID}.bed.gz > OUTPUT_BED.bed
gunzip -c PATH_TO_ZIPPED_UKBIOBANK_DATA_SPLIT_BY_CHROM_${SLURM_ARRAY_TASK_ID}_v2.bim.gz > OUTPUT_BIM.bim
plink1.9 --bed OUTPUT_BED.bed\
                          --bim OUTPUT_BIM.bim\
                          --fam PATH_TO_FAM_FILE.fam\
                          --pheno PATH_TO_ANY_OTHER_TARGET_VARS_TO_INCLUDE.txt\
                          --extract PATH_TO_RELEVANT_SNPS_${SLURM_ARRAY_TASK_ID}.txt\
                          --keep PATH_TO_INDIV_IN_TRAIN_SET.txt\
                          --make-bed\
                          --out PATH_TO_BED_TRAIN_${SLURM_ARRAY_TASK_ID}
rm OUTPUT_BED.bed
rm OUTPUT_BIM.bim

