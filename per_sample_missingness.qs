#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=65536
#SBATCH --partition=week
#SBATCH --time=7-00:00:00

#Load module
module load VCFtools

vcftools --gzvcf filename.vcf.gz --missing-indv --out missing_indv


