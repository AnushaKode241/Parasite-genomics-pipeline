#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32768
#SBATCH --partition=week
#SBATCH --time=120:00

module load GATK

## reference path
REF=Ov_nOv4_mtOv_wOv.fa
fileBaseName="Nov2023.Ov.WGA.comp"
ploidy=2
chrom=nuc

gatk CombineGVCFs -O $fileBaseName.$chrom.HC.vcf.gz -R $REF --variant variant.list
gatk IndexFeatureFile -I $fileBaseName.$chrom.HC.vcf.gz -O $fileBaseName.$chrom.HC.vcf.tbi
gatk GenotypeGVCFs -R $REF -V $fileBaseName.$chrom.HC.vcf.gz -O $fileBaseName.$chrom.GG.vcf.gz -ploidy $ploidy -all-sites
