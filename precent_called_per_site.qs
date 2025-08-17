#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=65536
#SBATCH --partition=week
#SBATCH --time=7-00:00:00

#--> graph % called per site per type [=1-F]: line graph, each line = kit 

#To get this, 

vcftools --gzvcf file --keep filename.list --missing-site --out filename 
#(% called per site)

#The output file will have, 
#Chrom chrom_pos NMISS FMISS Then do 1-F_MISS to obtain % called per site information. The plot for this is generated in JMP software.
