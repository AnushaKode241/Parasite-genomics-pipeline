#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=65536
#SBATCH --partition=week
#SBATCH --time=7-00:00:00

#distribution of coverage -- how far apart are contiguous regions of coverage? how much does this vary across the genome?


#Step1
#Generate separate output files for each kit
module load SAMtools

ID="DRC_1118_4_1_dop
DRC_1118_4_2_dop
DRC_1118_4_3_dop
DRC_1118_10_1_dop
DRC_1118_10_2_dop
DRC_1118_10_3_dop
DRC_1118_3_1_gphi
DRC_1118_3_2_gphi
DRC_1118_3_3_gphi
DRC_1118_9_1_gphi
DRC_1118_9_2_gphi
DRC_1118_9_3_gphi
DRC_1118_2_1_r2g
DRC_1118_2_2_r2g
DRC_1118_2_3_r2g
DRC_1118_8_1_r2g
DRC_1118_8_2_r2g
DRC_1118_8_3_r2g
DRC_1118_5_1_pico
DRC_1118_5_2_pico
DRC_1118_5_3_pico
DRC_1118_11_1_pico
DRC_1118_11_2_pico
DRC_1118_11_3_pico
DRC_1118_1_1_res
DRC_1118_1_2_res
DRC_1118_1_3_res
DRC_1118_7_1_res
DRC_1118_7_2_res
DRC_1118_7_3_res"

for i in $ID
do

samtools mpileup DRC_1118_4_1_dop.nOv4.mtOv.wOv.unique.bam DRC_1118_4_2_dop.nOv4.mtOv.wOv.unique.bam DRC_1118_4_3_dop.nOv4.mtOv.wOv.unique.bam DRC_1118_10_1_dop.nOv4.mtOv.wOv.unique.bam DRC_1118_10_2_dop.nOv4.mtOv.wOv.unique.bam DRC_1118_10_3_dop.nOv4.mtOv.wOv.unique.bam > output_filename.pileup

#Step2:Convert the output file into a txt file containing only start and end chromosome positions(bp)
cat output_filename.pileup | perl filter_mpileup.pl | awk '{if ($3>0) print $_}'| perl assess_ranges.pl > newfilename.txt

#See below on how to obtain the length (contiguity) and distance (how far apart each contiguous region is):
#For eg: 
#start 	end  endchrom-startchrom        length   distance
36    	50   50-36 = length      	14       4
54   	103  54-50 = distance  		49       0
