#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32768
#SBATCH --partition=week
#SBATCH --time=150:00:00
#SBATCH --mail-type=all

module load Trimmomatic

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
echo "$i"
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 2 -phred33 -trimlog $i\.trim.log $i\_R1.fastq.gz $i\_R2.fastq.gz $i\_1P.fq.gz $i\_1U.fq.gz $i\_2P.fq.gz $i\_2U.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:125
rm *trim.log
done