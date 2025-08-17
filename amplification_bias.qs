#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32768
#SBATCH --partition=week
#SBATCH --time=120:00

#To check if the amplification is biased towards mitochondria extract what percent reads are mapped to nuclear vs that of mitochondria

# Load the necessary modules
module load SAMtools

# Get the input BAM files
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

sample_name="$i"
output_file="${sample_name}.nu.bam.txt"

#For nuclear genome readcount
samtools view -c -F 4 $i\.nOv4.mtOv.wOv.unique.bam OVOC.OM1a_TELO_TELO OVOC.OM2_TELO_TELO OVOC.OM3_TELO_TELO OVOC.OM4_TELOL_LFR >> $output_file

#For mitochondrial genome readcount
#samtools view -c -F 4 $i\.nOv4.mtOv.wOv.unique.bam OVOC_MITOCHONDRIAL >> $output_file

done

