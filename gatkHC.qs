#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=65536
#SBATCH --partition=week
#SBATCH --time=7-00:00:00


module load SAMtools
module load GATK

REF=Ov_nOv4_mtOv_wOv.fa
prefix=nOv4.mtOv.wOv
genome="nu"
min=30
ploidy=2
dir=dir

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

samtools index $dir/$i\.$prefix\.unique.bam
gatk HaplotypeCaller -R $REF -I $dir/$i\.$prefix\.unique.bam -O $i\.$genome\.gatk.vcf.gz -ERC GVCF -ploidy $ploidy --linked-de-bruijn-graph -RF MappingQualityReadFilter --minimum-mapping-quality $min -L OVOC.OM1a_TELO_TELO -L OVOC.OM2_TELO_TELO -L OVOC.OM3_TELO_TELO -L OVOC.OM4_TELOL_LFR

done