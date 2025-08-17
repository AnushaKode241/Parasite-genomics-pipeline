#!/bin/bash
#SBATCH --partition=week
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL

module load SAMtools/1.15.1-GCC-11.2.0
module load BWA/0.7.17-GCCcore-11.2.0
module load BEDTools/2.30.0-GCC-11.2.0
module load GATK/4.2.5.0-GCCcore-11.2.0-Java-11

REF=/dir/homo_Ov.fa
dirA=/dir/OnchoMF_WGA_Seq_run13_4_23

CN="AGRF"
LIB="IlluminaDNA"
hprefix="homo_Ov"
prefix="nOv4.mtOv.wOv"

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
echo "mapping $i to $REF"
bwa mem $REF -R "@RG\tID:$i\tCN:$CN\tLB:IlluminaDNA\tSM:$i\tPL:ILLUMINA" $dirA/$i\_1P.fq.gz $dirA/$i\_2P.fq.gz | samtools view -bh - | samtools sort -T $i -> $i\.$hprefix\.sorted.bam

echo "$i host and parasite flagstat"
samtools flagstat $i\.$hprefix\.sorted.bam

echo "$i counting host contamination"
samtools view $i\.$hprefix\.sorted.bam | grep -c "NC_0"

echo "$i removing host contamination"
samtools view -h $i\.$hprefix\.sorted.bam | grep -v "NC_0" | samtools view -bh - > $i\.$prefix\.sorted.bam

echo "mapping stats $i parasite only sorted"
samtools flagstat $i\.$prefix\.sorted.bam

echo "removing duplicates $i"
gatk MarkDuplicatesWithMateCigar I=$i\.$prefix\.sorted.bam O=$i\.$prefix\.dedup.bam M=$i\.$prefix\.metrics REMOVE_DUPLICATES=true MINIMUM_DISTANCE=500
samtools view -bh -q 30 -F 2048 $i\.$prefix\.dedup.bam | samtools view -bh -F 256 - > $i\.$prefix\.unique.bam

echo "$i dedup"
samtools flagstat $i\.$prefix\.dedup.bam
echo "$i unique"
samtools flagstat $i\.$prefix\.unique.bam

echo "$i coverage"
bedtools genomecov -ibam $i\.$prefix\.unique.bam -g $REF -d | perl calc_cov_per_chrom.pl

done