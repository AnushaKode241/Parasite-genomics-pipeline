#!/bin/bash
#PBS -l ncpus=48
#PBS -l mem=40GB 
#PBS -l jobfs=40GB
#PBS -q expresssr
#PBS -P rl16
#PBS -l walltime=05:00:00
#PBS -l storage=gdata/rl16+scratch/rl16 
#PBS -l wd
#PBS -M N.Sirwani@latrobe.edu.au
#PBS -m abe 
#PBS -N kraken_big_samples


# --------------------------
# 1. Setup
# --------------------------

source /home/572/ns2132/miniconda3/etc/profile.d/conda.sh
conda activate microbiome

# Adjust this path to your FASTQ files

RAW_DIR=/g/data/rl16/oncho_mf


# Adjust to where you want outputs written
WORKDIR=/g/data/rl16/oncho_mf/kraken2

mkdir -p $WORKDIR

cd $WORKDIR

# Kraken2 database (untarred version!)
KRAKEN_DB=/g/data/rl16/nsirwani/Simulium_microbiome/k2_pluspf_16_GB

mkdir -p classified unclassified reports bracken logs

# --------------------------
# 2. Loop through samples
# --------------------------

# List of large samples to re-run
SAMPLES=(
  DRC_1118_10_1_dop
  DRC_1118_10_2_dop
  DRC_1118_10_3_dop
  DRC_1118_11_1_pico
  DRC_1118_11_2_pico
  DRC_1118_11_3_pico
  DRC_1118_1_1_res
  DRC_1118_1_2_res
  DRC_1118_1_3_res
  DRC_1118_2_1_r2g
  DRC_1118_2_2_r2g
  DRC_1118_2_3_r2g
  DRC_1118_3_1_gphi
  DRC_1118_3_2_gphi
  DRC_1118_3_3_gphi
  DRC_1118_4_1_dop
  DRC_1118_4_2_dop
  DRC_1118_4_3_dop
  DRC_1118_5_1_pico
  DRC_1118_5_2_pico
  DRC_1118_5_3_pico
  DRC_1118_7_1_res
  DRC_1118_7_2_res
  DRC_1118_7_3_res
  DRC_1118_8_1_r2g
  DRC_1118_8_2_r2g
  DRC_1118_8_3_r2g
  DRC_1118_9_1_gphi
  DRC_1118_9_2_gphi
  DRC_1118_9_3_gphi
)

for sample in "${SAMPLES[@]}"; do
  echo "=== Processing $sample ===" | tee logs/${sample}.log

  R1P=${RAW_DIR}/${sample}_1P.fq.gz
  R2P=${RAW_DIR}/${sample}_2P.fq.gz
  R1U=${RAW_DIR}/${sample}_1U.fq.gz
  R2U=${RAW_DIR}/${sample}_2U.fq.gz

  kraken2 \
    --db $KRAKEN_DB \
    --threads 16 \
    --paired $R1P $R2P \
    --unpaired $R1U $R2U \
    --report reports/${sample}.kraken2.report \
    --classified-out classified/${sample}.classified#.fq \
    --unclassified-out unclassified/${sample}.unclassified#.fq \
    > logs/${sample}.kraken2.output 2>&1

 # ---- BRACKEN abundance estimation ----
    # Adjust read length depending on your ONT/Porechop data if needed
    READLEN=100

    bracken \
        -d $KRAKEN_DB \
        -i reports/${sample}.kraken2.report \
        -o bracken/${sample}.bracken.txt \
        -l S \
        -r $READLEN \
        >> logs/${sample}.bracken.log 2>&1

done

echo "All samples completed."

