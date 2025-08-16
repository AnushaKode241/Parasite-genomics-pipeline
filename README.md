# Parasite-genomics-pipeline
Scripts and workflows for preprocessing NGS data and comparing whole genome amplification (WGA) methods, based on analyses from a PhD thesis focused on a parasitic nematode *Onchocerca volvulus*.

# WGA Comparison NGS Pipeline

This repository contains scripts and workflows developed during the first chapter of my PhD thesis, which focuses on comparing whole genome amplification (WGA) methods in *Onchocerca volvulus*. The project includes standard NGS preprocessing steps (trimming, mapping, variant calling) and downstream analyses specific to WGA performance evaluation.

This modular pipeline can be adapted for use with other nematode species, including *Wuchereria bancrofti*, *Brugia malayi*, and *Dirofilaria immitis*, as well as for broader research applications involving whole-genome amplification and next-generation Illumina sequencing data.

## Project Overview

The overall aim is to assess and compare different WGA methods by analyzing sequencing data generated from amplified DNA. The workflow includes:
   1. Quality trimming of Illumina FASTQ files (`trimmomatic v.0.32`)  
   2. Read mapping to the *O. volvulus* and human genomes (`bwa v.0.7.17`)  
   3. Filtering alignments (`samtools v.1.9`, `picard`, `bedtools v.2.26.0`)  
   4. Variant calling using GATK HaplotypeCaller  (`GATK v.4.2.6.1`, `VCFtools v.0.1.16`)
   5. Combining variants and genotype calling (`GATK v.4.2.6.1`, `VCFtools v.0.1.16`)
   6. Data filtering using (`VCFtools v.0.1.16`)
   7. Analysis of the filtered data using (`Rstudio v.4.3.2`) and (`JMP v.17.2.0`)

Environment setup and tool installation are not included here — users should install these tools via conda or system package managers.

#How to Use

Each script is modular. Run them step by step, editing paths and sample names as required.

Example:
```bash
bash scripts/trim.qs sample_R1.fastq.gz sample_R2.fastq.gz
bash scripts/map.qs trimmed_R1.fastq.gz trimmed_R2.fastq.gz reference.fasta
