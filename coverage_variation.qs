#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=65536
#SBATCH --partition=week
#SBATCH --time=7-00:00:00

#Check the variation in nuclear coverage between the whole-genome amplification kits used at 1x and 5x perhaps as a start
#report as a function of # reads mapped -- expect as increase # reads, will increase coverage. Better WGA method will show better curve. This is a type of rarefaction analysis

#To address this, we need to subsample the reads from the bam file and then collect the results for ("percentage of target genome with coverage") and generate the one single output file by combining all the coverage data from each sample

#NOTE: Better if each step is run separately i.e., create a .qs script file for each step

#Step1: collect number of sequences in bam file (see below code)
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

# To collect the no. of sequences in each bam file
for i in $ID

do

samtools stats $i.nOv4.mtOv.wOv.unique.bam > $i.stats

done

#Step2: subsample reads: Here I subsampled 500K to 10M reads in 500K increments for 5x coverage. Change this according to your data
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

for i in $(seq -w 5 5 100)
do
    for j in $ID
    do
echo "$j rep $i"
        f=$(grep "raw total sequences" $j.stats | cut -f3 | awk -v i=$i '{frac=i * 100000 / $1; if (frac > 1) {print "1.0"} else {print frac}}')
        mkdir -p t.subsampled_$i
        cd t.subsampled_$i
        samtools view -bs $f ../$j.nOv4.mtOv.wOv.unique.bam > $j.t.subsampled.bam && \
        samtools stats $j.t.subsampled.bam > $j.t.subsampled.stats && \
        samtools view -q 5 -b $j.t.subsampled.bam > $j.t.subsampled.q5.bam && \
        samtools stats -t ../target_regions.nuc $j.t.subsampled.q5.bam > $j.t.subsampled.stats.nuc && \
        samtools stats -t ../target_regions.mito $j.t.subsampled.q5.bam > $j.t.subsampled.stats.mito && \
        samtools stats -t ../target_regions.nuc -g 5 $j.t.subsampled.q5.bam > $j.t.subsampled.g5.stats.nuc && \
        samtools stats -t ../target_regions.mito -g 5 $j.t.subsampled.q5.bam > $j.t.subsampled.g5.stats.mito && \
        rm $j.t.subsampled.bam $j.t.subsampled.q5.bam
        cd ..
    done
done

#Step3: Collect the output/result
# Loop through the sample list
for i in $(cat sample.list); do

for subdir in t.subsampled_*; do

grep "percentage of target genome with coverage" $subdir/*.t.subsampled.g5.stats.mito | cut -f1,3 | sed s'|^/t.subsampled_||' | sed s'|/|\t|' | sed s'/.t.subsampled.g5.stats.mito:SN//' | $

done

done

# Paste temporary files into output file (run this command after obtaining the above stats)
#paste $(ls subsampled_*/tmp.*) > pct_cov_10x_nuc.txt

# Remove temporary files
#rm -f tmp.*

#Step4: For visualization generate the plots in R:
# Define data path
data_path <- "WGA_kits/pct_cov_10x_nuc.txt"

# Read data into a data frame
data <- read.table(data_path, header = TRUE)

# Define sample names and kits
sample_names <- colnames(data)[1:ncol(data)]

kit_colors <- c(rep("red", 6), rep("blue", 6), rep("green", 6), rep("black", 6), rep("orange", 6))

#kit_names = c("DOPlify", "GenomiPhi", "RTG_GenomiPhi",
              "PicoPLEX", "Resolve_DNA")

#kit_names <- colnames(data)[1:6(data)]

# Define kit colors
kit_c <- c(
  "DOPlify" = "red",
  "GenomiPhi" = "blue",
  "RTG_GenomiPhi" = "green",
  "PicoPLEX" = "black",
  "Resolve_DNA" = "orange"
)

# Extract coverage values
coverage_values <- data[, ]

par(mar = c(5, 4, 2, 10))
# Create a plot
plot(
  NULL,
  type = "l",
  xlab = "Sequencing effort (Millions)",
  ylab = "Percent Nuclear Genome Covered (10x)",
  main = "Coverage Comparison for Different WGA kits",
  xlim = c(0, 40),
  ylim = c(0, 100),
  xpd = TRUE
)

# Add lines for each sample, grouped by kit and colored accordingly
for (i in 1:length(sample_names)) {
  lines(
    x = as.numeric(rownames(coverage_values)),
    y = coverage_values[, i],
    col = kit_colors[[i]],
    lwd = 2
  )
}


# Add legend with kit names
legend(
  "right",
  inset=c(-0.45,0),
  kit_names,
  col = kit_c,
  title="WGA_kits",
  lwd = 2,
  bty = "n",
  xpd = TRUE,
  plot = TRUE
)

