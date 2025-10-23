#!/bin/bash

#############################################
##### Preprocessing and Quality control #####
#############################################

# ----------------------------
# ------ Preprocessing -------
# ----------------------------

# Set root directory
root_dir="/path/to/the/project"


base_dir="$root_dir/01_fastqs"
output_dir="$base_dir/merged"
mkdir -p "$output_dir"

# Sample list
samples=(
    "RPE_DMSO_1_S1"
    "RPE_DMSO_2_S2"
    "RPE_DMSO_3_S3"
    "RPE_p8_1_S4"
    "RPE_p8_2_S5"
    "RPE_p8_3_S6"
    "RPE_p15_1_S7"
    "RPE_p15_2_S8"
    "RPE_p15_3_S9"
)
# Loop through each sample and merge files
for sample in "${samples[@]}"; do
    # Find and merge R1 files
    cat $(find "$base_dir" -type f -name "${sample}_L00[1-4]_R1_001.fastq.gz" | sort) > "$output_dir/${sample}_merged_R1_001.fastq.gz"

    # Find and merge R2 files
    cat $(find "$base_dir" -type f -name "${sample}_L00[1-4]_R2_001.fastq.gz" | sort) > "$output_dir/${sample}_merged_R2_001.fastq.gz"
done


# ----------------------------
# ------ Quality Control -----
# ----------------------------

# Load fastqc and multiqc module
module load FastQC/0.12.1-Java-11
module load Miniconda3/202411
multiqc --version

# Rename input and output directory
base_dir="$root_dir/01_fastqs/merged"
output_dir="$root_dir/02_qc"

# Run FastQC for all files in the directory
fastqc -o "$output_dir" -t 8 "$base_dir"/*.fastq.gz

# Summarize all reports with Multiqc
multiqc $output_dir -o $output_dir

# -----------------------------
# ------ FastP processing -----
# -----------------------------

# Load fastp
module load fastp/0.24.0

# Create output directory for fastp
fastp_output_dir="$root_dir/03_fastp"
mkdir -p "$fastp_output_dir"

# Loop through each sample and process with fastp
for sample in "${sample[@]}"; do
    fastp \
        -i "$output_dir/${sample}_merged_R1_001.fastq.gz" \
        -I "$output_dir/${sample}_merged_R2_001.fastq.gz" \
        -o "$fastp_output_dir/${sample}_filtered_R1.fastq.gz" \
        -O "$fastp_output_dir/${sample}_filtered_R2.fastq.gz" \
        -h "$fastp_output_dir/${sample}_fastp.html" \
        -j "$fastp_output_dir/${sample}_fastp.json" \
        --thread 8
done

# Summarize all reports with Multiqc
multiqc $fastp_output_dir -o $fastp_output_dir
