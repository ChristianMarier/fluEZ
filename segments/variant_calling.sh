#!/bin/bash

# Script to perform variant calling for the FluEZ pipeline

# Usage check
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <project_directory> <sample_name> <threads>"
    exit 1
fi

# Assign variables from arguments
PROJECT_DIR=$1
SAMPLE_NAME=$2
THREADS=$3


#########################


# Establish paths
DEDUP_BAM="$PROJECT_DIR/bam/${SAMPLE_NAME}_dedup.bam"
REF_FASTA=$(grep "^REF_FASTA" "$PROJECT_DIR/settings.txt" | cut -d '=' -f2)

# Check for the existence of the project directory and input BAM file
if [ ! -d "$PROJECT_DIR" ]; then
    echo "Error: Project directory does not exist - $PROJECT_DIR"
    exit 1
fi

if [ ! -f "$DEDUP_BAM" ]; then
    echo "Error: Deduplicated BAM file does not exist - $DEDUP_BAM"
    exit 1
fi

echo "Input validation complete. Proceeding with variant calling for sample: $SAMPLE_NAME using $THREADS threads."
echo "Deduplicated BAM: $DEDUP_BAM"

# Ensure the module environment is clean
module purge

# Load the required modules
module load bcftools/1.15.1
module load tabix/0.2.6
module load python/cpu/3.6.5

# Confirm modules are loaded
echo "Loaded modules:"
module list

# Directory structure setup for outputs and logs
VARIANT_DIR="$PROJECT_DIR/variants"
LOG_DIR="$PROJECT_DIR/logs-vc"

# Create directories if they don't exist
mkdir -p "$VARIANT_DIR"
mkdir -p "$LOG_DIR"

echo "Directory structure set up at:"
echo "Variants: $VARIANT_DIR"
echo "Logs: $LOG_DIR"

# Setup paths for files
MPILEUP_OUT="$VARIANT_DIR/${SAMPLE_NAME}_mpileup.vcf"
VCF_OUT="$VARIANT_DIR/${SAMPLE_NAME}_raw.vcf"
VCF_NORM="$VARIANT_DIR/${SAMPLE_NAME}_norm.vcf"
FILTERED_VCF="$VARIANT_DIR/${SAMPLE_NAME}_filtered.vcf"
STATS="$LOG_DIR/${SAMPLE_NAME}_variant_stats.txt"


#########################


# Variant calling with bcftools
echo "Starting variant calling with bcftools..."
bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --redo-BAQ -d 10000 -Q 20 -C 50 -Ov -o "$MPILEUP_OUT" -f "$REF_FASTA" "$DEDUP_BAM" --threads $THREADS
bcftools call --ploidy 1 -A -mv -Ov --threads $THREADS -o "$VCF_OUT" "$MPILEUP_OUT"

# Check for successful variant calling
if [ ! -f "$VCF_OUT" ]; then
    echo "Variant calling failed to produce a VCF file."
    exit 1
fi

# Apply filters to remove artifacts and improve the quality of variant calls
echo "Applying filters to clean up the VCF..."
bcftools norm -m -any "$VCF_OUT" > "$VCF_NORM"
bcftools filter -i '(IMF > 0.5 || (DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)) & FORMAT/AD[0:1] > 50' -Ov -o "$FILTERED_VCF" "$VCF_NORM"

# Check for successful filtering
if [ ! -f "$FILTERED_VCF" ]; then
    echo "Filtering failed to produce a cleaned VCF file."
    exit 1
fi

# Validation of filtered variants
bcftools stats "$FILTERED_VCF" > "$STATS"

echo "Variant calling and filtering completed successfully."
echo "Raw VCF stored at: $VCF_OUT"
echo "Filtered VCF stored at: $FILTERED_VCF"

# Compress filtered variants for generating consensus
echo "Compressing and indexing the VCF file..."
bgzip -c "$FILTERED_VCF" > "${FILTERED_VCF}.gz"
tabix -p vcf "${FILTERED_VCF}.gz"

# Confirm compression and indexing
if [ -f "${FILTERED_VCF}.gz" ] && [ -f "${FILTERED_VCF}.gz.tbi" ]; then
    echo "VCF file compressed and indexed successfully."
else
    echo "Error compressing or indexing VCF file."
    exit 1
fi


#########################


# Low-frequency variant detection
echo "Detecting low-frequency variants..."

# Extract CHROM, POS, REF, ALT, and FORMAT columns from the mpileup VCF file and convert to TXT
echo -e "CHROM\tPOS\tREF\tALT\tFORMAT\tVALUES\tLENGTH" > "$VARIANT_DIR/${SAMPLE_NAME}_mpileup.txt"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FORMAT\n' "$MPILEUP_OUT" >> "$VARIANT_DIR/${SAMPLE_NAME}_mpileup.txt"

# Add LENGTH column with reference sequence lengths
echo "Adding LENGTH column to the TXT..."
awk 'BEGIN {OFS="\t"}
     NR==FNR {if ($1 ~ /^##contig/) {match($1, /ID=([^,]+),length=([0-9]+)/, a); len[a[1]] = a[2]}}
     NR!=FNR {print $0, len[$1]}' "$MPILEUP_OUT" "$VARIANT_DIR/${SAMPLE_NAME}_mpileup.txt" > "$VARIANT_DIR/${SAMPLE_NAME}_mpileup_lengths.txt"

# Run custom Python script to detect low-frequency variants
python3 "$PROJECT_DIR/fluEZ/scripts/detect-variants.py" --input_mp "$VARIANT_DIR/${SAMPLE_NAME}_mpileup_lengths.txt" --output_hf "$VARIANT_DIR/${SAMPLE_NAME}_high_freq_variants.txt" --output_lf "$VARIANT_DIR/${SAMPLE_NAME}_low_freq_variants.txt"

# Check if low-frequency variant detection was successful
if [ $? -ne 0 ]; then
    echo "Low-frequency variant detection failed."
    exit 1
fi

echo "Low-frequency variant detection completed successfully."
echo "Low-frequency variants stored at: $VARIANT_DIR/${SAMPLE_NAME}_low_freq_variants.txt"


#########################

# end


