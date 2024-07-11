#!/bin/bash

# Usage and argument validation
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <project_directory> <sample_name> <threads>"
    exit 1
fi

# Set up environment variables
PROJECT_DIR=$1
SAMPLE_NAME=$2
THREADS=$3
SETTINGS_FILE="$PROJECT_DIR/settings.txt"


#########################


# Check existence of necessary directories and files
if [ ! -d "$PROJECT_DIR/trimmed_fastq" ]; then
    echo "Error: Trimmed FASTQ directory does not exist."
    exit 1
fi

# Load reference genome path from settings file
REF_GENOME=$(grep "^GENOME_DIR=" $SETTINGS_FILE | cut -d '=' -f2)
if [ -z "$REF_GENOME" ]; then
    echo "Reference genome path not set in settings."
    exit 1
fi

# Module environment setup
module purge
module load bowtie2/2.5.3 sambamba/0.6.8

# Prepare directories
BAM_DIR="$PROJECT_DIR/bam"
LOG_DIR="$PROJECT_DIR/logs-align"
mkdir -p $BAM_DIR $LOG_DIR

# File paths
R1="$PROJECT_DIR/trimmed_fastq/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
R2="$PROJECT_DIR/trimmed_fastq/${SAMPLE_NAME}_trimmed_R2.fastq.gz"
TMP_BAM="$BAM_DIR/${SAMPLE_NAME}.tmp.bam"
FINAL_BAM="$BAM_DIR/${SAMPLE_NAME}.bam"
BAM_LOG="$LOG_DIR/${SAMPLE_NAME}_bowtie2.log"

# Skip if final BAM file already exists and is not corrupt
if [ -e "$FINAL_BAM" ] && [ -s "$FINAL_BAM" ]; then
    echo "BAM file already exists and is complete: $FINAL_BAM"
    exit 0
fi

# Remove potentially corrupt BAM files
rm -f "$TMP_BAM"

# Alignment command construction
CMD="bowtie2 -x ${REF_GENOME}genome -1 $R1 -2 $R2 -p $THREADS --no-unal --rg-id $SAMPLE_NAME --rg LB:lib1 --rg PL:illumina --rg SM:$SAMPLE_NAME | sambamba-0.6.8 view -S -f bam -o $TMP_BAM /dev/stdin"

# Execute alignment and check for errors
echo "Running Bowtie2 alignment..."
if ! eval $CMD 2>$BAM_LOG; then
    echo "Alignment failed for $SAMPLE_NAME, see log: $BAM_LOG"
    exit 1
fi

# Verify BAM file creation
if [ ! -s "$TMP_BAM" ]; then
    echo "Failed to create BAM file: $TMP_BAM"
    exit 1
fi


#########################


# BAM post-processing with Sambamba
echo "Sorting and indexing BAM file..."
if ! sambamba-0.6.8 sort -o $FINAL_BAM --tmpdir=$BAM_DIR $TMP_BAM && sambamba-0.6.8 index $FINAL_BAM; then
    echo "Failed to sort and index BAM file: $FINAL_BAM"
    exit 1
fi

# Calculate alignment statistics
STATS_LOG="$LOG_DIR/${SAMPLE_NAME}_stats.log"
if ! sambamba-0.6.8 flagstat $FINAL_BAM > $STATS_LOG; then
    echo "Failed to compute alignment statistics."
    exit 1
fi

# Update summary file
echo "$SAMPLE_NAME,$FINAL_BAM" >> $PROJECT_DIR/summary/align_summary.csv

# Cleanup temporary files
rm -f "$TMP_BAM"

echo "Alignment and processing complete for $SAMPLE_NAME."
echo "BAM file located at: $FINAL_BAM"
echo "Alignment stats available at: $STATS_LOG"


#########################


# end
