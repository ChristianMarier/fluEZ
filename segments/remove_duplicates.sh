#!/bin/bash

# Usage and argument validation
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <project_directory> <sample_name> <threads>"
    exit 1
fi

# Variables setup
PROJECT_DIR=$1
SAMPLE_NAME=$2
THREADS=$3


#########################


# Paths for input and output files
INPUT_BAM="$PROJECT_DIR/bam/${SAMPLE_NAME}.bam"
DEDUP_BAM="$PROJECT_DIR/bam/${SAMPLE_NAME}_dedup.bam"
DEDUP_BAM_TEMP="$PROJECT_DIR/bam/${SAMPLE_NAME}_dedup_temp.bam"
LOG_DIR="$PROJECT_DIR/logs-dedup"
DEDUP_LOG="$LOG_DIR/${SAMPLE_NAME}_dedup.log"
FLAGSTAT_LOG="$LOG_DIR/${SAMPLE_NAME}_flagstat_post_dedup.log"
BAM_DIR="$PROJECT_DIR/bam"

# Check existence of necessary directories and files
if [ ! -d "$PROJECT_DIR/bam" ] || [ ! -f "$INPUT_BAM" ]; then
    echo "Error: BAM directory or input BAM file does not exist."
    exit 1
fi

# Create necessary directories
mkdir -p $LOG_DIR

# Module environment setup
module purge
module load sambamba/0.6.8

# Remove potentially corrupt deduplicated BAM files
rm -f "$DEDUP_BAM_TEMP"

# Sambamba deduplication command
echo "Removing duplicates with Sambamba..."
if ! sambamba-0.6.8 markdup -r --tmpdir=$BAM_DIR --nthreads=$THREADS "$INPUT_BAM" "$DEDUP_BAM_TEMP" 2>$DEDUP_LOG; then
    echo "Deduplication failed, see log: $DEDUP_LOG"
    rm -f "$DEDUP_BAM_TEMP"
    exit 1
fi

# Verify creation of the deduplicated BAM file
if [ ! -s "$DEDUP_BAM_TEMP" ]; then
    echo "Deduplicated BAM file was not created: $DEDUP_BAM_TEMP"
    rm -f "$DEDUP_BAM_TEMP"
    exit 1
fi


#########################


# Run Sambamba flagstat for quality control
echo "Generating post-deduplication statistics..."
if ! sambamba-0.6.8 flagstat "$DEDUP_BAM_TEMP" > $FLAGSTAT_LOG; then
    echo "Failed to generate flagstat report post-deduplication."
    exit 1
fi

# Move temporary deduplicated BAM to final location
mv "$DEDUP_BAM_TEMP" "$DEDUP_BAM"
sambamba-0.6.8 index "$DEDUP_BAM"

# Update summary file with deduplication statistics
echo "$SAMPLE_NAME,$DEDUP_BAM" >> $PROJECT_DIR/summary/dedup_summary.csv

echo "Deduplication and post-processing complete for $SAMPLE_NAME."
echo "Deduplicated BAM located at: $DEDUP_BAM"
echo "Statistics logged in: $FLAGSTAT_LOG"


#########################


# end
