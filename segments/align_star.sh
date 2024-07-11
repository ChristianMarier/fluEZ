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
module load star/2.7.3a

# Prepare directories
BAM_DIR="$PROJECT_DIR/bam"
LOG_DIR="$PROJECT_DIR/logs-align"
mkdir -p $BAM_DIR $LOG_DIR

# File paths
R1="$PROJECT_DIR/trimmed_fastq/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
R2="$PROJECT_DIR/trimmed_fastq/${SAMPLE_NAME}_trimmed_R2.fastq.gz"
TMP_BAM="$BAM_DIR/${SAMPLE_NAME}.tmp.bam"
FINAL_BAM="$BAM_DIR/${SAMPLE_NAME}.bam"
BAM_LOG="$LOG_DIR/${SAMPLE_NAME}_star.log"

# Skip if final BAM file already exists and is not corrupt
if [ -e "$FINAL_BAM" ] && [ -s "$FINAL_BAM" ]; then
    echo "BAM file already exists and is complete: $FINAL_BAM"
    exit 0
fi

# Remove potentially corrupt BAM files
rm -f "$TMP_BAM"

# Alignment command construction
CMD="
STAR --runThreadN $THREADS --genomeDir ${REF_GENOME}STAR --genomeLoad NoSharedMemory --outFilterMismatchNoverLmax 0.2 --outFilterMultimapNmax 1 \
--outFilterType BySJout --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM MD XS \
--outSAMattrRGline ID:${SAMPLE_NAME} SM:${SAMPLE_NAME} LB:${SAMPLE_NAME} PL:ILLUMINA --outSAMmapqUnique 60 --twopassMode Basic --readFilesCommand zcat \
--readFilesIn $R1 $R2 --outFileNamePrefix ${LOG_DIR}/${SAMPLE_NAME}_ --quantMode GeneCounts --outSAMtype BAM Unsorted --outStd BAM_Unsorted | \
samtools sort -m 32G -T ${SAMPLE_NAME}.samtools -o $FINAL_BAM -
"

# Execute alignment and check for errors
echo "Running STAR alignment..."
eval "$CMD"

# Verify BAM file creation
if [ ! -s "$FINAL_BAM" ]; then
    echo "Failed to create BAM file: $FINAL_BAM"
    exit 1
fi


#########################


# BAM post-processing with Samtools
echo "Indexing BAM file..."
CMD="samtools index $FINAL_BAM"
eval "$CMD"

# Calculate alignment statistics
STATS_LOG="$LOG_DIR/${SAMPLE_NAME}_stats.log"
if ! samtools flagstat $FINAL_BAM > $STATS_LOG; then
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
