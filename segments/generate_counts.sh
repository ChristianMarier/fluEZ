#!/bin/bash

# Script to quantify gene expression from BAM files using featureCounts

## Ensure the script is called with the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <project_directory> <sample_name> <number_of_threads>"
    exit 1
fi

## Variables setup
PROJECT_DIR=$1
SAMPLE_NAME=$2
THREADS=$3


#######################


echo "Starting gene expression quantification for sample: $SAMPLE_NAME"
echo "Project Directory: $PROJECT_DIR"
echo "Sample Name: $SAMPLE_NAME"
echo "Number of Threads: $THREADS"

# Input BAM file path
BAM_FILE="$PROJECT_DIR/bam/${SAMPLE_NAME}.bam"

# Check for the existence of the project directory and input BAM file
if [ ! -d "$PROJECT_DIR" ]; then
    echo "Error: Project directory does not exist - $PROJECT_DIR"
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file does not exist - $BAM_FILE"
    exit 1
fi

# Check for the existence of the GTF file
REF_GENOME=$(grep "^GENOME_DIR=" "$PROJECT_DIR/settings.txt" | cut -d '=' -f2)
GTF_FILE="${REF_GENOME}genes.gtf"
if [ ! -f "$GTF_FILE" ]; then
    echo "Error: GTF file does not exist - $GTF_FILE"
    exit 1
fi

echo "Input validation complete. Proceeding with featureCounts quantification."


#######################


# Ensure the module environment is clean
module purge

# Load the required module for featureCounts
module load subread/1.6.3

# Confirm modules are loaded
echo "Loaded modules:"
module list

# Directory structure setup for outputs and logs
COUNTS_DIR="$PROJECT_DIR/counts"
#SUMMARY_DIR="$COUNTS_DIR/summary"
#LOG_DIR="$COUNTS_DIR/logs-counts"

# Create directories if they don't exist
mkdir -p "$COUNTS_DIR"
#mkdir -p "$SUMMARY_DIR"
#mkdir -p "$LOG_DIR"

echo "Directory structure set up at:"
echo "Counts: $COUNTS_DIR"
#echo "Summary: $SUMMARY_DIR"
#echo "Logs: $LOG_DIR"

# Setup paths for files
RAW_COUNTS="$COUNTS_DIR/${SAMPLE_NAME}_raw_counts.txt"
CLEAN_COUNTS="$COUNTS_DIR/${SAMPLE_NAME}_clean_counts.txt"
#SUMMARY_FILE="$SUMMARY_DIR/${SAMPLE_NAME}_counts_summary.txt"

# Run Type and Strand Flags
RUN_TYPE="-p"
STRAND_FLAG="-s 2"

# Construct featureCounts command
echo "Constructing featureCounts command..."
FC_CMD="featureCounts -T $THREADS -a $GTF_FILE -o $RAW_COUNTS $RUN_TYPE -g gene_name $STRAND_FLAG $BAM_FILE"

# Execute featureCounts command
echo "Running featureCounts..."
echo $FC_CMD
eval $FC_CMD

# Check for successful featureCounts execution
if [ ! -f "$RAW_COUNTS" ]; then
    echo "featureCounts failed to produce a raw counts file."
    exit 1
fi


#######################


# Extract gene counts and generate clean counts file
echo "Cleaning raw counts data..."
awk 'BEGIN {OFS="\t"} NR>2 {print $1, $7}' "$RAW_COUNTS" > "$CLEAN_COUNTS"


#######################


echo "Gene expression quantification complete for sample: $SAMPLE_NAME."
echo "Date and Time: $(date)"


# end



