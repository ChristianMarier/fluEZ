#!/bin/bash

# Script to perform post-processing quality control in the FluEZ pipeline

# Ensure the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <project_directory> <sample_name> <threads>"
    echo "Error: Incorrect number of arguments. Expected 3, got $#."
    exit 1
fi

# Set variables from input arguments
PROJECT_DIR=$1
SAMPLE_NAME=$2
THREADS=$3


#########################


# Define path for the deduplicated BAM file
DEDUP_BAM="$PROJECT_DIR/bam/${SAMPLE_NAME}_dedup.bam"

# Check for the existence of the project directory
if [ ! -d "$PROJECT_DIR" ]; then
    echo "Error: Project directory does not exist: $PROJECT_DIR"
    exit 1
fi

# Check for the existence of the deduplicated BAM file
if [ ! -f "$DEDUP_BAM" ]; then
    echo "Error: Deduplicated BAM file does not exist: $DEDUP_BAM"
    exit 1
fi

echo "QC Input validation and setup completed successfully."
echo "Proceeding with quality control post-processing for sample: $SAMPLE_NAME"
echo "Using deduplicated BAM file located at: $DEDUP_BAM"

echo "Setting up the environment for quality control post-processing..."

# Clear any previously loaded modules to ensure a clean environment
module purge

# Load required modules
module load samtools/1.9
module load picard-tools/2.18.20
#module load python/cpu/3.6.5
#python -m pip install altair

PICARD='java -jar /gpfs/share/apps/picard-tools/2.18.20/raw/picard/build/libs/picard.jar'

echo "Loaded modules:"
module list

# Define directories for output files and logs
QC_OUTPUT_DIR="$PROJECT_DIR/post-processing_qc"
#QC_REPORT_DIR="$QC_OUTPUT_DIR/reports"
#QC_LOG_DIR="$QC_OUTPUT_DIR/logs"
#QC_VIS_DIR="$QC_OUTPUT_DIR/visualizations"

# Create directories if they do not exist
#mkdir -p $QC_REPORT_DIR $QC_LOG_DIR $QC_VIS_DIR
mkdir -p $QC_OUTPUT_DIR

echo "Directories for storing quality control outputs have been set up at: $QC_OUTPUT_DIR"
#echo "Reports will be stored in: $QC_REPORT_DIR"
#echo "Logs will be stored in: $QC_LOG_DIR"
#echo "Visualizations will be stored in: $QC_VIS_DIR"

echo "Environment setup completed successfully. Ready to proceed with quality control processes."


#########################


echo "Calculating coverage statistics using Samtools..."

# Load the reference genome and regions from settings.txt
REF_FASTA=$(grep "^REF_FASTA=" "$PROJECT_DIR/settings.txt" | cut -d '=' -f2)
GENOME_BED=$(grep "^GENOME_DIR=" "$PROJECT_DIR/settings.txt" | cut -d '=' -f2)genome.bed

# Check if reference files are available
if [[ ! -f "$REF_FASTA" || ! -f "$GENOME_BED" ]]; then
    echo "Error: Reference fasta or BED file not found. Please check the settings."
    exit 1
fi

# Define output file path for coverage statistics
COVERAGE_STATS="$QC_OUTPUT_DIR/${SAMPLE_NAME}_coverage_stats.txt"

# Using Samtools to calculate depth of coverage across specified regions
samtools depth -a -b "$GENOME_BED" "$DEDUP_BAM" > "$COVERAGE_STATS"

# Check if coverage statistics were generated
if [ ! -s "$COVERAGE_STATS" ]; then
    echo "Error: Coverage statistics could not be generated."
    exit 1
fi

echo "Coverage statistics calculated and stored at: $COVERAGE_STATS"

# Generate visualization (under development)
#echo "Generating visualization for coverage statistics..."
#python $CODE_DIR/scripts/visualize_coverage.py --input "$COVERAGE_STATS" --output "$QC_VIS_DIR/${SAMPLE_NAME}_coverage.html"

echo "Coverage statistics complete for $SAMPLE_NAME."
echo "Date and Time: $(date)"


#########################


echo "Gathering alignment statistics for $SAMPLE_NAME using Samtools..."

# Define output files for alignment statistics
FLAGSTAT_OUTPUT="$QC_OUTPUT_DIR/${SAMPLE_NAME}_flagstat.txt"
STATS_OUTPUT="$QC_OUTPUT_DIR/${SAMPLE_NAME}_stats.txt"

# Run samtools flagstat to calculate basic alignment statistics
samtools flagstat "$DEDUP_BAM" > "$FLAGSTAT_OUTPUT"

# Check if the flagstat output was generated
if [ ! -s "$FLAGSTAT_OUTPUT" ]; then
    echo "Error: Failed to generate alignment statistics with flagstat."
    exit 1
fi

echo "Basic alignment statistics generated: $FLAGSTAT_OUTPUT"

# Run samtools stats for more detailed statistics
samtools stats "$DEDUP_BAM" > "$STATS_OUTPUT"

# Check if the detailed stats output was generated
if [ ! -s "$STATS_OUTPUT" ]; then
    echo "Error: Failed to generate detailed alignment statistics with samtools stats."
    exit 1
fi

echo "Detailed alignment statistics available at: $STATS_OUTPUT"

# Optional: Generate visualizations for alignment statistics (under development)
#echo "Visualizing alignment statistics..."
#python $CODE_DIR/scripts/visualize_alignment_stats.py --input "$STATS_OUTPUT" --output "$QC_VIS_DIR/${SAMPLE_NAME}_alignment_stats.html"

echo "Alignment statistics processed for $SAMPLE_NAME."
echo "Date and Time: $(date)"


#########################


echo "Calculating insert size metrics for $SAMPLE_NAME using Picard..."

# Define output paths for insert size metrics
INSERT_METRICS_FILE="$QC_OUTPUT_DIR/${SAMPLE_NAME}_insert_metrics.txt"
INSERT_HISTOGRAM_FILE="$QC_OUTPUT_DIR/${SAMPLE_NAME}_insert_histogram.pdf"

# Run Picard CollectInsertSizeMetrics
$PICARD CollectInsertSizeMetrics \
    I="$DEDUP_BAM" \
    O="$INSERT_METRICS_FILE" \
    H="$INSERT_HISTOGRAM_FILE" \
    M=0.5

# Check if the insert size metrics were generated successfully
if [ ! -s "$INSERT_METRICS_FILE" ]; then
    echo "Error: Failed to generate insert size metrics."
    exit 1
fi

echo "Insert size metrics generated: $INSERT_METRICS_FILE"
#echo "Insert size histogram plotted: $INSERT_HISTOGRAM_FILE"

# Optional: Visualize the histogram with Python for enhanced presentation (under development)
#echo "Visualizing insert size histogram..."
#python $CODE_DIR/scripts/visualize_insert_size.py --input "$INSERT_HISTOGRAM_FILE" --output "$QC_VIS_DIR/${SAMPLE_NAME}_insert_size_plot.html"

echo "Insert size metrics processed for $SAMPLE_NAME."
echo "Date and Time: $(date)"


#########################


echo "Analyzing GC content distribution for $SAMPLE_NAME..."

# Define output paths for GC content metrics
GC_CONTENT_OUTPUT="$QC_OUTPUT_DIR/${SAMPLE_NAME}_gc_content.txt"

# Calculate GC content using samtools
samtools depth -a "$DEDUP_BAM" | \
awk '{sum+=$3} END {print "Average Depth:", sum/NR}' > "$GC_CONTENT_OUTPUT"

# Optional: Enhance GC content visualization (under development)
#echo "Visualizing GC content distribution..."
#python $CODE_DIR/scripts/visualize_gc_content.py --input "$GC_CONTENT_OUTPUT" --output "$QC_VIS_DIR/${SAMPLE_NAME}_gc_content_plot.html"

echo "GC content analysis completed for $SAMPLE_NAME."
echo "Date and Time: $(date)"


# end



