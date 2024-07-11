#!/bin/bash


# Verify the correct number of arguments
if [ "$#" -ne 3 ]; then
	echo "Usage: $0 <project_directory> <sample_name> <threads>"
	exit 1
fi

# Assign input arguments
PROJECT_DIR=$1
SAMPLE_NAME=$2
THREADS=$3


#########################


# Setup variables for file paths
CLEAN_DIR="$PROJECT_DIR/clean_fastq"
TRIMMED_DIR="$PROJECT_DIR/trimmed_fastq"
LOGS_DIR="$PROJECT_DIR/logs-trim"

# Ensure input FASTQ files exist
R1="$CLEAN_DIR/${SAMPLE_NAME}_R1.fastq.gz"
R2="$CLEAN_DIR/${SAMPLE_NAME}_R2.fastq.gz"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
	echo "Error: Cleaned FASTQ files do not exist for $SAMPLE_NAME."
	exit 1
fi

# Create directories for outputs and logs
mkdir -p "$TRIMMED_DIR"
mkdir -p "$LOGS_DIR"

# Set up paths for output files
OUTPUT_R1="$TRIMMED_DIR/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
OUTPUT_R2="$TRIMMED_DIR/${SAMPLE_NAME}_trimmed_R2.fastq.gz"
TRIM_LOG="$LOGS_DIR/${SAMPLE_NAME}_trimming.log"


#########################


# Load Trimmomatic module
module purge
module add default-environment
module load trimmomatic/0.36

TRIMMOMATIC_JAR="${TRIMMOMATIC_ROOT}/trimmomatic-0.36.jar"

# Construct Trimmomatic command
CMD="java -jar $TRIMMOMATIC_JAR PE -threads $THREADS \
	 $R1 $R2 \
	 $OUTPUT_R1 $TRIMMED_DIR/${SAMPLE_NAME}_unpaired_R1.fastq.gz \
	 $OUTPUT_R2 $TRIMMED_DIR/${SAMPLE_NAME}_unpaired_R2.fastq.gz \
	 ILLUMINACLIP:/gpfs/data/sequence/references/contaminants/trimmomatic.fa:2:30:10:1:true SLIDINGWINDOW:4:20 MINLEN:36"

# Run Trimmomatic and redirect output to log
echo "Running Trimmomatic for $SAMPLE_NAME..."
$CMD &> $TRIM_LOG


#########################


# Check for successful trimming and validate output files
if [ ! -f "$OUTPUT_R1" ] || [ ! -f "$OUTPUT_R2" ]; then
	echo "Error: Trimming failed for $SAMPLE_NAME. Check $TRIM_LOG for details."
	exit 1
fi


#########################


# Log summary information
echo "Trimming completed for $SAMPLE_NAME."
echo "Log file: $TRIM_LOG"

# Extract and log trimming efficiency
echo "Generating summary for $SAMPLE_NAME..."
reads_before=$(zcat $R1 | echo $((`wc -l`/4)))
reads_after=$(zcat $OUTPUT_R1 | echo $((`wc -l`/4)))
echo "Reads before trimming: $reads_before"
echo "Reads after trimming: $reads_after"
percentage=$(echo "scale=2; $reads_after/$reads_before*100" | bc)
echo "Trimming efficiency: $percentage%" >> $LOGS_DIR/${SAMPLE_NAME}_trim_summary.csv

# Update summary CSV file
echo "$SAMPLE_NAME,$reads_before,$reads_after,$percentage" >> $PROJECT_DIR/summary/trim_summary.csv

echo "All trimming processes complete for $SAMPLE_NAME."


#########################



# end
