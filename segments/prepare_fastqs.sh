#!/bin/bash


# Set up basic variables
SEGMENT="FASTQ Preparation"
echo "Starting $SEGMENT segment..."

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
	echo "Error: Incorrect number of arguments."
	echo "Usage: $0 <project_directory> <sample_name>"
	exit 1
fi

PROJECT_DIR=$1
SAMPLE_NAME=$2

# Confirm the existence of the project directory
if [ ! -d "$PROJECT_DIR" ]; then
	echo "Error: Project directory does not exist: $PROJECT_DIR"
	exit 1
fi

# Define directory for cleaned FASTQ files and summary output
CLEAN_DIR="$PROJECT_DIR/clean_fastq"
SUMMARY_DIR="$PROJECT_DIR/summary"
mkdir -p "$CLEAN_DIR" "$SUMMARY_DIR"

# Paths for cleaned FASTQ files
CLEAN_R1="$CLEAN_DIR/${SAMPLE_NAME}_R1.fastq.gz"
CLEAN_R2="$CLEAN_DIR/${SAMPLE_NAME}_R2.fastq.gz"

# Check if clean files already exist
if [ -e "$CLEAN_R1" ] && [ -e "$CLEAN_R2" ]; then
	echo "Cleaned FASTQ files already exist for $SAMPLE_NAME, skipping..."
	exit 0
fi

# Find raw FASTQ files from CSV
RAW_CSV="$PROJECT_DIR/samples_fastq.csv"
if [ ! -f "$RAW_CSV" ]; then
	echo "Error: CSV file listing raw samples does not exist: $RAW_CSV"
	exit 1
fi

# Read CSV and process FASTQ files
while IFS=',' read -r sample r1 r2; do
	if [ "$sample" == "$SAMPLE_NAME" ]; then
		#Validate FASTQ files
		if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
			echo "Error: FASTQ files for $sample are missing or corrupt."
			exit 1
		fi

		# Concatenate multiple R1 and R2 files or link single files
		if grep -q "$SAMPLE_NAME" "$RAW_CSV"; then
			cat $(grep "$SAMPLE_NAME" "$RAW_CSV" | cut -d ',' -f2) > "$CLEAN_R1"
			cat $(grep "$SAMPLE_NAME" "$RAW_CSV" | cut -d ',' -f3) > "$CLEAN_R2"
		else
			ln -s "$r1" "$CLEAN_R1"
			ln -s "$r2" "$CLEAN_R2"
		fi

		# Calculate read counts
		reads_r1=$(zcat "$CLEAN_R1" | echo $((`wc -l`/4)))
		reads_r2=$(zcat "$CLEAN_R2" | echo $((`wc -l`/4)))

		# Validate read counts and output to summary CSV
		if [ "$reads_r1" -ne "$reads_r2" ] || [ "$reads_r1" -lt 100000 ]; then
			echo "Discrepancy or insufficient reads in $sample: R1=$reads_r1, R2=$reads_r2" >> "$SUMMARY_DIR/error_log.csv"
			rm "$CLEAN_R1" "$CLEAN_R2"  # Remove problematic files
		else
			echo "$sample,$reads_r1,$reads_r2" >> "$SUMMARY_DIR/fastq_summary.csv"
		fi
	fi
done < "$RAW_CSV"

echo "All FASTQ files processed and summaries created."



# end
