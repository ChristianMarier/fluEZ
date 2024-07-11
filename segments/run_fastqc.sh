#!/bin/bash


# Define script info
SEGMENT_NAME="FastQC Analysis"
SCRIPT_NAME="run_fastqc.sh"
SCRIPT_PATH=$(dirname $(readlink -f $0))

echo "Script: $SCRIPT_NAME"
echo "Segment: $SEGMENT_NAME"
echo "Path: $SCRIPT_PATH"

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
	echo "Error: Incorrect number of arguments."
	echo "Usage: $0 <project_directory> <sample_name> <threads>"
	exit 1
fi

# Arguments
PROJECT_DIR=$1
SAMPLE_NAME=$2
THREADS=$3


#########################


# Paths for the cleaned FASTQ files
R1="$PROJECT_DIR/clean_fastq/${SAMPLE_NAME}_R1.fastq.gz"
R2="$PROJECT_DIR/clean_fastq/${SAMPLE_NAME}_R2.fastq.gz"

# Check if project directory exists
if [ ! -d "$PROJECT_DIR" ]; then
	echo "Error: Project directory does not exist."
	exit 1
fi

# Validate FASTQ files
if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
	echo "Error: One or both FASTQ files are missing or empty."
	exit 1
fi

# FastQC results directory
FASTQC_RESULTS_DIR="$PROJECT_DIR/FastQC_results"
mkdir -p "$FASTQC_RESULTS_DIR"

# Output file names
HTML_R1="$FASTQC_RESULTS_DIR/${SAMPLE_NAME}_R1_fastqc.html"
ZIP_R1="$FASTQC_RESULTS_DIR/${SAMPLE_NAME}_R1_fastqc.zip"
HTML_R2="$FASTQC_RESULTS_DIR/${SAMPLE_NAME}_R2_fastqc.html"
ZIP_R2="$FASTQC_RESULTS_DIR/${SAMPLE_NAME}_R2_fastqc.zip"

# Check if output already exists
if [ -e "$HTML_R1" ] && [ -e "$HTML_R2" ]; then
	echo "FastQC has already been completed for $SAMPLE_NAME."
	exit 0
fi


#########################


# Load FastQC module
module purge
module load fastqc/0.11.7

# Display the FastQC version
echo "FastQC version: $(fastqc --version)"

# Construct FastQC command
FASTQC_CMD="fastqc --threads $THREADS --outdir $FASTQC_RESULTS_DIR --extract --quiet $R1 $R2"
echo "Executing: $FASTQC_CMD"
$FASTQC_CMD


#########################


# Check results
if [ ! -s "$HTML_R1" ] || [ ! -s "$HTML_R2" ]; then
	echo "Error: FastQC did not generate expected output."
	exit 1
fi

echo "FastQC analysis completed successfully for $SAMPLE_NAME. Results are in $FASTQC_RESULTS_DIR"


#########################


# Define paths for R1 and R2 FastQC output directories
FASTQC_R1_DIR="$FASTQC_RESULTS_DIR/${SAMPLE_NAME}_R1_fastqc"
FASTQC_R2_DIR="$FASTQC_RESULTS_DIR/${SAMPLE_NAME}_R2_fastqc"

# Function to parse FastQC data for a single sample
parse_fastqc_summary() {
	local fastqc_dir=$1
	local sample_name=$2

	# Define the path to the FastQC data file
	local fastqc_data_file="${fastqc_dir}/fastqc_data.txt"

	# Check if the FastQC data file exists
	if [[ ! -f "$fastqc_data_file" ]]; then
		echo "Error: FastQC data file does not exist for $sample_name at $fastqc_data_file" >&2
		return 1
	fi

	# Extract summary statistics from FastQC data
	local total_sequences=$(grep -m 1 "^Total Sequences" "$fastqc_data_file" | cut -f 2)
	local gc_content=$(grep -m 1 "^%GC" "$fastqc_data_file" | cut -f 2)
	local sequence_length=$(grep -m 1 "^Sequence length" "$fastqc_data_file" | cut -f 2)
	local mean_q_score=$(awk '/>>Per base sequence quality/{flag=1; next} />>END_MODULE/{flag=0} flag' "$fastqc_data_file" | awk -F'\t' 'NR == 2 {print $2}')

	# Write summary to CSV
	echo "$sample_name,$total_sequences,$gc_content,$sequence_length,$mean_q_score" >> "${PROJECT_DIR}/summary/fastqc_summary.csv"
}

# Generate summary for R1 and R2
parse_fastqc_summary "$FASTQC_R1_DIR" "${SAMPLE_NAME}_R1"
parse_fastqc_summary "$FASTQC_R2_DIR" "${SAMPLE_NAME}_R2"

# Output the CSV header if the file is newly created
if [[ ! -s "${PROJECT_DIR}/summary/fastqc_summary.csv" ]]; then
	echo "Sample,Sequences Checked,GC Content,Sequence Length,Mean Quality Score" > "${PROJECT_DIR}/summary/fastqc_summary.csv"
fi

# Append summary data
echo "Finished summarizing FastQC results for ${SAMPLE_NAME}"


#########################


# Function to evaluate FastQC data and determine quality and adapter trimming needs
evaluate_fastqc() {
	local fastqc_dir=$1
	local sample_name=$2

	# Define the path to the FastQC data file
	local fastqc_data_file="${fastqc_dir}/fastqc_data.txt"

	# Check if the FastQC data file exists
	if [[ ! -f "$fastqc_data_file" ]]; then
		echo "Error: FastQC data file does not exist for $sample_name at $fastqc_data_file" >&2
		return 1
	fi

	# Parse FastQC data to determine quality and adapter presence
	local mean_q_score=$(awk '/>>Per base sequence quality/{flag=1; next} />>END_MODULE/{flag=0} flag' "$fastqc_data_file" | awk -F'\t' '{total+=$2; count++} END {print total/count}')
	local adapter_content=$(awk '/>>Adapter Content/{flag=1; next} />>END_MODULE/{flag=0} flag' "$fastqc_data_file" | awk 'NR > 1 {if ($2 > 5) exit 1}')

	# Determine QC Status based on average quality score
	local qc_status="PASS"
	if (( $(echo "$mean_q_score < 30" | bc -l) )); then
		qc_status="FAIL"
	fi

	# Determine Adapter Status
	local adapter_status="PASS"
	if [ $? -eq 1 ]; then
		adapter_status="FAIL"
	fi

	# Update the CSV file with QC and Adapter status
	echo "$sample_name,$qc_status,$adapter_status" >> "${PROJECT_DIR}/summary/fastqc_quality_adapters.csv"
}

# Update CSV headers and process summaries for R1 and R2
if [[ ! -s "${PROJECT_DIR}/summary/fastqc_quality_adapters.csv" ]]; then
	echo "Sample,QC Status,Adapter Status" > "${PROJECT_DIR}/summary/fastqc_quality_adapters.csv"
fi

evaluate_fastqc "$FASTQC_R1_DIR" "${SAMPLE_NAME}_R1"
evaluate_fastqc "$FASTQC_R2_DIR" "${SAMPLE_NAME}_R2"

echo "Finished evaluating FastQC results for trimming needs of ${SAMPLE_NAME}."



#########################



# end
