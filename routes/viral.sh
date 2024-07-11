#!/bin/bash


##
## Influenza A virus whole genome sequencing variant discovery using Bowtie2, bcftools, samtools, and custom scripts
##


# Ensure the script is called with the correct number of arguments
if [ "$#" -ne 2 ]; then
	echo "Usage: $0 <project_directory> <sample_name>"
	exit 1
fi

# Variables setup
PROJECT_DIR=$1
SAMPLE_NAME=$2
THREADS=$SLURM_CPUS_PER_TASK
THREADS=$(( THREADS - 1 )) # Reserve a thread for overhead
CODE_DIR=${CODE_DIR:-$(dirname $(dirname $(readlink -f $0)))} # Directory of the script

echo "Starting Influenza A virus sequencing analysis for sample: $SAMPLE_NAME"
echo "Project Directory: $PROJECT_DIR"
echo "Sample Name: $SAMPLE_NAME"
echo "Number of Threads: $THREADS"
echo "Code Directory: $CODE_DIR"


#########################


# Segment 1: FASTQ Preparation
echo "Preparing FASTQ files..."
bash $CODE_DIR/segments/prepare_fastqs.sh $PROJECT_DIR $SAMPLE_NAME
if [ $? -ne 0 ]; then
	echo "FASTQ preparation failed for $SAMPLE_NAME."
	exit 1
fi

# Segment 2: Quality Control with FastQC
echo "Running FastQC..."
bash $CODE_DIR/segments/run_fastqc.sh $PROJECT_DIR $SAMPLE_NAME $THREADS
if [ $? -ne 0 ]; then
	echo "FastQC analysis failed for $SAMPLE_NAME."
	exit 1
fi

# Segment 3: Trimming
echo "Trimming FASTQ files..."
bash $CODE_DIR/segments/trim_fastqs.sh $PROJECT_DIR $SAMPLE_NAME $THREADS
if [ $? -ne 0 ]; then
	echo "Trimming failed for $SAMPLE_NAME."
	exit 1
fi

# Segment 4: Trimming Quality Control with FastQC
echo "Running FastQC on trimmed FASTQ files..."
bash $CODE_DIR/segments/run_fastqc_trimmed.sh $PROJECT_DIR $SAMPLE_NAME $THREADS
if [ $? -ne 0 ]; then
	echo "Trimming FastQC analysis failed for $SAMPLE_NAME."
	exit 1
fi

# Segment 5: Alignment using Bowtie2
echo "Aligning trimmed FASTQ files to reference genome..."
bash $CODE_DIR/segments/align_bowtie2.sh $PROJECT_DIR $SAMPLE_NAME $THREADS
if [ $? -ne 0 ]; then
	echo "Alignment failed for $SAMPLE_NAME."
	exit 1
fi

# Segment 6: Duplicate Removal using Sambamba
echo "Removing duplicates from aligned BAM files..."
bash $CODE_DIR/segments/remove_duplicates.sh $PROJECT_DIR $SAMPLE_NAME $THREADS
if [ $? -ne 0 ]; then
	echo "Duplicate removal failed for $SAMPLE_NAME."
	exit 1
fi

# Segment 7: Quality Control Post-processing
echo "Performing post-processing quality control..."
bash $CODE_DIR/segments/qc_post_processing.sh $PROJECT_DIR $SAMPLE_NAME $THREADS
if [ $? -ne 0 ]; then
	echo "Post-processing QC failed for $SAMPLE_NAME."
	exit 1
fi

# Segment 8: Variant Calling
echo "Performing variant calling..."
bash $CODE_DIR/segments/variant_calling.sh $PROJECT_DIR $SAMPLE_NAME $THREADS
if [ $? -ne 0 ]; then
	echo "Variant calling failed for $SAMPLE_NAME."
	exit 1
fi

# Segment 9: Consensus Sequence Generation
echo "Generating consensus sequence..."
bash $CODE_DIR/segments/generate_consensus.sh $PROJECT_DIR $SAMPLE_NAME $THREADS
if [ $? -ne 0 ]; then
	echo "Consensus sequence generation failed for $SAMPLE_NAME."
	exit 1
fi


#########################


echo "Analysis complete for $SAMPLE_NAME."
echo "Date and Time: $(date)"


# end
