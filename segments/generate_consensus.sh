#!/bin/bash

# Script to generate a consensus sequence from VCF files for the FluEZ pipeline.

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


## Define paths for input files
VCF_PATH="${PROJECT_DIR}/variants/${SAMPLE_NAME}_filtered.vcf.gz"
BAM_PATH="${PROJECT_DIR}/bam/${SAMPLE_NAME}_dedup.bam"
REF_FASTA=$(grep "^REF_FASTA" "$PROJECT_DIR/settings.txt" | cut -d '=' -f2)

echo "Setting up environment for consensus generation..."
echo "Checking input files..."

## Check if the VCF file exists
if [ ! -f "$VCF_PATH" ]; then
    echo "Error: VCF file does not exist at ${VCF_PATH}"
    exit 1
fi

## Check if the reference genome file exists
if [ ! -f "$REF_FASTA" ]; then
    echo "Error: Reference genome file does not exist at ${REF_FASTA}"
    exit 1
fi

## Check if the BAM file exists
if [ ! -f "$BAM_PATH" ]; then
    echo "Error: BAM file does not exist at ${BAM_PATH}"
    exit 1
fi

echo "Inputs verified: VCF and reference genome files are in place."
echo "Project Directory: $PROJECT_DIR"
echo "Sample Name: $SAMPLE_NAME"
echo "Number of Threads: $THREADS"
echo "VCF Path: $VCF_PATH"
echo "BAM Path: $BAM_PATH"
echo "Reference Genome Path: $REF_FASTA"


#######################


echo "Configuring environment for consensus sequence generation..."

## Clear any previously loaded modules to avoid conflicts
module purge

## Load necessary modules
module load bedtools/2.30.0
module load bcftools/1.15.1

echo "Modules loaded: bedtools and bcftools."

## Directory structure setup
# Create directory for storing final consensus sequences and logs
CONS_SEQ_DIR="${PROJECT_DIR}/consensus_sequences"

echo "Setting up directory..."

# Check and create directories if they do not exist
mkdir -p $CONS_SEQ_DIR

echo "Directories set up:"
echo "Consensus Sequences Directory: $CONS_SEQ_DIR"

echo "Starting consensus sequence generation for sample: $SAMPLE_NAME..."

## Path setup
OUTPUT_CONSENSUS="${CONS_SEQ_DIR}/${SAMPLE_NAME}_unmasked_consensus.fa"
MASK_BED="${CONS_SEQ_DIR}/${SAMPLE_NAME}_mask.bed"
MASKED_CONSENSUS="${CONS_SEQ_DIR}/${SAMPLE_NAME}_masked_consensus.fa"

## Generate consensus sequence
echo "Generating consensus sequence using bedtools and bcftools..."
bedtools genomecov -bga -ibam "$BAM_PATH" -g "$REF_FASTA" | awk '$4 < 100' | bedtools merge > "$MASK_BED"
sleep 2
bcftools consensus -f "$REF_FASTA" -o "$OUTPUT_CONSENSUS" "$VCF_PATH"
sleep 2
bcftools consensus -f "$REF_FASTA" -m "$MASK_BED" -o "$MASKED_CONSENSUS" "$VCF_PATH"

sleep 2

if [ $? -ne 0 ]; then
    echo "Error: Consensus sequence generation failed for $SAMPLE_NAME."
    exit 1
fi

echo "Consensus sequence generated at: $MASKED_CONSENSUS"


#######################


## Quality checks
echo "Performing quality checks on the generated consensus sequence..."

# Validate sequence length and absence of undefined bases
SEQ_LENGTH=$(grep -v '>' "$MASKED_CONSENSUS" | tr -d '\n' | wc -c)
UNDEFINED_BASES=$(grep -v '>' "$MASKED_CONSENSUS" | tr -d '\n' | grep -o 'N' | wc -l)

echo "Total sequence length: $SEQ_LENGTH"
echo "Number of undefined bases (N's): $UNDEFINED_BASES"

if [ "$UNDEFINED_BASES" -gt 0 ]; then
    echo "Warning: There are undefined bases in the consensus sequence."
    # Additional steps might include re-running with adjusted parameters or further investigation
fi

echo "Consensus sequence quality checks completed. Results stored in log files."


#######################

# end



