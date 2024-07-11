#!/bin/bash


##
## Script to locate and return paths to specific reference data for a given genome
##


# Check if exactly two arguments are provided
if [ "$#" -ne 2 ]; then
	echo "Usage: $0 GENOME_DIR SETTING_NAME"
	echo "GENOME_DIR: The root directory where genome references are located."
	echo "SETTING_NAME: The type of reference data sought (e.g., FASTA, GTF, etc.)."
	exit 1
fi

GENOME_DIR="$1"
SETTING_NAME="$2"


#########################



# Define functions to find files, directories, and basenames
find_file() {
	local name="$1"
	# Prefer the shortest match if multiple are found
	find "$GENOME_DIR" -type f -name "$name" | sort | head -n 1
}

find_dir() {
	local name="$1"
	# Case-insensitive search
	find "$GENOME_DIR" -type d -iname "$name" | sort | head -n 1
}

find_basename() {
	local suffix="$1"
	local file=$(find "$GENOME_DIR" -type f -name "*$suffix" | sort | head -n 1)
	echo "${file%$suffix}"
}



#########################


# Main logic to locate appropriate reference data based on SETTING_NAME
case "$SETTING_NAME" in
	REF_FASTA)
		find_file "*genome.fa" ;;
	REF_DICT)
		find_file "*genome.dict" ;;
	REF_GTF)
		find_file "*genes.gtf" ;;
	REF_REFFLAT)
		find_file "*refFlat.txt.gz" ;;
	REF_CHROM_SIZES)
		find_file "*chrom.sizes" ;;
	REF_FASTQ_SCREEN_CONF)
		find_file "*fastq_screen.conf" ;;
	REF_STAR_DIR)
		find_dir "STAR" ;;
	REF_BOWTIE1_BASENAME)
		find_basename ".1.ebwt" ;;
	REF_BOWTIE2_BASENAME)
		find_basename ".1.bt2" ;;
	REF_BWA_BASENAME)
		find_basename ".bwt" ;;
	*)
		echo "Error: Unknown setting name '$SETTING_NAME'."
		exit 1
esac


#########################


# If the directory or file does not exist, exit with an error
if [ $? -ne 0 ]; then
	echo "Error: The specified GENOME_DIR does not exist or the required reference cannot be found."
	exit 1
fi


#########################



# end
