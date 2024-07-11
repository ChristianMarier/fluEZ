#!/bin/bash


##
## Script to modify CSV files to enhance their compatibility with bioinformatics tools
##


# Check if exactly one argument (the CSV file name) is provided
if [ "$#" -ne 1 ]; then
	echo "Error: Incorrect usage."
	echo "Usage: $0 <CSV-file>"
	exit 1
fi

# Variable for the CSV file
CSV_FILE="$1"

# Check if the specified CSV file exists and is not empty
if [ ! -s "$CSV_FILE" ]; then
	echo "Error: Specified CSV file does not exist or is empty: $CSV_FILE"
	exit 1
fi

# Load the necessary module for dos2unix
module load dos2unix/7.4.0

# Normalize text file line endings
dos2unix "$CSV_FILE"
mac2unix "$CSV_FILE"

# Only modify commas within quotes and ensure not to disrupt entire data structure
awk -F'"' '{
	for (i=2; i<NF; i+=2) gsub(",", "-", $i)
} {OFS="\""; print}' "$CSV_FILE" > temp.csv && mv temp.csv "$CSV_FILE"

# Remove all double-quote characters 
#sed -i -e 's/"//g' "$CSV_FILE"

# Remove empty lines
grep -vE '^,+$' "$CSV_FILE" > temp.csv && mv temp.csv "$CSV_FILE"

# Ensure the CSV file ends with a newline character
sed -i -e '$a\' "$CSV_FILE"

echo "CSV file has been cleaned and standardized: $CSV_FILE"



# end
