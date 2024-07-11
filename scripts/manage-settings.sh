#!/bin/bash


##
## Script to manage settings within the specified "settings.txt" file
##

# Usage instructions
usage() {
	echo "Usage: $0 <settings-file> <setting-name> [new-value]"
	echo " - <settings-file>: Path to settings.txt file"
	echo " - <setting-name>: Name of the setting to retrieve or set"
	echo " - [new-value]: Optional. New value to set for the given setting"
	exit 1
}

# Verify the correct number of arguments
if [ $# -lt 2 ]; then
	echo "Error: Incorrect number of arguments supplied."
	usage
fi

# Variables
SETTINGS_FILE="$1"
SETTING_NAME="$2"
NEW_VALUE="$3"


#########################


# Check if the specified settings file exists and is not empty
if [ ! -s "$SETTINGS_FILE" ]; then
	echo "Error: Specified settings file does not exist or is empty."
	exit 1
fi


#########################


# Attempt to retrieve the current value of the specified setting
CURRENT_VALUE=$(grep "^$SETTING_NAME=" "$SETTINGS_FILE" | cut -d '=' -f2)

if [ -n "$CURRENT_VALUE" ]; then
	echo "$CURRENT_VALUE"
	exit 0
elif [ -z "$NEW_VALUE" ] && [[ "$SETTING_NAME" == REF_* ]]; then
	# Compute the setting value dynamically for reference settings
	GENOME_DIR=$(grep "^GENOME_DIR=" "$SETTINGS_FILE" | cut -d '=' -f2)
	if [ -z "$GENOME_DIR" ]; then
		echo "Error: GENOME_DIR setting is not set."
		exit 1
	fi
	NEW_VALUE=$(bash "${BASH_SOURCE%/*}/locate-reference.sh" "$GENOME_DIR" "$SETTING_NAME")
	echo "$SETTING_NAME=$NEW_VALUE" >> "$SETTINGS_FILE"
elif [ -n "$NEW_VALUE" ]; then
	# Use the provided new value
	echo "$SETTING_NAME=$NEW_VALUE" >> "$SETTINGS_FILE"
else
	echo "Error: No new value provided for $SETTING_NAME and it cannot be computed dynamically."
	exit 1
fi


#########################


# Attempt to retrieve the setting again to confirm it was updated
CONFIRMED_VALUE=$(grep "^$SETTING_NAME=" "$SETTINGS_FILE" | cut -d "=" -f2)
if [ -n "$CONFIRMED_VALUE" ]; then
	echo "$CONFIRMED_VALUE"
else
	echo "Error: Failed to set the value for $SETTING_NAME."
fi


#########################



# end
