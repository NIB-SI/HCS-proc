#!/bin/bash

# The path to the CSV file with filenames to move
CSV_FILE="/PATH/TO/preprocessing/qc/remove/remove_image_list.csv"

# The source directory where the images currently exist
SOURCE_DIR="/PATH/TO/IMAGES/"

# The destination directory where the images will be moved
DEST_DIR="/PATH/TO/input_removed"

# Create the destination directory if it does not exist
mkdir -p "$DEST_DIR"

# Read each filename from the CSV file, skipping the header line
tail -n +2 "$CSV_FILE" | while read -r FILENAME; do
    # Construct the full path to the file
    FULL_PATH="$SOURCE_DIR/$FILENAME"

    # Debug: Print the full path
    echo "Looking for file at: $FULL_PATH"

    # Check if the file exists
    if [[ -f "$FULL_PATH" ]]; then
        # Move the file to the destination directory
        mv "$FULL_PATH" "$DEST_DIR"
        echo "Moved file: $FILENAME"
    else
        # Print an error message to stderr
        >&2 echo "File not found: $FILENAME"
    fi
done
