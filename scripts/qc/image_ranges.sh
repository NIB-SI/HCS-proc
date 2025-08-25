#!/bin/bash

# Define variables
IMAGE_DIR="/PATH/TO/IMAGES/"
NUM_INSTANCES=20  # Number of parallel instances
OUTPUT_FILE="/PATH/TO/IMAGES/preprocessing/image_ranges_qc.txt"

# Count TIFF files in all subdirectories and calculate ranges
TOTAL_FILES=$(find $IMAGE_DIR -type f -name '*.tif' | wc -l)
TOTAL_IMAGES=$((TOTAL_FILES / 3))  # Every 3 files constitute one image
IMAGES_PER_INSTANCE=$((TOTAL_IMAGES / NUM_INSTANCES))

# Create ranges
> $OUTPUT_FILE
for i in $(seq 1 $NUM_INSTANCES); do
    FIRST_IMAGE=$(( (i - 1) * IMAGES_PER_INSTANCE + 1 ))
    LAST_IMAGE=$(( i * IMAGES_PER_INSTANCE ))

    # Adjust for the last instance to cover all images
    if [ $i -eq $NUM_INSTANCES ]; then
        LAST_IMAGE=$TOTAL_IMAGES
    fi

    echo "$FIRST_IMAGE $LAST_IMAGE" >> $OUTPUT_FILE
done

echo "Image ranges written to $OUTPUT_FILE"
