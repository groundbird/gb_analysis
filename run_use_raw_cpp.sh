#!/bin/bash

# Define an array of measurement IDs
MEAS_IDS=(

    8448
    
)

# Path to the Python script
PYTHON_SCRIPT="./use_raw_cpp.py"

# Check if the Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Python script not found at $PYTHON_SCRIPT"
    exit 1
fi

# Process each measurement ID
for meas_id in "${MEAS_IDS[@]}"; do
    echo "Processing measurement ID: $meas_id"
    python "$PYTHON_SCRIPT" "$meas_id"

    # Check if the script ran successfully
    if [ $? -eq 0 ]; then
        echo "Successfully processed measurement ID: $meas_id"
    else
        echo "Error occurred while processing measurement ID: $meas_id"
    fi

    echo "----------------------------------------"
done

echo "All measurements have been processed."
