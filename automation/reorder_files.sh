#!/bin/bash

# Directory containing the files
gen_data_dir="/home/murali/Documents/rass/automation/gen_data"

# Ensure all month directories exist
for month_num in {01..12}; do
    month_name=$(date -d "2021-$month_num-01" +%b)
    mkdir -p "$gen_data_dir/$month_name"
done

# Process files in the gen_data directory
for file_path in "$gen_data_dir"/*_wf.csv "$gen_data_dir"/*_cleaned.csv "$gen_data_dir"/*_N-S.png "$gen_data_dir"/*_E-W.png; do  
    filename=$(basename "$file_path")
    date_str=${filename:1:8}  # Assumes date is at position 1 to 8 in the filename
    if [[ $date_str =~ ^[0-9]{8}$ ]]; then
        year=${date_str:0:4}
        month=${date_str:4:2}
        day=${date_str:6:2}
        month_name=$(date -d "$year-$month-$day" +%b)
        
        # Move the file to the target directory
        target_dir="$gen_data_dir/$month_name"
        mv "$file_path" "$target_dir"
        echo "Moved $filename to $target_dir"
    else
        echo "Skipping $filename, date extraction failed."
    fi
done
