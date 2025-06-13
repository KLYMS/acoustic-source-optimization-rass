#!/bin/bash

DIRECTORY="/home/murali/Documents/rass/data/2021_weather_data"
DIRECTORY_GEN_DATA="/home/murali/Documents/rass/automation/gen_data"
CONDA_PYTHON="/home/murali/anaconda3/bin/python"

# Ensure the DIRECTORY_GEN_DATA exists
mkdir -p "$DIRECTORY_GEN_DATA"

# Compile the C++ code
g++ /home/murali/Documents/rass/automation/pde_cor_edg_auto.cpp -o /home/murali/Documents/rass/automation/pde_cor_edg_auto.o -I /usr/include/alglib -L /usr/lib -lalglib

for SUB_DIRECTORY in "$DIRECTORY"/*; do
    if [ -d "$SUB_DIRECTORY" ]; then
        for FILE in "$SUB_DIRECTORY"/*.CSV; do
            if [ -f "$FILE" ]; then
                FILENAME=$(basename "$FILE")
                FILEPATH=$(realpath "$FILE")

                # Extract the base name without extension
                BASENAME="${FILENAME%.CSV}"

                # Run the Python script to clean data
                python data_analysis_clean.py "$FILEPATH"

                # Run the compiled C++ program
                /home/murali/Documents/rass/automation/pde_cor_edg_auto.o "$DIRECTORY_GEN_DATA/${BASENAME}_cleaned.csv"

                # Run the Python parser script using Anaconda interpreter
                $CONDA_PYTHON "/home/murali/Documents/rass/automation/backscattring_cluster_parser.py" "$DIRECTORY_GEN_DATA/${BASENAME}_wf.csv" x
                $CONDA_PYTHON "/home/murali/Documents/rass/automation/backscattring_cluster_parser.py" "$DIRECTORY_GEN_DATA/${BASENAME}_wf.csv" y

                # Remove the cleaned CSV file
                # rm "$DIRECTORY_GEN_DATA/${BASENAME}_cleaned.csv"
            fi
        done
    fi
done
