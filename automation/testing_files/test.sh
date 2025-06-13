#!/bin/bash

# Configurable variables
DIRECTORY="/home/murali/Documents/rass/data/2021_weather_data"
DIRECTORY_GEN_DATA="/home/murali/Documents/rass/automation/gen_data"
CONDA_PYTHON="/home/murali/anaconda3/bin/python"
CPP_SOURCE="/home/murali/Documents/rass/automation/pde_cor_edg_auto.cpp"
CPP_OUTPUT="/home/murali/Documents/rass/automation/pde_cor_edg_auto.o"
INCLUDE_PATH="/usr/include/alglib"
LIB_PATH="/usr/lib"
LIB_NAME="alglib"
DATA_CLEAN_SCRIPT="data_analysis_clean.py"
PARSER_SCRIPT="/home/murali/Documents/rass/automation/backscattring_cluster_parser.py"
INPUT_EXTENSION=".CSV"
CLEANED_SUFFIX="_cleaned.csv"
WF_SUFFIX="_wf.csv"

# Ensure the DIRECTORY_GEN_DATA exists
mkdir -p "$DIRECTORY_GEN_DATA"

# Compile the C++ code
g++ "$CPP_SOURCE" -o "$CPP_OUTPUT" -I "$INCLUDE_PATH" -L "$LIB_PATH" -l "$LIB_NAME"

# Iterate through subdirectories and files
for SUB_DIRECTORY in "$DIRECTORY"/*; do
    if [ -d "$SUB_DIRECTORY" ]; then
        for FILE in "$SUB_DIRECTORY"/*"$INPUT_EXTENSION"; do
            if [ -f "$FILE" ]; then
                FILENAME=$(basename "$FILE")
                FILEPATH=$(realpath "$FILE")

                # Extract the base name without extension
                BASENAME="${FILENAME%$INPUT_EXTENSION}"

                # Run the Python script to clean data
                python "$DATA_CLEAN_SCRIPT" "$FILEPATH"

                # Run the compiled C++ program
                "$CPP_OUTPUT" "$DIRECTORY_GEN_DATA/${BASENAME}$CLEANED_SUFFIX"

                # Run the Python parser script using Anaconda interpreter
                "$CONDA_PYTHON" "$PARSER_SCRIPT" "$DIRECTORY_GEN_DATA/${BASENAME}$WF_SUFFIX" x
                "$CONDA_PYTHON" "$PARSER_SCRIPT" "$DIRECTORY_GEN_DATA/${BASENAME}$WF_SUFFIX" y

                # Optionally remove the cleaned CSV file
                # rm "$DIRECTORY_GEN_DATA/${BASENAME}$CLEANED_SUFFIX"
            fi
        done
    fi
done
