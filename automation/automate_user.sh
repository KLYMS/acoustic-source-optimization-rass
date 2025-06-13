#!/bin/bash

# making it user friendly

# Default values
DEFAULT_DIRECTORY="./data/2021_weather_data"
DEFAULT_DIRECTORY_GEN_DATA="./gen_data"
DEFAULT_CONDA_PYTHON="python"
DEFAULT_CPP_SOURCE="./pde_cor_edg_auto.cpp"
DEFAULT_CPP_OUTPUT="./pde_cor_edg_auto.o"
DEFAULT_INCLUDE_PATH="/usr/include/alglib"
DEFAULT_LIB_PATH="/usr/lib"
DEFAULT_LIB_NAME="alglib"
DEFAULT_DATA_CLEAN_SCRIPT="./data_analysis_clean.py"
DEFAULT_PARSER_SCRIPT="./backscattring_cluster_parser.py"
DEFAULT_INPUT_EXTENSION=".CSV"
DEFAULT_CLEANED_SUFFIX="_cleaned.csv"
DEFAULT_WF_SUFFIX="_wf.csv"

# Get command-line arguments or use default values
DIRECTORY="${1:-$DEFAULT_DIRECTORY}"
DIRECTORY_GEN_DATA="${2:-$DEFAULT_DIRECTORY_GEN_DATA}"
CONDA_PYTHON="${3:-$DEFAULT_CONDA_PYTHON}"
CPP_SOURCE="${4:-$DEFAULT_CPP_SOURCE}"
CPP_OUTPUT="${5:-$DEFAULT_CPP_OUTPUT}"
INCLUDE_PATH="${6:-$DEFAULT_INCLUDE_PATH}"
LIB_PATH="${7:-$DEFAULT_LIB_PATH}"
LIB_NAME="${8:-$DEFAULT_LIB_NAME}"
DATA_CLEAN_SCRIPT="${9:-$DEFAULT_DATA_CLEAN_SCRIPT}"
PARSER_SCRIPT="${10:-$DEFAULT_PARSER_SCRIPT}"
INPUT_EXTENSION="${11:-$DEFAULT_INPUT_EXTENSION}"
CLEANED_SUFFIX="${12:-$DEFAULT_CLEANED_SUFFIX}"
WF_SUFFIX="${13:-$DEFAULT_WF_SUFFIX}"

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
