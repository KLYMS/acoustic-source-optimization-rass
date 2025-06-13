import os
import subprocess
import sys

# Default values
DEFAULT_DIRECTORY = "./data/2021_weather_data"
DEFAULT_DIRECTORY_GEN_DATA = "./gen_data"
DEFAULT_CONDA_PYTHON = "python"
DEFAULT_CPP_SOURCE = "./pde_cor_edg_auto.cpp"
DEFAULT_CPP_OUTPUT = "./pde_cor_edg_auto.o"
DEFAULT_INCLUDE_PATH = "/usr/include/alglib"
DEFAULT_LIB_PATH = "/usr/lib"
DEFAULT_LIB_NAME = "alglib"
DEFAULT_DATA_CLEAN_SCRIPT = "./data_analysis_clean.py"
DEFAULT_PARSER_SCRIPT = "./backscattring_cluster_parser.py"
DEFAULT_INPUT_EXTENSION = ".CSV"
DEFAULT_CLEANED_SUFFIX = "_cleaned.csv"
DEFAULT_WF_SUFFIX = "_wf.csv"

# Get command-line arguments or use default values
DIRECTORY = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_DIRECTORY
DIRECTORY_GEN_DATA = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_DIRECTORY_GEN_DATA
CONDA_PYTHON = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_CONDA_PYTHON
CPP_SOURCE = sys.argv[4] if len(sys.argv) > 4 else DEFAULT_CPP_SOURCE
CPP_OUTPUT = sys.argv[5] if len(sys.argv) > 5 else DEFAULT_CPP_OUTPUT
INCLUDE_PATH = sys.argv[6] if len(sys.argv) > 6 else DEFAULT_INCLUDE_PATH
LIB_PATH = sys.argv[7] if len(sys.argv) > 7 else DEFAULT_LIB_PATH
LIB_NAME = sys.argv[8] if len(sys.argv) > 8 else DEFAULT_LIB_NAME
DATA_CLEAN_SCRIPT = sys.argv[9] if len(sys.argv) > 9 else DEFAULT_DATA_CLEAN_SCRIPT
PARSER_SCRIPT = sys.argv[10] if len(sys.argv) > 10 else DEFAULT_PARSER_SCRIPT
INPUT_EXTENSION = sys.argv[11] if len(sys.argv) > 11 else DEFAULT_INPUT_EXTENSION
CLEANED_SUFFIX = sys.argv[12] if len(sys.argv) > 12 else DEFAULT_CLEANED_SUFFIX
WF_SUFFIX = sys.argv[13] if len(sys.argv) > 13 else DEFAULT_WF_SUFFIX

# Ensure the DIRECTORY_GEN_DATA exists
os.makedirs(DIRECTORY_GEN_DATA, exist_ok=True)

# Compile the C++ code
subprocess.run(["g++", CPP_SOURCE, "-o", CPP_OUTPUT, "-I", INCLUDE_PATH, "-L", LIB_PATH, "-l", LIB_NAME])

# Iterate through subdirectories and files
for sub_directory in os.listdir(DIRECTORY):
    sub_dir_path = os.path.join(DIRECTORY, sub_directory)
    if os.path.isdir(sub_dir_path):
        for file in os.listdir(sub_dir_path):
            if file.endswith(INPUT_EXTENSION):
                file_path = os.path.join(sub_dir_path, file)
                filename = os.path.basename(file_path)
                basename = filename[:-len(INPUT_EXTENSION)]

                # Run the Python script to clean data
                subprocess.run([CONDA_PYTHON, DATA_CLEAN_SCRIPT, file_path])

                # Run the compiled C++ program
                subprocess.run([CPP_OUTPUT, f"{DIRECTORY_GEN_DATA}/{basename}{CLEANED_SUFFIX}"])

                # Run the Python parser script using Anaconda interpreter
                subprocess.run([CONDA_PYTHON, PARSER_SCRIPT, f"{DIRECTORY_GEN_DATA}/{basename}{WF_SUFFIX}", "x"])
                subprocess.run([CONDA_PYTHON, PARSER_SCRIPT, f"{DIRECTORY_GEN_DATA}/{basename}{WF_SUFFIX}", "y"])

                # Optionally remove the cleaned CSV file
                # os.remove(f"{DIRECTORY_GEN_DATA}/{basename}{CLEANED_SUFFIX}")
