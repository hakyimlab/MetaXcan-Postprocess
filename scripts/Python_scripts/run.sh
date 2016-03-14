#!/bin/bash
# Helper bash script that launches run_MetaXcanPostprocessing file 

if [ -z $1 ]; then
    echo "Usage: ./run.sh <input_file_name>"
    echo "(but do not include the file extension such as .csv"

else
    python run_MetaXcanPostprocessing.py $1 $1.log
fi