#!/bin/bash
# Helper bash script that launches run_MetaXcanPostprocessing file 

if [ -z $1 ]; then
    echo "Usage: ./run.sh <project name>"
    echo "(e.g ./run.sh ovarican_cancer)"

else
    python run_MetaXcanPostprocessing.py $1
fi