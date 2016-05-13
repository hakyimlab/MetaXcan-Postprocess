#!/bin/bash
# Helper bash script that launches run_MetaXcanPostprocessing file 

if [ -z $1 ]; then
    echo "Usage: ./run.sh <project name>"
    echo "(only single word, such as breastcancer, diabetase, ovarycancer)"

else
    python run_MetaXcanPostprocessing.py $1
fi