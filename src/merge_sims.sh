#!/bin/bash
# Accepts a filename as an argument for the first command

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <experiment_folders_list>"
    exit 1
fi

experiment_folders_list=$1
all_files="all_files.txt"
summary_csv="summary.csv"
analysis="analysis"
tar_file="analysis.tar.gz"

# First command with filename argument
sh ~/sim2/tumorEvoSim/src/create_file_list.sh $experiment_folders_list

python ~/sim2/tumorEvoSim/src/merge_sims3.py $all_files

python ~/sim2/tumorEvoSim/src/process_timeseries.py $summary_csv

tar -czvf $tar_file $analysis
tar -czvf "$summary_csv".tar.gz $summary_csv  
tar -czvf comp.tar.gz analysis/comp.csv 
