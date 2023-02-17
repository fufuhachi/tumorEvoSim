#!/bin/python 
import numpy as np
import pandas as pd
import sys
import classes
import simulation
from multiprocessing import Pool, cpu_count


def process_file(file):
    print(f'processing {file}')
    if file.split('_')[-1].startswith('time'):
        sim = classes.load_object(file)
        rep = file.split('rep=')[-1].split('_')[0]
        summary = simulation.tumor_summary(sim.tumor, rep = rep)
        summary['t'] = sim.tumor.t
        return summary

if __name__ == "__main__":
    with open(sys.argv[1],'r') as f:
        files = f.splitlines()
    num_workers = int(sys.argv[2]) if len(sys.argv) > 3 else cpu_count()
    print(f'using {num_workers} workers')
    with Pool(num_workers) as pool:
        dataframes = pool.map(process_file, files)
    dataframes = [df for df in dataframes if df is not None] 
    print('concatenating all dataframes')
    timedata = pd.concat(dataframes)
    print('writing to file')
    timedata.to_csv(f'summary.csv')
    print('done')