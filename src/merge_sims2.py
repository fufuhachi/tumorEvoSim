#!/bin/python
import numpy as np
import pandas as pd
import sys
import classes
import simulation
import os
import multiprocessing
import sys
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count

def process_folder(folder):
    try:
        rep = int(folder.split('_')[-1])
        for file in os.listdir(folder):
            if file.split('_')[-1].startswith('time'):
                sim = classes.load_object(os.path.join(folder,file))
                summary = simulation.tumor_summary(sim.tumor, rep = rep)
                summary['t'] = sim.tumor.t
                return summary
    except(FileNotFoundError):
        print(f'{folder} not found. Trying next folder')

if __name__ == '__main__':
    folder_list = sys.argv[1]
    try: 
        max_rep = int(sys.argv[2])
    except(IndexError):
        max_rep = np.infty #change to desired number of replicates if less than total 
    print(f'max rep: {max_rep}')

    num_workers = int(sys.argv[3]) if len(sys.argv) > 3 else cpu_count()
    print(f'using {num_workers} workers')

    completed = 0
    with open(folder_list, 'r') as f:
        folders = f.readlines()
        folders = [folder.strip('\n') for folder in folders]
        folders = folders[:max_rep] # limit to max_rep number of folders
    print(f'processing {len(folders)} folders')
    
    with Pool(num_workers) as pool:
        dataframes = pool.map(process_folder, folders)
        dataframes = [df for df in dataframes if df is not None] # remove None values
    
    print('concatenating all dataframes')
    timedata = pd.concat(dataframes)
    print('writing to file')
    timedata.to_csv(f'experiment_nreps_{len(dataframes)}.csv')
    print('done')
