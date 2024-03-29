#script to convert .pkl files into tumor summaries 
import numpy as np
import pandas as pd
import sys
import classes
import simulation
import os
import multiprocessing
#functions
#depricated, use version in simulation.py
'''def tumor_summary(tumor, rep = 0):
    """turn tumor into dataframe with all the information: 
        ncells x 5 matrix where columns are 'cell_ID' 'x','y' 'r' 'angle' 'genotype' 'n_drivers' 'drivers' """
    #turn list of cells into columns of information 
    b = tumor.params['init_birth_rate']
    dr_params = tumor.params['dr_params']
    decay = 33
    mat = tumor.graph
    cell_ID = mat[mat > 0]
    x, y = np.indices(mat.shape)
    x = x[mat >0] - tumor.center[0]
    y = y[mat>0] - tumor.center[1]
    r = np.sqrt(x**2 + y**2)
    angle = (360/2/np.pi)*np.arctan2(y,x)
    
    genotype = np.array([tumor.cells.get_item(id).gen.ID for id in cell_ID], dtype = int)
    n_drivers = np.array([tumor.cells.get_item(id).gen.n_drivers for id in cell_ID], dtype = int)
    drivers = [tuple(tumor.cells.get_item(id).gen.drivers) for id in cell_ID]
    
    death_rate = [tumor.cells.get_item(id).get_death_rate(**dr_params) for id in cell_ID]
    
    df = pd.DataFrame({'cell_ID' : cell_ID, 'x':x,'y':y,'r':r,'angle':angle, 
    'genotype':genotype,'n_drivers': n_drivers, 'drivers':drivers, 'death_rate': death_rate})
    df['cell_hge'] = df['death_rate']/(b - df['death_rate'] + decay)
    df['rep'] = rep
    return df'''


if __name__ == '__main__':
    folder_list = sys.argv[1]
    try: 
        max_rep = int(sys.argv[2])
    except(IndexError):
        max_rep = np.infty #change to desired number of replicates if less than total 
    print(f'max rep: {max_rep}')

    completed = 0
    with open(folder_list, 'r') as f:
        dataframes = []
        for folder in f.readlines():
            folder = folder.strip('\n')
            print(folder)
            rep = int(folder.split('_')[-1])
            if completed >= max_rep:
                break
            else: 
                try:
                    for file in os.listdir(folder):
                        if file.split('_')[-1].startswith('time'):
                            
                            sim = classes.load_object(os.path.join(folder,file))
                            summary = simulation.tumor_summary(sim.tumor, rep = rep)
                            summary['t'] = sim.tumor.t
                            dataframes.append(summary)
                    completed+=1
                    print(f'{folder} dataframes created')
                except(FileNotFoundError):
                    print(f'{folder} not found. Trying next folder')

        print('concatenating all dataframes')
        timedata = pd.concat(dataframes)
        print('writing to file')
        timedata.to_csv(f'experiment_nreps_{completed}.csv')
        print('done')

    
    



    
#make text file with file names, 

