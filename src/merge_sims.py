#script to convert .pkl files into tumor summaries 
import numpy as np
import pandas as pd
import sys
import classes
from main import *
import os

#functions
def pkl_to_npz(path, param_names, params):
    """Accepts path to .pkl file and simulation params that generated it. Returns pandas dataframe with 
    columns of all information """
    if path.split('.')[-1] !='pkl':
        print('path must be to .pkl\nexiting...')
        sys.exit()
    tumor = classes.load_object(path)
    df = tumor_summary(tumor, param_names, params)
    df.to_csv(path.split('.')[-1]+'.npz')
def tumor_summary(tumor, param_names, params):
    """turn tumor into dataframe with all the information: 
        ncells x 5 matrix where columns are 'cell_ID' 'x','y' 'r' 'angle' 'genotype' 'n_drivers' 'drivers'"""
    #turn list of cells into columns of information 
    mat = tumor.graph
    cell_ID = mat[mat > 0]
    x, y = np.indices(mat.shape)
    x = x[mat > 0] - tumor.center[0]
    y = y[mat > 0] - tumor.center[1]
    r = np.sqrt(x**2 + y**2)
    angle = (360/2/np.pi)*np.arctan2(y,x)
    
    genotype = np.array([tumor.cells.get_item(id).gen.ID for id in cell_ID], dtype = int)
    n_drivers = np.array([tumor.cells.get_item(id).gen.n_drivers for id in cell_ID], dtype = int)
    drivers = np.array([tuple(tumor.cells.get_item(id).gen.drivers) for id in cell_ID],dtype = object)
   
    df = pd.DataFrame({'x':x,'y':y,'r':r,'angle':angle, 'genotype':genotype,'n_drivers': n_drivers, 'drivers':drivers})
    for name, par in zip(param_names, params):
        df[name] = par
    return df

if __name__ == '__main__':
    to_concat = []
    folder_list = sys.argv[1]
    catted_file = sys.argv[2]
    for folder in np.loadtxt(folder_list, dtype = 'str'):
        for i in range(100):
            try:
                tumor = classes.load_object(f'{folder}/rep={i}_ncells=5000.pkl').tumor
            except(FileNotFoundError):
                print(f'{folder}/rep={i}_ncells=5000.pkl not found')
                continue
            sf = folder.split('_')
            param_names = (sf[0],sf[2],sf[4], 'rep')
            params = (sf[1],sf[3],sf[5], i)
            df = tumor_summary(tumor, param_names, params)
            to_concat.append(df)
            
        print(f'{folder} done')

    catted = pd.concat(objs = to_concat, ignore_index=True)
    catted.to_csv(f'{catted_file}.npy')



    
#make text file with file names, 

