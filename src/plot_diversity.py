#!/bin/python 
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt 
import numpy as np 
import sys
import os 

def get_div_old(data):
    tiss = data.groupby(['norm_t_binned','rep','is_outer'])['tissue'].transform(lambda x: inv_simpson_ind(x)).reset_index()['tissue']
    bl = data.groupby(['norm_t_binned','rep','is_outer'])['blood'].transform(lambda x: inv_simpson_ind(x)).reset_index()['blood']
    data['tissue_div'] = tiss
    data['blood_div'] = bl
    return data
def get_div(data):
    data['tiss_inner_div'] = data.groupby(['rep','norm_t_binned'])['tissue_inner_count'].transform(inv_simpson_ind)
    data['tiss_outer_div'] = data.groupby(['rep','norm_t_binned'])['tissue_outer_count'].transform(inv_simpson_ind)
    data['bl_inner_div'] = data.groupby(['rep','norm_t_binned'])['blood_inner_count'].transform(inv_simpson_ind)
    data['bl_outer_div'] = data.groupby(['rep','norm_t_binned'])['blood_outer_count'].transform(inv_simpson_ind)
    return data

def plot_div_old(div_data, outdir):
    """plot showing mean entropy by deme in blood and tissue"""
    sns.boxplot(data = div_data, x = 'norm_t_binned', y = 'tissue_div', hue = 'is_outer')
    #sns.boxplot(data = div_data, x = 'norm_t_binned', y = 'blood_div',color = 'red',hue = 'is_outer')
    plt.xticks(rotation = 45)
    plt.savefig(os.path.join(outdir,'div.png'))
    plt.show(block = False)
    plt.close()

def plot_div(data, outdir):
    sns.lineplot(data = data, x = 'norm_t_binned',y = 'tiss_inner_div')#,label = 'tiss inner')
    sns.lineplot(data = data, x = 'norm_t_binned',y = 'tiss_outer_div')#,label = 'tiss outer')
    sns.lineplot(data = data, x = 'norm_t_binned',y = 'bl_inner_div')#,label = 'blood inner')
    sns.lineplot(data = data, x = 'norm_t_binned',y = 'bl_outer_div')#,label = 'blood outer')
    plt.xticks(rotation = 45)
    plt.savefig(os.path.join(outdir,'div.png'))
    plt.show(block = False)
    plt.close()

def inv_simpson_ind(vec):
    """input: vector of frequencies.
    output: inverse simpson diversity index"""
    if vec is None or len(vec) == 0 or vec.sum() == 0:
        return 0
    vec = vec/vec.sum()
    return 1/(np.power(vec, 2).sum())

if __name__ == '__main__':
    data = pd.read_csv(sys.argv[1])
    outdir = '' if len(sys.argv) < 3 else sys.argv[2]
    out = get_div(data)
    plot_div(out, outdir)
    print('done')