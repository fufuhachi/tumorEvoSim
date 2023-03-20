"""script to make age v mean r on time-normalized runs"""
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns 
import sys
import os
import numpy as np

def plot_agevr_diff(data, outdir):
    t_all = data['norm_t_binned'].unique()
    t_all.sort()
    for t in t_all:
        t
        data = data[data['norm_t_binned']==t]
        for cutoff in [0,.01,.1]:
            visible = data[(data['blood'] > cutoff) | (data['tissue'] > cutoff)]
            sns.scatterplot(data = visible, x = 'centroid_r', y = 'norm_age', hue = 'diff', size = 'tissue', alpha = .7)
            plt.show(block = False)
            plt.savefig(os.path.join(outdir, f'agevr_diff_t_{t}_{cutoff}.png'))
            plt.close()

def plot_agevr_pcterr(data, outdir):
    t_all = data['norm_t_binned'].unique()
    t_all.sort()
    for t in t_all:
        for cutoff in [0,.01,.1]:
            visible = data[((data['blood'] > cutoff) | (data['tissue'] > cutoff)) & (data['norm_t_binned']==t)]
            sns.scatterplot(data = visible, x = 'centroid_r', y = 'norm_age', hue = 'pcterr', size = 'tissue', alpha = .7)
            plt.show(block = False)
            plt.savefig(os.path.join(outdir, f'agevr_pcterr_t_{t}_{cutoff}.png'))
            plt.close()
            
if __name__ == '__main__':
    print('loading data')
    data = pd.read_csv(sys.argv[1])
    outdir = '' if len(sys.argv) < 3 else sys.argv[2]

    data = pd.read_csv(sys.argv[1])
    #make more bins 
    binned_time = pd.cut(data['norm_t'],20).apply(lambda x: np.round(x.mid,2))
    data['norm_t_binned'] = binned_time.astype(float)
    #choose median sample per rep per bin 
    med_vals = data.groupby(['norm_t_binned','rep'])['t'].median().reset_index()[['rep','t']] #reset_index to convert groupby object to DataFrame
    # get the index of median column in the original dataframe
    med_df = pd.merge(data, med_vals, on  = ['rep','t'], how = 'inner')
    data = med_df
    plot_agevr_pcterr(data, outdir)