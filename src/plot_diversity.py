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
    data['tiss_whole_div'] = data.groupby(['rep','norm_t_binned'])['tissue'].transform(inv_simpson_ind)
    data['bl_whole_div'] = data.groupby(['rep','norm_t_binned'])['blood'].transform(inv_simpson_ind)
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
    print('plotting div boxplot')
    #sns.lineplot(data = data, x = 'norm_t_binned',y = 'tiss_inner_div',label = 'tiss inner')
    #sns.lineplot(data = data, x = 'norm_t_binned',y = 'tiss_outer_div',label = 'tiss outer')
    #sns.lineplot(data = data, x = 'norm_t_binned',y = 'bl_inner_div',label = 'blood inner')
    #sns.lineplot(data = data, x = 'norm_t_binned',y = 'bl_outer_div',label = 'blood outer')

    df = pd.melt(data, id_vars = ['norm_t_binned'], value_vars = ['tiss_whole_div','bl_whole_div']) #hue_kws = {'label': {'tiss_whole_div'}})
    sns.boxplot(data = df, x = 'norm_t_binned', y = 'value', hue = 'variable')
    plt.xticks(rotation = 45)
    plt.title('Inverse Simpson Diversity')
    plt.xlabel(['time (normalized)'])
    
    plt.savefig(os.path.join(outdir,'div.png'))
    plt.show(block = False)
    plt.close()
def plot_div_line(data, outdir):
    print('plotting div lineplot')
    data = data[['norm_t_binned', 'tiss_whole_div', 'bl_whole_div']]
    data.columns = ['norm_t_binned', 'tissue', 'blood']
    df = pd.melt(data, id_vars = ['norm_t_binned'], value_vars = ['tissue','blood']) #hue_kws = {'label': {'tiss_whole_div'}})
    
    sns.lineplot(data = df, x = 'norm_t_binned', y = 'value', hue = 'variable', err_style = 'band', errorbar = 'sd', palette = ['b', 'r'], )
    #sns.lineplot(data, x = 'norm_t_binned', y = 'tissue', color = 'blue', err_style = 'band', errorbar = 'sd')
    #sns.lineplot(data, x = 'norm_t_binned', y = 'blood', color = 'red', err_style = 'band', errorbar = 'sd')
    plt.xticks(rotation = 45)
    plt.title('Blood vs Tissue ITH')
    plt.ylabel('Inv. Simpson D')
    plt.xlabel('normalized time')
    #plt.legend(labels=['tissue','blood'])
    plt.savefig(os.path.join(outdir,'div_line.png'))
    plt.show(block = False)
    plt.close()
def inv_simpson_ind(vec):
    cutoff  = 0
    """input: vector of frequencies.
    output: inverse simpson diversity index"""
    vec = vec/vec.sum()
    vec[vec < cutoff] = 0
    if vec is None or len(vec) == 0 or vec.sum() == 0:
        return 0
    return 1/(np.power(vec, 2).sum())
def cf_hist(data, outdir,cutoff = 0):
    """create violin plot of clone fractions in blood and tissue over normed time
    --data: median-sampled time-normalized clone frequencies"""
    
    tisslab = 'log10_cf_tissue'
    bllab = 'log10_cf_blood'
    data[tisslab] = data['tissue'][data['tissue']>cutoff]
    data[bllab] = data['blood'][data['blood']>cutoff]
    df = pd.melt(data, id_vars = ['norm_t_binned'], value_vars = [tisslab, bllab]) 
    sns.violinplot(data = df, x = 'norm_t_binned', y = 'value', hue = 'variable',cut = 0, split = False)
    plt.xlabel('normed time')
    plt.ylabel('clone fraction')
    if cutoff >0:
        plt.title(f'cutoff = {cutoff}')
    plt.savefig(os.path.join(outdir,'cfhist.png'))
    plt.show(block = False)
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
    
    #sample = data[for data[['rep','t']] in med_vals]
    #print(sample)
    #print(data.shape, med_vals.shape)
    outdir = sys.argv[2].strip() if len(sys.argv) > 2 else ''
    out = get_div(med_df)
    #plot_div(out, outdir)
    plot_div_line(out, outdir)
    #cf_hist(med_df, outdir, cutoff = .01)
    print('done')