#!/bin/python 
"""script to make difference plot on time-normalized runs"""
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns 
import sys
import os
import numpy as np

def diff_scatter(data, outdir):
    """plot scatter plot of blood cf - tissue cf"""
    print('plotting diff scatter')
    for cutoff in [0,.01,.1,.2,.3]:
        visible = data[(data['blood'] > cutoff) | (data['tissue'] > cutoff)]
        sns.scatterplot(data = visible, x = 'norm_t_binned',y = 'diff',hue = 'n_drivers')
        plt.show(block = False)
        plt.savefig(os.path.join(outdir, f'diff_cutoff_{cutoff}.png'))
        plt.close()

def diff_box(data, outdir):
    print('plotting diff box')
    for cutoff in [0,.01,.1,.2,.3]:
        visible = data[(data['blood'] > cutoff) | (data['tissue'] > cutoff)]
        sns.boxplot(data = visible[visible['diff']>0], x = 'norm_t_binned',y = 'diff')
        sns.boxplot(data = visible[visible['diff']<0],x = 'norm_t_binned',y = 'diff')
        plt.xticks(rotation = 45)
        plt.xlabel('normalized time')
        plt.ylabel('blood - tissue')
        plt.show(block = False)
        plt.savefig(os.path.join(outdir, f'diffbox_cutoff_{cutoff}.png'))
        plt.close()
def diff_age(data, outdir):
    print('plotting diff_age')
    palette = 'viridis_r'
    edgecolor = 'white'
    linewidth = 3
    s = 40
    for cutoff in [0,.01,.05, .1,.2,.3]:
        visible = data[(data['blood'] > cutoff) | (data['tissue'] > cutoff)]
        #visible['norm_t_binned_scaled'] = 20*visible['norm_t_binned']
        sns.scatterplot(data = visible[visible['diff']>0], x = 
        'norm_t_binned',y = 'diff', hue = 'norm_age', hue_norm = (0,1), 
        legend = False, alpha = .9, palette= palette, edgecolor = edgecolor, s = s)
        sns.scatterplot(data = visible[visible['diff']<0],x = 
        'norm_t_binned',y = 'diff', hue = 'norm_age', hue_norm = (0,1), 
        alpha = .9,palette= palette, edgecolor = edgecolor, s = s)
        sns.lineplot(data = visible[visible['diff']>0], x = 
        'norm_t_binned',y = 'diff', err_style = 'band', color = 'C0', alpha = .5,linewidth = linewidth)
        sns.lineplot(data = visible[visible['diff']<0],x = 
        'norm_t_binned',y = 'diff',err_style = 'band', color = 'C0', alpha = .5,linewidth = linewidth)
        #plt.xticks(np.arange(20), np.linspace(0,1,20),rotation = 45)
        plt.ylim([-.5,.5])
        plt.xlabel('normalized time')
        plt.ylabel('blood - tissue')
        plt.title('Dominant Clone Distortion')
        plt.legend(title = 'normalized age')
        plt.show(block = False)
        plt.savefig(os.path.join(outdir, f'diff_age_cutoff_{cutoff}.png'))
        plt.close()
def diff_cent(data, outdir):
    print('plotting diff_cent')
    palette = 'viridis_r'
    edgecolor = 'white'
    linewidth = 3
    s = 40
    
    for cutoff in [0,.01,.05, .1,.2,.3]:
        visible = data[(data['blood'] > cutoff) | (data['tissue'] > cutoff)]
        #visible['norm_t_binned_scaled'] = 20*visible['norm_t_binned']
        sns.scatterplot(data = visible[visible['diff']>0], x = 
        'norm_t_binned',y = 'diff', hue = 'centroid_r', 
        legend = False, alpha = .9, palette= palette, edgecolor = edgecolor, s = s)
        sns.scatterplot(data = visible[visible['diff']<0],x = 
        'norm_t_binned',y = 'diff', hue = 'centroid_r', 
        alpha = .9,palette= palette, edgecolor = edgecolor, s = s)
        sns.lineplot(data = visible[visible['diff']>0], x = 
        'norm_t_binned',y = 'diff', err_style = 'band', color = 'C0', alpha = .5,linewidth = linewidth)
        sns.lineplot(data = visible[visible['diff']<0],x = 
        'norm_t_binned',y = 'diff',err_style = 'band', color = 'C0', alpha = .5,linewidth = linewidth)
        #plt.xticks(np.arange(20), np.linspace(0,1,20),rotation = 45)
        plt.ylim([-.5,.5])
        plt.xlabel('normalized time')
        plt.ylabel('blood - tissue')
        plt.title('Dominant Clone Distortion')
        plt.legend(title = 'clone dispersal distance')
        plt.show(block = False)
        plt.savefig(os.path.join(outdir, f'diff_cent_cutoff_{cutoff}.png'))
        plt.close()
def logprob_scatter(data, outdir):
    print('plotting logprob scatter')
    for cutoff in [0,.01,.1,.2,.3]:
        visible = data[(data['blood'] > cutoff) | (data['tissue'] > cutoff)]
        sns.scatterplot(data = visible[visible['logprob']>0], x = 'norm_t_binned',y = 'logprob',hue = 'n_drivers')
        sns.scatterplot(data = visible[visible['logprob']<0],x = 'norm_t_binned',y = 'logprob',hue = 'n_drivers')
        plt.xticks(rotation = 45)
        plt.xlabel('normalized time')
        plt.ylabel('log(blood) - log(tissue)')
        plt.savefig(os.path.join(outdir, f'logprob_cutoff_{cutoff}.png'))
        plt.show(block = False)
        plt.close()

def pcterr_scatter(data, outdir):
    print('plotting pcterr scatter')
    for cutoff in [0,.01,.1,.2,.3]:
        visible = data[(data['blood'] > cutoff) | (data['tissue'] > cutoff)]
        sns.scatterplot(data = visible[visible['pcterr']>0], x = 'norm_t_binned',y = 'pcterr',
        size = 'tissue', hue = 'norm_age',hue_norm = (0,1), legend = False, alpha = .5)
        sns.scatterplot(data = visible[visible['pcterr']<0],x = 'norm_t_binned',y = 'pcterr',
        size = 'tissue', hue = 'norm_age',hue_norm = (0,1), alpha = .5)
        plt.xticks(rotation = 45)
        plt.xlabel('normalized time')
        plt.ylabel('pcterr')
        plt.savefig(os.path.join(outdir, f'pcterr_cutoff_{cutoff}.png'))
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
    data = med_df
    #get the weighted average of 
    #data['diff_wtd'] = data['diff']*data['tissue']/data['tissue'].sum()
    #data['norm_age_wtd'] = data['norm_age']*data['tissue']/data['tissue'].sum() 

    #diff_scatter(data, outdir)
    #diff_box(data, outdir)
    #logprob_scatter(data, outdir)
    #pcterr_scatter(data, outdir)
    diff_age(data, outdir)
    diff_cent(data, outdir)
    print('done')


    
    
    
    
    
