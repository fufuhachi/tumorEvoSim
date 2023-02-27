#!/bin/python 
"""script to make difference plot on time-normalized runs"""
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns 
import sys
import os

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
        sns.scatterplot(data = visible[visible['pcterr']>=0], x = 'norm_t_binned',y = 'pcterr',size = 'tissue',legend = False, alpha = .5)
        sns.scatterplot(data = visible[visible['pcterr']<=0],x = 'norm_t_binned',y = 'pcterr',size = 'tissue',alpha = .5)
        plt.xticks(rotation = 45)
        plt.xlabel('normalized time')
        plt.ylabel('pcterr')
        plt.savefig(os.path.join(outdir, f'pcterr_cutoff_{cutoff}.png'))
        plt.show(block = False)
        plt.close()

if __name__ == '__main__':
    print('loading data')
    data = pd.read_csv(sys.argv[1])
    outdir = sys.argv[2] if len(sys.argv) > 2 else ''
    diff_scatter(data, outdir)
    diff_box(data, outdir)
    logprob_scatter(data, outdir)
    pcterr_scatter(data, outdir)
    print('done')


    
    
    
    
    
