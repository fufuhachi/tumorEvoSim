#!/bin/python 
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import sys
import os
import numpy as np
def plot_corr_all(data, outdir = ''):
    """plot correlation of concatenated runs """
    print('plotting corr with all data')
    for cutoff in [.01,.1,.2,.3]:
        vis = data[(data['blood'] > cutoff) | (data['tissue'] > cutoff)]
        corr = vis.groupby(['norm_t_binned'])[['blood','tissue']].corr().iloc[0::2,-1].reset_index()
        #p = comp_vis.groupby(['norm_t_binned']).apply(lambda x: pearsonr(x['blood'],x['tissue']))
        sns.lineplot(data = corr, x = 'norm_t_binned',y = 'tissue',label = f'cutoff = {cutoff}',linestyle = '--')
        plt.ylabel('pearson r')
        plt.xlabel('norm. t')
    plt.title('blood-tissue correlation')
    #plt.ylim((0,1.1))
    plt.savefig(os.path.join(outdir, 'corr.png'))
    plt.show(block = False)
    plt.close()
    print('done')

def plot_corr_posneg(data, outdir = ''):
    """plot corelation of concatenated runs split by negative and positive cf difference (blood - tissue)"""
    print('plotting corr with pos and neg diff separately')
    fig, ax = plt.subplots()
    for cutoff in [.01,.1,.2,.3]:
        vis = data[(data['blood'] > cutoff) | (data['tissue'] > cutoff)]
        corrpos = vis[vis['diff']>=0].groupby(['norm_t_binned'])[['blood','tissue']].corr().iloc[0::2,-1].reset_index()
        corrneg = vis[vis['diff']<=0].groupby(['norm_t_binned'])[['blood','tissue']].corr().iloc[0::2,-1].reset_index()
        #p = comp_vis.groupby(['norm_t_binned']).apply(lambda x: pearsonr(x['blood'],x['tissue']))
        sns.lineplot(data = corrpos, x = 'norm_t_binned',y = 'tissue',label = f'cutoff = {cutoff} +',linestyle = '-', ax = ax)
        sns.lineplot(data = corrneg, x = 'norm_t_binned',y = 'tissue',label = f'cutoff = {cutoff} -',linestyle = '--',color = ax.lines[-1].get_color(),ax = ax) 
    plt.title('blood-tissue correlation')
    #plt.ylim((0,1.1))
    plt.ylabel('pearson r')
    plt.xlabel('norm. t')
    plt.savefig(os.path.join(outdir, 'corr-posneg.png'))
    plt.show(block = False)


    plt.close()
    print('done')
def plot_bloodvtiss_scatterplots(data, outdir):
    print('making blood tissue scatterplots...')
    tbins = data['norm_t_binned'].unique()
    tbins.sort()
    sns.scatterplot(data = data, x = 'tissue', y = 'blood',hue = 'norm_age')
    plt.plot(np.linspace(0,1),np.linspace(0,1))
    plt.savefig(os.path.join(outdir,'bloodvtiss_all_time.png'))
    plt.show(block = False)
    plt.close()
    for t in tbins:
        data = data[data['norm_t_binned']==t]
        sns.scatterplot(data = data, x = 'tissue',y = 'blood', size = 'n_drivers',hue = 'norm_age')
        plt.title(f't = {t}')
        plt.plot(np.linspace(0,1),np.linspace(0,1))
        plt.savefig(os.path.join(outdir,f'bloodvtiss_t_{t}.png'))
        plt.show(block = False)
        plt.close()
print('done')
if __name__ == '__main__':
    data = pd.read_csv(sys.argv[1])
    outdir = sys.argv[2] if len(sys.argv) > 2 else ''
    plot_corr_posneg(data,outdir)
    plot_bloodvtiss_scatterplots(data, outdir)
    