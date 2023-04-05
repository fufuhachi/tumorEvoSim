import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import sys
import simulation
import os 
def tumor_timeplot(data, label = '', size = None, hue = None, size_norm = None, hue_norm = None, outdir = '', dim = 150, cmap = None):
    t = data.name
    print(f'plotting t = {t}')
    
    fig, ax = plt.subplots()
    ax = sns.scatterplot(data = data,x = 'x',y = 'y', size = size, hue = hue,
    size_norm = size_norm, hue_norm = hue_norm, legend = False, palette = cmap,s =1, marker = 's')
    #simulation.tumor_scatter(x = data['x'], y = data['y'], c = data['genotype'], cmap_name = 'viridis', dim = 90)
    plt.title(f't = {t}')
    ax.set_xlim((-dim, dim))
    ax.set_ylim((-dim, dim))
    ax.set_box_aspect(1)
    plt.savefig(os.path.join(outdir, f'{label}_t_{t}.png'))
    plt.show(block = False)
    plt.close()

def timeplot2(data, label, vmin, vmax, xmin,xmax, ymin, ymax, cmap = None):
    t = data.name
    print(f'plotting t = {t}')
    pivoted = data.pivot(index = 'x', columns = 'y', values = 'genotype')
    sns.heatmap(pivoted,vmin = vmin, vmax = vmax, cmap = cmap, square = True, cbar = False)
    plt.title(f't = {t}')
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.savefig(os.path.join(outdir, f'{label}_t_{t}.png'))
    plt.show(block = False)
    plt.close()
def clone_timeplot(clones, colors):
    """accepts a matrix of time-series data for blood and tissue clone fractions
    and a list of the colors used for each clone"""
if __name__ == "__main__":
    data = pd.read_csv(sys.argv[1])
    label = sys.argv[2] if len(sys.argv) > 2 else 'tumor-timeplot'
    outdir = sys.argv[3] if len(sys.argv) > 3 else ''
    #scale the time 
    data['norm_t'] = data['t']/data['t'].max() 
    binned_time = pd.cut(data['norm_t'],60).apply(lambda x: np.round(x.mid,2))
    data['norm_t_binned'] = binned_time.astype(float)
    med_vals = data.groupby(['norm_t_binned','rep'])['t'].median().reset_index()[['rep','t']] #reset_index to convert groupby object to DataFrame
    # get the index of median column in the original dataframe
    med_df = pd.merge(data, med_vals, on  = ['rep','t'], how = 'inner')
    data = med_df


    hgemax = data['cell_hge'].max()
    data['gen_categ'] = data['genotype']#.astype(str)
    genmax = data['genotype'].max()
    hue_norm=(0, genmax)
    dim = data[['x','y']].abs().values.flatten().max()+2
    vals = np.linspace(0,1,data['genotype'].max())
    rng = np.random.default_rng(seed = 1234)
    rng.shuffle(vals)
    cmap = plt.cm.colors.ListedColormap(plt.cm.tab20(vals))
    xmin = data['x'].min()
    xmax = data['x'].max()
    ymin = data['y'].min()
    ymax = data['y'].max()
    #data.groupby(['norm_t']).apply(lambda x: tumor_timeplot(x, outdir = outdir, label = label, hue_norm = hue_norm, hue = 'genotype', dim = dim, cmap = cmap))
    #data.groupby(['norm_t']).apply(lambda x: timeplot2(x, label = label, cmap = cmap, vmin = 0, vmax = genmax, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
    #use same cmap to make clone plot
    #use only clones that occupy at least 10% of tumor at any one point 

    data['above_cutoff'] = (data['blood'] > .1) | (data['tissue'] > .1)
    included = data.groupby(['rep','genotype'])['above_cutoff'].transform(any)
    to_plot = data[included]
    clones = to_plot['genotype'].unique().sort()
    print(clones)