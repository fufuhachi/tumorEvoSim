"""Script to take in dataframe with summary of tumor information for multiple replicates at multiple timepoints. 
See simulation.tumor_summary for details on the format of the csv. 
    Inputs: a file string for the csv 
    Outputs: plots aggregated by replicate measuring the distortion of clonal fraction over time in ctDNA compared to tissue. The estimates of ctDNA fraction are based on the crude assumption 
    that the clone frequencies in the blood are the weighted average of the tissue frequencies, scaled by the cell death rate. """
import pandas as pd
import numpy as np
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import simulation
import ast
import json
#read file 
datapath = sys.argv[-1]

print('loading data')

#setting globals
FOLDER = '/home/trachman/sim2/hybrid_bdg/netneg/birthsel'
PATH = f'{FOLDER}/rep0.csv'
RADIUS = 40 #radius of inner region (get from experiment parameter file)
CUTOFF = .9 #cutoff for fraction of a clone over time that must live in outer region for it to count as a range expander
VAF_CUTOFF = .01
SAVEFIGFOLDER = f'../tumorEvoSim/{FOLDER}/figs'
Path(SAVEFIGFOLDER).mkdir(parents=True, exist_ok=True)
#titlestring = 'expon. growth death-based selection'
TITLESTRING = 'bdg moving edge birth-based selection'
MIN_FREQUENCY = 0
N_CLONES_TO_PLOT = 100

timedata = pd.read_csv(datapath)
timedata['cell_hge'] = timedata['death_rate']/(timedata['birth_rate']-timedata['death_rate']+33)
timedata['is_outer'] = timedata['r'] > RADIUS
outer_count = timedata.groupby(['rep','t'])['is_outer'].sum()
pop_size = timedata.groupby(['rep','t'])['is_outer'].count()
pop_size = pop_size.reset_index()
per_time_count = outer_count.reset_index() #get fraction outside of boundary per time as column 'is_outer' 
per_time_count['pop_size'] = pop_size['is_outer']

#make comparison file
print('making comparison files')
tissue = timedata.groupby(['rep','t','genotype'])['genotype'].count()/timedata.groupby('t')['genotype'].size()
blood = timedata.groupby(['rep','t','genotype'])['cell_hge'].sum()/timedata.groupby(['t'])['cell_hge'].sum()

comparison = tissue.compare(blood)
comparison.columns = ['tissue','blood']
comparison['diff'] = (comparison['blood'] - comparison['tissue'])
comparison['pcterr'] = 100*(comparison['blood'] - comparison['tissue'])/comparison['tissue']
comp_reset = comparison.reset_index()
max_freqs = comp_reset.groupby(['rep','genotype'])['blood'].max()
top_gens = max_freqs.sort_values()[-N_CLONES_TO_PLOT:]

clone_outer_fraction = (timedata.groupby(['rep','genotype'])['is_outer'].sum()/timedata.groupby(['rep','genotype'])['is_outer'].count()).reset_index()
clone_outer_fraction['range_expander'] = clone_outer_fraction['is_outer'] > CUTOFF

#plot clone frequencies for each replicate
print('plotting clone frequencies for each replicate')
replist = list(set(timedata['rep']))
replist.sort()
for rep in replist:
    current = comparison.reset_index() #get comparison as 2d dataframe
    current = current[current['rep']==rep] #slice by replicate
    comp_reset_rep = current.reset_index() 
    max_freqs = comp_reset_rep.groupby('genotype')['blood'].max()
    top_gens = max_freqs.sort_values()[-N_CLONES_TO_PLOT:]
    outer_count = per_time_count[per_time_count['rep']==rep] #
    range_expanders = clone_outer_fraction[clone_outer_fraction['rep']==rep]

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax2.set_ylabel('population size')
    ax1.set_ylabel('clonal fraction')
    for gen in top_gens.index:
        data = comp_reset_rep[comp_reset_rep['genotype']==gen]
        is_expander = range_expanders[range_expanders['genotype'] == gen]['range_expander'].values
        time = data['t']
        #plot fraction difference 
        #ax1.plot(data['t'],data['diff'], color = 'green' if is_expander else 'blue', label = 'expanding clone' if is_expander else 'non-expanding clone')
        line,  = ax1.plot(time, data['blood'], label = f'clone # {gen} (blood)',linestyle = '--')
        ax1.plot(time, data['tissue'], label = f'clone # {gen} (tissue)', color = line.get_color())
    ax2.plot(outer_count['t'], outer_count['pop_size'], color = 'black', label = 'pop. size')
    plt.title(f'rep {rep}')
    plt.legend()
    ax1.set_xlabel('time')
    plt.savefig(f'{SAVEFIGFOLDER}/timeplot_rep_{rep}.png')
    plt.show(block = False)

def plot_wrapper(toplot, label):
    simulation.tumor_scatter(toplot['x'],toplot['y'], toplot['n_drivers'], dim = 160, show = False)
    plt.title(label)
    plt.colorbar()
    plt.savefig(f'{SAVEFIGFOLDER}/driverPlot_{label}.png')
    plt.show(block = False)

timedata[timedata['rep']==0].groupby(['rep','t']).apply(lambda x: plot_wrapper(x, f't = {x.name[1]:.2f}'))

diffsum = comp_reset.groupby(['rep','genotype'])['diff'].sum().reset_index()
clone_outer_fraction['diffsum'] = diffsum['diff']
marginal_by_genotype = clone_outer_fraction.groupby(['rep', 'range_expander'])['diffsum'].sum().reset_index()
sns.boxplot(data = marginal_by_genotype,x = 'range_expander', y = 'diffsum' )
plt.ylabel('cumulative difference in blood')
plt.xlabel('clone location')
plt.xticks([0,1], ['inner region', 'outer region'])
plt.title(TITLESTRING)
plt.savefig(f'{SAVEFIGFOLDER}/diffsum.png')
plt.ylim([-21,21])
plt.show(block = False)