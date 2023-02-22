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
from scipy.stats import wilcoxon 
from scipy.stats import ks_2samp
import os
#read file 

DATAPATH = sys.argv[-1]
print('loading data')
#setting globals
#FOLDER = '/home/trachman/sim2/hybrid_bdg/netneg/birthsel'
#PATH = f'{FOLDER}/rep0.csv'
RADIUS = 40 #radius of inner region (get from experiment parameter file)
CUTOFF = .9 #cutoff for fraction of a clone over time that must live in outer region for it to count as a range expander
VAF_CUTOFF = .01

OUTDIR = 'analysis'
SAVEFIGFOLDER = os.path.join(OUTDIR,'figs')
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
Path(SAVEFIGFOLDER).mkdir(parents=True, exist_ok=True)

#titlestring = 'expon. growth death-based selection'
TITLESTRING = 'bdg moving edge birth-based selection'
MIN_FREQUENCY = 0
N_CLONES_TO_PLOT = 100

timedata = pd.read_csv(DATAPATH)
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
comparison['logprob'] = np.log10((1e-10+comparison['blood'])/(1e-10+comparison['tissue']))
comp_reset = comparison.reset_index()
comp_reset['cell_ID'] = timedata['cell_ID']
comp_reset['cell_hge'] = timedata['cell_hge']
comp_reset['is_outer'] = timedata['is_outer']
comp_reset['n_drivers'] = timedata['n_drivers']
comp_reset['death_rate'] = timedata['death_rate']
max_freqs = comp_reset.groupby(['rep','genotype'])['blood'].max()
top_gens = max_freqs.sort_values()[-N_CLONES_TO_PLOT:]

clone_outer_fraction = (timedata.groupby(['rep','genotype'])['is_outer'].sum()/timedata.groupby(['rep','genotype'])['is_outer'].count()).reset_index()
clone_outer_fraction['range_expander'] = clone_outer_fraction['is_outer'] > CUTOFF

scaled_time = comp_reset.groupby('rep')['t'].transform(lambda x: x/x.max())
binned_time = pd.cut(scaled_time, bins = 20).apply(lambda x: np.round(x.mid,2))
comp_reset['norm_t'] = scaled_time
comp_reset['norm_t_binned'] = binned_time.astype(float)


#wilcoxon signed-rank test
print('computing wilcoxon test')
wilcox= comp_reset.groupby(['rep','t']).apply(lambda x: wilcoxon(x['blood'], x['tissue'],zero_method = "pratt")[1]).reset_index()
wilcox['binned'] = binned_time
wilcox['-log10p'] = -np.log10(wilcox[0])

sns.boxplot(data = wilcox, x = 'binned', y ='-log10p')
plt.hlines(y = -np.log10(.05), xmin = 0, xmax = 20, linestyles = ['--'],colors = ['r'],label = 'p = .05')
plt.xticks(rotation = 45)
plt.savefig(f'{SAVEFIGFOLDER}/wilcoxon.png')
plt.title(TITLESTRING)
plt.legend()
plt.show(block  = False)
plt.close()
#kolmogorov-smirnov 2 sample test 
print('computing ks test')
ks = comp_reset.groupby(['rep','t']).apply(lambda x: ks_2samp(x['blood'], x['tissue'])[1]).reset_index()
ks['binned'] = binned_time
ks['-log10p'] = -np.log10(ks[0])

sns.boxplot(data = ks, x = 'binned', y ='-log10p')
plt.hlines(y = -np.log10(.05), xmin = 0, xmax = 20, linestyles = ['--'],colors = ['r'],label = 'p = .05')
plt.xticks(rotation = 45)
plt.savefig(f'{SAVEFIGFOLDER}/ks_2samp.png')
plt.title(TITLESTRING)
plt.legend()
plt.show(block  = False)
plt.close()


#plot clone frequencies for each replicate
print('plotting clone frequencies for each replicate')
replist = comp_reset['rep'].unique()
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
    plt.close()
print('done')
################################################
print('making blood tissue scatterplots...')
tbins = comp_reset['norm_t_binned'].unique()
tbins.sort()
sns.scatterplot(data = comp_reset, x = 'tissue', y = 'blood',hue = 'norm_age')
plt.plot(np.linspace(0,1),np.linspace(0,1))
plt.savefig(f'{SAVEFIGFOLDER}/bloodvtiss_all_time.png')
plt.show(block = False)
for t in tbins:
    data = comp_reset[comp_reset['norm_t_binned']==t]
    sns.scatterplot(data = data, x = 'tissue',y = 'blood', size = 'n_drivers',hue = 'norm_age')
    plt.title(f't = {t}')
    plt.plot(np.linspace(0,1),np.linspace(0,1))
    plt.savefig(f'{SAVEFIGFOLDER}/bloodvtiss_t_{t}.png')
    plt.show(block = False)
print('done')
##########################################
print('making correlation plot...')
for cutoff in [.01,.1,.2,.3]:
    comp_vis = comp_reset[(comp_reset['blood'] > cutoff) | (comp_reset['tissue'] > cutoff)]
    corr = comp_vis.groupby(['norm_t_binned'])[['blood','tissue']].corr().iloc[0::2,-1].reset_index()
    #p = comp_vis.groupby(['norm_t_binned']).apply(lambda x: pearsonr(x['blood'],x['tissue']))
    sns.lineplot(data = corr, x = 'norm_t_binned',y = 'tissue',label = f'cutoff = {cutoff}',linestyle = '--')
    plt.ylabel('pearson r')
    plt.xlabel('norm. t')
plt.title('blood-tissue correlation')
plt.ylim((0,1.1))
plt.savefig(f'{SAVEFIGFOLDER}/corr.png')
plt.show()
print('done')
#save comparison file
comp_reset.to_csv(os.path.join(OUTDIR,'comp.csv'))
#save analysis files
wilcox.to_csv(os.path.join(OUTDIR,'wilcox.csv'))
ks.to_csv(os.path.join(OUTDIR,'ks.csv'))
print('done with everything')
