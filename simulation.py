#code for simulation
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import classes
import time
from scipy.stats import binned_statistic
import argparse 
import sys
#simulation code
#simulation overview (following Chkaidze et al. 2019):
#set up 2d or 3d grid

def stopping_condition(tumor):
    nmax = classes.params['n_cells']
    imax = classes.params['max_iter']
    return tumor.N == nmax or tumor.iter > imax or tumor.hit_bound
def run(tumor = None):
    t = []
    n = []
    if tumor is None:
        tumor = classes.Tumor()
    while not stopping_condition(tumor):
        tumor.iterate()
        if tumor.iter%40000 ==0:
            print(f'size = {tumor.N} at {int(tumor.t/365)} years {tumor.t%365} days')

    return tumor
"""save tumor snapshots and further analysis"""
def save_tumor(tumor):
    #TODO
    pass
"""convert int to rbg tuple"""
def int_to_rgb(RGBint):
    Blue =  RGBint & 255
    Green = (RGBint >> 8) & 255
    Red =   (RGBint >> 16) & 255
    return (Red, Blue, Green)
def plot_tumor(tumor,drivers=False):
    graph = tumor.graph.copy()
    pos = np.array([np.array(item.pos) for item in tumor.cells.items])
    pos = tuple([tuple(row) for row in pos.T])
    if drivers:
        gens = [item.gen.n_drivers for item in tumor.cells.items]
    else:
        gens = [item.gen.ID for item in tumor.cells.items]
    graph[pos] = gens
    if classes.params['dim'] ==2:
        graph[graph==0]=-np.max(graph.flatten())
        sns.heatmap(graph)
        plt.show()
    else:
        print('not yet implemented!')
        raise(NotImplementedError)
def plot_slice(tumor, ax = 0):
    if classes.params['dim']==2:
        plot_tumor(tumor)
    else:
        graph = tumor.graph.copy()
        pos = np.array([np.array(item.pos) for item in tumor.cells.items])
        pos = tuple([tuple(row) for row in pos.T])
        gens = [item.gen.ID for item in tumor.cells.items]
        graph[pos] = gens
        graph[graph==0] = -np.max(graph.flatten())
        if ax=='x' or ax==0:
            sns.heatmap(graph[tumor.center[0]])
        elif ax=='y' or ax==1:
            sns.heatmap(graph[:,tumor.center[0],:])
        elif ax=='z' or ax==2:
            sns.heatmap(graph[:,:,tumor.center[0]])
        else:
            print('invalid axis')
        plt.title(f'Slice along {ax} axis')
        plt.show()
        return graph

"""Given a tumor genotype ListDict, construct tree
    This tree includes possibly extinct genotypes and represents the true evolutionary history 
    of the simulation 
"""
def get_tumor_phylogeny(gens):
    
    root = gens.get_item(0)
    return make_newick(root)
    
"""Given a list of genotypes and root, return newick tree"""
def make_newick(root):
    #base case
    if len(root.children)==0:
        return f'{root.ID}'
    else:
        return f'{tuple([make_newick(child) for child in root.children])}{root.ID}'
"""given a list of genotypes, return genotypes sorted by pop and drivers"""
def get_genotype_dist(gens):
    genlist = gens.items
    df = pd.DataFrame([[g.ID,g.number,g.n_drivers] for g in genlist])
    df.columns = ['clone','number','n_drivers']
    bypop = df.sort_values(by = 'number')
    bydriv = df.sort_values(by='n_drivers')
    return bypop,bydriv
"""plot top <topn> genotypes by population and numebr of driver mutations, return correlation between groups """
def plot_genotype_dist(bypop, bydriv, topn = 10):
    bypop.iloc[-topn:].plot.bar(x = 'clone', y='number')
   
    plt.show()
    bypop.iloc[-topn:].plot.bar(x='clone', y = 'n_drivers')
    plt.show()
    plt.plot(bypop.iloc[-topn:]['clone'],(bypop['number'].sum()/bypop.iloc[-topn:]['number'])**2,'.')
    plt.show()
def plot_growth(tumor):
    plt.scatter(tumor.t_traj, tumor.N_traj)
    plt.xlabel('time (days)')
    plt.ylabel('size (cells)')
    plt.show()
"""function to save entire tumor object and all related objects to folder"""
def save_tumor(tumor):
    #TODO
    pass
"""function to read tumor object from folder so that simulation can be started"""
def load_tumor(tumor):
    #TODO
    pass
"""Take bulk sample at specified location
inputs:
    tumor: tumor object to be sampled
    pos: position, 2d or 3d tuple defining center of cube to be sampled. must match dimensions of tumor
    length: side length of cube to be sampled 
    depth: read depth
    cutoff: minimal detectable VAF 
outputs:
    vector of VAFs 
"""
def bulk_sample(tumor, pos: tuple, length: int, depth: int = 0, cutoff: float=0):
    #TODO 
    try:
        assert(len(tumor.graph.shape)==len(pos))
    except(AssertionError):
        print(f'expected dimension {tumor.graph.shape} but got {len(pos)}')
        print('exiting...')
        sys.exit()
    try:
        assert(length%2==1)
    except(AssertionError):
        print('length must be odd')
        print('exiting...')
        sys.exit()
    genmat = get_gen_graph(tumor)
    l = (length-1)/2
    try:
        sample = genmat[pos[0]-l:pos[0]+l,pos[1]-l:pos[1]+l]
        if len(pos)==3:
            sample = sample[:,:,pos[2]-l:pos[2]+l]
    except(IndexError):
        print('sample out of bounds')
        print('exiting...')
        sys.exit()
    #get frequencies of genotypes
    gens, counts = np.unique(sample.flatten(),return_counts = True)
    

    
"""return lattice populated by genotypes"""
def get_gen_graph(tumor):
    graph = tumor.graph.copy()
    pos = np.array([np.array(item.pos) for item in tumor.cells.items])
    pos = tuple([tuple(row) for row in pos.T])
    gens = [item.gen.ID for item in tumor.cells.items]
    graph[pos] = gens
    return graph


   
if __name__=="__main__":
    t0 = time.time()
    out = run()
    t1 = time.time()
    print(t1-t0)
    #plot_tumor(out)
    gens = out.gens
    pop, dr = get_genotype_dist(gens)
    plot_genotype_dist(pop,dr)

        
