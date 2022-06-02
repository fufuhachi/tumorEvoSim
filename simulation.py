#code for simulation
from multiprocessing.dummy import Array
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
    maxxed_cells = tumor.N == nmax
    maxxed_iter = tumor.iter > imax
    return maxxed_cells or maxxed_iter or tumor.hit_bound
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
def plot_tumor(tumor,drivers=False,trim = 0):
    graph = tumor.graph.copy()
    pos = np.array([np.array(item.pos) for item in tumor.cells.items])
    pos = tuple([tuple(row) for row in pos.T])
    if drivers:
        gens = [item.gen.n_drivers for item in tumor.cells.items]
    else:
        gens = [item.gen.ID for item in tumor.cells.items]
    graph[pos] = gens
    if len(tumor.graph.shape) ==2:
        graph = graph[trim:-(1+trim),trim:-(1+trim)]
        graph[graph==0]=-np.max(graph.flatten())
        sns.heatmap(graph,cbar= False,square = True)
        plt.show()
    else:
        print('not yet implemented!')
        raise(NotImplementedError)
def plot_slice(tumor, ax = 0,trim=0):
    if len(tumor.graph.shape)==2:
        plot_tumor(tumor,trim)
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
def genotype_plot(bypop, bydriv, topn = 10):
    bypop.iloc[-topn:].plot.bar(x = 'clone', y='number')
   
    plt.show()
    bypop.iloc[-topn:].plot.bar(x='clone', y = 'n_drivers')
    plt.show()
    plt.plot(bypop.iloc[-topn:]['clone'],(bypop['number'].sum()/bypop.iloc[-topn:]['number'])**2,'.')
    plt.show()
def plot_genotype_dist(**kwargs):
    bypop, bydriv = get_genotype_dist(gens)
    genotype_plot(bypop = bypop, bydriv = bydriv, topn=10)
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
def bulk_sample(tumor, pos: tuple, length = None, depth: int = 0):
    #TODO 
    try:
        assert(len(tumor.graph.shape)==len(pos))
    except(AssertionError):
        print(f'expected dimension {tumor.graph.shape} but got {len(pos)}')
        print('exiting...')
        sys.exit()
   
    genmat = get_gen_graph(tumor)
    if length is not None: 
        try:
            l = int((length-1)/2)
            assert(length%2==1)
        except(AssertionError):
            print('length must be odd')
            print('exiting...')
            sys.exit()
        try:
            sample = genmat[pos[0]-l:pos[0]+l+1,pos[1]-l:pos[1]+l+1]
            
            if len(pos)==3:
                sample = sample[:,:,pos[2]-l:pos[2]+l+1]
            assert(sample.shape == (length, length))
            
        except(IndexError):
            print('sample out of bounds')
            print('exiting...')
            sys.exit()
        except(AssertionError):
            print(f'dimensions of sample are wrong, expected ___ but got {sample.shape}')
            print('exiting...')
    else:
        sample = genmat[genmat>0] #sample whole tumor
    #get frequencies of genotypes
    #gens, counts = np.unique(sample.flatten(),return_counts = True)

   # for g in sample.flatten():
   #     if g==1:
   #         tumor.gens.get_item(g).neut = np.array([-1])
   #     if g==2: 
   #         print(tumor.gens.get_item(g).neut)
   #     all_muts = np.append(all_muts, tumor.gens.get_item(g).neut)
   #     all_muts = np.append(all_muts, tumor.gens.get_item(g).drivers)
    
    all_muts = map(lambda genid: get_muts(genid, tumor), sample.flatten())

    all_muts = np.hstack(np.array(list(all_muts),dtype=object)).flatten()
    muts, counts = np.unique(all_muts,return_counts = True)
    vafs = counts/counts.sum() if counts.sum() > 0 else counts
    #turn gens into VAFs
    #get all mutations in set, 
    
    if depth <1:
        m, f = muts, vafs
    else:
        #simulate sequence reads
        m, f = sequence(muts, vafs, depth)
    return m, f
def genotype_sample(tumor, pos, length):
    try:
        assert(len(tumor.graph.shape)==len(pos))
    except(AssertionError):
        print(f'expected dimension {tumor.graph.shape} but got {len(pos)}')
        print('exiting...')
        sys.exit()
   
    genmat = get_gen_graph(tumor)
    if length is not None: 
        try:
            l = int((length-1)/2)
            assert(length%2==1)
        except(AssertionError):
            print('length must be odd')
            print('exiting...')
            sys.exit()
        try:
            sample = genmat[pos[0]-l:pos[0]+l+1,pos[1]-l:pos[1]+l+1]
            
            if len(pos)==3:
                sample = sample[:,:,pos[2]-l:pos[2]+l+1]
            assert(sample.shape == (length, length))
            
        except(IndexError):
            print('sample out of bounds')
            print('exiting...')
            sys.exit()
        except(AssertionError):
            print(f'dimensions of sample are wrong, expected ___ but got {sample.shape}')
            print('exiting...')
    else:
        sample = genmat[genmat>0]
    print(sample)
    gens, counts = np.unique(sample.flatten(),return_counts = True)
    return gens, counts/counts.sum()
def get_muts(genID,tumor):
    if genID==0:
        return np.array([])
    gen = tumor.gens.get_item(genID)
    #if genID==1:
       #neut = np.array([-1])
    return np.append(gen.neut, gen.drivers)
"""Given a list of mutations and their frequencies, simulate bulk sequencing in the manner of
Chkhaidze et al: 
Get coverage with distribution pois(depth) for each mutation
Get get frequency as binom(vaf, coverage)
return only those greater than cutoff
"""
def sequence(muts, vafs,depth):
    coverage = np.random.poisson(depth, muts.shape[0])
    sampled_vafs = np.array([np.random.binomial(cov, f) for f, cov in zip(vafs, coverage)])
    sampled_vafs = sampled_vafs/coverage
    return muts, sampled_vafs
    
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

        
