
#from ast import Index

import numpy as np

#import pandas as pd
import pickle
#import sys
#import time

#classes and object methods
class Simulation():
    def __init__(self, params):
        self.params = params
        self.graph = set_graph(params)
        self.tumor = Tumor(self)
        self.t_traj = []
        self.N_traj = []
    def stopping_condition(self):
        nmax = self.params['n_cells']
        imax = self.params['max_iter']
        return self.tumor.N == nmax or self.tumor.iter > imax or self.tumor.hit_bound or self.tumor.N ==0
    
    def run(self,rep=0):
        prev = 0
        while not self.stopping_condition():
            self.tumor.iterate()
            
            if self.tumor.iter%40000==0:
                print(f'size = {self.tumor.N} at {int(self.tumor.t/365)} years {self.tumor.t%365} days')
            if self.tumor.N%self.params['save_interval']==0 and self.tumor.N!=prev:
                #print(f'saving with driviers {self.tumor.drivers}') #DEBUG
                save_object(self, f'{self.params["exp_path"]}/rep={rep}_ncells={self.tumor.N}.pkl')
                prev = self.tumor.N #only save same cell size once 
        #save final 
        save_object(self, f'{self.params["exp_path"]}/rep={rep}_ncells={self.tumor.N}.pkl')
        return self.tumor
    
#classes and object methods
        
class Cell():
    """a cell defined by a unique ID with mutable attribute position (optional) and immutable attribute genotype
    gen: Genotype of the cell. 
    """
    def __init__(self,sim,ID, gen, pos) -> None:
        self.sim = sim
        self.ID = ID
        self.pos = pos #nonneg DIM-tuple 
        self.gen = gen #Genotype 
        #self.death_rate = self.sim.params['init_death_rate'] if self.sim.params['fixed_death_rate'] else self.sim.params['init_death_rate']*np.power(self.sim.params['driver_dr_factor'], self.gen.n_drivers)
        
    def get_death_rate(self):
        return self.sim.params['dr_function'](self,**self.sim.params['dr_params'])
        
    """mutate function depends only on tumor-level mutation order (infinite allele assumption)
        Samples poisson-distributed number of neutral mutations, driver mutations, and mutator mutations 
        If total is nonzero, creates a new genotype with new mutations added. Each mutation assigned a unique integer based on the order it appeared. Genotype's parent is current genotype. 
    """
    def default_mutate(self):
            #print('starting mutation')
            #print('current cells are:')
            #print(tumor.cells)
            #sample number of mutations
            n_neuts = np.random.poisson(self.gen.neut_rate)
            n_drivers = np.random.poisson(self.gen.driver_rate)
            n_mutators = np.random.poisson(self.gen.mutator_rate)
            n = self.sim.tumor.next_snp
            if n_neuts + n_drivers + n_mutators >0:
                neutral = None
                drivers = None
                mutators = None
                if n_neuts>0:
                    neutral = np.append(self.gen.neut, np.arange(n,n+n_neuts))
                    n+=n_neuts
                if n_drivers >0:
                    print(f"drivers {np.arange(n,n+n_drivers)} at size {self.sim.tumor.N}")
                    drivers = np.append(self.gen.drivers, np.arange(n,n+n_drivers))
                    n+=n_drivers
                    
                if n_mutators >0:
                    print(f"mutator at size {self.sim.tumor.N}")
                    mutators = np.append(self.gen.mutators, np.arange(n,n+n_mutators))
                    n+=n_mutators
                new_gen = Genotype(sim = self.sim,ID = self.sim.tumor.gens.len()+1, parent = self.gen, neutral = neutral, drivers =drivers, mutators = mutators)
                
                self.gen.children = np.append(self.gen.children, new_gen)
                self.gen = new_gen
                self.sim.tumor.add_genotype(new_gen)
            self.sim.tumor.next_snp = n
            #update max death rate
            
    def progression_mutate(self):
        #TODO
        raise(NotImplementedError)
    def get_empty_nbrs(self):
        nbrs = (self.sim.params['neighbors'] + np.array(self.pos))
        tups = tuple([tuple(col) for col in nbrs.T])
        ids = self.sim.graph[tups]
        return nbrs[ids==0]
    
    def __repr__(self):
        return str(f'Cell# {self.ID} with gen {self.gen.ID} at {self.pos}')
    
        
"""object to store mutations of a particular lineage. Mutations are either passengers or drivers 
Drivers are a subset of mutations that increase the birth rate of the cell 
Attributes:
    id: int representing the order in which the Genotype appeared during the simulation (a genotype id)
    muts: np array of mutations (ints)
    drivers: np array of driver mutations (ints), a subset of muts 
    n_drivers: the number of driver mutations 
    br: the birth rate (a function of n_drivers)
    parent: the genotype id of the parent Genotype 
    children: np int array of child Genotypes 
    number: the number of individuals with this Genotype 

"""
class Genotype():
    def __init__(self, sim, ID=1, parent=None, neutral=None, drivers=None, mutators = None) -> None:
        self.sim = sim
        self.ID = ID
        self.neut = np.array([], dtype = int) if neutral is None else neutral
        self.n_neut = self.neut.shape[0]
        self.drivers = np.array([], dtype = int) if drivers is None else drivers
        self.n_drivers = self.drivers.shape[0]
        self.mutators = np.array([], dtype = int) if mutators is None else mutators
        self.n_mutators = self.mutators.shape[0] #TODO: justify mutator 
        self.neut_rate = self.sim.params['init_mut_prob']*np.power(self.sim.params['mutator_factor'], self.n_mutators)
        self.driver_rate = self.sim.params['driver_prob']*np.power(self.sim.params['mutator_factor'], self.n_mutators)
        self.mutator_rate = self.sim.params['mutator_prob']*np.power(self.sim.params['mutator_factor'], self.n_mutators)
       
        self.parent = parent
        self.children = np.array([],dtype = int)
        self.number = 1
        
        #add genotype to genotype list
    def add(self):
        #self.sim.tumor.add_driver_count(self.n_drivers) #not necessary if death rate is not divided by dmax 
        self.number+=1
    def remove(self):
        if self.number ==0:
            print("call to remove extinct genotype!!")
            raise(Exception)
        else:
            #if self.number==1:
            #    self.sim.tumor.remove_driver_count(self.n_drivers)
            self.number-=1

    def __repr__(self):
        parID = 'none' if self.parent is None else self.parent.ID
        return str(f'Genotype {self.ID} from {parID} with {self.n_drivers} drivers')
    
"""Class containing Tumor object
    attributes:
    graph: the graph (lattice for now) containing genotype ids 
    gens: Array of genotypes in the population 
    cells: ListDict of cells in population 
    center: coordinate of center 
    N: population size (size of cells)
    next_snp: next mutation when progression is None, None otherwise 
    progression: matrix defining mutation tree to be sampled from if defined. Sampling depends on current genotype 
    bmax: maximum population birth rate 
    dmax: maximum population death rate
    t: current time 
    iter: current iteration 
    hit_bound: boolean of whether tumor has reached boundary 
    drivers: total driver mutations in tumor
    mutators: total mutators in tumor 
"""
class Tumor():
    def __init__(self,sim) -> None:
        self.sim = sim
        self.graph = sim.graph
        self.params = sim.params
        coord = int(self.params['boundary']/2)
        self.center = tuple(coord*np.ones(self.params['dim'],dtype = int))
        self.gens = ListDict()
        first_gen = Genotype(self.sim)
        self.gens.add_item(first_gen)
        self.cells = ListDict()
        first_cell = Cell(self.sim,ID = 1,gen = first_gen, pos = self.center)
        self.cells.add_item(first_cell)
        
        self.graph[self.center] = first_cell.ID
        self.N = 1
        if self.params['progression'] is None:
            #assume infinite sites 
            self.next_snp = 0
        else:
            self.next_snp = None
           
        self.progression = self.params['progression']
        self.bmax = self.params['init_birth_rate']
        self.dmax = self.params['max_death_rate']
        
        
        self.hit_bound = False 
        self.drivers = np.array([],dtype = int) 
        self.driver_counts = np.array([], dtype = int) #ascending sorted list of number of drivers in each genotype
        self.max_drivers = 0
        self.mutators = np.array([], dtype = int)
        self.iter = 1
        self.t = 0

    """add a number of drivers: Rationale: Rather than use priority queue to track highest death rate, 
    use fact that death rate is function of number of driver mutations to index in constant time
    """
    def add_driver_count(self,n):
        if n+1 > self.driver_counts.shape[0]:
            print('added driver count')
            self.driver_counts = np.append(self.driver_counts, np.zeros(n+1-self.driver_counts.shape[0]))
            self.driver_counts[-1]=1
        else:
            if self.driver_counts[n-1]==0:
                self.driver_counts[n-1] = 1
                self.dmax = self.get_dmax()
      
    """remove a number of drivers"""
    def remove_driver_count(self,n):
        if n+1 > self.driver_counts.shape[0]:
            print('driver count should have been added')
            raise(IndexError)
        else:
            if self.driver_counts[n-1]==1:
                print('removing driver count')
                self.driver_counts[n-1] = 0
                self.dmax = self.get_dmax()
    def get_dmax(self):
        min_driv =  self.driver_counts[self.driver_counts ==1][0]
        return self.params['init_death_rate']*np.power(self.params['driver_dr_factor'], min_driv)
    """update the time of the simulation after a birth/death event """
    def update_time(self):
        self.t+=1/(self.bmax*self.N)
        self.iter+=1
    """function to perform one iteration of growth"""
    def iterate(self):
        #choose random cell
        #rejection sample birth and death 
        self.rand_birth_death_mutate()
        self.update_time()
        if self.iter%40000 ==0:
            self.sim.t_traj = np.append(self.sim.t_traj, self.t)
            self.sim.N_traj = np.append(self.sim.N_traj, self.N)
        pass
    """Function to simulate birth, death, and mutation. """
    def rand_birth_death_mutate(self):
        #print('starting birth_death_mutate')
        #print(f'current cells: {self.cells}')
        #print('current graph:')
        #print(self)

        cell = self.cells.choose_random_item()
        #print(f'chose {cell}')
        gen = cell.gen
        nbrs = cell.get_empty_nbrs()
        
        if nbrs.shape[0] >0:
            nbr = nbrs[np.random.randint(0,nbrs.shape[0])]
            #print(f'chose nbr {nbr}')
            self.hit_bound = check_bound(nbr, self.params['boundary'])
            #birth rate held constant, always gives birth
            new_cell = Cell(self.sim, ID = self.iter+1, gen = gen, pos = tuple(nbr))
            self.add_cell(new_cell)
            #print(f'cell pop is now: {self.cells}')
            #mutate daughter cell, sample to determine whether to kill parent cell. If no, mutate parent cell 
            new_cell.default_mutate()
           # if np.random.random() <gen.death_rate/self.dmax:
            if np.random.random() < cell.get_death_rate():
                #kill cell
                #print(f'{cell} died')
                self.remove_cell(cell)
            else:
                cell.default_mutate()
                
        return

    """function to add cell to cell list, change relevant fields"""        
    def add_cell(self, born_cell):
        self.graph[born_cell.pos] = born_cell.ID
        self.cells.add_item(born_cell)
        self.N+=1
        born_cell.gen.add()
    """function to remove cell from cell list change relevant fields"""
    def remove_cell(self, doomed_cell):
        #print('removing cell from')
        #print(self.cells)
        #print('and dictionary is')
        #print(self.cells.item_to_position)
        self.graph[doomed_cell.pos]=0
        self.cells.remove_item(doomed_cell)
        if self.N ==0:
            print('tumor pop is empty')
            raise(Exception)
        else:
            self.N-=1
        doomed_cell.gen.remove()
        #print('cells are now')
        #print(self.cells)
        #print('and dictionary is')
        #print(self.cells.item_to_position)

    def add_genotype(self, gen):
        self.gens.add_item(gen)
        self.drivers = np.append(self.drivers, gen.drivers)
        self.mutators = np.append(self.mutators, gen.mutators)
    def __repr__(self):
        return str(self.graph)

def check_bound(arr, boundary):
        return (arr==boundary-1).any() or (arr==0).any()
  
def set_graph(params):
    l = int(params['boundary'])
    d = int(params['dim'])
    t = tuple(l*np.ones(d,dtype = int))
    return np.zeros(t,dtype=int)

def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)
def load_object(filename):
    with open(filename,'rb') as outp:
        obj = pickle.load(outp)
        return obj
"""Class to allow for cells to be added, removed, and randomly selected
Caveat: does not allow for efficient movement of cells unless key is cell ID
"""
class ListDict(object):
    def __init__(self):
        self.item_to_position = {}
        self.items = np.array([])
    def add_item(self, item):
        if item in self.item_to_position:
            return
        self.items = np.append(self.items, item)
        self.item_to_position[item.ID] = self.items.shape[0]-1
    def remove_item(self, item):
        position = self.item_to_position.pop(item.ID)
        last_item = self.items[-1]
        self.items = self.items[:-1]
        if position != self.items.shape[0]:
            self.items[position] = last_item
            self.item_to_position[last_item.ID] = position

    def choose_random_item(self):
        return np.random.choice(self.items)
    def get_item(self, index):
        return self.items[self.item_to_position[index]]
    def len(self):
        return self.items.shape[0]
    def __repr__(self):
        return f'{str([item.__repr__() for item in self.items])}'
    
"""result of treatment: periodic changes in cell death rate"""
def set_periodic_death_rate(by = 'time',interval = None):
    if by=='time':
        if interval is None:
            interval = 30 #set to every 30 days
    else:
        #interval is n_cells
        #TODO
        pass
        return interval


"""Using params file, get cell death rate for a given configuration"""

def one_fixed_death_rate(cell,death_rate):
    return death_rate
    
def one_changing_death_rate(cell,init_death_rate,driver_dr_factor):
    return init_death_rate*np.power(driver_dr_factor, cell.gen.n_drivers)
    
def radial_death_rate(cell,radius, inner_rate, outer_rate, driver_dr_factor):
    a = np.array(cell.pos)
    b = np.array(cell.sim.tumor.center)
    if np.linalg.norm(a-b) < radius:
        return inner_rate*np.power(driver_dr_factor, cell.gen.n_drivers)
    return outer_rate*np.power(driver_dr_factor, cell.gen.n_drivers)

def programmed_death_rate(time, start=60,end=130,dr1 = .1, dr2 = .65, ):
    if time < start or time > end:
        return dr1
    return dr2
    
if __name__ == "__main__":
    
   pass