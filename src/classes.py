
import numpy as np
import pickle
import sys

#classes and object methods
class Simulation():
    def __init__(self, params):
        self.params = params
        self.graph = set_graph(params)
        self.tumor = Tumor(self)
        self.t_traj = []
        self.N_traj = []
        self.nbrs = set_nbrs(params['neighborhood'], params['dim'])

    def stopping_condition(self):
        nmax = self.params['n_cells']
        imax = self.params['max_iter']
        return self.tumor.N == nmax or self.tumor.iter > imax or self.tumor.hit_bound or self.tumor.N ==0
    
    def run(self,rep=0):
        prev = 0
        animate_it = False
        if 'frame_rate' in self.params:
            animate_it = True
        while not self.stopping_condition():
            self.tumor.iterate()
            if self.tumor.iter%40000==0:
                print(f'size = {self.tumor.N} at {int(self.tumor.t/365)} years {self.tumor.t%365} days')
                if self.params['save_by'] == 'time':
                    save_object(self, f'{self.params["exp_path"]}/rep={rep}_ncells={self.tumor.N}_time={self.tumor.t:.2f}.pkl')
                    
            #    
            #else: 
            #    counter = self.tumor.N
            if self.tumor.N%self.params['save_interval']==0 and self.tumor.N!=prev:
                #print(f'saving with driviers {self.tumor.drivers}') #DEBUG
                save_object(self, f'{self.params["exp_path"]}/rep={rep}_ncells={self.tumor.N}_time={self.tumor.t}.pkl')
                prev = self.tumor.N #only save same cell size once 
            
            if animate_it:
                self.animate(rep)
        #save final 
        save_object(self, f'{self.params["exp_path"]}/rep={rep}_ncells={self.tumor.N}.pkl')
        return self.tumor

    def animate(self,rep):
        if self.tumor.N%self.params['frame_rate']==0:
            np.savetxt(f'{self.params["exp_path"]}/anim/rep={rep}_{self.tumor.N}.npz',self.tumor.graph)

       
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
        BIRTH_FUNCTIONS = {'default': self.one_fixed_rate, 'fixed': self.one_fixed_rate, 'one_changing': self.one_changing_rate,
        'radial': self.radial_rate}

        DEATH_FUNCTIONS = {'default': self.one_changing_rate,
        'radial':self.radial_rate,'one_changing':self.one_changing_rate, 
        'radial_prop': self.radial_prop_rate, 'radial_bern': self.radial_bern_rate, 
        'nbr_based': self.nbr_based_rate, 'resistance_model': self.resistance_model_death, 'radial_nbr_hybrid':self.radial_nbr_hybrid}

        MUTATE_FUNCTIONS = {'default': self.default_mutate, 'fixed_number': self.fixed_number_mutate, 
        'progression': self.progression_mutate,'resistance_model': self.resistance_model_mutate}

        PUSH_FUNCTIONS = {'default': self.default_push, 'prop': self.prop_push}

         
        
        try:
            self.br_function = BIRTH_FUNCTIONS[self.sim.params['br_function']]
        except(KeyError):
            print(f'birth rate function not defined. Must be one of {[k for k in BIRTH_FUNCTIONS.keys()]}')
            print('exiting...')
            sys.exit()

        try:
            self.dr_function = DEATH_FUNCTIONS[self.sim.params['dr_function']]
        except(KeyError):
            print(f'death rate function not defined. Must be one of {[k for k in DEATH_FUNCTIONS.keys()]}')
            print('exiting...')
            sys.exit()

        try: 
            self.mutate = MUTATE_FUNCTIONS[self.sim.params['mutate_function']]
        except(KeyError):
            print(f'mutate function not defined. Must be one of {[k for k in MUTATE_FUNCTIONS.keys()]}')
            print('exiting...')
            sys.exit()
        
        try:
            self.get_push_rate = PUSH_FUNCTIONS[self.sim.params['push_function']]
        except(KeyError):
            print(f'birth rate function not defined. Must be one of {[k for k in BIRTH_FUNCTIONS.keys()]}')
            print('exiting...')
            sys.exit()

    def get_birth_rate(self):
        """wrapper function that returns the cell birth rate"""
        return max(0,min(1,self.br_function(**self.sim.params['br_params'], is_birth = True)))
    def get_death_rate(self):
        """wrapper function that returns the cell deathr rate"""
        return max(0, min(1,self.dr_function(**self.sim.params['dr_params'], is_birth = False)))

    
    #Growth rate functions: either birth or death 
    def one_fixed_rate(self, init_rate, is_birth = False):
        """default birth rate, just return the given birth rate"""
        return init_rate
   
    """Using params file, get cell death rate for a given configuration"""
    
        
    def one_changing_rate(self,init_rate, is_birth = False):
        s = self.sim.params['driver_advantage']
        select_birth = self.sim.params['select_birth']
        if is_birth and select_birth:
            s = -s
        elif is_birth or select_birth:
            s = 0
        return init_rate*np.power(1-s, self.gen.n_drivers)
        
    def radial_rate(self,radius, inner_rate, outer_rate, is_birth = False):
        s = self.sim.params['driver_advantage']
        select_birth = self.sim.params['select_birth']
        if is_birth and select_birth:
            s = -s
        elif is_birth or select_birth:
            s = 0
        a = np.array(self.pos)
        b = np.array(self.sim.tumor.center)
        if np.linalg.norm(a-b) < radius:
            #print('in')
            return inner_rate*np.power(1-s, self.gen.n_drivers)
        #print('out')
        return outer_rate*np.power(1-s, self.gen.n_drivers)

    def radial_prop_rate(self, prop, inner_rate, outer_rate, is_birth = False):
        """ inner death rate if within prop% of radius, otherwise inner"""
        s = self.sim.params['driver_advantage']
        select_birth = self.sim.params['select_birth']
        if is_birth and select_birth:
            s = -s
        elif is_birth or select_birth:
            s = 0
        n_cells = self.sim.tumor.N 
        n_driv = self.gen.n_drivers
        r = calc_radius(n_cells, self.sim.params['dim']) 
        a = np.array(self.pos)
        b = np.array(self.sim.tumor.center)
        #print(f'r = {r}')
        #print(f'dist = {np.linalg.norm(a-b)} and prop*r = {prop*r}')
        if np.linalg.norm(a-b) < prop*r:
            #print('inside radius')
            return inner_rate*np.power(1-s, n_driv)
        #print('outside radius')
        return outer_rate*np.power(1-s,n_driv)

    def radial_bern_rate(self, prop, inner_rate, outer_rate, is_birth = False):
        """death rate outer_rate with prob p = min(d(cell, center)/radius,1) and inner_rate with 
        prob 1-p
        """
        s = self.sim.params['driver_advantage']
        select_birth = self.sim.params['select_birth']
        if is_birth and select_birth:
            s = -s
        elif is_birth or select_birth:
            s = 0
        n_cells = self.sim.tumor.N 
        n_driv = self.gen.n_drivers
        r = calc_radius(n_cells, self.sim.params['dim']) 
        a = np.array(self.pos)
        b = np.array(self.sim.tumor.center)
        pcell = np.linalg.norm(a-b)/r
        probability = pcell/2/prop if pcell <=prop else (pcell+1-2*prop)/(2*(1-prop))
        if np.random.random() < probability:
            return outer_rate*np.power(1-s, n_driv)
        return inner_rate*np.power(1-s, n_driv)

    def programmed_death_rate(time, start=60,end=130,dr1 = .1, dr2 = .65, ):
        if time < start or time > end:
            return dr1
        return dr2

    def nbr_based_rate(self, inner_rate, outer_rate, is_birth = False):
        s = self.sim.params['driver_advantage']
        select_birth = self.sim.params['select_birth']
        if is_birth and select_birth:
            s = -s
        elif is_birth or select_birth:
            s = 0
        empty_nbrs = self.get_empty_nbrs()
        n_driv= self.gen.n_drivers
        if len(empty_nbrs) > 0:
            return outer_rate*np.power(1-s, n_driv)
        return inner_rate*np.power(1-s, n_driv)

    def resistance_model_death(self, **dr_params):
        """if cell has fewer mutations than needed for resistance, specified by the mutation
        function parameter 'muts_to_res', then death rate is radial_death rate. Otherwise
        return the inner death rate. Assumed that inner rate is lower, as this is a model 
        of acquired resistance to treatment"""

        if self.gen.n_drivers < self.sim.params['mutate_params']['muts_to_res']:
            return self.radial_rate(**dr_params)
        else: 
            return dr_params['inner_rate']
    def radial_nbr_hybrid(self, init_radius, inner_rate, outer_rate, is_birth = False):
        """assume as in radial_rate that the outer death rate is experienced past a fixed radius
        However, beyond that radius, only cells with empty neighbors experience the outer death rate"""
        s = self.sim.params['driver_advantage']
        select_birth = self.sim.params['select_birth']
        if is_birth and select_birth:
            s = -s
        elif is_birth or select_birth:
            s = 0
        a = np.array(self.pos)
        b = np.array(self.sim.tumor.center)
        if np.linalg.norm(a-b) < init_radius: #if cell is in the inner region
            #print('in')
            return inner_rate*np.power(1-s, self.gen.n_drivers)
        else:
            return self.nbr_based_rate(inner_rate, outer_rate, is_birth)

    
    #MUTATE FUNCTIONS
        
    def default_mutate(self):
        """mutate function depends only on tumor-level mutation order (infinite allele assumption)
        Samples poisson-distributed number of passenger mutations, driver mutations, and mutator mutations 
        If total is nonzero, creates a new genotype with new mutations added. Each mutation assigned a unique integer based on the order it appeared. Genotype's parent is current genotype. 
        """
        #print('starting mutation')
        #print('current cells are:')
        #print(tumor.cells)
        #sample number of mutations
        n_passen = np.random.poisson(self.gen.passen_rate)
        n_drivers = np.random.poisson(self.gen.driver_rate)
        n_mutators = np.random.poisson(self.gen.mutator_rate)
        n = self.sim.tumor.next_snp
        if n_passen + n_drivers + n_mutators >0:
            passen = None
            drivers = None
            mutators = None
            if n_passen>0:
                passen = np.append(self.gen.passen, np.arange(n,n+n_passen))
                n+=n_passen
            if n_drivers >0:
                #print(f"drivers {np.arange(n,n+n_drivers)} at size {self.sim.tumor.N}")
                drivers = np.append(self.gen.drivers, np.arange(n,n+n_drivers))
                n+=n_drivers
                self.sim.tumor.total_drivers+=n_drivers
            if n_mutators >0:
                #print(f"mutator at size {self.sim.tumor.N}")
                mutators = np.append(self.gen.mutators, np.arange(n,n+n_mutators))
                n+=n_mutators
            new_gen = Genotype(sim = self.sim,ID = self.sim.tumor.gens.len()+1, parent = self.gen, passen = passen, drivers = drivers, mutators = mutators)
            self.gen.children = np.append(self.gen.children, new_gen)
            self.gen.remove()
            self.gen = new_gen
            self.sim.tumor.add_genotype(new_gen)
        self.sim.tumor.next_snp = n
    
    def passen_mutate(self):
        """function to produce mutations in default mutate but only for passenger mutations
        returns a mutations list containing the passenger mutations this round 
        """
        n = self.sim.tumor.next_snp
        n_new_passen  = np.random.poisson(self.gen.passen_rate)
        total_passen  = self.gen.passen
        if n_new_passen > 0:
            total_passen = np.append(total_passen, np.arange(n,n+n_new_passen))
        
        self.sim.tumor.next_snp += n_new_passen
        return total_passen

    def driver_mutate(self):
        """function to produce mutations in default mutate but only for driver mutations
        returns a mutations list containing the driver mutations this round 
        """
        n = self.sim.tumor.next_snp
        n_new_driv  = np.random.poisson(self.gen.driver_rate)
        total_drivers  = self.gen.drivers
        if n_new_driv > 0:
            total_drivers = np.append(total_drivers, np.arange(n,n+n_new_driv))
        n += n_new_driv
        self.sim.tumor.next_snp += n_new_driv
        self.sim.tumor.total_drivers += n_new_driv
        return total_drivers

    def update_genome(self, passengers = None, drivers = None, mutators = None):
        new_gen = Genotype(sim = self.sim,ID = self.sim.tumor.gens.len()+1, parent = self.gen, passen = passengers, drivers = drivers, mutators = mutators)
        self.gen.children = np.append(self.gen.children, new_gen)
        self.gen.remove()
        self.gen = new_gen
        self.sim.tumor.add_genotype(new_gen)
        

    def fixed_number_mutate(self, number = 1, by = 'random', value = 100):
        """allow a fixed number of mutations to appear in the whole population, and then forbid further mutations. Accepts a certain number of mutations, does not allow them afterwards"""
        #needs work 
        new_drivers = []
        if by == 'size': #once tumor is 
            # a certain size, give out number of mutations to a group of cells all at the same time (meant to ensure some survive drift) until there are the desired number of mutations in the population
            if self.sim.tumor.N >= value and len(self.sim.tumor.tracked_variants) < number:
                snp = self.sim.tumor.next_snp
                new_drivers.append(snp)
                self.sim.tumor.next_snp += 1
                print(f'added driver to cell {self}')
                self.sim.tumor.tracked_variants.add(snp)
                print(f'added, set is now {self.sim.tumor.tracked_variants}')

        
        n_passen = np.random.poisson(self.gen.passen_rate)
        n = self.sim.tumor.next_snp
        if n_passen + len(new_drivers)>0:
            passen = np.append(self.gen.passen, np.arange(n,n+n_passen))
            n+=n_passen
            new_gen = Genotype(sim = self.sim,ID = self.sim.tumor.gens.len()+1, parent = self.gen, passen = passen, drivers = np.append(self.gen.drivers, new_drivers))
            self.gen.children = np.append(self.gen.children, new_gen)
            self.gen.remove()
            self.gen = new_gen
            self.sim.tumor.add_genotype(new_gen)
        self.sim.tumor.next_snp = n

    def resistance_model_mutate(self,muts_to_res = 1):
        """model of resistance to treatment. cell mutation determines death rate
        if specified, selective advantage should be 0, as death rate is determined
        not by the advantage of a driver but by the existence of a certain number of 
        mutations 
        """
        #ndrivbefore = self.gen.n_drivers
        total_passen = self.passen_mutate()
        #assert(self.gen.n_drivers == ndrivbefore)
        total_drivers = self.gen.drivers
        if self.gen.n_drivers < muts_to_res:
           total_drivers = self.driver_mutate()
        if len(total_drivers) > self.n_drivers or len(total_passen) > self.n_passen: 
            self.update_genome(passengers = total_passen, drivers = total_drivers)
        #print(f'before was {ndrivbefore} now is {self.gen.n_drivers}')

        
        
        
    def progression_mutate(self, progression):
        """given a list of mutations that are allowed, only sample mutations if a cell has fewer than the allowed mutations"""
        
        raise(NotImplementedError)

    #PUSH FUNCTIONS

    def default_push(self):
        """return guaranteed probability of pushing"""
        return 0

    def prop_push(self, prop):
        """return push based on proportion of radius"""
        n_cells = self.sim.tumor.N 
        r = calc_radius(n_cells, self.sim.params['dim']) 
        a = np.array(self.pos)
        b = np.array(self.sim.tumor.center)
        return int(np.linalg.norm(a-b) >= prop*r)
            

    def get_empty_nbrs(self):
        nbrs = (self.sim.nbrs + np.array(self.pos))
        tups = tuple([tuple(col) for col in nbrs.T])
        ids = self.sim.graph[tups]
        return nbrs[ids==0]

    def get_all_nbrs(self):
        nbrs = (self.sim.nbrs + np.array(self.pos))
        return nbrs

    

    def __repr__(self):
        return str(f'Cell# {self.ID} with gen {self.gen.ID} at {self.pos}')
    
        

class Genotype():
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
    def __init__(self, sim, ID=1, parent=None, passen=None, drivers=None, mutators = None) -> None:
        self.sim = sim
        self.ID = ID
        self.passen = np.array([], dtype = int) if passen is None else passen
        self.n_passen = self.passen.shape[0]
        self.drivers = np.array([], dtype = int) if drivers is None else drivers
        self.n_drivers = self.drivers.shape[0]
        self.mutators = np.array([], dtype = int) if mutators is None else mutators
        self.n_mutators = self.mutators.shape[0] #TODO: justify mutator 
        self.passen_rate = self.sim.params['passen_rate']*np.power(self.sim.params['mutator_factor'], self.n_mutators)
        self.driver_rate = self.sim.params['driver_rate']*np.power(self.sim.params['mutator_factor'], self.n_mutators)
        self.mutator_rate = self.sim.params['mutator_rate']*np.power(self.sim.params['mutator_factor'], self.n_mutators)
        
        self.parent = parent
        self.children = np.array([],dtype = int)
        self.number = 1
        
        #add genotype to genotype list
    def add(self):
        #self.sim.tumor.add_driver_count(self.n_drivers) #not necessary if death rate is not divided by dmax #comment out! 
        self.number+=1
    def remove(self):
        if self.number ==0:
            print("call to remove extinct genotype!!")
            raise(Exception)
        else:
            #if self.number==1:
                #self.sim.tumor.remove_driver_count(self.n_drivers) 
            self.number-=1
            if self.ID in self.sim.tumor.tracked_variants:
                self.sim.tumor.tracked_variants.pop(self.ID)
                print(f'popped, now is {self.sim.tumor.tracked_variants}')

    def __repr__(self):
        parID = 'none' if self.parent is None else self.parent.ID
        return str(f'Genotype {self.ID} from {parID} with {self.n_drivers} drivers')
    

class Tumor():
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
        self.hit_bound = False 
        self.tracked_variants = set()
        self.bmax = self.params['br_params']['init_rate']
        self.iter = 0
        self.t = 0
        self.next_snp = 0
        self.drivers = np.array([])
        self.mutators = np.array([])
        self.total_drivers = 0
       
        if self.params['model'] == 'default' or self.params['model'] == 'bdg_nonSpatialDeath':
            self.model = self.bdg_nonSpatialDeath
        elif self.params['model'] == 'bdg_spatialDeath':
            self.model = self.bdg_spatialDeath
        elif self.params['model'] == 'nonSpatial':
            self.model = self.nonSpatialGrowth
        else:
            print(f'model {self.params["model"]} not defined!\nexiting...')
            sys.exit()

    def run_growth_model(self):
        """wrapper to run growth model"""
        return self.model(**self.params['model_params'])
    def update_time(self):
        """update the time of the simulation after a birth/death event """
        self.t+=np.log(2)/(self.bmax*self.N)
        self.iter+=1

    def iterate(self):
        self.update_time()
        self.run_growth_model() 
        migrate_event = 0
        if self.N > 0:
            if 0 < self.params['n_migrations'] < 1:
                if np.random.random() < self.params['n_migrations']:
                    self.nbr_swap()

            else:
                while migrate_event < self.params['n_migrations']:
                    self.nbr_swap()
                    migrate_event +=1
            

        if self.iter%400 == 0:
            self.sim.t_traj = np.append(self.sim.t_traj, self.t)
            self.sim.N_traj = np.append(self.sim.N_traj, self.N)


    def bdg_spatialDeath(self):
        """Function to simulate birth, death, and mutation. Note: Death also spatially determined here. Bug in old code. This represents Waclaw 2015's 'quiescent core' model"""
        
        cell = self.cells.choose_random_item()
        #print(f'chose {cell}')
        gen = cell.gen
        nbrs = cell.get_empty_nbrs()
        
        if nbrs.shape[0] >0:
            nbr = nbrs[np.random.randint(0,nbrs.shape[0])]
            #print(f'chose nbr {nbr}')
            self.hit_bound = check_bound(nbr, self.params['boundary'])
            br = cell.get_birth_rate()
            #print(br)
            if np.random.random() < br:
                new_cell = Cell(self.sim, ID = self.iter+1, gen = gen, pos = tuple(nbr))
                self.add_cell(new_cell)
            #print(f'cell pop is now: {self.cells}')
            #mutate daughter cell, sample to determine whether to kill parent cell. If no, mutate parent cell 
                new_cell.mutate(**self.params['mutate_params'])
            #if np.random.random() <cell.get_death_rate()/self.dmax: #comment out! 
            dr = cell.get_death_rate()
            #print(dr)
            if np.random.random() < dr:
                #kill cell
                #print(f'{cell} died')
                self.remove_cell(cell)
            else:
                cell.mutate(**self.params['mutate_params'])
        return

    def bdg_nonSpatialDeath(self):
        """growth where cells only give birth when there is empty space but die regardless. Note that the death rate function could kill cells dependent on spatial factors,
        but this function would be evaluated on all cells regardless of location or surrounding cells
         """

        cell = self.cells.choose_random_item()
        #print(f'chose {cell}')
        gen = cell.gen
        nbrs = cell.get_empty_nbrs()

        br = cell.get_birth_rate()
       
        #can_push = cell.get_push_rate()
        if nbrs.shape[0] >0:
            nbr = nbrs[np.random.randint(0,nbrs.shape[0])]
            #print(f'chose nbr {nbr}')
            self.hit_bound = check_bound(nbr, self.params['boundary'])
            
            #print(br)
            if np.random.random() < br:
                new_cell = Cell(self.sim, ID = self.iter+1, gen = gen, pos = tuple(nbr))
                self.add_cell(new_cell)
            #print(f'cell pop is now: {self.cells}')
            #mutate daughter cell, sample to determine whether to kill parent cell. If no, mutate parent cell 
                new_cell.mutate(**self.params['mutate_params'])
            #if np.random.random() <cell.get_death_rate()/self.dmax: #comment out! 
       

        dr = cell.get_death_rate()
        if np.random.random() < dr:
            self.remove_cell(cell)
        return

    def bdg_nonSpatialSeparateDeath(self):
        """version of bdg_nonSpatialDeath where the cell chosen for death is different from the cell chosen for birth and independently chosen.
        This is a control to ensure the former scheme does not bias the simulation results (statistically it shouldn't)"""
        birth_cell = self.cells.choose_random_item()
        #print(f'chose {cell}')
        gen = birth_cell.gen
        nbrs = birth_cell.get_empty_nbrs()

        br = birth_cell.get_birth_rate()
       
        #can_push = cell.get_push_rate()
        if nbrs.shape[0] >0:
            nbr = nbrs[np.random.randint(0,nbrs.shape[0])]
            #print(f'chose nbr {nbr}')
            self.hit_bound = check_bound(nbr, self.params['boundary'])
            
            #print(br)
            if np.random.random() < br:
                new_cell = Cell(self.sim, ID = self.iter+1, gen = gen, pos = tuple(nbr))
                self.add_cell(new_cell)
            #print(f'cell pop is now: {self.cells}')
            #mutate daughter cell, sample to determine whether to kill parent cell. If no, mutate parent cell 
                new_cell.mutate(**self.params['mutate_params'])
            #if np.random.random() <cell.get_death_rate()/self.dmax: #comment out! 
       
        death_cell = self.cells.choose_random_item()
        dr = death_cell.get_death_rate()
        if np.random.random() < dr:
            self.remove_cell(death_cell)
        return


    def nonSpatialGrowth(self,push_forward= False):
        """growth where cells are randomly pushed to make room for offspring at every timestep, not BDG"""
        cell = self.cells.choose_random_item()
        #print(f'chose {cell}')
        gen = cell.gen
        br = cell.get_birth_rate()
        
        if np.random.random() < br:
            if not push_forward:
                direction = choose_direction(self.graph, cell.pos, self.sim.nbrs)
            else:
                #choose lattice direction closest to radial vector
                direction = get_closest_vector(np.array(cell.pos) - np.array(self.center), self.sim.nbrs)
                
                
            birth_pos = cell.pos
            pushed = self.simple_pushing(cell = cell, direction = direction)
            if pushed:
                new_cell = Cell(self.sim, ID = self.iter+1, gen = gen, pos = birth_pos)
                self.add_cell(new_cell)
                new_cell.mutate(**self.params['mutate_params'])

        dr = cell.get_death_rate()

        if np.random.random() < dr:
            self.remove_cell(cell)
        return

    def simple_pushing(self, cell, direction):
        """simplest pushing scheme: shift cells along specified lattice direction, update their positions, starting position should now be empty, 
        direction is one of the possible directions specified by the model (6, 8, 26)
        """
        cur_id = cell.ID
        #print(f'cur_id = {cur_id}')
        cur_pos = cell.pos
        #print(f'cur_pos = {cur_pos}')
        traceback = []
        
        while cur_id > 0:
            traceback.append(cur_pos)
            cur_pos = tuple(np.array(cur_pos)+direction)
            #print(f'cur_pos = {cur_pos}')
            try:
                cur_id = self.graph[cur_pos]
                #assert(cur_id >=0)
            except(IndexError):
                self.hit_bound = True
                return False
            #except(AssertionError):
                ##forbidden zone
                #return False 
        #print(traceback)
        #have list of cells to move, move them
        empty_pos = cur_pos
        for pos in traceback[::-1]:
            #print(pos)
            cur_cell = self.cells.get_item(self.graph[pos])
            cur_cell.pos = empty_pos
            self.graph[empty_pos] = cur_cell.ID
            self.graph[pos] = 0
            empty_pos = pos
        return True

    
    def nbr_swap(self):
        """method to choose a random cell from the tumor, choose a random neighbor, swap locations
        input: none
        output: none 
        """
        cell = self.cells.choose_random_item()
        #print(f'chose {cell}')
        new_pos = cell.pos+self.sim.nbrs[np.random.randint(self.sim.nbrs.shape[0])]
        new_pos = tuple(new_pos)
        #print(f'moving to {new_pos}')
        try:
            occupant_id = self.graph[new_pos]
        except(IndexError):
            self.hit_bound = True
            return
        #print(f'occupied by {occupant_id}')
        if occupant_id>0:
                #swap
                other = self.cells.get_item(occupant_id)
                temp_pos = cell.pos
                cell.pos = other.pos
                other.pos = temp_pos
                self.graph[cell.pos] = cell.ID
                self.graph[other.pos] = other.ID
                pass
        else:
                self.graph[cell.pos] = 0 
                cell.pos = new_pos
                self.graph[cell.pos] = cell.ID
        return 
            
            #if unoccupied, move there, otherwise swap
            
    def add_cell(self, born_cell):
        """function to add cell to cell list, change relevant fields"""        
        self.graph[born_cell.pos] = born_cell.ID
        self.cells.add_item(born_cell)
        self.N+=1
        born_cell.gen.add()
    
    def remove_cell(self, doomed_cell):
        """function to remove cell from cell list change relevant fields"""
        self.graph[doomed_cell.pos]=0
        self.cells.remove_item(doomed_cell)
        #if self.N ==0:
            #print('tumor pop is empty')
            #raise(Exception)
        self.N-=1
        doomed_cell.gen.remove()
        

    def add_genotype(self, gen):
        self.gens.add_item(gen)
        self.drivers = np.append(self.drivers, gen.drivers)
        self.mutators = np.append(self.mutators, gen.mutators)
    def __repr__(self):
        return str(self.graph)


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


#helper functions
"""return spherical radius of tumor"""
def calc_radius(n_cells, dim):
    if dim==2:
        return np.sqrt(n_cells/np.pi)
    else:
        return np.cbrt(3*n_cells/4/np.pi)
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

def choose_direction(graph, start_pos, directions):
    """given a list of possible directions, choose randomly from non-infinite distances to 0 element"""
    distances = []
    for v in directions:
        cur_pos = np.array(start_pos)
        cur_id = graph[start_pos]
        dist = 0
        while cur_id > 0:
            cur_pos += v
            dist+=1
            try:
                cur_id = graph[tuple(cur_pos)]
            except(IndexError):
                dist = np.infty
                break
        distances = np.append(distances, dist)
    
    options = directions[distances < np.infty]
    direction = options[np.random.randint(len(options))]
    
    return direction

def set_nbrs(neighborhood, dim):
    """configure the neighboring cells on the lattice"""
    nbrs = []
    if neighborhood == 'moore':
            for i in [0,-1,1]:
                for j in [0,-1,1]:
                    if dim ==2:
                        nbrs.append([i,j])
                    else:
                        for k in [0,-1,1]:
                            nbrs.append([i,j,k])
            nbrs = np.array(nbrs)[1:]
    else:   
        if dim==3:
            nbrs = np.array([[1,0,0], [-1,0,0], [0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
        else:
            nbrs = np.array([[1,0],[-1,0],[0,1],[0,-1]])
    return nbrs




def get_closest_vector(v, vectors):
    """returns a vector from vectors that is the closest to v in angle
    inputs: v: 1d numpy array as row vector
            vectors: 2d numpy array where rows are the vectors of interest
    """
    if (v == 0).all():
        return vectors[np.random.randint(low = 0, high = len(vectors))]
    u = v/np.linalg.norm(v)
    return vectors[np.argmax(vectors@u.T)]


if __name__ == "__main__":
    import seaborn as sns
    import matplotlib.pyplot as plt
    #import simulation
    #np.random.seed(6)
    params = {'dim':2}#load_object("params.pkl")
    params['init_birth_rate'] = .6
    #params['dr_params']['init_death_rate'] = 0.6
    params['dr_params'] = {}
    params['br_params'] = {'init_rate':params['init_birth_rate']}
    params['dr_params']['inner_rate'] = .1
    params['dr_params']['outer_rate'] = .9
    params['dr_params']['radius'] = calc_radius(5000, 2)
    params['boundary'] = 300
    params['n_cells'] = 10000
    #params['n_migrations'] = 5
   
    params['driver_rate'] = .01
    params['dr_function'] = 'radial'
    params['br_function'] = 'default'
    params['model'] = 'nonSpatial'
    params['model_params'] = {'push_forward':True}
    #params['mutate_function'] = 'fixed_number'
    #params['mutate_params'] = {'number': 1, 'by': 'size', 'value': 5000}
    
    sys.exit()
    N = 0
    attempts = 0
    while N == 0 and attempts < 10000:
        test_sim = Simulation(params)
        test_sim.run()
        tum = test_sim.tumor
        N = tum.N
        attempts +=1
        print('there was an attempt')
    plt.plot(test_sim.t_traj, test_sim.N_traj)
    plt.show()
    #cell = tum.cells.get_item(1)
    g = simulation.get_fitness_graph(tum)
    sns.heatmap(g)
    plt.show()
    sns.heatmap(tum.graph)
    plt.show()

    g = np.zeros((4,4))
    idx = np.indices((4,4))
    g[idx[0][1]] = 1
    g[idx[0][2]] = 1
    g[idx[0][3]] = 1
    g[(3,1)] = 0
    print(g)
    dir = choose_direction(g, (2,2), test_sim.nbrs)
    print(dir)
    #print(tum)
    #for cell in tum.cells.items:
        #print(cell)
    #tum.simple_pushing(cell, [0,1])
    #for cell in tum.cells.items:
        #print(cell)
    #print(tum)

    
   