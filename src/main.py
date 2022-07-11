#file with end-user functions to set parameters and run simulation 
import numpy as np
import pickle
import classes
import os
from pathlib import Path
import sys
import inspect
import json

#absolute path to parameter file 
PARAMS_PATH = '' #change to program directory 
#defualt parameters: 
DIM = 2
INIT_BIRTH_RATE = 1#np.log(2)
FIXED_DEATH_RATE = False
DRIVER_ADVANTAGE = .1
INIT_DEATH_RATE = .95*INIT_BIRTH_RATE
MAX_DEATH_RATE = INIT_BIRTH_RATE

PAD = 5 #fraction of radius to pad boundary 
MUTATOR_FACTOR = .01
BOUNDARY = 300
MAX_ITER = int(1e9)
MAX_POP = 50000
PUSH = 0
INIT_MUT_PROB = 0#.02 #taken from Waclaw et al. 2015 
DRIVER_PROB = 4e-5

MUTATOR_PROB = 0#4e-6 #need to implement way of increasing mutation probability in individual cell lines 
PROGRAMMED_DEATH_RATE = False
SEED = 123
#SAVE_INTERVAL = 1000 #number of cells to save after 
"""return spherical radius of tumor"""
def calc_radius(n_cells, dim):
    if dim==2:
        return np.sqrt(n_cells/np.pi)
    else:
        return np.cbrt(3*n_cells/4/np.pi)
        
"""Using params file, get cell death rate for a given configuration"""
def one_fixed_death_rate(cell,death_rate):
    return death_rate
    
def one_changing_death_rate(cell,init_death_rate,driver_dr_factor):
    return init_death_rate*np.power(driver_dr_factor, cell.gen.n_drivers)
    
def radial_death_rate(cell,radius, inner_rate, outer_rate, driver_dr_factor):
    a = np.array(cell.pos)
    b = np.array(cell.sim.tumor.center)
    if np.linalg.norm(a-b) < radius:
        #print('in')
        return inner_rate*np.power(driver_dr_factor, cell.gen.n_drivers)
    #print('out')
    return outer_rate*np.power(driver_dr_factor, cell.gen.n_drivers)
def radial_prop_death_rate(cell, prop, inner_rate, outer_rate, driver_dr_factor):
    """ inner death rate if within prop% of radius, otherwise inner"""
    n_cells = cell.sim.tumor.N 
    n_driv = cell.gen.n_drivers
    r = calc_radius(n_cells, cell.sim.params['dim']) #assuming 2 dimensions out of laziness
    a = np.array(cell.pos)
    b = np.array(cell.sim.tumor.center)
    #print(f'r = {r}')
    #print(f'dist = {np.linalg.norm(a-b)} and prop*r = {prop*r}')
    if np.linalg.norm(a-b) < prop*r:
        #print('inside radius')
        return inner_rate*np.power(driver_dr_factor, n_driv)
    #print('outside radius')
    return outer_rate*np.power(driver_dr_factor,n_driv)


def radial_gradient_death_rate(cell, radius, inner_rate, outer_rate, driver_dr_factor):
    """Rather than the death rate changing at a fixed radius, death rate changes as a gradient """

def radial_random_death_rate(cell, radius, inner_rate, outer_rate, dirver_dr_factor):
    """"""
    

def programmed_death_rate(time, start=60,end=130,dr1 = .1, dr2 = .65, ):
    if time < start or time > end:
        return dr1
    return dr2

DR_FUNCTIONS = {'default': one_fixed_death_rate,
'radial':radial_death_rate,'one_changing_death_rate':programmed_death_rate, 
'radial_prop': radial_prop_death_rate}

DR_PARAMS = dict() 

def check_fcn_args(fcn,params,kwargs):
    args = inspect.getfullargspec(fcn)[0]
    for a in args:
        if a != 'cell':
            if a not in params:
                try:
                    assert(a in kwargs)
                    params[a] = kwargs[a]
                except(AssertionError):
                    print(f'{a} must be defined when using {fcn.__name__}')
                    print('exiting...')
                    sys.exit() 
def set_nbrs(neighborhood = 'moore', dim = DIM):
    NBRS = []
    if neighborhood == 'moore':
            for i in [0,-1,1]:
                for j in [0,-1,1]:
                    if dim ==2:
                        NBRS.append([i,j])
                    else:
                        for k in [0,-1,1]:
                            NBRS.append([i,j,k])
            NBRS = np.array(NBRS)[1:]
    else:   
        if dim==3:
            NBRS = np.array([[1,0,0], [-1,0,0], [0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
        else:
            NBRS = np.array([[1,0],[-1,0],[0,1],[0,-1]])
    return NBRS

    
NBRS = set_nbrs()

def get_kwargs_from_file(path):
    if path is None:
        return {}

    ext = path.split('.')[-1]
    obj = None
    if ext == 'pkl':
        obj =  classes.load_object(path)
    if ext == 'txt':
        with open(path,'r') as f:
            obj = json.load(f)
        f.close()
    try:
        assert(type(obj) is dict)
    except(AssertionError):
        print(f'expected dict but got {type(obj)}')
        print('exiting...')
        sys.exit() 
    return obj
def calc_radius(n_cells, dim):
    if dim==2:
        return np.sqrt(n_cells/np.pi)
    else:
        return np.cbrt(3*n_cells/4/np.pi)

"""wrapper to handle user parameters. Accepts a file or a set of keyword arguments """
def config_params(kwargs):
   
    if 'dim' not in kwargs:
            kwargs['dim']=DIM
    if 'n_cells' not in kwargs:
            kwargs['n_cells'] = MAX_POP
    if 'neighborhood' not in kwargs:
            kwargs['neighborhood'] = 'moore'
    if kwargs['neighborhood']!= 'moore':
        kwargs['neighborhood'] = 'von_neumann'
    if 'fixed_death_rate' not in kwargs:
            kwargs['fixed_death_rate'] = False
    if 'init_death_rate' not in kwargs:
            kwargs['init_death_rate'] = INIT_DEATH_RATE
    kwargs['death_rate'] = kwargs['init_death_rate']
    if 'init_birth_rate' not in kwargs:
            kwargs['init_birth_rate'] = INIT_BIRTH_RATE
    if 'driver_advantage' not in kwargs:
            kwargs['driver_advantage'] = DRIVER_ADVANTAGE
    if 'driver_prob' not in kwargs:
            kwargs['driver_prob'] = DRIVER_PROB
    if 'max_death_rate' not in kwargs:
        kwargs['max_death_rate'] = MAX_DEATH_RATE
    if 'mutator_factor' not in kwargs:
        kwargs['mutator_factor']= MUTATOR_FACTOR
    if 'max_iter' not in kwargs:
        kwargs['max_iter'] = MAX_ITER
    if 'init_mut_prob' not in kwargs:
        kwargs['init_mut_prob'] = INIT_MUT_PROB
    if 'mutator_prob' not in kwargs:
        kwargs['mutator_prob'] = MUTATOR_PROB
    if 'progression' not in kwargs:
        kwargs['progression'] = None
    if 'exp_path' not in kwargs:
        kwargs['exp_path'] = PARAMS_PATH
    if 'reps' not in kwargs:
        kwargs['reps'] = 1
    if 'save_interval' not in kwargs:
        kwargs['save_interval'] = MAX_POP #only save one 
    if 'rep_start' not in kwargs:
        kwargs['rep_start'] = 0
    if 'dr_params' not in kwargs:
        kwargs['dr_params'] = DR_PARAMS
    

    print('starting sanity checks...')
    if 'dr_function' not in kwargs:
        kwargs['dr_function'] = DR_FUNCTIONS['default']
    else:
        try:
            kwargs['dr_function'] = DR_FUNCTIONS[kwargs['dr_function']]
        except(KeyError):
            print('death rate function not defined')
            print('exiting...')
            sys.exit()
   
    #set dependent parameters, sanity checks
    kwargs['driver_dr_factor'] = 1-kwargs['driver_advantage'] #death rate factor
    r = calc_radius(kwargs['n_cells'],kwargs['dim'])
    kwargs['boundary'] = int(r*(1+PAD))
    if kwargs['fixed_death_rate']: kwargs['driver_prob'] = 0 
    kwargs['neighbors'] = set_nbrs(kwargs['neighborhood'],kwargs['dim'])

    check_fcn_args(kwargs['dr_function'],kwargs['dr_params'],kwargs)
    probs = ['init_mut_prob', 'init_birth_rate', 'init_death_rate', 'driver_prob', 'driver_dr_factor','max_death_rate']
    for p in probs:
        try:
            assert(1>=kwargs[p]>=0)
        except(AssertionError):
            print(f'parameter {p} must be between 0 and 1 inclusive')
            print('exiting...')
            sys.exit()
    try:
        if 'driver_dr_factor' in kwargs['dr_params']:
            assert(kwargs['driver_dr_factor']==kwargs['dr_params']['driver_dr_factor'])
    except(AssertionError):
        print('driver advantage must match')
        print('exiting...')
        sys.exit() 
    #write params to pkl file for rest of simulation to use
    
    with open(os.path.join(PARAMS_PATH,'params.pkl'), 'wb') as outfile:
        pickle.dump(kwargs, outfile, protocol=pickle.HIGHEST_PROTOCOL)    
    #write params to experiment location 
    Path(kwargs['exp_path']).mkdir(parents=True, exist_ok=True)
    with open(os.path.join(kwargs['exp_path'],'params.pkl'), 'wb') as outfile:
        pickle.dump(kwargs, outfile, protocol=pickle.HIGHEST_PROTOCOL)    
    return kwargs
  
    
def simulateTumor(**kwargs):
    params = config_params(kwargs)
    rep = params['rep_start']
    sim_list = []
    print('starting simulation...')
    while rep<params['reps']:
        sim = classes.Simulation(params)
        sim.run(rep) 
        sim_list.append(sim)
        rep+=1
    print('done!')
    return sim_list[0] if len(sim_list)==1 else sim_list

if __name__ == '__main__': 
    try:
        config_file = sys.argv[1]
    except(IndexError):
        config_file = None
    kwargs = get_kwargs_from_file(config_file)
    out = simulateTumor(**kwargs)
    
    #import simulation
    #import matplotlib.pyplot as plt
    #exp_path = 'Documents/current/tumor_stuff/hi_s_2d'
    #out = simulateTumor(dim =2, driver_advantage = 4,n_cells = 100000,exp_path = exp_path)
    #print(out.tumor.hit_bound)
    #simulation.plot_slice(out.tumor,ax = 'x')
    #simulation.plot_slice(out.tumor,ax = 'y')
    #simulation.plot_slice(out.tumor,ax = 'z')
    
    #dr_params = {'radius':10 ,'inner_rate': .1 ,'outer_rate': .9}
    #for r in [0,10,15,20,30,40,45,50]:
        #exp_path = f'../test/r={r}_di=.9_do=.1'
    
        #dr_params = {'radius':r ,'inner_rate': .9 ,'outer_rate': .1}
        #out = simulateTumor(dim =2, driver_advantage = .5,n_cells = 20000,exp_path = exp_path, driver_prob = 1e-2, dr_function = 'radial', dr_params = dr_params, save_interval= 5000)
        #ax = simulation.plot_drivers(out.tumor, by_fitness = True)[2]
        #ax = simulation.plot_circle(ax, out.tumor.center, r)
        #ax.set_title(f'r = {r}, inner dr = {dr_params["inner_rate"]}, outer dr = {dr_params["outer_rate"]}')
        #plt.show()
    #simulation.plot_tumor(out.tumor)
    
    