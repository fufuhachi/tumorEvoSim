#file with end-user functions to set parameters and run simulation 
from lib2to3.pgen2 import driver
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

DRIVER_ADVANTAGE = .1
INIT_DEATH_RATE = .95*INIT_BIRTH_RATE
MAX_DEATH_RATE = INIT_DEATH_RATE

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

def get_kwargs_from_file(path):
    if path is None:
        return {}

    ext = path.split('.')[-1]
    obj = None
    if ext == 'pkl':
        obj =  classes.load_object(path)
    if ext == 'txt' or ext == 'json':
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

    if 'first_rep' not in kwargs and 'last_rep' not in kwargs:
        kwargs['first_rep'] = 0
        kwargs['last_rep'] = 0
    elif 'last_rep' not in kwargs and 'first_rep' in kwargs:
        kwargs['last_rep'] = kwargs['first_rep'] + kwargs['reps']-1
    else:
        kwargs['reps'] = kwargs['last_rep'] - kwargs['first_rep']+1

    if 'n_migrations' not in kwargs:
        kwargs['n_migrations'] = 0
    
    if 'save_interval' not in kwargs:
        kwargs['save_interval'] = MAX_POP #only save one 
        
    if 'dr_params' not in kwargs:
        kwargs['dr_function'] = 'default'
        kwargs['dr_params'] = {'init_death_rate': kwargs['init_death_rate'], 'driver_dr_factor': 1-kwargs['driver_advantage']}

    if 'br_params' not in kwargs:
        kwargs['br_function'] = 'default'
        kwargs['br_params'] = {'init_birth_rate': kwargs['init_birth_rate']}

    if 'mutate_params' not in kwargs:
        kwargs['mutate_function'] = 'default'
        kwargs['mutate_params'] = dict()

    print('starting sanity checks...')
    #check replicate number validity
    if kwargs['first_rep'] < 0 or kwargs['last_rep'] < 0 or kwargs['reps'] < 0:
        print('first_rep, last_rep, and reps must be nonnegative \n exiting...')
        sys.exit() 
    if kwargs['last_rep'] - kwargs['first_rep'] +1 != kwargs['reps']:
        print(f'last_rep ({kwargs["last_rep"]}) - first_rep ({kwargs["first_rep"]}) +1 must equal number of reps ({kwargs["reps"]})  \n exiting...')
        sys.exit()
    
        
   
    #set dependent parameters, sanity checks
    r = calc_radius(kwargs['n_cells'],kwargs['dim'])
    kwargs['boundary'] = int(r*(1+PAD))
    
    

    #check_fcn_args(kwargs['dr_function'],kwargs['dr_params'],kwargs) #depricated
    
    #probs = ['init_mut_prob', 'init_birth_rate', 'init_death_rate', 'driver_prob', 'driver_dr_factor','max_death_rate']

    #for p in probs:
    #    try:
    #        assert(1>=kwargs[p]>=0)
     #   except(AssertionError):
     #       print(f'parameter {p} must be between 0 and 1 inclusive')
     #       print('exiting...')
     #       sys.exit()
    try:
        if 'driver_dr_factor' in kwargs['dr_params']:
            assert(kwargs['driver_dr_factor']==kwargs['dr_params']['driver_dr_factor'])
    except(AssertionError):
        print(f'driver advantage must match, dr_params has {kwargs["dr_params"]["driver_dr_factor"]} but params has {kwargs["driver_dr_factor"]}')
        print('exiting...')
        sys.exit() 
    #check animation params
    if 'frame_rate' in kwargs:
        try:
            assert(type(kwargs['frame_rate'])==int)
            assert(kwargs['n_cells'] >= kwargs['frame_rate'] >=0)
        except(AssertionError):
            print('invalid animation params\nexiting...')
            sys.exit()
    
    #write params to pkl file for rest of simulation to use
        
    with open(os.path.join(PARAMS_PATH,'params.pkl'), 'wb') as outfile:
        pickle.dump(kwargs, outfile, protocol=pickle.HIGHEST_PROTOCOL)    
    #write params to experiment location 
    Path(kwargs['exp_path']).mkdir(parents=True, exist_ok=True)
    with open(os.path.join(kwargs['exp_path'],'params.pkl'), 'wb') as outfile:
        pickle.dump(kwargs, outfile, protocol=pickle.HIGHEST_PROTOCOL)    
    if 'frame_rate' in kwargs:
        Path(os.path.join(kwargs['exp_path'], 'anim')).mkdir(parents=True, exist_ok=True)
    return kwargs
   
  
    
def simulateTumor(**kwargs):
    max_attempts = 100 #maximum number of times to try a particular rep before moving on
    params = config_params(kwargs)
    first_rep = params['first_rep']
    cur_rep = first_rep
    sim_list = []
    last_rep = params['last_rep']
    """print(f'params are:')
    for k in params.keys():
        print(f'{k}:\n\t{params[k]}\n')
        
    print('\n\n')"""
    
    print('starting simulation...')
    while cur_rep < last_rep +1: 
        attempts = 0
        sim = classes.Simulation(params)
        sim.run(cur_rep) 
        sim_list.append(sim)
        print(f'rep {cur_rep} complete')
        if sim.tumor.N > 0 or attempts == max_attempts:
           cur_rep+=1
        attempts +=1
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
    
    