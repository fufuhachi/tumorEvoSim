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

DRIVER_ADVANTAGE = 0
INIT_DEATH_RATE = .95*INIT_BIRTH_RATE

PAD = 5 #fraction of radius to pad boundary 
MUTATOR_FACTOR = .01
BOUNDARY = 300
MAX_ITER = int(1e9)
MAX_POP = 50000
PUSH = 0
PASSEN_RATE = 0#.02 #taken from Waclaw et al. 2015 
DRIVER_RATE = 1e-2#4e-5

MUTATOR_RATE = 0#4e-6 #need to implement way of increasing mutation probability in individual cell lines 
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
            kwargs['dim']=DIM #dimension of simulation can be 2 or 3
    if 'n_cells' not in kwargs:
            kwargs['n_cells'] = MAX_POP 
    if 'neighborhood' not in kwargs:
            kwargs['neighborhood'] = 'moore'
    if kwargs['neighborhood']!= 'moore':
        kwargs['neighborhood'] = 'von_neumann'
    
    if 'driver_advantage' not in kwargs:
            kwargs['driver_advantage'] = DRIVER_ADVANTAGE
    if 'select_birth' not in kwargs:
            kwargs['select_birth'] = False #seleciton actions on death by default and can only act on 1 
    if 'driver_rate' not in kwargs:
            kwargs['driver_rate'] = DRIVER_RATE
    if 'mutator_factor' not in kwargs:
        kwargs['mutator_factor']= MUTATOR_FACTOR
    if 'max_iter' not in kwargs:
        kwargs['max_iter'] = MAX_ITER
    if 'passen_rate' not in kwargs:
        kwargs['passen_rate'] = PASSEN_RATE
    if 'mutator_rate' not in kwargs:
        kwargs['mutator_rate'] = MUTATOR_RATE
    if 'progression' not in kwargs:
        kwargs['progression'] = None
    if 'exp_path' not in kwargs:
        kwargs['exp_path'] = PARAMS_PATH
    if 'reps' not in kwargs:
        kwargs['reps'] = 1

    if 'first_rep' not in kwargs and 'last_rep' not in kwargs:
        kwargs['first_rep'] = 0
        kwargs['last_rep'] = kwargs['reps']-1
    elif 'last_rep' not in kwargs and 'first_rep' in kwargs:
        kwargs['last_rep'] = kwargs['first_rep'] + kwargs['reps']-1
    else:
        kwargs['reps'] = kwargs['last_rep'] - kwargs['first_rep']+1

    if 'n_migrations' not in kwargs:
        kwargs['n_migrations'] = 0
    
    if 'save_interval' not in kwargs:
        kwargs['save_interval'] = MAX_POP #only save one 
    if 'save_by' not in kwargs:
        kwargs['save_by'] = 'cells'#default to saving by cells
    elif kwargs['save_by'] == 'time':
        kwargs['save_by']
    
    #check growth params    
    if 'dr_params' not in kwargs and 'dr_function' not in kwargs:
        kwargs['dr_function'] = 'default'
        kwargs['dr_params'] = {'init_rate': kwargs['init_death_rate']}

    elif 'dr_params' in kwargs and 'dr_function' not in kwargs:
        print('must specify dr function if dr params are specified\nexiting...')
        sys.exit()
    elif 'dr_params' not in kwargs:
        kwargs['dr_params'] = {}
    
    if kwargs['dr_function'] == 'resistance_model':
        assert(kwargs['mutate_function'] == 'resistance_model')

    if 'br_params' not in kwargs and 'br_function' not in kwargs:
        kwargs['br_function'] = 'default'
        kwargs['br_params'] = {'init_rate': kwargs['init_birth_rate']}

    elif 'br_params' in kwargs and 'br_function' not in kwargs:
        print('must specify dr function if br params are specified\nexiting...')
        sys.exit()

    elif 'br_params' not in kwargs:
        kwargs['br_params'] = {}
    
    if 'mutate_params' not in kwargs:
        kwargs['mutate_function'] = 'default'
        kwargs['mutate_params'] = dict()

    if kwargs['mutate_function'] == 'resistance_model':
        assert(kwargs['dr_function'] == 'resistance_model')
        assert(kwargs['driver_advantage']==0) 
        print(f'resistance model selected. Death rate params are: {kwargs["dr_params"]}')
        
    if 'push_params' not in kwargs:
        kwargs['push_function'] = 'default'
        kwargs['push_params'] = dict()

    if 'model_params' not in kwargs:
        kwargs['model_params'] = {}

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
    
    
    #for p in probs:
    #    try:
    #        assert(1>=kwargs[p]>=0)
     #   except(AssertionError):
     #       print(f'parameter {p} must be between 0 and 1 inclusive')
     #       print('exiting...')
     #       sys.exit()
    
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
    #write params to text file for convenience
    with open(os.path.join(kwargs['exp_path'],'params.txt'), 'w') as outfile:
        outfile.write(json.dumps(kwargs))
    
    
    if 'frame_rate' in kwargs:
        Path(os.path.join(kwargs['exp_path'], 'anim')).mkdir(parents=True, exist_ok=True)
    return kwargs
   
  
    
def simulateTumor(**kwargs):
    max_attempts = 10000 #maximum number of times to try a particular rep before moving on
    params = config_params(kwargs)
    first_rep = params['first_rep']
    cur_rep = first_rep
    sim_list = []
    last_rep = params['last_rep']
    attempts = 0
    """print(f'params are:')
    for k in params.keys():
        print(f'{k}:\n\t{params[k]}\n')
        
    print('\n\n')"""
    
    print('starting simulation...')
    print(cur_rep)
    while cur_rep < last_rep +1: 
        sim = classes.Simulation(params)
        while sim.tumor.N < 2 and attempts < max_attempts:
           sim = classes.Simulation(params)
           print(f'trying rep {cur_rep}')
           sim.run(cur_rep) 
           attempts +=1
        cur_rep+=1
        sim_list.append(sim)
    print('done!')

    return sim_list[0] if len(sim_list)==1 else sim_list

if __name__ == '__main__': 
    try:
        config_file = sys.argv[1]
    except(IndexError):
        config_file = None
    try: 
        display = sys.argv[2]
    except(IndexError):
        display = False
    kwargs = get_kwargs_from_file(config_file)
    out = simulateTumor(**kwargs)
    if display:
        print('plotting tumor...')
        import simulation
        import seaborn as sns
        import matplotlib.pyplot as plt
        g = simulation.get_fitness_graph(out.tumor)
        sns.heatmap(g)
        plt.show()
        simulation.plot_growth(out.tumor)
        plt.show()
        summary = simulation.tumor_summary(out.tumor)
        simulation.tumor_scatter(summary['x'], summary['y'], summary['death_rate'])
        plt.title('death_rate distribution')
        plt.show()
        simulation.tumor_scatter(summary['x'], summary['y'], summary['birth_rate'])
        plt.title('birth_rate distribution')
        plt.show()

    
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
        #out = simulateTumor(dim =2, driver_advantage = .5,n_cells = 20000,exp_path = exp_path, driver_rate = 1e-2, dr_function = 'radial', dr_params = dr_params, save_interval= 5000)
        #ax = simulation.plot_drivers(out.tumor, by_fitness = True)[2]
        #ax = simulation.plot_circle(ax, out.tumor.center, r)
        #ax.set_title(f'r = {r}, inner dr = {dr_params["inner_rate"]}, outer dr = {dr_params["outer_rate"]}')
        #plt.show()
    #simulation.plot_tumor(out.tumor)
    
    