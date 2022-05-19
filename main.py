#file with end-user functions to set parameters and run simulation 
import numpy as np
import pickle
#absolute path to parameter file 
PARAMS_PATH = '~/Documents/current/tumor_stuff/tumorEvoSim'
#defualt parameters: 
DIM = 2
INIT_BIRTH_RATE = 0.69
FIXED_DEATH_RATE = False
DRIVER_ADVANTAGE = .05
INIT_DEATH_RATE = .95*INIT_BIRTH_RATE
MAX_DEATH_RATE = INIT_BIRTH_RATE

PAD = .1 #fraction of radius to pad boundary 
MUTATOR_FACTOR = .01
BOUNDARY = 300
MAX_ITER = int(1e9)
MAX_POP = 50000
PUSH = 0
INIT_MUT_PROB = .02 #taken from Waclaw et al. 2015 
DRIVER_PROB = 4e-5

MUTATOR_PROB = 0#4e-6 #need to implement way of increasing mutation probability in individual cell lines 
PROGRAMMED_DEATH_RATE = False
SEED = 123


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
def set_params(*args):
    return args
    
NBRS = set_nbrs()
def get_kwargs_from_file(path):
    #TODO
    raise(NotImplementedError)
def calc_radius(n_cells, dim):
    if dim==2:
        return np.sqrt(n_cells/np.pi)
    else:
        return np.cbrt(3*n_cells/4/np.pi)
"""wrapper to take user parameters. Accepts a file or a set of keyword arguments """
def config_params(kwargs):
    if 'filePath' in kwargs:
        kwargs = get_kwargs_from_file(kwargs['filePath'])
    if 'dim' not in kwargs:
            kwargs['dim']=DIM
    if 'n_cells' not in kwargs:
            kwargs['n_cells'] = MAX_POP
    if 'neighborhood' not in kwargs:
            kwargs['neighborhood'] = 'moore'
    if 'fixed_death_rate' not in kwargs:
            kwargs['fixed_death_rate'] = False
    if 'death_rate' not in kwargs:
            kwargs['death_rate'] = INIT_DEATH_RATE
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
    if 'exp_path' not in kwargs:
        kwargs['exp_path'] = PARAMS_PATH

    
    #set dependent parameters, sanity checks
    kwargs['driver_dr_factor'] = 1-kwargs['driver_advantage'] #death rate factor
    r = calc_radius(kwargs['n_cells'],kwargs['dim'])
    kwargs['boundary'] = int(r*(1+PAD))
    if kwargs['fixed_death_rate']: kwargs['driver_prob'] = 0 
    kwargs['neighbors'] = set_nbrs(kwargs['neighborhood'],kwargs['dim'])

    #write params to pkl file for rest of simulation to use
    with open(PARAMS_PATH+'/params.pkl', 'wb') as outfile:
        pickle.dump(kwargs, outfile, protocol=pickle.HIGHEST_PROTOCOL)    
    #write params to experiment location 
    with open(kwargs['exp_path']+'/params.pkl', 'wb') as outfile:
        pickle.dump(kwargs, outfile, protocol=pickle.HIGHEST_PROTOCOL)    
    return kwargs
  
    
import simulation


if __name__ == '__main__': 
    pass
