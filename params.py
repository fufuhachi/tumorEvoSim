#constants
import numpy as np
import main
#user-defined params

DIM = 2
INIT_BIRTH_RATE = 0.69
FIXED_DEATH_RATE = False
MOORE = True
DRIVER_FACTOR = .95
INIT_DEATH_RATE = .95*INIT_BIRTH_RATE
MAX_DEATH_RATE = INIT_BIRTH_RATE


MUTATOR_FACTOR = .01
BOUNDARY = 300
MAX_ITER = int(1e9)
MAX_POP = 50000
PUSH = 0
INIT_MUT_PROB = .02 #taken from Waclaw et al. 2015 
DRIVER_PROB = 4e-5
if FIXED_DEATH_RATE: DRIVER_PROB = 0
MUTATOR_PROB = 0#4e-6 #need to implement way of increasing mutation probability in individual cell lines 
PROGRAMMED_DEATH_RATE = False
SEED = 123
def set_nbrs(moore = MOORE, dim = DIM):
    NBRS = []
    if moore:
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


        

