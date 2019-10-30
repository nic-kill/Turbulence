import schwimmbad
import numpy as np

def func(i):
    '''
    A useless function
    '''
    print(str(i+1))
    return i

# Use multipool - same as multiprocessing
with schwimmbad.MultiPool() as pool:
    inputs = [i for i in np.arange(0,10,2)]
    out1 = list(pool.map(func, inputs))

# Use serial pool
with schwimmbad.SerialPool() as pool:
    inputs = [i for i in np.arange(10,20,2)]
    out2 = list(pool.map(func, inputs))


print(out1, out2)
