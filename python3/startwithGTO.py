##file format,  

##list in order component-particle information, ( n , gaussian-exponent, position ) 
##each row gets an 'index' under a column of that name.



import os
from numpy import *
import pandas as pd
from andromedaPy import galaxy
if os.getenv('LAUNCH')==None :
    print('warning,  bootup Jupyter after ignition.sh')
from scipy import linalg
import itertools

g = galaxy()
g.read_record('found/found').fieldInputs(nStates = 1, OpIndex = 0)
g.i()
dims = g.dims()
comps = g.comps()
B = len(g.dims()[0])
N = g.comps()[0][0]
D = len(comps)

gto = pd.read_csv('gto.csv')
gto.index = gto['index']
gto = gto.drop('index',axis=1)

print(B,D,N)
reset = 0
for i in gto.index:
    g.setCurrent(Current = 0)
    gto1 = gto.loc[i]
    
    g.GTO(gammaPy = reshape(gto1.values,(D,B,3)))
    g.to_file(reset = reset, index = i )
    reset = 0
