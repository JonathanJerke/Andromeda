##
##file format, across component-particles, ( index-space, index-momentum ),
##adding an 'index' on each row.
##

##explicitly delete the found/found.vector file to reset.



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

gto = pd.read_csv('sg.csv')
gto.index = gto['index']
gto = gto.drop('index',axis=1)

print(B,D,N)
for i in gto.index:
    g.setCurrent(Current = 0)
    gto1 = gto.loc[i]
    
    g.SG(gammaPy = reshape(gto1.values,(D,B,2)))
    g.to_file(reset = i==0, index = i )
