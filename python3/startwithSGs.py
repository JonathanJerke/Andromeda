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

print(B,D,N)
reset = 0
index = 0
for i in gto.index:
    g.setCurrent(Current = 0)
    gto1 = gto.loc[i]
    
    g.SG(gammaPy = reshape(gto1.values,(D,B,2)))
    g.to_file(reset = reset, index = index )
    index += 1
    reset = 0


