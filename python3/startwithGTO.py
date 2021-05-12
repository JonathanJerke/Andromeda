#!/usr/bin/env python3


##file format,  

##list in order component-particle information, ( n , gaussian-exponent, position ) 
##each row gets an 'index' under a column of that name.


##explicitly delete the found/found.vector file to reset.

from sys import argv

##first argv is opinional -- Ecut 

if len(argv) >= 2:
    Ecut = float(argv[1])
else:
    Ecut = 1000000000000000000000

import os
from numpy import *
import pandas as pd
from andromedaPy import galaxy
if os.getenv('LAUNCH')==None :
    print('warning,  bootup Jupyter after ignition.sh')
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

try:
    gto = gto[gto['pot'] < Ecut]
    gto = gto.drop('pot',axis=1)
    print('Ecut to ', len(gto))
except:
    1



print(B,D,N)
for i in gto.index:
    g.setCurrent(Current = 0)
    gto1 = gto.loc[i]
    
    g.GTO(gammaPy = reshape(gto1.values,(D,B,3)))
    g.to_file(reset = i==0, index = i )
