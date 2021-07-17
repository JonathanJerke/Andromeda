#!/usr/bin/env python3


##edit these
#angular momentum of GTO
l = [0,1]
#exponent of GTO
b = [0.125]
#position of GTO
x = [0]
##edit these

from numpy import *
import pandas as pd
import os
from andromedaPy import galaxy
if os.getenv('LAUNCH')==None :
    print('warning,  bootup Jupyter after ignition.sh')
from scipy import linalg
import itertools

g = galaxy()
g.read_record('found/found').fieldInputs(nStates = 1, OpIndex = -1)
g.i()
dims = g.dims()
comps = g.comps()
B = len(g.dims()[0])
N = g.comps()[0][0]
D = len(comps)



def epot( go ):
	"""evaluate expectation value
	Parameters
	----------
	go
		vectorized coordinate in proper format
	
	Returns
	-------
	float"""
	
	g.setCurrent( Current = 0 ) 
	g.GTO(gammaPy = reshape(go,(D,B,3)))
	su = 0.
	for term in g.terms():
		su += g.dot(vector = 100 , term = term , vector2 = 100 )
	return su



gto = pd.DataFrame(data = list(itertools.product(l,b,x,l,b,x,l,b,x,l,b,x,l,b,x,l,b,x,l,b,x,l,b,x,l,b,x,l,b,x,l,b,x,l,b,x)))
gto['index'] = gto.index
gto['pot'] = 0.

for index in gto.index:
    gto.loc[index,'pot'] = epot(gto.drop(['pot','index'],axis =1 ).loc[index].values)
    
gto.to_csv('gto.csv',index = False)

print(histogram(gto['pot'].values))