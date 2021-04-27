#!/usr/bin/python3

##make superior sg.csv file.

#J. Phys. Chem. Lett. 2020, 11, 15, 6468–6474 
##https://doi.org/10.1021/acs.jpclett.0c01435
#J. Chem. Phys. 152, 214102 (2020)  
##https://doi.org/10.1063/5.0005681

import os
from sys import argv

varTOT = 0.
if len(argv)!= 2:
	print(" Ecut ")
	exit()
	
Ecut  = float(argv[1])

from numpy import *
import pandas as pd
from andromedaPy import galaxy
if os.getenv('LAUNCH')==None :
    print('warning,  bootup Jupyter after ignition.sh')
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
	g.SG(gammaPy = reshape(go,(D,B,2)))
	su = 0.
	for term in g.terms():
		su += g.dot(vector = 100 , term = term , vector2 = 100 )
	return su

def leth ( qq, go ):
	"""least than?
	Parameters
	----------
	qq  : [2*d]
	gog : [2*d] 
	
	Returns
	-------
	True if lesser, else False"""
	
	for z in range(2*d):
		if ( qq[z] < go[z] ):
			return True
	return False
	
def grth ( qq, go):
	"""greater than?
	Parameters
	----------
	qq  : [2*d]
	gog : [2*d] 
	
	Returns
	-------
	True if greater, else False"""
	
	for z in range(2*d):
		if ( qq[z] > go[z] ):
			return True
	return False	
	
def addTo (ri,gogo,pot,qq,vbartemp):
	"""add a vector to list
	Parameters
	----------
	ri 
		order
	gogo
		vectors
	pot
		potential evaluations
	--additions--
		qq
		vbartemp
	
	Returns 
	-------
	ri,gogo,pot
	"""
	ri += [len(pot)]
	gogo += [qq]
	pot  += [vbartemp]
	return ri,gogo,pot
	
def sort (le,ri,gogo,qq):
	"""determine if unique vector coordinates
	Parameters
	----------
	le
	ri
	gogo
	qq
	
	Returns 
	-------
	le,ri,un
	"""
		
	un = 0
	for iii in range(len(gogo)):	
		if grth(qq,gogo[iii]) : 
			if (ri[iii]>0):		        
				iii = ri[iii]			
			else:							
				un=1                            
		elif leth(qq,gogo[iii]):	
			if (le[iii]>0):				
				iii = le[iii]				
			else:
				un=-1
	return le,ri,un


d = B * D
gogo =  [[0.,0.5] * d ]

arrae = []
for dd in range(d):
	##spatial +
	v = [ 0.] * 2 * d
	v[2*dd] = 1
	arrae += [v]
	
	#spatial - -
	v = [ 0.] * 2 * d
	v[2*dd] = -1
	arrae += [v]

	#momentum + 
	v = [ 0.] * 2 * d
	v[2*dd+1] = 1
	arrae += [v]

qq  = gogo[0]
pot = [epot(gogo[0])]

ci = 0
cf = 1

le = [-1]
ri = [-1]

while ( cf >= ci ):
	for i in range( ci,cf ):
		for dr in arrae:
			qq = list(array(gogo[i]) + array(dr))
			
			vbartemp = epot(qq)
		
			if ( vbartemp < Ecut ):
				le,ri,un = sort ( le,ri,gogo,qq)
			
				if ( un == 1 ):
					ri,gogo,pot = addTo(ri,gogo,pot, qq, vbartemp )
					le += [-1]
				elif ( un == -1 ):
					le,gogo,pot = addTo(le,gogo,pot, qq, vbartemp )
					ri += [-1]
					vbarTOT += vbartemp	
	ci = cf + 1
	cf = len(pot)

dat = pd.DataFrame( data = gogo )
dat['index'] = dat.index
dat.to_csv('sg.csv',index = False)
