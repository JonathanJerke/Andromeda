##BUILDING GTOs and taking inner products

import andromedaPy as apy
import numpy as np

l =   list(map( lambda x : float(x), open('GTOs','r').read().split() ) )
A = np.reshape(l,(int(len(l)/9),3,1,3))

##lists by dimension, (n,alpha,y)

g = apy.galaxy()
ls = [1,1,1]
cs = [g.comps(100),g.comps(100),g.comps(100)]
ds = [[g.dims(lattice = 0.1)],[g.dims(lattice = 0.1)],[g.dims(lattice = 0.1)]]
bs = [g.bases('Sinc'),g.bases('Sinc'),g.bases('Sinc')]
g.spaces(ls,cs,ds,bs)
g.block(['component','diagonal','total-parallel','matrixElement-parallel','multiply-parallel','permute','permute-parallel','transfer'])

g.calculationInputs(RAMmax = 1,numVectors = 0,Lambda = 1,numNames = 0)
g.fieldInputs(canonRank = 1,nStates = len(A),OpIndex = 0)
g.i()
reset = 1
index = 0
for (ii,a) in enumerate(A):
        g.GTO(gammaPy = a, vector = g.vectors(ii) )
        for iii in range(ii+1):
                print( ii, iii, g.dot(vector = g.vectors(ii) , vector2 = g.vectors(iii) ) ) 




