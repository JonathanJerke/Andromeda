##BUILDING A foundation from a list of SG coordinates

import andromedaPy as apy
import numpy as np

#l =   list(map( lambda x : int(2*(float(x))), open('lithium','r').read().split() ) )
#A = np.reshape(l,(913,3,3,2))

g = apy.galaxy()
ls = [1]
cs = [g.comps(10)]
ds = [[g.dims(lattice = 1.0)]]
bs = [g.bases('Sinc')]
g.spaces(ls,cs,ds,bs)
g.block(['component','diagonal','total-parallel','matrixElement-parallel','multiply-parallel','permute','permute-parallel','transfer'])

g.calculationInputs(RAMmax = 1,numVectors = 0,numNames = 0)
g.fieldInputs(canonRank = 1,nStates = 1,OpIndex = 0)
g.i()
g.rename('found')
g.SG(gammaPy = [[[0,1]]])
print(g.Current())
print(max(g.streams()),min(g.streams()))
