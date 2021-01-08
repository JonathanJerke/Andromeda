##BUILDING A foundation from a list of SG coordinates

import andromedaPy as apy
import numpy as np

l =   list(map( lambda x : int(2*(float(x))), open('lithium','r').read().split() ) )
A = np.reshape(l,(913,3,3,2))

g = apy.galaxy()
ls = [1,1,1]
cs = [g.comps(25),g.comps(25),g.comps(25)]
ds = [[g.dims(lattice = 0.4),g.dims(lattice = 0.4),g.dims(lattice = 0.4)],[g.dims(lattice = 0.4),g.dims(lattice = 0.4),g.dims(lattice = 0.4)],[g.dims(lattice = 0.4),g.dims(lattice = 0.4),g.dims(lattice = 0.4)]]
bs = [g.bases('Sinc'),g.bases('Sinc'),g.bases('Sinc')]
g.spaces(ls,cs,ds,bs)
g.block(['component','diagonal','total-parallel','matrixElement-parallel','multiply-parallel','permute','permute-parallel','transfer'])

g.calculationInputs(RAMmax = 1,numVectors = 0,numNames = 0)
g.fieldInputs(canonRank = 1,nStates = 1,OpIndex = 0)
g.i()
g.rename('found/found')
reset = 1
index = 0
for a in A:
	g.SG(gammaPy = a)
	g.to_file(reset = reset, index = index )
	index += 1
	reset = 0
