##SHOWS I CAN CREATE A SOP vector , expand it, and compress it again

import andromedaPy as apy

g = apy.galaxy()
cs = [g.comps(5),g.comps(5),g.comps(5)]
ds = [[g.dims(lattice = 0.333)],[g.dims(lattice = 0.333)],[g.dims(lattice = 0.333)]]
ls = [1,1,1]
bs = [g.bases('Sinc'),g.bases('Sinc'),g.bases('Sinc')]
g.spaces(ls,cs,ds,bs)
#g.block(['total','train','component','diagonal','total-parallel','matrixElement-parallel','multiply-parall>
g.block(['component','diagonal','total-parallel','matrixElement-parallel','multiply-parallel','permute','permute-parallel','transfer'])

g.calculationInputs(RAMmax = 1,numVectors = 0,numNames = 0)
g.fieldInputs(canonRank = 3,nStates = 1,OpIndex = 0)
g.i()
g.rename('here')
g.gaussian(width = 1.0)

f = apy.galaxy()
g.full(f)
print(f.isbooted())


f.compress(g)
