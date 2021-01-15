import andromedaPy as apy

g = apy.galaxy()
ls = [1,1,1]
cs = [g.comps(40),g.comps(40),g.comps(40)]
ds = [[g.dims(lattice = 0.333),g.dims(lattice = 0.333)],[g.dims(lattice = 0.333),g.dims(lattice = 0.333)],[g.dims(lattice = 0.333),g.dims(lattice = 0.333)]]
bs = [g.bases('Sinc'),g.bases('Sinc'),g.bases('Sinc')]
g.spaces(ls,cs,ds,bs)
g.block(['copy','component','diagonal','total-parallel','matrixElement-parallel','multiply-parallel','permute','permute-parallel','transfer'])
g.calculationInputs(RAMmax = 4,numVectors = 0,numNames = 0,Lambda = 1)
g.fieldInputs(canonRank = 4,nStates = 2,OpIndex = 0)

g.i()
g.gaussian(vector = g.vectors(2),width = 1.0)
g.gaussian(vector = g.vectors(2),width = 1.5)
g.gaussian(vector = g.vectors(2),width = 2.0)

g.decompose(origin = g.vectors(2),canonRank = 2)
del g
