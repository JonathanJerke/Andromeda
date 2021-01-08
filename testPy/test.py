import andromedaPy as apy
import h5py

g = apy.galaxy()
ls = [1,1,1]
cs = [g.comps(40),g.comps(40),g.comps(40)]
ds = [[g.dims(lattice = 0.333),g.dims(lattice = 0.333)],[g.dims(lattice = 0.333),g.dims(lattice = 0.333)],[g.dims(lattice = 0.333),g.dims(lattice = 0.333)]]
bs = [g.bases('Sinc'),g.bases('Sinc'),g.bases('Sinc')]
g.spaces(ls,cs,ds,bs)
g.block(['total','train','copy','component','diagonal','total-parallel','matrixElement-parallel','multiply-parallel','permute','permute-parallel','transfer'])
g.calculationInputs(RAMmax = 4,numVectors = 0,numNames = 0)
g.fieldInputs(canonRank = 1,nStates = 1,OpIndex = 0)


g.i()
g.rename('here')
g.gaussian(width = 1.0)
g.to_file()

print(h5py.File('here.1.0_mac','r'))
print(list(h5py.File('here.1.0_mac','r')['  0-0'])[:10])

g.read_file('here.vector')
print(g.vectors())
print(g.auxVectors())
del g
