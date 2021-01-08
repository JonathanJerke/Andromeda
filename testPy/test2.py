##create a galaxy then copy its definitions
import andromedaPy as apy

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
print(g.labels())
print(g.comps())
print(g.dims())
print(g.bases())

gg = apy.galaxy()
gg.spaces(g.labels(),g.comps(),g.dims(),g.bases())
gg.i()

print(gg.labels())
print(gg.comps())
print(gg.dims())
print(gg.bases())
