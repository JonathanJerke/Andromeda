###-------------###
####andromeda.pyx##
###-------------###


#*  Copyright 2021 Jonathan Jerke and Bill Poirier.
#*  We acknowledge the generous support of Texas Tech University,
#*  the Robert A. Welch Foundation, and the Army Research Office.


#*   *   This file is part of Andromeda.

#*   *   Andromeda is free software: you can redistribute it and/or modify
#*   *   it under the terms of the GNU General Public License as published by
#*   *   the Free Software Foundation, either version 3 of the License.

#*   *   Andromeda is distributed in the hope that it will be useful,
#*   *   but WITHOUT ANY WARRANTY; without even the implied warranty of
#*   *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#*   *   GNU General Public License for more details.

#*   *   You should have received a copy of the GNU General Public License
#*   *   along with Andromeda.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
include "system.pxi"

from constants cimport inta
from constants cimport floata
from constants cimport field
from constants cimport calculation
from constants cimport division

from Model cimport initCal
from Model cimport initField
from Model cimport iModel
from Model cimport fModel
from input cimport blockA
from input cimport allowQ
from input cimport resetA
from input cimport readShell
from ioPrint cimport tLoadEigenWeights
from ioPrint cimport printOut
from coreUtil cimport vectorLen
from coreUtil cimport tBoot
from coreUtil cimport printExpectationValues
from coreUtil cimport streams
from coreUtil cimport tMatrixElements
from coreUtil cimport SG
from coreUtil cimport GTO
from coreUtil cimport CanonicalRank
from coreUtil cimport defSpiralMatrix

from Decompose cimport CanonicalRankDecomposition
from Decompose cimport canonicalRankDecomposition

from constants cimport dimensions_label
from constants cimport metric_label
from constants cimport function_label

from constants cimport basisElementType
from constants cimport componentType
from constants cimport functionType
from constants cimport metricType
from constants cimport bodyType
from constants cimport blockMemoryType

from Compression cimport canonicalRankCompression

from libc.string cimport strcpy


cdef class galaxy:
	cdef calculation calculation
	cdef field field

	def __cinit__(self):
		self.calculation = initCal()
		self.field = initField()
		self.calculation.rt.NLanes = 1
		self.calculation.rt.NSlot = 1

	def __dealloc__(self):
		fModel(&self.field.f)
	
	def isbooted(self):
		return self.field.f.bootedMemory == 1
	
	def read_record(self,filepy):
		"""Give this member function a file in a directory structure of a linux computation,
		it will load up the same initial state.
		
		Parameters
		----------
		filepy : str
		
		Returns
		-------
		self	
		"""
		if self.isbooted() :
			print("warning, already booted")
			return self

		readShell(1, [str(filepy).encode('utf-8')],&self.calculation,&self.field)
		return self
        
	def dims(self, lattice:floata = 1, attack:floata =0.5, origin:floata =0.0,
												 anchor:floata =0.5 ):
		"""if not booted Returns new floata parameters for a linear dimension.
		else returns dims in current galaxy.
		
		Parameters
		----------
		lattice : floata 
		attack  : floata
		origin  : floata
		anchor  : floata		
		
		Returns
		-------
		either one dimensions_label or all
		"""
		if not self.isbooted():
			return dimensions_label(lattice = lattice , attack = attack, 
									origin = origin, anchor = anchor )
		else :
			all = []
			for space in range(SPACE):
				if self.field.f.canon[space].body != bodyType.nada:
					some = []
					for body in range(1,self.field.f.canon[space].body+1):
						particle = self.field.f.canon[space].particle[body]
						particle.origin += ( 
		particle.lattice*(self.field.f.canon[space].count1Basis-1)*particle.anchor)	
						some += [particle]
					all += [some]
			return all
			
	def comps(self, N = 1 , inc=1, i=1,periodic = False):
		"""Returns new component enum for a linear dimension.
		"""
		if not self.isbooted():
			if i <= 0 :
				return [N,inc,componentType.nullComponent]
			if not periodic:
				if i == 1 :
					return [N,inc,componentType.spatialComponent1]
				elif i == 2 :
					return [N,inc,componentType.spatialComponent2]
				elif i == 3 :
					return [N,inc,componentType.spatialComponent3]
			if periodic:
				if i == 1 :
					return [N,inc,componentType.periodicComponent1]
				elif i == 2 :
					return [N,inc,componentType.periodicComponent2]
				elif i == 3 :
					return [N,inc,componentType.periodicComponent3]
		else:
			all = []
			for space in range(SPACE):
				if self.field.f.canon[space].body != bodyType.nada:
					all += [[self.field.f.canon[space].count1Basis,self.field.f.canon[space].count1Inc,self.field.f.canon[space].component]]
			if len ( all ) > 0 :
				return all
			else :
				return [0,inc,componentType.nullComponent]

	def bases(self, desc : str = 'Sinc' ):
		"""Returns enumeration of basis types
		"""
		if not self.isbooted():
			names = dict({'Sinc':basisElementType.SincBasisElement,
					'Gaussian':basisElementType.GaussianBasisElement,
					'Dirac':basisElementType.DiracDeltaElement,
					'State':basisElementType.StateBasisElement,
					'overlap':basisElementType.overlapBasisElement})
			return names[desc]
		else :
			all = []
			for space in range(SPACE):
				if self.field.f.canon[space].body != bodyType.nada:
					all += [self.field.f.canon[space].basis]
			return all

	def labels(self):
		"""Returns enumeration of basis types
		"""
		all = []
		for space in range(SPACE):
			if self.field.f.canon[space].body != bodyType.nada:
				all += [self.field.f.canon[space].label]
		return all

	def spaces(self,labels:[inta], comps, dims, bases ):
		"""Absolute definition of space
        Three equal length vectors are conjoined.
        
        Parameters
        ----------
        labels : [int]
        		Per Space , bundle names
        comps  : [self.comps]
        		Per Space , component lengths
        dims   : [[self.dims]]
        		Per Space and Particle , float character
        bases  : [self.bases]
        		Per Space , type
        Sets internal clock of definitions
        
        Returns 
        -------
        self
		"""
		if not self.isbooted():
			if len(comps) == len(dims) :
				if len(comps) == len(bases):
					if len(comps) <= SPACE:
						if self.field.f.bootedMemory == 1 :
							print("warning, already booted")
							return self
						"""Define 
						"""		
					
						for (space,comp) in enumerate(comps):
							count = 0
							for labelish in labels:
								if labelish == labels[space]:
									count += 1
							self.field.f.canon[space].component = count				
							self.field.f.canon[space].space = comp[2]
							self.field.f.canon[space].count1Basis = comp[0]
							self.field.f.canon[space].count1Inc = comp[2]
						for (space,label) in enumerate(labels):
							self.field.f.canon[space].label = label
						for (space,dim) in enumerate(dims):
							for (body, particle) in enumerate(dim):
								self.field.f.canon[space].particle[body+1] = particle
								self.field.f.canon[space].particle[body+1].origin -= ( 
		particle['lattice']*(self.field.f.canon[space].count1Basis-1)*particle['anchor'])	
							if len(dim)==1:
								self.field.f.canon[space].body = bodyType.one
							elif len(dim)==2:
								self.field.f.canon[space].body = bodyType.two
							elif len(dim)==3:
								self.field.f.canon[space].body = bodyType.three
						for (space,base) in enumerate(bases):
							self.field.f.canon[space].basis = base
		return self
	
	def block(self, blockDescs:[str]):
		"""block memory allocations and some controls
		
		Parameters
		----------
		blockDesc : [str]
		
		Returns
		-------
		self
		"""
		blockNames = dict({'total':blockMemoryType.blockTotalVectorBlock,
		'train':blockMemoryType.blockTrainVectorsblock,
		'copy':blockMemoryType.blockCopyBlock,
		'transfer':blockMemoryType.blockTransferBasisblock,
		'matrixElements':blockMemoryType.blockMatrixElementsblock,
		'permute':blockMemoryType.blockPermutationsblock,
		'multiply-parallel':blockMemoryType.blockParallelMultiplyblock,
		'matrixElement-parallel':blockMemoryType.blockParallelMatrixElementblock,
		'permute-parallel':blockMemoryType.blockParallelPermuteblock,
		'print':blockMemoryType.blockPrintStuffblock,
		'total-parallel':blockMemoryType.blockTotalVectorParallelBlock,
		'component':blockMemoryType.blockComponentblock,
		'diagonal':blockMemoryType.blockDiagonalMatrixblock})
		resetA(&self.calculation.rt)
		for bl in blockDescs:
			blockA(&self.calculation.rt,blockNames[bl])
		return self
			
	def i(self):
		"""Initiate allocation, do not overwrite
		"""
		if self.isbooted() :
			print("warning, already booted")
			return self
		iModel(&self.calculation, &self.field)
		return self
						
	def calculationInputs ( self, numNames:inta=-1, numVectors:inta=-1, shiftFlag:inta=-1,
	Lambda:inta=-1
		,RAMmax=-1 , shiftFlag : inta = -1):
		"""Relevant calculation.input 's
		
		Parameters
		----------
		numNames : int
		numVectors : int
		shiftFlag : int
		Lambda : int
		RAMmax : int
		maxGB : int
		
		Returns
		-------
		self
		"""
		if self.isbooted() :
			print("warning, already booted")
			return self
		if RAMmax >= 0 :
			self.calculation.i.RAMmax = RAMmax
		if Lambda >= 0 :
			self.calculation.i.Lambda = Lambda
		if shiftFlag >= 0 :
			self.calculation.i.shiftFlag = shiftFlag
		if numVectors >= 0 :
			self.calculation.i.numVectors = numVectors
		if numNames >= 0 :
			self.calculation.i.numNames = numNames
		if shiftFlag >= 0 :
			self.calculation.i.shiftFlag = shiftFlag
		#print(self.calculation.i)
		return self
		
	def runTimeInputs ( self,boost : inta = -1, dynamic:inta=-1, tolerance:floata=-1.0,relativeTolerance : floata=-1.0,
	threshold:floata =-1.0, maxCondition : floata = -1.0, condition : floata = -1.0 , 
	maxCycle : inta = -1 ):
		"""Relevant calculation.input 's
		
		Parameters
		----------
		boost : inta
		dynamic : inta
		tolerance : floata
		relativeTolerance : floata
		threshold : floata
		Xcondition : floata
		alpha : floata
		maxCycle : inta
				
		Returns
		-------
		self
		"""
		if boost >= 0 :
			self.calculation.rt.boost = boost
		if dynamic >= 0 :
			self.calculation.rt.dynamic = dynamic
		if tolerance >= 0 :
			self.calculation.rt.TOLERANCE = tolerance
		if relativeTolerance >= 0 :
			self.calculation.rt.relativeTOLERANCE = relativeTolerance
		if threshold >= 0 :
			self.calculation.rt.THRESHOLD = threshold
		if Xcondition >= 0 :
			self.calculation.rt.XCONDITION = maxCondition
		if alpha >= 0 :
			self.calculation.rt.ALPHA = condition
		if maxCycle >= 0 :
			self.calculation.rt.MAX_CYCLE = maxCycle
		return self
		
	def fieldInputs( self, flex :inta = -1, OpIndex:inta  = -2 , body:inta  =-1,
		irrep:inta = -1, Iterations:inta = -1
	,nStates:inta = -1,iRank:inta = -1,canonRank:inta= -1,xRank:inta = -1,
	qFloor:inta = -1,filter:inta = -1,collect:inta=-1):
		"""Relevant field.input 's
		
		Parameters
		----------
		flex  : int
		OpIndex  :int
		body  : int
		irrep : int
		Iterations : int
		nStates: int
		iRank  : int
		canonRank : int
		xRank  : int
		qFloor : int
		filter :int
		collect:int
		
		Returns
		-------
		self
		"""
		if self.isbooted() :
			print("warning, already booted")
			return self
		if flex >= 0 :
			self.field.i.flex =flex
		if OpIndex >= -1 :
			self.field.i.OpIndex =OpIndex
		if body >= 1:
			self.field.i.body = body
		if irrep >= 0 :
			self.field.i.irrep = irrep
		if Iterations >= 0 :
			self.field.i.Iterations = Iterations
		if nStates >= 0:
			self.field.i.nStates = nStates
		if iRank >=0:
			self.field.i.iRank = iRank
		if canonRank >= 0:
			self.field.i.canonRank = canonRank
		if xRank >= 0:
			self.field.i.xRank = xRank
		if qFloor >= 0:
			self.field.i.qFloor = qFloor
		if filter >= 0:
			self.field.i.filter = filter
		if collect >= 0:
			self.field.i.collect = collect
		#print(self.field.i)
		return self
		
	def vectors(self, index :inta = 0):
		"""Vectors are addressed via these enumations.
		
		The number of them is by nStates.
		
		Returns
		-------
		division.eigenVectors
		"""
		return (int(division.eigenVectors)+index)
		
	def auxVectors (self,index : inta = 0):
		"""auxiliary Vectors are addressed via these enumations.
		
		Requires a booted galaxy.
		The number of them is by qFloor.
		
		Returns
		-------
		division arrayed after allocated eigenVectors
		"""

		if self.isbooted() :
			return (int(self.field.f.user)+index)
		else:
			return 0
		
	def read_file ( self, filename : str ,vector : division = division.eigenVectors, 
	collect : inta = 0 ):
		"""Standard Input procedure
		
		Parameters
		----------
		filename : str
		vector : division
		collect : inta
		
		Returns 
		-------
		inta 
			Number of vectors loaded
		"""
		cdef inta count = 0
		tLoadEigenWeights (  &self.calculation, self.field ,filename.encode('utf-8'), 
				&count,  vector, collect)
		return count
		
	def rename(self , name:str ):
		"""Name calculation.
		
		Parameters
		----------
		name : str
		
		Returns 
		-------
		self
		"""
		strcpy(self.calculation.name , name.encode('utf-8'))
		return self
				
	def to_file ( self, vector : division = division.eigenVectors, reset : inta = 1, index : inta = 0 ):
		"""Standard Input procedure
		
		Writes to calculation name.
		
		Parameters
		----------
		vector : division
		reset : inta
			Will overwrite the .vector file
		index : inta
			Meant for ease of indexing
		
		Returns 
		-------
		self
		"""
		printOut(  &self.calculation, self.field,reset, index, vector)
		return self
		
	def gaussian ( self, vector : division = division.eigenVectors, spin : inta = 0, 
	width : floata = 1.0):
		"""Places a correctly band-limited gaussian.
		
		Parameters
		----------
		vector : division
		spin : inta
		width : floata
		
		Returns
		-------
		self
		"""
		tBoot(self.field.f, vector, spin, width)
		return self
		
	def SG ( self, gammaPy : [[[int]]] ,amplitude :floata = 1.0, vector : division = division.eigenVectors, 
	spin : inta = 0):
		"""Places a correctly band-limited Symmetrized Gaussian.
		
		Parameters
		----------
		vector : division
		spin : inta
		gammaPy : [[[int]]]
		
		Returns
		-------
		self
		"""
		cdef inta gamma[SPACE*MAXBODY*2]
		index : inta  = 0
		for space in range(SPACE):
			if self.field.f.canon[space].body != bodyType.nada:
				for b in gammaPy[space]:	
					gamma[index] = b[0]
					gamma[index+1] = b[1]
					index += 2
		
		SG(self.field.f, vector, spin,amplitude, gamma)
		return self
	
	
	def GTO ( self, gammaPy : [[[float]]] , amplitude :floata = 1.0,
	vector : division = division.eigenVectors, spin : inta = 0):
		"""Places a correctly band-limited Gaussian Type Orbital.
		
		Parameters
		----------
		vector : division
		spin : inta
		gammaPy : [[[float]]]
		
		Returns
		-------
		self
		"""
		cdef inta   gamma[SPACE*MAXBODY]
		cdef floata delta[SPACE*MAXBODY*2]
		index : inta  = 0
		for space in range(SPACE):
			if self.field.f.canon[space].body != bodyType.nada:
				for b in gammaPy[space]:
					gamma[index] = int(b[0])
					delta[2*index] = b[1]
					delta[2*index+1] = b[2]
					index += 1
		
		GTO(self.field.f, vector, spin,amplitude, gamma,delta)
		return self
	
	
	def Current ( self, vector : division = division.eigenVectors , spin : inta = 0):
		return self.field.f.name[int(vector)].Current[spin]		
	
	def setCurrent ( self, vector : division = division.eigenVectors , spin : inta = 0, Current : inta = 0):
		self.field.f.name[int(vector)].Current[spin] = Current
		return self.field.f.name[int(vector)].Current[spin]			
	
	def streams( self, space : inta =0,vector: division = division.eigenVectors, 
	 index : inta = 0, spin : inta = 0, inputStream : [floata]= [] ):
		"""Streams will input/output the Andromeda structures.
		inputStream empty will lead to accessing Andromeda structures,
		otherwise Andromeda structures will be written to...in either case, Relevant
		[floata] will be outputted.	
		
		Parameters
		----------
		vector : division
		space  : inta
		index  : inta
		spin   : inta
		inputStream : [floata]
		
		Returns
		-------
		[floata]
		"""
		cdef double * pt = streams(self.field.f, vector, spin , space ) + vectorLen(self.field.f,space)*index
		if inputStream == []:
			outStream = []
			for c in range(vectorLen(self.field.f,space)):
				outStream += [pt[c]]
			return outStream
		else:
			for c in range(vectorLen(self.field.f,space)):
				pt[c] = inputStream[c]
			return inputStream
				
	def expectation(self, vector : division = division.eigenVectors):
		"""Print Expectation values.
		Resets blockMemory 10
		
		Returns
		-------
		self
		"""
		for b in range(BLOCK_COUNT):
			if self.calculation.rt.memBlock[b] == blockMemoryType.blockPrintStuffblock:
				self.calculation.rt.memBlock[b] = blockMemoryType.passBlock
		
		printExpectationValues (  &self.calculation,   self.field.f ,  division.Ha  , vector)
		return self

	def dot(self, vector : division = division.eigenVectors, matrix = division.nullOverlap,
	 vector2: division = division.eigenVectors):
		"""Print dot.
		Returns
		-------
		floata
		"""
		return tMatrixElements ( 0, self.field.f ,  vector, 0 , matrix, 0,vector2,0)

	def matmul(self, vectorIn : division = division.eigenVectors, matrix = division.Iterator,
	 vectorOut: division = division.eigenVectors, canonRank : inta):
		"""
		Parameters
		----------
		vectorIn : division
		matrix : division
		vectorOut : division
		canonRank :inta
				 
		Returns
		-------
		self
		"""
		tHXpY(self.field.f, vectorOut, defSpiralMatrix(&f1.f, matrix), 
		self.calculation.i.shiftFlag, vectorIn, 
		self.calculation.rt.TOLERANCE,
		self.calculation.rt.relativeTOLERANCE,
		self.calculation.rt.ALPHA,
		self.calculation.rt.THRESHOLD,
		self.calculation.rt.MAX_CYCLE,
		self.calculation.rt.XCONDITION,
		canonRank,
		self.calculation.rt.dynamic )
		return self

	def full( self, g: galaxy ,  vector : division  = division.eigenVectors ):
		"""Create a new galaxy with all elements explicitly written down.
		
		Limitations currently allow for 3D data Cubes (not rectangles).
		
		Parameters
		----------
		vector : division 
			The vector in self to analyze
			
		Returns
		-------
		galaxy
		"""
		cdef floata *cp[SPACE] 
		cdef floata *pt
		print(tMatrixElements(0,self.field.f,vector,0,
			division.nullOverlap,0,vector,0))
		for space in range(SPACE):
			if self.field.f.canon[space].body != bodyType.nada:
				cp[space] = streams(self.field.f,division.copyVector,0,space)
		
		blocks = ['component','diagonal','total-parallel',
		'matrixElement-parallel','multiply-parallel','permute','permute-parallel',
		'transfer']
		if allowQ(&self.calculation.rt,blockMemoryType.blockCopyBlock)==0:
			print('need copy block')
			return self
		
		comps = self.comps()
		spaces = 1
		cs = []
		xc = 0
		for comp in comps:
				c1 = comp[0]
				if xc < c1 :
					xc = c1
		cs = [g.comps(xc)]
		ds = []
		dims = self.dims()
		for dim in dims:
			for d in dim:
				ds += [d]
		ds = [ds]
		ls = [1]
		bs = [g.bases('Sinc')]
		g.spaces(ls,cs,ds,bs).block(blocks)
		g.calculationInputs(Lambda = 4,RAMmax = 1,numVectors = 0,numNames = 0)
		g.fieldInputs(canonRank = 1,nStates = 1,OpIndex = 0)
		g.i()
		pt = streams(g.field.f,division.eigenVectors,0,0)
		g.field.f.name[int(division.eigenVectors)].Current[0] = 2
		self.field.f.name[int(division.copyVector)].Current[0] = 1
		for ii in range(pow(cs[0][0],len(ds[0]))):
			iv = 1
			for space in range(SPACE):
				if self.field.f.canon[space].body != bodyType.nada:
					c1 = pow(comps[space][0],len(dims[space]))
					for c in range(c1):
						if c == (int(ii/iv)%c1):
							cp[space][c] = 1.0
						else:
							cp[space][c] = 0.0
					iv *= c1
			pt[ii] = tMatrixElements(0,self.field.f,division.copyVector,0,
			division.nullOverlap,0,vector,0)
		return g

	def compress ( self, spat: [[int]], g : galaxy , vector : division = division.eigenVectors,canonRank :inta = 1):
		"""self-> g
		testing...
		"""
		cdef inta spatial[SPACE][SPACE]
		for (ss,s) in enumerate(spat):
			for s2 in s:
				spatial[ss][s2] = 1
	

		canonicalRankCompression(spatial,NULL,self.field.f,0,NULL,(division.eigenVectors),0,
		1,0,g.field.f,1,vector,0,canonRank,0,
		self.calculation.rt.TOLERANCE,
		self.calculation.rt.relativeTOLERANCE,
		self.calculation.rt.ALPHA,
		self.calculation.rt.XCONDITION,
		self.calculation.rt.MAX_CYCLE)
		return self

	def decompose ( self, origin : division , ospin : inta = 0, 
	alloy : division = division.eigenVectors , spin = 0, canonRank : inta = 1 ):
		"""origin -> alloy
		"""
		
		if False:
			canonicalRankDecomposition( self.field.f , NULL,0,NULL, origin,0,
			CanonicalRank(self.field.f,origin,ospin),ospin, 1,alloy,0 , canonRank,  spin ,
		self.calculation.rt.TOLERANCE,
		self.calculation.rt.relativeTOLERANCE,
		self.calculation.rt.ALPHA,
		self.calculation.rt.XCONDITION,
		self.calculation.rt.MAX_CYCLE)
		
		
		return CanonicalRankDecomposition( self.field.f, NULL, origin, ospin , alloy, spin , 
		self.calculation.rt.TOLERANCE,
		self.calculation.rt.relativeTOLERANCE,
		self.calculation.rt.ALPHA,
		self.calculation.rt.THRESHOLD,
		self.calculation.rt.MAX_CYCLE,
		self.calculation.rt.XCONDITION,
		canonRank,
		self.calculation.rt.dynamic )
