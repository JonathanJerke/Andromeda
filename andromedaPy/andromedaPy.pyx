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
cimport andromedaPy
include "system.pxi"

from constants cimport inta
from constants cimport floata
from constants cimport field
from constants cimport calculation
from constants cimport division

ctypedef inta** inta2


from Model cimport initCal
from Model cimport initField
from Model cimport iModel
from Model cimport fModel
from input cimport blockA
from input cimport allowQ
from input cimport resetA
from input cimport bootShell
from ioPrint cimport tLoadEigenWeights
from ioPrint cimport printOut
from coreUtil cimport tBoot
from coreUtil cimport printExpectationValues
from coreUtil cimport streams
from coreUtil cimport tMatrixElements

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
		if self.field.f.bootedMemory == 1 :
			print("warning, already booted")
			return self

		bootShell(1, [str(filepy).encode('utf-8')],&self.calculation,&self.field)
		return self
        
	def dims(self, lattice:floata, attack:floata =0.5, origin:floata =0.0,
												 anchor:floata =0.5 ):
		"""Returns new floata parameters for a linear dimension.
		"""
		return dimensions_label(lattice = lattice , attack = attack, 
									origin = origin, anchor = anchor )

	def comps(self, N , inc=1, i=1,periodic = False):
		"""Returns new component enum for a linear dimension.
		"""
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

		return [0,inc,componentType.nullComponent]


	def bases(self, desc : str ):
		"""Returns enumeration of basis types
		"""
		names = dict({'Sinc':basisElementType.SincBasisElement,
					'Gaussian':basisElementType.GaussianBasisElement,
					'Dirac':basisElementType.DiracDeltaElement,
					'State':basisElementType.StateBasisElement,
					'overlap':basisElementType.overlapBasisElement})
		return names[desc]

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
		if self.field.f.bootedMemory == 1 :
			print("warning, already booted")
			return self
		iModel(&self.calculation, &self.field)
		return self
						
	def calculationInputs ( self, numNames:inta=-1, numVectors:inta=-1, shiftFlag:inta=-1,
	Lambda:inta=-1
		,RAMmax=-1 ):
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
		if self.field.f.bootedMemory == 1 :
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
		print(self.calculation.i)
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
		if self.field.f.bootedMemory == 1 :
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
		print(self.field.i)
		return self
		
	def vectors(self):
		"""Vectors are addressed via these enumations.
		
		The number of them is by nStates.
		
		Returns
		-------
		division.eigenVectors
		"""
		return division.eigenVectors
				
	def auxVectors (self):
		"""auxiliary Vectors are addressed via these enumations.
		
		Requires a booted galaxy.
		The number of them is by qFloor.
		
		Returns
		-------
		division arrayed after allocated eigenVectors
		"""

		if self.field.f.bootedMemory == 1 :
			return self.field.f.user	
		else:
			return division.nullName
		
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
		
	def metric(self, funcDesc: str = 'Coulomb', intervalDesc : str = 'interval',
										 betas : [inta,inta] =[0,1],interval : inta = 7,
										  contr : inta = 2):
		"""Metric definition by description
		
		Parameters
		----------
		funcDesc : str
		intervalDesc : str
			Type of interval or Dirac
		betas : [floata,floata]
			interval span or first float only
		interval : int
			CanonRank of operator
		contr : int
			Off diagonals
			
		Returns
		-------
		metric_label
		"""
		funcNames = dict(
			{'null':functionType.nullFunction,'Pseudo':functionType.Pseudo,
			'Yukawa':functionType.Yukawa,'Coulomb':functionType.Coulomb,
			'Morse':functionType.Morse,'LennardJones':functionType.LennardJones,
			'LDA':functionType.LDA,'BLYP':functionType.BLYP,
			'Gaussian':functionType.Gaussian}
		)
		
		intervalName = dict({'dirac':metricType.dirac,
							'separateDirac':metricType.separateDirac,
							'interval':metricType.separateDirac,
							'semiDefinite':metricType.semiIndefinite,
							'pureInterval':metricType.pureInterval,
							'pureSemiIndefinite':metricType.pureSemiIndefinite}
		)
		
		zs = np.zeros(SPACE)
		return metric_label(pow = zs,powB = zs,deriv = zs,
			fn =function_label(interval = interval, contr = contr,
						fn = funcNames[funcDesc],param = np.zeros(MAX_PARAM_FUNC)) ,
			metric = intervalName[intervalDesc],beta = betas)
		
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
		
		
		blocks = ['total','copy','component','diagonal','total-parallel',
		'matrixElement-parallel','multiply-parallel','permute','permute-parallel',
		'transfer']
		if allowQ(&self.calculation.rt,blockMemoryType.blockCopyBlock)==0:
			print('need copy block')
			return self
		
		spaces = 1
		dims = 0
		cs = []
		xc = 0
		for space in range(SPACE):
			if self.field.f.canon[space].body != bodyType.nada:
				cp[space]= streams(self.field.f,division.copyVector,0,space)
				c1 = self.field.f.canon[space].count1Basis
				spaces *= c1	
				if xc < c1 :
					xc = c1
				dims += 1
		cs = [g.comps(xc)]
		ds = []
		for d in range(dims):
			ds += [g.dims(lattice = 1)]
		ds = [ds]
		ls = [1]
		bs = [g.bases('Sinc')]
		g.spaces(ls,cs,ds,bs).block(blocks)
		g.calculationInputs(RAMmax = 4,numVectors = 0,numNames = 0)
		g.fieldInputs(canonRank = 1,nStates = 1,OpIndex = 0)
		g.i()
		pt = streams(self.field.f,division.eigenVectors,0,0)
		self.field.f.name[int(division.copyVector)].Current[0] = 1
		for ii in range(spaces):
			iv = 1
			for space in range(SPACE):
				if self.field.f.canon[space].body != bodyType.nada:
					c1 = self.field.f.canon[space].count1Basis
					for c in range(c1):
						if c == (int(ii/iv)%c1):
							cp[space][c] = 1.0
						else:
							cp[space][c] = 0.0
					iv *= c1
			pt[ii] = tMatrixElements(0,self.field.f,division.copyVector,0,division.nullOverlap,0,vector,0)		
		return g
		
	def compress ( self, g : galaxy , vector : division = division.eigenVectors):
		"""self-> g
		
		testing...
		"""
		cdef inta spatial[SPACE][SPACE]
		for s in range(SPACE):
			for s2 in range(SPACE):
				spatial[s][s2] = 0
	
		spatial[0][0] = 1
		spatial[0][1] = 1
		spatial[0][2] = 1


		canonicalRankCompression(spatial,NULL,self.field.f,0,NULL,division.eigenVectors,0,1,0,g.field.f,1,vector,0,1,0,
		self.calculation.rt.TOLERANCE,
		self.calculation.rt.relativeTOLERANCE,
		self.calculation.rt.ALPHA,
		self.calculation.rt.XCONDITION,
		self.calculation.rt.MAX_CYCLE)
		return self