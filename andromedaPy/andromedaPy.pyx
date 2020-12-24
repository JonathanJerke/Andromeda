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

from Model cimport initCal
from Model cimport initField
from Model cimport iModel
from Model cimport fModel
from input cimport blockA
from input cimport resetA
from ioPrint cimport tLoadEigenWeights
from ioPrint cimport printOut
from coreUtil cimport tBoot

from constants cimport dimensions_label
from constants cimport metric_label
from constants cimport function_label

from constants cimport basisElementType
from constants cimport componentType
from constants cimport functionType
from constants cimport metricType
from constants cimport bodyType
from constants cimport blockMemoryType


cdef class galaxy:
	cdef calculation calculation
	cdef field field

	def __cinit__(self):
		self.calculation = initCal()
		self.field = initField()

	def __dealloc__(self):
		fModel(&self.field.f)
    
	def dims(self, floata lattice, floata attack=0.5, floata origin =0.0,
													floata anchor =0.5 ):
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


	def bases(self, desc ):
		"""Returns enumeration of basis types
		"""
		names = dict({'Sinc':basisElementType.SincBasisElement,
					'Gaussian':basisElementType.GaussianBasisElement,
					'Dirac':basisElementType.DiracDeltaElement,
					'State':basisElementType.StateBasisElement,
					'overlap':basisElementType.overlapBasisElement})
		return names[desc]

	def spaces(self,labels, comps, dims, bases ):
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
			
	def block(self, blockDescs):
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
						
	def calculationInputs ( self, numNames=-1, numVectors=-1, shiftFlag=-1,Lambda=-1
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
		return self
		
	def fieldInputs( self, flex = -1, OpIndex = -2 , body =-1,irrep = -1, Iterations = -1
	,nStates = -1,iRank = -1,canonRank= -1,xRank = -1,qFloor = -1,filter = -1,collect=-1):
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
		
	def vectors( self ):
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
		
	def read_file ( self, filename,vector = division.eigenVectors, collect = 0 ):
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
	 	count = inta(0)
		tLoadEigenWeights (  &self.calculation, self.field ,filename.encode('utf-8'), 
				&count,  vector, collect)
		return count
		
		
	def to_file ( self, vector = division.eigenVectors, reset = 1, index = 1 ):
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
		
	def gaussian ( self, vector = division.eigenVectors, spin = 0, width = 1.0):
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
		
	def metric(self, funcDesc = 'Coulomb', intervalDesc = 'interval',
										 betas =[0,1],interval = 7, contr = 2):
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
		