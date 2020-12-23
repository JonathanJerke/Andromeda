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

from constants cimport floata
from constants cimport field
from constants cimport calculation

from Model cimport initCal
from Model cimport initField
from Model cimport iModel
from Model cimport fModel

from constants cimport dimensions_label
from constants cimport metric_label
from constants cimport function_label

from constants cimport basisElementType
from constants cimport componentType
from constants cimport functionType
from constants cimport metricType
from constants cimport bodyType

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
							self.field.f.canon[space].particle[body+1].origin -= particle.lattice *(self.field.f.canon[space].count1Basis-1)*particle.anchor					
						if len(dim)==1:
							self.field.f.canon[space].body = bodyType.one
						elif len(dim)==2:
							self.field.f.canon[space].body = bodyType.two
						elif len(dim)==3:
							self.field.f.canon[space].body = bodyType.three
						
						
					for (space,base) in enumerate(bases):
						self.field.f.canon[space].basis = base
						
		return self	
					
	def i(self):
		"""Initiate allocation, do not overwrite
		"""
		if self.field.f.bootedMemory == 1 :
			print("warning, already booted")
			return self
		iModel(&self.calculation, &self.field)
		return self
					
	def metric(self, funcDesc = 'Coulomb', intervalDesc = 'interval',
										 betas =[0,1],interval = 7, contr = 2):
		"""Metric definition by description
		
		Parameters
		----------
		funcDesc : str
		intervalDesc : str
		betas : [floata,floata]
			interval span or first float only
		interval : int
			CanonRank of operator
		contr : int
			Off diagonal measure
			
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
		