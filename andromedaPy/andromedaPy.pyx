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

cimport andromedaPy
from constants cimport floata
from constants cimport field
from constants cimport calculation
from Model cimport initCal
from Model cimport initField
from Model cimport iModel
from Model cimport fModel
from constants cimport dimensions_label
from constants cimport basisElementType
from constants cimport componentType

cdef enum componentType:
    nullComponent,
    spatialComponent1,
    spatialComponent2,
    spatialComponent3,

cdef class base:
    cdef calculation calculation
    cdef field field

    def __cinit__(self):
        self.calculation = initCal()
        self.field = initField()

    def __dealloc__(self):
        fModel(&self.field.f)
    
    def dims(self, floata lattice, floata attack=0.5, floata origin =0.0,floata anchor =0.5 ):
        """Returns new floata parameters for a linear dimension.
        """
        return dimensions_label(lattice = lattice , attack = attack, origin = origin, anchor = anchor )

    def comps(self, i=1,periodic = False):
        """Returns new component enum a linear dimension.
        """
        if i <= 0 :
            return componentType.nullComponent
        if not periodic:
            if i == 1 :
                return componentType.spatialComponent1
            elif i == 2 :
                return componentType.spatialComponent2
            elif i == 3 :
                return componentType.spatialComponent3
                
        if periodic:
            if i == 1 :
                return componentType.periodicComponent1
            elif i == 2 :
                return componentType.periodicComponent2
            elif i == 3 :
                return componentType.periodicComponent3

        else :
            return componentType.nullComponent

    def comps(self,component, dims, count1Inc=1, basis=basisElementType.SincBasisElement):
        """Returns new info for parameters for a component bundle.
        """
        nbody = len(dims)
        if (nbody > 3) :
            print('warning: nbody > 3')
            return self
    
