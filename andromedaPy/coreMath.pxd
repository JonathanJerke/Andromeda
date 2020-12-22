###-------------###
###coreMath.pxd####
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

include "system.pxi"


from constants cimport inta
from constants cimport floata
from constants cimport genusType
from constants cimport bodyType
from constants cimport division
from constants cimport sinc_label
from constants cimport field
from constants cimport calculation

cdef extern from "../Object/coreMath.h":
    floata SymmetrizedGaussianInSinc( floata K, inta n , inta m , floata X )
    floata Power ( floata b, inta n )
    #extern floata complex Faddeeva_erfcx(floata complex z, floata relerr)
    #DCOMPLEX expErf ( DCOMPLEX z )
    floata momentumIntegralInTrain ( floata beta, floata kl , floata d,  genusType hidden,   bodyType body )
    floata delta ( inta n )
    inta sign( inta n )
    inta imax( inta x1, inta x2 )
    inta imin( inta x1, inta x2 )
    floata max( floata x1, floata x2 )
    floata min( floata x1, floata x2 )
