###-------------###
#Decompose.pxd#####
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
from constants cimport division
from constants cimport sinc_label
from constants cimport field
from constants cimport calculation

cdef extern from "../Object/Compression.h":
    floata canonicalRankCompression( inta ** spatial, floata * cofact,sinc_label  f1 ,inta G,floata *GG, division origin,inta l1,inta l2,inta os, sinc_label  f2 ,inta neo,division alloy ,inta l3 , inta l4, inta spin ,floata tolerance,floata relativeTolerance, floata condition,floata maxCondition, inta maxCycle);

    floata CanonicalRankCompression ( inta ** spatial,  floata * cofact, sinc_label  f1 , division origin,inta os, sinc_label  f2 , division alloy,inta spin,  floata tolerance ,  floata relativeTolerance, floata condition,floata threshold, inta maxCycle ,floata maxCondition, inta canon,inta X1 );
