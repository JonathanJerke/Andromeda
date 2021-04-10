###-------------###
######input.pxd####
###-------------###


#*  Copyright 2021 Jonathan Jerke and Bill Poirier.
#*  Ongoing support for this program is coordinated through quantumgalaxies.org.
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
from constants cimport blockMemoryType
from constants cimport runTime
from constants cimport calculation
from constants cimport field

cdef extern from "../Object/input.h":
    inta resetA(   runTime *rt);
    inta allowQ(   runTime *  f1,   blockMemoryType a );
    inta blockA(   runTime *  f1,   blockMemoryType a );
    inta bootShell (inta argc , char * argv[],  calculation * c1,   field *f)
    inta readShell (inta argc , char * argv[],  calculation * c1,   field *f)
