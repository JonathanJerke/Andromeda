###-------------###
#interfaceMath.pxd#
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
from constants cimport blockType
from constants cimport division
from constants cimport basisElement_label
from constants cimport sinc_label
from constants cimport input
from constants cimport field
from constants cimport calculation

cdef extern from "../Object/interfaceMath.h":
    inta tdsyev( inta rank,   char job , inta n, floata * ar, inta ns , floata * w )
    inta tdsygv( inta rank,   sinc_label f1, char job , inta n, floata * sr, floata * ar, inta ns , floata * w )
    inta tdgeqr( inta rank,   sinc_label f1,inta len, inta n, floata * ar, inta ns ,floata *w, floata *xr , inta xs )
    #inta tzheev( inta rank,   sinc_label f1, char job , inta n,DCOMPLEX * ar, inta ns , floata * w )
    #inta tzhegv( inta rank,   sinc_label f1, char job , inta n,DCOMPLEX * sr, DCOMPLEX * ar, inta ns , floata * w )
    void transpose(inta N, inta M, floata * orig, floata* targ)
    floata Sinc( floata d , floata x)
    #DCOMPLEX periodicSinc ( floata d , floata x, floata momentum, inta N1 )
    floata SS( floata d1 , floata x , floata d2, floata y )
    inta tdpotrf ( inta L1, floata * array ,inta LS1) ;
    inta tdpotrs ( inta L1, inta M2, floata * array,inta LS1, floata * arrayo,inta inc)
    floata tdpocon (inta rank,  sinc_label f1,  inta L1 , floata * Matrix ,inta stride)
    inta tdgels ( inta rank,  sinc_label f1 , inta L1, inta M2, floata * array, floata * arrayo ,inta inc)
    inta tdgesvd ( inta rank,   sinc_label f1 ,  inta M1, inta M2, floata * ge, floata * m1, floata* m2 )
    inta tInverse(   sinc_label f1, inta n, floata * ar)
    #inta tzInverse(   sinc_label f1, inta n, DCOMPLEX * ar)
    inta tLowdin( inta n , floata *ar, floata *lowdinVec, floata * lowdinMatrix )
