###-------------###
######mAls.pxd#####
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

cdef extern from "../Object/mAls.h":
    floata canonicalRankDecomposition( sinc_label  f1 , floata * cofact,inta G,floata *GG, division origin,inta l1,inta l2,inta os, inta neo,division alloy ,inta l3 , inta l4,  inta spin ,floata tolerance,floata relativeTolerance, floata condition,floata maxCondition, inta maxCycle)
    floata CanonicalRankDecomposition ( sinc_label  f1 , floata * cofact, division origin,inta os, division alloy,inta spin,  floata tolerance ,  floata relativeTolerance, floata condition,floata threshold, inta maxCycle ,floata maxCondition, inta canon,inta X1 )
    floata printExpectationValues ( calculation *c,   sinc_label  f1 , division ha  , division vector)
    floata tMatrixElements ( inta rank,  sinc_label  f1 ,   division bra, inta bspin,  division mat, inta mspin,  division ket, inta kspin )
    inta tOuterProductSu( sinc_label  f1,  division vector , inta a, division vector2,inta b, division proj, inta c)
    floata pMatrixElement (   sinc_label  f1 ,   division alloy1 , inta spin1,division op, inta ospin, division alloy2 , inta spin2)
    inta tOuterProductSuOne( sinc_label  f1,inta space, division vector , inta a,   division vector2,inta b,   division proj, inta c)
    inta tGEMV (inta rank, sinc_label  f1,inta dim, division equals,inta e, inta espin, division left,inta l,inta lspin, division right, inta r,inta rspin )
    inta tGEVV (inta rank,  sinc_label  f1,inta dim, division equals,inta e, inta espin, division left,inta l,inta lspin, division right, inta r,inta rspin )
    floata tDOT (inta rank, sinc_label  f1,inta dim,char leftChar, division left,inta l,inta lspin, char rightChar, division right, inta r,inta rspin )
    inta tHX(  inta rank, sinc_label f1 ,division left, inta l, inta im, floata prod, division ket , inta k, inta sp2,   division oket, inta o,inta ospin )
    void tHXpY (   sinc_label f1 ,  division bra,   division left,inta shiftFlag,  division right ,  floata tolerance ,floata relativeTolerance,floata condition,floata threshold, inta maxCycle, floata maxCondition, inta canon,inta X1)
