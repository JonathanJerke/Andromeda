###-------------###
#####saUtil.pxd####
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
from constants cimport genusType
from constants cimport bodyType
from constants cimport blockType
from constants cimport division
from constants cimport basisElement_label
from constants cimport sinc_label
from constants cimport input
from constants cimport field
from constants cimport calculation

cdef extern from "../Object/saUtil.h":
    void testSAAgain (   sinc_label f1 ,   division vector )
    inta tPerms(  bodyType bd)
    void tTestSA (  bodyType bd, inta n)
    inta tSize(  bodyType bd)
    inta tPaths(  bodyType bd , inta irrep )
    inta nEqua(  bodyType bd, inta *a )
    void tIntr (   bodyType bd , inta eq , floata * a)
    inta tInnerTest(   sinc_label f1,   division A ,  division B)
    inta tSA (  bodyType bd, inta X, inta Y, inta Z, inta T )
    floata tGetProjection(   bodyType bd, inta type , inta perm )
    floata tGetVector(  bodyType bd , inta type , inta perm )
    inta tClassifyComponents(   sinc_label  , floata * up, floata * entropy )
    inta tClassify( sinc_label  ,   division label)
    inta tBuildIrr ( inta rank,   sinc_label , inta meta ,   division origin, inta ospin,   division targ , inta tspin)
    inta tPermuteOne(inta rank,   sinc_label , inta dim, inta leftChar ,   division left, inta l, inta lspin,   division equals,inta e, inta espin)
    inta tPermute(inta rank,   sinc_label , inta leftChar ,   division left, inta lspin,   division equals, inta espin)
    inta tAllCompPermMultiplyMP( sinc_label  f1 ,   division left ,inta lspin,   division right ,inta rspin, floata * sequ)
    bodyType commandSA(  bodyType bd, inta act,   blockType cl,   blockType bl,inta perm[], inta op[])
    inta tAddUpComponents( inta rank,   sinc_label  f1 ,   division left ,   division right ,  floata *up)
    inta tTabulateInnerProjection( sinc_label  f1 ,   division vec, floata *up)
    floata tTestInner (   bodyType bd, inta i , inta j )
    inta tTest (   bodyType bd )
    inta tSizeUp(inta rank,   sinc_label  f1 , inta type,   division label)
    inta testSA (   sinc_label f1 ,   division vector )
    inta irreps (   bodyType bd, inta type )
