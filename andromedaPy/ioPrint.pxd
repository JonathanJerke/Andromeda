###-------------###
####ioPrint.pxd####
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
from constants cimport mea
from constants cimport genusType
from constants cimport bodyType
from constants cimport blockType
from constants cimport division
from constants cimport basisElement_label
from constants cimport sinc_label
from constants cimport input
from constants cimport field
from constants cimport calculation
from libc.stdio cimport FILE

cdef extern from "../Object/ioPrint.h":
    inta writeFast( sinc_label f1,char * filename, inta space, division label ,inta spin)
    inta readFast ( sinc_label f1,char * filename, inta command ,inta space, division label ,inta spin, inta space2)
    void tFilename (char * cycleName, inta count, inta body ,inta IRREP, inta cmpl, char * filename)
    inta printOut(  calculation *c,  field f1, inta reset,inta lv,  division vector )
    inta printVector (  calculation *c,  sinc_label f1,char * name,char * vectorName,inta iv, inta irrep, mea * vector)
    double evaluateDensityBracket( double x [], size_t dim , void * params )
    double evaluateVectorBracket( double x [], size_t dim , void * params )
    inta tSymmetryClass (char * cycleName, char * read , char * filename,inta cmplFlag, inta cmpl)
    void tFromReadToToken (char * read , char * token)
    void outputFormat(  sinc_label   f1, FILE * out,   division output ,inta spin)
    inta inputFormat(  sinc_label f1,char * name,    division buffer, inta input)
    inta printOutput (   sinc_label f1,inta number)
    inta printVectorOutput (   sinc_label f1,inta number)
    inta printFaceOutput (   sinc_label f1,inta number)
    inta ioArray(  calculation *c1,   field f,char * name,inta N1, floata * matrix, inta ioIn)
    inta tLoadEigenWeights (   calculation * c1,  field *f1, char * filename ,inta *ct,   division input, inta collect)
    void tFilename (char * cycleName, inta count, inta body ,inta IRREP, inta cmpl, char * filename)
    mea tFromReadToFilename (char * cycleName, char * read , char * filename,inta cmplFlag, inta cmpl,char * title,inta *number)
