###-------------###
###coreUtil.pxd####
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

cdef extern from "../Object/coreUtil.h":
    inta vector1Len(  sinc_label f1, inta space)
    void length1(  sinc_label f1, inta *len)
    division name (   sinc_label   f1,   division label)
    inta part (   sinc_label f1 ,   division label )
    inta species (   sinc_label  f1 ,   division label )
    bodyType bodies (   sinc_label  f1 ,   division label)
    bodyType Bodies (   sinc_label  f1 ,   division label,inta space)
    inta header (   sinc_label f1 ,   division label )
    inta sizeofDivision(  sinc_label f1,   division head, inta space )
    inta vectorLen(  sinc_label f1, inta space)
    inta matrixLen(  sinc_label f1,   bodyType body,inta space)
    inta length (   sinc_label  f1 ,   division label, inta *lens )
    void assignParticle(  sinc_label  f1,   division ma, inta label ,   bodyType ba )
    void assignOneWithPointers(   sinc_label f1,   division oneMat, inta label )
    void assignTwoWithPointers(   sinc_label f1,   division twoMat )
    inta outerVectorLen(  sinc_label f1,   bodyType bd, inta space)
    inta alloc (   sinc_label  f1 ,   division label,inta space )
    inta pZero (   sinc_label  *f1 ,   division label, inta spin )
    inta zero (   sinc_label  f1 ,   division label, inta spin )
    inta myZero (   sinc_label  f1 ,   division label, inta spin )
    floata tTrace(   sinc_label f1 ,   division label )
    inta pClear (   sinc_label * f1 ,   division label )
    inta tClear (   sinc_label  f1 ,   division label )
    inta CanonicalRank(   sinc_label  f1 ,   division label , inta spin )
    inta CanonicalOperator(   sinc_label f1,   division label, inta spin )
    inta spins (   sinc_label  f1 ,   division label )
    inta tReplace(   sinc_label f1 ,   division label,inta spin,inta space,inta l )
    void  fromBeginning(   sinc_label  f1 ,  division new,   division head )
    inta sumTo2(  sinc_label f1,floata scalar, inta space,   blockType bl,   division mat,inta ms,   division sum,inta spin)
    inta sumTo3(  sinc_label f1,floata scalar, inta space,   blockType bl,   division mat,inta ms,   division sum,inta spin)
    inta sumTo4(  sinc_label f1,floata scalar, inta space,   blockType bl,   division mat,inta ms,   division sum,inta spin)
    floata* myStreams (   sinc_label f1,   division label ,inta spin )
    floata* pMyStreams (   sinc_label *f1,   division label ,inta spin )
    floata* streams (   sinc_label f1,   division label ,inta spin, inta space )
    floata* pStreams (   sinc_label *f1,   division label ,inta spin, inta space )
    inta diagonalOp(  bodyType bd,  inta act,   blockType op,   blockType bl, inta N1,floata * vector, floata * toep, floata* vectorOut)
    division anotherLabel(  sinc_label *f1,   inta particle,  bodyType body)
    void assignSplit(   sinc_label f1,   division twoMat ,inta len,   division oneMat,   division bufcp )
    inta topezOp(floata origin, floata lattice,  bodyType bd,inta act,   blockType tv,   blockType bl,  inta N1,floata * vector , inta pw, floata * vectorOut)
    void assignView(inta lane,   sinc_label  f1,   division A,inta part )
    void assignViewBlock(inta lane,   sinc_label  f1,    division A )
    division ocean(inta lane,   sinc_label f1,  inta l, inta spin)
    floata tEqua (   sinc_label f1 ,   division targ ,inta tspin,   division orig,inta ospin )
    inta pScaleOne(   sinc_label *f1,   division label,inta spin, floata scalar )
    inta tScaleOne(   sinc_label  f1,   division label,inta spin, floata scalar )
    #inta tScale(   sinc_label  f1,   division label, DCOMPLEX scalar )
    inta tAddTw(   sinc_label f1 ,   division left, inta lspin,   division right , inta rspin)
    inta tEquals(   sinc_label  f1 ,   division left ,   division right)
    inta tAddTwo(   sinc_label f1 ,   division left ,   division right)
    inta tSumMatrices(  sinc_label f1,   division sum ,inta spin,   division mat  )
    inta tPartialSumMatrices(  sinc_label f1,   division sum ,   division mat ,inta tGetType )
    inta tAlt(  sinc_label  f1 ,   division label, inta spin , inta space1,   bodyType body)
    inta tEnd(  sinc_label f1 ,   division label, inta spin , inta space1,   bodyType body)
    inta tPauli (   sinc_label f1  )
    inta tId (   sinc_label f1 ,   division label,inta spin )
    inta pBoot (   sinc_label *f1 ,   division label,inta spin )
    inta tBoot (   sinc_label f1 ,   division label,inta spin ,floata scale)
    floata vectorElement (  sinc_label f1,   division state, inta l1,inta l2 , inta l3 )
    floata matrixElement (  sinc_label  f1,   division label, inta i , inta i2, inta j,inta j2, inta k , inta k2 )
    inta assignCores(  sinc_label  f1, inta parallel )
    inta defineCores(  calculation * c,   field * f)
    inta Rank(   sinc_label  f1 ,   division label )
    floata volume (   input * f1 )
    inta xAddTw(   sinc_label f1 ,   division left, inta lspin,  sinc_label f2 ,    division right , inta rspin)
    void xsAdd (floata scalar ,  inta dim ,  sinc_label  f1 ,   division targ ,inta tspin,  sinc_label  f2 ,   division orig,inta o,inta ospin )
    void xsEqu (floata scalar ,  inta dim ,  sinc_label  f1 ,   division targ ,inta t,inta tspin,inta dim2,  sinc_label  f2 ,   division orig,inta o,inta ospin )
    inta ready (   sinc_label f1 )
    inta bootedQ (   sinc_label f1)
    floata traceOne(   sinc_label  f1 ,   division label , inta spin )
    division defSpiralVector(   sinc_label *f1, inta term,   division ket)
    division defSpiralMatrix(   sinc_label *f1,   division H)
    division defSpiralGrid(   sinc_label *f1,   division bra, inta term, floata diagonalPreference)
    division defRefVector(   sinc_label *f1, inta spiralOp,   division ket)
    inta zeroSpiraly(   sinc_label f1,   division spiral)
    floata xEqua (   sinc_label  f1 ,   division targ ,inta tspin,  sinc_label  f2 ,   division orig,inta ospin )
    floata xOneBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label  f2,   division out,inta s2, inta periodic)
    floata xTwoBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label  f2,   division out,inta s2, inta periodic)
    floata xThreeBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label  f2,   division out,inta s2, inta periodic)
    floata xFourBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label  f2,   division out,inta s2, inta periodic)
    basisElement_label grabBasis (  sinc_label  f1, inta space, inta particle, inta elementIndex)
    basisElement_label transformBasis( inta flip,floata scale,   basisElement_label ba )
    inta  countLinesFromFile(   calculation *c1,  field f1,inta location, inta * ir,inta *ix)
    inta defineTerms(  calculation * c,   field *f,   division head, inta memory)
    inta InvertOp(  bodyType bd,inta invert, inta N1,floata * vector, floata* vectorOut)
    inta balance (  sinc_label f1,    division alloy, inta spin)
    void linkDetails(  sinc_label f1,   division linkHeader)
    void chainDetails(  sinc_label f1,   division chainHeader)
    void loopDetails(  sinc_label f1,   division loopHeader)


    void SG ( sinc_label f1, division vector , inta spin, floata amplitude, inta * gamma )
    void GTO( sinc_label f1, division vector , inta spin, floata amplitude, inta * gamma, floata *delta )



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
