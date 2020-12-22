/**
 *  coreUtil.h
 *
 *
 *  Copyright 2021 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University,
 *  the Robert A. Welch Foundation, and the Army Research Office.
 *
 
 *   *   This file is part of Andromeda.
 
 *   *   Andromeda is free software: you can redistribute it and/or modify
 *   *   it under the terms of the GNU General Public License as published by
 *   *   the Free Software Foundation, either version 3 of the License.
 
 *   *   Andromeda is distributed in the hope that it will be useful,
 *   *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   *   GNU General Public License for more details.
 
 *   *   You should have received a copy of the GNU General Public License
 *   *   along with Andromeda.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef coreUtil_h
#define coreUtil_h

#include "constants.h"
#include "coreMath.h"
#include "coreForce.h"
#include "ioPrint.h"
inta vector1Len(  sinc_label f1, inta space);
void length1(  sinc_label f1, inta *len);
division name (   sinc_label   f1,   division label);
inta part (   sinc_label f1 ,   division label );
inta species (   sinc_label  f1 ,   division label );
bodyType bodies (   sinc_label  f1 ,   division label);
bodyType Bodies (   sinc_label  f1 ,   division label,inta space);
inta header (   sinc_label f1 ,   division label );
inta sizeofDivision(  sinc_label f1,   division head, inta space );
inta vectorLen(  sinc_label f1, inta space);
inta matrixLen(  sinc_label f1,   bodyType body,inta space);
inta length (   sinc_label  f1 ,   division label, inta *lens );
void assignParticle(  sinc_label  f1,   division ma, inta label ,   bodyType ba );
void assignOneWithPointers(   sinc_label f1,   division oneMat, inta label );
void assignTwoWithPointers(   sinc_label f1,   division twoMat );
inta outerVectorLen(  sinc_label f1,   bodyType bd, inta space);
inta alloc (   sinc_label  f1 ,   division label,inta space );
inta pZero (   sinc_label  *f1 ,   division label, inta spin );
inta zero (   sinc_label  f1 ,   division label, inta spin );
inta myZero (   sinc_label  f1 ,   division label, inta spin );
double tTrace(   sinc_label f1 ,   division label );
inta pClear (   sinc_label * f1 ,   division label );
inta tClear (   sinc_label  f1 ,   division label );
inta CanonicalRank(   sinc_label  f1 ,   division label , inta spin );
inta CanonicalOperator(   sinc_label f1,   division label, inta spin );
inta spins (   sinc_label  f1 ,   division label );
inta tReplace(   sinc_label f1 ,   division label,inta spin,inta space,inta l );
void  fromBeginning(   sinc_label  f1 ,  division new,   division head );
inta sumTo2(  sinc_label f1,double scalar, inta space,   blockType bl,   division mat,inta ms,   division sum,inta spin);
inta sumTo3(  sinc_label f1,double scalar, inta space,   blockType bl,   division mat,inta ms,   division sum,inta spin);
inta sumTo4(  sinc_label f1,double scalar, inta space,   blockType bl,   division mat,inta ms,   division sum,inta spin);
floata* myStreams (   sinc_label f1,   division label ,inta spin );
floata* pMyStreams (   sinc_label *f1,   division label ,inta spin );
floata* streams (   sinc_label f1,   division label ,inta spin, inta space );
floata* pStreams (   sinc_label *f1,   division label ,inta spin, inta space );
inta diagonalOp(  bodyType bd,  inta act,   blockType op,   blockType bl, inta N1,floata * vector, floata * toep, floata* vectorOut);
division anotherLabel(  sinc_label *f1,   inta particle,  bodyType body);
void assignSplit(   sinc_label f1,   division twoMat ,inta len,   division oneMat,   division bufcp );
inta topezOp(double origin, double lattice,  bodyType bd,inta act,   blockType tv,   blockType bl,  inta N1,floata * vector , inta pw, floata * vectorOut);
void assignView(inta lane,   sinc_label  f1,   division A,inta part );
void assignViewBlock(inta lane,   sinc_label  f1,    division A );
division ocean(inta lane,   sinc_label f1,  inta l, inta spin);
double tEqua (   sinc_label f1 ,   division targ ,inta tspin,   division orig,inta ospin );
inta pScaleOne(   sinc_label *f1,   division label,inta spin, double scalar );
inta tScaleOne(   sinc_label  f1,   division label,inta spin, double scalar );
inta tScale(   sinc_label  f1,   division label, DCOMPLEX scalar );
inta tAddTw(   sinc_label f1 ,   division left, inta lspin,   division right , inta rspin);
inta tEquals(   sinc_label  f1 ,   division left ,   division right);
inta tAddTwo(   sinc_label f1 ,   division left ,   division right);
inta tSumMatrices(  sinc_label f1,   division sum ,inta spin,   division mat  );
inta tPartialSumMatrices(  sinc_label f1,   division sum ,   division mat ,inta tGetType );
inta tAlt(  sinc_label  f1 ,   division label, inta spin , inta space1,   bodyType body);
inta tEnd(  sinc_label f1 ,   division label, inta spin , inta space1,   bodyType body);
inta tPauli (   sinc_label f1  );
inta tId (   sinc_label f1 ,   division label,inta spin );
inta pBoot (   sinc_label *f1 ,   division label,inta spin );
inta tBoot (   sinc_label f1 ,   division label,inta spin ,floata scale);
double vectorElement (  sinc_label f1,   division state, inta l1,inta l2 , inta l3 );
double matrixElement (  sinc_label  f1,   division label, inta i , inta i2, inta j,inta j2, inta k , inta k2 );
inta assignCores(  sinc_label  f1, inta parallel );
inta defineCores(  calculation * c,   field * f);
inta Rank(   sinc_label  f1 ,   division label );
double volume (   input * f1 );
inta xAddTw(   sinc_label f1 ,   division left, inta lspin,  sinc_label f2 ,    division right , inta rspin);
void xsAdd (double scalar ,  inta dim ,  sinc_label  f1 ,   division targ ,inta tspin,  sinc_label  f2 ,   division orig,inta o,inta ospin );
void xsEqu (double scalar ,  inta dim ,  sinc_label  f1 ,   division targ ,inta t,inta tspin,inta dim2,  sinc_label  f2 ,   division orig,inta o,inta ospin );
inta ready (   sinc_label f1 );
inta bootedQ (   sinc_label f1);
double traceOne(   sinc_label  f1 ,   division label , inta spin );
division defSpiralVector(   sinc_label *f1, inta term,   division ket);
division defSpiralMatrix(   sinc_label *f1,   division H);
division defSpiralGrid(   sinc_label *f1,   division bra, inta term, double diagonalPreference);
division defRefVector(   sinc_label *f1, inta spiralOp,   division ket);
inta zeroSpiraly(   sinc_label f1,   division spiral);
double xEqua (   sinc_label  f1 ,   division targ ,inta tspin,  sinc_label  f2 ,   division orig,inta ospin );
double xOneBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label  f2,   division out,inta s2, inta periodic);
double xTwoBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label  f2,   division out,inta s2, inta periodic);
double xThreeBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label  f2,   division out,inta s2, inta periodic);
double xFourBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label  f2,   division out,inta s2, inta periodic);
  basisElement_label grabBasis (  sinc_label  f1, inta space, inta particle, inta elementIndex);
  basisElement_label transformBasis( inta flip,double scale,   basisElement_label ba );
inta  countLinesFromFile(   calculation *c1,  field f1,inta location, inta * ir,inta *ix);
inta defineTerms(  calculation * c,   field *f,   division head, inta memory);
inta InvertOp(  bodyType bd,inta invert, inta N1,floata * vector, floata* vectorOut);
inta balance (  sinc_label f1,    division alloy, inta spin);
void linkDetails(  sinc_label f1,   division linkHeader);
void chainDetails(  sinc_label f1,   division chainHeader);
void loopDetails(  sinc_label f1,   division loopHeader);
#endif /* coreUtil_h */
