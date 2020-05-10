/*
 *  coreUtil.h
 *
 *
 *  Copyright 2020 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University,
 *  the Robert A. Welch Foundation, and Army Research Office.
 *
 
 *   *   This file is part of Andromeda.
 
 *   *   Andromeda is free software: you can redistribute it and/or modify
 *   *   it under the terms of the GNU General Public License as published by
 *   *   the Free Software Foundation, either version 3 of the License, or
 *   *   (at your option) any later version.
 
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
double lattice ( struct input * f1, INT_TYPE space );
INT_TYPE vector1Len(struct sinc_label f1, INT_TYPE space);
void length1(struct sinc_label f1, INT_TYPE *len);
INT_TYPE spaces( struct sinc_label  f1, enum division label);
enum division name ( struct sinc_label   f1, enum division label);
INT_TYPE pPart ( struct sinc_label *f1 , enum division label );

INT_TYPE part ( struct sinc_label f1 , enum division label );
INT_TYPE species ( struct sinc_label  f1 , enum division label );
enum bodyType bodies ( struct sinc_label  f1 , enum division label);
enum bodyType Bodies ( struct sinc_label  f1 , enum division label,INT_TYPE space);
enum particleType particle ( struct sinc_label  f1 , enum division label ,INT_TYPE space);
INT_TYPE header ( struct sinc_label f1 , enum division label );
INT_TYPE sizeofDivision(struct sinc_label f1, enum division head, INT_TYPE space );
INT_TYPE vectorLen(struct sinc_label f1, INT_TYPE space);
INT_TYPE pVectorLen(struct sinc_label *f1, INT_TYPE space);

INT_TYPE matrixLen(struct sinc_label f1, enum bodyType body,INT_TYPE space);
INT_TYPE length ( struct sinc_label  f1 , enum division label, INT_TYPE *lens );
void assignParticle(struct sinc_label  f1, enum division ma, enum particleType pa , enum bodyType ba );
void assignOneWithPointers( struct sinc_label f1, enum division oneMat, enum particleType particle );
void assignTwoWithPointers( struct sinc_label f1, enum division twoMat );
INT_TYPE outerVectorLen(struct sinc_label f1, enum bodyType bd, INT_TYPE space);
INT_TYPE alloc ( struct sinc_label  f1 , enum division label,INT_TYPE space );
INT_TYPE pZero ( struct sinc_label  *f1 , enum division label, INT_TYPE spin );
INT_TYPE zero ( struct sinc_label  f1 , enum division label, INT_TYPE spin );
INT_TYPE myZero ( struct sinc_label  f1 , enum division label, INT_TYPE spin );
double tTrace( struct sinc_label f1 , enum division label );
INT_TYPE pClear ( struct sinc_label * f1 , enum division label );

INT_TYPE tClear ( struct sinc_label  f1 , enum division label );
INT_TYPE CanonicalRank( struct sinc_label  f1 , enum division label , INT_TYPE spin );
INT_TYPE CanonicalOperator( struct sinc_label f1, enum division label, INT_TYPE spin );
INT_TYPE spins ( struct sinc_label  f1 , enum division label );
double sumSquare (struct sinc_label  f1,  enum division alloy);
INT_TYPE tReplace( struct sinc_label f1 , enum division label,INT_TYPE spin,INT_TYPE space,INT_TYPE l );
void  fromBeginning( struct sinc_label  f1 ,enum division new, enum division head );

Stream_Type* myStreams ( struct sinc_label f1, enum division label ,INT_TYPE spin );
Stream_Type* pMyStreams ( struct sinc_label *f1, enum division label ,INT_TYPE spin );
Stream_Type* streams ( struct sinc_label f1, enum division label ,INT_TYPE spin, INT_TYPE space );
Stream_Type* pStreams ( struct sinc_label *f1, enum division label ,INT_TYPE spin, INT_TYPE space );
INT_TYPE diagonalOp(enum bodyType bd,  INT_TYPE act, enum block op, enum block bl, INT_TYPE N1,Stream_Type * vector, Stream_Type * toep, Stream_Type* vectorOut);
enum division anotherLabel(struct sinc_label *f1, enum particleType particle,enum bodyType body);
void assignSplit( struct sinc_label f1, enum division twoMat ,INT_TYPE len, enum division oneMat, enum division bufcp );
INT_TYPE topezOp(enum bodyType bd,INT_TYPE act, enum block tv, enum block bl,  INT_TYPE N1,Stream_Type * vector , INT_TYPE pw, Stream_Type * vectorOut);
double levelDetermine ( INT_TYPE M1 , double * array ,double level);
void assignView(INT_TYPE lane, struct sinc_label  f1, enum division A,INT_TYPE part );
void assignViewBlock(INT_TYPE lane, struct sinc_label  f1,  enum division A );
enum division ocean(INT_TYPE lane, struct sinc_label f1,  INT_TYPE l, INT_TYPE spin);
double tEqua ( struct sinc_label f1 , enum division targ ,INT_TYPE tspin, enum division orig,INT_TYPE ospin );
INT_TYPE pScaleOne( struct sinc_label *f1, enum division label,INT_TYPE spin, double scalar );
INT_TYPE tScaleOne( struct sinc_label  f1, enum division label,INT_TYPE spin, double scalar );
INT_TYPE tScale( struct sinc_label  f1, enum division label, DCOMPLEX scalar );
INT_TYPE tAddTw( struct sinc_label f1 , enum division left, INT_TYPE lspin, enum division right , INT_TYPE rspin);
INT_TYPE tEquals( struct sinc_label  f1 , enum division left , enum division right);
INT_TYPE tAddTwo( struct sinc_label f1 , enum division left , enum division right);
INT_TYPE tSumMatrices(struct sinc_label f1, enum division sum ,INT_TYPE spin, enum division mat  );
INT_TYPE tPartialSumMatrices(struct sinc_label f1, enum division sum , enum division mat ,INT_TYPE tGetType );
INT_TYPE tAlt(struct sinc_label  f1 , enum division label, INT_TYPE spin , INT_TYPE space1);
INT_TYPE tEnd(struct sinc_label f1 , enum division label, INT_TYPE spin , INT_TYPE space1);
INT_TYPE tPauli ( struct sinc_label f1  );
INT_TYPE tId ( struct sinc_label f1 , enum division label,INT_TYPE spin );
INT_TYPE pBoot ( struct sinc_label *f1 , enum division label,INT_TYPE spin );
INT_TYPE tBoot ( struct sinc_label f1 , enum division label,INT_TYPE spin );
double vectorElement (struct sinc_label f1, enum division state, INT_TYPE l1,INT_TYPE l2 , INT_TYPE l3 );
double matrixElement (struct sinc_label  f1, enum division label, INT_TYPE i , INT_TYPE i2, INT_TYPE j,INT_TYPE j2, INT_TYPE k , INT_TYPE k2 );
void vectorArray (struct sinc_label f1, enum division oneVector, enum division array,INT_TYPE M1);
void pVectorArray (struct sinc_label *f1, enum division oneVector, enum division array,INT_TYPE M1);
void pNuclearArray (struct input c, struct field *f1,  enum division array,INT_TYPE M1);

void nuclearArray (struct input c, struct field f1,  enum division array,INT_TYPE M1);
INT_TYPE assignCores(struct sinc_label  f1, INT_TYPE parallel );
INT_TYPE defineCores(struct calculation * c, struct field * f);
INT_TYPE Rank( struct sinc_label  f1 , enum division label );
double volume ( struct input * f1 );
INT_TYPE xAddTw( struct sinc_label f1 , enum division left, INT_TYPE lspin,struct sinc_label f2 ,  enum division right , INT_TYPE rspin);
void xsAdd (double scalar ,  INT_TYPE dim ,struct sinc_label  f1 , enum division targ ,INT_TYPE tspin,struct sinc_label  f2 , enum division orig,INT_TYPE o,INT_TYPE ospin );
void xsEqu (double scalar ,  INT_TYPE dim ,struct sinc_label  f1 , enum division targ ,INT_TYPE t,INT_TYPE tspin,INT_TYPE dim2,struct sinc_label  f2 , enum division orig,INT_TYPE o,INT_TYPE ospin );
enum division defSpiralVector( struct sinc_label *f1, INT_TYPE term, enum division ket);
enum division defSpiralMatrix( struct sinc_label *f1, enum division H);
enum division defSpiralGrid( struct sinc_label *f1, enum division bra, INT_TYPE term, double diagonalPreference);
enum division defRefVector( struct sinc_label *f1, INT_TYPE spiralOp, enum division ket);
INT_TYPE zeroSpiraly( struct sinc_label f1, enum division spiral);
double xEqua ( struct sinc_label  f1 , enum division targ ,INT_TYPE tspin,struct sinc_label  f2 , enum division orig,INT_TYPE ospin );
double xTwoBand (struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct sinc_label  f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
double xThreeBand (struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct sinc_label  f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
double xFourBand (struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct sinc_label  f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
double xOneBand (struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct sinc_label  f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
void printVectorAllocations(struct sinc_label f1);
struct basisElement grabBasis (struct sinc_label  f1, INT_TYPE space, INT_TYPE particle, INT_TYPE elementIndex);
struct basisElement defineSincBasis (enum noteType note, enum componentType space, enum basisElementType basis, double lattice , double origin, INT_TYPE count1, INT_TYPE elementIndex );
struct basisElement defineGaussBasis (enum noteType note, enum componentType space, enum basisElementType basis, double lattice , double origin, INT_TYPE count1, INT_TYPE elementIndex );
struct basisElement defineSpinorBasis (enum noteType note, enum componentType space,INT_TYPE total, INT_TYPE elementIndex );
struct basisElement transformBasis( INT_TYPE flip,double scale, struct basisElement ba );
INT_TYPE  countLinesFromFile( struct calculation *c1,struct field f1,INT_TYPE location, INT_TYPE * ir,INT_TYPE *ix);
INT_TYPE completeInverse (INT_TYPE rank, struct sinc_label  f1, INT_TYPE dim,enum division vector,INT_TYPE v,INT_TYPE spin, enum division ov , INT_TYPE v2,INT_TYPE sp2);
INT_TYPE defineTerms(struct calculation * c, struct sinc_label *f1, enum division head, INT_TYPE memory);
INT_TYPE InvertOp(enum bodyType bd,INT_TYPE invert, INT_TYPE N1,Stream_Type * vector, Stream_Type* vectorOut);
#endif /* coreUtil_h */
