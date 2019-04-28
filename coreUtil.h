/*
 *  coreUtil.h
 *
 *
 *  Copyright 2019 Jonathan Jerke and Bill Poirier.
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
double lattice ( struct field * f1, INT_TYPE space );
INT_TYPE vector1Len(struct field *f1, INT_TYPE space);
void length1(struct field *f1, INT_TYPE *len);
INT_TYPE spaces( struct field * f1, enum division label);
enum division name ( struct field * f1, enum division label);
INT_TYPE part ( struct field * f1 , enum division label );
INT_TYPE species ( struct field * f1 , enum division label );
enum bodyType bodies ( struct field * f1 , enum division label);
enum bodyType Bodies ( struct field * f1 , enum division label,INT_TYPE space);
enum particleType particle ( struct field * f1 , enum division label ,INT_TYPE space);
INT_TYPE header ( struct field * f1 , enum division label );
INT_TYPE sizeofDivision(struct field *f1, enum division head, INT_TYPE space );
INT_TYPE defineSpinors (struct field *f1 );
INT_TYPE vectorLen(struct field *f1, INT_TYPE space);
INT_TYPE matrixLen(struct field *f1, enum bodyType body,INT_TYPE space);
INT_TYPE length ( struct field * f1 , enum division label, INT_TYPE *lens );
void assignParticle(struct field * f1, enum division ma, enum particleType pa , enum bodyType ba );
void assignOneWithPointers( struct field *f1, enum division oneMat, enum particleType particle );
void assignTwoWithPointers( struct field *f1, enum division twoMat );
INT_TYPE outerVectorLen(struct field *f1, enum bodyType bd, INT_TYPE space);
INT_TYPE alloc ( struct field * f1 , enum division label,INT_TYPE space );
INT_TYPE zero ( struct field * f1 , enum division label, INT_TYPE spin );
INT_TYPE myZero ( struct field * f1 , enum division label, INT_TYPE spin );
double traceOne( struct field * f1 , enum division label , INT_TYPE spin );
double tTrace( struct field * f1 , enum division label );
INT_TYPE tClear ( struct field * f1 , enum division label );
INT_TYPE CanonicalRank( struct field * f1 , enum division label , INT_TYPE spin );
INT_TYPE spins ( struct field * f1 , enum division label );
double getPosition(struct field * f1, INT_TYPE at , INT_TYPE space );
double sumSquare (struct field * f1,  enum division alloy);
INT_TYPE tReplace( struct field *f1 , enum division label,INT_TYPE spin,INT_TYPE space,INT_TYPE l );
void  fromBeginning( struct field * f1 ,enum division new, enum division head );
Stream_Type* myStreams ( struct field * f1, enum division label ,INT_TYPE spin );
Stream_Type* streams ( struct field * f1, enum division label ,INT_TYPE spin, INT_TYPE space );
double levelDetermine ( INT_TYPE M1 , double * array ,double level);
void assignView(INT_TYPE lane, struct field * f1, enum division A,INT_TYPE part );
void assignViewBlock(INT_TYPE lane, struct field * f1,  enum division A );
enum division ocean(INT_TYPE lane, struct field * f1,  INT_TYPE l, INT_TYPE spin);
double tEqua ( struct field * f1 , enum division targ ,INT_TYPE tspin, enum division orig,INT_TYPE ospin );
INT_TYPE tScaleOne( struct field * f1, enum division label,INT_TYPE spin, double scalar );
INT_TYPE tScale( struct field * f1, enum division label, double scalar );
INT_TYPE tAddTw( struct field* f1 , enum division left, INT_TYPE lspin, enum division right , INT_TYPE rspin);
INT_TYPE tEquals( struct field * f1 , enum division left , enum division right);
INT_TYPE tAddTwo( struct field * f1 , enum division left , enum division right);
INT_TYPE tSumMatrices(struct field *f1, enum division sum ,INT_TYPE spin, enum division mat  );
INT_TYPE tPartialSumMatrices(struct field *f1, enum division sum , enum division mat ,INT_TYPE tGetType );
INT_TYPE tAlt(struct field * f1 , enum division label, INT_TYPE spin , INT_TYPE space1);
INT_TYPE tEnd(struct field * f1 , enum division label, INT_TYPE spin , INT_TYPE space1);
INT_TYPE tPauli ( struct field * f1  );
INT_TYPE tId ( struct field *f1 , enum division label,INT_TYPE spin );
INT_TYPE tBoot ( struct field *f1 , enum division label,INT_TYPE spin );
double vectorElement (struct field * f1, enum division state, INT_TYPE l1,INT_TYPE l2 , INT_TYPE l3 );
double matrixElement (struct field * f1, enum division label, INT_TYPE i , INT_TYPE i2, INT_TYPE j,INT_TYPE j2, INT_TYPE k , INT_TYPE k2 );
void vectorArray (struct field * f1, enum division oneVector, enum division array,INT_TYPE M1);
void nuclearArray (struct field * f1,  enum division array,INT_TYPE M1);
INT_TYPE assignCores(struct field * f1, INT_TYPE parallel );
INT_TYPE Rank( struct field * f1 , enum division label );
double volume ( struct field * f1 );
INT_TYPE xAddTw( struct field* f1 , enum division left, INT_TYPE lspin,struct field* f2 ,  enum division right , INT_TYPE rspin);
void xsAdd (double scalar ,  INT_TYPE dim ,struct field * f1 , enum division targ ,INT_TYPE tspin,struct field * f2 , enum division orig,INT_TYPE o,INT_TYPE ospin );
double xEqua ( struct field * f1 , enum division targ ,INT_TYPE tspin,struct field * f2 , enum division orig,INT_TYPE ospin );
double xTwoBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
double xThreeBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
double xFourBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
double xOneBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
void printVectorAllocations(struct field *f1);
struct basisElement grabBasis (struct field * f1, INT_TYPE component, INT_TYPE particle, INT_TYPE elementIndex);
INT_TYPE  countLinesFromFile(struct calculation *c1, INT_TYPE location);
#endif /* coreUtil_h */
