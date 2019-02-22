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
INT_TYPE spaces( struct field * f1, enum division label);
INT_TYPE name ( struct field * f1, enum division label);
INT_TYPE part ( struct field * f1 , enum division label );
INT_TYPE property ( struct field * f1, enum division label);
INT_TYPE species ( struct field * f1 , enum division label );
enum body bodies ( struct field * f1 , enum division label );
enum bodyType bodyType ( struct field * f1 , enum division label );
INT_TYPE header ( struct field * f1 , enum division label );
INT_TYPE purpose ( struct field * f1, enum division label);
INT_TYPE memory ( struct field * f1, enum division label);
INT_TYPE tPath ( struct field * f1, enum division label);

void initPointerTensors(struct field *f1);
INT_TYPE defineSpinors (struct field *f1 );
INT_TYPE * vectorLen ( struct field * f1 , enum division label );
INT_TYPE length ( struct field * f1 , enum division label, INT_TYPE *lens );
INT_TYPE alloc ( struct field * f1 , enum division label );
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
ADDRESS_TYPE  fromBegining( struct field * f1 , enum division label );
ADDRESS_TYPE  fromMyBegining( struct field * f1 , enum division label  );
Stream_Type* myStreams ( struct field * f1, enum division label ,INT_TYPE spin );
Stream_Type* streams ( struct field * f1, enum division label ,INT_TYPE spin, INT_TYPE space );
double levelDetermine ( INT_TYPE M1 , double * array ,double level);
enum division rivers(INT_TYPE rank, struct field * f1, enum division A, INT_TYPE category);
INT_TYPE riversBranch(INT_TYPE rank, struct field * f1, enum division A, INT_TYPE spin, INT_TYPE category);
enum division ocean(INT_TYPE rank, struct field * f1, enum division A, INT_TYPE l, INT_TYPE spin);
double tEqua ( struct field * f1 , enum division targ ,INT_TYPE tspin, enum division orig,INT_TYPE ospin );
INT_TYPE tScaleOne( struct field * f1, enum division label,INT_TYPE spin, double scalar );
INT_TYPE tScale( struct field * f1, enum division label, double scalar );
INT_TYPE tAddTw( struct field* f1 , enum division left, INT_TYPE lspin, enum division right , INT_TYPE rspin);
INT_TYPE tEquals( struct field * f1 , enum division left , enum division right);
INT_TYPE tAddTwo( struct field * f1 , enum division left , enum division right);
INT_TYPE tSumMatrices(struct field *f1, enum division sum , enum division mat  );
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
INT_TYPE xAddTw( struct field* f1 , enum division left, INT_TYPE lspin,struct field* f2 ,  enum division right , INT_TYPE rspin);
double xEqua ( struct field * f1 , enum division targ ,INT_TYPE tspin,struct field * f2 , enum division orig,INT_TYPE ospin );
double xTwoBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
double xThreeBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
double xFourBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
double xOneBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic);
void printVectorAllocations(struct field *f1);
struct basisElement grabBasis (struct field * f1, INT_TYPE component, INT_TYPE particle, INT_TYPE elementIndex);
#endif /* coreUtil_h */
