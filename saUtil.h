/*
 *  saUtil.h
 *
 *
 *  Copyright 2018 Jonathan Jerke and Bill Poirier.
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


#ifndef saUtil_h
#define saUtil_h

#include "constants.h"
#include "coreUtil.h"
#include "mAls.h"
#define sr7 0.37796447301
#define sr6 0.40824829046
#define sr5 0.44721359550
#define hf  0.50000000000
#define sr3 0.57735026918
#define sr2 0.70710678118


enum grp {
    A1,
    A2,
    EE,
    T1,
    T2
};
INT_TYPE tPerms(enum body bd);
void tTestSA (enum body bd, INT_TYPE n);
INT_TYPE tSize(enum body bd);
INT_TYPE tPaths(enum body bd , INT_TYPE type );
void gm ( enum body bd, double *b, double *m , double *a);
INT_TYPE nEqua(enum body bd, INT_TYPE *a );
void tIntr ( enum body bd , INT_TYPE eq , double * a);
INT_TYPE tInnerTest( struct field * f1, enum division A ,enum division B);
double deg(struct field *f1, INT_TYPE cl );
double get(enum body bd , INT_TYPE type , INT_TYPE i );;
double getter(enum body bd , INT_TYPE type , INT_TYPE i );;

INT_TYPE tSA (enum body bd, INT_TYPE X, INT_TYPE Y, INT_TYPE Z, INT_TYPE T );
INT_TYPE tClassifyComponents( struct field * f1 , double * up, double * entropy);
INT_TYPE tClassify(INT_TYPE rank, struct field * f1 , enum division label);
INT_TYPE tBuildIrr ( INT_TYPE rank, struct field * f1, char type/*1..6 (3b)*/ , enum division origin, INT_TYPE ospin, enum division targ , INT_TYPE tspin);
INT_TYPE t1Permute( INT_TYPE rank, struct field * f1, char leftChar, enum division left, INT_TYPE lspin , enum division equals, INT_TYPE espin, INT_TYPE space );
INT_TYPE tPermute(INT_TYPE rank, struct field * f1, char leftChar , enum division left, INT_TYPE lspin, enum division equals, INT_TYPE espin);
double tPermMultiply( INT_TYPE rank, struct field * f1 , char type , enum division left ,INT_TYPE lspin, enum division right ,INT_TYPE rspin);
INT_TYPE tAllCompPermMultiplyMP( INT_TYPE rank, struct field * f1 , enum division left ,INT_TYPE lspin, enum division right ,INT_TYPE rspin, double * sequ);
INT_TYPE tAddUpComponents( INT_TYPE rank, struct field * f1 , enum division left , enum division right ,  double *up);
INT_TYPE tTabulateProjection( INT_TYPE rank, struct field * f1 , enum division left , enum division right ,  double *up);
INT_TYPE tTabulateComponentProjection( INT_TYPE rank, struct field * f1 , enum division left , enum division right ,  double *up);
INT_TYPE tVeto ( enum body bd , INT_TYPE type , INT_TYPE eq);
double tTestInner ( enum body bd, INT_TYPE i , INT_TYPE j );
INT_TYPE tTest ( enum body bd );
INT_TYPE tIR (enum body bd, INT_TYPE  X, INT_TYPE Y, INT_TYPE Z,INT_TYPE T );
double tGetGet ( enum body bd, INT_TYPE i , INT_TYPE j );
INT_TYPE tSizeUp(INT_TYPE rank, struct field * f1 , INT_TYPE type, enum division label);
INT_TYPE tIterator(struct field *f1, INT_TYPE * cSA ,INT_TYPE type, INT_TYPE index);

#endif /* saUtil_h */
