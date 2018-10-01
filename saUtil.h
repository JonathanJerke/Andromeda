/*
 *  saUtil.h
 *
 *
 *  Copyright 2018 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University
 *  and the Robert A. Welch Foundation.
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

INT_TYPE tInnerTest( struct field * f1, enum division A ,enum division B);
double deg(struct field *f1, INT_TYPE cl );
double get(enum body bd , INT_TYPE type , INT_TYPE i );;

INT_TYPE tClassifyComponents( struct field * f1 , double * up, double * entropy);
INT_TYPE tClassify(INT_TYPE rank, struct field * f1 , enum division label);
INT_TYPE tBuildIrr ( INT_TYPE rank, struct field * f1, char type/*1..6 (3b)*/ , enum division origin, INT_TYPE ospin, enum division targ , INT_TYPE tspin);
INT_TYPE tPermute(INT_TYPE rank, struct field * f1, char leftChar , enum division left, INT_TYPE lspin, enum division equals, INT_TYPE espin);
double tPermMultiply( INT_TYPE rank, struct field * f1 , char type , enum division left ,INT_TYPE lspin, enum division right ,INT_TYPE rspin);
INT_TYPE tAllCompPermMultiplyMP( INT_TYPE rank, struct field * f1 , enum division left ,INT_TYPE lspin, enum division right ,INT_TYPE rspin, double * sequ);
INT_TYPE tAddUpComponents( INT_TYPE rank, struct field * f1 , enum division left , enum division right ,  double *up);
INT_TYPE xConstructFoundation (struct calculation * calc , enum division usr, INT_TYPE UR, struct calculation * calc2, enum division usz, INT_TYPE UZ ,INT_TYPE mx);

#endif /* saUtil_h */
