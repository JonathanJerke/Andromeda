/*
 *  saUtil.h
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


#ifndef saUtil_h
#define saUtil_h

#include "constants.h"
#include "coreUtil.h"
#include "mAls.h"
#define sr7 0.3779644730092272272
#define sr6 0.4082482904638630164
#define sr5 0.4472135954999579393
#define hf  0.5000000000000000000
#define sr3 0.5773502691896257645
#define sr2 0.7071067811865475244

enum grp {
    nullGroup,
    A1,
    A2,
    EE,
    T1,
    T2
};

void testSAAgain ( struct sinc_label f1 , enum division vector );
INT_TYPE tPerms(enum bodyType bd);
void tTestSA (enum bodyType bd, INT_TYPE n);
INT_TYPE tSize(enum bodyType bd);
INT_TYPE tPaths(enum bodyType bd , INT_TYPE irrep );
INT_TYPE nEqua(enum bodyType bd, INT_TYPE *a );
void tIntr ( enum bodyType bd , INT_TYPE eq , double * a);
INT_TYPE tInnerTest( struct sinc_label f1, enum division A ,enum division B);
double deg(struct sinc_label , INT_TYPE cl );
INT_TYPE tSA (enum bodyType bd, INT_TYPE X, INT_TYPE Y, INT_TYPE Z, INT_TYPE T );
double tGetProjection( enum bodyType bd, INT_TYPE type , INT_TYPE perm );
double tGetVector(enum bodyType bd , INT_TYPE type , INT_TYPE perm );
INT_TYPE tClassifyComponents( struct sinc_label  , double * up, double * entropy );
INT_TYPE tClassify(INT_TYPE rank, struct sinc_label  , enum division label);
INT_TYPE tBuildIrr ( INT_TYPE rank, struct sinc_label , INT_TYPE meta , enum division origin, INT_TYPE ospin, enum division targ , INT_TYPE tspin);
INT_TYPE matrixAction ( enum bodyType bd,INT_TYPE act, enum block bk, INT_TYPE direction);
INT_TYPE tPermuteOne(INT_TYPE rank, struct sinc_label , INT_TYPE dim, INT_TYPE leftChar , enum division left, INT_TYPE l, INT_TYPE lspin, enum division equals,INT_TYPE e, INT_TYPE espin);
INT_TYPE tCat3(enum bodyType bd ,  INT_TYPE irrep,INT_TYPE cat, INT_TYPE space);
INT_TYPE tBuild3IrrOne ( INT_TYPE rank, struct sinc_label  f1,INT_TYPE space, INT_TYPE meta , enum division origin, INT_TYPE ospin, enum division targ , INT_TYPE tspin);
INT_TYPE tBuild3Irr ( INT_TYPE rank, struct sinc_label  f1, INT_TYPE meta , enum division origin, INT_TYPE ospin, enum division targ , INT_TYPE tspin);
INT_TYPE tPermute(INT_TYPE rank, struct sinc_label , INT_TYPE leftChar , enum division left, INT_TYPE lspin, enum division equals, INT_TYPE espin);
INT_TYPE tAllCompPermMultiplyMP( INT_TYPE rank, struct sinc_label  f1 , enum division left ,INT_TYPE lspin, enum division right ,INT_TYPE rspin, double * sequ);
enum bodyType commandSA(enum bodyType bd, INT_TYPE act, enum block cl, enum block bl,INT_TYPE perm[], INT_TYPE op[]);
INT_TYPE tAddUpComponents( INT_TYPE rank, struct sinc_label  f1 , enum division left , enum division right ,  double *up);
INT_TYPE tTabulateInnerProjection( INT_TYPE rank, struct sinc_label  f1 , enum division vec, double *up);
double tTestInner ( enum bodyType bd, INT_TYPE i , INT_TYPE j );
INT_TYPE tTest ( enum bodyType bd );
INT_TYPE tSizeUp(INT_TYPE rank, struct sinc_label  f1 , INT_TYPE type, enum division label);
INT_TYPE testSA ( struct sinc_label f1 , enum division vector );
INT_TYPE irreps ( enum bodyType bd, INT_TYPE type );
INT_TYPE mapir(enum bodyType bd , INT_TYPE class );
#endif /* saUtil_h */