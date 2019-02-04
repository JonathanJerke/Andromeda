/*
 *  eigen.h
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


#ifndef eigen_h
#define eigen_h
#include "mAls.h"
#include "saUtil.h"
#include "interfaceMath.h"

struct sortClass {
    double * str[3];
    INT_TYPE *mmm;
    INT_TYPE i;
    INT_TYPE nG;
    INT_TYPE n1[3];
};

int sortComp (const void * elem1, const void * elem2);
int sort2Comp (const void * elem1, const void * elem2);

INT_TYPE tNBodyConstruction (struct calculation * c1, enum division build, enum division eigen);
INT_TYPE t1BodyConstruction(struct calculation * c1, enum division eigen);

INT_TYPE tSASplit ( struct field * f1, INT_TYPE type , INT_TYPE Ve ,INT_TYPE target, enum division usz, enum division vector);

INT_TYPE tFoundationLevel( struct field * f1, enum division A , double lvlm,double lvlx, INT_TYPE ops,enum division build,INT_TYPE xB, double lvl1, double lvl2, double lvl3,INT_TYPE * mmm, INT_TYPE type,double seekPower);
//void tDFTChallenge(struct field * f1, INT_TYPE index);
INT_TYPE tOCSB (struct calculation * c1 , enum division usz);
INT_TYPE tSquareVectors(struct calculation * c1, INT_TYPE EV2, enum division usz,enum division usr );
INT_TYPE tGreatDivideIteration ( struct field * f1, enum division A , INT_TYPE I1, INT_TYPE I2, enum division usz, INT_TYPE foundation,INT_TYPE nMult, INT_TYPE shift);
INT_TYPE tMinorDivideIteration ( struct field * f1, enum division A , INT_TYPE I1, INT_TYPE I2, enum division usz, INT_TYPE foundation,INT_TYPE nMult, double shift);

INT_TYPE tEigenCycle (struct field * f1, enum division A ,char permutation,  INT_TYPE Ne, enum division usz, INT_TYPE quantumBasisSize ,INT_TYPE iterations,INT_TYPE foundation, INT_TYPE type,INT_TYPE flag,  enum division outputSpace, enum division outputValues);
INT_TYPE tSelect(struct field * f1, INT_TYPE Ve, INT_TYPE type, enum division usr, enum division usa, INT_TYPE testFlag);
INT_TYPE tCollect (struct field * f1, INT_TYPE type,enum division usz, INT_TYPE target,double seekPower);
INT_TYPE tFilter(struct field * f1, INT_TYPE Ve, INT_TYPE type, enum division usr);
INT_TYPE tSAboot(struct calculation *c1);
INT_TYPE tEdges(struct calculation *c1);
#endif /* eigen_h */
