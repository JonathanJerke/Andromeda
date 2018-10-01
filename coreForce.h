/*
 *  coreForce.h
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

#ifndef coreForce_h
#define coreForce_h

#include "constants.h"
#include "coreUtil.h"
#include "coreMath.h"
#include "mAls.h"


DCOMPLEX ei ( double arg );
double test ( double p , struct general_index * pa );
DCOMPLEX FSS ( double p , struct general_index * pa );
DCOMPLEX FDD ( double p , struct general_index * pa );
double periodicSSGSSa( double alpha ,struct general_2index * pa);
double periodicSGSa( double alpha ,struct general_2index * pa);
double periodicSGS ( double a, double x0, INT_TYPE power, double d , INT_TYPE n, INT_TYPE m, INT_TYPE N1  );
double periodicSSGSS ( double a, double d , INT_TYPE n  , INT_TYPE m ,INT_TYPE n2, INT_TYPE m2 , INT_TYPE N1);

double eeFunc ( double k , double alpha , struct general_2index * pa );
double potFunc ( double k, double alpha, struct general_2index * pa );

void cubicFunc(void * arg,size_t n,const double * x,double * y);
double elementCal (double a, double b,struct general_2index * aAf );
void myBuildOne (struct field * f1,double beta ,   enum division transfer, INT_TYPE periodic);
void myBuildTwo (struct field * f1,double beta ,   enum division transfer, INT_TYPE periodic);
void getDescription ( struct function_label *fn ,double scalar,FILE * outString);

double findMyInterval1 ( struct field * f1  , double gamma , struct interaction_label lab,double lvl, INT_TYPE periodic);
double findMyInterval2 ( struct field * f1  , double gamma , struct interaction_label lab,double lvl,INT_TYPE periodic);
void mySeparateExactOne (struct field * f1, double scalar, enum division basis);
void mySeparateExactTwo (struct field * f1, INT_TYPE periodic, double scalar,  enum division basis);
INT_TYPE separateExternal( struct calculation * c1,INT_TYPE periodic, INT_TYPE atom,double scalar, INT_TYPE dim, enum division basis );
INT_TYPE separateBackground( struct calculation * c1,INT_TYPE periodic, INT_TYPE Ns,INT_TYPE background, INT_TYPE dim, enum division basis );
INT_TYPE separateKinetic( struct field * f1, INT_TYPE periodic,double gamma );
INT_TYPE separateBoost( struct field * f1, enum division in,INT_TYPE dim ,double vectorMomentum );
INT_TYPE buildElectronProtonInteraction ( struct field * f1, enum division mat);
INT_TYPE separateHarmonicExternal( struct calculation * c1,INT_TYPE atom, double scalar, INT_TYPE dim, enum division basis );
void separateX ( struct field * f1, double vectorDipole );
double tTestTwoBody( struct field * f1, enum division mat,INT_TYPE periodic, INT_TYPE * p);
double tRMSDevRandom( struct field * f1, enum division mat, INT_TYPE periodic ,INT_TYPE Nc);
#endif /* coreForce_h */
