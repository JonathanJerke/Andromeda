/*
 *  coreForce.h
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

#ifndef coreForce_h
#define coreForce_h

#include "constants.h"
#include "coreMath.h"
#include "coreUtil.h"
#include "mAls.h"

double c10c00 ( double a, double bd, double xd, double b1, double x1, double b2, double x2, double b3 ,double x3 );
double c11c00 ( double a, double bd1, double xd1, double bd2, double xd2, double b2, double x2, double b3 ,double x3 );
double c11c00 ( double a, double bd1, double xd1, double bd2, double xd2, double b2, double x2, double b3 ,double x3 );
double c10c10 (  double a, double bd1, double xd1,  double b2, double x2, double bd2, double xd2,double b3 ,double x3 );
double c11c10 ( double a , double bd1, double  xd1,double  bd2,double  xd2, double bd3,double xd3,double  b1,double x1);
double c11c11 (double a, double bd1,double xd1,  double bd2, double xd2, double bd3, double xd3,double bd4, double xd4);
double aGGCGG(double a , struct general_2index * pa);

double collective( double beta ,struct general_2index * pa);
DCOMPLEX ei ( double arg );
double No(double beta1);
double GoG( double beta1, double beta2 , double x );
double test ( double p , struct general_index * pa );
DCOMPLEX FSS ( double p , struct general_index * pa );
DCOMPLEX FDD ( double p , struct general_index * pa );
DCOMPLEX FGG ( double p , struct general_index * pa);
double periodicSSGSSa( double alpha ,struct general_2index * pa);
double periodicSGSa( double alpha ,struct general_2index * pa);
double periodicSGS ( double a, double x0, INT_TYPE power, double d , INT_TYPE n, INT_TYPE m, INT_TYPE N1  );
double periodicSSGSS ( double a, double d , INT_TYPE n  , INT_TYPE m ,INT_TYPE n2, INT_TYPE m2 , INT_TYPE N1);

double eeFunc ( double k , double alpha , struct general_2index * pa );
double potFunc ( double k, double alpha, struct general_2index * pa );
double collectives (double beta , struct general_2index * pa );
void cubicFunc(void * arg,size_t n,const double * x,double * y);
double elementCal (double a, double b,struct general_2index * aAf );
void myBuildOne (struct field * f1,double beta ,   enum division transfer, INT_TYPE periodic);
void myBuildTwo (struct field * f1,double beta ,   enum division transfer, INT_TYPE periodic);
void getDescription ( struct function_label *fn ,double scalar,FILE * outString);
double monteCarloElementCal (double beta, struct general_2index *aAf  );
double findMyInterval1 ( struct field * f1  , double gamma , struct interaction_label lab,double lvl, INT_TYPE periodic);
double findMyInterval2 ( struct field * f1  , double gamma , struct interaction_label lab,double lvl,INT_TYPE periodic);
void mySeparateExactOne (struct field * f1, double scalar, enum division basis);
void mySeparateExactTwo (struct field * f1, INT_TYPE periodic, double scalar,  enum division basis,INT_TYPE plus, INT_TYPE particle1,INT_TYPE particle2);
INT_TYPE separateExternal( struct calculation * c1,INT_TYPE periodic, INT_TYPE atom,double scalar, INT_TYPE dim, enum division basis , INT_TYPE particle1 );
INT_TYPE separateBackground( struct calculation * c1,INT_TYPE periodic, INT_TYPE Ns,INT_TYPE background, INT_TYPE dim, enum division basis );
INT_TYPE separateKinetic( struct field * f1, INT_TYPE periodic,enum division akinetic, double mass , INT_TYPE particle1);
INT_TYPE separateBoost( struct field * f1, enum division in,INT_TYPE dim ,double vectorMomentum );
INT_TYPE buildElectronProtonInteraction ( struct field * f1, enum division mat,INT_TYPE spin);
INT_TYPE separateVector( struct field * f1, INT_TYPE periodic,enum division aVector,  double amass, double vs[], INT_TYPE particle1 );
INT_TYPE separateHarmonicExternal( struct calculation * c1,INT_TYPE periodic, double scalar, double vs[], enum division basis, INT_TYPE particle1 );
void separateX ( struct field * f1, double vectorDipole );
double tTestTwoBody( struct field * f1, enum division mat,INT_TYPE periodic, INT_TYPE * p);
double tRMSDevRandom( struct field * f1, enum division mat, INT_TYPE periodic ,INT_TYPE Nc);
INT_TYPE tZeroSum ( struct field * f1, enum division mat,INT_TYPE spin);
double BoB (struct basisElement b1, struct basisElement b2);
double BdB (struct basisElement b1, struct basisElement b2);
double Bd2B (struct basisElement b1, struct basisElement b2);
double Bx2B (struct basisElement b1, struct basisElement b2);
#endif /* coreForce_h */
