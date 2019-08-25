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
#include "ioPrint.h"

double HeavisideTheta ( double x );
DCOMPLEX Complex ( double r, double i );
DCOMPLEX aaGGCGG( double beta, struct general_2index * pa);
//DCOMPLEX aabGdnGdm( INT_TYPE n,INT_TYPE m, double boost, struct general_index * pa);
DCOMPLEX aaGdnGdm( INT_TYPE n, INT_TYPE m, struct general_index * pa);
DCOMPLEX aaGGCD( double invbeta,  double position,struct general_index * pa);
double aaGetGamma (  double b1,INT_TYPE l1, double o1,double b2,INT_TYPE l2,double o2);

DCOMPLEX BgB (double beta, struct basisElement b1, INT_TYPE action , INT_TYPE powSpace,double origin, struct basisElement b2);
double GTOnorm ( struct basisElement ba );

double collective( double beta ,struct general_2index * pa);
DCOMPLEX ei ( double arg );
double No(double beta1);
double GoG( double beta1, double beta2 , double x );
//double test2 ();
DCOMPLEX FS ( double p , struct general_index * pa );
DCOMPLEX FSS ( double p , struct general_index * pa );
DCOMPLEX FSSp ( double p , struct general_index * pa );
DCOMPLEX FDD ( double p , struct general_index * pa );
DCOMPLEX FGG( double k, struct general_index * pa);
DCOMPLEX FGS( double p ,INT_TYPE dn, struct general_index * pa );
double periodicSSGSSa( double alpha ,struct general_2index * pa);
double periodicSGSa( double alpha ,struct general_2index * pa);
double periodicSGS ( double a, double x0, INT_TYPE power, double d , INT_TYPE n, INT_TYPE m, INT_TYPE N1  );
double periodicSSGSS ( double a, double d , INT_TYPE n  , INT_TYPE m ,INT_TYPE n2, INT_TYPE m2 , INT_TYPE N1);
INT_TYPE separateOverlap( struct sinc_label f1, INT_TYPE periodic,enum division akinetic,  double amass, INT_TYPE particle1 );
double eeFunc ( double k , double alpha , struct general_2index * pa );
double potFunc ( double k, double alpha, struct general_2index * pa );
double collectives (double beta , struct general_2index * pa );
void cubicFunc(void * arg,size_t n,const double * x,double * y);
double elementCal (double a, double b,struct general_2index * aAf );
void myBuildOne (struct sinc_label f1,double beta ,   enum division transfer, INT_TYPE periodic);
void myBuildTwo (struct sinc_label f1,double beta ,   enum division transfer, INT_TYPE periodic);
void getDescription ( struct function_label *fn ,double scalar,FILE * outString);
double monteCarloElementCal (double beta, struct general_2index *aAf  );
double findMyInterval1 ( struct sinc_label f1  , double gamma , struct interaction_label lab,double lvl, INT_TYPE periodic);
double findMyInterval2 ( struct sinc_label f1  , double gamma , struct interaction_label lab,double lvl,INT_TYPE periodic);
void mySeparateExactOne (struct sinc_label f1, double scalar, enum division basis);
void mySeparateExactTwo (struct sinc_label f1, struct interaction_label twoBody,enum division interactionExchange, double scalar,  enum division basis,INT_TYPE periodicOverRide, INT_TYPE particle1,INT_TYPE diagonal);
void mySeparateEwaldCoulomb1(struct sinc_label f1,INT_TYPE nVec,double *  occupy, enum division vectors,INT_TYPE part1,enum division interaction,enum division shorten, double scalar,INT_TYPE plus,double rescale, enum particleType particle);
void mySeparateExactOneByOne (struct sinc_label f1, struct interaction_label twoBody,INT_TYPE periodic, enum division interactionExchangePlus,enum division shorten ,double scalar,  INT_TYPE plus,double rescale, enum particleType particle1,enum particleType particle2);
INT_TYPE separateExternal( struct calculation * c1,struct sinc_label f1,enum division linear, INT_TYPE periodic, INT_TYPE atom,double scalar, INT_TYPE dim, enum division basis , INT_TYPE particle1 );
INT_TYPE separateBackground( struct calculation * c1,INT_TYPE periodic, INT_TYPE Ns,INT_TYPE background, INT_TYPE dim, enum division basis );
INT_TYPE separateKinetic( struct sinc_label f1, INT_TYPE periodic,enum division akinetic, double mass , INT_TYPE particle1);
INT_TYPE separateBoost( struct sinc_label f1, enum division in,INT_TYPE dim ,double vectorMomentum );
INT_TYPE buildElectronProtonInteraction ( struct sinc_label f1, enum division mat,INT_TYPE spin);
INT_TYPE separateVector( struct sinc_label f1, INT_TYPE periodic,enum division aVector,  double amass, double vs[], INT_TYPE particle1 );
void separateDerivatives( struct sinc_label f1, INT_TYPE periodic,enum division mat, INT_TYPE *x, INT_TYPE *grad,double mag,INT_TYPE particle1 );
void separateX ( struct sinc_label f1, double vectorDipole );
double tTestTwoBody( struct sinc_label f1, enum division mat,INT_TYPE periodic, INT_TYPE * p);
double tRMSDevRandom( struct sinc_label f1, enum division mat, INT_TYPE periodic ,INT_TYPE Nc);
INT_TYPE tZeroSum ( struct sinc_label f1, enum division mat,INT_TYPE spin);
DCOMPLEX BoB (struct basisElement b1, struct basisElement b2);
DCOMPLEX BdB (struct basisElement b1, struct basisElement b2);
DCOMPLEX Bd2B (struct basisElement b1, struct basisElement b2);
DCOMPLEX Bx2B (struct basisElement b1, struct basisElement b2);
DCOMPLEX periodicBoostOverlapBasis2 ( INT_TYPE N1, double d, INT_TYPE b, INT_TYPE m ,double d2, INT_TYPE bb, INT_TYPE mm );
double quadCal(struct general_2index * aAf );
INT_TYPE buildPairWisePotential(struct calculation *c1, struct sinc_label f1, enum division pair, enum particleType particle1 , INT_TYPE overline,enum spinType cmpl);
INT_TYPE buildExternalPotential(struct calculation *c1, struct sinc_label f1, enum division single, enum particleType particle1, INT_TYPE overline,enum spinType cmpl);
double testPotential (struct calculation c, struct field f , enum division wavey, enum division potential );
#endif /* coreForce_h */
