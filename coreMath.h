/*
 *  coreMath.h
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

#ifndef coreMath_h
#define coreMath_h
#include "constants.h"
#include "coreUtil.h"
#include "interfaceMath.h"
double Power ( double b, INT_TYPE n );

double sqr(double arg);
double cube ( double arg );
double delta ( INT_TYPE n );
INT_TYPE sign( INT_TYPE n );
INT_TYPE imax( INT_TYPE x1, INT_TYPE x2 );
INT_TYPE imin( INT_TYPE x1, INT_TYPE x2 );
double max( double x1, double x2 );
double min( double x1, double x2 );
INT_TYPE maxZ( struct input * f1 );
INT_TYPE sumZ( struct input * f1 );
DCOMPLEX hyperGeometric (double gamma, INT_TYPE lambda, double delta);
DCOMPLEX periodicBoostOverlapBasis (  INT_TYPE N1, double d, INT_TYPE b, double k ,double d2, INT_TYPE bb, INT_TYPE kki );
DCOMPLEX periodicBoostMomentumBasis (  INT_TYPE N1, double d, INT_TYPE b, double k ,double d2, INT_TYPE bb, INT_TYPE kki );
DCOMPLEX periodicBoostKineticBasis ( INT_TYPE N1, double d, INT_TYPE b, double k ,double d2, INT_TYPE bb, INT_TYPE kki );
DCOMPLEX periodicBoostOverlap0 ( INT_TYPE N1, double d, INT_TYPE b, double k ,double d2, INT_TYPE bb, double kk );
DCOMPLEX periodicBoostMomentum0 ( INT_TYPE N1, double d, INT_TYPE b, double k ,double d2, INT_TYPE bb, double kk );
DCOMPLEX periodicBoostKinetic0 ( INT_TYPE N1, double d, INT_TYPE b, double k ,double d2, INT_TYPE bb, double kk );
DCOMPLEX periodicBoostOverlapBasisBasis (  double p, INT_TYPE N1, double d, INT_TYPE b, INT_TYPE ki ,double d2, INT_TYPE bb, INT_TYPE kki );
DCOMPLEX periodicBoostMomentumBasisBasis (double p ,  INT_TYPE N1, double d, INT_TYPE b, INT_TYPE ki ,double d2, INT_TYPE bb, INT_TYPE kki );
DCOMPLEX periodicBoostKineticBasisBasis (  INT_TYPE N1, double d, INT_TYPE b, INT_TYPE ki ,double d2, INT_TYPE bb, INT_TYPE kki );


#endif /* coreMath_h */
