/*
 *  coreMath.h
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

#ifndef coreMath_h
#define coreMath_h
#include "constants.h"
#include "coreUtil.h"
#include "interfaceMath.h"
#include "Faddeeva.h"
double Power ( double b, INT_TYPE n );

extern double complex Faddeeva_erfcx(double complex z, double relerr);
DCOMPLEX expErf ( DCOMPLEX z );
double momentumIntegralInTrain ( double beta, double kl , double d,enum genus hidden, enum bodyType body );
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
DCOMPLEX periodicBoost0 (INT_TYPE N, double P, INT_TYPE bb, double kk ,double dd, INT_TYPE LL,INT_TYPE b, double k ,double d,INT_TYPE L);
DCOMPLEX periodicBoostBasisBasis( INT_TYPE N , double P, INT_TYPE bb, INT_TYPE kki ,double dd, INT_TYPE LL,  INT_TYPE b, INT_TYPE ki,double d,INT_TYPE L );
#endif /* coreMath_h */
