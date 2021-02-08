/**
 *  coreMath.h
 *
 *
 *  Copyright 2021 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University,
 *  the Robert A. Welch Foundation, and the Army Research Office.
 *
 
 *   *   This file is part of Andromeda.
 
 *   *   Andromeda is free software: you can redistribute it and/or modify
 *   *   it under the terms of the GNU General Public License as published by
 *   *   the Free Software Foundation, either version 3 of the License.
 
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

floata SymmetrizedGaussianInSinc( floata K, inta n , inta m , floata X );
floata GaussianInSinc( floata K, inta n, floata alpha, floata y, floata X );
double Power ( double b, inta n );
extern double complex Faddeeva_erfcx(double complex z, double relerr);
DCOMPLEX expErf ( DCOMPLEX z );
double delta ( inta n );
inta sign( inta n );
inta imax( inta x1, inta x2 );
inta imin( inta x1, inta x2 );
double max( double x1, double x2 );
double min( double x1, double x2 );
double momentumIntegralInTrain ( double beta, double kl , double d,  genusType hidden,   bodyType body );
double momentumSumInPeriodicTrain ( double k , double l ,inta N1, inta Q );
#endif /* coreMath_h */
