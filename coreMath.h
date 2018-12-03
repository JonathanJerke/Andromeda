/*
 *  coreMath.h
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

#ifndef coreMath_h
#define coreMath_h
#include "constants.h"
#include "coreUtil.h"
#include "interfaceMath.h"

double sqr(double arg);
double cube ( double arg );
double delta ( INT_TYPE n );
INT_TYPE sign( INT_TYPE n );
INT_TYPE imax( INT_TYPE x1, INT_TYPE x2 );
INT_TYPE imin( INT_TYPE x1, INT_TYPE x2 );
double max( double x1, double x2 );
double min( double x1, double x2 );
INT_TYPE maxZ( struct field * f1 );
INT_TYPE sumZ( struct field * f1 );
double traceOne( struct field * f1 , enum division label , INT_TYPE spin );

#endif /* coreMath_h */
