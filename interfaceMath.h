/*
 *  interfaceMath.h
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

#ifndef interfaceMath_h
#define interfaceMath_h
#include "constants.h"
#include "coreUtil.h"
INT_TYPE tdsyev( INT_TYPE rank, struct field * f1, char job , INT_TYPE n, double * ar, INT_TYPE ns , double * w );
INT_TYPE tdsygv( INT_TYPE rank, struct field * f1, char job , INT_TYPE n, double * sr, double * ar, INT_TYPE ns , double * w );
INT_TYPE tdgeqr( INT_TYPE rank, struct field * f1,INT_TYPE len, INT_TYPE n, double * ar, INT_TYPE ns ,double *w);
INT_TYPE tzheev( INT_TYPE rank, struct field * f1, char job , INT_TYPE n,DCOMPLEX * ar, INT_TYPE ns , double * w );
INT_TYPE tzhegv( INT_TYPE rank, struct field * f1, char job , INT_TYPE n,DCOMPLEX * sr, DCOMPLEX * ar, INT_TYPE ns , double * w );
void transpose(INT_TYPE N, INT_TYPE M, Stream_Type * orig, Stream_Type* targ);
double Sinc( double d , double x);
double periodicSinc ( double d , double x, INT_TYPE N1 );
double SS( double d1 , double x , double d2, double y ) ;
double periodicSS( double d1 , double x ,INT_TYPE N1, double d2, double y,INT_TYPE N2 )    ;
INT_TYPE tdpotrf ( INT_TYPE L1, double * array ) ;
INT_TYPE tdpotrs ( INT_TYPE L1, INT_TYPE M2, double * array, double * arrayo );
double tdpocon (INT_TYPE rank,struct field * f1,  INT_TYPE L1 , double * Matrix );
#endif /* interfaceMath_h */
