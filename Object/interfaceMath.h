/**
 *  interfaceMath.h
 *
 *
 *  Copyright 2020 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University,
 *  the Robert A. Welch Foundation, and the Army Research Office.
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
inta tdsyev( inta rank,   char job , inta n, double * ar, inta ns , double * w );
inta tdsygv( inta rank,   sinc_label f1, char job , inta n, double * sr, double * ar, inta ns , double * w );
inta tdgeqr( inta rank,   sinc_label f1,inta len, inta n, double * ar, inta ns ,double *w, double *xr , inta xs );
inta tzheev( inta rank,   sinc_label f1, char job , inta n,DCOMPLEX * ar, inta ns , double * w );
inta tzhegv( inta rank,   sinc_label f1, char job , inta n,DCOMPLEX * sr, DCOMPLEX * ar, inta ns , double * w );
void transpose(inta N, inta M, floata * orig, floata* targ);
double Sinc( double d , double x);
DCOMPLEX periodicSinc ( double d , double x, double momentum, inta N1 );
double SS( double d1 , double x , double d2, double y ) ;
inta tdpotrf ( inta L1, double * array ,inta LS1) ;
inta tdpotrs ( inta L1, inta M2, double * array,inta LS1, double * arrayo,inta inc);
double tdpocon (inta rank,  sinc_label f1,  inta L1 , double * Matrix ,inta stride);
inta tdgels ( inta rank,  sinc_label f1 , inta L1, inta M2, double * array, double * arrayo ,inta inc);
inta tdgesvd ( inta rank,   sinc_label f1 ,  inta M1, inta M2, floata * ge, floata * m1, floata* m2 );
inta tInverse(   sinc_label f1, inta n, double * ar);
inta tzInverse(   sinc_label f1, inta n, DCOMPLEX * ar);
inta tLowdin( inta n , double *ar, double *lowdinVec, double * lowdinMatrix );
#endif /* interfaceMath_h */
