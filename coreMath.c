/*
 *  coreMath.c
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

#include "coreMath.h"


double sqr(double arg){
    return arg*arg;
}

double cube ( double arg ){
    return arg*arg*arg;
}

double delta ( INT_TYPE n ){
    if ( n == 0 )
        return 1.;
    else
        return 0.;
}

INT_TYPE sign( INT_TYPE n ){
    if ( llabs(n) % 2 == 0 )
        return 1;
    else
        return -1;
};

INT_TYPE imax( INT_TYPE x1, INT_TYPE x2 ){
    if ( x1 > x2 )
        return x1;
    else {
        return x2;
    }
}
INT_TYPE imin( INT_TYPE x1, INT_TYPE x2 ){
    if ( x1 < x2 )
        return x1;
    else {
        return x2;
    }
}


double max( double x1, double x2 ){
    if ( x1 > x2 )
        return x1;
    else {
        return x2;
    }
}

double min( double x1, double x2 ){
    if ( x1 < x2 )
        return x1;
    else {
        return x2;
    }
}



double traceOne( struct field * f1 , enum division label , INT_TYPE spin ){
    Stream_Type * base;
    double sum,sum2,product;
    INT_TYPE l,i,space;
    INT_TYPE *N1 = vectorLen(f1, label);
    if ( species(f1, label) != matrix ){
        printf("oops trace %d\n",label);
        return 0;
    }
    sum2 = 0.;
    //printf("trace %d (%d), %d-body \n", label, spin, bodies(f1,label)  );
    for ( l = 0 ; l < CanonicalRank(f1,label,spin); l++)
    {
        product = 1.;
        for ( space = 0; space < SPACE ; space++){
            sum = 0.;
            base = streams(f1, label, spin, space );
            
            for ( i = 0; i < N1[space] ; i++ ){
                sum += base[ l*N1[space]*N1[space] + i*N1[space]+i];
                
            }
            product *= sum;
         //   printf("%lld %lld -- %f\n", l, space+1,sum);
        }
        sum2 += product;//*coeff(f1, label ,spin, l);
       // printf("++  %f \n", product);
    }
    
    
    return sum2;
}


INT_TYPE maxZ( struct field * f1 ){
    INT_TYPE maxz = 0;
    INT_TYPE a,ai;
    for ( a= 1 ;a <= f1->Na ; a++)
        if ( maxz < f1->atoms[a].label.Z){
            maxz=f1->atoms[a].label.Z;
            ai = a;
        }
    return maxz;
}

INT_TYPE sumZ( struct field * f1 ){
    INT_TYPE a,ai=0;
    for ( a= 1 ;a <= f1->Na ; a++)
        ai +=  f1->atoms[a].label.Z;
    return ai;
}
