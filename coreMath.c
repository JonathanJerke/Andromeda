/**
 *  coreMath.c
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

#include "coreMath.h"
///have not bothered to put this in a better spot because its double loading
#undef __cplusplus
#ifndef FADDEEVA
#include "../Faddeeva/Faddeeva.cc"
#define FADDEVA
#endif

/**
 *Kronecker delta
 */
double delta ( inta n ){
    if ( n == 0 )
        return 1.;
    else
        return 0.;
}

/**
 *even/odd sign
 */
inta sign( inta n ){
    if ( llabs(n) % 2 == 0 )
        return 1;
    else
        return -1;
};

/**
 * integer max
 */
inta imax( inta x1, inta x2 ){
    if ( x1 > x2 )
        return x1;
    else {
        return x2;
    }
}
/**
 * integer min
*/
inta imin( inta x1, inta x2 ){
    if ( x1 < x2 )
        return x1;
    else {
        return x2;
    }
}

/**
 *float max
*/
double max( double x1, double x2 ){
    if ( x1 > x2 )
        return x1;
    else {
        return x2;
    }
}

/**
 *float max
*/
double min( double x1, double x2 ){
    if ( x1 < x2 )
        return x1;
    else {
        return x2;
    }
}

/**
 *b^n
*/

double Power ( double b, inta n ){
    inta i;
    double va = 1.;
    if ( n > 0 ){
        for ( i = 0; i < n ; i++)
            va *= b;
    } else if ( n < 0 ){
        for ( i = 0; i < -n ; i++)
            va /= b;

    }
    return va;
}

/**
 * Center of the magic, Faddeeva routines
 */
DCOMPLEX expErf ( DCOMPLEX z ){
    DCOMPLEX w;
    w = Faddeeva_w(z*I,1e-15);
    return cexp(-(cimag(z)*cimag(z))) - cexp(-(creal(z)*creal(z))-2.*I*creal(z)*cimag(z)) *w;
};

/**
 *Center of the magic, the implementation of the Sinc integrals
 */
double momentumIntegralInTrain ( double beta, double kl , double d,   genusType hidden,   bodyType body ){
    switch(hidden){
        case eikonDiagonal:
                switch ( body ){
                        case one :
                        //the nominal exp(-2pi*sinc-n) is erased
                        return creal(sqrt(pi)*(
                                               2.*(pi+I*beta*beta*kl*d*d)*expErf(pi/beta/d+I*kl*beta*d)
                        /*+(pi-I*beta*beta*kl*d)*expErf(pi/beta/d-I*kl*beta*d)*/
                                               -I* 2.* d*d*kl *beta*beta * expErf(I * beta * kl*d)
                                               )
                                     + 2.*exp(-(pi/beta/d)*(pi/beta/d))*cos(2.*pi*kl)*beta - 2*beta)/sqrt(4.*pi*pi*pi);
                    case two:
                        return creal(
            sqrt(pi)*(2.*(2.*pi*pi + beta*beta + I * 4. *pi *kl*beta*beta -2.*kl*kl*beta*beta*beta*beta)*expErf(pi/beta+I * kl * beta)
             + 8.*pi*beta*beta*kl*expErf(I*kl*beta)/I)
                        -8.*pi*beta
            +
                                     (4. *exp(-pi*pi/beta/beta)*(pi*beta*cos(2.*pi*kl)+beta*beta*beta*kl*sin(2.*pi*kl)))
                                        
                                     )/sqrt(16.*pi*pi*pi*pi*pi);
                default:
                    break;

                
                }
        case eikonSemiDiagonal:
            switch ( body ){
                case one :
                    printf("error\n");
                    return 0;
                case two:
                    return (-4.*sqrt(pi)*creal((-I *pi + beta*beta * kl )*expErf(pi/beta + I * kl*beta) + pi* I*expErf( I * kl * beta)) + 0*4. *exp(-pi*pi/beta/beta) * sin(2.*pi*kl) * beta)/sqrt(64.*pi*pi*pi*pi*pi);
            default:
                break;

            }
        case eikonOffDiagonal:
            switch (body){
                case one:
                    return creal(-2*I*expErf(pi/beta/d-I*beta*kl*d) /* -I*expErf(pi/beta/d+I*beta*kl*d) */-2 * I * expErf(I * kl * beta*d))/(4*pi);
                case two:
                    return 2.*creal(expErf(pi/beta + I * kl * beta ))/(8.*pi*pi);
                    default:
                        break;

            }
        default:
            break;

    }
    return 0.;
}
    
