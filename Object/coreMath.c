/**
 *  coreMath.c
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

/**
 * Periodic Boundary Conditions
 * its interesting, without boundary conditions, the operator looks like Derivatives on non-diagonal terms with a vector core
 * for PBC, each body gets an operator, this forms a split operator for 2-bodies.
 */
double momentumSumInPeriodicTrain ( double k, double l , inta N1, inta Q ){
    double su = 0.;
    inta n,m, N12 = (N1-1)/2;
    for ( n = -N12; n <= N12 ; n++){
        for ( m = -N12 ; m <= N12 ;m++){
            if ( n + m == Q ){
                su += cos( 2.0/N1 * pi * ( n * k + m * l ) );
            }
        }
    }
    return 2.0*pi*su/(N1);
}
  
/**
 *Bill's Magic
 * n=1 -> 1/2, n = 3 -> 3/2 ...
 */
floata SymmetrizedGaussianInSinc( floata K, inta n , inta m , floata X ){
    double spi = sqrt(pi);
    if ( fabs(X-m*spi)> 5*spi )
        return 0.;
    return creal(-sqrt(spi)*I/2. * cexp(-(1./4.)* I *n *(pi + 2. *spi * (m *spi + X)))*
                 (
                  cexp(0.5*I*(1.+2.*m)* pi * n ) *(
                                         +expErf((- K +  I *m *spi + 0.5*n *spi - I* X))
                                         -expErf((  K +  I *m *spi + 0.5*n *spi - I* X))
                                         )
                  +cexp(I * n * spi * X )*(
                                           -expErf(( - K -  I *m *spi + 0.5*n *spi  + I* X))
                                           +expErf((   K -  I *m *spi + 0.5*n *spi  + I* X))
                                           )
                  )
                 )/sqrt(2.*K/sqrt(2.)) ;
}



/*
 *Old fashion made new
 * GTOs
 * GaussianInSinc = < (x-y)^n Exp[-alpha (x-y)^2] *normalization | X > < X | Sinc >
 *
 */
floata GaussianInSinc( floata K, inta n, floata alpha, floata y, floata X ){
    floata sspi = sqrt(sqrt(pi/2./alpha));
    floata spi = sqrt(pi);
    floata sa   = sqrt(alpha);
    X -= y;
    floata erfBase = sa * ( expErf(( K - 2. * I * X * alpha )/(2.*sa))+expErf(( K + 2. * I * X * alpha )/(2.*sa)) );
    if ( n == 0 )
        ///S
        return (erfBase)*sspi/sqrt(2.*K);
    else if ( n == 1 ){
        ///P
        floata func = (
                                -sin(K*X)
                                );
        return (4.*exp(-K*K/alpha/4.)*func/spi + erfBase * 2. * X * sa)*sspi /sqrt(2.*K);
    }
    else if ( n == 2 ){
        ///D
        floata func = (
                                      cos(K*X) * K
                               - 2. * sin(K*X) * X * alpha
                               )/sqrt(3.*alpha);
        return (4.*exp(-K*K/alpha/4.)*func/spi + erfBase * 4. * X * X * alpha / sqrt(3.))*sspi/sqrt(2.*K);
    }
    else if ( n == 3 ){
        ///F
        floata func = (
                                 sin(K*X) * K * K
                            + 2.*cos(K*X) * K * X *alpha
                            - 2.*sin(K*X) * alpha* (1. + 2. * X*X*alpha)
                               )/sqrt(15.*alpha*alpha);
        return (4.*exp(-K*K/alpha/4.)*func/spi + erfBase * 8. * X * X * X * alpha * sa/ sqrt(15.))*sspi/sqrt(2.*K);
    }
    else if ( n == 4 ){
        ///G
        floata func = (
                               -    cos(K*X) * K * K * K
                               + 2.*sin(K*X) * K * K * X * alpha
                               + 2.*cos(K*X) * K * alpha * (3.+2.*X*X*alpha)
                               - 4.*sin(K*X) * X * alpha * alpha *(1.+2.*X*X*alpha)
                          )/sqrt(105.*alpha*alpha*alpha);
        return (4.*exp(-K*K/alpha/4.)*func/spi + erfBase * 16. * X * X * X * X * alpha * alpha / sqrt(105.))*sspi/sqrt(2.*K);
    }
    else if ( n == 5 ){
        ///H
        floata func = (
                                -    sin(K*X) * K * K * K * K
                                - 2.*cos(K*X) * K * K * K * X * alpha
                                + 4.*sin(K*X) * K * K * alpha * (3.+X*X*alpha)
                                + 4.*cos(K*X) * K * X * alpha * alpha*(3.+2.*X*X*alpha)
                                - 4.*sin(K*X) * alpha * alpha * (3.+2*X*X*alpha+4.*alpha*alpha*X*X*X*X)
                               )/3./sqrt(105.*alpha*alpha*alpha*alpha);
        return (4.*exp(-K*K/alpha/4.)*func/spi + erfBase * 32. * X * X * X * X * X * alpha*alpha * sa /3./ sqrt(105.))*sspi /sqrt(2.*K);
    }
    else if ( n == 6 ){
        ///I
        floata func = (
                                    cos(K*X) * K * K * K * K * K
                               - 2.*sin(K*X) * K * K * K * K * X * alpha
                               - 4.*cos(K*X) * K * K * K * alpha * alpha *(5.+X*X*alpha)
                               + 8.*sin(K*X) * K * K * X * alpha * alpha *(3.+X*X*alpha)
                               + 4.*cos(K*X) * K * alpha * alpha * alpha *(15.+6.*X*X*alpha+4.*X*X*alpha*alpha)
                               - 8.*sin(K*X) * X * alpha * alpha * alpha * alpha*(3.+2*X*X*alpha+4.*alpha*alpha*X*X*X*X)
                          )/3./sqrt(1155.*alpha*alpha*alpha*alpha*alpha);
        return (4.*exp(-K*K/alpha/4.)*func/spi + erfBase * 64. * X * X * X * X * X * X * alpha * alpha * alpha /3./ sqrt(1155.))*sspi/sqrt(2.*K);
    }

    return 0;
}
