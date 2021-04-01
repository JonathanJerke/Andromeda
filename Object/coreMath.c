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
 *Center of the magic Gaussian-Sinc novelity, the implementation of the GaussianSinc integrals
 */
DCOMPLEX gaussianSincfourierIntegralInTrain (double d, double gammax, double x, double gammay, double y, double kl, double momentum, genusType hidden ){
    if (hidden == eikonDiagonal){
            double gg = ( gammax + gammay);

            DCOMPLEX exp1 = -(
                            (1./(4.*d*d*(gg)))*
                            (d*d*(4.*gammax*gammay*x*x
                                  - 4.*gammax*x*(2.*gammay*y + I*momentum)
                                  + 4.*gammax*gammay*y*y
                                  + momentum*momentum
                                  - 4*I*gammay*y*momentum) +
                           4*pi*d*(2*I*d*kl*(gg) - 2*I*gammax*x -
                                   2*I*gammay*y + momentum) + 4*pi*pi));
            DCOMPLEX exp2 = -((1./(4.*d*d)*(gg)))*(4*pi*pi -
                                                              4*d*pi*(-2*I*x*gammax - 2*I*y*gammay +
                                                                       2*I*d*kl*(gg) + momentum) +
                                                              d*d*(4*x*x*gammax*gammay + 4*y*y*gammax*gammay -
                                                                       4*x*gammax*(2*y*gammay + I*momentum) -
                                                                   4*I*y*gammay*momentum + momentum*momentum));
            DCOMPLEX exp3 = -((4*x*x*gammax*gammay + 4*y*y*gammax*gammay -
                             4*x*gammax*(2*y*gammay + I*momentum) -
                                 4*I*
                             y*gammay*momentum + momentum*momentum)/(4*gg));
            
            DCOMPLEX expterm = 1./4*sqrt(gg/pi)/pi*d * (2* exp1 + 2*exp2 - 4 * exp3);
            
            DCOMPLEX preamble = 1./4./pi*exp(-gammax*gammay * ( x-y)/gg)*cexp(I *kl*momentum*d );
            
            DCOMPLEX experf1 = preamble * I*expErf((-((2*pi)/d) - momentum - 2*I*d*kl*(gg) +
                              2*I*(x*gammax + y*gammay))/(2*sqrt(gg)));
            DCOMPLEX experf2 = preamble * I*expErf((((2*pi)/d) - momentum - 2*I*d*kl*(gg) +
                                       2*I*(x*gammax + y*gammay))/(2*sqrt(gg)));
            DCOMPLEX experf3 = preamble * I*expErf(( - momentum - 2*I*d*kl*(gg) +
                                       2*I*(x*gammax + y*gammay))/(2*sqrt(gg)));

            DCOMPLEX term1 = -pi*experf1;
            DCOMPLEX term2 = pi*experf2;
            DCOMPLEX term3 = (I*pi + d*(2*x*gammax + 2*y*gammay - 2*d*kl*(gg) +
                                        I*momentum))*experf1;
            DCOMPLEX term4 = -2*d*d*kl*gammax*experf2;
            DCOMPLEX term5 = 2*x*d*gammax*experf2;
            DCOMPLEX term6 =- 2 * kl* d * d* gammay *experf2;
            DCOMPLEX term7 = 2* y * gammay*d * experf2;
            DCOMPLEX term8 = momentum * d * I * experf2;
            DCOMPLEX term9 = -pi * I * experf2;
            DCOMPLEX term10 =-pi * I * experf3 ;
            DCOMPLEX term11 = 4 * gammax*kl*d*d* experf3;
            DCOMPLEX term12 = -4 *x * d * gammax *experf3;
            DCOMPLEX term13 = 4 * gammay * kl*d*d*experf3;
            DCOMPLEX term14 = -4 *gammay*d*y*experf3;
            DCOMPLEX term15 = - 2*momentum*I*d*experf3;
            DCOMPLEX term16 = I * pi * experf3 ;
            return expterm+ term1+term2+term3+term4+term5+term6+term7+term8+term9+term10+term11+term12+term13+term14+term15+term16;
            
    } else if ( hidden == eikonOffDiagonal) {
        double gg = ( gammax + gammay);
        DCOMPLEX preamble = 1./4./pi*exp(-gammax*gammay * ( x-y)/gg)*cexp(I *kl*momentum*d );
        DCOMPLEX experf1 = preamble * I*expErf((-((2*pi)/d) - momentum - 2*I*d*kl*(gg) +
                          2*I*(x*gammax + y*gammay))/(2*sqrt(gg)));
        DCOMPLEX experf2 = preamble * I*expErf((((2*pi)/d) - momentum - 2*I*d*kl*(gg) +
                                   2*I*(x*gammax + y*gammay))/(2*sqrt(gg)));
        DCOMPLEX experf3 = preamble * I*expErf(( - momentum - 2*I*d*kl*(gg) +
                                   2*I*(x*gammax + y*gammay))/(2*sqrt(gg)));

        return experf1 + experf2 - 2*experf3;
            
    }
    return 0.;
}


/**
 * Spatial one-body Sinc elements
 * This can be further reduced, may be helpful for quantum computing...
 */

DCOMPLEX spatialSincfourierIntegralInTrain ( inta m1, floata lattice1, inta m2, floata lattice2, floata origin,floata momentum ){
    floata d ;
    floata b ;
    inta m,n;
    if ( lattice1 > lattice2 ){
        d = lattice1;
        b = lattice2;
        m = m1;
        n = m2;
    }else {
        d = lattice2;
        b = lattice1;
        m = m2;
        n = m1;
    }
    DCOMPLEX phase = cexp(I*momentum*origin);
    DCOMPLEX base  = sqrt(b*d)/(2. *pi);
    if ( m1 == m2 ){
        if ( fabs(momentum) < (pi/b -pi/d)){
            ///m
            return phase*base* cexp(I*m*b*momentum)*(2*pi/d) ;
        }else if ( momentum >= 0 ){
            ///r
            if ( momentum < pi/b + pi/d )
                return phase*base* cexp(I*m*b*momentum)* ( pi/d + pi/b - momentum );
            else
                return 0.;
        }else if ( momentum < 0 ){
            ///l
            if ( momentum > -pi/b - pi/d )
                return phase*base * cexp(I*m*b*momentum)* ( pi/d + pi/b + momentum );
            else
                return 0.;
        }
    }
    else {
        if ( fabs(momentum) < fabs(pi/b -pi/d)){
            ///m
            return 0.;
        }else if ( momentum >= 0 ){
            ///r
            if ( momentum < pi/d + pi/b )
                return phase*I/(-b*m+d*n)* base * (- cexp(I*m*b*(momentum-pi/d)+I*n*pi)+cexp(I*n*d*(momentum-pi/b)+I*m*pi));
            else
                return 0.;
        }else if ( momentum < 0 ){
            ///l
            if ( momentum > -pi/b - pi/d )
                return -phase*I/(-b*m+d*n)*base* (- cexp(I*m*b*(momentum+pi/d)+I*n*pi)+cexp(I*n*d*(momentum+pi/b)+I*m*pi));
            else
                return 0.;
        }

        
        
    }
    return 0.;
}



/**
 * Periodic Boundary Conditions
 * its interesting, without boundary conditions, the operator looks like Derivatives on non-diagonal terms with a vector core
 * for PBC, each body gets an operator, this forms a split operator for 2-bodies.
 *
 * MAY NEED A BLOCH K
 */
DCOMPLEX periodicSincfourierIntegralInTrain ( inta m1, inta m2,floata lattice , floata origin,inta N1, inta momentumIndex , inta order){
    DCOMPLEX su = 0.,s;
    DCOMPLEX phase = 1;//cexp(I*momentumIndex*2.*pi/N1*origin/lattice);

    inta n,m, N12 = (N1-1)/2;
    for ( n = -N12; n <= N12 ; n++){
        for ( m = -N12 ; m <= N12 ;m++){
            if ( n + m == momentumIndex ){
                s = cexp( I * 2.0/(N1) * pi* ( n * m1 + m * m2 ) );
                if ( order > 0 )
                    s *= cpow( I * 2.* pi *m/lattice/N1,order );
                su += s;
                ///all variables are in ratio against total length, so no lattice required.
            }
        }
    }
    return su*phase/N1;
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
