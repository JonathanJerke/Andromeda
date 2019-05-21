/*
 *  coreForce.c
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

#include "coreForce.h"


void getDescription ( struct function_label *fn ,double scalar,FILE * outString){
    if ( fn->fn == nullFunction){
        fprintf(outString,"\tnullFunction\n");
    }else if ( fn->fn == Pseudo ){
        fprintf(outString,"\tPseudo = %1.3f Erf(r/%1.3f)/r \n",scalar*fn->param[0],fn->param[2]);//, fn->param[1]);
    }else if ( fn->fn == Yukawa ){
        fprintf(outString,"\tYukawa = %1.3f exp(- r %1.3f)/r \n",scalar*fn->param[0],fn->param[2]);//, fn->param[1]);
    }else if ( fn->fn == Coulomb ){
        fprintf(outString,"\tCoulomb = %1.3f /r \n",scalar*fn->param[0]);//,fn->param[1]);
    }else if ( fn->fn == Morse ){
        fprintf(outString,"\tMorse = %1.3f (1-exp[-%1.3f *( r - %1.3f )] )^2 -1 \n",scalar*fn->param[0],fn->param[3],fn->param[2]);//,fn->param[1]);
    }else if ( fn->fn == LennardJones ){
        fprintf(outString,"\tLennardJones = %1.3f  rm = %f\n",fn->param[0],fn->param[1]);
    }
    fflush(outString);
}

DCOMPLEX Complex ( double r, double i ){
    return r + i*I;
    
}

double HeavisideTheta ( double x ){
    if ( x >= 0 )
        return 1.;
    else
        return 0.;
}

DCOMPLEX ei ( double arg ){
    return cexp(I*arg);
}

double No(double beta1){
    return 1./sqrt(sqrt(pi / 2. / beta1 ) ) ;
};//single dimensions of guassian

double Nol ( double beta1, INT_TYPE l ){
    double base =  No(beta1)/pow(beta1,0.500*l) ;
    switch (l){
        case 0:
            return base;
        case 1:
            return base;
        case 2:
            return base/sqrt(3.);
        case 3:
            return base/sqrt(15.);
        case 4:
            return base/sqrt(105.);
        case 5:
            return base/sqrt(945.);
        case 6:
            return base/sqrt(10395.);
        case 7:
            return base/sqrt(135135.);
        
    }
    exit(0);
    
}


double GoG( double beta1, double beta2 , double x ){
    double va = sqrt(pi/(beta1+beta2))*No(beta1)*No(beta2)*exp(-(beta1*beta2)/(beta1+beta2)*x*x);
    return va;
}



//DCOMPLEX FGG ( double p , struct general_index * pa){
//    double xc1,xc2;
//    DCOMPLEX fgg;
//    double b0,b1,x0,x1;
//    INT_TYPE l0,l1;
//    double x;
//
//
//    if ( pa->l0 < pa->l1 ){
//        b0 = pa->b1;
//        b1 = pa->b0;
//        l0 = pa->l1;
//        l1 = pa->l0;
//        x0 = pa->x1;
//        x1 = pa->x0;
//
//    }else {
//        b1 = pa->b1;
//        b0 = pa->b0;
//        l1 = pa->l1;
//        l0 = pa->l0;
//        x1 = pa->x1;
//        x0 = pa->x0;
//
//    }
//    x = x1-x0;
//    xc1 =GoG(b0, b1, x)*exp(-sqr(p /2.)/(b0+b1))/No(b0)/No(b1);
//
//    xc2 = ((b0*x0+ b1*x1)/(b0+b1));
//    fgg = xc1 * ei( p * xc2 );
//
//    if ( l0 == 0 && l1 == 0 ){
//        return fgg;
//    }
//    else if ( l0 == 1 && l1 == 0 ){
//        return fgg*(b0)/(b0+b1) * (  2. * b1 * x +I* p  );
//    } else if ( l0 == 1 && l1 == 1 ){
//        return -fgg*(b0*b1)/sqr(b0+b1) * (( p*p - 2.*(b1+b0)+4.*b0*b1*x*x -I*( 2. * ( b1 - b0 ) * x * p)));
//
//
//    }  else if ( l0 == 2 && l1 ==    0 ){
//        return fgg*(b0*(-2*b1*b1 - b0*(p*p - 4*b1*b1*x*x + b1*(2 - 4*I*p*x))))/(b0 + b1)/(b0 + b1);
//    }  else if ( l0 == 2 && l1 ==    1 ){
//        return (fgg/(b0 + b1)*(b0 + b1)*(b0 + b1))*(b0*b1*(-2*I*b1*b1*p + b0*((-I)*p*p*p + 2*b1*p*(I - 2*p*x) +
//                                                                              4*b1*b1*x*(3 + I*p*x)) + 2*b0*b0*(p*p*x + 2*b1*x*(3 - 2*b1*x*x) -
//                                                                                                                2*I*p*(-1 + 2*b1*x*x))));
//    }  else if ( l0 == 2 && l1 ==    2 ){
//
//        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b1*(2*b1*b1*b1*p*p+ b0*b1*(p*p*p*p + 2*b1*p*p*(-3 - 2*I*p*x) -
//                                                           4*b1*b1*(-3 - 6*I*p*x + p*p*x*x)) +
//                                       2*b0*b0*b0*(p*p + 8*b1*b1*b1*x*x*x*x + 8*I*b1*b1*x*x*(3*I + p*x) -
//                                               2*b1*(-3 + 6*I*p*x + p*p*x*x)) + 2*b0*b0*b1*(8*b1*b1*x*x*(-3 - I*p*x) +
//                                                                                            p*p*(-3 + 2*I*p*x) + 4*b1*(3 + 2*p*p*x*x))));
//
//
//    }  else if ( l0 == 3 && l1 ==    0 ){
//        return fgg*(b0*b0*(I*p + 2*b1*x)*(-6*b1*b1 -
//                                      b0*(p*p - 4*b1*b1*x*x + b1*(6 - 4*I*p*x))))/(b0 + b1)/(b0 + b1)/(b0 + b1);
//
//    }  else if ( l0 == 3 && l1 ==    1 ){
//       return  (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*(6*b1*b1*(p*p + b1*(-2 - 2*I*p*x)) +
//                                  b0*(p*p*p*p - 6*I*b1*p*p*p*x + 8*b1*b1*b1*x*x*(6 + I*p*x) -
//                                      12*b1*b1*(2 - 2*I*p*x + p*p*x*x)) - 2*b0*b0*(8*b1*b1*b1*x*x*x*x + p*p*(3 - I*p*x) +
//                                                                                   12*I*b1*b1*x*x*(2*I + p*x) - 6*b1*(-1 + 3*I*p*x + p*p*x*x))));
//    }  else if ( l0 == 3 && l1 ==    2 ){
//        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*(6*b1*b1*b1*p*(I*p*p + 2*b1*(-2*I + p*x)) +
//                                             b0*b1*(I*p*p*p*p*p + 6*b1*p*p*p*(-I + p*x) + 8*b1*b1*b1*x*(15 + 12*I*p*x - p*p*x*x) -
//                                                    12*I*b1*b1*p*(1 - 5*I*p*x + p*p*x*x)) +
//                                             2*b0*b0*b0*(I*p*p*p + 16*b1*b1*b1*b1*x*x*x*x*x + 8*b1*b1*b1*x*x*x*(-10 + 3*I*p*x) +
//                                                     2*b1*p*(9*I + 9*p*x - I*p*p*x*x) - 12*b1*b1*x*(-5 + 6*I*p*x + p*p*x*x)) +
//                                             2*b0*b0*b1*(16*b1*b1*b1*x*x*x*(-5 - I*p*x) - p*p*p*(5*I + 2*p*x) +
//                                                        24*b1*b1*x*(5 - I*p*x + p*p*x*x) + 6*I*b1*p*(4 + 3*I*p*x + 2*p*p*x*x))));
//
//    }  else if ( l0 == 3 && l1 ==    3 ){
//        return -((fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*b1*(6*b1*b1*b1*p*p*(p*p+ b1*(-6 - 2*I*p*x)) -
//                                            6*b0*b0*b1*(p*p*p*p*(2 - I*p*x) - 6*b1*p*p*(4 + p*p*x*x) +
//                                                       8*b1*b1*b1*x*x*(-15 - 10*I*p*x + p*p*x*x) +
//                                                       12*b1*b1*(5 + 5*I*p*x + 4*p*p*x*x +
//                                                                I*p*p*p*x*x*x)) +
//                                            b0*b1*(p*p*p*p*p*p + 6*b1*p*p*p*p*(-2 - I*p*x) -
//                                                   12*b1*b1*p*p*(-3 - 8*I*p*x + p*p*x*x) +
//                                                   8*b1*b1*b1*(-15 - 45*I*p*x + 18*p*p*x*x +
//                                                           I*p*p*p*x*x*x)) +
//                                            4*b0*b0*b0*b0*(16*b1*b1*b1*b1*x*x*x*x*x*x + 3*I*p*p*(3*I + p*x) +
//                                                    24*I*b1*b1*b1*x*x*x*x*(5*I + p*x) -
//                                                    12*b1*b1*x*x*(-15 + 10*I*p*x + p*p*x*x) +
//                                                    b1*(-30 + 90*I*p*x + 36*p*p*x*x - 2*I*p*p*p*x*x*x)) +
//
//                                            6*b0*b0*b0*(p*p*p*p + 16*b1*b1*b1*b1*x*x*x*x*(-5 - I*p*x) +
//                                                    24*b1*b1*b1*x*x*(10 + p*p*x*x) -
//                                                    2*b1*p*p*(-3 + 8*I*p*x + p*p*x*x) +
//                                                    12*I*b1*b1*(5*I + 5*p*x + 4*I*p*p*x*x +
//                                                                p*p*p*x*x*x)))));
//
//    }  else if ( l0 == 4 && l1 ==    0 ){
//        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*(12*b1*b1*b1*b1 -
//                                      12*b0*b1*b1*(-p*p + 4*b1*b1*x*x + b1*(-2 + 4*I*p*x)) +
//                                      b0*b0*(p*p*p*p + 16*b1*b1*b1*b1*x*x*x*x + 4*b1*p*p*(3 - 2*I*p*x) +
//                                            16*b1*b1*b1*x*x*(-3 + 2*I*p*x) -
//                                             12*b1*b1*(-1 + 4*I*p*x + 2*p*p*x*x))));
//
//    }  else if ( l0 == 4 && l1 ==    1 ){
//        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*
//                                b1*(12*I*b1*b1*b1*b1*p + 12*b0*b1*b1*(I*p*p*p + 2*b1*b1*x*(-5 - 2*I*p*x) +
//                                                              2*b1*p*(-I + 2*p*x)) +
//                                    I*b0*b0*(p*p*p*p*p + 4*b1*p*p*p*(1 - 2*I*p*x) +
//                                            16*b1*b1*b1*b1*x*x*x*(-10*I + p*x) +
//                                            16*b1*b1*b1*x*(15*I + 9*p*x + 2*I*p*p*x*x) -
//                                            12*b1*b1*p*(7 - 2*I*p*x + 2*p*p*x*x)) -
//                                    2*
//                                    b0*b0*b0*(16*b1*b1*b1*b1*x*x*x*x*x + 16*b1*b1*b1*x*x*x*(-5 + 2*I*p*x) +
//                                          p*p*p*(4*I + p*x) +
//                                          4*b1*p*(6*I + 9*p*x - 2*I*p*p*x*x) -
//                                          12*b1*b1*x*(-5 + 8*I*p*x + 2*p*p*x*x))));
//
//    }  else if ( l0 == 4 && l1 ==    2 ){
//        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*(-12*b1*b1*b1*b1*b1*p*p + 12*b0*b1*b1*b1*(-p*p*p*p + 2*b1*p*p*(3 + 2*I*p*x) +
//                                                                    2*b1*b1*(-5 - 10*I*p*x + 2*p*p*x*x)) +
//                                         2*b0*b0*b0*b0*(-p*p*p*p + 32*b1*b1*b1*b1*b1*x*x*x*x*x*x + 16*b1*b1*b1*b1*x*x*x*x*(-15 + 4*I*p*x) +
//                                                 2*b1*p*p*(-18 + 12*I*p*x + p*p*x*x) - 8*b1*b1*b1*x*x*(-45 + 40*I*p*x + 6*p*p*x*x) +
//                                                 4*b1*b1*(-15 + 60*I*p*x + 36*p*p*x*x - 4*I*p*p*p*x*x*x)) -
//                                         b0*b0*b1*(p*p*p*p*p*p + 4*b1*p*p*p*p*(-1 - 2*I*p*x) + 16*b1*b1*b1*b1*x*x*(-45 - 20*I*p*x + p*p*x*x) -
//                                                  12*b1*b1*p*p*(9 - 8*I*p*x + 2*p*p*x*x) + 8*b1*b1*b1*(45 + 42*p*p*x*x + 4*I*p*p*p*x*x*x)) +
//                                         2*b0*b0*b0*b1*(16*b1*b1*b1*b1*x*x*x*x*(-15 - 2*I*p*x) + p*p*p*p*(7 - 2*I*p*x) -
//                                                    8*b1*p*p*(3 + 6*I*p*x + 2*p*p*x*x) + 16*b1*b1*b1*x*x*(45 - 10*I*p*x + 4*p*p*x*x) +
//                                                        12*I*b1*b1*(15*I + 30*p*x + 4*I*p*p*x*x+ 4*p*p*p*x*x*x))));
//    }  else if ( l0 == 4 && l1 ==    3 ){
//
//        return -((fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*b1*(12*I*b1*b1*b1*b1*b1*p*p*p + 12*b0*b1*b1*b1*p*
//                                      (I*p*p*p*p + 2*b1*p*p*(-5*I + 2*p*x) +
//                                       b1*b1*(30*I - 30*p*x - 4*I*p*p*x*x)) +
//
//                                      4*b0*b0*b0*b0*b0*(32*b1*b1*b1*b1*b1*x*x*x*x*x*x*x + 16*b1*b1*b1*b1*x*x*x*x*x*(-21 + 4*I*p*x) -
//                                              3*p*p*p*(4*I + p*x) -
//                                              24*b1*b1*b1*x*x*x*(-35 + 20*I*p*x + 2*p*p*x*x) +
//
//                                              4*b1*b1*x*(-105 + 180*I*p*x + 60*p*p*x*x - 4*I*p*p*p*x*x*x) +
//                                              2*b1*p*(-60*I - 90*p*x + 24*I*p*p*x*x + p*p*p*x*x*x)) +
//                                      I*b0*b0*b1*(p*p*p*p*p*p*p + 4*b1*p*p*p*p*p*(-3 - 2*I*p*x) - 12*b1*b1*p*p*p*
//                                                 (5 - 14*I*p*x + 2*p*p*x*x) +
//                                                 8*b1*b1*b1*p*(75 - 90*I*p*x + 66*p*p*x*x +
//                                                           4*I*p*p*p*x*x*x) +
//                                                 16*b1*b1*b1*b1*x*(105*I - 135*p*x - 30*I*p*p*x*x + p*p*p*x*x*x)) -
//
//                                      6*b0*b0*b0*b1*(p*p*p*p*p*(3*I + p*x) +
//                                                 4*b1*p*p*p*(-10*I + 3*p*x - 2*I*p*p*x*x) +
//                                                 16*b1*b1*b1*b1*x*x*x*(-35 - 15*I*p*x + p*p*x*x) +
//                                                 8*b1*b1*b1*x*(105 + 30*I*p*x + 30*p*p*x*x +
//                                                           4*I*p*p*p*x*x*x) -
//                                                 12*b1*b1*p*(-5*I + 25*p*x - 4*I*p*p*x*x + 2*p*p*p*x*x*x)) +
//                                      6*b0*b0*b0*b0*(I*p*p*p*p*p - 32*I*b1*b1*b1*b1*b1*x*x*x*x*x*(-7*I + p*x) +
//                                              2*b1*p*p*p*(10*I + 11*p*x - I*p*p*x*x) +
//                                              16*b1*b1*b1*b1*x*x*x*(70 - 5*I*p*x + 4*p*p*x*x) -
//                                              4*b1*b1*p*(45*I + 28*I*p*p*x*x + 4*p*p*p*x*x*x) +
//
//                                                     8*I*b1*b1*b1*x*(105*I + 75*p*x + 20*I*p*p*x*x + 6*p*p*p*x*x*x)))));
//    }  else if ( l0 == 4 && l1 ==    4 ){
//
//        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*b1*(12*b1*b1*b1*b1*b1*b1*p*p*p*p - 12*b0*b1*b1*b1*b1*p*p*
//                                           (-p*p*p*p + 2*b1*p*p*(7 + 2*I*p*x) +
//                                            4*b1*b1*(-15 - 10*I*p*x + p*p*x*x)) +
//                                           4*b0*b0*b0*b1*b1*(p*p*p*p*p*p*(-5 + 2*I*p*x) + 8*b1*p*p*p*p*(15 + 2*p*p*x*x) +
//
//                                                        16*b1*b1*b1*b1*x*x*(-210 - 210*I*p*x + 45*p*p*x*x + 2*I*p*p*p*x*x*x) -
//
//                                                        24*I*b1*b1*p*p*(-30*I + 25*p*x - 10*I*p*p*x*x + 2*p*p*p*x*x*x) +
//
//                                                        16*b1*b1*b1*(105 + 210*I*p*x + 45*p*p*x*x + 50*I*p*p*p*x*x*x -
//                                                                 4*p*p*p*p*x*x*x*x)) +
//                                           b0*b0*
//                                           b1*b1*(p*p*p*p*p*p*p*p + 4*b1*p*p*p*p*p*p*(-5 - 2*I*p*x) +
//                                                 12*b1*b1*p*p*p*p*(5 + 20*I*p*x - 2*p*p*x*x) +
//
//                                                 16*b1*b1*b1*p*p*(15 - 120*I*p*x + 45*p*p*x*x + 2*I*p*p*p*x*x*x) +
//
//                                                 16*b1*b1*b1*b1*(105 + 420*I*p*x - 270*p*p*x*x - 40*I*p*p*p*x*x*x +
//                                                          p*p*p*p*x*x*x*x)) +
//                                           4*b0*b0*b0*b0*b0*b0*(3*p*p*p*p+ 64*b1*b1*b1*b1*b1*b1*x*x*x*x*x*x*x*x + 128*I*b1*b1*b1*b1*b1*x*x*x*x*x*x*(7*I + p*x) -
//                                                   12*b1*p*p*(-15 + 10*I*p*x + p*p*x*x) -
//                                                   96*b1*b1*b1*b1*x*x*x*x*(-35 + 14*I*p*x + p*p*x*x)) +
//
//                                                   16*b1*b1*b1*x*x*(-210 + 210*I*p*x + 45*p*p*x*x - 2*I*p*p*p*x*x*x) +
//
//                                                   4*b1*b1*(105 - 420*I*p*x - 270*p*p*x*x + 40*I*p*p*p*x*x*x +
//                                                           p*p*p*p*x*x*x*x)) -
//                                           12*b0*b0*b0*b0*b1*(-p*p*p*p*p*p + 32*b1*b1*b1*b1*b1*x*x*x*x*(-35 - 14*I*p*x + p*p*x*x) +
//                                                       b1*p*p*p*p*(-5 + 20*I*p*x + 2*p*p*x*x) +
//                                                       8*b1*b1*p*p*(30 - 25*I*p*x + 10*p*p*x*x -
//                                                                   2*I*p*p*p*x*x*x) +
//                                                       16*b1*b1*b1*b1*x*x*(210 + 70*I*p*x + 25*p*p*x*x + 4*I*p*p*p*x*x*x) -
//                                                       24*b1*b1*b1*(35 + 50*p*p*x*x + 2*p*p*p*p*x*x*x*x)) +
//                                           8*b0*b0*b0*b0*b0*
//                                           b1*(3*p*p*p*p*(-7 + 2*I*p*x) - 64*I*b1*b1*b1*b1*b1*x*x*x*x*x*x*(-7*I + p*x) +
//                                               32*b1*b1*b1*b1*x*x*x*x*(105 + 4*p*p*x*x) +
//                                               2*b1*p*p*(15 + 120*I*p*x + 45*p*p*x*x -
//                                                         2*I*p*p*p*x*x*x +
//                                               24*I*b1*b1*b1*x*x*(210*I + 70*p*x + 25*I*p*p*x*x + 4*p*p*p*x*x*x) -
//
//                                               8*b1*b1*(-105 + 210*I*p*x - 45*p*p*x*x + 50*I*p*p*p*x*x*x +
//                                                        4*p*p*p*p*x*x*x*x))));
//    }
//
//        return 0;
//}


//DCOMPLEX FSSprev ( double p , struct general_index * pa ){
//    INT_TYPE n = pa->n;
//    INT_TYPE m = pa->m;
//    double d = pa->d;
//    if ( fabs( p * d) >= 2 * pi )
//        return 0.;
//    if ( n == m )
//        return (2 * pi-fabs(p*d))*ei( p * d * n )/ (2*pi );
//    if ( p > 0 )
//        return (I/(-m+n)/ (2*pi)) *  ei( n * p*d ) *( -ei( (m-n)*(d*p-pi) ) +  sign(n-m)   ) ;
//    if ( p < 0 )
//        return (-I/(-m+n)/ (2*pi))  * ei( n * p*d ) *(  sign(n-m) - ei((m-n)*(d*p+pi) ) );
//    return 0.;
//};

DCOMPLEX FGS( double p ,INT_TYPE dn, struct general_index * pa ){
    
    // Integrate[[   D[Gaus( beta, (x-o) ),{o,l}] Exp[ i p x ] Sinc ( x/d- n )
    // one gaussian, one sinc, fourier transformed.
    struct general_2index g2;
    double beta,realpart,imagepart;
    g2.i[0].bra.basis = SincBasisElement;
    g2.i[0].ket.basis = nullBasisElement;
    g2.i[1].bra.basis = DiracDeltaElement;
    g2.i[1].ket.basis = nullBasisElement;
    
    g2.gaussianAccelerationFlag = 0;
    g2.momentumShift = p ;
    g2.i[0].pointer = 1;
    g2.i[1].pointer = 1;
    g2.i[0].action = 0;
    g2.i[1].action = dn;
    g2.momentumShift = 0;
    
    if ( pa->bra.basis == GaussianBasisElement && pa->ket.basis == SincBasisElement ){
        
        beta = pa->bra.length;
        g2.powSpace = pa->ket.index;
        g2.i[1].bra.origin = pa->bra.origin;
        g2.i[0].bra = pa->ket;
        g2.i[0].d = pa->ket.length;

    }
    else if ( pa->bra.basis == SincBasisElement && pa->ket.basis == GaussianBasisElement ){
        beta = pa->ket.length;
        g2.powSpace = pa->ket.index;
        g2.i[1].bra.origin = pa->ket.origin;
        g2.i[0].bra = pa->bra;
        g2.i[0].d = pa->bra.length;

    }else {
        
        
        printf("wrong mail stop\n");
        exit(0);
    }
    g2.realFlag = 1;
    realpart =  collective(sqrt(beta), &g2);
    g2.realFlag = 0;
    imagepart = collective(sqrt(beta), &g2);
    return /*Nol(beta,g2.i[1].action)*/ sqrt(pi*g2.i[0].d)*(realpart + I * imagepart);
};





DCOMPLEX FS ( double p , struct general_index * pa ){
    double n = pa->n;
    double m = pa->m;
    double d = pa->d;
    INT_TYPE pt = pa->pointer;
    
    
    
    if ( fabs( p * d) > imax(1,pt) * pi ){
        return 0.;
    }
    if ( pt < 2){
        if ( pt == 0 )
            return ei( n * p*d );
        else
            return ei( m * p*d);
    }else {
        printf("FS\n");
        exit(0);
    }
    
}

DCOMPLEX FSS( double p , struct general_index * pa ){
    double d = pa->bra.length;
    INT_TYPE n = pa->bra.index;
    double x1 = n*d + pa->bra.origin;
    
    double d2 = pa->ket.length;
    INT_TYPE m = pa->ket.index;
    double x2 = m*d2 + pa->ket.origin;
//    if ( pa->ket.origin != pa->bra.origin){
//        printf("origins\n");
//        exit(0);
//    }
    double cof = sqrt(d*d2)/(2.*pi);
    double Pi = pi;
    if ( n == m  ) {
        return cof*(cexp(I*p*x1)*HeavisideTheta(p + (1/d + 1/d2)*Pi)*((-(d*d2*p) + (d + d2)*Pi + (d*d2*p - d*Pi + d2*Pi)*HeavisideTheta(-p + (-(1/d) + 1/d2)*Pi))*
                                                                        HeavisideTheta(p + (-(1/d) + 1/d2)*Pi)*HeavisideTheta(-p + (1/d + 1/d2)*Pi) -
                                                                        (-1 + HeavisideTheta(p + (-(1/d) + 1/d2)*Pi))*((2*d*Pi + (d*d2*p - d*Pi + d2*Pi)*HeavisideTheta(-p + (-(1/d) + 1/d2)*Pi))*HeavisideTheta(p + (1/d + 1/d2)*Pi) +(d*d2*p + (d + d2)*Pi)*HeavisideTheta(-p + (-(1/d) + 1/d2)*Pi)*HeavisideTheta(-p - ((d + d2)*Pi)/(d*d2)))))/(d*d2);
        
        
    } else {
        return cof*(Complex(0,1)*HeavisideTheta(p + (1/d + 1/d2)*Pi)*((cexp(I*p*x1) - cexp(I*(((d + d2)*Pi*(x1 - x2))/(d*d2) + p*x2)) +
                                                                   (-cexp(I*p*x1) + cexp(I*((Pi*(x1 - x2))/d2 + p*x2 + (Pi*(-x1 + x2))/d)))*HeavisideTheta(-p + (-(1/d) + 1/d2)*Pi))*HeavisideTheta(p + (-(1/d) + 1/d2)*Pi)*
                                                                  HeavisideTheta(-p + (1/d + 1/d2)*Pi) + (-1 + HeavisideTheta(p + (-(1/d) + 1/d2)*Pi))*
                                                                  ((cexp(I*p*x1)*(-1 + cexp((2*I*Pi*(x1 - x2))/d2)) + (cexp(I*p*x1) - cexp(I*((Pi*(x1 - x2))/d2 + p*x2 + (Pi*(-x1 + x2))/d)))*
                                                                    HeavisideTheta(-p + (-(1/d) + 1/d2)*Pi))*HeavisideTheta(p + (1/d + 1/d2)*Pi) +
                                                                   cexp((I*Pi*(x1 - x2))/d2)*(cexp((I*(d2*p*x1 + Pi*x1 - Pi*x2))/d2) - cexp((I*(-(Pi*x1) + d*p*x2 + Pi*x2))/d))*HeavisideTheta(-p + (-(1/d) + 1/d2)*Pi)*
                                                                   HeavisideTheta(-p - ((d + d2)*Pi)/(d*d2)))))/(cexp((I*Pi*(x1 - x2))/d2)*(x1 - x2));
    }
}
DCOMPLEX FSSp ( double p , struct general_index * pa ){
    INT_TYPE n = pa->n;
    INT_TYPE m = pa->m;
    double d = pa->d;
    
    if ( fabs( p * d) > 2 * pi ){
        return 0.;
    }
    
    return ei((n+m)*(1./2.*p*d+pi)) *( delta(n-m) -fabs( p * d / 2./pi ) * Sinc(1., (p*d/2./pi) *(n-m)));

};

DCOMPLEX periodicFSSp ( double p , struct general_index * pa ){
    INT_TYPE n = pa->bra.index;
    INT_TYPE m = pa->ket.index;
    INT_TYPE N1 = pa->bra.grid,pp;
    double d = pa->bra.length;
    double sh = 2*pi/d/N1;
    double sp = 1./N1;
    DCOMPLEX sum = 0.;
    INT_TYPE N12 = (N1-1)/2;
    for ( pp = -N12; pp <= N12 ; pp++ ){
        if ( fabs(pp*sh-p) <= pi/d ){
            sum += sp*ei(m*p*d+(n-m)*d*(pp*sh));
        }
    }
    return sum;
};

//double test ( double p , struct general_index * pa ){
////    printf("%f %f\n", cabs(FSSnew(p,pa)),cabs(FSS(p,pa)));
////    printf("%f %f %f %f\n", creal(FSSnew(p,pa)),cimag(FSSnew(p,pa)),creal(FSS(p,pa)),cimag(FSS(p,pa)));
//    return cabs(FSS(p,pa)-FSSprev(p, pa));
//}




DCOMPLEX FDD ( double p , struct general_index * pa ){
    return ei( p*pa->x0 );
};



DCOMPLEX FB ( double p , struct general_index * pa ){
    
    INT_TYPE periodic , component;
    if ( pa->bra.type != nullComponent ){
        if ( pa->bra.type > 3 )
            periodic = 1;
        else
            periodic = 0;
    }
    else {
        printf("messed up definitions! %d\n",pa->bra.basis);
        exit(0);
    }
    component =  ( (pa->bra.type -1) % COMPONENT ) ;

    
    
    if ( pa->bra.basis == SincBasisElement && pa->ket.basis == SincBasisElement){
        if ( ! periodic ){
            if ( fabs(pa->bra.origin ) > 1e-6 ||  fabs( pa->bra.origin - pa->ket.origin ) > 1e-6 || fabs(pa->bra.length - pa->ket.length) > 1e-6 )
                return FSS(p,pa);
            else{
                pa->n = pa->bra.index;
                pa->m = pa->ket.index;
                pa->d = pa->bra.length;
                pa->pointer = 2;
                return FSSp(p, pa);
            }
        }else if ( periodic ){
            return periodicFSSp(p,pa);
        }
    }else if ( pa->bra.basis == DiracDeltaElement || pa->ket.basis == DiracDeltaElement){
        if ( pa->bra.basis == DiracDeltaElement && pa->ket.basis == nullBasisElement){
            pa->x0 = pa->bra.origin;
        }
        else  if ( pa->ket.basis == DiracDeltaElement && pa->bra.basis == nullBasisElement){
            pa->x0 = pa->ket.origin;
        }else {
            printf("parity correspondence\n");
            exit(0);
        }
        return FDD(p, pa);
    }else if ( pa->bra.basis == GaussianBasisElement || pa->ket.basis == GaussianBasisElement){
        return FGG(p,pa);
    }else     if ( pa->bra.basis == SincBasisElement && pa->ket.basis == nullBasisElement){
        pa->n = pa->bra.index+pa->bra.origin/pa->bra.length;
        pa->d = pa->bra.length;
        pa->pointer = 0;
        return FS(p, pa);
    }else     if ( pa->bra.basis == nullBasisElement && pa->ket.basis == SincBasisElement){
        pa->m = pa->ket.index+pa->ket.origin/pa->ket.length;
        pa->d = pa->ket.length;
        pa->pointer = 1;
        return FS(p, pa);
    }


    printf("off tracks\n");
    exit(0);
    return 0;
}


double SdS(INT_TYPE arg){
    if ( arg == 0 )
        return 0.;
    else {
        if ( (llabs(arg) % 2) == 0 )
            return 1./(double)(arg);
        else
            return -1./(double)(arg);
    }
    return 0;
}


double Sd2S(INT_TYPE arg){
    if ( arg == 0 )
        return -pi*pi/3.;
    else {
        if ( llabs(arg) % 2  == 0 )
            return -2./((double)arg*arg);
        else
            return 2./((double)arg*arg);
    }
    return 0;
}



double periodicSd2S ( double arg, INT_TYPE N ){
    
#ifdef GSL_LIB
    
    
    if ( arg == 0 ){
        return Sd2S(0)-4.0*gsl_sf_zeta(2.)/sqr(2.*N)+4.0*gsl_sf_hzeta(2.,0.500)/sqr(2.*N);
    }
    
    return -cos(pi* arg)/2./sqr(N)*(gsl_sf_hzeta(2., 1-(1.*fabs(arg))/(2.*N))+gsl_sf_hzeta(2., (1.*fabs(arg))/(2.*N))-gsl_sf_hzeta(2., 0.5+(1.*fabs(arg))/(2.*N))-gsl_sf_hzeta(2.,0.5- (1.*fabs(arg))/(2.*N)));
#else
    
    if ( N == 9 ){
        double cal9[] = {-3.249252,    1.957614,    -0.451818,    0.162463,    -0.043633,    -0.043633,    0.162463,    -0.451818,    1.957614};
      //  return cal9[abs(arg)];
        return 0;
    }else {

    
    printf("boot gsl\n");
    exit(0);
    }
#endif
    
}

double periodicSdS ( INT_TYPE arg, INT_TYPE N ){
    
    if ( ! arg )
        return 0;
    else {
        return sign(arg)* pi / N / sin( pi * arg *1./ N ) ;
    }
}


double aaGetGamma (  double b1,INT_TYPE l1, double o1,double b2,INT_TYPE l2,double o2){
    return 1./4./(b1+b2);
}

double aaGetDelta ( double b1,INT_TYPE l1, double o1,double b2,INT_TYPE l2,double o2){
    return (b1*o1+b2*o2)/(b1+b2);
}

double aaGetConst( double b1,INT_TYPE l1, double o1,double b2,INT_TYPE l2,double o2){
    return exp( - b1*b2/(b1+b2)*sqr(o1-o2))/ pow(b1+b2,1/2.)/sqrt(2.);
}

DCOMPLEX aaGetPoly( double b1,INT_TYPE l1, double o1,double b2,INT_TYPE l2,double o2, INT_TYPE lambda1 ){
    DCOMPLEX va = 0.;
    if ( l1 > l2 ){
        printf("switch");
        exit(0);
    }
    if ( l1 == 0 && l2 == 0 ){
        switch ( lambda1 ) {
            case 0:
                va =  1.;
                break;
        }
    }
    else if ( l1 == 0 && l2 == 1 ){
        switch ( lambda1 ) {
            case 0:
                va =  b1 *(o1 - o2);
                break;

            case 1:
                va =  I/2.;
                break;

        }
    }    else if ( l1 == 1 && l2 == 1 ){
        switch ( lambda1 ) {
            case 0:
                va =  (1/2.)*(b2 + b1*(1 - 2*b2*(o1 - o2)*(o1 - o2)));
                break;

            case 1:
                va =  (1/2.)*I*(b1 - b2)*(o1 - o2);
                break;

            case 2:
                va = -(1/4.);
                break;

        }
    }else if ( l1 == 0 && l2 == 2 ){
        switch ( lambda1 ) {
            case 0:
                va =  (1/2.)*(b1 + b2 + 2*b1*b1*(o1 - o2)*(o1 - o2));
                break;

            case 1:
                va =  I*b1*(o1 - o2);

                break;
            case 2:
                va = -(1/4.);

                break;
        }
    }else if ( l1 == 1 && l2 == 2 ){
        switch ( lambda1 ) {
            case 0:
                va =  (-(1/2.))*((-b1)*b2 + b2*b2 + 2*b1*b1*(-1 + b2*(o1 - o2)*(o1 - o2)))*(o1 - o2);
                break;
            case 1:
                va =  (1/4.)*I*(3*b2 + b1*(3 - 4*b2*(o1 - o2)*(o1 - o2)) + 2*b1*b1*(o1 - o2)*(o1 - o2));
                break;
            case 2:
                va = (-(1/4.))*(2*b1 - b2)*(o1 - o2);
                break;
            case 3:
                va = -(I/8.);
                break;

        }
    }else if ( l1 == 2 && l2 == 2 ){
        switch ( lambda1 ) {
            case 0:
                va =  (1/4.)*(-6*b1*b2*(-1 + b2*(o1 - o2)*(o1 - o2)) + b2*b2*(3 + 2*b2*(o1 - o2)*(o1 - o2)) +
                              b1*b1*(3 - 6*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 2*b1*b1*b1*(o1 - o2)*(o1 - o2));
                break;

            case 1:
                va =  (1/2.)*I*(-3*b2*b2 + b1*b1*(3 - 2*b2*(o1 - o2)*(o1 - o2)) + 2*b1*b2*b2*(o1 - o2)*(o1 - o2))*(o1 - o2);
                break;
            case 2:
                va = (1/4.)*((-b2)*(3 + b2*(o1 - o2)*(o1 - o2)) + b1*(-3 + 4*b2*(o1 - o2)*(o1 - o2)) - b1*b1*(o1 - o2)*(o1 - o2));
                break;

            case 3:
                va =  (-(1/4.))*I*(b1 - b2)*(o1 - o2);
                break;

            case 4:
                va = 1./16.;
                break;

        }
    }
        //{(1/2)*b1*(3*b1 + 3*b2 + 2*b1^2*(o1 - o2)^2)*(o1 - o2),
        //(3/4)*I*(b1 + b2 + 2*b1^2*(o1 - o2)^2), (-(3/4))*b1*(o1 - o2), -(I/8)}
        else if ( l1 == 0 && l2 == 3 ){
            switch ( lambda1 ) {
                case 0:
                    va =  (1/2.)*b1*(3*b1 + 3*b2 + 2*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                case 1:
                    va =  (3/4.)*I*(b1 + b2 + 2*b1*b1*(o1 - o2)*(o1 - o2));
                    break;
                case 2:
                    va =  (-(3/4.))*b1*(o1 - o2);
                    break;
                    
                case 3:
                    va =  -(I/8.);
                    break;
                    
                
                    
            }
        }
//    {(1/4)*(3*b1^2 + 3*b2^2 - 6*b1*b2*(-1 + b2*(o1 - o2)^2) -
//            2*b1^3*(-3 + 2*b2*(o1 - o2)^2)*(o1 - o2)^2),
//        (1/4)*I*(6*b1*b2 - 3*b2^2 + b1^2*(9 - 6*b2*(o1 - o2)^2) + 2*b1^3*(o1 - o2)^2)*
//        (o1 - o2), (-(3/4))*(b2 - b1*(-1 + b2*(o1 - o2)^2) + b1^2*(o1 - o2)^2),
//        (-(1/8))*I*(3*b1 - b2)*(o1 - o2), 1/16}
//
    
        else if ( l1 == 1 && l2 == 3 ){
            switch ( lambda1 ) {
                case 0:
                    va = (1/4.)*(3*b1*b1 + 3*b2*b2 - 6*b1*b2*(-1 + b2*(o1 - o2)*(o1-o2)) -
                                        2*b1*b1*b1*(-3 + 2*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 1:
                    va =  (1/4.)*I*(6*b1*b2 - 3*b2*b2 + b1*b1*(9 - 6*b2*(o1 - o2)*(o1-o2)) + 2*b1*b1*b1*(o1 - o2)*(o1 - o2))*
                            (o1 - o2);
                    break;
                case 2:
                    va =  (-(3/4.))*(b2 - b1*(-1 + b2*(o1 - o2)*(o1 - o2)) + b1*b1*(o1 - o2)*(o1 - o2)) ;
                    break;
                    
                case 3:
                    va = (-(1/8.))*I*(3*b1 - b2)*(o1 - o2);
                    break;
                    
                case 4:
                    va = 1/16.;
                    break;
                    
            }
        }
    
//    {(1/4)*(-6*b2*b2*b2 - 6*b1*b1*b2*(-2 + b2*(o1 - o2)*(o1 - o2)) + 3*b1*b2*b2*(-1 + 2*b2*(o1 - o2)*(o1 - o2)) +
//            b1*b1*b1*(9 - 10*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 2*b1*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2),
//        (-(1/8))*I*(30*b1*b2*(-1 + b2*(o1 - o2)*(o1 - o2)) - 3*b2*b2*(5 + 2*b2*(o1 - o2)*(o1 - o2)) -
//                    3*b1*b1*(5 - 6*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 2*b1*b1*b1*(-9 + 4*b2*(o1 - o2)*(o1 - o2))*
//                    (o1 - o2)*(o1 - o2)), (-(1/4))*(-6*b2*b2 + b1*b1*(9 - 6*b2*(o1 - o2)*(o1 - o2)) +
//                                            3*b1*b2*(1 + b2*(o1 - o2)*(o1 - o2)) + b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2),
//        (-(1/8))*I*(b1*(5 - 6*b2*(o1 - o2)*(o1 - o2)) + b2*(5 + b2*(o1 - o2)*(o1 - o2)) + 3*b1*b1*(o1 - o2)*(o1 - o2)),
//        (1/16)*(3*b1 - 2*b2)*(o1 - o2), I/32}

    
        else if ( l1 == 2 && l2 == 3 ){
            switch ( lambda1 ) {
                case 0:
                    va = (1/4.)*(-6*b2*b2*b2 - 6*b1*b1*b2*(-2 + b2*(o1 - o2)*(o1 - o2)) + 3*b1*b2*b2*(-1 + 2*b2*(o1 - o2)*(o1 - o2)) +
                                            b1*b1*b1*(9 - 10*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 2*b1*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                case 1:
                    va =  (-(1/8.))*I*(30*b1*b2*(-1 + b2*(o1 - o2)*(o1 - o2)) - 3*b2*b2*(5 + 2*b2*(o1 - o2)*(o1 - o2)) -
                                                         3*b1*b1*(5 - 6*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 2*b1*b1*b1*(-9 + 4*b2*(o1 - o2)*(o1 - o2))*
                                                          (o1 - o2)*(o1 - o2));
                    break;
                case 2:
                    va =  (-(1/4.))*(-6*b2*b2 + b1*b1*(9 - 6*b2*(o1 - o2)*(o1 - o2)) +
                                                                            3*b1*b2*(1 + b2*(o1 - o2)*(o1 - o2)) + b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                case 3:
                    va = (-(1/8.))*I*(b1*(5 - 6*b2*(o1 - o2)*(o1 - o2)) + b2*(5 + b2*(o1 - o2)*(o1 - o2)) + 3*b1*b1*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 4:
                    va = (1/16.)*(3*b1 - 2*b2)*(o1 - o2);
                    break;
                    
                case 5:
                    va = I/32.;
                    break;

            }
        }
    //    {(1/8)*(15*b2*b2*b2 + 18*b2*b2*b2*b2*o1*1 - 3*b1*b2*b2*(-15 + 6*b2*(o1 - o2)*(o1 - o2)+
    //                                                4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 3*b1*b1*b2*(15 - 24*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
    //            b1*b1*b1*(15 - 18*b2*(o1 - o2)*(o1 - o2) + 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 8*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) -
    //            6*b1*b1*b1*b1*(-3 + 2*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) - 36*b2*b2*b2*b2*o1*o2 + 18*b2*b2*b2*b2*o2*o2),
    //        (3/8)*I*((-b2*b2*b2)*(15 + 2*b2*(o1 - o2)*(o1 - o2)) + b1*b2*b2*(-15 + 16*b2*(o1 - o2)*(o1 - o2)) +
    //                 b1*b1*b1*(15 - 16*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
    //                 b1*b1*(15*b2 - 4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 2*b1*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2),
    //        (1/16)*(-9*b2*b2*(5 + 4*b2*(o1 - o2)*(o1 - o2)) + 6*b1*b2*(-15 + 12*b2*(o1 - o2)*(o1 - o2) +
    //                                                          2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) - 9*b1*b1*(5 - 8*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
    //                12*b1*b1*b1*(-3 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)), (-(1/8))*I*(b1 - b2)*
    //        (b1*(15 - 8*b2*(o1 - o2)*(o1 - o2)) + b2*(15 + b2*(o1 - o2)*(o1 - o2)) + b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2),
    //        (3/32)*(b1*(5 - 6*b2*(o1 - o2)*(o1 - o2)) + b2*(5 + 2*b2*(o1 - o2)*(o1 - o2)) + 2*b1*b1*(o1 - o2)*(o1 - o2)),
    //        (3/32)*I*(b1 - b2)*(o1 - o2), -(1/64)}
        else if ( l1 == 3 && l2 == 3 ){
            switch ( lambda1 ) {
                case 0:
                    va = (1/8.)*(15*b2*b2*b2 + 18*b2*b2*b2*b2*o1*1 - 3*b1*b2*b2*(-15 + 6*b2*(o1 - o2)*(o1 - o2)+
                                                                                                                                4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 3*b1*b1*b2*(15 - 24*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                                                                            b1*b1*b1*(15 - 18*b2*(o1 - o2)*(o1 - o2) + 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 8*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) -
                                                                                            6*b1*b1*b1*b1*(-3 + 2*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) - 36*b2*b2*b2*b2*o1*o2 + 18*b2*b2*b2*b2*o2*o2);
                    break;
                    
                case 1:
                    va = (3/8.)*I*((-b2*b2*b2)*(15 + 2*b2*(o1 - o2)*(o1 - o2)) + b1*b2*b2*(-15 + 16*b2*(o1 - o2)*(o1 - o2)) +
                                                   b1*b1*b1*(15 - 16*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                                   b1*b1*(15*b2 - 4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 2*b1*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2) ;
                    break;
                case 2:
                    va =  (1/16.)*(-9*b2*b2*(5 + 4*b2*(o1 - o2)*(o1 - o2)) + 6*b1*b2*(-15 + 12*b2*(o1 - o2)*(o1 - o2) +
                                                                                                                                              2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) - 9*b1*b1*(5 - 8*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                                                                                     12*b1*b1*b1*(-3 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 3:
                    va =(-(1/8.))*I*(b1 - b2)*
                            (b1*(15 - 8*b2*(o1 - o2)*(o1 - o2)) + b2*(15 + b2*(o1 - o2)*(o1 - o2)) + b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2) ;
                    break;
                    
                case 4:
                    va = (3/32.)*(b1*(5 - 6*b2*(o1 - o2)*(o1 - o2)) + b2*(5 + 2*b2*(o1 - o2)*(o1 - o2)) + 2*b1*b1*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 5:
                    va = (3/32.)*I*(b1 - b2)*(o1 - o2);
                    break;
                    
                case 6:
                    va = -(1/64.);
                    break;
                    
            }
        }

//    {(3*b1*b2)/2 + (3*b2*b2)/4 + (3/4)*b1*b1*(1 + 4*b2*(o1 - o2)*(o1 - o2)) + 3*b1*b1*b1*(o1 - o2)*(o1 - o2) +
//        b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2), I*b1*(3*b1 + 3*b2 + 2*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2),
//        (-(3/4))*(b1 + b2 + 2*b1*b1*(o1 - o2)*(o1 - o2)), (-(1/2))*I*b1*(o1 - o2), 1/16}
    
        else if ( l1 == 0 && l2 == 4 ){
            switch ( lambda1 ) {
                case 0:
                    va = (3*b1*b2)/2. + (3*b2*b2)/4. + (3/4.)*b1*b1*(1 + 4*b2*(o1 - o2)*(o1 - o2)) + 3*b1*b1*b1*(o1 - o2)*(o1 - o2) +
                    b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2);
                    break;
                    
                case 1:
                    va = I*b1*(3*b1 + 3*b2 + 2*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                case 2:
                    va = (-(3/4.))*(b1 + b2 + 2*b1*b1*(o1 - o2)*(o1 - o2)) ;
                    break;
                    
                case 3:
                    va =(-(1/2.))*I*b1*(o1 - o2);
                    break;
                    
                case 4:
                    va =1/16.;
                    break;
            }
        }
    
//    {(-(1/4))*(-6*b1*b2*b2 + 3*b2*b2*b2 + 4*b1*b1*b1*(-3 + b2*(o1 - o2)*(o1 - o2)) +
//               3*b1*b1*b2*(-7 + 4*b2*(o1 - o2)) + 4*b1*b1*b1*b1*(-2 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2))*
//        (o1 - o2), (1/8)*I*(15*b2*b2 - 6*b1*b2*(-5 + 4*b2*(o1 - o2)*(o1 - o2)) +
//                            3*b1*b1*(5 + 4*b2*(o1 - o2)*(o1 - o2)) - 4*b1*b1*b1*(-9 + 4*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) +
//                            4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)), (-(1/4))*(9*b1*b2 - 3*b2*b2 - 6*b1*b1*(-2 + b2*(o1 - o2)*(o1 - o2)) +
//                                                           4*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2), (-(1/8))*I*(5*b2 + b1*(5 - 4*b2*(o1 - o2)*(o1 - o2)) +
//                                                                                                      6*b1*b1*(o1 - o2)*(o1 - o2)), (1/16)*(4*b1 - b2)*(o1 - o2), I/32}
        else if ( l1 == 1 && l2 == 4 ){
            switch ( lambda1 ) {
                case 0:
                    va = (-(1/4.))*(-6*b1*b2*b2 + 3*b2*b2*b2 + 4*b1*b1*b1*(-3 + b2*(o1 - o2)*(o1 - o2)) +
                                   3*b1*b1*b2*(-7 + 4*b2*(o1 - o2)*(o1 - o2)) + 4*b1*b1*b1*b1*(-2 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2))*
                    (o1 - o2);
                    break;
                    
                case 1:
                    va = (1/8.)*I*(15*b2*b2 - 6*b1*b2*(-5 + 4*b2*(o1 - o2)*(o1 - o2)) +
                                  3*b1*b1*(5 + 4*b2*(o1 - o2)*(o1 - o2)) - 4*b1*b1*b1*(-9 + 4*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) +
                                  4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2));
                    break;
                case 2:
                    va =   (-(1/4.))*(9*b1*b2 - 3*b2*b2 - 6*b1*b1*(-2 + b2*(o1 - o2)*(o1 - o2)) +
                                     4*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                case 3:
                    va = (-(1/8.))*I*(5*b2 + b1*(5 - 4*b2*(o1 - o2)*(o1 - o2)) +
                                     6*b1*b1*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 4:
                    va =(1/16.)*(4*b1 - b2)*(o1 - o2);
                    break;
                    
                case 5:
                    va =I/32.;
                    break;

            }
        }
//    {(1/8)*(3*b2*b2*b2*(5 + 2*b2*(o1 - o2)*(o1 - o2)) - 9*b1*b2*b2*(-5 + 4*b2*(o1 - o2)*(o1 - o2)) +
//            b1*b1*b1*(15 + 24*b2*(o1 - o2)*(o1 - o2) - 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//            3*b1*b1*b2*(15 - 18*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//            4*b1*b1*b1*b1*(9 - 7*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) +
//            4*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)), (-(1/4))*I*(15*b2*b2*b2 + 3*b1*b1*b2*(-15 + 8*b2*(o1 - o2)*(o1 - o2)) +
//                                             b1*b1*b1*(-30 + 24*b2*(o1 - o2)*(o1 - o2) - 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) - 12*b1*b2*b2*b2*(o1 - o2)*(o1 - o2) +
//                                             4*b1*b1*b1*b1*(-3 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2))*(o1 - o2),
//        (1/16)*(-3*b2*b2*(15 + 4*b2*(o1 - o2)*(o1 - o2)) + 6*b1*b2*(-15 + 14*b2*(o1 - o2)*(o1 - o2)) -
//                3*b1*b1*(15 - 8*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                8*b1*b1*b1*(-9 + 4*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) - 4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)),
//        (-(1/4))*I*(-5*b2*b2 + b1*b2*(5 + 2*b2*(o1 - o2)*(o1 - o2)) - 2*b1*b1*(-5 + 3*b2*(o1 - o2)*(o1 - o2)) +
//                    2*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2), (1/32)*(b1*(15 - 16*b2*(o1 - o2)*(o1 - o2)) +
//                                                           b2*(15 + 2*b2*(o1 - o2)*(o1 - o2)) + 12*b1*b1*(o1 - o2)*(o1 - o2)), (1/16)*I*(2*b1 - b2)*(o1 - o2),
//        -(1/64)}
    
        else if ( l1 == 2 && l2 == 4 ){
            switch ( lambda1 ) {
                case 0:
                    va = (1/8.)*(3*b2*b2*b2*(5 + 2*b2*(o1 - o2)*(o1 - o2)) - 9*b1*b2*b2*(-5 + 4*b2*(o1 - o2)*(o1 - o2)) +
                                b1*b1*b1*(15 + 24*b2*(o1 - o2)*(o1 - o2) - 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                3*b1*b1*b2*(15 - 18*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                4*b1*b1*b1*b1*(9 - 7*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) +
                                4*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 1:
                    va = (-(1/4.))*I*(15*b2*b2*b2 + 3*b1*b1*b2*(-15 + 8*b2*(o1 - o2)*(o1 - o2)) +
                                     b1*b1*b1*(-30 + 24*b2*(o1 - o2)*(o1 - o2) - 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) - 12*b1*b2*b2*b2*(o1 - o2)*(o1 - o2) +
                                     4*b1*b1*b1*b1*(-3 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                case 2:
                    va = (1/16.)*(-3*b2*b2*(15 + 4*b2*(o1 - o2)*(o1 - o2)) + 6*b1*b2*(-15 + 14*b2*(o1 - o2)*(o1 - o2)) -
                                 3*b1*b1*(15 - 8*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                 8*b1*b1*b1*(-9 + 4*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) - 4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) ;
                    break;
                    
                case 3:
                    va = (-(1/4.))*I*(-5*b2*b2 + b1*b2*(5 + 2*b2*(o1 - o2)*(o1 - o2)) - 2*b1*b1*(-5 + 3*b2*(o1 - o2)*(o1 - o2)) +
                                     2*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                case 4:
                    va = (1/32.)*(b1*(15 - 16*b2*(o1 - o2)*(o1 - o2)) +
                                 b2*(15 + 2*b2*(o1 - o2)*(o1 - o2)) + 12*b1*b1*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 5:
                    va = (1/16.)*I*(2*b1 - b2)*(o1 - o2);
                    break;
                    
                case 6:
                    va =  -(1/64.);
                    break;
            }
        }
//    {(1/8)*(-45*b2*b2*b2*b2*o1 + 15*b1*b2*b2*b2*(-5 + 4*b2*(o1 - o2)*(o1 - o2))*(o1 - o2) +
//            3*b1*b1*b1*b2*(45 - 40*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) -
//            4*b1*b1*b1*b1*(-15 + 15*b2*(o1 - o2)*(o1 - o2) - 9*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*
//            (o1 - o2) - 6*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) - 12*b1*b1*b1*b1*b1*(-2 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//            45*b2*b2*b2*b2*o2 - 3*b1*b1*b2*b2*(-15*o1 - 10*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                                      15*o2)), (1/16)*I*(105*b2*b2*b2 - 3*b1*b2*b2*(-105 + 60*b2*(o1 - o2)*(o1 - o2) +
//                                                                               16*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 6*b1*b1*((105*b2)/2 - 75*b2*b2*(o1 - o2)*(o1 - o2) +
//                                                                                                              28*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) - b1*b1*b1*(-105 - 72*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 32*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                                                         90*b2*b2*b2*b2*(o1 - o2)*(o1 - o2) + 12*b1*b1*b1*b1*(15 - 11*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*
//                                                         (o1 - o2)*(o1 - o2) + 12*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)),
//        (3/16)*(2*b1*b2*b2*(15 - 22*b2*(o1 - o2)*(o1 - o2)) + b2*b2*b2*(45 + 4*b2*(o1 - o2)*(o1 - o2)) -
//                4*b1*b1*b1*(15 - 14*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                b1*b1*b2*(-75 + 24*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                4*b1*b1*b1*b1*(-4 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2))*(o1 - o2),
//        (-(1/32))*I*(15*b2*b2*(7 + 4*b2*(o1 - o2)*(o1 - o2)) -
//                     2*b1*b2*(-105 + 90*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                     3*b1*b1*(35 - 40*b2*(o1 - o2)*(o1 - o2) + 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) -
//                     24*b1*b1*b1*(-5 + 2*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) + 4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)),
//        (1/32)*((-b2*b2)*(45 + 2*b2*(o1 - o2)*(o1 - o2)) - 12*b1*b1*(-5 + 3*b2*(o1 - o2)*(o1 - o2)) +
//                3*b1*b2*(5 + 8*b2*(o1 - o2)*(o1 - o2)) + 8*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2),
//        (3/64)*I*(b1*(7 - 8*b2*(o1 - o2)*(o1 - o2)) + b2*(7 + 2*b2*(o1 - o2)*(o1 - o2)) +
//                  4*b1*b1*(o1 - o2)*(o1 - o2)), (-(1/64))*(4*b1 - 3*b2)*(o1 - o2), -(I/128)}
        else if ( l1 == 3 && l2 == 4 ){
            switch ( lambda1 ) {
                case 0:
                    va = (1/8.)*(-45*b2*b2*b2*b2*o1 + 15*b1*b2*b2*b2*(-5 + 4*b2*(o1 - o2)*(o1 - o2))*(o1 - o2) +
                                3*b1*b1*b1*b2*(45 - 40*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) -
                                4*b1*b1*b1*b1*(-15 + 15*b2*(o1 - o2)*(o1 - o2) - 9*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*
                                (o1 - o2) - 6*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) - 12*b1*b1*b1*b1*b1*(-2 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2) +
                                45*b2*b2*b2*b2*o2 - 3*b1*b1*b2*b2*(-15*o1 - 10*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
                                                                   15*o2));
                    break;
                    
                case 1:
                    va = (1/16.)*I*(105*b2*b2*b2 - 3*b1*b2*b2*(-105 + 60*b2*(o1 - o2)*(o1 - o2) +
                                                              16*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 6*b1*b1*((105*b2)/2 - 75*b2*b2*(o1 - o2)*(o1 - o2) +
                                                                                                                           28*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) - b1*b1*b1*(-105 - 72*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 32*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                   90*b2*b2*b2*b2*(o1 - o2)*(o1 - o2) + 12*b1*b1*b1*b1*(15 - 11*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*
                                   (o1 - o2)*(o1 - o2) + 12*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2));
                    break;
                case 2:
                    va =
                    (3/16.)*(2*b1*b2*b2*(15 - 22*b2*(o1 - o2)*(o1 - o2)) + b2*b2*b2*(45 + 4*b2*(o1 - o2)*(o1 - o2)) -
                            4*b1*b1*b1*(15 - 14*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                            b1*b1*b2*(-75 + 24*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                            4*b1*b1*b1*b1*(-4 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                    
                case 3:
                   va = (-(1/32.))*I*(15*b2*b2*(7 + 4*b2*(o1 - o2)*(o1 - o2)) -
                                  2*b1*b2*(-105 + 90*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                  3*b1*b1*(35 - 40*b2*(o1 - o2)*(o1 - o2) + 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) -
                                  24*b1*b1*b1*(-5 + 2*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) + 4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2));
                    break;
                    
                    
                case 4:
                    va = (1/32.)*((-b2*b2)*(45 + 2*b2*(o1 - o2)*(o1 - o2)) - 12*b1*b1*(-5 + 3*b2*(o1 - o2)*(o1 - o2)) +
                                 3*b1*b2*(5 + 8*b2*(o1 - o2)*(o1 - o2)) + 8*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    
                break;
                    
                case 5:
                    va =(3/64.)*I*(b1*(7 - 8*b2*(o1 - o2)*(o1 - o2)) + b2*(7 + 2*b2*(o1 - o2)*(o1 - o2)) +
                                  4*b1*b1*(o1 - o2)*(o1 - o2)) ;
                    break;
                    
                case 6:
                    va = (-(1/64.))*(4*b1 - 3*b2)*(o1 - o2);
                    break;
                    
                case 7:
                    va =  -(I/128.);
                    break;
            }
        }
    
    
//    {(1/32)*(210*b2*b2*b2*b2 + 24*b1*b2*b2*b2*(35 + 5*b2*(o1 - o2)*(o1 - o2) - 14*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//             40*b1*b1*b1*b2*(21 - 36*b2*(o1 - o2)*(o1 - o2) + 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//             12*b1*b1*b2*b2*(105 - 120*b2*(o1 - o2)*(o1 - o2) + 10*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                           8*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 2*b1*b1*b1*b1*(105 + 60*b2*(o1 - o2)*(o1 - o2) + 60*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                                                         80*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 16*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 360*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2) +
//             24*b1*b1*b1*b1*b1*(15 - 14*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) +
//             24*b1*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 24*b2*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)),
//        (-(1/4))*I*(105*b2*b2*b2*b2*(o1 - o2) + 2*b1*b1*b1*b2*(-105 + 75*b2*(o1 - o2)*(o1 - o2) -
//                                                    4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) + b1*b1*b1*b1*(-105 + 120*b2*(o1 - o2)*(o1 - o2) -
//                                                                                          60*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 8*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) + 30*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                    60*b1*b1*b2*b2*b2*b2*(5/2 - b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2) + 6*b1*b1*b1*b1*b1*(-5 + 2*b2*(o1 - o2)*(o1 - o2))*
//                    (o1 - o2)*(o1 - o2)*(o1 - o2) - 6*b1*b2*b2*b2*(20*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 35*(-o1 + o2))),
//        (1/16)*(-210*b2*b2*b2 - 270*b2*b2*b2*b2*o1^2 - 6*b1*b1*b1*b1*(45 - 30*b2*(o1 - o2)*(o1 - o2) +
//                                                    4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) - 12*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 12*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                540*b2*b2*b2*b2*o1*o2 - 270*b2*b2*b2*b2*o2^2 + 90*b1*b2*b2*(-7 + 2*b2*o1^2 + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                                                             4*b2*o1*o2 + 2*b2*o2^2) + 6*b1*b1*b2*(-105 + 150*b2*o1^2 - 40*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                                                                                                  4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 300*b2*o1*o2 + 150*b2*o2^2) +
//                b1*b1*b1*(-210 + 180*b2*o1^2 - 240*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 64*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                      360*b2*o1*o2 + 180*b2*o2^2)), (1/8)*I*(5*b2*b2*b2*(21 + 4*b2*(o1 - o2)*(o1 - o2)) +
//                                                             b1*b1*b1*(-105 + 100*b2*(o1 - o2)*(o1 - o2) - 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                                                             b1*b2*b2*(105 - 100*b2*(o1 - o2)*(o1 - o2) - 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                                                             3*b1*b1*b2*(-35 + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 4*b1*b1*b1*b1*(-5 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2))*
//        (o1 - o2), (1/32)*(b2*b2*(105 + 90*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) -
//                           2*b1*b2*(-105 + 75*b2*(o1 - o2)*(o1 - o2) + 16*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                           3*b1*b1*(35 - 50*b2*(o1 - o2)*(o1 - o2) + 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) -
//                           2*b1*b1*b1*(-45 + 16*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) + 2*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)),
//        (1/16)*I*(b1 - b2)*(b1*(21 - 10*b2*(o1 - o2)*(o1 - o2)) + b2*(21 + 2*b2*(o1 - o2)*(o1 - o2)) +
//                            2*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2), (1/32)*((-b2)*(7 + 3*b2*(o1 - o2)*(o1 - o2)) +
//                                                                   b1*(-7 + 8*b2*(o1 - o2)*(o1 - o2)) - 3*b1*b1*(o1 - o2)*(o1 - o2)), (-(1/32))*I*(b1 - b2)*(o1 - o2),
//        1/256}
        else if ( l1 == 4 && l2 == 4 ){
            switch ( lambda1 ) {
                case 0:
                    va = (1/32.)*(210*b2*b2*b2*b2 + 24*b1*b2*b2*b2*(35 + 5*b2*(o1 - o2)*(o1 - o2) - 14*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                 40*b1*b1*b1*b2*(21 - 36*b2*(o1 - o2)*(o1 - o2) + 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                 12*b1*b1*b2*b2*(105 - 120*b2*(o1 - o2)*(o1 - o2) + 10*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
                                                 8*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 2*b1*b1*b1*b1*(105 + 60*b2*(o1 - o2)*(o1 - o2) + 60*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
                                                                                                                                          80*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 16*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 360*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2) +
                                 24*b1*b1*b1*b1*b1*(15 - 14*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) +
                                 24*b1*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 24*b2*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 1:
                    va =  (-(1/4.))*I*(105*b2*b2*b2*b2*(o1 - o2) + 2*b1*b1*b1*b2*(-105 + 75*b2*(o1 - o2)*(o1 - o2) -
                                                                                 4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) + b1*b1*b1*b1*(-105 + 120*b2*(o1 - o2)*(o1 - o2) -
                                                                                                                                                                                  60*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 8*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) + 30*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) -
                                      60*b1*b1*b2*b2*b2*b2*(5/2. - b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2) + 6*b1*b1*b1*b1*b1*(-5 + 2*b2*(o1 - o2)*(o1 - o2))*
                                      (o1 - o2)*(o1 - o2)*(o1 - o2) - 6*b1*b2*b2*b2*(20*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 35*(-o1 + o2)));
                    break;
                case 2:
                    va = (1/16.)*(-210*b2*b2*b2 - 270*b2*b2*b2*b2*o1*o1 - 6*b1*b1*b1*b1*(45 - 30*b2*(o1 - o2)*(o1 - o2) +
                                                                                       4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) - 12*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 12*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
                                 540*b2*b2*b2*b2*o1*o2 - 270*b2*b2*b2*b2*o2*o2 + 90*b1*b2*b2*(-7 + 2*b2*o1*o1 + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
                                                                                             4*b2*o1*o2 + 2*b2*o2*o2) + 6*b1*b1*b2*(-105 + 150*b2*o1*o1 - 40*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
                                                                                                                                   4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 300*b2*o1*o2 + 150*b2*o2*o2) +
                                 b1*b1*b1*(-210 + 180*b2*o1*o1 - 240*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 64*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
                                           360*b2*o1*o2 + 180*b2*o2*o2)) ;
                    break;
                    
                case 3:
                    va =(1/8.)*I*(5*b2*b2*b2*(21 + 4*b2*(o1 - o2)*(o1 - o2)) +
                                 b1*b1*b1*(-105 + 100*b2*(o1 - o2)*(o1 - o2) - 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                 b1*b2*b2*(105 - 100*b2*(o1 - o2)*(o1 - o2) - 4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                 3*b1*b1*b2*(-35 + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) + 4*b1*b1*b1*b1*(-5 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2))*
                    (o1 - o2) ;
                    break;
                    
                case 4:
                    va = (1/32.)*(b2*b2*(105 + 90*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) -
                                 2*b1*b2*(-105 + 75*b2*(o1 - o2)*(o1 - o2) + 16*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                 3*b1*b1*(35 - 50*b2*(o1 - o2)*(o1 - o2) + 24*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) -
                                 2*b1*b1*b1*(-45 + 16*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) + 2*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 5:
                    va = (1/16.)*I*(b1 - b2)*(b1*(21 - 10*b2*(o1 - o2)*(o1 - o2)) + b2*(21 + 2*b2*(o1 - o2)*(o1 - o2)) +
                                             2*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                case 6:
                    va = (1/32.)*((-b2)*(7 + 3*b2*(o1 - o2)*(o1 - o2)) +
                                 b1*(-7 + 8*b2*(o1 - o2)*(o1 - o2)) - 3*b1*b1*(o1 - o2)*(o1 - o2));
                    break;
                case 7:
                    va = (-(1/32.))*I*(b1 - b2)*(o1 - o2);
                    break;
                    
                case 8:
                    va = 1/256.;
                    break;
            
            }
        }
//    {(1/4)*b1*(30*b1*b2 + 15*b2*b2 + 5*b1*b1*(3 + 4*b2*(o1 - o2)*(o1 - o2)) + 20*b1*b1*b1*(o1 - o2)*(o1 - o2) +
//               4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2),
//        (5/8)*I*(6*b1*b2 + 3*b2*b2 + 3*b1*b1*(1 + 4*b2*(o1 - o2)*(o1 - o2)) + 12*b1*b1*b1*(o1 - o2*(o1 - o2) +
//                 4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)), (-(5/4))*b1*(3*b1 + 3*b2 + 2*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2),
//        (-(5/8))*I*(b1 + b2 + 2*b1*b1*(o1 - o2)*(o1 - o2)), (5/16)*b1*(o1 - o2), I/32}
        else if ( l1 == 0 && l2 == 5 ){
            switch ( lambda1 ) {
                case 0:
                    va = (1/4.)*b1*(30*b1*b2 + 15*b2*b2 + 5*b1*b1*(3 + 4*b2*(o1 - o2)*(o1 - o2)) + 20*b1*b1*b1*(o1 - o2)*(o1 - o2) +
                                   4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                case 1:
                    va = (5/8.)*I*(6*b1*b2 + 3*b2*b2 + 3*b1*b1*(1 + 4*b2*(o1 - o2)*(o1 - o2)) + 12*b1*b1*b1*(o1 - o2)*(o1 - o2) +4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2));
                    break;
                case 2:
                    va = (-(5/4.))*b1*(3*b1 + 3*b2 + 2*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                case 3:
                    va =(-(5/8.))*I*(b1 + b2 + 2*b1*b1*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 4:
                    va =(5/16.)*b1*(o1 - o2);
                    break;
                    
                case 5:
                    va =I/32.;
                    break;
                    
            }
        }
                                  
////    List((45*Power(b1,2)*b2 + 15*Power(b2,3) -
////          15*b1*Power(b2,2)*(-3 + 2*b2*Power(o1 - o2,2)) -
////          5*Power(b1,3)*(-3 - 18*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4)) -
////          20*Power(b1,4)*(-3 + b2*Power(o1 - o2,2))*Power(o1 - o2,2) -
////          4*Power(b1,5)*(-5 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,4))/8.
//
////    ,
//         Complex(0,0.125)*(75*Power(b1,3) + 45*b1*Power(b2,2) - 15*Power(b2,3) -
//                           15*Power(b1,2)*b2*(-9 + 4*b2*Power(o1 - o2,2)) -
//                           20*Power(b1,4)*(-3 + b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                           4*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2)
//
////    ,
//         (-5*(9*Power(b2,2) - 6*b1*b2*(-3 + 2*b2*Power(o1 - o2,2)) +
//              3*Power(b1,2)*(3 + 4*b2*Power(o1 - o2,2)) -
//              8*Power(b1,3)*(-3 + b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//              4*Power(b1,4)*Power(o1 - o2,4)))/16.
//
////    ,
//         Complex(0,-0.625)*(4*b1*b2 - Power(b2,2) + Power(b1,2)*(5 - 2*b2*Power(o1 - o2,2)) +
//                            2*Power(b1,3)*Power(o1 - o2,2))*(o1 - o2)
//
////    ,
//         (5*(3*b2 + b1*(3 - 2*b2*Power(o1 - o2,2)) + 4*Power(b1,2)*Power(o1 - o2,2)))/32.
//
////    ,
//         Complex(0,0.03125)*(5*b1 - b2)*(o1 - o2),
//
// //  -0.015625
    
                                  
    
        else if ( l1 == 1 && l2 == 5 ){
            switch ( lambda1 ) {
                case 0:
                    va = (45*Power(b1,2)*b2 + 15*Power(b2,3) -
                                    15*b1*Power(b2,2)*(-3 + 2*b2*Power(o1 - o2,2)) -
                                    5*Power(b1,3)*(-3 - 18*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4)) -
                                   20*Power(b1,4)*(-3 + b2*Power(o1 - o2,2))*Power(o1 - o2,2) -
                                    4*Power(b1,5)*(-5 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,4))/8.;
                    break;
                    
                case 1:
                    va = Complex(0,0.125)*(75*Power(b1,3) + 45*b1*Power(b2,2) - 15*Power(b2,3) -
                                           15*Power(b1,2)*b2*(-9 + 4*b2*Power(o1 - o2,2)) -
                                           20*Power(b1,4)*(-3 + b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
                                           4*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2);
                    break;
                case 2:
                    va =(-5*(9*Power(b2,2) - 6*b1*b2*(-3 + 2*b2*Power(o1 - o2,2)) +
                                           3*Power(b1,2)*(3 + 4*b2*Power(o1 - o2,2)) -
                                    8*Power(b1,3)*(-3 + b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
                                           4*Power(b1,4)*Power(o1 - o2,4)))/16;
                    break;
                    
                case 3:
                    va = Complex(0,-0.625)*(4*b1*b2 - Power(b2,2) + Power(b1,2)*(5 - 2*b2*Power(o1 - o2,2)) +
                                            2*Power(b1,3)*Power(o1 - o2,2))*(o1 - o2);
                    break;
                    
                case 4:
                    va =(5*(3*b2 + b1*(3 - 2*b2*Power(o1 - o2,2)) + 4*Power(b1,2)*Power(o1 - o2,2)))/32.;
                    break;
                    
                case 5:
                    va =Complex(0,0.03125)*(5*b1 - b2)*(o1 - o2);
                    break;
                   
                case 6:
                    va =  -0.015625;
                    break;

            }
        }
//    {(1/8)*(-15*b1*b2*b2*b2*o1 - 30*b2*b2*b2*b2*o1 + 30*b1*b2*b2*b2*b2*o1*o1*o1 +
//            15*b1*b1*b2*b2*(9*o1 - 4*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) - 9*o2) + 75*b1*b1*b1*b1*(o1 - o2) +
//            5*b1*b1*b1*b2*(39 - 30*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) +
//            4*b1*b1*b1*b1*b1*(15 - 9*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//            4*b1*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 15*b1*b2*b2*b2*o2 + 30*b2*b2*b2*b2*o2 - 90*b1*b2*b2*b2*b2*o1*o1*o2 +
//            90*b1*b2*b2*b2*b2*o1*o2*o2 - 30*b1*b2*b2*b2*b2*o2*o2*o2),
//        (1/16)*I*(15*b2*b2*b2*(7 + 2*b2*(o1 - o2)*(o1 - o2)) - 15*b1*b2*b2*(-21 + 16*b2*(o1 - o2)*(o1 - o2)) -
//                  15*b1*b1*b1*(-7 - 20*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                  30*b1*b1*((21*b2)/2 - 9*b2*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
//                  20*b1*b1*b1*b1*(15 - 9*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) -
//                  4*b1*b1*b1*b1*b1*(-15 + 4*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)),
//        (1/16)*(-45*b1*b2*b2*o1 + 90*b2*b2*b2*o1 - 60*b1*b2*b2*b2*o1*o1*o1 +
//                180*b1*b1*b2*(-2 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2) +
//                5*b1*b1*b1*(-45 + 24*b2*(o1 - o2)*(o1 - o2) - 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) +
//                40*b1*b1*b1*b1*(-3 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2) - 4*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 45*b1*b2*b2*o2 -
//                90*b2*b2*b2*o2 + 180*b1*b2*b2*b2*o1*o1*o2 - 180*b1*b2*b2*b2*o1*o2*o2 + 60*b1*b2*b2*b2*o2*o2*o2),
//        (-(5/32))*I*(b2*b2*(21 + 4*b2*(o1 - o2)*(o1 - o2)) - 6*b1*b2*(-7 + 6*b2*(o1 - o2)*(o1 - o2)) +
//                     b1*b1*(21 + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) - 8*b1*b1*b1*(-5 + 2*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) +
//                     4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)), (5/32)*(-6*b2*b2 + b1*b1*(15 - 8*b2*(o1 - o2)*(o1 - o2)) +
//                                                  b1*b2*(9 + 2*b2*(o1 - o2)*(o1 - o2)) + 4*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2),
//        (1/64)*I*(b1*(21 - 20*b2*(o1 - o2)*(o1 - o2)) + b2*(21 + 2*b2*(o1 - o2)*(o1 - o2)) +
//                  20*b1*b1*(o1 - o2)*(o1 - o2)), (-(1/64))*(5*b1 - 2*b2)*(o1 - o2), -(I/128)}
        else if ( l1 == 2 && l2 == 5 ){
            switch ( lambda1 ) {
                case 0:
                    va = (1/8.)*(-15*b1*b2*b2*b2*o1 - 30*b2*b2*b2*b2*o1 + 30*b1*b2*b2*b2*b2*o1*o1*o1 +
                                15*b1*b1*b2*b2*(9*o1 - 4*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) - 9*o2) + 75*b1*b1*b1*b1*(o1 - o2) +
                                5*b1*b1*b1*b2*(39 - 30*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) +
                                4*b1*b1*b1*b1*b1*(15 - 9*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2) +
                                4*b1*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 15*b1*b2*b2*b2*o2 + 30*b2*b2*b2*b2*o2 - 90*b1*b2*b2*b2*b2*o1*o1*o2 +
                                90*b1*b2*b2*b2*b2*o1*o2*o2 - 30*b1*b2*b2*b2*b2*o2*o2*o2);
                    break;
                    
                case 1:
                    va =  (1/16.)*I*(15*b2*b2*b2*(7 + 2*b2*(o1 - o2)*(o1 - o2)) - 15*b1*b2*b2*(-21 + 16*b2*(o1 - o2)*(o1 - o2)) -
                                    15*b1*b1*b1*(-7 - 20*b2*(o1 - o2)*(o1 - o2) + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                    30*b1*b1*((21*b2)/2 - 9*b2*b2*(o1 - o2)*(o1 - o2) + 4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) +
                                    20*b1*b1*b1*b1*(15 - 9*b2*(o1 - o2)*(o1 - o2) + 2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) -
                                    4*b1*b1*b1*b1*b1*(-15 + 4*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2));
                    break;
                case 2:
                    va =(1/16.)*(-45*b1*b2*b2*o1 + 90*b2*b2*b2*o1 - 60*b1*b2*b2*b2*o1*o1*o1 +
                                180*b1*b1*b2*(-2 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2) +
                                5*b1*b1*b1*(-45 + 24*b2*(o1 - o2)*(o1 - o2) - 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2))*(o1 - o2) +
                                40*b1*b1*b1*b1*(-3 + b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2)*(o1 - o2) - 4*b1*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 45*b1*b2*b2*o2 -
                                90*b2*b2*b2*o2 + 180*b1*b2*b2*b2*o1*o1*o2 - 180*b1*b2*b2*b2*o1*o2*o2 + 60*b1*b2*b2*b2*o2*o2*o2);
                    break;
                    
                case 3:
                    va =(-(5/32.))*I*(b2*b2*(21 + 4*b2*(o1 - o2)*(o1 - o2)) - 6*b1*b2*(-7 + 6*b2*(o1 - o2)*(o1 - o2)) +
                                     b1*b1*(21 + 8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)) - 8*b1*b1*b1*(-5 + 2*b2*(o1 - o2)*(o1 - o2))*(o1 - o2)*(o1 - o2) +
                                     4*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2));
                    break;
                    
                case 4:
                    va =(5/32.)*(-6*b2*b2 + b1*b1*(15 - 8*b2*(o1 - o2)*(o1 - o2)) +
                                b1*b2*(9 + 2*b2*(o1 - o2)*(o1 - o2)) + 4*b1*b1*b1*(o1 - o2)*(o1 - o2))*(o1 - o2);
                    break;
                    
                case 5:
                    va =(1/64.)*I*(b1*(21 - 20*b2*(o1 - o2)*(o1 - o2)) + b2*(21 + 2*b2*(o1 - o2)*(o1 - o2)) +
                                  20*b1*b1*(o1 - o2)*(o1 - o2));
                    break;
                case 6:
                    va =(-(1/64.))*(5*b1 - 2*b2)*(o1 - o2);
                    break;
                    
                case 7:
                    va =-(I/128.);
                    break;

            }
        }
//    List(
//
//         (105*Power(b2,4) + 90*Power(b2,5)*Power(o1,2) -
//          60*b1*Power(b2,3)*(-7 + 3*b2*Power(o1 - o2,2) + Power(b2,2)*Power(o1 - o2,4)) +
//          30*Power(b1,2)*Power(b2,2)*(21 - 26*b2*Power(o1 - o2,2) +
//                                      8*Power(b2,2)*Power(o1 - o2,4)) -
//          20*Power(b1,3)*b2*(-21 + 18*b2*Power(o1 - o2,2) -
//                             15*Power(b2,2)*Power(o1 - o2,4) + 4*Power(b2,3)*Power(o1 - o2,6)) +
//          5*Power(b1,4)*(21 + 90*b2*Power(o1 - o2,2) - 60*Power(b2,2)*Power(o1 - o2,4) +
//                         8*Power(b2,3)*Power(o1 - o2,6)) -
//          4*Power(b1,5)*(-75 + 60*b2*Power(o1 - o2,2) - 24*Power(b2,2)*Power(o1 - o2,4) +
//                         4*Power(b2,3)*Power(o1 - o2,6))*Power(o1 - o2,2) -
//          12*Power(b1,6)*(-5 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,4) -
//          180*Power(b2,5)*o1*o2 + 90*Power(b2,5)*Power(o2,2))/16.
//
//         ,
//
//
//         Complex(0,-0.0625)*(315*Power(b2,4)*(o1 - o2) +
//                             5*Power(b1,4)*(-105 + 60*b2*Power(o1 - o2,2) - 36*Power(b2,2)*Power(o1 - o2,4) +
//                                            8*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) + 30*Power(b2,5)*Power(o1 - o2,3) -
//                             12*Power(b1,5)*(25 - 14*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4))*
//                             Power(o1 - o2,3) - 12*Power(b1,6)*Power(o1 - o2,5) +
//                             30*Power(b1,2)*Power(b2,2)*(b2*Power(o1 - o2,3) + 4*Power(b2,2)*Power(o1 - o2,5) +
//                                                         21*(-o1 + o2)) + 60*Power(b1,3)*
//                             ((35*Power(b2,2)*Power(o1 - o2,3))/2. - 4*Power(b2,3)*Power(o1 - o2,5) +
//                              21*b2*(-o1 + o2)) - 60*b1*((13*Power(b2,4)*Power(o1 - o2,3))/2. +
//                                                         7*Power(b2,3)*(-o1 + o2)))
//
//         ,
//
//
//         (-420*Power(b2,3) - 270*Power(b2,4)*Power(o1,2) +
//          30*Power(b1,2)*b2*(-42 + 51*b2*Power(o1 - o2,2) -
//                             20*Power(b2,2)*Power(o1 - o2,4)) -
//          30*b1*Power(b2,2)*(42 - 27*b2*Power(o1 - o2,2) - 4*Power(b2,2)*Power(o1 - o2,4)) -
//          60*Power(b1,4)*(15 - 10*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4))*
//          Power(o1 - o2,2) + 24*Power(b1,5)*(-5 + b2*Power(o1 - o2,2))*Power(o1 - o2,4) +
//          540*Power(b2,4)*o1*o2 - 270*Power(b2,4)*Power(o2,2) -
//          5*Power(b1,3)*(84 + 90*b2*Power(o1,2) - 16*Power(b2,3)*Power(o1 - o2,6) -
//                         180*b2*o1*o2 + 90*b2*Power(o2,2)))/32.
//
//         ,
//
//
//         Complex(0,-0.03125)*(-5*Power(b2,3)*(63 + 4*b2*Power(o1 - o2,2)) +
//                              35*b1*Power(b2,2)*(-3 + 8*b2*Power(o1 - o2,2)) -
//                              5*Power(b1,2)*b2*(-147 + 60*b2*Power(o1 - o2,2) +
//                                                8*Power(b2,2)*Power(o1 - o2,4)) +
//                              5*Power(b1,3)*(105 - 80*b2*Power(o1 - o2,2) + 24*Power(b2,2)*Power(o1 - o2,4)) -
//                              20*Power(b1,4)*(-10 + 3*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                              4*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2)
//
//         ,
//
//
//         (5*(3*Power(b2,2)*(7 + 3*b2*Power(o1 - o2,2)) -
//             2*b1*b2*(-21 + 18*b2*Power(o1 - o2,2) + Power(b2,2)*Power(o1 - o2,4)) +
//             3*Power(b1,2)*(7 - 5*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4)) -
//             6*Power(b1,3)*(-5 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//             2*Power(b1,4)*Power(o1 - o2,4)))/32.
//
//         ,
//         Complex(0,0.015625)*(-(Power(b2,2)*(63 + 2*b2*Power(o1 - o2,2))) -
//                              15*Power(b1,2)*(-7 + 4*b2*Power(o1 - o2,2)) +
//                              6*b1*b2*(7 + 5*b2*Power(o1 - o2,2)) + 20*Power(b1,3)*Power(o1 - o2,2))*(o1 - o2)
//
//         ,
//         (-(b2*(14 + 3*b2*Power(o1 - o2,2))) + b1*(-14 + 15*b2*Power(o1 - o2,2)) -
//          10*Power(b1,2)*Power(o1 - o2,2))/64.
//
//         ,
//
//         Complex(0,-0.0078125)*(5*b1 - 3*b2)*(o1 - o2)
//
//         ,
//         0.00390625)
        else if ( l1 == 3 && l2 == 5 ){
            switch ( lambda1 ) {
                case 0:
                    va = (105*Power(b2,4) + 90*Power(b2,5)*Power(o1,2) -
                          60*b1*Power(b2,3)*(-7 + 3*b2*Power(o1 - o2,2) + Power(b2,2)*Power(o1 - o2,4)) +
                          30*Power(b1,2)*Power(b2,2)*(21 - 26*b2*Power(o1 - o2,2) +
                                                      8*Power(b2,2)*Power(o1 - o2,4)) -
                          20*Power(b1,3)*b2*(-21 + 18*b2*Power(o1 - o2,2) -
                                             15*Power(b2,2)*Power(o1 - o2,4) + 4*Power(b2,3)*Power(o1 - o2,6)) +
                          5*Power(b1,4)*(21 + 90*b2*Power(o1 - o2,2) - 60*Power(b2,2)*Power(o1 - o2,4) +
                                         8*Power(b2,3)*Power(o1 - o2,6)) -
                          4*Power(b1,5)*(-75 + 60*b2*Power(o1 - o2,2) - 24*Power(b2,2)*Power(o1 - o2,4) +
                                         4*Power(b2,3)*Power(o1 - o2,6))*Power(o1 - o2,2) -
                          12*Power(b1,6)*(-5 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,4) -
                          180*Power(b2,5)*o1*o2 + 90*Power(b2,5)*Power(o2,2))/16.;
                    break;
                    
                case 1:
                    va = Complex(0,-0.0625)*(315*Power(b2,4)*(o1 - o2) +
                                             5*Power(b1,4)*(-105 + 60*b2*Power(o1 - o2,2) - 36*Power(b2,2)*Power(o1 - o2,4) +
                                                            8*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) + 30*Power(b2,5)*Power(o1 - o2,3) -
                                             12*Power(b1,5)*(25 - 14*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4))*
                                             Power(o1 - o2,3) - 12*Power(b1,6)*Power(o1 - o2,5) +
                                             30*Power(b1,2)*Power(b2,2)*(b2*Power(o1 - o2,3) + 4*Power(b2,2)*Power(o1 - o2,5) +
                                                                         21*(-o1 + o2)) + 60*Power(b1,3)*
                                             ((35*Power(b2,2)*Power(o1 - o2,3))/2. - 4*Power(b2,3)*Power(o1 - o2,5) +
                                              21*b2*(-o1 + o2)) - 60*b1*((13*Power(b2,4)*Power(o1 - o2,3))/2. +
                                                                         7*Power(b2,3)*(-o1 + o2)));
                    break;
                case 2:
                    va =(-420*Power(b2,3) - 270*Power(b2,4)*Power(o1,2) +
                         30*Power(b1,2)*b2*(-42 + 51*b2*Power(o1 - o2,2) -
                                            20*Power(b2,2)*Power(o1 - o2,4)) -
                         30*b1*Power(b2,2)*(42 - 27*b2*Power(o1 - o2,2) - 4*Power(b2,2)*Power(o1 - o2,4)) -
                         60*Power(b1,4)*(15 - 10*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4))*
                         Power(o1 - o2,2) + 24*Power(b1,5)*(-5 + b2*Power(o1 - o2,2))*Power(o1 - o2,4) +
                         540*Power(b2,4)*o1*o2 - 270*Power(b2,4)*Power(o2,2) -
                         5*Power(b1,3)*(84 + 90*b2*Power(o1,2) - 16*Power(b2,3)*Power(o1 - o2,6) -
                                        180*b2*o1*o2 + 90*b2*Power(o2,2)))/32.;
                    break;
                    
                case 3:
                    va =         Complex(0,-0.03125)*(-5*Power(b2,3)*(63 + 4*b2*Power(o1 - o2,2)) +
                                                      35*b1*Power(b2,2)*(-3 + 8*b2*Power(o1 - o2,2)) -
                                                      5*Power(b1,2)*b2*(-147 + 60*b2*Power(o1 - o2,2) +
                                                                        8*Power(b2,2)*Power(o1 - o2,4)) +
                                                      5*Power(b1,3)*(105 - 80*b2*Power(o1 - o2,2) + 24*Power(b2,2)*Power(o1 - o2,4)) -
                                                      20*Power(b1,4)*(-10 + 3*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
                                                      4*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2);
;
                    break;
                    
                case 4:
                    va =  (5*(3*Power(b2,2)*(7 + 3*b2*Power(o1 - o2,2)) -
                              2*b1*b2*(-21 + 18*b2*Power(o1 - o2,2) + Power(b2,2)*Power(o1 - o2,4)) +
                              3*Power(b1,2)*(7 - 5*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4)) -
                              6*Power(b1,3)*(-5 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
                              2*Power(b1,4)*Power(o1 - o2,4)))/32.;
                    break;
                    
                case 5:
                    va =         Complex(0,0.015625)*(-(Power(b2,2)*(63 + 2*b2*Power(o1 - o2,2))) -
                                                      15*Power(b1,2)*(-7 + 4*b2*Power(o1 - o2,2)) +
                                                      6*b1*b2*(7 + 5*b2*Power(o1 - o2,2)) + 20*Power(b1,3)*Power(o1 - o2,2))*(o1 - o2);
                    break;
                case 6:
                    va =  (-(b2*(14 + 3*b2*Power(o1 - o2,2))) + b1*(-14 + 15*b2*Power(o1 - o2,2)) -
                           10*Power(b1,2)*Power(o1 - o2,2))/64.;
                    break;
                    
                case 7:
                    va =Complex(0,-0.0078125)*(5*b1 - 3*b2)*(o1 - o2);
                    break;
                case 8:
                    va =0.00390625;
                    break;
            }
        }

//    List((
    
//          60*Power(b1,2)*Power(b2,3)*(-7 + 19*b2*Power(o1 - o2,2) -
//                                      6*Power(b2,2)*Power(o1 - o2,4))*(o1 - o2) +
//          15*b1*Power(b2,4)*(-77 + 36*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
//          (o1 - o2) - 40*Power(b1,4)*b2*
//          (-42 + 45*b2*Power(o1 - o2,2) - 18*Power(b2,2)*Power(o1 - o2,4) +
//           2*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) +
//          10*Power(b1,3)*Power(b2,2)*(147 - 72*b2*Power(o1 - o2,2) -
//                                      18*Power(b2,2)*Power(o1 - o2,4) + 8*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) +
//          Power(b1,5)*(525 - 300*b2*Power(o1 - o2,2) + 252*Power(b2,2)*Power(o1 - o2,4) -
//                       112*Power(b2,3)*Power(o1 - o2,6) + 16*Power(b2,4)*Power(o1 - o2,8))*(o1 - o2) +
//          12*Power(b1,6)*(25 - 18*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
//          Power(o1 - o2,3) + 12*Power(b1,7)*Power(o1 - o2,5) +
//          60*Power(b2,5)*(-2*b2*Power(o1 - o2,3) + 7*(-o1 + o2)))/16.
//
//         ,
//         Complex(0,0.03125)*(60*b1*Power(b2,3)*
//                             (63 - 7*b2*Power(o1 - o2,2) - 18*Power(b2,2)*Power(o1 - o2,4)) +
//                             15*Power(b2,4)*(63 + 84*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4)) +
//                             720*Power(b1,3)*b2*(5.25 - 7*b2*Power(o1 - o2,2) +
//                                                 5*Power(b2,2)*Power(o1 - o2,4) - Power(b2,3)*Power(o1 - o2,6)) +
//                             60*Power(b1,2)*Power(b2,2)*(94.5 - 112*b2*Power(o1 - o2,2) +
//                                                         21*Power(b2,2)*Power(o1 - o2,4) + 4*Power(b2,3)*Power(o1 - o2,6)) +
//                             5*Power(b1,4)*(189 + 420*b2*Power(o1 - o2,2) - 180*Power(b2,2)*Power(o1 - o2,4) -
//                                            48*Power(b2,3)*Power(o1 - o2,6) + 16*Power(b2,4)*Power(o1 - o2,8)) -
//                             4*Power(b1,5)*(-525 + 450*b2*Power(o1 - o2,2) - 156*Power(b2,2)*Power(o1 - o2,4) +
//                                            16*Power(b2,3)*Power(o1 - o2,6))*Power(o1 - o2,2) -
//                             12*Power(b1,6)*(-25 + 8*b2*Power(o1 - o2,2))*Power(o1 - o2,4))
//
//         ,
//         (15*b1*Power(b2,3)*(49*o1 - 33*b2*Power(o1 - o2,3) - 2*Power(b2,2)*Power(o1 - o2,5) -
//                             49*o2) + 30*Power(b2,4)*(14 + 3*b2*Power(o1 - o2,2))*(o1 - o2) +
//          5*Power(b1,4)*(-105 + 90*b2*Power(o1 - o2,2) - 48*Power(b2,2)*Power(o1 - o2,4) +
//                         8*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) -
//          3*Power(b1,5)*(75 - 38*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
//          Power(o1 - o2,3) - 6*Power(b1,6)*Power(o1 - o2,5) +
//          15*Power(b1,2)*Power(b2,2)*(-21*o1 - 24*b2*Power(o1 - o2,3) +
//                                      14*Power(b2,2)*Power(o1 - o2,5) + 21*o2) +
//          5*Power(b1,3)*b2*(-231*o1 + 180*b2*Power(o1 - o2,3) -
//                            24*Power(b2,2)*Power(o1 - o2,5) - 4*Power(b2,3)*Power(o1 - o2,7) + 231*o2))/8.
         
//         ,
//         Complex(0,0.0625)*(-5*Power(b2,3)*
//                            (63 + 63*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) +
//                            5*b1*Power(b2,2)*(-189 + 84*b2*Power(o1 - o2,2) +
//                                              38*Power(b2,2)*Power(o1 - o2,4)) -
//                            5*Power(b1,2)*b2*(189 - 252*b2*Power(o1 - o2,2) +
//                                              80*Power(b2,2)*Power(o1 - o2,4) + 4*Power(b2,3)*Power(o1 - o2,6)) +
//                            5*Power(b1,3)*(-63 - 40*Power(b2,2)*Power(o1 - o2,4) +
//                                           16*Power(b2,3)*Power(o1 - o2,6)) -
//                            5*Power(b1,4)*(105 - 70*b2*Power(o1 - o2,2) + 12*Power(b2,2)*Power(o1 - o2,4))*
//                            Power(o1 - o2,2) + 2*Power(b1,5)*(-25 + 4*b2*Power(o1 - o2,2))*Power(o1 - o2,4))
         
//         ,
//         ((-60*Power(b2,3)*(7 + b2*Power(o1 - o2,2)) +
//           5*b1*Power(b2,2)*(-63 + 78*b2*Power(o1 - o2,2) +
//                             2*Power(b2,2)*Power(o1 - o2,4)) +
//           15*Power(b1,3)*(35 - 30*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4)) -
//           10*Power(b1,2)*b2*(-63 + 15*b2*Power(o1 - o2,2) +
//                              8*Power(b2,2)*Power(o1 - o2,4)) -
//           10*Power(b1,4)*(-15 + 4*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//           2*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2))/32.
//
//         ,
//         Complex(0,0.015625)*(Power(b2,2)*
//                              (189 + 126*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) -
//                              2*b1*b2*(-189 + 147*b2*Power(o1 - o2,2) + 20*Power(b2,2)*Power(o1 - o2,4)) +
//                              3*Power(b1,2)*(63 - 70*b2*Power(o1 - o2,2) + 40*Power(b2,2)*Power(o1 - o2,4)) -
//                              10*Power(b1,3)*(-21 + 8*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                              10*Power(b1,4)*Power(o1 - o2,4))
//
//         ,
//         -((-2*Power(b2,2)*(14 + b2*Power(o1 - o2,2)) -
//            5*Power(b1,2)*(-7 + 4*b2*Power(o1 - o2,2)) +
//            b1*b2*(7 + 15*b2*Power(o1 - o2,2)) + 5*Power(b1,3)*Power(o1 - o2,2))*(o1 - o2))/
//         32.,
    
//    Complex(0,-0.015625)*(b1*(9 - 10*b2*Power(o1 - o2,2)) +
//                                   3*b2*(3 + b2*Power(o1 - o2,2)) + 5*Power(b1,2)*Power(o1 - o2,2))
//
//         ,
//         ((5*b1 - 4*b2)*(o1 - o2))/256.
//
//         ,
//
//         Complex(0,0.001953125))
//
    
        else if ( l1 == 4 && l2 == 5 ){
            switch ( lambda1 ) {
                case 0:
                    va = (60*Power(b1,2)*Power(b2,3)*(-7 + 19*b2*Power(o1 - o2,2) -
                                                6*Power(b2,2)*Power(o1 - o2,4))*(o1 - o2) +
                    15*b1*Power(b2,4)*(-77 + 36*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
                    (o1 - o2) - 40*Power(b1,4)*b2*
                    (-42 + 45*b2*Power(o1 - o2,2) - 18*Power(b2,2)*Power(o1 - o2,4) +
                     2*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) +
                    10*Power(b1,3)*Power(b2,2)*(147 - 72*b2*Power(o1 - o2,2) -
                                                18*Power(b2,2)*Power(o1 - o2,4) + 8*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) +
                    Power(b1,5)*(525 - 300*b2*Power(o1 - o2,2) + 252*Power(b2,2)*Power(o1 - o2,4) -
                                 112*Power(b2,3)*Power(o1 - o2,6) + 16*Power(b2,4)*Power(o1 - o2,8))*(o1 - o2) +
                    12*Power(b1,6)*(25 - 18*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
                    Power(o1 - o2,3) + 12*Power(b1,7)*Power(o1 - o2,5) +
                    60*Power(b2,5)*(-2*b2*Power(o1 - o2,3) + 7*(-o1 + o2)))/16.;
                    break;
                    
                case 1:
                    va =   Complex(0,0.03125)*(60*b1*Power(b2,3)*
                                               (63 - 7*b2*Power(o1 - o2,2) - 18*Power(b2,2)*Power(o1 - o2,4)) +
                                               15*Power(b2,4)*(63 + 84*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4)) +
                                               720*Power(b1,3)*b2*(5.25 - 7*b2*Power(o1 - o2,2) +
                                                                   5*Power(b2,2)*Power(o1 - o2,4) - Power(b2,3)*Power(o1 - o2,6)) +
                                               60*Power(b1,2)*Power(b2,2)*(94.5 - 112*b2*Power(o1 - o2,2) +
                                                                           21*Power(b2,2)*Power(o1 - o2,4) + 4*Power(b2,3)*Power(o1 - o2,6)) +
                                               5*Power(b1,4)*(189 + 420*b2*Power(o1 - o2,2) - 180*Power(b2,2)*Power(o1 - o2,4) -
                                                              48*Power(b2,3)*Power(o1 - o2,6) + 16*Power(b2,4)*Power(o1 - o2,8)) -
                                               4*Power(b1,5)*(-525 + 450*b2*Power(o1 - o2,2) - 156*Power(b2,2)*Power(o1 - o2,4) +
                                                              16*Power(b2,3)*Power(o1 - o2,6))*Power(o1 - o2,2) -
                                               12*Power(b1,6)*(-25 + 8*b2*Power(o1 - o2,2))*Power(o1 - o2,4));
                    break;
                case 2:
                    va =(15*b1*Power(b2,3)*(49*o1 - 33*b2*Power(o1 - o2,3) - 2*Power(b2,2)*Power(o1 - o2,5) -
                                            49*o2) + 30*Power(b2,4)*(14 + 3*b2*Power(o1 - o2,2))*(o1 - o2) +
                         5*Power(b1,4)*(-105 + 90*b2*Power(o1 - o2,2) - 48*Power(b2,2)*Power(o1 - o2,4) +
                                        8*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) -
                         3*Power(b1,5)*(75 - 38*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
                         Power(o1 - o2,3) - 6*Power(b1,6)*Power(o1 - o2,5) +
                         15*Power(b1,2)*Power(b2,2)*(-21*o1 - 24*b2*Power(o1 - o2,3) +
                                                     14*Power(b2,2)*Power(o1 - o2,5) + 21*o2) +
                         5*Power(b1,3)*b2*(-231*o1 + 180*b2*Power(o1 - o2,3) -
                                           24*Power(b2,2)*Power(o1 - o2,5) - 4*Power(b2,3)*Power(o1 - o2,7) + 231*o2))/8.;
                    break;
                
                    
                case 3:
                    va =         Complex(0,0.0625)*(-5*Power(b2,3)*
                                                    (63 + 63*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) +
                                                    5*b1*Power(b2,2)*(-189 + 84*b2*Power(o1 - o2,2) +
                                                                      38*Power(b2,2)*Power(o1 - o2,4)) -
                                                    5*Power(b1,2)*b2*(189 - 252*b2*Power(o1 - o2,2) +
                                                                      80*Power(b2,2)*Power(o1 - o2,4) + 4*Power(b2,3)*Power(o1 - o2,6)) +
                                                    5*Power(b1,3)*(-63 - 40*Power(b2,2)*Power(o1 - o2,4) +
                                                                   16*Power(b2,3)*Power(o1 - o2,6)) -
                                                    5*Power(b1,4)*(105 - 70*b2*Power(o1 - o2,2) + 12*Power(b2,2)*Power(o1 - o2,4))*
                                                    Power(o1 - o2,2) + 2*Power(b1,5)*(-25 + 4*b2*Power(o1 - o2,2))*Power(o1 - o2,4));
                    break;
                    
                case 4:
                    va =  ((-60*Power(b2,3)*(7 + b2*Power(o1 - o2,2)) +
                                   5*b1*Power(b2,2)*(-63 + 78*b2*Power(o1 - o2,2) +
                                                     2*Power(b2,2)*Power(o1 - o2,4)) +
                                   15*Power(b1,3)*(35 - 30*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4)) -
                                   10*Power(b1,2)*b2*(-63 + 15*b2*Power(o1 - o2,2) +
                                                      8*Power(b2,2)*Power(o1 - o2,4)) -
                                   10*Power(b1,4)*(-15 + 4*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
                                   2*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2))/32.;
                    break;
                case 5:
                    va =   Complex(0,0.015625)*(Power(b2,2)*
                                               (189 + 126*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) -
                                               2*b1*b2*(-189 + 147*b2*Power(o1 - o2,2) + 20*Power(b2,2)*Power(o1 - o2,4)) +
                                               3*Power(b1,2)*(63 - 70*b2*Power(o1 - o2,2) + 40*Power(b2,2)*Power(o1 - o2,4)) -
                                               10*Power(b1,3)*(-21 + 8*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
                                               10*Power(b1,4)*Power(o1 - o2,4));
                    break;
                    
                case 6:
                    va = -((-2*Power(b2,2)*(14 + b2*Power(o1 - o2,2)) -
                           5*Power(b1,2)*(-7 + 4*b2*Power(o1 - o2,2)) +
                           b1*b2*(7 + 15*b2*Power(o1 - o2,2)) + 5*Power(b1,3)*Power(o1 - o2,2))*(o1 - o2))/
                    32.;
                    
                case 7:
                    va = Complex(0,-0.015625)*(b1*(9 - 10*b2*Power(o1 - o2,2)) +
                                              3*b2*(3 + b2*Power(o1 - o2,2)) + 5*Power(b1,2)*Power(o1 - o2,2));
                    break;
                case 8:
                    va =((5*b1 - 4*b2)*(o1 - o2))/256.;
                    break;
                case 9:
                    va =Complex(0,0.001953125);
                    break;
            }
        }
    
    
//    {(1/32)*(-60*(2*b2*(o1 - o2)*(o1 - o2)- 5)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*b1*b1*b1*b1*b1*b1*b1 -
//             20*(8*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 48*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                 105*b2*(o1 - o2)*(o1 - o2) - 105)*(o1 - o2)*(o1 - o2)*b1*b1*b1*b1*b1*b1 +
//             (-32*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 240*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//              120*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 2100*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//              3150*b2*(o1 - o2)*(o1 - o2) + 945)*b1*b1*b1*b1*b1 +
//             15*b2*(16*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 160*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                    420*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 420*b2*(o1 - o2)*(o1 - o2) + 315)*b1*b1*b1*b1 -
//             10*b2*b2*(16*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 12*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                      630*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 1470*b2*(o1 - o2)*(o1 - o2) - 945)*
//             b1*b1*b1 +
//             30*b2*b2*b2*(32*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 70*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                      210*b2*(o1 - o2)*(o1 - o2) + 315)*b1*b1 -
//             15*b2*b2*b2*b2*(8*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 140*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                      210*b2*(o1 - o2)*(o1 - o2) - 315)*b1 +
//             300*b2*b2*b2*b2*b2*(b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 7*b2*(o1 - o2)*(o1 - o2) + 63/20)),
//        (1/32)*-5*I*(-12*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*b1*b1*b1*b1*b1*b1*b1 -
//                     12*(4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 23*b2*(o1 - o2)*(o1 - o2) + 35)*(o1 - o2)*(o1 - o2)*(o1 - o2)*
//                     b1*b1*b1*b1*b1*b1 - (16*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 192*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                             612*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 840*b2*(o1 - o2)*(o1 - o2) + 945)*(o1 -
//                                                                               o2)*
//                     b1*b1*b1*b1*b1 + b2*(16*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 900*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                                2940*b2*(o1 - o2)*(o1 - o2) - 2835)*(o1 - o2)*b1*b1*b1*b1 +
//                     6*(-32*b2*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 150*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                        315*b2*b2*(o1 - o2))*b1*b1*b1 -
//                     6*b2*b2*b2*(-8*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                             102*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 490*b2*(o1 - o2)*(o1 - o2) - 315)*(o1 -
//                                                                               o2)*
//                     b1*b1 -
//                     3*b2*b2*b2*b2*(92*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 280*b2*(o1 - o2)*(o1 - o2) - 945)*
//                     (o1 - o2)*b1 + 3*b2*b2*b2*b2*b2*(4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                                            140*b2*(o1 - o2)*(o1 - o2)*(o1 - o2) + 315*o1 - 315*o2)),
//        (5/64)*(12*(4*b2*(o1 - o2)*(o1 - o2) - 15)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*b1*b1*b1*b1*b1*b1 +
//                16*(2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 27*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                    90*b2*(o1 - o2)*(o1 - o2) - 105)*(o1 - o2)*(o1 - o2)*b1*b1*b1*b1*b1 +
//                5*(-16*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 96*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                   36*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 168*b2*o1*o1 - 168*b2*o2*o2 +
//                   336*b2*o1*o2 - 189)*b1*b1*b1*b1 - 4*b2*(-8*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                                                    120*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 900*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                                                    1470*b2*o1*o1 -
//                                                    1470*b2*o2*o2 + 2940*b2*o1*o2 + 945)*b1*b1*b1 +
//                6*b2*(-72*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 30*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                        980*b2*o1*o1 + 980*b2*o2*o2 - 1960*b2*o1*o2 - 945)*b1*b1 -
//                12*b2*b2*b2*(-4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 120*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                         70*b2*o1*o1 + 70*b2*o2*o2 - 140*b2*o1*o2 + 315)*b1 +
//                15*b2*b2*b2*b2*(-12*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 112*b2*(o1 - o2)*(o1 - o2) - 63)),
//        (1/16)*-5*I*(2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*b1*b1*b1*b1*b1*b1 +
//                     (4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 48*b2*(o1 - o2)*(o1 - o2) + 105)*(o1 - o2)*
//                     b1*b1*b1*b1*b1 - 5*(4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 30*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                               63*b2*(o1 - o2)*(o1 - o2) - 63)*b1*b1*b1*b1 +
//                     10*b2*(2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 42*b2*(o1 - o2)*(o1 - o2) + 63)*b1*b1*b1 -
//                     2*b2*b2*b2*(2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 75*b2*(o1 - o2)*(o1 - o2) - 210)*
//                     (o1 - o2)*(o1 - o2)*b1*b1 + 3*b2*b2*b2*(16*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                                                105*b2*(o1 - o2)*(o1 - o2) - 210)*b1 -
//                     b2*b2*b2*b2*(2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 105*b2*(o1 - o2)*(o1 - o2) + 315))*
//        (o1 - o2), (1/64)*-5*(2*(2*b2*(o1 - o2)*(o1 - o2) - 15)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*
//                              b1*b1*b1*b1*b1 - 10*(4*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 27*b2*(o1 - o2)*(o1 - o2) + 42)*
//                              (o1 - o2)*(o1 - o2)*b1*b1*b1*b1 + 5*(16*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                                                    60*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 42*b2*(o1 - o2)*(o1 - o2) - 63)*b1*b1*b1 -
//                              5*b2*(8*b2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 60*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                                    252*b2*(o1 - o2)*(o1 - o2) + 189)*b1*b1 +
//                              b2*b2*(4*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 270*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) +
//                                    210*b2*(o1 - o2)*(o1 - o2) - 945)*b1 -
//                              15*b2*b2*b2*(2*b2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 28*b2*(o1 - o2)*(o1 - o2) + 21)),
//        (1/64)*I*(b1 - b2)*(2*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                            6*b1*b1*b1*(8*b2*(o1 - o2)*(o1 - o2) - 35)*(o1 - o2)*(o1 - o2) +
//                            b2*b2*(2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 210*b2*(o1 - o2)*(o1 - o2) + 945) -
//                            6*b1*b2*(8*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 105*b2*(o1 - o2)*(o1 - o2) - 315) +
//                            b1*b1*(152*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 630*b2*(o1 - o2)*(o1 - o2) + 945))*
//        (o1 - o2), (1/128)*-5*(2*b1*b1*b1*b1*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) -
//                               4*b1*b1*b1*(5*b2*(o1 - o2)*(o1 - o2) - 14)*(o1 - o2)*(o1 - o2) +
//                               b2*b2*(2*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) + 56*b2*(o1 - o2)*(o1 - o2) + 63) -
//                               2*b1*b2*(10*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2)+ 42*b2*(o1 - o2)*(o1 - o2) - 63) +
//                               b1*b1*(40*b2*b2*(o1 - o2)*(o1 - o2)*(o1 - o2)*(o1 - o2) - 84*b2*(o1 - o2)*(o1 - o2) + 63)),
//        (1/64)*-5*I*(b1 - b2)*(b1*b1*(o1 - o2)*(o1 - o2) +
//                               b2*(b2*(o1 - o2)*(o1 - o2) + 9) + b1*(9 - 4*b2*(o1 - o2)*(o1 - o2)))*(o1 -
//                                                                                     o2),
//        (5/512)*(4*b1*b1*(o1 - o2)*(o1 - o2) + b2*(4*b2*(o1 - o2)*(o1 - o2) + 9) +
//                 b1*(9 - 10*b2*(o1 - o2)*(o1 - o2))), (5/512)*I*(b1 - b2)*(o1 - o2),
//        -(1/1024)}
    
        else if ( l1 == 5 && l2 == 5 ){
            switch ( lambda1 ) {
                case 0:
                    va = (300*Power(b2,5)*(3.15 + 7*b2*Power(o1 - o2,2) + Power(b2,2)*Power(o1 - o2,4)) -
                          15*b1*Power(b2,4)*(-315 - 210*b2*Power(o1 - o2,2) +
                                             140*Power(b2,2)*Power(o1 - o2,4) + 8*Power(b2,3)*Power(o1 - o2,6)) +
                          30*Power(b1,2)*Power(b2,3)*(315 - 210*b2*Power(o1 - o2,2) -
                                                      70*Power(b2,2)*Power(o1 - o2,4) + 32*Power(b2,3)*Power(o1 - o2,6)) +
                          15*Power(b1,4)*b2*(315 - 420*b2*Power(o1 - o2,2) +
                                             420*Power(b2,2)*Power(o1 - o2,4) - 160*Power(b2,3)*Power(o1 - o2,6) +
                                             16*Power(b2,4)*Power(o1 - o2,8)) -
                          10*Power(b1,3)*Power(b2,2)*(-945 + 1470*b2*Power(o1 - o2,2) -
                                                      630*Power(b2,2)*Power(o1 - o2,4) + 12*Power(b2,3)*Power(o1 - o2,6) +
                                                      16*Power(b2,4)*Power(o1 - o2,8)) +
                          Power(b1,5)*(945 + 3150*b2*Power(o1 - o2,2) - 2100*Power(b2,2)*Power(o1 - o2,4) -
                                       120*Power(b2,3)*Power(o1 - o2,6) + 240*Power(b2,4)*Power(o1 - o2,8) -
                                       32*Power(b2,5)*Power(o1 - o2,10)) -
                          20*Power(b1,6)*(-105 + 105*b2*Power(o1 - o2,2) - 48*Power(b2,2)*Power(o1 - o2,4) +
                                          8*Power(b2,3)*Power(o1 - o2,6))*Power(o1 - o2,2) -
                          60*Power(b1,7)*(-5 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,4))/32.;
                    break;
                    
                case 1:
                    va = Complex(0,-0.15625)*(6*Power(b1,3)*(-315*Power(b2,2)*(o1 - o2) +
                                                             150*Power(b2,4)*Power(o1 - o2,5) - 32*Power(b2,5)*Power(o1 - o2,7)) +
                                              3*Power(b2,5)*(315*o1 + 140*b2*Power(o1 - o2,3) + 4*Power(b2,2)*Power(o1 - o2,5) -
                                                             315*o2) - 3*b1*Power(b2,4)*(-945 + 280*b2*Power(o1 - o2,2) +
                                                                                         92*Power(b2,2)*Power(o1 - o2,4))*(o1 - o2) -
                                              6*Power(b1,2)*Power(b2,3)*(-315 + 490*b2*Power(o1 - o2,2) -
                                                                         102*Power(b2,2)*Power(o1 - o2,4) - 8*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) +
                                              Power(b1,4)*b2*(-2835 + 2940*b2*Power(o1 - o2,2) -
                                                              900*Power(b2,2)*Power(o1 - o2,4) + 16*Power(b2,4)*Power(o1 - o2,8))*(o1 - o2) -
                                              Power(b1,5)*(945 - 840*b2*Power(o1 - o2,2) + 612*Power(b2,2)*Power(o1 - o2,4) -
                                                           192*Power(b2,3)*Power(o1 - o2,6) + 16*Power(b2,4)*Power(o1 - o2,8))*(o1 - o2) -
                                              12*Power(b1,6)*(35 - 23*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
                                              Power(o1 - o2,3) - 12*Power(b1,7)*Power(o1 - o2,5));
                    break;
                case 2:
                    va =(5*(15*Power(b2,4)*(-63 - 112*b2*Power(o1 - o2,2) - 12*Power(b2,2)*Power(o1 - o2,4)) +
                            16*Power(b1,5)*(-105 + 90*b2*Power(o1 - o2,2) - 27*Power(b2,2)*Power(o1 - o2,4) +
                                            2*Power(b2,3)*Power(o1 - o2,6))*Power(o1 - o2,2) +
                            12*Power(b1,6)*(-15 + 4*b2*Power(o1 - o2,2))*Power(o1 - o2,4) -
                            4*Power(b1,3)*b2*(945 - 1470*b2*Power(o1,2) + 900*Power(b2,2)*Power(o1 - o2,4) -
                                              120*Power(b2,3)*Power(o1 - o2,6) - 8*Power(b2,4)*Power(o1 - o2,8) +
                                              2940*b2*o1*o2 - 1470*b2*Power(o2,2)) +
                            5*Power(b1,4)*(-189 - 168*b2*Power(o1,2) - 36*Power(b2,2)*Power(o1 - o2,4) +
                                           96*Power(b2,3)*Power(o1 - o2,6) - 16*Power(b2,4)*Power(o1 - o2,8) +
                                           336*b2*o1*o2 - 168*b2*Power(o2,2)) -
                            12*b1*Power(b2,3)*(315 + 70*b2*Power(o1,2) - 120*Power(b2,2)*Power(o1 - o2,4) -
                                               4*Power(b2,3)*Power(o1 - o2,6) - 140*b2*o1*o2 + 70*b2*Power(o2,2)) +
                            6*Power(b1,2)*Power(b2,2)*(-945 + 980*b2*Power(o1,2) -
                                                       30*Power(b2,2)*Power(o1 - o2,4) - 72*Power(b2,3)*Power(o1 - o2,6) -
                                                       1960*b2*o1*o2 + 980*b2*Power(o2,2))))/64.;
                    break;
                    
                case 3:
                    va =Complex(0,-0.3125)*(-(Power(b2,4)*(315 + 105*b2*Power(o1 - o2,2) +
                                                           2*Power(b2,2)*Power(o1 - o2,4))) +
                                            3*b1*Power(b2,3)*(-210 + 105*b2*Power(o1 - o2,2) +
                                                              16*Power(b2,2)*Power(o1 - o2,4)) +
                                            10*Power(b1,3)*b2*(63 - 42*b2*Power(o1 - o2,2) + 2*Power(b2,3)*Power(o1 - o2,6)) -
                                            5*Power(b1,4)*(-63 + 63*b2*Power(o1 - o2,2) - 30*Power(b2,2)*Power(o1 - o2,4) +
                                                           4*Power(b2,3)*Power(o1 - o2,6)) -
                                            2*Power(b1,2)*Power(b2,3)*(-210 + 75*b2*Power(o1 - o2,2) +
                                                                       2*Power(b2,2)*Power(o1 - o2,4))*Power(o1 - o2,2) +
                                            Power(b1,5)*(105 - 48*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
                                            Power(o1 - o2,2) + 2*Power(b1,6)*Power(o1 - o2,4))*(o1 - o2);
                    break;
                    
                case 4:
                    va =(-5*(-15*Power(b2,3)*(21 + 28*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) +
                             b1*Power(b2,2)*(-945 + 210*b2*Power(o1 - o2,2) +
                                             270*Power(b2,2)*Power(o1 - o2,4) + 4*Power(b2,3)*Power(o1 - o2,6)) -
                             5*Power(b1,2)*b2*(189 - 252*b2*Power(o1 - o2,2) +
                                               60*Power(b2,2)*Power(o1 - o2,4) + 8*Power(b2,3)*Power(o1 - o2,6)) +
                             5*Power(b1,3)*(-63 + 42*b2*Power(o1 - o2,2) - 60*Power(b2,2)*Power(o1 - o2,4) +
                                            16*Power(b2,3)*Power(o1 - o2,6)) -
                             10*Power(b1,4)*(42 - 27*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
                             Power(o1 - o2,2) + 2*Power(b1,5)*(-15 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,4)))
                    /64.;
                    break;
                    
                case 5:
                    va =Complex(0,0.015625)*(b1 - b2)*(Power(b2,2)*
                                                       (945 + 210*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) -
                                                       6*b1*b2*(-315 + 105*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4)) +
                                                       Power(b1,2)*(945 - 630*b2*Power(o1 - o2,2) + 152*Power(b2,2)*Power(o1 - o2,4)) -
                                                       6*Power(b1,3)*(-35 + 8*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
                                                       2*Power(b1,4)*Power(o1 - o2,4))*(o1 - o2);
                    break;
                case 6:
                    va =(-5*(Power(b2,2)*(63 + 56*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) -
                             2*b1*b2*(-63 + 42*b2*Power(o1 - o2,2) + 10*Power(b2,2)*Power(o1 - o2,4)) +
                             Power(b1,2)*(63 - 84*b2*Power(o1 - o2,2) + 40*Power(b2,2)*Power(o1 - o2,4)) -
                             4*Power(b1,3)*(-14 + 5*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
                             2*Power(b1,4)*Power(o1 - o2,4)))/128.;
                    break;
                    
                case 7:
                    va = Complex(0,-0.078125)*(b1 - b2)*(b1*(9 - 4*b2*Power(o1 - o2,2)) +
                                                         b2*(9 + b2*Power(o1 - o2,2)) + Power(b1,2)*Power(o1 - o2,2))*(o1 - o2);
                    break;
                case 8:
                    va =(5*(b1*(9 - 10*b2*Power(o1 - o2,2)) + b2*(9 + 4*b2*Power(o1 - o2,2)) +
                            4*Power(b1,2)*Power(o1 - o2,2)))/512.;
                    break;
                case 9:
                    va =Complex(0,0.009765625)*(b1 - b2)*(o1 - o2);
                    break;
                    
                case 10:
                    va = -0.0009765625;
                    break;
            }
       }
//            else if ( l1 == 0 && l2 == 6 ){
//            switch ( lambda1 ) {
//                case 0:
//                    va = (45*b1*Power(b2,2) + 15*Power(b2,3) + 45*Power(b1,2)*b2*(1 + 2*b2*Power(o1 - o2,2)) +
//                          15*Power(b1,3)*(1 + 12*b2*Power(o1 - o2,2)) +
//                          30*Power(b1,4)*(3 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                          60*Power(b1,5)*Power(o1 - o2,4) + 8*Power(b1,6)*Power(o1 - o2,6))/8.;
//                    break;
//
//                case 1:
//                    va = Complex(0,0.75)*b1*(30*b1*b2 + 15*Power(b2,2) +
//                                             5*Power(b1,2)*(3 + 4*b2*Power(o1 - o2,2)) + 20*Power(b1,3)*Power(o1 - o2,2) +
//                                             4*Power(b1,4)*Power(o1 - o2,4))*(o1 - o2);
//                    break;
//                case 2:
//                    va =(-15*(6*b1*b2 + 3*Power(b2,2) + 3*Power(b1,2)*(1 + 4*b2*Power(o1 - o2,2)) +
//                              12*Power(b1,3)*Power(o1 - o2,2) + 4*Power(b1,4)*Power(o1 - o2,4)))/16.;
//                    break;
//
//                case 3:
//                    va =Complex(0,-1.25)*b1*(3*b1 + 3*b2 + 2*Power(b1,2)*Power(o1 - o2,2))*(o1 - o2);
//                    break;
//
//                case 4:
//                    va =(15*(b1 + b2 + 2*Power(b1,2)*Power(o1 - o2,2)))/32.;
//                    break;
//
//                case 5:
//                    va =Complex(0,0.1875)*b1*(o1 - o2);
//                    break;
//                case 6:
//                    va =-0.015625;
//                    break;
//            }
//        }
//
//            else if ( l1 == 1 && l2 == 6 ){
//                switch ( lambda1 ) {
//                    case 0:
//                        va = -((-45*b1*Power(b2,3) + 15*Power(b2,4) +
//                                45*Power(b1,2)*Power(b2,2)*(-5 + 2*b2*Power(o1 - o2,2)) +
//                                15*Power(b1,3)*b2*(-17 + 4*b2*Power(o1 - o2,2)) +
//                                30*Power(b1,4)*(-3 - 5*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) +
//                                12*Power(b1,5)*(-10 + 3*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                                8*Power(b1,6)*(-3 + b2*Power(o1 - o2,2))*Power(o1 - o2,4))*(o1 - o2))/8.;
//                        break;
//
//                    case 1:
//                        va = Complex(0,0.0625)*(105*Power(b2,3) + 45*Power(b1,2)*b2*(7 + 2*b2*Power(o1 - o2,2)) -
//                                                45*b1*Power(b2,2)*(-7 + 4*b2*Power(o1 - o2,2)) -
//                                                15*Power(b1,3)*(-7 - 48*b2*Power(o1 - o2,2) + 16*Power(b2,2)*Power(o1 - o2,4)) -
//                                                30*Power(b1,4)*(-15 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,2) -
//                                                12*Power(b1,5)*(-15 + 4*b2*Power(o1 - o2,2))*Power(o1 - o2,4) +
//                                                8*Power(b1,6)*Power(o1 - o2,6));
//                        break;
//                    case 2:
//                        va =(-3*(60*b1*Power(b2,2) - 15*Power(b2,3) + 10*Power(b1,3)*(9 + 2*b2*Power(o1 - o2,2)) -
//                                 15*Power(b1,2)*b2*(-11 + 4*b2*Power(o1 - o2,2)) -
//                                 20*Power(b1,4)*(-4 + b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                                 8*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2))/16.;
//                        break;
//
//                    case 3:
//                        va =Complex(0,-0.15625)*(21*Power(b2,2) - 6*b1*b2*(-7 + 4*b2*Power(o1 - o2,2)) +
//                                                 3*Power(b1,2)*(7 + 12*b2*Power(o1 - o2,2)) -
//                                                 4*Power(b1,3)*(-15 + 4*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                                                 12*Power(b1,4)*Power(o1 - o2,4));
//                        break;
//
//                    case 4:
//                        va =(5*(15*b1*b2 - 3*Power(b2,2) - 6*Power(b1,2)*(-3 + b2*Power(o1 - o2,2)) +
//                                8*Power(b1,3)*Power(o1 - o2,2))*(o1 - o2))/32.;
//                        break;
//
//                    case 5:
//                        va =Complex(0,0.046875)*(7*b2 + b1*(7 - 4*b2*Power(o1 - o2,2)) +
//                                                 10*Power(b1,2)*Power(o1 - o2,2));
//                        break;
//                    case 6:
//                        va =-((6*b1 - b2)*(o1 - o2))/64.;
//                        break;
//
//                    case 7:
//                        va =Complex(0,-0.0078125);
//                        break;
//
//                }
//            }
//
//
//            else if ( l1 == 2 && l2 == 6 ){
//                switch ( lambda1 ) {
//                    case 0:
//                        va = (420*b1*Power(b2,3) + 105*Power(b2,4) - 270*b1*Power(b2,4)*Power(o1,2) +
//                              60*Power(b1,3)*b2*(7 + 5*b2*Power(o1 - o2,2) - 2*Power(b2,2)*Power(o1 - o2,4)) +
//                              90*Power(b1,2)*Power(b2,2)*(7 - 6*b2*Power(o1 - o2,2) +
//                                                          2*Power(b2,2)*Power(o1 - o2,4)) +
//                              15*Power(b1,4)*(7 + 66*b2*Power(o1 - o2,2) - 40*Power(b2,2)*Power(o1 - o2,4) +
//                                              8*Power(b2,3)*Power(o1 - o2,6)) + 30*Power(b2,5)*Power(o1 - o2,2) +
//                              6*Power(b1,5)*(75 - 20*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
//                              Power(o1 - o2,2) + 4*Power(b1,6)*
//                              (45 - 22*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*Power(o1 - o2,4) +
//                              8*Power(b1,7)*Power(o1 - o2,6) + 540*b1*Power(b2,4)*o1*o2 -
//                              270*b1*Power(b2,4)*Power(o2,2))/16.;
//                        break;
//
//                    case 1:
//                        va = Complex(0,-0.125)*(105*Power(b2,4)*(o1 - o2) +
//                                                15*Power(b1,4)*(-21 - 10*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
//                                                (o1 - o2) - 90*b1*Power(b2,4)*Power(o1 - o2,3) -
//                                                12*Power(b1,5)*(25 - 12*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4))*
//                                                Power(o1 - o2,3) + 4*Power(b1,6)*(-9 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,5) +
//                                                90*Power(b1,2)*Power(b2,2)*(3*b2*Power(o1 - o2,3) + 7*(-o1 + o2)) +
//                                                120*Power(b1,3)*((17*Power(b2,2)*Power(o1 - o2,3))/4. -
//                                                                 Power(b2,3)*Power(o1 - o2,5) + 7*b2*(-o1 + o2)));
//                        break;
//                    case 2:
//                        va =(-15*Power(b2,3)*(14 + 3*b2*Power(o1 - o2,2)) +
//                             90*b1*Power(b2,2)*(-7 + 5*b2*Power(o1 - o2,2)) -
//                             90*Power(b1,2)*b2*(7 - 4*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) +
//                             30*Power(b1,3)*(-7 - 27*b2*Power(o1 - o2,2) + 10*Power(b2,2)*Power(o1 - o2,4)) -
//                             15*Power(b1,4)*(45 - 20*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
//                             Power(o1 - o2,2) + 12*Power(b1,5)*(-15 + 4*b2*Power(o1 - o2,2))*Power(o1 - o2,4) -
//                             4*Power(b1,6)*Power(o1 - o2,6))/16.;
//                        break;
//
//                    case 3:
//                        va =Complex(0,-0.0625)*(-105*Power(b2,3) + 15*b1*Power(b2,2)*(7 + 4*b2*Power(o1 - o2,2)) -
//                                                15*Power(b1,2)*b2*(-35 + 16*b2*Power(o1 - o2,2)) +
//                                                5*Power(b1,3)*(63 - 20*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4)) -
//                                                20*Power(b1,4)*(-10 + 3*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                                                12*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2);
//                        break;
//
//                    case 4:
//                        va =(5*(3*b1*b2*(14 - 11*b2*Power(o1 - o2,2)) + 3*Power(b2,2)*(7 + b2*Power(o1 - o2,2)) +
//                                3*Power(b1,2)*(7 + 3*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) -
//                                Power(b1,3)*(-45 + 16*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                                6*Power(b1,4)*Power(o1 - o2,4)))/32.;
//                        break;
//
//                    case 5:
//                        va =Complex(0,0.03125)*(-21*Power(b2,2) + Power(b1,2)*(63 - 30*b2*Power(o1 - o2,2)) +
//                                                6*b1*b2*(7 + b2*Power(o1 - o2,2)) + 20*Power(b1,3)*Power(o1 - o2,2))*(o1 - o2);
//                        break;
//                    case 6:
//                        va =(-(b2*(14 + b2*Power(o1 - o2,2))) + 2*b1*(-7 + 6*b2*Power(o1 - o2,2)) -
//                             15*Power(b1,2)*Power(o1 - o2,2))/64.;
//                        break;
//
//                    case 7:
//                        va =Complex(0,-0.015625)*(3*b1 - b2)*(o1 - o2);
//                        break;
//                    case 8:
//                        va =0.00390625;
//                        break;
//                }
//            }
//            else if ( l1 == 3 && l2 == 6 ){
//                switch ( lambda1 ) {
//                    case 0:
//                        va =(-315*Power(b2,5)*o1 + 2520*Power(b1,3)*Power(b2,2)*(o1 - o2) +
//                             630*Power(b1,2)*Power(b2,3)*(o1 - o2) - 630*b1*Power(b2,4)*(o1 - o2) +
//                             6*Power(b1,5)*(105 + 75*b2*Power(o1 - o2,2) - 48*Power(b2,2)*Power(o1 - o2,4) +
//                                            4*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) -
//                             15*Power(b1,4)*b2*(-147 + 114*b2*Power(o1 - o2,2) -
//                                                48*Power(b2,2)*Power(o1 - o2,4) + 8*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) -
//                             1860*Power(b1,3)*Power(b2,3)*Power(o1 - o2,3) +
//                             180*Power(b1,2)*Power(b2,4)*Power(o1 - o2,3) +
//                             450*b1*Power(b2,5)*Power(o1 - o2,3) - 30*Power(b2,6)*Power(o1 - o2,3) -
//                             4*Power(b1,6)*(-150 + 99*b2*Power(o1 - o2,2) - 30*Power(b2,2)*Power(o1 - o2,4) +
//                                            4*Power(b2,3)*Power(o1 - o2,6))*Power(o1 - o2,3) +
//                             360*Power(b1,3)*Power(b2,4)*Power(o1 - o2,5) -
//                             180*Power(b1,2)*Power(b2,5)*Power(o1 - o2,5) -
//                             24*Power(b1,7)*(-3 + b2*Power(o1 - o2,2))*Power(o1 - o2,5) + 315*Power(b2,5)*o2)/16. ;
//                        break;
//
//                    case 1:
//                        va = Complex(0,0.09375)*(315*Power(b2,4) -
//                                                 30*b1*Power(b2,3)*(-42 + 21*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4)) +
//                                                 60*Power(b1,2)*Power(b2,2)*(31.5 - 35*b2*Power(o1 - o2,2) +
//                                                                             11*Power(b2,2)*Power(o1 - o2,4)) +
//                                                 20*Power(b1,3)*b2*(63 - 21*b2*Power(o1 - o2,2) + 24*Power(b2,2)*Power(o1 - o2,4) -
//                                                                    8*Power(b2,3)*Power(o1 - o2,6)) +
//                                                 5*Power(b1,4)*(63 + 378*b2*Power(o1 - o2,2) - 240*Power(b2,2)*Power(o1 - o2,4) +
//                                                                40*Power(b2,3)*Power(o1 - o2,6)) + 210*Power(b2,5)*Power(o1 - o2,2) -
//                                                 2*Power(b1,5)*(-525 + 300*b2*Power(o1 - o2,2) - 108*Power(b2,2)*Power(o1 - o2,4) +
//                                                                16*Power(b2,3)*Power(o1 - o2,6))*Power(o1 - o2,2) +
//                                                 4*Power(b1,6)*(75 - 34*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
//                                                 Power(o1 - o2,4) + 8*Power(b1,7)*Power(o1 - o2,6));
//                        break;
//                    case 2:
//                        va =(3*(210*Power(b2,4)*o1 + 15*Power(b2,5)*Power(o1,3) -
//                                30*Power(b1,2)*Power(b2,2)*(21*o1 - 5*b2*Power(o1 - o2,3) -
//                                                            2*Power(b2,2)*Power(o1 - o2,5) - 21*o2) -
//                                30*b1*Power(b2,3)*(-7 + 8*b2*Power(o1 - o2,2))*(o1 - o2) +
//                                5*Power(b1,4)*(-84 + 15*b2*Power(o1 - o2,2) - 12*Power(b2,2)*Power(o1 - o2,4) +
//                                               4*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) -
//                                12*Power(b1,5)*(25 - 13*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4))*
//                                Power(o1 - o2,3) + 4*Power(b1,6)*(-6 + b2*Power(o1 - o2,2))*Power(o1 - o2,5) -
//                                210*Power(b2,4)*o2 - 45*Power(b2,5)*Power(o1,2)*o2 +
//                                45*Power(b2,5)*o1*Power(o2,2) - 15*Power(b2,5)*Power(o2,3) +
//                                30*Power(b1,3)*b2*(-35*o1 + 26*b2*Power(o1 - o2,3) -
//                                                   6*Power(b2,2)*Power(o1 - o2,5) + 35*o2)))/16.;
//                        break;
//
//                    case 3:
//                        va =Complex(0,-0.03125)*(315*Power(b2,3)*(2 + b2*Power(o1 - o2,2)) -
//                                                 30*b1*Power(b2,2)*(-63 + 42*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4)) +
//                                                 30*Power(b1,2)*b2*(63 - 63*b2*Power(o1 - o2,2) + 26*Power(b2,2)*Power(o1 - o2,4)) -
//                                                 10*Power(b1,3)*(-63 - 126*b2*Power(o1 - o2,2) + 30*Power(b2,2)*Power(o1 - o2,4) +
//                                                                 8*Power(b2,3)*Power(o1 - o2,6)) +
//                                                 45*Power(b1,4)*(35 - 20*b2*Power(o1 - o2,2) + 4*Power(b2,2)*Power(o1 - o2,4))*
//                                                 Power(o1 - o2,2) - 12*Power(b1,5)*(-25 + 6*b2*Power(o1 - o2,2))*Power(o1 - o2,4) +
//                                                 4*Power(b1,6)*Power(o1 - o2,6));
//                        break;
//
//                    case 4:
//                        va =(3*(-5*Power(b2,3)*(21 + b2*Power(o1 - o2,2)) -
//                                5*Power(b1,2)*b2*(-63 + 27*b2*Power(o1 - o2,2) +
//                                                  2*Power(b2,2)*Power(o1 - o2,4)) +
//                                5*Power(b1,3)*(42 - 25*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4)) +
//                                85*b1*Power(b2,3)*Power(o1 - o2,2) -
//                                10*Power(b1,4)*(-10 + 3*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                                4*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2))/32.;
//                        break;
//
//                    case 5:
//                        va =Complex(0,0.046875)*(21*Power(b2,2)*(3 + b2*Power(o1 - o2,2)) +
//                                                 b1*b2*(126 - 105*b2*Power(o1 - o2,2) - 4*Power(b2,2)*Power(o1 - o2,4)) +
//                                                 3*Power(b1,2)*(21 - 7*b2*Power(o1 - o2,2) + 10*Power(b2,2)*Power(o1 - o2,4)) -
//                                                 5*Power(b1,3)*(-21 + 8*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                                                 10*Power(b1,4)*Power(o1 - o2,4));
//                        break;
//                    case 6:
//                        va =-((Power(b1,2)*(84 - 45*b2*Power(o1 - o2,2)) - Power(b2,2)*(42 + b2*Power(o1 - o2,2)) +
//                               6*b1*b2*(7 + 3*b2*Power(o1 - o2,2)) + 20*Power(b1,3)*Power(o1 - o2,2))*(o1 - o2))
//                        /64.;
//                        break;
//
//                    case 7:
//                        va =Complex(0,-0.0234375)*(-6*b1*(-1 + b2*Power(o1 - o2,2)) + b2*(6 + b2*Power(o1 - o2,2)) +
//                                                   5*Power(b1,2)*Power(o1 - o2,2));
//                        break;
//                    case 8:
//                        va =(3*(2*b1 - b2)*(o1 - o2))/256.;
//                        break;
//                    case 9:
//                        va =Complex(0,0.001953125);
//                        break;
//
//                }
//            }
//            else if ( l1 == 4 && l2 == 6 ){
//                switch ( lambda1 ) {
//                    case 0:
//                        va = (945*Power(b2,5) + 1260*b1*Power(b2,4)*(3.75 - Power(b2,2)*Power(o1 - o2,4)) +
//                              30*Power(b1,3)*Power(b2,2)*(315 - 420*b2*Power(o1 - o2,2) +
//                                                          238*Power(b2,2)*Power(o1 - o2,4) - 40*Power(b2,3)*Power(o1 - o2,6)) +
//                              90*Power(b1,2)*Power(b2,3)*(105 - 105*b2*Power(o1 - o2,2) +
//                                                          14*Power(b2,2)*Power(o1 - o2,4) + 4*Power(b2,3)*Power(o1 - o2,6)) +
//                              3*Power(b1,5)*(315 + 2520*b2*Power(o1 - o2,2) - 2100*Power(b2,2)*Power(o1 - o2,4) +
//                                             608*Power(b2,3)*Power(o1 - o2,6) - 48*Power(b2,4)*Power(o1 - o2,8)) +
//                              15*Power(b1,4)*b2*(315 + 84*Power(b2,2)*Power(o1 - o2,4) -
//                                                 88*Power(b2,3)*Power(o1 - o2,6) + 16*Power(b2,4)*Power(o1 - o2,8)) +
//                              1260*Power(b2,6)*Power(o1 - o2,2) +
//                              2*Power(b1,6)*(1575 - 1050*b2*Power(o1 - o2,2) + 516*Power(b2,2)*Power(o1 - o2,4) -
//                                             144*Power(b2,3)*Power(o1 - o2,6) + 16*Power(b2,4)*Power(o1 - o2,8))*
//                              Power(o1 - o2,2) + 60*Power(b2,7)*Power(o1 - o2,4) +
//                              12*Power(b1,7)*(75 - 44*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4))*
//                              Power(o1 - o2,4) + 24*Power(b1,8)*Power(o1 - o2,6))/32.;
//                        break;
//
//                    case 1:
//                        va = Complex(0,-0.0625)*(1890*Power(b2,5)*(o1 - o2) +
//                                                 30*Power(b1,3)*Power(b2,2)*(-315 + 196*b2*Power(o1 - o2,2) -
//                                                                             6*Power(b2,2)*Power(o1 - o2,4) - 8*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) +
//                                                 30*Power(b1,4)*b2*(-315 + 294*b2*Power(o1 - o2,2) -
//                                                                    120*Power(b2,2)*Power(o1 - o2,4) + 16*Power(b2,3)*Power(o1 - o2,6))*(o1 - o2) -
//                                                 3*Power(b1,5)*(945 + 60*Power(b2,2)*Power(o1 - o2,4) -
//                                                                96*Power(b2,3)*Power(o1 - o2,6) + 16*Power(b2,4)*Power(o1 - o2,8))*(o1 - o2) +
//                                                 420*Power(b2,6)*Power(o1 - o2,3) -
//                                                 1440*Power(b1,2)*Power(b2,4)*(2.625 - b2*Power(o1 - o2,2))*Power(o1 - o2,3) +
//                                                 4*Power(b1,6)*(-525 + 360*b2*Power(o1 - o2,2) - 96*Power(b2,2)*Power(o1 - o2,4) +
//                                                                8*Power(b2,3)*Power(o1 - o2,6))*Power(o1 - o2,3) +
//                                                 12*Power(b1,7)*(-15 + 4*b2*Power(o1 - o2,2))*Power(o1 - o2,5) -
//                                                 45*b1*Power(b2,4)*(-105*o1 + 56*b2*Power(o1 - o2,3) +
//                                                                    4*Power(b2,2)*Power(o1 - o2,5) + 105*o2));
//                        break;
//                    case 2:
//                        va =(3*(-1575*Power(b2,4) - 1680*Power(b2,5)*Power(o1,2) +
//                                60*b1*Power(b2,3)*(-105 + 28*b2*Power(o1 - o2,2) +
//                                                   22*Power(b2,2)*Power(o1 - o2,4)) +
//                                8*Power(b1,5)*(-525 + 375*b2*Power(o1 - o2,2) -
//                                               132*Power(b2,2)*Power(o1 - o2,4) + 16*Power(b2,3)*Power(o1 - o2,6))*
//                                Power(o1 - o2,2) - 60*Power(b2,6)*Power(o1 - o2,4) -
//                                4*Power(b1,6)*(225 - 92*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4))*
//                                Power(o1 - o2,4) - 16*Power(b1,7)*Power(o1 - o2,6) + 3360*Power(b2,5)*o1*o2 -
//                                1680*Power(b2,5)*Power(o2,2) +
//                                5*Power(b1,4)*(-315 - 1176*b2*Power(o1,2) + 660*Power(b2,2)*Power(o1 - o2,4) -
//                                               32*Power(b2,3)*Power(o1 - o2,6) - 16*Power(b2,4)*Power(o1 - o2,8) +
//                                               2352*b2*o1*o2 - 1176*b2*Power(o2,2)) +
//                                20*Power(b1,3)*b2*(-315 + 294*b2*Power(o1,2) - 228*Power(b2,2)*Power(o1 - o2,4) +
//                                                   52*Power(b2,3)*Power(o1 - o2,6) - 588*b2*o1*o2 + 294*b2*Power(o2,2)) +
//                                30*Power(b1,2)*Power(b2,2)*(-315 + 364*b2*Power(o1,2) -
//                                                            86*Power(b2,2)*Power(o1 - o2,4) - 8*Power(b2,3)*Power(o1 - o2,6) -
//                                                            728*b2*o1*o2 + 364*b2*Power(o2,2))))/64.;
//                        break;
//
//                    case 3:
//                        va =Complex(0,-0.125)*(15*Power(b1,4)*(63*(o1 - o2) - 35*b2*Power(o1 - o2,3) +
//                                                               20*Power(b2,2)*Power(o1 - o2,5) - 4*Power(b2,3)*Power(o1 - o2,7)) +
//                                               20*Power(b1,3)*((441*b2*(o1 - o2))/4. - 84*Power(b2,2)*Power(o1 - o2,3) +
//                                                               15*Power(b2,3)*Power(o1 - o2,5) + Power(b2,4)*Power(o1 - o2,7)) -
//                                               630*Power(b2,4)*(o1 - o2) - 15*Power(b1,2)*Power(b2,2)*
//                                               (-63 - 14*b2*Power(o1 - o2,2) + 18*Power(b2,2)*Power(o1 - o2,4))*(o1 - o2) -
//                                               105*Power(b2,5)*Power(o1 - o2,3) +
//                                               3*Power(b1,5)*(175 - 90*b2*Power(o1 - o2,2) + 12*Power(b2,2)*Power(o1 - o2,4))*
//                                               Power(o1 - o2,3) - 2*Power(b1,6)*(-15 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,5) -
//                                               15*b1*Power(b2,3)*(63*o1 - 49*b2*Power(o1,3) - 2*Power(b2,2)*Power(o1 - o2,5) -
//                                                                  63*o2 + 147*b2*Power(o1,2)*o2 - 147*b2*o1*Power(o2,2) + 49*b2*Power(o2,3)));
//                        break;
//
//                    case 4:
//                        va =(15*Power(b2,3)*(105 + 84*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) -
//                             15*b1*Power(b2,2)*(-315 + 168*b2*Power(o1 - o2,2) +
//                                                46*Power(b2,2)*Power(o1 - o2,4)) +
//                             15*Power(b1,2)*b2*(315 - 378*b2*Power(o1 - o2,2) +
//                                                132*Power(b2,2)*Power(o1 - o2,4) + 4*Power(b2,3)*Power(o1 - o2,6)) -
//                             5*Power(b1,3)*(-315 - 252*b2*Power(o1 - o2,2) - 60*Power(b2,2)*Power(o1 - o2,4) +
//                                            64*Power(b2,3)*Power(o1 - o2,6)) +
//                             30*Power(b1,4)*(105 - 65*b2*Power(o1 - o2,2) + 12*Power(b2,2)*Power(o1 - o2,4))*
//                             Power(o1 - o2,2) - 6*Power(b1,5)*(-75 + 16*b2*Power(o1 - o2,2))*Power(o1 - o2,4) +
//                             4*Power(b1,6)*Power(o1 - o2,6))/64.;
//                        break;
//
//                    case 5:
//                        va =Complex(0,0.09375)*(-14*Power(b2,3)*(9 + b2*Power(o1 - o2,2)) +
//                                                b1*Power(b2,2)*(-63 + 112*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4)) -
//                                                4*Power(b1,2)*b2*(-63 + 21*b2*Power(o1 - o2,2) + 5*Power(b2,2)*Power(o1 - o2,4)) +
//                                                Power(b1,3)*(189 - 140*b2*Power(o1 - o2,2) + 40*Power(b2,2)*Power(o1 - o2,4)) -
//                                                10*Power(b1,4)*(-7 + 2*b2*Power(o1 - o2,2))*Power(o1 - o2,2) +
//                                                2*Power(b1,5)*Power(o1 - o2,4))*(o1 - o2);
//                        break;
//                    case 6:
//                        va =(-(Power(b2,2)*(315 + 168*b2*Power(o1 - o2,2) + 2*Power(b2,2)*Power(o1 - o2,4))) +
//                             6*b1*b2*(-105 + 84*b2*Power(o1 - o2,2) + 8*Power(b2,2)*Power(o1 - o2,4)) -
//                             9*Power(b1,2)*(35 - 28*b2*Power(o1 - o2,2) + 20*Power(b2,2)*Power(o1 - o2,4)) +
//                             20*Power(b1,3)*(-21 + 8*b2*Power(o1 - o2,2))*Power(o1 - o2,2) -
//                             30*Power(b1,4)*Power(o1 - o2,4))/128.;
//                        break;
//
//                    case 7:
//                        va =Complex(0,-0.03125)*(9*b1*b2*(1 + b2*Power(o1 - o2,2)) -
//                                                 Power(b2,2)*(18 + b2*Power(o1 - o2,2)) -
//                                                 3*Power(b1,2)*(-9 + 5*b2*Power(o1 - o2,2)) + 5*Power(b1,3)*Power(o1 - o2,2))*
//                        (o1 - o2);
//                        break;
//                    case 8:
//                        va =(3*(b1*(15 - 16*b2*Power(o1 - o2,2)) + b2*(15 + 4*b2*Power(o1 - o2,2)) +
//                                10*Power(b1,2)*Power(o1 - o2,2)))/512.;
//                        break;
//                    case 9:
//                        va =Complex(0,0.00390625)*(3*b1 - 2*b2)*(o1 - o2);
//                        break;
//
//                    case 10:
//                        va =-0.0009765625;
//                        break;
//                }
//            }
//            else if ( l1 == 5 && l2 == 6 ){
//                switch ( lambda1 ) {
//                    case 0:
//                        va = ;
//                        break;
//
//                    case 1:
//                        va = ;
//                        break;
//                    case 2:
//                        va =;
//                        break;
//
//                    case 3:
//                        va =;
//                        break;
//
//                    case 4:
//                        va =;
//                        break;
//
//                    case 5:
//                        va =;
//                        break;
//                    case 6:
//                        va =;
//                        break;
//
//                    case 7:
//                        va =;
//                        break;
//                    case 8:
//                        va =;
//                        break;
//                    case 9:
//                        va =;
//                        break;
//
//                    case 10:
//                        va =;
//                        break;
//                    case 11:
//                        va =;
//                        break;
//                }
//            }
//            else if ( l1 == 6 && l2 == 6 ){
//                switch ( lambda1 ) {
//                    case 0:
//                        va = ;
//                        break;
//
//                    case 1:
//                        va = ;
//                        break;
//                    case 2:
//                        va =;
//                        break;
//
//                    case 3:
//                        va =;
//                        break;
//
//                    case 4:
//                        va =;
//                        break;
//
//                    case 5:
//                        va =;
//                        break;
//                    case 6:
//                        va =;
//                        break;
//
//                    case 7:
//                        va =;
//                        break;
//                    case 8:
//                        va =;
//                        break;
//                    case 9:
//                        va =;
//                        break;
//
//                    case 10:
//                        va =;
//                        break;
//
//                    case 11:
//                        va =;
//                        break;
//                    case 12:
//                        va =;
//                        break;
//                }
//            }
//
//            else if ( l1 == 0 && l2 == 7 ){
//                switch ( lambda1 ) {
//                    case 0:
//                        va = ;
//                        break;
//
//                    case 1:
//                        va = ;
//                        break;
//                    case 2:
//                        va =;
//                        break;
//
//                    case 3:
//                        va =;
//                        break;
//
//                    case 4:
//                        va =;
//                        break;
//
//                    case 5:
//                        va =;
//                        break;
//                    case 6:
//                        va =;
//                        break;
//
//                    case 7:
//                        va =;
//                        break;
//                    case 8:
//                        va =;
//                        break;
//                    case 9:
//                        va =;
//                        break;
//
//                    case 10:
//                        va =;
//                        break;
//
//                    case 11:
//                        va =;
//                        break;
//                    case 12:
//                        va =;
//                        break;
//                }
//            }
    return va/pow(b1+b2,l1+l2);
}

DCOMPLEX gauGetPoly( double b1,INT_TYPE l1, double o1, INT_TYPE lambda1 ){
    DCOMPLEX va = 0.;
    if ( l1 == 0 ){
        switch ( lambda1 ) {
            case 0:
                va =  1;
                break;
        }
    }
    else if ( l1 == 1 ){
        switch ( lambda1 ) {
            case 1:
                va =  I/2.;
                break;
                
        }
    }    else if ( l1 == 2 ){
        switch ( lambda1 ) {
            case 0:
                va =  1/2.;
                break;
                
            case 2:
                va = -1/4. ;
                break;
                
        }
    }else if ( l1 == 3 ){
        switch ( lambda1 ) {
            case 1:
                va = 3*I/4.    ;
                break;
                
            case 3:
                va =  -I/8.     ;
                
                break;
        }
    }else if ( l1 == 4 ){
        switch ( lambda1 ) {
            case 0:
                va = 3/4.      ;
                break;
                
            case 2:
                va =  -3/4.      ;
                break;
            case 4:
                va =  1/16.    ;
                break;
                
        }
    }else if ( l1 == 5 ){
        switch ( lambda1 ) {
            case 1:
                va = I*15/8.   ;
                break;
                
            case 3:
                va =  -I*5/8.    ;
                break;
            case 5:
                va =  I/32.    ;
                break;
                
        }
    }else if ( l1 == 6 ){
        switch ( lambda1 ) {
            case 0:
                va = 15*1/8.     ;
                break;
                
            case 2:
                va =  -1*45/16.    ;
                break;
            case 4:
                va =  1*15/32.  ;
                break;
            case 6:
                va =  -1/64.     ;
                break;
                
        }
    }else if ( l1 == 7 ){
        switch ( lambda1 ) {
            case 1:
                va = I*105/16. ;
                break;
                
            case 3:
                va =  -I*105/32. ;
                break;
            case 5:
                va =  I*21/64. ;
                break;
            case 7:
                va =  -I/128.   ;
                break;
                
        }
    }else if ( l1 == 8 ){
        switch ( lambda1 ) {
            case 0:
                va = 105/16.         ;
                break;
                
            case 2:
                va = - 105/8.          ;
                break;
            case 4:
                va =  105/32.        ;
                break;
            case 6:
                va =  -7/32.          ;
                break;
            case 8:
                va =  1/256.           ;
                break;
                
        }
    }
    
    
    return va*pow(b1,-l1-0.5+(l1-lambda1)/2.)/sqrt(2.);
}


DCOMPLEX aaGGCGG( double invbeta, struct general_2index * pa){
    struct basisElement b1,  k1, b2,  k2;
    
    if ( pa->i[0].bra.index < pa->i[0].ket.index){
        b1 = pa->i[0].bra;
        k1 = pa->i[0].ket;
    }else {
        k1 = pa->i[0].bra;
        b1 = pa->i[0].ket;

    }
    
    if ( pa->i[1].bra.index < pa->i[1].ket.index){
        b2 = pa->i[1].bra;
        k2 = pa->i[1].ket;
    }else {
        k2 = pa->i[1].bra;
        b2 = pa->i[1].ket;
    }
    
    INT_TYPE l,ll,l2,xl = b1.index + k1.index , xl2 = b2.index + k2.index;
    double gamma , delta;
    DCOMPLEX value=0.;
    DCOMPLEX hg[1+xl+xl2],p1[1+xl],p2[1+xl2];
    gamma = sqr(invbeta/2.);
    gamma += aaGetGamma(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin);
    gamma += aaGetGamma(b2.length, b2.index, b2.origin, k2.length, k2.index, k2.origin);
    delta = 0.;
    delta += aaGetDelta(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin);
    delta -= aaGetDelta(b2.length, b2.index, b2.origin, k2.length, k2.index, k2.origin);

    for ( ll = 0; ll <= xl+xl2; ll++){
        hg[ll] =  hyperGeometric(sqrt(gamma), ll, delta);
    }
    for ( ll = 0; ll <= xl; ll++)
        p1[ll] = aaGetPoly(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin,ll);
    
    for ( ll = 0; ll <= xl2; ll++)
        p2[ll] = sign(ll)*aaGetPoly(b2.length, b2.index, b2.origin, k2.length, k2.index, k2.origin,ll);
    
    for ( l = 0; l <= xl ; l++)
        for ( l2 = 0; l2 <= xl2 ; l2++){
            value += p1[l]* p2[l2]* hg[l+l2];
        }
//    for ( ll = 0; ll <= xl+xl2; ll++){
//       printf("%d %f i%f , %f i%f , %f i%f\n", ll,xg[ll] ,hg[ll], hg[ll]*xg[ll]);
//    }
    return value*aaGetConst(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin)*aaGetConst(b2.length, b2.index, b2.origin, k2.length, k2.index, k2.origin);//invbeta
}

DCOMPLEX aaGdnGdm( INT_TYPE n,INT_TYPE m,  struct general_index * pa){
    struct basisElement b1,  k1;
    
    b1 = pa->bra;
    k1 = pa->ket;
    
    INT_TYPE i,l,ll,ll2,l1 = b1.index ,l2 = k1.index;;
    double gamma , delta;
    DCOMPLEX in,im,value=0.;
    DCOMPLEX xv[1+l2+l1+n+m],hg[1+l2+n+m+l1],dn[1+n+l1],dm[1+m+l2];
    
    for ( ll = 0; ll <= l1+n+l2+m ; ll++)
        xv[ll] = 0.;
    for ( ll = 0; ll <= l1+n ; ll++)
        dn[ll] = 0.;
    for ( ll = 0; ll <= l2+m ; ll++)
        dm[ll] = 0.;

    gamma = 0.;
    gamma += 1/4.*(1/b1.length + 1/k1.length);
    delta = 0.;
    delta +=  b1.origin - k1.origin;
    
    for ( ll = 0; ll <= m+n+l1+l2; ll++){
        hg[ll] =  hyperGeometric(sqrt(gamma), ll, delta);
    }
    
    in = 1.;
    for ( i = 0 ; i < n ; i++)
        in *= I;
    
    im = 1.;
    for ( i = 0 ; i < m ; i++)
        im *= -I;

    for ( ll = 0; ll <= l1; ll++){
        dn[ll+n] = conj(gauGetPoly(b1.length, b1.index, b1.origin,ll));
    }
    for ( ll2 = 0; ll2 <= l2; ll2++){
        dm[ll2+m] = (gauGetPoly(k1.length, k1.index,k1.origin,ll2));
    }

    for ( ll = 0; ll <= l1+n ; ll++)
        for ( ll2= 0; ll2 <= l2+m ; ll2++){
            xv[ll+ll2] += dn[ll]*hg[ll+ll2]*dm[ll2];
        }
    for ( ll = 0; ll <= l1+n+l2+m ; ll++){
        value += xv[ll];
    }
    return in*im*value;
}

double GTOnorm ( struct basisElement ba ){
    struct general_index ga ;
    ga.bra = ba;
    ga.ket = ba;
 //   printf("..(%f %f) \n",aaGdnGdm(0, 0, &ga) );
    return sqrt ( 1./aaGdnGdm(0, 0, &ga));
    
}



INT_TYPE nCp ( INT_TYPE N , INT_TYPE P ){
    double va = 1;
    if ( N == 2 )
        switch ( P ){
            case 1 :
                va = 2;
                break;
                
        }
    if ( N == 3 )
        switch ( P ){
            case 1 :
                va = 3;
                break;
            case 2 :
                va = 3;
                break;
                
        }
    if ( N == 4 )
        switch ( P ){
            case 1 :
                va = 4;
                break;
            case 2 :
                va = 6;
                break;
            case 3 :
                va = 4;
                break;
        }
    if ( N == 5 )
        switch ( P ){
            case 1 :
                va = 5;
                break;
            case 2 :
                va = 10;
                break;
            case 3 :
                va = 10;
                break;
            case 4 :
                va = 5;
                break;
        }
    if ( N == 6 )
        switch ( P ){
            case 1 :
                va = 6;
                break;
            case 2 :
                va = 15;
                break;
            case 3 :
                va = 20;
                break;
            case 4 :
                va = 15;
                break;
            case 5 :
                va = 6;
                break;
        }
    if ( N == 7 )
        switch ( P ){
            case 1 :
                va = 7;
                break;
            case 2 :
                va = 21;
                break;
            case 3 :
                va = 35;
                break;
            case 4 :
                va = 35;
                break;
            case 5 :
                va = 21;
                break;
            case 6 :
                va = 7;
                break;
        }
    if ( N == 8 )
        switch ( P ){
            case 1 :
                va = 8;
                break;
            case 2 :
                va = 28;
                break;
            case 3 :
                va = 56;
                break;
            case 4 :
                va = 70;
                break;
            case 5 :
                va = 56;
                break;
            case 6 :
                va = 28;
                break;
            case 7 :
                va = 8;
                break;
        }
    if ( N == 9 )
        switch ( P ){
            case 1 :
                va = 9;
                break;
            case 2 :
                va = 36;
                break;
            case 3 :
                va = 84;
                break;
            case 4 :
                va = 126;
                break;
            case 5 :
                va = 126;
                break;
            case 6 :
                va = 84;
                break;
            case 7 :
                va = 36;
                break;
            case 8 :
                va = 9;
                break;

        }
    if ( N == 10 )
        switch ( P ){
            case 1 :
                va = 10;
                break;
            case 2 :
                va = 45;
                break;
            case 3 :
                va = 120;
                break;
            case 4 :
                va = 210;
                break;
            case 5 :
                va = 252;
                break;
            case 6 :
                va = 210;
                break;
            case 7 :
                va = 120;
                break;
            case 8 :
                va = 45;
                break;
            case 9 :
                va = 10;
                break;
        }
    else if ( P > 10 ){
        exit(0);
    }
    return va;
}

//DCOMPLEX aabGdnGdm( INT_TYPE n,INT_TYPE m, double boost, struct general_index * pa){
//    struct basisElement b1,  k1;
//
//    b1 = pa->bra;
//    k1 = pa->ket;
//
//    if ( boost != 0.)
//        exit(0);
//
//    INT_TYPE p,i,l,ll,lll,ll2,l1 = b1.index ,l2 = k1.index;;
//    double va,gamma , delta;
//    DCOMPLEX in,im,value=0.;
//    DCOMPLEX xv[1+l2+l1+n+m],hg[1+l2+n+m+l1],dn[1+n+l1],dm[1+m+l2];
//
//    for ( ll = 0; ll <= l1+n+l2+m ; ll++)
//        xv[ll] = 0.;
//    for ( ll = 0; ll <= l1+n ; ll++)
//        dn[ll] = 0.;
//    for ( ll = 0; ll <= l2+m ; ll++)
//        dm[ll] = 0.;
//
//    gamma = 0.;
//    gamma += 1/4.*(1/b1.length + 1/k1.length);
//    delta = 0.5*boost/b1.length;//Needs be complex!!! cannot do that!
//    delta +=  b1.origin - k1.origin;
//
//    for ( ll = 0; ll <= m+n+l1+l2; ll++){
//        hg[ll] =  hyperGeometric(sqrt(gamma), ll, delta);
//    }
//
//    in = 1.;
//    for ( i = 0 ; i < n ; i++)
//        in *= I;
//
//    im = 1.;
//    for ( i = 0 ; i < m ; i++)
//        im *= -I;
//
//    for ( ll = 0; ll <= l1; ll++){
//        double ga = gauGetPoly(b1.length, b1.index, b1.origin,ll);
//        //(k+boost)^(ll+n)
//        for ( lll = 0 ; lll <= ll + n ; lll++){
//            double ba = 1.;
//            for ( p = 0; p < ll+n-lll; p++)
//                ba *= boost ;
//            dn[lll] += ga * nCp( ll+n,lll) * ba;
//        }
//    }
//    for ( ll2 = 0; ll2 <= l2; ll2++){
//        dm[ll2+m] = gauGetPoly(k1.length, k1.index,k1.origin,ll2);
//    }
//
//    for ( ll = 0; ll <= l1+n ; ll++)
//        for ( ll2= 0; ll2 <= l2+m ; ll2++){
//            xv[ll+ll2] += dn[ll]*hg[ll+ll2]*dm[ll2];
//        }
//    for ( ll = 0; ll <= l1+n+l2+m ; ll++){
//        value += xv[ll];
//        printf("%d + %f\n", ll, xv[ll]);
//    }
//    return in*im*value* exp(-sqr(boost/2.)/b1.length);
//}



DCOMPLEX aaGGCD( double invbeta,   double position, struct general_index * pa){
    struct basisElement b1,  k1;
    
    
    if ( pa->bra.index < pa->ket.index){
        b1 = pa->bra;
        k1 = pa->ket;
    }else {
        k1 = pa->bra;
        b1 = pa->ket;
    }
    
    INT_TYPE l,ll,xl = b1.index + k1.index;
    double gamma , delta;
    DCOMPLEX value=0.;
    DCOMPLEX hg[1+xl],p1[1+xl];
    gamma = sqr(invbeta/2.);
    gamma += aaGetGamma(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin);
    delta =  - position;
    delta += aaGetDelta(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin);
    
    for ( ll = 0; ll <= xl; ll++){
        hg[ll] =  hyperGeometric(sqrt(gamma), ll, delta);
    }
    
    for ( ll = 0; ll <= xl; ll++){
        p1[ll] = aaGetPoly(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin,ll);
    }
    
    for ( l = 0; l <= xl ; l++)
    {
        value += p1[l]*hg[l];
    }
    return value*aaGetConst(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin);//invbeta
}


DCOMPLEX FGG( double k, struct general_index * pa){
    struct basisElement b1,  k1;
        if ( pa->bra.index < pa->ket.index){
            b1 = pa->bra;
            k1 = pa->ket;
        }else {
            k1 = pa->bra;
            b1 = pa->ket;
        }

    INT_TYPE l,ll,xl = b1.index + k1.index ;
    double gamma , delta;
    DCOMPLEX value=0.;
    DCOMPLEX hg[1+xl],p1[1+xl];
    gamma = 0;
    gamma += aaGetGamma(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin);
    delta =  0;
    delta += aaGetDelta(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin);
    
    double kp = 1;
    for ( ll = 0; ll <= xl; ll++){
        hg[ll] =  exp(- sqr(gamma * k) ) * kp * cexp(I * delta * k ) ;
        kp *= k;
    }
    
    for ( ll = 0; ll <= xl; ll++){
        p1[ll] = aaGetPoly(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin,ll);
    }
    
    for ( l = 0; l <= xl ; l++)
        {
            value += p1[l]*hg[l];
        }
    return value*aaGetConst(b1.length, b1.index, b1.origin, k1.length, k1.index, k1.origin);
}


double BoB (struct basisElement b1, struct basisElement b2 ){
    
    if ( b1.type == nullComponent || b2.type == nullComponent    )
    {
        printf("null\n");
        exit(0);
    }
    if ( b1.type != b2.type ){
        printf("compo");
        exit(0);
    }
    
    if ( b1.basis == SincBasisElement && b2.basis == SincBasisElement ){
        double va ;
        if ( b2.type <= 3 ){
           va = SS ( b1.length,b1.length*(b1.index)+ b1.origin, b2.length, b2.length*( b2.index ) + b2.origin);
        }
        else {//periodic
            va =  periodicSS ( b1.length,b1.length*(b1.index) + b1.origin, b1.grid, b2.length, b2.length*( b2.index ) + b2.origin,b2.grid);
        }
        
        return va;
    }else if ( b1.basis == GaussianBasisElement && b2.basis == GaussianBasisElement ){
        
        if ( b2.type <= 3 ){
            
            struct general_index pa ;
            pa.bra = b1;
            pa.ket = b2;
            
            return aaGdnGdm(0, 0, &pa);
            
        }else{
            printf("p-GG?");
            exit(0);
        }
    }else if ( b1.basis == GaussianBasisElement && b2.basis == SincBasisElement ){
        if ( b2.type <= 3 ){

            struct general_index pa ;
            pa.bra = b1;
            pa.ket = b2;



            return FGS(0,0,&pa);
        }else{
            printf("p-GS?");
            exit(0);
        }
    }else if ( b1.basis == SincBasisElement && b2.basis == GaussianBasisElement ){
        if ( b2.type <= 3 ){
            struct general_index pa ;
            pa.bra = b1;
            pa.ket = b2;
            return FGS(0,0,&pa);
        }else{
            printf("p-GS?");
            exit(0);
        }
    }

    printf("over rails\n");
    exit(0);
    return 0;
}

DCOMPLEX BdB (struct basisElement b1, struct basisElement b2){
    if ( b1.type == nullComponent || b2.type == nullComponent    )
    {
        printf("null\n");
        exit(0);
    }
    if ( b1.type != b2.type ){
        printf("compo");
        exit(0);
    }
    
    if ( b1.basis == SincBasisElement && b2.basis == SincBasisElement ){
        if ( b2.type <= 3 ){
            double arg = b1.length*b1.index + b1.origin - (b2.length*b2.index + b2.origin);
            if ( fabs(arg) < 1e-6 ){
                if ( b1.length > b2.length ){
                    arg /= b1.length;
                    return sqrt(b2.length/b1.length)/sqr(b1.length) * ( arg*(pi*pi/3. - pi * pi * pi * pi/30.*arg*arg + pi * pi * pi * pi*pi*pi/840.*arg*arg*arg*arg)  );
                }else {
                    arg /= b2.length;
                    
                    return sqrt(b1.length/b2.length)/sqr(b2.length) * ( arg*(pi*pi/3. - pi * pi * pi * pi/30.*arg*arg + pi * pi * pi * pi*pi*pi/840.*arg*arg*arg*arg) );
                }
            } else {
                if ( b1.length > b2.length ){
                    arg /= b1.length;
                    
                    return sqrt(b2.length/b1.length)/(b1.length) * ( ( Sinc(1,arg)- cos(pi*arg))/arg );
                }else {
                    arg /= b2.length;
                    
                    return sqrt(b1.length/b2.length)/(b2.length) * ( ( Sinc(1,arg)- cos(pi*arg))/arg );
                }
                
            }
        }
        else {
            
            double arg = b1.length*b1.index + b1.origin - (b2.length*b2.index + b2.origin);
            
            if ( b1.length > b2.length ){
                arg /= b1.length;
                
                return periodicSdS(arg,b1.grid) * sqrt(b2.length/b1.length);
            }else {
                arg /= b2.length;
                
                return periodicSdS(arg,b2.grid) * sqrt(b1.length/b2.length);
            }
        }
    }else if ( b1.basis == GaussianBasisElement && b2.basis == GaussianBasisElement ){
        if ( b2.type <= 3 ){
            
            struct general_index pa ;
            pa.bra = b1;
            pa.ket = b2;
            
            return aaGdnGdm(0, 1, &pa);
        }else{
            printf("p-GdG?");
            exit(0);
        }
        
        
    }else if ( b1.basis == GaussianBasisElement && b2.basis == SincBasisElement ){
        if ( b2.type <= 3 ){
            
            struct general_index pa ;
            pa.bra = b1;
            pa.ket = b2;
            
            return (FGS(0,1,&pa));//signs
        }else{
            printf("p-GdG?");
            exit(0);
        }
        
    }else if ( b1.basis == SincBasisElement && b2.basis == GaussianBasisElement ){
        if ( b2.type <= 3 ){
            
            struct general_index pa ;
            pa.bra = b1;
            pa.ket = b2;
            
            return -(FGS(0,1,&pa));//signs
        }else{
            printf("p-GdG?");
            exit(0);
        }
        
    }
    
    printf("over rails\n");
    exit(0);
    return 0;
}

DCOMPLEX BgB (double beta, struct basisElement b1, INT_TYPE action , INT_TYPE powSpace,double origin, struct basisElement b2){
    
        struct general_2index g2;
        double realpart,imagepart;
        g2.i[0].bra = b1;
        g2.i[0].ket = b2;
        g2.i[0].bra.type = 1;
        g2.i[1].bra.type = 1;

        g2.i[1].bra.basis = DiracDeltaElement;
        g2.i[1].ket.basis = nullBasisElement;
        g2.i[1].bra.origin = origin;
    
        g2.gaussianAccelerationFlag = 0;
        g2.momentumShift = 0;
        g2.i[0].pointer = 2;
        g2.i[1].pointer = 2;
        g2.i[0].action = 0;
        g2.i[1].action = action;
        g2.powSpace = powSpace;
        g2.momentumShift = 0;
        g2.realFlag = 1;
        realpart =  collective(sqrt(beta), &g2);
        g2.realFlag = 0;
        imagepart = collective(sqrt(beta), &g2);
        return (realpart + I * imagepart);
};


double Bd2B (struct basisElement b1, struct basisElement b2){
    if ( b1.type == nullComponent || b2.type == nullComponent    )
    {
        printf("null\n");
        exit(0);
    }
    if ( b1.type != b2.type ){
        printf("compo");
        exit(0);
    }
    
    if ( b1.basis == SincBasisElement && b2.basis == SincBasisElement ){
        if ( b2.type <= 3 ){
            double arg = b1.length*b1.index + b1.origin - (b2.length*b2.index + b2.origin);
            if ( fabs(arg) < 1e-6 ){
                if ( b1.length > b2.length ){
                    arg /= b1.length;
                    return sqrt(b2.length/b1.length)/sqr(b1.length) * ( - sqr(pi)* Sinc(1,arg) + 2 *(pi*pi/3. - pi * pi * pi * pi/30.*arg*arg + pi * pi * pi * pi*pi*pi/840.*arg*arg*arg*arg)  );
                }else {
                    arg /= b2.length;
                    
                    return sqrt(b1.length/b2.length)/sqr(b2.length) * ( - sqr(pi)* Sinc(1,arg) + 2 *(pi*pi/3. - pi * pi * pi * pi/30.*arg*arg + pi * pi * pi * pi*pi*pi/840.*arg*arg*arg*arg) );
                }
            } else {
                if ( b1.length > b2.length ){
                    arg /= b1.length;
                    
                    return sqrt(b2.length/b1.length)/sqr(b1.length) * ( - sqr(pi)* Sinc(1,arg) + 2 *( Sinc(1,arg)- cos(pi*arg))/arg/arg );
                }else {
                    arg /= b2.length;
                    
                    return sqrt(b1.length/b2.length)/sqr(b2.length) * ( - sqr(pi)* Sinc(1,arg) + 2 *( Sinc(1,arg)- cos(pi*arg))/arg/arg );
                }
                
            }
        }
        else {
            //            printf("1: %f %d %d\n", b1.length,b1.index, b1.auxIndex );
            //            printf("2: %f %d %d\n", b2.length,b2.index, b2.auxIndex );
            
            double arg = b1.length*b1.index + b1.origin - (b2.length*b2.index + b2.origin);
            
            if ( b1.length > b2.length ){
                arg /= b1.length;
                //                printf("arg: %f\n", arg);
                return periodicSd2S(arg,b1.grid) * sqrt(b2.length/b1.length)/sqr(b1.length);
            }else {
                arg /= b2.length;
                //                printf("arg: %f\n", arg);
                return periodicSd2S(arg,b2.grid) * sqrt(b1.length/b2.length)/sqr(b2.length);
            }
        }
    }else if ( b1.basis == GaussianBasisElement && b2.basis == GaussianBasisElement ){
        if ( b2.type <= 3 ){
            
            struct general_index pa ;
            pa.bra = b1;
            pa.ket = b2;
            
            return aaGdnGdm(0, 2, &pa);
        }else{
            printf("p-Gd2G?");
            exit(0);
        }
        
    }else if ( b1.basis == GaussianBasisElement && b2.basis == SincBasisElement ){
        if ( b2.type <= 3 ){
            
            struct general_index pa ;
            pa.bra = b1;
            pa.ket = b2;
            return FGS(0,2,&pa);
        }else{
            printf("p-Gd2S?");
            exit(0);
        }
    }else if ( b1.basis == SincBasisElement && b2.basis == GaussianBasisElement ){
        if ( b2.type <= 3 ){
            
            struct general_index pa ;
            pa.bra = b1;
            pa.ket = b2;
            return FGS(0,2,&pa);
        }else{
            printf("p-Gd2S?");
            exit(0);
        }
    }
    
    printf("over rails\n");
    exit(0);
    return 0;
}


DCOMPLEX poly( double k ,double beta, INT_TYPE powSpace ){
    if ( ! powSpace )
        return 1.;
    else if ( powSpace == 1 )
        return 0.50*(I * k)/ beta;//here
    else if ( powSpace == 2 )
        return -0.250*(k*k - 2. *beta )/ beta/beta;
    else if ( powSpace == 3 )
        return 0.5*0.250*(-I * k)*(-6. *beta + k*k )/ beta/beta/beta;
    else if ( powSpace == 4 )
        return 0.250*0.250*(k*k*k*k -12. * k*k * beta  + 12. * beta * beta ) / beta/beta/beta/beta;
    else if ( powSpace == 5 )
        return 0.5*0.250*0.250*(I * k)*(k*k*k*k -20. * k*k * beta  + 60. * beta * beta ) / beta/beta/beta/beta/beta;
    else if ( powSpace == 6 )
        return -0.250*0.250*0.250*( k*k*k*k*k*k - 30. * k*k  *k*k * beta + 180. * k*k*beta*beta-120 * beta*beta*beta ) / beta/beta/beta/beta/beta/beta;
    else if ( powSpace == 7 )
        return 0.5*0.250*0.250*0.250*(-I * k)*( k*k*k *k*k*k - 42. * k*k*k*k * beta + 420. * k*k*beta*beta-840. * beta*beta*beta ) / beta/beta/beta/beta/beta/beta/beta;
    else if ( powSpace == 8 )
        return 0.250*0.250*0.250*0.250*(k*k*k*k *k*k*k*k - 56.* k*k*k *k*k*k*beta + 840. *k*k*k*k *beta*beta -
                                        3360.* k*k*beta*beta*beta + 1680.*beta*beta*beta*beta)/ beta/beta/ beta/beta/ beta/beta/ beta/beta;
    else
    {
        printf("pow> 8!\n");
        exit(0);
    }
}

double gaussianSinc ( double k, void * arg ){
    struct general_2index *pa = (struct general_2index *) arg;
    INT_TYPE ai = 0;

    double beta = pa->beta ;
    DCOMPLEX act = 1.;
    for ( ai = 0; ai < pa->i[1].action ; ai++)
        act *= /*I*/ k ;

    if ( pa->realFlag != 1 )
        act *= -I ;
    double va =   exp(-sqr(k /2./ beta))*0.28209479177387814/*1/(2*sqrt(pi))*//beta * creal(act* poly(k,beta*beta,pa->powSpace) *FB(-k,&pa->i[0])*FB(k,&pa->i[1]));
    if ( isnan(va) || isinf (va)){
        printf("nan or inf\n");
        exit(0);
    }
    return va;
}

void gaussianSincFunc(void * arg,size_t n,const double * x,double * y)
{
    struct general_2index *af = (struct general_2index *) arg;
    size_t i ;
    for ( i=0;i<n;i++)
    {
        y[i] = gaussianSinc(x[i],af);
    }
    
}

double mcGS ( double x [], size_t dim , void * params ){
    struct general_2index* ga = (struct general_2index*)(params);
    ga[0].beta = x[0];
    ga[1].beta = x[0];
    ga[2].beta = x[0];
    
    return gaussianSinc(x[1], ga)*gaussianSinc(x[2], ga+1)*gaussianSinc(x[3], ga+2);
    
    
}



double collective( double beta ,struct general_2index * pa){
    double value= 0.;
    pa->beta = beta;
    INT_TYPE periodic , component;
    if ( pa->i[0].bra.type != nullComponent ){
        if ( pa->i[0].bra.type > 3 )
            periodic = 1;
        else
            periodic = 0;
    }
    else {
        printf("messed up definitions\n");
        exit(0);
    }
   component =  ( (pa->i[0].bra.type -1) % COMPONENT ) ;
    
    if ( periodic && (( pa->i[0].bra.note == interactionCell && pa->i[0].ket.note == interactionCell ) || ( pa->i[1].bra.note == interactionCell && pa->i[1].ket.note == interactionCell ) )){
        if ( pa->i[0].bra.basis == SincBasisElement&&pa->i[1].bra.basis == SincBasisElement){
            INT_TYPE targParticle = 0;
            if (( pa->i[0].bra.note == interactionCell && pa->i[0].ket.note == interactionCell ) && !( pa->i[1].bra.note == interactionCell || pa->i[1].ket.note == interactionCell )){
                targParticle = 0;
            }else if (!( pa->i[0].bra.note == interactionCell || pa->i[0].ket.note == interactionCell ) && ( pa->i[1].bra.note == interactionCell && pa->i[1].ket.note == interactionCell )){
                targParticle = 1;
            }else {
                printf("interaction cell?");
                exit(0);
            }
            
            
            INT_TYPE N1 = pa->i[targParticle].bra.grid;//
            INT_TYPE N12 = (N1-1)/2;
            double d = pa->i[targParticle].bra.length;
            double kSmall = 2*pi / N1 / d;
            INT_TYPE k;
            for ( k = - N12 ; k <= N12 ; k++)
            {
                value += gaussianSinc(kSmall*k+pa->momentumShift,pa)*kSmall;
            }
            
        }else {
            printf("not completed basis\n");
            exit(0);
        }
    }else//else not periodic
        
    {
        
            if ( pa->gaussianAccelerationFlag == 0){
                double l1=0.;
                double l2=0.;
                double pint;
                double abs_error;

                
                //INIT
                
#ifdef APPLE
                quadrature_integrate_function g;
                g.fun = gaussianSincFunc;                               // Called to evaluate the function to integrate
                g.fun_arg = pa;                                  // Passed as first argument to the callback
                
                quadrature_integrate_options options;
                bzero(&options,sizeof(options));
                options.integrator = QUADRATURE_INTEGRATE_QAGS;    // Use QAGS: adaptive subdivision with convergence acceleration
                
                options.abs_tolerance = 1.0e-9;                    // Requested absolute tolerance on result
                options.rel_tolerance = 1.0e-9;                        // Requested absolute tolerance on result
                options.max_intervals = 100;                        // Max number of intervals
                options.qag_points_per_interval = 25;
                quadrature_status status;
                
#else
#ifdef GSL_LIB
                
                gsl_function F;
                F.function = &gaussianSinc;
                F.params = pa;

                gsl_integration_workspace * workspace= gsl_integration_workspace_alloc (1000);
#endif
#endif
            //INIT

                if ( (pa->i[0].bra.basis == SincBasisElement && pa->i[0].ket.basis == SincBasisElement ) || (pa->i[1].bra.basis == SincBasisElement && pa->i[1].ket.basis == SincBasisElement )  || (pa->i[0].bra.basis == SincBasisElement && pa->i[0].ket.basis == nullFunction )|| (pa->i[0].bra.basis == nullFunction && pa->i[0].ket.basis == SincBasisElement ) || (pa->i[1].bra.basis == SincBasisElement && pa->i[1].ket.basis == nullFunction )|| (pa->i[1].bra.basis == nullFunction && pa->i[1].ket.basis == SincBasisElement ) ){
                    //finite integral
                    if (pa->i[0].bra.basis == SincBasisElement)
                        l1 += pi/pa->i[0].bra.length;
                    if (pa->i[0].bra.basis == SincBasisElement)
                        l1 += pi/pa->i[0].ket.length;
                    
                    if (pa->i[1].ket.basis == SincBasisElement)
                        l2 += pi/pa->i[1].bra.length;
                    if (pa->i[1].ket.basis == SincBasisElement)
                        l2 += pi/pa->i[1].ket.length;
                    
                    if ( l1 > 0 && l2 > 0 )
                        pint = min(l1,l2);
                    else if ( l1 > 0 )
                        pint = l1;
                    else
                        pint = l2;

#ifdef APPLE
                    
                    value =  quadrature_integrate(&g, -pint+pa->momentumShift,pint+pa->momentumShift, &options, &status, &abs_error, 0, NULL);
                    if ( status != QUADRATURE_SUCCESS)
                    {
                        printf("error\n");
                        exit(0);
                    }
#else
#ifdef GSL_LIB
                    gsl_integration_qag (&F,  -pint+pa->momentumShift,  pint+pa->momentumShift, 1e-9, 1e-6,1000,6,workspace, &value, &abs_error);
                    gsl_integration_workspace_free(workspace);
#endif
#endif
                    
                }else//infinite integral!
                {
#ifdef APPLE
                    options.integrator = QUADRATURE_INTEGRATE_QAGS;    // Use QAGS: adaptive subdivision with convergence acceleration
                    value =  quadrature_integrate(&g, -INFINITY,+INFINITY, &options, &status, &abs_error, 0, NULL);
#else
                    
#ifdef GSL_LIB
                    gsl_integration_qagi (&F, 1e-9, 1e-9,1000,workspace, &value, &abs_error);
		gsl_integration_workspace_free(workspace);    
#endif
#endif
                }
            }
                
                
        else//accelerated interal
        {
            if (  pa->i[0].bra.basis == GaussianBasisElement && pa->i[1].bra.basis == GaussianBasisElement && pa->i[0].ket.basis == GaussianBasisElement && pa->i[1].ket.basis == GaussianBasisElement)
                value =  aaGGCGG(1/beta,pa)/(beta);
            
            else if (  pa->i[0].bra.basis == GaussianBasisElement && pa->i[0].ket.basis == GaussianBasisElement && pa->i[1].bra.basis == DiracDeltaElement  )
                value =  aaGGCD(1/beta,pa->i[1].bra.origin,&pa->i[0])/(beta);
            
            else if (  pa->i[0].bra.basis == GaussianBasisElement && pa->i[0].ket.basis == GaussianBasisElement && pa->i[1].ket.basis == DiracDeltaElement  )
                value =  aaGGCD(1/beta,pa->i[1].ket.origin,&pa->i[0])/(beta);
            else {
                printf("not accelerated\n");
                exit(0);
            }
            //add OTHER PRE-integrated terms here...
        }
    }//not periodic
    return value;
}

double inverseLaplaceTransform(double beta, struct function_label * fl){
    double value2;
    value2 = 0.0;
    if ( fl->fn == Yukawa ){
        double m = fl->param[2];
        value2  += 2/sqrt(pi)*exp(-sqr(m/2./beta));
    }
    else if ( fl->fn == Morse ){
        double R = fl->param[2];
        double a = fl->param[3];
        value2  += - a * exp (     R * a - sqr(a/beta /2.))/sqr(beta);
        value2  +=   a * exp ( 2 * R * a - sqr(a/beta    ))/sqr(beta);
        value2 *= 2/sqrt(pi);
    }
    else if ( fl->fn == Coulomb || fl->fn == Pseudo || fl->fn == nullFunction ){
        if ( beta > fl->param[2] && fl->fn == Pseudo )
            value2 = 0.;
        else
            value2 += 2./sqrt(pi);//
    } else if ( fl->fn == LennardJones ){
        double rm = fl->param[1];
        if ( beta > fl->param[2] && fl->fn == Pseudo )
            value2 = 0.;
        else
            value2 += 2./sqrt(pi) * ( 1./120 * pow(rm,12) * pow(beta,11) - 2*pow(rm,6) * pow(beta,5) ) ;
    }
    return value2*fl->param[0];
}

double collectives (double beta , struct general_2index * pa ){
    
     if (  pa->point == 1){
        double l,r;
        pa->i[0].pointer = 0;
        l = collective(beta, pa);
        pa->i[0].pointer = 1;
        r = collective(beta, pa);
        return l*r;
    }
    else if ( pa->point == 2 ){
        pa->i[0].pointer = 2;
        pa->i[1].pointer = 2;

        return collective ( beta , pa );
    }
    else{
        printf("what no point?");
    
        exit(0);
    }
    }



double element ( double beta, void * aAf){
    struct general_2index *af = (struct general_2index *) aAf;
    
    
    
    return inverseLaplaceTransform(beta,af->fl)*collectives(beta,af)*collectives(beta,af+1)*collectives(beta,af+2);
}


void elementFunc(void * arg,size_t n,const double * x,double * y)
{
    struct general_2index *af = (struct general_2index *) arg;
    size_t i;
    for ( i=0;i<n;i++)
    {
        y[i] = element(x[i], af);
    }
    
}

//double monteCarloElementCal (double beta, struct general_2index *aAf  ){
//    double res= 0.;
//
//    INT_TYPE dim;
//    double xl[4],xu[4];
//
//    xl[0] = 1e-6;
//    xu[0] = beta;
//
//    xl[1] = -aAf[0].boundary;
//    xu[1] = aAf[0].boundary;
//
//    xl[2] = -aAf[1].boundary;
//    xu[2] = aAf[1].boundary;
//
//    xl[3] = -aAf[2].boundary;
//    xu[3] = aAf[2].boundary;
//
//    dim = 4;
//#ifndef APPLE
//    rank = 0;
//    gsl_monte_function G;
//    G.f = &mcGS;
//    G.dim = dim;
//    G.params = aAf;
//
//    INT_TYPE iter = 0;
//    const gsl_rng_type * T;
//    gsl_rng  * r;
//    gsl_monte_vegas_state * s;
//    gsl_rng_env_setup();
//    T = gsl_rng_default;
//    r = gsl_rng_alloc(T);
//    s = gsl_monte_vegas_alloc(G.dim);
//
//    gsl_monte_vegas_integrate ( & G , xl, xu, dim , 1000, r,s,&res, &err );
//    do
//    {
//        gsl_monte_vegas_integrate ( & G , xl, xu, dim , 1000, r,s,&res, &err );
//        printf("%f %f %f\n", res,err, gsl_monte_vegas_chisq (s));
//    } while (fabs ((gsl_monte_vegas_chisq (s) - 1.0) > 0.1 || err/max(1.,res) > 1e-9)&& num++ < 100);
//
//    gsl_monte_vegas_free(s);
//    gsl_rng_free(r);
//#else
//
//    //    {
//    //        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
//    //
//    //        gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
//    //                                   &res, &err);
//    //        display_results ("vegas warm-up", res, err);
//    //
//    //        printf ("converging...\n");
//    //
//    //        do
//    //        {
//    //            gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
//    //                                       &res, &err);
//    //            printf ("result = % .6f sigma = % .6f "
//    //                    "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
//    //        }
//    //        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//    //
//    //        display_results ("vegas final", res, err);
//    //
//    //        gsl_monte_vegas_free (s);
//    //    }
//#endif
//
//    return res;
//}

double elementCal (double a, double b,struct general_2index * aAf ){
#ifdef APPLE
    quadrature_integrate_function g;
    g.fun = elementFunc;                               // Called to evaluate the function to integrate
    g.fun_arg = aAf;                                  // Passed as first argument to the callback
    
    quadrature_integrate_options options;
    bzero(&options,sizeof(options));
    options.integrator = QUADRATURE_INTEGRATE_QAGS;    // Use QAGS: adaptive subdivision with convergence acceleration
    
    options.abs_tolerance = 1.0e-12;                    // Requested absolute tolerance on result
    options.max_intervals = 100;                        // Max number of intervals
    //options.qag_points_per_interval = 25;
    quadrature_status status;
    double  value,abs_error;
    if ( b < a)
        value =  quadrature_integrate(&g, 1e-15, INFINITY, &options, &status, &abs_error, 0, NULL);
    else
        value =  quadrature_integrate(&g, a, b, &options, &status, &abs_error, 0, NULL);
#else
    double value, abs_error;

    gsl_function F;
    F.function = &element;
    F.params = aAf;
    size_t neval;
    gsl_integration_workspace * workspace= gsl_integration_workspace_alloc (1000);

    if ( b < a)
        gsl_integration_qagiu(&F, a, 1e-9, 1e-9,1000,workspace, &value, &abs_error);
    else
        gsl_integration_qag (&F,  a,  b, 1e-9, 1e-9,1000,4,workspace, &value, &abs_error);
    gsl_integration_workspace_free(workspace);

#endif
    return value;

}






//void mySeparateExactOne (struct field * f1, double scalar, enum division basis){
//
//
//    INT_TYPE i,I1,I2,I3,I4,alpha,sp,space;
//    INT_TYPE N1[SPACE];
//    length1(f1,N1);
//    INT_TYPE N2 = N1*N1;
//    INT_TYPE N12 = (N1-1)/2;
//    double value,d = f1->sinc.d;
//    if ( scalar < 0 )
//        scalar = -pow(fabs(scalar), 1./SPACE);
//    else
//        scalar =  pow(fabs(scalar), 1./SPACE);
//
//    f1->sinc.tulip[diagonalCube].header = Cube;
//
//
//    tClear(f1, quadCube);
//    f1->sinc.tulip[quadCube].Current[0] = 1;
//
//    for ( alpha = 0 ; alpha < CanonicalRank(f1, interactionExchange, 0) ; alpha++){
//        zero(f1, quadCube,0);
//
//#pragma omp parallel for private (I1,I2,I3,I4,value,sp,space) schedule(dynamic,1)
//        for ( I1 = 0 ;I1 < N1 ;I1++){// body 0 -1
//            for ( I2 = 0 ;I2 < N1 ;I2++)// body 0 -2
//                for ( I3= 0 ;I3 < N1 ;I3++)// body 1 -1
//                    for ( I4 = 0 ;I4 < N1 ;I4++)// body 1 -2
//                        for ( space = 0; space < SPACE ; space++){
//                            value  = streams(f1, interactionExchange, 0,space)[ N2*N2 * alpha + (N2*(N1*I1+I3)+(N1*I2+I4))];
//                            streams(f1, quadCube,0,space)[(N2*(N1*I1+I2)+(N1*I3+I4))]=value*scalar;
//                        }
//        }
//
//        if ( basis ){
//           // tTransformTo(f1, quadCube, quad,basis );
//            tAddTwo(f1, interactionDirect, quad);
//        } else {
//            tAddTwo(f1, interactionDirect, quadCube);
//        }
//    }
//}

double gaussQuad(INT_TYPE pt , INT_TYPE nm, INT_TYPE which ){
    double gk3X [] = {
        -0.7745966692414833770359,
        0,
        0.7745966692414833770359
    };
    
    double gk3W [] = {
        0.5555555555555555555556,
        0.8888888888888888888889,
        0.555555555555555555556
    };
    
    double gk7X [] = {
        0.949107912342759
        ,0.741531185599394
        ,0.405845151377397
        ,0.0
        ,-0.405845151377397
        ,-0.741531185599394
        ,-0.949107912342759
    };//7
    
    double gk7W [] = {
        0.129484966168870
        ,0.279705391489277
        ,0.381830050505119
        ,0.417959183673469
        ,0.381830050505119
        ,0.279705391489277
        ,0.129484966168870};
    
    
    double gk10X[] = {
        -0.9840853600948424644962,
        -0.9061798459386639927976,
        -0.7541667265708492204408,
        -0.5384693101056830910363,
        -0.2796304131617831934135,
        0,
        0.2796304131617831934135,
        0.5384693101056830910363,
        0.7541667265708492204408,
        0.9061798459386639927976
    };
    
    double gk10W[]= {
        0.04258203675108183286451,
        0.1152333166224733940246,
        0.186800796556492657468,
        0.2410403392286475866999,
        0.272849801912558922341,
        0.2829874178574912132043,
        0.272849801912558922341,
        0.2410403392286475866999,
        0.1868007965564926574678,
        0.1152333166224733940246
    };
    
    
    
    
    double gk15X [] = {
        0.991455371120813
        ,0.949107912342759
        ,0.864864423359769
        ,0.741531185599394
        ,0.586087235467691
        ,0.405845151377397
        ,0.207784955007898
        ,0.0
        ,-0.207784955007898
        ,-0.405845151377397
        ,-0.586087235467691
        ,-0.741531185599394
        ,-0.864864423359769
        ,-0.949107912342759
        ,-0.991455371120813
    };//8
    
    
    double gk15W [] = {
        0.022935322010529
        ,0.063092092629979
        ,0.104790010322250
        ,0.140653259715525
        ,0.169004726639267
        ,0.190350578064785
        ,0.204432940075298
        ,0.209482141084728
        ,0.204432940075298
        ,0.190350578064785
        ,0.169004726639267
        ,0.140653259715525
        ,0.104790010322250
        ,0.063092092629979
        ,0.022935322010529
    };
    
    double gk35X [] = {-0.9984329706060580765167,
        -0.990575475314417335675,
        -0.9746592569674310674486,
        -0.9506755217687677612227,
        -0.9190961368038916732426,
        -0.880239153726985902123,
        -0.834274092850134363076,
        -0.7815140038968014069252,
        -0.722472287372409906023,
        -0.6576711592166907658503,
        -0.5875692123340352515726,
        -0.512690537086476967886,
        -0.4336872952097993710771,
        -0.3512317634538763152972,
        -0.2659465074516820191091,
        -0.1784841814958478558507,
        -0.0895856394252266354625,
        0.00,
        0.0895856394252266354625,
        0.1784841814958478558507,
        0.2659465074516820191091,
        0.3512317634538763152972,
        0.4336872952097993710771,
        0.5126905370864769678863,
        0.5875692123340352515726,
        0.6576711592166907658503,
        0.7224722873724099060232,
        0.7815140038968014069252,
        0.834274092850134363076,
        0.880239153726985902123,
        0.9190961368038916732426,
        0.9506755217687677612227,
        0.9746592569674310674486,
        0.9905754753144173356754,
        0.9984329706060580765167};
    
    double gk35W[] = {0.0042189757937769386559,
        0.0117858375622890860476,
        0.020022233953295123733,
        0.02785672245786342694433,
        0.0352497467518800321001,
        0.04244263020500089117947,
        0.0494361418239590678556,
        0.0559940445300935170357,
        0.062000915268229960367,
        0.067528016718131089923,
        0.072589890114190143846,
        0.0770562230462021636789,
        0.0808366784408178713044,
        0.0839725719912347754276,
        0.0864905322056519368263,
        0.0883087822976044066468,
        0.08936184586778888574733,
        0.0896964219439813653622,
        0.0893618458677888857473,
        0.0883087822976044066468,
        0.086490532205651936826,
        0.08397257199123477542758,
        0.080836678440817871304,
        0.0770562230462021636789,
        0.0725898901141901438458,
        0.067528016718131089923,
        0.062000915268229960367,
        0.05599404453009351703567,
        0.0494361418239590678556,
        0.0424426302050008911795,
        0.0352497467518800321001,
        0.0278567224578634269443,
        0.020022233953295123733,
        0.0117858375622890860476,
        0.00421897579377693865586};
    
    double gk99X [] ={-0.999804199685638727995,
        -0.9988201506066353793618,
        -0.996819814299826469401,
        -0.9937886619441677907601,
        -0.9897642384140644710632,
        -0.984757895914213004359,
        -0.9787573030312066845299,
        -0.97176220090155538014,
        -0.9637900181363554282611,
        -0.954853658674137233555,
        -0.944955059221329473156,
        -0.934100294755810149059,
        -0.922305475453936168755,
        -0.909585655828073285213,
        -0.8959496245835934207352,
        -0.881408445573008910037,
        -0.8659800043151457264443,
        -0.8496821198441657010349,
        -0.832528501487460168517,
        -0.8145344273598554315395,
        -0.7957203207625248381361,
        -0.7761068943454466350181,
        -0.7557118903695143082457,
        -0.734554254237402696214,
        -0.7126570667050088308474,
        -0.6900438244251321135048,
        -0.666735700841571667866,
        -0.642754832419237664057,
        -0.6181268178008792991612,
        -0.5928776941089007124559,
        -0.5670315494953917328631,
        -0.5406132469917260665582,
        -0.5136506284017343823111,
        -0.486171941452492042177,
        -0.4582036915390298548496,
        -0.4297729933415765246586,
        -0.4009095739292798809373,
        -0.3716435012622848888637,
        -0.3420031959918559601078,
        -0.3120175321197487622079,
        -0.2817177082410294178775,
        -0.2511351786125772735072,
        -0.2202997655098053353243,
        -0.1892415924618135864853,
        -0.157992877664368358666,
        -0.126585997269672051068,
        -0.09505164998612337566,
        -0.0634206849826867860288,
        -0.03172586345907315363082,
        0.0,
        0.03172586345907315363082,
        0.06342068498268678602884,
        0.09505164998612337566,
        0.126585997269672051068,
        0.157992877664368358666,
        0.189241592461813586485,
        0.2202997655098053353243,
        0.2511351786125772735072,
        0.2817177082410294178775,
        0.3120175321197487622079,
        0.3420031959918559601078,
        0.3716435012622848888637,
        0.4009095739292798809373,
        0.4297729933415765246586,
        0.4582036915390298548496,
        0.486171941452492042177,
        0.5136506284017343823111,
        0.540613246991726066558,
        0.567031549495391732863,
        0.5928776941089007124559,
        0.618126817800879299161,
        0.6427548324192376640569,
        0.666735700841571667866,
        0.6900438244251321135048,
        0.712657066705008830847,
        0.7345542542374026962137,
        0.7557118903695143082457,
        0.7761068943454466350181,
        0.7957203207625248381361,
        0.81453442735985543154,
        0.8325285014874601685173,
        0.8496821198441657010349,
        0.8659800043151457264443,
        0.881408445573008910037,
        0.8959496245835934207352,
        0.909585655828073285213,
        0.9223054754539361687554,
        0.934100294755810149059,
        0.944955059221329473156,
        0.9548536586741372335552,
        0.9637900181363554282611,
        0.97176220090155538014,
        0.97875730303120668453,
        0.9847578959142130043593,
        0.9897642384140644710632,
        0.9937886619441677907601,
        0.9968198142998264694013,
        0.9988201506066353793618,
        0.999804199685638727995};
    
    double gk99W[] = {0.0005274769683783323143,
        0.0014779639281743620209,
        0.0025218765731496496845,
        0.0035329557014832599892,
        0.0045141331625998620836,
        0.005501420399338078204,
        0.0064998725321664845168,
        0.0074869102764140445198,
        0.0084551856397750456701,
        0.0094175722296862066762,
        0.010378732924116607707,
        0.0113278438578780228795,
        0.0122591748799473589077,
        0.0131792068212079366057,
        0.0140911110472705440377,
        0.0149880989346802956593,
        0.01586572369887289313,
        0.016727898553777318682,
        0.017576872641244826238,
        0.0184077538519825820281,
        0.0192169329826556442117,
        0.0200070643859274929265,
        0.0207798535198561693337,
        0.0215314818077816867538,
        0.022258913930129946636,
        0.0229641216311680166298,
        0.0236484943738631350276,
        0.02430890292981194970869,
        0.0249427312116380137806,
        0.0255515676298983967962,
        0.0261366295727170561997,
        0.0266952731071580377427,
        0.027225205808796711996,
        0.0277278075508688728381,
        0.0282042205700292005214,
        0.0286521667109019998496,
        0.0290696140121504756931,
        0.0294578454056420734076,
        0.0298179973237107521896,
        0.030148081066799933145,
        0.0304462797396858877386,
        0.030713855857739083392,
        0.030951992615566392528,
        0.03115893749659741587,
        0.0313330525912310597589,
        0.03147563582699890351852,
        0.031587959458685973161,
        0.03166846605661171584462,
        0.031715664503973344041,
        0.0317309313985218944091,
        0.031715664503973344041,
        0.0316684660566117158446,
        0.031587959458685973161,
        0.031475635826998903519,
        0.0313330525912310597589,
        0.03115893749659741587,
        0.03095199261556639252752,
        0.030713855857739083392,
        0.0304462797396858877386,
        0.030148081066799933145,
        0.0298179973237107521896,
        0.029457845405642073408,
        0.0290696140121504756931,
        0.0286521667109019998496,
        0.02820422057002920052136,
        0.027727807550868872838,
        0.027225205808796711996,
        0.026695273107158037743,
        0.0261366295727170561997,
        0.025551567629898396796,
        0.0249427312116380137806,
        0.0243089029298119497087,
        0.0236484943738631350276,
        0.02296412163116801663,
        0.02225891393012994663599,
        0.0215314818077816867538,
        0.020779853519856169334,
        0.0200070643859274929265,
        0.0192169329826556442117,
        0.0184077538519825820281,
        0.017576872641244826238,
        0.016727898553777318682,
        0.01586572369887289313,
        0.014988098934680295659,
        0.0140911110472705440377,
        0.0131792068212079366057,
        0.012259174879947358908,
        0.0113278438578780228795,
        0.010378732924116607707,
        0.0094175722296862066762,
        0.0084551856397750456701,
        0.0074869102764140445198,
        0.0064998725321664845168,
        0.005501420399338078204,
        0.0045141331625998620836,
        0.0035329557014832599892,
        0.002521876573149649685,
        0.0014779639281743620209,
        0.00052747696837833231426};
    double *gkX,*gkW;
    INT_TYPE ngk;
    if ( pt == 3 ){
        gkX = gk3X;
        gkW = gk3W;
        ngk = 3;
    }else
    if ( pt == 7 ) {
        gkX = gk7X;
        gkW = gk7W;
        ngk = 7;
    }
        else if ( pt == 10 ){
            gkX = gk10X;
            gkW = gk10W;
            ngk = 10;
        
    }else if ( pt == 15 ){
        gkX = gk15X;
        gkW = gk15W;
        ngk = 15;
    }else if ( pt == 35 ){
        gkX = gk35X;
        gkW = gk35W;
        ngk = 35;
    }else if ( pt == 99 ){
        gkX = gk99X;
        gkW = gk99W;
        ngk = 99;
    }else {
        exit(0);
    }
    
    if ( which )
        return 0.5*(1+gkX[nm]);
    else
        return 0.5*gkW[nm];//shift to interval [0,1]
}

void mySeparateExactTwo (struct field * f1, enum division interactionExchange, double scalar,  enum division basis,INT_TYPE overline, INT_TYPE particle1){//overline = ##particle number ...  = intercellular interaction
    //https://keisan.casio.com/exec/system/1329114617
    zero(f1,interactionExchange,0);
    
    if ( f1->twoBody.func.fn == nullFunction)
        return ;
    
    if(0){
        INT_TYPE comp[COMPONENT+1],c,space,flagc=1,dim=0;
        
        for ( c = 0; c <= COMPONENT; c++)
            comp[c]= 0;
        for ( space = 0; space < SPACE ; space++)
            if ( f1->sinc.rose[space].body >= one && f1->sinc.rose[space].particle == particle1 ){
                dim++;
        
            if ( f1->sinc.tulip[interactionExchange].space[space].body == two ){
                comp[f1->sinc.rose[space].component] += f1->sinc.tulip[interactionExchange].space[space].body;
                if (f1->sinc.rose[space].body < two )
                    flagc = 0;
            }else if ( f1->sinc.tulip[interactionExchange].space[space].body != nada ){
                flagc = 0;
            }
        }
        for ( c = 1 ; c <= dim ; c++)
            if ( 2 != comp[c] )
                flagc = 0;
        
        if ( ! flagc ){
            printf("error in two setup \n");
            exit(1);
            
        }
    }
    INT_TYPE flagPow;
    long double cpow ;
    INT_TYPE spaces;

    INT_TYPE I1,I2,I3,I4,beta,section;
    INT_TYPE n1[SPACE];
    INT_TYPE n2[SPACE];
    length(f1, interactionExchange,n2);
    length1(f1,n1);
    INT_TYPE N1,space;
    spaces = 0;
    for ( space = 0 ;space < SPACE  ; space++)
        if ( f1->sinc.rose[space].body != nada )
            if ( f1->sinc.rose[space].component != nullComponent )
                if ( f1->sinc.rose[space].particle == particle1 )

                if ( f1->sinc.tulip[interactionExchange].space[space].body == two )
            spaces++;

    
    
    enum functionType fn = f1->twoBody.func.fn;
    if ( fn == nullFunction )
        return ;
    
    struct general_2index g2;
    g2.realFlag = 1;
    double value,g,x;
    double constant;
    struct function_label fl = f1->twoBody.func;
    INT_TYPE si,interval = fl.interval;
    double * param   = fl.param;
    getDescription(&fl, scalar, stdout);
    
    INT_TYPE ngk;
    f1->sinc.tulip[quadCube].header = Cube;
    
    tClear(f1, quadCube);
    tId(f1,quadCube,0 );
    
    for ( section = 0 ; section < 2 ; section++){
        if ( imax(0,interval - section) == 0 ) {
                ngk = 7;
            }else if ( imax(0,interval - section) == 1 ){
                ngk = 15;
            }else if ( imax(0,interval - section) == 2 ){
                ngk = 35;
            }else if ( imax(0,interval - section) == 3 ){
                ngk = 99;
            }else {
                printf("oops\n");
                exit(0);
            }

            
            
        
            for ( beta = 0; beta < ngk ; beta++){
                
                flagPow = 1;

                g = gaussQuad(ngk,beta,1);
                constant = gaussQuad(ngk, beta, 0);
                
                if ( section == 0 ){
                    x = g ;
                }
                else {
                    x = ( g ) / (1. - g)+1 ;
                    constant /= sqr(1.-g);
                }
                x *=  param[1];
                constant *= param[1];

                //divide 0->infinity
                
                constant *= scalar*inverseLaplaceTransform(x,&fl);//for purposes of scaling tw-body interactions by Ne...
                cpow = powl(fabsl(constant),1./spaces);
                for ( space = 0; space < SPACE ; space++)
                    if ( f1->sinc.rose[space].body != nada )
                    if ( f1->sinc.rose[space].component != nullComponent )
                    if ( f1->sinc.tulip[interactionExchange].space[space].body == two && f1->sinc.rose[space].particle == particle1 )
                    {
      					//printf("%d %d\n", beta,space);
				//fflush(stdout);
		                  N1 = n1[space];
#ifdef OMP
#pragma omp parallel for private (si,I1,I2,I3,I4,g2,value)
#endif
                        for ( si = 0 ; si < N1*N1*N1*N1; si++)
                            
                        {
                            I1 = ( si ) % N1;
                            
                            I3 = ( si / N1) % N1;
                            I2 = ( si /  (N1*N1) ) % N1;
                            I4 = ( si / ( N1*N1*N1) ) % N1;
                            
                            g2.realFlag = 1;
                            g2.momentumShift = 0;
                            g2.gaussianAccelerationFlag = 0;
                            g2.point = 2;
                            g2.powSpace = 0;
                            
                            g2.i[0].action = 0;
                            g2.i[1].action = 0;
                            
                            
                            g2.i[0].bra = grabBasis ( f1, space, particle1, I1);
                            g2.i[0].ket = grabBasis ( f1, space, particle1, I2);
                            g2.i[1].bra = grabBasis ( f1, space, particle1, I3);
                            g2.i[1].ket = grabBasis ( f1, space, particle1, I4);
                            if ( overline ){
                                g2.i[overline].bra.note = interactionCell ;
                                g2.i[overline].ket.note = interactionCell ;
                            }

                            
                            value = collectives(x, &g2)*cpow;
                            if (flagPow && constant < 0 )
                                value *= -1.;
                            streams(f1, quadCube,0,space)[si]=value;
                        }
                        flagPow = 0;

                    }
                tAddTw(f1, interactionExchange, 0,quadCube,0);
            }
    }
    return ;
}


void mySeparateEwaldCoulomb1(struct field * f1,INT_TYPE nVec, double *  occupy,enum division vectors, INT_TYPE part1, enum division interactionExchange, enum division interactionEwald,enum division shorten, double scalar,INT_TYPE plus,double rescale, enum particleType particle){
    INT_TYPE info,dim,rank=0,l,vol,vos=0,x,r,rr = 0,n1[SPACE],space,j1,j2,i1,i2,vor,vor2, vox,lll;
    length1(f1,n1);
    double sumDis =0 ;
    Stream_Type * streamIn, *streamOut;
    enum division vo = vectors ,current;
    tClear(f1,shorten);
    
    for ( lll = 0; lll < 2 ; lll++){
        for ( vo = vectors ; vo < vectors+nVec ; vo++){
            vox = CanonicalRank(f1, vo, 0);
            current = interactionExchange;
            if ( lll == 1 )
                current = interactionEwald;

            for ( r = 0 ; r < CanonicalRank(f1, current, 0); r++){


                tClear(f1,copy);
                tClear(f1, diagonalCube );
                tId(f1, diagonalCube, 0);//set protons to 1
                
                for ( vor = 0 ; vor < vox; vor++){
                    for ( vor2 = 0 ; vor2 < vox; vor2++){
                        tId(f1, copy,0);
                    }
                }
                
                
                for ( dim = 0 ;dim < SPACE ; dim++)
                    if ( f1->sinc.rose[dim].body != nada && f1->sinc.rose[dim].particle == particle){
                        streamIn = streams(f1,current,0,dim)+r*n1[dim]*n1[dim]*n1[dim]*n1[dim];
                        f1->sinc.tulip[diagonalCube].space[dim].block = tv1;

#ifdef OMP
#pragma omp parallel for private (j1,j2,i1,i2,vor,vor2,streamOut) schedule(dynamic,1)
#endif
                        for ( j1 = 0 ; j1 < n1[dim] ; j1++){
#ifdef OMP
                            rank = omp_get_thread_num();
#else
                            rank = 0;
#endif
                            streamOut = streams(f1,diagonalCube,rank,dim);


                            
                            for ( j2 = 0; j2 < n1[dim];j2++){

                                for ( i1 = 0 ; i1 < n1[dim]; i1++)
                                    for ( i2 = 0; i2 < n1[dim]; i2++){
                                        streamOut[i2*n1[dim]+i1] = streamIn[(j2*n1[dim]+i2)*n1[dim]*n1[dim] + j1*n1[dim]+i1];
                                    }
                                for ( vor = 0 ; vor < vox; vor++){
                                    f1->sinc.tulip[canonicalmeVector].Current[rank] =0;
                                    tGEMV(rank, f1, dim, canonicalmeVector, rank, diagonalCube, 0, rank, vo, vor, 0);
                                    f1->sinc.tulip[canonicalmeVector].Current[rank] =1;
                                    for ( vor2 = 0 ; vor2 < vox; vor2++){
                                        ( streams(f1, copy,0, dim) + n1[dim]*n1[dim]*(vox * vor2 + vor))[j2*n1[dim]+j1] = occupy[vo]*tDOT(rank, f1, dim, 'T', vo, vor2, 0, 'N', canonicalmeVector, 0, rank);
                                    }
                                }
                            }
                        }
                    }
                    tClear(f1,copyTwo);
                    if ( lll == 0 )
                        tScale(f1, copy, -1.);

                    tCycleDecompostionListOneMP(-1, f1, copy, 0,NULL, copyTwo, 0, f1->mem1->rt->CANON,part1, 1);
                tMultiplyMP(0, &info, f1, 1., -1, copy, 0, 'T', copyTwo, 0, 'N', copyTwo, 0);
                if (fabs( traceOne(f1, copy, 0)) > 1e-9){
                    tAddTw(f1, shorten, 0, copyTwo, 0);
                    //sumDis += sqrt(distance(f1, copy, copyTwo));
                }
                printf("Ewald ++%f = %f %d\n", traceOne(f1, copy, 0),traceOne(f1, shorten,0),CanonicalRank(f1, shorten, 0));

            }
        }
    }
}




void mySeparateExactOneByOne (struct field * f1, INT_TYPE part1, enum division interactionExchangePlus,enum division shorten, double scalar,INT_TYPE plus,double rescale, enum particleType particle1,enum particleType particle2){
    //https://keisan.casio.com/exec/system/1329114617
    zero(f1,interactionExchangePlus,0);
    struct name_label n;
    INT_TYPE Mt=0,mm,xx,dim;
    INT_TYPE cra;
    double sumDis = 0.;
    if ( f1->twoBody.func.fn == nullFunction)
        return ;

    if(0){
        INT_TYPE comp[COMPONENT+1],c,space,flagc=1,dim=0;
        
        for ( c = 0; c <= COMPONENT; c++)
            comp[c]= 0;
        for ( space = 0; space < SPACE ; space++)
            if ( f1->sinc.rose[space].body >= one && f1->sinc.rose[space].particle == particle1 ){
                dim++;
                if ( f1->sinc.tulip[interactionExchangePlus ].space[space].body == one ){
                    comp[f1->sinc.rose[space].component] += f1->sinc.tulip[interactionExchangePlus].space[space].body;
                    if (f1->sinc.rose[space].body < one )
                        flagc = 0;
                }else if ( f1->sinc.tulip[interactionExchangePlus].space[space].body != nada ){
                    flagc = 0;
                }
            }
        for ( c = 1 ; c <= dim ; c++)
            if ( comp[c] > 2  )
                flagc = 0;
        
        if ( ! flagc ){
            printf("error in two setup \n");
            exit(1);
            
        }
    }
    INT_TYPE flagPow;;
    long double cpow ;
    INT_TYPE spaces;
    spaces = 0;

    INT_TYPE I1,I2 = 0,I3,I4 = 0,beta,section,space2;
    INT_TYPE n1[SPACE];
    INT_TYPE N2;
    length1(f1,n1);
    INT_TYPE N1,space;
    
    for ( space = 0 ;space < SPACE  ; space++)
        if ( f1->sinc.rose[space].body != nada )
            if ( f1->sinc.rose[space].component != nullComponent )
            if ( f1->sinc.tulip[interactionExchangePlus].space[space].body >= one )
                if ( f1->sinc.rose[space].particle == particle1 || f1->sinc.rose[space].particle == particle2 )

            spaces++;

    
    enum functionType fn = f1->twoBody.func.fn;
    if ( fn == nullFunction )
        return ;
    
    struct general_2index g2;
    
    g2.realFlag = 1;
    double value,g,x;
    double constant;
    INT_TYPE si,interval = f1->twoBody.func.interval;
    double * param   = f1->twoBody.func.param;
    getDescription(&f1->twoBody.func, scalar, stdout);
    struct function_label fl = f1->twoBody.func;
    
    INT_TYPE ngk,co;
    tClear(f1, interactionExchangePlus);
    for ( section = 0 ; section < 2 ; section++){
        if ( imax(0,interval - section) == 0 ) {
            ngk = 7;
        }else if ( imax(0,interval - section) == 1 ){
            ngk = 15;
        }else if ( imax(0,interval - section) == 2 ){
            ngk = 35;
        }else if ( imax(0,interval - section) == 3 ){
            ngk = 99;
        }else {
            printf("oops\n");
            exit(0);
        }
        
        
        

        for ( beta = 0; beta < ngk ; beta++){
            tClear(f1, copy);

            flagPow = 1;
            
            
            g = gaussQuad(ngk,beta,1);
            constant = gaussQuad(ngk, beta, 0);
            
            if ( section == 0 ){
                x = g ;
            }
            else {
                x = ( g ) / (1. - g)+1 ;
                constant /= sqr(1.-g);
            }
            x *=  param[1];
            constant *= param[1];
            
            
            myZero(f1, oneByOneBuffer,0);
            constant *= scalar*inverseLaplaceTransform(x,&fl);//for purposes of scaling tw-body interactions by Ne...
            Mt= 0;
            cpow = powl(fabsl(constant),1./COMPONENT);

            for ( co = 1 ; co <= COMPONENT ; co++)
                for ( space = 0 ; space < SPACE ; space++)
                    if ( f1->sinc.rose[space].component == co && f1->sinc.rose[space].particle == particle1)
                        for (space2 = space+1 ; space2 < SPACE ; space2++)
                            if (  co == f1->sinc.rose[space2].component && f1->sinc.rose[space2].particle == particle2 )
                                if (space != space2 ){
                                    if (f1->sinc.rose[space].body != nada &&  f1->sinc.rose[space2].body != nada && f1->sinc.tulip[interactionExchangePlus].space[space].body == one && f1->sinc.tulip[interactionExchangePlus].space[space2].body == one  )
                                        {
                                            N1 = n1[space];
                                            N2 = n1[space2];
#ifdef OMP
#pragma omp parallel for private (si,I1,I2,I3,I4,g2,value) schedule(dynamic,1)
#endif
                                            for ( si = 0 ; si < N1*N1*N2*N2 ; si++)
                                            {
                                                
                                                I1 = ( si ) % N1;
                                                I2 = ( si / N1) % N2;
                                                I3 = ( si /  (N2*N1) ) % N1;
                                                I4 = ( si / ( N2*N1*N1) ) % N2;
                                                g2.realFlag = 1;
                                                g2.momentumShift = 0;
                                                g2.gaussianAccelerationFlag = 0;
                                                g2.point = 2;
                                                g2.powSpace = 0;
                                                
                                                g2.i[0].action = 0;
                                                g2.i[1].action = 0;
                                                
                                                g2.i[0].bra = grabBasis ( f1, space, f1->sinc.rose[space].particle, I1);
                                                g2.i[0].ket = grabBasis ( f1, space, f1->sinc.rose[space].particle, I3);
                                                g2.i[1].bra = grabBasis( f1, space2, f1->sinc.rose[space2].particle, I2);
                                                g2.i[1].ket = grabBasis( f1, space2, f1->sinc.rose[space2].particle, I4);
                                                g2.i[1].bra.length *= rescale;
                                                g2.i[1].ket.length *= rescale;
                                                if ( plus == -1 ){
                                                    g2.i[1].bra.origin *= -1.;
                                                    g2.i[1].ket.origin *= -1.;
                                                    g2.i[1].bra.index *= -1;
                                                    g2.i[1].ket.index *= -1;
                                                }
                                                value = collectives(x, &g2)*cpow;
                                                if (flagPow && constant < 0 )
                                                    value *= -1.;
                                                myStreams(f1, oneByOneBuffer,0)[N1*N1*(N2*I4+I2)+(N1*I3+I1)]=value;
                                            }
                                            
                                            Mt = tdgesvd(0, f1, N1*N1,N2*N2, myStreams(f1, oneByOneBuffer,0), streams(f1, copy,0,space), streams(f1, copy,0,space2));
                                        
                                            if (Mt !=  imin(N1*N1,N2*N2)){
                                                printf("congruity\n");
                                                exit(0);
                                            }
                                            
                                        }else
                                        
                                        
                                        {//default to space 1.
                                            if (f1->sinc.rose[space].body == nada || f1->sinc.tulip[interactionExchangePlus].space[space].body == nada  )
                                            {
                                                printf("oops\n");//try changing order of particle inputs
                                                exit(0);
                                            }

                                            
                                            N1 = n1[space];
#ifdef OMP
#pragma omp parallel for private (si,I1,I3,g2,value) schedule(dynamic,1)
#endif

                                            for ( si = 0 ; si < N1*N1 ; si++)
                                            {
                                                
                                                I1 = ( si ) % N1;
                                                I3 = ( si / N1) % N1;
                                                g2.realFlag = 1;
                                                g2.momentumShift = 0;
                                                g2.gaussianAccelerationFlag = 0;
                                                g2.point = 2;
                                                g2.powSpace = 0;
                                                
                                                g2.i[0].action = 0;
                                                g2.i[1].action = 0;
                                                
                                                g2.i[1].bra.type = f1->sinc.rose[space].component;
                                                g2.i[1].ket.type = f1->sinc.rose[space].component;

                                                g2.i[0].bra = grabBasis ( f1, space, f1->sinc.rose[space].particle, I1);
                                                g2.i[0].ket = grabBasis ( f1, space, f1->sinc.rose[space].particle, I3);
                                                g2.i[1].bra.basis = DiracDeltaElement ;
                                                g2.i[1].bra.origin = 0;
                                                g2.i[1].ket.basis = nullBasisElement;
                                                value = collectives(x, &g2)*cpow;
                                                if (flagPow && constant < 0 )
                                                    value *= -1.;
                                                myStreams(f1, oneByOneBuffer,0)[(N1*I3+I1)]=value;
                                            }
                                            if ( Mt == 0 ){
                                                printf("order");
                                                exit(1);
                                            }
                                            for ( mm = 0 ; mm < Mt ;mm++){
                                                cblas_dcopy(N1*N1, myStreams(f1, oneByOneBuffer,0), 1, streams(f1, copy,0,space)+mm*N1*N1, 1);
                                            }

                                        }
                                    flagPow = 0;

                                    }
            
            f1->sinc.tulip[copy].Current[0] = Mt;
            if ( part1 ){
                tClear(f1,copyTwo);
                tCycleDecompostionListOneMP(-1, f1, copy, 0,NULL, copyTwo, 0, f1->mem1->rt->CANON,part1, 1);
                tAddTw(f1, shorten, 0, copyTwo, 0);
                sumDis += sqrt(distance(f1, copy, copyTwo));
            }else {
                tAddTw(f1, interactionExchangePlus, 0, copy, 0);
            }
        }
        
    }
    if ( part(f1,interactionExchangePlus) < CanonicalRank(f1, interactionExchangePlus, 0)){
        printf("whatgives?\n");
        exit(9);
    }
    printf("max error %f percent\n", 100*sumDis/fabs(traceOne(f1,shorten,0)) );
    return ;
}



INT_TYPE separateExternal( struct calculation * c1,enum division linear, INT_TYPE periodic, INT_TYPE atom,double scalar, INT_TYPE dim, enum division basis , INT_TYPE particle1){
    //https://keisan.casio.com/exec/system/1329114617
    struct field * f1 = &c1->i.c;
    zero(f1,linear,0);
    
    if ( f1->oneBody.func.fn == nullFunction)
        return 0;
    
    INT_TYPE flagPow;
    long double cpow ;
    INT_TYPE spaces;
   
    if(0){
        INT_TYPE comp[COMPONENT+1],c,space,flagc=1,dim= 0;
        
        for ( c = 0; c <= COMPONENT ; c++)
            comp[c]= 0;
        for ( space = 0; space < SPACE ; space++)
            if ( f1->sinc.rose[space].body >= one && f1->sinc.rose[space].particle == particle1 ){
                dim++;
                if ( f1->sinc.tulip[linear].space[space].body == one ){
                    comp[f1->sinc.rose[space].component] += f1->sinc.tulip[linear].space[space].body;
                    if ( f1->sinc.rose[space].body < one )
                        flagc = 0;
                }
                else if ( f1->sinc.tulip[linear].space[space].body != nada )
                    flagc = 0;
            }
        for ( c = 1 ; c <= dim ; c++)
            if ( 1 != comp[c] )
                flagc = 0;
        
        if ( ! flagc ){
            printf("error in linear setup \n");
            exit(1);
            
        }
    }
    

    const INT_TYPE psu_stride = 17;
    const INT_TYPE psu_length = 11;
    double lda[] = {/*test*/
       // 2, 1, 0, 1,      10, 1.2, 0,0,0,1,     0, 1, 1, 0, 0, 1.5, 1,/*test*/
        1, 1, 0, 2, 2, 0.2, -4.06633, 0.677832, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        3, 3, 0, 4, 4, 0.4, -14.0094, 9.50991, -1.75327, 0.0834586, 0, 0, 0,
        0, 0, 0, 0, 4, 4, 0, 4, 4, 0.325, -23.991, 17.1718, -3.31896,
        0.165083, 0, 0, 0, 0, 0, 0, 0, 5, 3, 0, 3, 2, 0.4325, -5.60048,
        0.806284, 0, 0, 1, 0.373882, 6.23522, 0, 0, 0, 0, 6, 4, 0, 3, 2,
        0.346473, -8.57533, 1.23413, 0, 0, 1, 0.304523, 9.53419, 0, 0, 0, 0,
        7, 5, 0, 3, 2, 0.288905, -12.2046, 1.75582, 0, 0, 1, 0.256912,
        0.256912, 0, 0, 0, 0, 8, 6, 0, 3, 2, 0.247754, -16.4822, 2.37014, 0,
        0, 1, 0.222203, 18.1996, 0, 0, 0, 0, 13, 3, 0, 6, 1, 0.45, -6.83406,
        0, 0, 0, 2, 0.465436, 2.81408, 1.93952, 3, 0.546243, 1.91601, 14, 4,
        0, 6, 1, 0.44, -6.91363, 0, 0, 0, 2, 0.424334, 3.20813, 2.58888, 3,
        0.485359, 2.65622, 15, 5, 0, 6, 1, 0.43, -6.64097, 0, 0, 0, 2,
        0.390738, 3.65826, 3.15066, 3, 0.440846, 3.28594, 16, 6, 0, 6, 1,
        0.42, -6.59607, 0, 0, 0, 2, 0.362614, 4.22284, 3.66966, 3, 0.405311,
        3.88535};
    
    double blyp[] = {1, 1, 0, 2, 2, 0.2, -4.10561, 0.692787, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        3, 3, 0, 4, 4, 0.4, -14.1026, 9.65027, -1.79063, 0.0857313, 0, 0, 0,
        0, 0, 0, 0, 4, 4, 0, 4, 4, 0.325, -24.0586, 17.2529, -3.33239,
        0.165305, 0, 0, 0, 0, 0, 0, 0, 5, 3, 0, 3, 2, 0.424087, -6.08744,
        0.980916, 0, 0, 1, 0.371141, 6.32735, 0, 0, 0, 0, 6, 4, 0, 3, 2,
        0.337633, -9.12847, 1.42513, 0, 0, 1, 0.302528, 9.65073, 0, 0, 0, 0,
        7, 5, 0, 3, 2, 0.281959, -12.7548, 1.94859, 0, 0, 1, 0.255444,
        13.6594, 0, 0, 0, 0, 8, 6, 0, 3, 2, 0.24245, -17.0171, 2.56133, 0, 0,
        1, 0.221084, 18.3556, 0, 0, 0, 0, 13, 3, 0, 6, 1, 0.45, -5.54822, 0,
        0, 0, 2, 0.505838, 3.02008, 1.06418, 3, 0.577572, 1.53528, 14, 4, 0,
        6, 1, 0.44, -5.97966, 0, 0, 0, 2, 0.444927, 3.4402, 1.88129, 3,
        0.503637, 2.28821, 15, 5, 0, 6, 1, 0.43, -5.87283, 0, 0, 0, 2,
        0.403545, 3.8762, 2.54131, 3, 0.452751, 2.9405, 16, 6, 0, 6, 1, 0.42,
        -6.0083, 0, 0, 0, 2, 0.370401, 4.37362, 3.19573, 3, 0.413079, 3.5911};

    
    
    
    INT_TYPE section;
    INT_TYPE dims1[SPACE];
    length1(f1,dims1);

    INT_TYPE gPoint = 2,i,beta,I1,space,I2,a,ii,spin;
    spaces = 0;
    for ( space = 0 ;space < SPACE  ; space++)
        if ( f1->sinc.rose[space].body != nada )
            if ( f1->sinc.rose[space].component != nullComponent )
                if ( f1->sinc.tulip[linear].space[space].body == one )
                    if ( f1->sinc.rose[space].particle == particle1 )
                        spaces++;

    double constant,value,x,g;
    Stream_Type  *stream[SPACE];
    struct function_label fl;
    struct general_2index g2;
    struct name_label u = f1->sinc.tulip[diagonalCube];
    f1->sinc.tulip[diagonalCube].header = Cube;
    for ( i = 0; i < SPACE ; i++)
        if ( f1->sinc.rose[i].body != nada )
        stream[i] =  streams( f1, diagonalCube,0,i  );
    zero(f1, diagonalCube , 0);
    f1->sinc.tulip[diagonalCube].Current[0] = 1;
    
    
    /////IGNORE TRAIN COMMAND
    INT_TYPE *n1,*m1;
    if ( basis ){
        f1->sinc.tulip[linear].header = f1->sinc.tulip[basis].header-1;
        n1 = dims1;
        m1 = n1;
        f1->sinc.tulip[diagonalCube].header = Cube;
    }
    //    else if ( header (f1, linear ) == SincSub  ){
    //        n1 = f1->sinc.n1;
    //        m1 = f1->sinc.n1;
    //        f1->sinc.tulip[diagonalCube].header = SincSub;
    //        f1->sinc.tulip[linear].header = SincSub;
    //    }
    else {
        n1 = dims1;
        m1 = n1 ;
        f1->sinc.tulip[diagonalCube].header = Cube;
        f1->sinc.tulip[linear].header = Cube;
        
    }
    spin = 0;
   // f1->sinc.tulip[linear].stop[0][0] = CanonicalRank(f1, linear, 0);
    enum functionType fn = f1->oneBody.func.fn;
    if ( fn == nullFunction )
        return 0;
    
    
    double va,*psu=NULL,*PSU=NULL;
    if ( fn == LDA )
        psu = lda;
    else if ( fn == BLYP )
        psu = blyp;
    
    INT_TYPE interval = f1->oneBody.func.interval;
    double * param   = f1->oneBody.func.param;

    double *gkX, *gkW,Z, gGX[10],gGW[10];
    INT_TYPE si,ngk,ai,powSpace[SPACE][10];
    
    struct basisElement boa;
    INT_TYPE mA,xA;
    if ( ! atom ){
        xA = f1->Na;
        mA = 1;
    } else if ( atom > 0 ) {
        mA = atom;
        xA = atom;
    } else {
        printf("separateExternal\n");

        exit(0);
    }
    printf("\t  Load %d Nuclear Fields |\n",xA);
    for ( a = mA ; a <= xA; a++){
        ii = 0;
        printf("\t  %1.3f\t %1.3f\t %1.3f \t|\t",getPosition(f1, a, 1),getPosition(f1, a, 2),getPosition(f1, a, 3));
        
        //section 0    long range
        //section 1    SHORT RANGE
        
        
        //GTH pseudo-potentials!!!
        //section 2     Gaussians-straight
        //section 3     Gaussians-r^2
        //section 4     Gaussians-r^4
        //section 5     Gaussians-r^8
        //section 6     s-wave--straight
        //section 7     s-wave--r^2
        //section 8     p-wave--straight
        if ( fn == LDA || fn == BLYP ){
            if ( psu == NULL ){
                printf("separateExternal 2\n");
                
                exit(0);
            }
            PSU = NULL;
            for ( ai = 0; ai < psu_length ; ai++){
                if ( psu[psu_stride*ai] ==  f1->atoms[a].label.Z)
                    PSU = psu+psu_stride * ai;
            }
            if ( PSU == NULL )
            {
                printf ("not coded\n");
                exit(0);
            }
            
            Z = PSU[1];
            
            fl.interval = interval;
            fl.fn = Pseudo;
            //fprintf(outString,"Pseudo = %1.3f Erf(r/%1.3f)/r @ %1.3f\n",fn->param[0],fn->param[2], fn->param[1]);
            fl.param[0] = 1;
            fl.param[2] = sqrt(2.)*PSU[5];
            fl.param[1] = 1./(sqrt(2.)*PSU[5]);
            
            
            
        }
        else{
            fl.interval = interval;
            fl.param[0] = param[0];
            fl.param[1] = param[1];
            fl.param[2] = param[2];
            fl.param[3] = param[3];
            fl.fn = fn;
            Z = f1->atoms[a].label.Z;
        }
        getDescription(&fl, Z * scalar, stdout);
        
        
        for ( section = 0 ; section < 9 ; section++){
            ngk = 0;
            if ( (( fn == LDA || fn == BLYP )&& section == 0) || ( ( fn != LDA && fn != BLYP ) && section < 2)  ){
                gPoint = 2;
                
                
                
                
                if ( imax(0,interval - section) == 0 ) {
                    ngk = 7;
                }else if ( imax(0,interval - section) == 1 ){
                    ngk = 15;
                }else if ( imax(0,interval - section) == 2 ){
                    ngk = 35;
                }else if ( imax(0,interval - section) == 3 ){
                    ngk = 99;
                }else {
                    printf("oops\n");
                    exit(0);
                }
            } else if ( fn == LDA || fn == BLYP ){
                gkW = gGW;
                gkX = gGX;
                if ( 2 <= section && section <= 5 ){
                    if ( ! ( PSU[4] > section - 2 ) )
                        continue;
                    gPoint = 2;
                    
                    if ( section == 2 ){
                        ngk = 1;
                        for ( ai = 0 ; ai < ngk ; ai++){
                            gGX[ai] = 1./sqrt(2.)/PSU[5];
                            gGW[ai] = -PSU[section + 4]/Z*pow(sqr(sqrt(2.)*gGX[ai]), section - 2);
                        }
                        for ( space = 0 ; space < SPACE ; space++)
                            powSpace[space][0] = 0;
                        
                    }
                    else if ( section == 3 ){
                        ngk = 3;
                        for ( ai = 0 ; ai < ngk ; ai++){
                            gGX[ai] = 1./sqrt(2.)/PSU[5];
                            gGW[ai] = -PSU[section + 4]/Z*pow(sqr(sqrt(2.)*gGX[ai]), section - 2);
                        }
                        //a^2 + b^2 + c^2
                        powSpace[0][0] = 2;
                        powSpace[1][0] = 0;
                        powSpace[2][0] = 0;
                        
                        powSpace[0][1] = 0;
                        powSpace[1][1] = 2;
                        powSpace[2][1] = 0;
                        
                        powSpace[0][2] = 0;
                        powSpace[1][2] = 0;
                        powSpace[2][2] = 2;
                        
                    }
                    else if ( section == 4 ){
                        ngk = 6;
                        for ( ai = 0 ; ai < ngk ; ai++){
                            gGX[ai] = 1./sqrt(2.)/PSU[5];
                            gGW[ai] = -PSU[section + 4]/Z*pow(sqr(sqrt(2.)*gGX[ai]), section - 2);
                        }
                        //a^2 + 2 a b + b^2 + 2 a c + 2 b c + c^2
                        gGW[0]  *= 1;
                        powSpace[0][0] = 4;
                        powSpace[1][0] = 0;
                        powSpace[2][0] = 0;
                        
                        gGW[1]  *= 2;
                        powSpace[0][1] = 2;
                        powSpace[1][1] = 2;
                        powSpace[2][1] = 0;
                        
                        gGW[2]  *= 1;
                        powSpace[0][2] = 0;
                        powSpace[1][2] = 4;
                        powSpace[2][2] = 0;
                        
                        gGW[3]  *= 2;
                        powSpace[0][3] = 0;
                        powSpace[1][3] = 2;
                        powSpace[2][3] = 2;
                        
                        gGW[4]  *= 1;
                        powSpace[0][4] = 0;
                        powSpace[1][4] = 0;
                        powSpace[2][4] = 4;
                        
                        gGW[5]  *= 2;
                        powSpace[0][5] = 2;
                        powSpace[1][5] = 0;
                        powSpace[2][5] = 2;
                        
                        
                    }
                    
                    else if ( section == 5 ){
                        ngk = 10;
                        for ( ai = 0 ; ai < ngk ; ai++){
                            gGX[ai] = 1./sqrt(2.)/PSU[5];
                            gGW[ai] = -PSU[section + 4]/Z*pow(sqr(sqrt(2.)*gGX[ai]), section - 2);
                        }
                        //a^3 + 3 a^2 b + 3 a b^2 + b^3 + 3 a^2 c + 6 a b c + 3 b^2 c +
                        //3 a c^2 + 3 b c^2 + c^3
                        
                        
                        
                        gGW[0]  *= 1;
                        powSpace[0][0] = 6;
                        powSpace[1][0] = 0;
                        powSpace[2][0] = 0;
                        
                        gGW[1]  *= 3;
                        powSpace[0][1] = 4;
                        powSpace[1][1] = 2;
                        powSpace[2][1] = 0;
                        
                        gGW[2]  *= 3;
                        powSpace[0][2] = 2;
                        powSpace[1][2] = 4;
                        powSpace[2][2] = 0;
                        
                        gGW[3]  *= 1;
                        powSpace[0][3] = 0;
                        powSpace[1][3] = 6;
                        powSpace[2][3] = 0;
                        
                        gGW[4]  *= 3;
                        powSpace[0][4] = 4;
                        powSpace[1][4] = 0;
                        powSpace[2][4] = 2;
                        
                        gGW[5]  *= 6;
                        powSpace[0][5] = 2;
                        powSpace[1][5] = 2;
                        powSpace[2][5] = 2;
                        
                        gGW[6]  *= 3;
                        powSpace[0][6] = 0;
                        powSpace[1][6] = 4;
                        powSpace[2][6] = 2;
                        
                        gGW[7]  *= 3;
                        powSpace[0][7] = 2;
                        powSpace[1][7] = 0;
                        powSpace[2][7] = 4;
                        
                        gGW[8]  *= 3;
                        powSpace[0][8] = 0;
                        powSpace[1][8] = 2;
                        powSpace[2][8] = 4;
                        
                        gGW[9]  *= 1;
                        powSpace[0][9] = 0;
                        powSpace[1][9] = 0;
                        powSpace[2][9] = 6;
                        
                        
                    }
                    
                }else if ( 6 <= section && section <= 7 ){
                    if ( ! ( PSU[10] > section - 6 ) )
                        continue;
                    
                    gPoint = 1;
                    if ( section == 6 ){
                        ngk = 1;
                        for ( ai = 0 ; ai < ngk ; ai++){
                            gGW[ai] = -4.*PSU[section + 6]/Z/sqrt(pi)/cube(PSU[11])/2.5464827020350906;
                            gGX[ai] = 1./sqrt(2.)/PSU[11];
                        }
                        powSpace[0][0] = 0;
                        powSpace[1][0] = 0;
                        powSpace[2][0] = 0;
                    }
                    else if ( section == 7 ){
                        ngk = 3;
                        for ( ai = 0 ; ai < ngk ; ai++){
                            gGW[ai] = -16./15.*PSU[section + 6]/Z/(sqrt(pi))/pow(PSU[11],7.)/2.5464827020350906;
                            gGX[ai] = 1./sqrt(2.)/PSU[11];
                        }
                        powSpace[0][0] = 2;
                        powSpace[1][0] = 0;
                        powSpace[2][0] = 0;
                        
                        powSpace[0][1] = 0;
                        powSpace[1][1] = 2;
                        powSpace[2][1] = 0;
                        
                        powSpace[0][2] = 0;
                        powSpace[1][2] = 0;
                        powSpace[2][2] = 2;
                        
                    }
                } else if ( section == 8 ){
                    if ( ! ( PSU[14] ) )
                        continue;
                    
                    gPoint = 1;
                    if ( section == 8 ){
                        ngk = 3;
                        for ( ai = 0 ; ai < ngk ; ai++){
                            gGW[ai] = -8/3.*PSU[16]/Z/sqrt(pi)/pow(PSU[15],5.)/2.5464827020350906;
                            gGX[ai] = 1./sqrt(2.)/PSU[15];
                        }
                        powSpace[0][0] = 1;
                        powSpace[1][0] = 0;
                        powSpace[2][0] = 0;
                        
                        powSpace[0][1] = 0;
                        powSpace[1][1] = 1;
                        powSpace[2][1] = 0;
                        
                        powSpace[0][2] = 0;
                        powSpace[1][2] = 0;
                        powSpace[2][2] = 1;
                        
                    }
                }
            } else {
                continue;
            }
            
                    
            for ( beta = 0; beta < ngk ; beta++){
                
                flagPow = 1;
                
                if ( section == 0 ){
                    g = gaussQuad(ngk,beta,1);
                    constant = gaussQuad(ngk, beta, 0);
                    
                    x = g ;
                }
                else if ( section == 1 ){
                    g = gaussQuad(ngk,beta,1);
                    constant = gaussQuad(ngk, beta, 0);
                    
                    x = ( g ) / (1. - g)+1 ;
                    constant /= sqr(1.-g);
                }else {
                    x = gkX[beta];
                    constant = gkW[beta];
                }
                
                x *=  param[1];
                constant *= param[1];
                
                constant *= Z*scalar;
                //this is odd but,
                if (section > 2 )
                    constant *= 2*pi;//old definition
                else
                    constant *= inverseLaplaceTransform(x,&fl);
                //end weirdness

                tClear(f1,diagonalCube);
                zero(f1,diagonalCube,0);
                tId(f1,diagonalCube,0);
                
                cpow = powl(fabsl(constant),1./spaces);

                for ( space = 0 ;space < SPACE  ; space++)
                    if ( f1->sinc.rose[space].body != nada )

                    if ( f1->sinc.rose[space].component != nullComponent )
                        
                        
                    if ( f1->sinc.tulip[linear].space[space].body == one )
                        if(f1->sinc.rose[space].particle == particle1  )
                {
#ifdef OMP
#pragma omp parallel for private (si,I1,I2,g2,value) schedule(dynamic,1)
#endif
                    for ( si = 0; si < dims1[space]*dims1[space]; si++)
                    {
                        I1 = si % dims1[space];
                        I2 = (si/dims1[space])% dims1[space];
                        
                        g2.realFlag = 1;
                        g2.momentumShift = 0;
                        g2.gaussianAccelerationFlag = 0;
                        g2.point = gPoint;
                        if ( section >= 2 )
                            g2.powSpace = powSpace[space][beta];
                        else
                            g2.powSpace = 0;
                        
                        
                        
                        g2.i[0].bra = grabBasis(f1, space, particle1, I1);
                        g2.i[0].ket = grabBasis ( f1, space, particle1, I2);
                        g2.i[1].bra.type = f1->sinc.rose[space].component;
                        g2.i[1].ket.type = f1->sinc.rose[space].component;

                        g2.i[1].bra.basis = DiracDeltaElement;
                        g2.i[1].bra.origin = f1->atoms[a].position[space%COMPONENT + 1 ];
                        g2.i[1].ket.basis = nullBasisElement;
                        g2.i[1].action = ( dim == space ) ;
                        
                        value = collectives(x, &g2)*cpow;
               //         printf("%1.15f\n", value);
                        if (flagPow && constant < 0 )
                            value *= -1.;

                        (stream[space])[I1*dims1[space]+I2] = value;
                    }
                    flagPow = 0;

                }
               // printf("trace %f\n", traceOne(f1, diagonalCube, 0));
               // exit(0);

                tAddTw(f1, linear,0, diagonalCube,0);

            }
        }
    }
    //outputFormat(f1, stdout, protonRepulsion, 0);
    return 0;
}

INT_TYPE separateKinetic( struct field * f1, INT_TYPE periodic,enum division akinetic,  double amass, INT_TYPE particle1 ){
    INT_TYPE space,dim,I1,I2;
    INT_TYPE dims1[SPACE];
    struct general_index o1;
    length1(f1,dims1);
    Stream_Type * stream;
    f1->sinc.tulip[diagonalCube].Current[0] = 1;

    for ( dim = 0 ; dim < SPACE; dim++){
        if ( f1->sinc.rose[dim].body != nada )
            if (( f1->sinc.tulip[akinetic].space[dim].body == one && f1->sinc.rose[dim].particle == particle1  )){
            for ( space = 0 ;space < SPACE; space++)
                if ( f1->sinc.rose[space].body != nada )

            {
                stream =  streams( f1, diagonalCube, 0 , space );
                for ( I1 = 0 ; I1 < dims1[space] ; I1++)
                    for( I2 = 0; I2 < dims1[space] ; I2++){
                        o1.bra = grabBasis(f1, space, particle1, I1);
                        o1.ket = grabBasis ( f1, space, particle1, I2);

                        if ( dim == space  ){
                            (stream )[dims1[space]*I1+I2] = - 0.5/amass*Bd2B(o1.bra,o1.ket);
//                            printf ("%d %d %d %f \n", space,I1,I2,Bd2B(grabBasis(f1, space, particle1, I1),grabBasis(f1, space, particle1, I2))/ Sd2S(I1-I2));

                        }
                        else{
                            (stream )[dims1[space]*I1+I2] = BoB(o1.bra,o1.ket);
//                            printf("%f\n",BoB(o1.bra,o1.ket));
                        }
                    }
            }
            tAddTw(f1, akinetic, 0, diagonalCube, 0);
        }
    }
    
    struct name_label k = f1->sinc.tulip[akinetic];
    double va = traceOne(f1, akinetic, 0);
    return 0;
}


void separateDerivatives( struct field * f1, INT_TYPE periodic,enum division mat, INT_TYPE *x, INT_TYPE *grad,double mag,INT_TYPE particle1 ){
    DCOMPLEX bca;
    double b0 = 0.0001,powSpace,spaces=0;
    INT_TYPE space,I1,I2,flagCMPL=0;
    INT_TYPE dims[SPACE],signFlag= 1;
    length1(f1, dims);
    
    Stream_Type * stream[2];
    
    
    INT_TYPE cur[2];
    cur[0] = CanonicalRank(f1, mat, 0);
    cur[1] = CanonicalRank(f1, mat, 1);
    
    tId(f1,mat,0);
    for (space = 0; space < SPACE ; space++)
        if ( f1->sinc.rose[space].body != nada )
            if ( grad[space] % 2 == 1 )
                flagCMPL = 1;
    
    
    if (flagCMPL)
        tId(f1,mat,1);
    for (space = 0; space < SPACE ; space++)
        if ( f1->sinc.rose[space].body != nada )
            spaces++;
    
        powSpace = pow( fabs(mag),1./spaces);

   // printf("mag %f\n", mag);
    for (space = 0; space < SPACE ; space++)
        if ( f1->sinc.rose[space].body != nada )
        {
            
            stream[0] =  streams( f1, mat, 0 , space ) + dims[space]*dims[space]*cur[0];//last one is COMPLEX:: cmpl!
            stream[1] =  streams( f1, mat, 1 , space ) + dims[space]*dims[space]*cur[1];//last one is COMPLEX:: cmpl!
            for ( I1 = 0 ; I1 < dims[space] ; I1++){
                for( I2 = 0; I2 < dims[space] ; I2++)
                {
                    
                    bca =  powSpace * BgB(b0,grabBasis(f1, space, particle1, I1),grad[space],x[space],0,grabBasis(f1, space, particle1, I2))/pow(b0,grad[space]);
                    if ( mag < 0  && signFlag )
                        mag *= -1;

                    (stream[0])[dims[space]*I1+I2] = creal(bca);
                    if ( flagCMPL )
                        (stream[1])[dims[space]*I1+I2] = cimag(bca);
                    
//                    if ( x[space] == 2 ){
//                        if ( I1 == I2 )
//                            (stream[0])[dims[space]*I1+I2] = mag*sqr((I1-(dims[space]-1)/2) * f1->sinc.rose[space].lattice);
//                        else
//                            (stream[0])[dims[space]*I1+I2] = 0.;
//                    }else {
//                        if ( I1 == I2 )
//                            (stream[0])[dims[space]*I1+I2] = 1;
//                        else
//                            (stream[0])[dims[space]*I1+I2] = 0.;
//
//                    }
                    
//                    if ( I1 == I2 )
//                        printf("%d %d %1.12f %1.12f\n",space,I1, creal(bca),cimag(bca));
                }
            
            }
            signFlag = 0;
        }
//    double u = traceOne(f1, vectorMomentum, 0);
//    double u2 = traceOne(f1, vectorMomentum, 1);
    return;
}





//double tTestTwoBody( struct field * f1, enum division mat,INT_TYPE periodic, INT_TYPE * p){
//    INT_TYPE r,space,N1 = f1->sinc.N1;
//
//
//    double sum = 0.,product;
//    for ( r = 0; r < CanonicalRank(f1, mat, 0 );r++){
//        product = 1;
//        for ( space = 0 ; space < SPACE ; space++)
//            product *= streams(f1, mat, 0,space)[p[4*space] + N1*p[4*space+1]+ N1*N1*p[4*space+2]+N1*N1*N1*p[4*space+3]+r*N1*N1*N1*N1];
//        sum += product;
//    }
//    INT_TYPE i ;
//    struct general_2index g3[3];
//    for ( i = 0; i < 3; i++){
//        g3[i].realFlag = 1;
//        g3[i].i[0].d = f1->sinc.d;
//        g3[i].i[1].d = f1->sinc.d;
//
//        g3[i].i[0].n = p[i*4];
//        g3[i].i[0].m = p[i*4+2];
//        g3[i].i[1].n = p[i*4+1];
//        g3[i].i[1].m = p[i*4+3];
//        g3[i].point = 2;
//        g3[i].fl = & f1->twoBody.func;
//        g3[i].powSpace = 0;
//        if ( i == 0 )
//            g3[i].periodic = ( periodic ) % 2;
//        else if ( i == 1 )
//            g3[i].periodic = ( periodic / 2 ) % 2;
//        else if (i == 2 )
//            g3[i].periodic = ( periodic / 4) % 2;
//        else
//        {
//            printf("here\n");
//            exit(0);
//        }
//        g3[i].body = 2;
//        g3[i].N1 = N1;
//    }
//   // printf("%1.15f \t %1.15f\n", sum, elementCal(1e-3,-1, g3));
//    {
//        double value;
//        value = sum - elementCal(1e-3,-1, g3);
//        if ( isnan(value)){
//            printf(".");
//            return 0.;
//        } else
//            return sqr(value);
//    }
//}

//double tRMSDevRandom( struct field * f1, enum division mat, INT_TYPE periodic ,INT_TYPE Nc){
//    INT_TYPE p[12],i,j;
//    double sum = 0.;
//    for ( i = 0; i < Nc ; i++){
//        for ( j = 0 ; j < 12 ; j++)
//            p[j] = rand()% f1->sinc.N1;
//        sum += tTestTwoBody(f1, mat, periodic, p);
//    
//    }
//    return sqrt(sum/Nc);
//}

#if 0
//CEG paper

INT_TYPE buildElectronProtonInteraction ( struct field * f1, enum division mat){
    INT_TYPE space,r,i,j,n,m,N1 = f1->sinc.N1, N2 = f1->sinc.N1*f1->sinc.N1;
    double value;
    if ( ! CanonicalRank(f1, interactionExchange, 0))
        return 0;
    for ( r = 0; r < CanonicalRank(f1, interactionExchange, 0); r++){
        
        for ( i = 0; i < N1 ; i++)
            for ( j = 0; j < N1; j++)
                for (space = 0; space < SPACE ;space++)

            {
                value = 0.;
                for ( n = 0 ; n < N1 ; n++)
                    for ( m = 0 ; m < N1 ; m++)
                        value += streams(f1, interactionExchange, 0,space)[(N2*(N1*i+n)+(N1*j+m))+r*N2*N2];
                    streams(f1, mat,0,space)[N1*i+j+r*N2] =value/N1;
            }
    }
    f1->sinc.tulip[mat].Current[0] = CanonicalRank(f1, interactionExchange, 0);
    tScaleOne(f1, mat, 0,  -f1->Ne);
    
    
    if ( bodies(f1, eigenVectors)==two){
        tClear(f1, copy);
        tId(f1, copy,0);
        tScale(f1, copy,-3.*traceOne(f1, mat, 0)/(pow(N1,SPACE))/sqr(f1->Ne));
        //
        tAddTw(f1, mat, 0,copy,0);
    }
    else     if ( bodies(f1, eigenVectors)==three){
        tClear(f1, copy);
        tId(f1, copy,0);
        tScale(f1, copy,-6.*traceOne(f1, mat, 0)/(pow(N1,SPACE))/sqr(f1->Ne));
        //
        tAddTw(f1, mat,0, copy,0);
    }
    else     if ( bodies(f1, eigenVectors)==four){
        tClear(f1, copy);
        tId(f1, copy,0);
        tScale(f1, copy,-10.*traceOne(f1, mat, 0)/(pow(N1,SPACE))/sqr(f1->Ne));
        //
        tAddTw(f1, mat, 0,copy,0);
    }
    
    
    f1->sinc.tulip[mat].stop[0][0] = CanonicalRank(f1, mat,0);
    
    return 0;
}
#endif

INT_TYPE buildElectronProtonInteraction ( struct field * f1, enum division mat, INT_TYPE spin){
    INT_TYPE space,r,i,j,n,m,N1,N2 ;
    INT_TYPE n1[space];
    length1(f1,n1);
    double value;
    if ( ! CanonicalRank(f1, interactionExchange, 0))
        return 0;
    for ( r = 0; r < CanonicalRank(f1, interactionExchange, 0); r++){
        
        for ( i = 0; i < N1 ; i++)
            for ( j = 0; j < N1; j++)
                for (space = 0; space < SPACE ;space++)
                    
                {
                    N1 = n1[space] ;
                    N2 = N1*N1;
                    value = 0.;
                    for ( n = 0 ; n < N1 ; n++)
                        for ( m = 0 ; m < N1 ; m++)
                            value += streams(f1, interactionExchange, 0,space)[(N2*(N1*i+n)+(N1*j+m))+r*N2*N2];
                    streams(f1, mat,spin,space)[N1*i+j+r*N2] =value/N1/*NORM OF UNIFORM WAVEFUNCTION*/;
                }
    }
    f1->sinc.tulip[mat].Current[spin] = CanonicalRank(f1, interactionExchange, 0);
    tScaleOne(f1, mat, spin,  -f1->Ne);

    return 0;
}

INT_TYPE tZeroSum ( struct field * f1, enum division mat,INT_TYPE spin){
    INT_TYPE space,r,n,m,N1 , N2 ,n2,m2,flag = 1;
    INT_TYPE n1[SPACE];
    length1(f1,n1);
    double value,prod = 1,psum,sum = 0.;
    for ( r = 0; r < CanonicalRank(f1, mat, spin); r++){
        prod = 1.;
        for (space = 0; space < SPACE ;space++)
            if ( f1->sinc.rose[space].body != nada ){
              //  printf("%d %d\n",space,f1->sinc.rose[space].component );
                if ( f1->sinc.rose[space].component > 3 )
        
                {
                    flag = 0;
                    N1= n1[space];
                    N2 = N1*N1;
                    if ( bodies(f1, mat ) == one)
                        
                        
                        
                    {
                        value = 0.;
                        for ( n = 0 ; n < N1 ; n++)
                            for ( m = 0 ; m < N1 ; m++)
                                value += streams(f1, mat, spin,space)[(N1*(n)+(m))+r*N2];
                        prod *= value;
                    }
                    else
                        if ( bodies(f1, mat ) == two )
                            
                        {
                            value = 0.;
                            for ( n = 0 ; n < N1 ; n++)
                                for ( m = 0 ; m < N1 ; m++)
                                    for ( n2 = 0 ; n2 < N1 ; n2++)
                                        for ( m2 = 0 ; m2 < N1 ; m2++)
                                            value += streams(f1, mat, spin,space)[(N2*(N1*n2+n)+(N1*m2+m))+r*N2*N2];
                            prod *= value;
                            
                        }
                    sum += prod;

                }
            }
    }
    if (! flag ) {
        enum division label;
        if ( bodies(f1, mat) == one )
            label = diagonalCube;
        else if ( bodies(f1, mat) == two )
            label = quadCube;
        else
        {
            printf("zerosum \n");
            exit(0);
        }
        if ( mat == label ){
            printf("overlap\n");
            exit(0);
        }
        tClear(f1, label);
        tId(f1, label, 0);
        tScaleOne(f1, label, spin, -sum/traceOne(f1, label, 0));
        tAddTw(f1, mat, spin, label, 0);
        printf("%d-tr %f\n", mat,traceOne(f1, mat, spin));
    }else {
      //  printf("failed %d\n", flag);
    }
    
    return 0;
}


void separateX ( struct field * f1,  double vectorDipole ){
    INT_TYPE space, i;
    INT_TYPE N1;
    INT_TYPE n1[SPACE];
    length1(f1,n1);
    tClear(f1, X);
    zero(f1, X,0);
    if ( fabs(vectorDipole) < 1e-6 )
        return;
    tId(f1, X,0);
    for ( space = 0; space < 1 ; space++){
        N1 = n1[space];
        for ( i = 0; i < N1 ; i++){
            //streams(f1,X,0,space)[i*N1+i+space*N1*N1] = f1->sinc.d*(i-(N1-1)/2)*vectorDipole;
        }
    }
    return;
}//build X final
