/*
 *  coreForce.c
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

#include "coreForce.h"


void getDescription ( struct function_label *fn ,double scalar,FILE * outString){
    if ( fn->fn == nullFunction){
        fprintf(outString,"nullFunction\n");
    }else if ( fn->fn == Pseudo ){
        fprintf(outString,"Pseudo = %1.3f Erf(r/%1.3f)/r @ %1.3f\n",fn->param[0],fn->param[2], fn->param[1]);
    }else if ( fn->fn == Yukawa ){
        fprintf(outString,"Yukawa = %1.3f exp(- r %1.3f)/r @ %1.3f\n",fn->param[0],fn->param[2], fn->param[1]);
    }else if ( fn->fn == Coulomb ){
        fprintf(outString,"Coulomb = %1.3f /r @ %1.3f\n",fn->param[0],fn->param[1]);
    }else if ( fn->fn == Morse ){
        fprintf(outString,"Morse = %1.3f (1-exp[-%1.3f *( r - %1.3f )] )^2 -1 @ %1.3f\n",fn->param[0],fn->param[3],fn->param[2],fn->param[1]);
    }else if ( fn->fn == LJ ){
        //fprintf(outString,"LJ = %1.3f (..+1./r^6) @ %1.3f\n",fn->param[0],fn->param[2],fn->param[1],fn->param[1]);
    }
}



DCOMPLEX ei ( double arg ){
    return cexp(I*arg);
}

double No(double beta1){
    return 1./sqrt(sqrt(pi / 2. / beta1 ) ) ;
};//single dimensions of guassian



double GoG( double beta1, double beta2 , double x ){
    double va = sqrt(pi/(beta1+beta2))*No(beta1)*No(beta2)*exp(-(beta1*beta2)/(beta1+beta2)*x*x);
    return va;
}


DCOMPLEX FGG ( double p , struct general_index * pa){
    double xc1,xc2;
    DCOMPLEX fgg;
    double b0,b1,x0,x1;
    INT_TYPE l0,l1;
    double x;
    
    
    if ( pa->l0 < pa->l1 ){
        b0 = pa->b1;
        b1 = pa->b0;
        l0 = pa->l1;
        l1 = pa->l0;
        x0 = pa->x1;
        x1 = pa->x0;
        
    }else {
        b1 = pa->b1;
        b0 = pa->b0;
        l1 = pa->l1;
        l0 = pa->l0;
        x1 = pa->x1;
        x0 = pa->x0;

    }
    x = x1-x0;
    xc1 =GoG(b0, b1, x)*exp(-sqr(p /2.)/(b0+b1));
        
    xc2 = ((b0*x0+ b1*x1)/(b0+b1));
    fgg = xc1 * ei( p * xc2 );

    if ( l0 == 1 ){
        fgg /= pow(b0,0.5*l0);
    }else if ( l0 == 2 ){
        fgg /= sqrt(3.)*pow(b0,0.5*l0);
    }else if ( l0 == 3 ){
        fgg /= sqrt(15.)*pow(b0,0.5*l0);
    }else if ( l0 == 4 ){
        fgg /= sqrt(105.)*pow(b0,0.5*l0);
    }
    
    if ( l1 == 1 ){
        fgg /= pow(b1,0.5*l1);
    }else if ( l1 == 2 ){
        fgg /= sqrt(3.)*pow(b1,0.5*l1);
    }else if ( l1 == 3 ){
        fgg /= sqrt(15.)*pow(b1,0.5*l1);
    }else if ( l1 == 4 ){
        fgg /= sqrt(105.)*pow(b1,0.5*l1);
    }

    if ( l0 == 0 && l1 == 0 ){
        return fgg;
    }
    else if ( l0 == 1 && l1 == 0 ){
        return fgg*(b0)/(b0+b1) * (  2. * b1 * x +I* p  );
    } else if ( l0 == 1 && l1 == 1 ){
        return -fgg*(b0*b1)/sqr(b0+b1) * (( p*p - 2.*(b1+b0)+4.*b0*b1*x*x -I*( 2. * ( b1 - b0 ) * x * p)));
        
        
    }  else if ( l0 == 2 && l1 ==    0 ){
        return fgg*(b0*(-2*b1*b1 - b0*(p*p - 4*b1*b1*x*x + b1*(2 - 4*I*p*x))))/(b0 + b1)/(b0 + b1);
    }  else if ( l0 == 2 && l1 ==    1 ){
        return (fgg/(b0 + b1)*(b0 + b1)*(b0 + b1))*(b0*b1*(-2*I*b1*b1*p + b0*((-I)*p*p*p + 2*b1*p*(I - 2*p*x) +
                                                                              4*b1*b1*x*(3 + I*p*x)) + 2*b0*b0*(p*p*x + 2*b1*x*(3 - 2*b1*x*x) -
                                                                                                                2*I*p*(-1 + 2*b1*x*x))));
    }  else if ( l0 == 2 && l1 ==    2 ){

        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b1*(2*b1*b1*b1*p*p+ b0*b1*(p*p*p*p + 2*b1*p*p*(-3 - 2*I*p*x) -
                                                           4*b1*b1*(-3 - 6*I*p*x + p*p*x*x)) +
                                       2*b0*b0*b0*(p*p + 8*b1*b1*b1*x*x*x*x + 8*I*b1*b1*x*x*(3*I + p*x) -
                                               2*b1*(-3 + 6*I*p*x + p*p*x*x)) + 2*b0*b0*b1*(8*b1*b1*x*x*(-3 - I*p*x) +
                                                                                            p*p*(-3 + 2*I*p*x) + 4*b1*(3 + 2*p*p*x*x))));
        
        
    }  else if ( l0 == 3 && l1 ==    0 ){
        return fgg*(b0*b0*(I*p + 2*b1*x)*(-6*b1*b1 -
                                      b0*(p*p - 4*b1*b1*x*x + b1*(6 - 4*I*p*x))))/(b0 + b1)/(b0 + b1)/(b0 + b1);
        
    }  else if ( l0 == 3 && l1 ==    1 ){
       return  (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*(6*b1*b1*(p*p + b1*(-2 - 2*I*p*x)) +
                                  b0*(p*p*p*p - 6*I*b1*p*p*p*x + 8*b1*b1*b1*x*x*(6 + I*p*x) -
                                      12*b1*b1*(2 - 2*I*p*x + p*p*x*x)) - 2*b0*b0*(8*b1*b1*b1*x*x*x*x + p*p*(3 - I*p*x) +
                                                                                   12*I*b1*b1*x*x*(2*I + p*x) - 6*b1*(-1 + 3*I*p*x + p*p*x*x))));
    }  else if ( l0 == 3 && l1 ==    2 ){
        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*(6*b1*b1*b1*p*(I*p*p + 2*b1*(-2*I + p*x)) +
                                             b0*b1*(I*p*p*p*p*p + 6*b1*p*p*p*(-I + p*x) + 8*b1*b1*b1*x*(15 + 12*I*p*x - p*p*x*x) -
                                                    12*I*b1*b1*p*(1 - 5*I*p*x + p*p*x*x)) +
                                             2*b0*b0*b0*(I*p*p*p + 16*b1*b1*b1*b1*x*x*x*x*x + 8*b1*b1*b1*x*x*x*(-10 + 3*I*p*x) +
                                                     2*b1*p*(9*I + 9*p*x - I*p*p*x*x) - 12*b1*b1*x*(-5 + 6*I*p*x + p*p*x*x)) +
                                             2*b0*b0*b1*(16*b1*b1*b1*x*x*x*(-5 - I*p*x) - p*p*p*(5*I + 2*p*x) +
                                                        24*b1*b1*x*(5 - I*p*x + p*p*x*x) + 6*I*b1*p*(4 + 3*I*p*x + 2*p*p*x*x))));
        
    }  else if ( l0 == 3 && l1 ==    3 ){
        return -((fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*b1*(6*b1*b1*b1*p*p*(p*p+ b1*(-6 - 2*I*p*x)) -
                                            6*b0*b0*b1*(p*p*p*p*(2 - I*p*x) - 6*b1*p*p*(4 + p*p*x*x) +
                                                       8*b1*b1*b1*x*x*(-15 - 10*I*p*x + p*p*x*x) +
                                                       12*b1*b1*(5 + 5*I*p*x + 4*p*p*x*x +
                                                                I*p*p*p*x*x*x)) +
                                            b0*b1*(p*p*p*p*p*p + 6*b1*p*p*p*p*(-2 - I*p*x) -
                                                   12*b1*b1*p*p*(-3 - 8*I*p*x + p*p*x*x) +
                                                   8*b1*b1*b1*(-15 - 45*I*p*x + 18*p*p*x*x +
                                                           I*p*p*p*x*x*x)) +
                                            4*b0*b0*b0*b0*(16*b1*b1*b1*b1*x*x*x*x*x*x + 3*I*p*p*(3*I + p*x) +
                                                    24*I*b1*b1*b1*x*x*x*x*(5*I + p*x) -
                                                    12*b1*b1*x*x*(-15 + 10*I*p*x + p*p*x*x) +
                                                    b1*(-30 + 90*I*p*x + 36*p*p*x*x - 2*I*p*p*p*x*x*x)) +
                                            
                                            6*b0*b0*b0*(p*p*p*p + 16*b1*b1*b1*b1*x*x*x*x*(-5 - I*p*x) +
                                                    24*b1*b1*b1*x*x*(10 + p*p*x*x) -
                                                    2*b1*p*p*(-3 + 8*I*p*x + p*p*x*x) +
                                                    12*I*b1*b1*(5*I + 5*p*x + 4*I*p*p*x*x +
                                                                p*p*p*x*x*x)))));

    }  else if ( l0 == 4 && l1 ==    0 ){
        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*(12*b1*b1*b1*b1 -
                                      12*b0*b1*b1*(-p*p + 4*b1*b1*x*x + b1*(-2 + 4*I*p*x)) +
                                      b0*b0*(p*p*p*p + 16*b1*b1*b1*b1*x*x*x*x + 4*b1*p*p*(3 - 2*I*p*x) +
                                            16*b1*b1*b1*x*x*(-3 + 2*I*p*x) -
                                             12*b1*b1*(-1 + 4*I*p*x + 2*p*p*x*x))));
        
    }  else if ( l0 == 4 && l1 ==    1 ){
        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*
                                b1*(12*I*b1*b1*b1*b1*p + 12*b0*b1*b1*(I*p*p*p + 2*b1*b1*x*(-5 - 2*I*p*x) +
                                                              2*b1*p*(-I + 2*p*x)) +
                                    I*b0*b0*(p*p*p*p*p + 4*b1*p*p*p*(1 - 2*I*p*x) +
                                            16*b1*b1*b1*b1*x*x*x*(-10*I + p*x) +
                                            16*b1*b1*b1*x*(15*I + 9*p*x + 2*I*p*p*x*x) -
                                            12*b1*b1*p*(7 - 2*I*p*x + 2*p*p*x*x)) -
                                    2*
                                    b0*b0*b0*(16*b1*b1*b1*b1*x*x*x*x*x + 16*b1*b1*b1*x*x*x*(-5 + 2*I*p*x) +
                                          p*p*p*(4*I + p*x) +
                                          4*b1*p*(6*I + 9*p*x - 2*I*p*p*x*x) -
                                          12*b1*b1*x*(-5 + 8*I*p*x + 2*p*p*x*x))));
        
    }  else if ( l0 == 4 && l1 ==    2 ){
        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*(-12*b1*b1*b1*b1*b1*p*p + 12*b0*b1*b1*b1*(-p*p*p*p + 2*b1*p*p*(3 + 2*I*p*x) +
                                                                    2*b1*b1*(-5 - 10*I*p*x + 2*p*p*x*x)) +
                                         2*b0*b0*b0*b0*(-p*p*p*p + 32*b1*b1*b1*b1*b1*x*x*x*x*x*x + 16*b1*b1*b1*b1*x*x*x*x*(-15 + 4*I*p*x) +
                                                 2*b1*p*p*(-18 + 12*I*p*x + p*p*x*x) - 8*b1*b1*b1*x*x*(-45 + 40*I*p*x + 6*p*p*x*x) +
                                                 4*b1*b1*(-15 + 60*I*p*x + 36*p*p*x*x - 4*I*p*p*p*x*x*x)) -
                                         b0*b0*b1*(p*p*p*p*p*p + 4*b1*p*p*p*p*(-1 - 2*I*p*x) + 16*b1*b1*b1*b1*x*x*(-45 - 20*I*p*x + p*p*x*x) -
                                                  12*b1*b1*p*p*(9 - 8*I*p*x + 2*p*p*x*x) + 8*b1*b1*b1*(45 + 42*p*p*x*x + 4*I*p*p*p*x*x*x)) +
                                         2*b0*b0*b0*b1*(16*b1*b1*b1*b1*x*x*x*x*(-15 - 2*I*p*x) + p*p*p*p*(7 - 2*I*p*x) -
                                                    8*b1*p*p*(3 + 6*I*p*x + 2*p*p*x*x) + 16*b1*b1*b1*x*x*(45 - 10*I*p*x + 4*p*p*x*x) +
                                                        12*I*b1*b1*(15*I + 30*p*x + 4*I*p*p*x*x+ 4*p*p*p*x*x*x))));
    }  else if ( l0 == 4 && l1 ==    3 ){

        return -((fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*b1*(12*I*b1*b1*b1*b1*b1*p*p*p + 12*b0*b1*b1*b1*p*
                                      (I*p*p*p*p + 2*b1*p*p*(-5*I + 2*p*x) +
                                       b1*b1*(30*I - 30*p*x - 4*I*p*p*x*x)) +
                                      
                                      4*b0*b0*b0*b0*b0*(32*b1*b1*b1*b1*b1*x*x*x*x*x*x*x + 16*b1*b1*b1*b1*x*x*x*x*x*(-21 + 4*I*p*x) -
                                              3*p*p*p*(4*I + p*x) -
                                              24*b1*b1*b1*x*x*x*(-35 + 20*I*p*x + 2*p*p*x*x) +
                                              
                                              4*b1*b1*x*(-105 + 180*I*p*x + 60*p*p*x*x - 4*I*p*p*p*x*x*x) +
                                              2*b1*p*(-60*I - 90*p*x + 24*I*p*p*x*x + p*p*p*x*x*x)) +
                                      I*b0*b0*b1*(p*p*p*p*p*p*p + 4*b1*p*p*p*p*p*(-3 - 2*I*p*x) - 12*b1*b1*p*p*p*
                                                 (5 - 14*I*p*x + 2*p*p*x*x) +
                                                 8*b1*b1*b1*p*(75 - 90*I*p*x + 66*p*p*x*x +
                                                           4*I*p*p*p*x*x*x) +
                                                 16*b1*b1*b1*b1*x*(105*I - 135*p*x - 30*I*p*p*x*x + p*p*p*x*x*x)) -
                                      
                                      6*b0*b0*b0*b1*(p*p*p*p*p*(3*I + p*x) +
                                                 4*b1*p*p*p*(-10*I + 3*p*x - 2*I*p*p*x*x) +
                                                 16*b1*b1*b1*b1*x*x*x*(-35 - 15*I*p*x + p*p*x*x) +
                                                 8*b1*b1*b1*x*(105 + 30*I*p*x + 30*p*p*x*x +
                                                           4*I*p*p*p*x*x*x) -
                                                 12*b1*b1*p*(-5*I + 25*p*x - 4*I*p*p*x*x + 2*p*p*p*x*x*x)) +
                                      6*b0*b0*b0*b0*(I*p*p*p*p*p - 32*I*b1*b1*b1*b1*b1*x*x*x*x*x*(-7*I + p*x) +
                                              2*b1*p*p*p*(10*I + 11*p*x - I*p*p*x*x) +
                                              16*b1*b1*b1*b1*x*x*x*(70 - 5*I*p*x + 4*p*p*x*x) -
                                              4*b1*b1*p*(45*I + 28*I*p*p*x*x + 4*p*p*p*x*x*x) +
                                              
                                                     8*I*b1*b1*b1*x*(105*I + 75*p*x + 20*I*p*p*x*x + 6*p*p*p*x*x*x)))));
    }  else if ( l0 == 4 && l1 ==    4 ){

        return (fgg/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1)/(b0 + b1))*(b0*b0*b1*b1*(12*b1*b1*b1*b1*b1*b1*p*p*p*p - 12*b0*b1*b1*b1*b1*p*p*
                                           (-p*p*p*p + 2*b1*p*p*(7 + 2*I*p*x) +
                                            4*b1*b1*(-15 - 10*I*p*x + p*p*x*x)) +
                                           4*b0*b0*b0*b1*b1*(p*p*p*p*p*p*(-5 + 2*I*p*x) + 8*b1*p*p*p*p*(15 + 2*p*p*x*x) +
                                                        
                                                        16*b1*b1*b1*b1*x*x*(-210 - 210*I*p*x + 45*p*p*x*x + 2*I*p*p*p*x*x*x) -
                                                        
                                                        24*I*b1*b1*p*p*(-30*I + 25*p*x - 10*I*p*p*x*x + 2*p*p*p*x*x*x) +
                                                        
                                                        16*b1*b1*b1*(105 + 210*I*p*x + 45*p*p*x*x + 50*I*p*p*p*x*x*x -
                                                                 4*p*p*p*p*x*x*x*x)) +
                                           b0*b0*
                                           b1*b1*(p*p*p*p*p*p*p*p + 4*b1*p*p*p*p*p*p*(-5 - 2*I*p*x) +
                                                 12*b1*b1*p*p*p*p*(5 + 20*I*p*x - 2*p*p*x*x) +
                                                 
                                                 16*b1*b1*b1*p*p*(15 - 120*I*p*x + 45*p*p*x*x + 2*I*p*p*p*x*x*x) +
                                                 
                                                 16*b1*b1*b1*b1*(105 + 420*I*p*x - 270*p*p*x*x - 40*I*p*p*p*x*x*x +
                                                          p*p*p*p*x*x*x*x)) +
                                           4*b0*b0*b0*b0*b0*b0*(3*p*p*p*p+ 64*b1*b1*b1*b1*b1*b1*x*x*x*x*x*x*x*x + 128*I*b1*b1*b1*b1*b1*x*x*x*x*x*x*(7*I + p*x) -
                                                   12*b1*p*p*(-15 + 10*I*p*x + p*p*x*x) -
                                                   96*b1*b1*b1*b1*x*x*x*x*(-35 + 14*I*p*x + p*p*x*x)) +
                                                   
                                                   16*b1*b1*b1*x*x*(-210 + 210*I*p*x + 45*p*p*x*x - 2*I*p*p*p*x*x*x) +
                                                   
                                                   4*b1*b1*(105 - 420*I*p*x - 270*p*p*x*x + 40*I*p*p*p*x*x*x +
                                                           p*p*p*p*x*x*x*x)) -
                                           12*b0*b0*b0*b0*b1*(-p*p*p*p*p*p + 32*b1*b1*b1*b1*b1*x*x*x*x*(-35 - 14*I*p*x + p*p*x*x) +
                                                       b1*p*p*p*p*(-5 + 20*I*p*x + 2*p*p*x*x) +
                                                       8*b1*b1*p*p*(30 - 25*I*p*x + 10*p*p*x*x -
                                                                   2*I*p*p*p*x*x*x) +
                                                       16*b1*b1*b1*b1*x*x*(210 + 70*I*p*x + 25*p*p*x*x + 4*I*p*p*p*x*x*x) -
                                                       24*b1*b1*b1*(35 + 50*p*p*x*x + 2*p*p*p*p*x*x*x*x)) +
                                           8*b0*b0*b0*b0*b0*
                                           b1*(3*p*p*p*p*(-7 + 2*I*p*x) - 64*I*b1*b1*b1*b1*b1*x*x*x*x*x*x*(-7*I + p*x) +
                                               32*b1*b1*b1*b1*x*x*x*x*(105 + 4*p*p*x*x) +
                                               2*b1*p*p*(15 + 120*I*p*x + 45*p*p*x*x -
                                                         2*I*p*p*p*x*x*x +
                                               24*I*b1*b1*b1*x*x*(210*I + 70*p*x + 25*I*p*p*x*x + 4*p*p*p*x*x*x) -
                                               
                                               8*b1*b1*(-105 + 210*I*p*x - 45*p*p*x*x + 50*I*p*p*p*x*x*x +
                                                        4*p*p*p*p*x*x*x*x))));
    }
        
        return 0;
}


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

DCOMPLEX FGS( double p , struct general_index * pa ){
    
    // Integrate[[   D[Gaus( beta, (x-o) ),{o,l}] Exp[ i p x ] Sinc ( x/d- n )
    // one gaussian, one sinc, fourier transformed.
    struct general_2index g2;
    double beta,realpart,imagepart;
    g2.periodic = 0;
    g2.i[0].d = pa->d;
    g2.i[0].bra.basis = SincBasisElement;
    g2.i[0].ket.basis = nullBasisElement;
    g2.i[1].bra.basis = DiracDelta;
    g2.i[1].ket.basis = nullBasisElement;
    g2.fl->fn = nullFunction;
    g2.gaussianAccelerationFlag = 0;
    g2.fl->param[0] = 1;
    g2.point = 1;
    g2.momentumShift = p ;
    g2.i[0].pointer = 1;
    g2.momentumShift = 0;

#if 0
    if ( pa->bra.basis == GaussianBasisElement && pa->ket.basis == SincBasisElement ){
        
        beta = pa->bra.length;
        g2.powSpace = g2.i[0].bra.index;
        g2.i[1].bra.origin = g2.i[0].bra.origin;
        g2.i[0].bra = pa->ket;


    }
    else if ( pa->ket.basis == GaussianBasisElement && pa->bra.basis == SincBasisElement ){
        beta = pa->ket.length;
        g2.powSpace = g2.i[0].ket.index;
        g2.i[1].bra.origin = g2.i[0].ket.origin;
        g2.i[0].bra = pa->bra;


    }else {
        printf("wrong mail stop\n");
        exit(0);
    }
#else
    beta = pa->b0;
    g2.powSpace = pa->l0;
    g2.i[1].bra.origin = pa->x0;
    g2.i[0].bra.index = pa->n;
#endif
    g2.realFlag = 1;
    realpart =  collective(sqrt(beta), &g2);
    g2.realFlag = 0;
    imagepart = collective(sqrt(beta), &g2);
    return realpart + I * imagepart;
};



DCOMPLEX FSS ( double p , struct general_index * pa ){
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
    }else if ( pt == 2 ){
        //return ei((n+m)*(1./2.*p*d+pi)) *( delta(n-m) -fabs( p * d / 2./pi ) * Sinc(1., (p*d/2./pi) *(n-m)));
        return ei((n+m)*(1./2.*p*d+pi)) *( Sinc(1. ,(n-m)) -fabs( p * d / 2./pi ) * Sinc(1., (p*d/2./pi) *(n-m)));

    }else
    {
        printf("pointy");
        exit(0);
    }
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
    if ( pa->bra.basis == SincBasisElement && pa->ket.basis == SincBasisElement){
        pa->n = pa->bra.index;
        pa->m = pa->ket.index;
        pa->d = pa->bra.length;
        pa->pointer = 2;
        if ( (pa->bra.length - pa->ket.length ) > 1e-3 ){
            printf("parties correspondence\n");
            exit(0);
        }
        return FSS(p, pa);
    }else if ( pa->bra.basis == DiracDelta || pa->ket.basis == DiracDelta){
        if ( pa->bra.basis == DiracDelta && pa->ket.basis == nullBasisElement){
            pa->x0 = pa->bra.origin;
        }
        else  if ( pa->ket.basis == DiracDelta && pa->bra.basis == nullBasisElement){
            pa->x0 = pa->ket.origin;
        }else {
            printf("parity correspondence\n");
            exit(0);
        }
        return FDD(p, pa);
    }else if ( pa->bra.basis == GaussianBasisElement || pa->ket.basis == GaussianBasisElement){
        pa->b0 = pa->bra.length;
        pa->b1 = pa->ket.length;
        pa->l0 = pa->bra.index;
        pa->l1 = pa->ket.index;
        pa->x0 = pa->bra.origin;
        pa->x1 = pa->ket.origin;
        return FGG(p, pa);
    }else     if ( pa->bra.basis == SincBasisElement && pa->ket.basis == nullBasisElement){
        pa->n = pa->bra.index;
        pa->d = pa->bra.length;
        pa->pointer = 0;
        return FSS(p, pa);
    }else     if ( pa->bra.basis == nullBasisElement && pa->ket.basis == SincBasisElement){
        pa->m = pa->ket.index;
        pa->d = pa->ket.length;
        pa->pointer = 1;
        return FSS(p, pa);
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




double periodicSd2S ( INT_TYPE arg, INT_TYPE N ){
    
#ifdef GSL_LIB
    
    
    if ( arg == 0 ){
        return Sd2S(0)-4.0*gsl_sf_zeta(2.)/sqr(2.*N)+4.0*gsl_sf_hzeta(2.,0.500)/sqr(2.*N);
    }
    
    return -sign(arg)/2./sqr(N)*(gsl_sf_hzeta(2., 1-(1.*labs(arg))/(2.*N))+gsl_sf_hzeta(2., (1.*labs(arg))/(2.*N))-gsl_sf_hzeta(2., 0.5+(1.*labs(arg))/(2.*N))-gsl_sf_hzeta(2.,0.5- (1.*labs(arg))/(2.*N)));
#else
    
    if ( N == 9 ){
        double cal9[] = {-3.249252,    1.957614,    -0.451818,    0.162463,    -0.043633,    -0.043633,    0.162463,    -0.451818,    1.957614};
        return cal9[abs(arg)];
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



double BoB (struct basisElement b1, struct basisElement b2){
    
    if ( b1.basis == SincBasisElement && b2.basis == SincBasisElement ){
        if ( b1.periodic == 1 ){
            return periodicSS ( b1.length,b1.length*(b1.index) + b1.origin, b1.auxIndex, b2.length, b2.length*( b2.index ) + b2.origin,b2.auxIndex);
        }
        else{
            return SS ( b1.length,b1.length*(b1.index)+ b1.origin, b2.length, b2.length*( b2.index ) + b2.origin);
        }
    }else if ( b1.basis == GaussianBasisElement && b2.basis == GaussianBasisElement ){
        
        struct general_index pa ;
        pa.b0 = b1.length;
        pa.l0 = b1.index;
        pa.x0 = b1.origin;
        pa.b1 = b2.length;
        pa.l1 = b2.index;
        pa.x1 = b2.origin;
        
        return FGG(0,&pa);
        
    }else if ( b1.basis == GaussianBasisElement && b2.basis == SincBasisElement ){
        
            struct general_index pa ;
            pa.b0 = b1.length;
            pa.l0 = b1.index;
            pa.x0 = b1.origin;
            pa.d = b2.length;
            pa.n = b2.index;
        
            return FGS(0,&pa);
        
    }else if ( b1.basis == SincBasisElement && b2.basis == GaussianBasisElement ){
        
        struct general_index pa ;
        pa.b0 = b2.length;
        pa.l0 = b2.index;
        pa.x0 = b2.origin;
        pa.d = b1.length;
        pa.n = b1.index;
        
        return FGS(0,&pa);
    }

    printf("over rails\n");
    exit(0);
    return 0;
}

double BdB (struct basisElement b1, struct basisElement b2){
    if ( b1.basis == SincBasisElement && b2.basis == SincBasisElement ){
        if ( b1.periodic == 0 ){
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
                
                return periodicSd2S(arg,b1.auxIndex) * sqrt(b2.length/b1.length);
            }else {
                arg /= b2.length;
                
                return periodicSd2S(arg,b2.auxIndex) * sqrt(b1.length/b2.length);
            }
        }
    }else if ( b1.basis == GaussianBasisElement && b2.basis == GaussianBasisElement ){
        
        struct general_index pa ;
        pa.b0 = b1.length;
        pa.l0 = b1.index+1;
        pa.x0 = b1.origin;
        pa.b1 = b2.length;
        pa.l1 = b2.index;
        pa.x1 = b2.origin;
        
        return FGG(0,&pa);
        
    }else if ( b1.basis == GaussianBasisElement && b2.basis == SincBasisElement ){
        
        struct general_index pa ;
        pa.b0 = b1.length;
        pa.l0 = b1.index+1;
        pa.x0 = b1.origin;
        pa.d = b2.length;
        pa.n = b2.index;
        
        return FGS(0,&pa);
        
    }else if ( b1.basis == SincBasisElement && b2.basis == GaussianBasisElement ){
        
        struct general_index pa ;
        pa.b0 = b2.length;
        pa.l0 = b2.index+1;
        pa.x0 = b2.origin;
        pa.d = b1.length;
        pa.n = b1.index;
        
        return FGS(0,&pa);
    }
    
    printf("over rails\n");
    exit(0);
    return 0;
}


double Bd2B (struct basisElement b1, struct basisElement b2){
    if ( b1.basis == SincBasisElement && b2.basis == SincBasisElement ){
        if ( b1.periodic == 0 ){
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
            
            double arg = b1.length*b1.index + b1.origin - (b2.length*b2.index + b2.origin);
            
            if ( b1.length > b2.length ){
                arg /= b1.length;
                
                return periodicSd2S(arg,b1.auxIndex) * sqrt(b2.length/b1.length);
            }else {
                arg /= b2.length;
                
                return periodicSd2S(arg,b2.auxIndex) * sqrt(b1.length/b2.length);
            }
        }
    }else if ( b1.basis == GaussianBasisElement && b2.basis == GaussianBasisElement ){
        
        struct general_index pa ;
        pa.b0 = b1.length;
        pa.l0 = b1.index+2;
        pa.x0 = b1.origin;
        pa.b1 = b2.length;
        pa.l1 = b2.index;
        pa.x1 = b2.origin;
        
        return FGG(0,&pa);
        
    }else if ( b1.basis == GaussianBasisElement && b2.basis == SincBasisElement ){
        
        struct general_index pa ;
        pa.b0 = b1.length;
        pa.l0 = b1.index+2;
        pa.x0 = b1.origin;
        pa.d = b2.length;
        pa.n = b2.index;
        
        return FGS(0,&pa);
        
    }else if ( b1.basis == SincBasisElement && b2.basis == GaussianBasisElement ){
        
        struct general_index pa ;
        pa.b0 = b2.length;
        pa.l0 = b2.index+2;
        pa.x0 = b2.origin;
        pa.d = b1.length;
        pa.n = b1.index;
        
        return FGS(0,&pa);
    }
    
    printf("over rails\n");
    exit(0);
    return 0;
}


DCOMPLEX poly( double k ,double beta, INT_TYPE powSpace ){
    if ( ! powSpace )
        return 1.;
    else if ( powSpace == 1 )
        return 0.50*(-I * k)/ beta;
    else if ( powSpace == 2 )
        return -0.250*(k*k - 2. *beta )/ beta/beta;
    else if ( powSpace == 3 )
        return 0.5*0.250*(I * k)*(-6. *beta + k*k )/ beta/beta/beta;
    else if ( powSpace == 4 )
        return 0.250*0.250*(k*k*k*k -12. * k*k * beta  + 12. * beta * beta ) / beta/beta/beta/beta;
    else if ( powSpace == 5 )
        return 0.5*0.250*0.250*(-I * k)*(k*k*k*k -20. * k*k * beta  + 60. * beta * beta ) / beta/beta/beta/beta/beta;
    else if ( powSpace == 6 )
        return -0.250*0.250*0.250*( k*k*k*k*k*k - 30. * k*k  *k*k * beta + 180. * k*k*beta*beta-120 * beta*beta*beta ) / beta/beta/beta/beta/beta/beta;
    else if ( powSpace == 7 )
        return 0.5*0.250*0.250*0.250*(I * k)*( k*k*k *k*k*k - 42. * k*k*k*k * beta + 420. * k*k*beta*beta-840. * beta*beta*beta ) / beta/beta/beta/beta/beta/beta/beta;
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

    double alpha = pa->alpha ;
    DCOMPLEX act = 1.;
    if ( pa->body == 1 ){
        if ( pa->i[1].action > 0 )
            act *= I * k ;
        if ( pa->i[1].action > 1 )
            act *= I * k ;
        if ( pa->i[1].action > 2 ){
            printf("type more\n");
            fflush(stdout);
            exit(0);
        }
    }
    if ( pa->realFlag == 0 )
        act *= -I ;
    
    return  exp(-sqr(k /2./ alpha))/alpha * creal(act* poly(k,alpha*alpha,pa->powSpace) *FB(-k,&pa->i[0])*FB(k,&pa->i[1]))/(2.*pi);

}

void gaussianSincFunc(void * arg,size_t n,const double * x,double * y)
{
    struct general_2index *af = (struct general_2index *) arg;
    size_t i ;
    for ( i=0;i<n;i++)
    {
        y[i] = gaussianSinc(x[i],af);
       // printf("%d %f %f\n", i, x[i],y[i]);
    }
    
}

double mcGS ( double x [], size_t dim , void * params ){
    struct general_2index* ga = (struct general_2index*)(params);
    ga[0].alpha = x[0];
    ga[1].alpha = x[0];
    ga[2].alpha = x[0];
    
    return gaussianSinc(x[1], ga)*gaussianSinc(x[2], ga+1)*gaussianSinc(x[3], ga+2);
    
    
}


double c10c00 ( double a, double bd, double xd, double b1, double x1, double b2, double x2, double b3 ,double x3 ){
    return (2*sqrt(bd)*(b1*(b2 + b3)*(x1 - xd) + a*a*(b2*x2 + b3*x3 + b1*(x1 - xd) -
                                                      b2*xd - b3*xd)))/((b2 + b3)*(b1 + bd) + a*a*(b1 + b2 + b3 + bd));
};

double c11c00 ( double a, double bd1, double xd1, double bd2, double xd2, double b2, double x2, double b3 ,double x3 ){
    
    return -((2*sqrt(bd1*bd2)*
              (sqr(b2 + b3)*(-bd2 + bd1*(-1 + 2*bd2*sqr(xd1 - xd2))) -
               a*a*(b2 + b3)*(b3 + 2*(bd1 + bd2) +
                              b2*(1 + 2*(xd1 - xd2)*(bd1*(x2 - xd1) + bd2*(-x2 + xd2))) +
                              2*(xd1 - xd2)*(2*bd1*bd2*(-xd1 + xd2) +
                                             b3*(bd1*(x3 - xd1) + bd2*(-x3 + xd2)))) -
               a*a*a*a*(b3 + bd1 + bd2 + 2*b2*b2*(x2 - xd1)*(x2 - xd2) +
                        2*(b3*x3 + bd1*xd1 - (b3 + bd1)*xd2)*(b3*(x3 - xd1) +
                                                              bd2*(-xd1 + xd2)) + b2*(1 + 2*(xd1 - xd2)*(bd1*(x2 - xd1) +
                                                                                                         bd2*(-x2 + xd2)) - 2*b3*(-2*xd1*xd2 + x3*(xd1 + xd2) +
                                                                                                                                  x2*(-2*x3 + xd1 + xd2))))))/
             sqr((b2 + b3)*(bd1 + bd2) + a*a*(b2 + b3 + bd1 + bd2)));
    
    
    
    
}


double c10c10 (  double a, double bd1, double xd1,  double b2, double x2, double bd2, double xd2,double b3 ,double x3 ){
    double result= (2*sqrt(bd1*bd2)*(2*b2*b3*(b2 + bd1)*(b3 + bd2)*(x2 - xd1)*(x3 - xd2) +
                                     a*a*(2*b2*b2*(x2 - xd1)*(b3*(x2 + x3 - 2*xd2) + bd2*(x2 - xd2)) +
                                          b2*(bd2*(1 + 2*bd1*(x2 - xd1)*(xd1 - xd2)) + 2*b3*b3*(x2 + x3 - 2*xd1)*
                                              (x3 - xd2) + b3*(1 + 2*bd1*(x2 - xd1)*(x3 + xd1 - 2*xd2) +
                                                               2*bd2*(x3 - xd2)*(x2 - 2*xd1 + xd2))) +
                                          bd1*(bd2 + 2*b3*b3*(x3 - xd1)*(x3 - xd2) +
                                               b3*(1 + 2*bd2*(x3 - xd2)*(-xd1 + xd2)))) +
                                     a*a*a*a*(bd1 + bd2 - 2*bd1*bd2*xd1*xd1 + 2*b2*b2*(x2 - xd1)*(x2 - xd2) +
                                              2*b3*b3*(x3 - xd1)*(x3 - xd2) + 4*bd1*bd2*xd1*xd2 - 2*bd1*bd2*xd2*xd2 +
                                              b3*(1 + 2*bd1*(x3 - xd1)*(xd1 - xd2) + 2*bd2*(x3 - xd2)*(-xd1 + xd2)) +
                                              b2*(1 - 2*bd2*x2*xd1 + 2*bd1*(x2 - xd1)*(xd1 - xd2) + 2*bd2*x2*xd2 +
                                                  2*bd2*xd1*xd2 - 2*bd2*xd2*xd2 + b3*(x2*(4*x3 - 2*(xd1 + xd2)) -
                                                                                      2*(-2*xd1*xd2 + x3*(xd1 + xd2)))))))/
    sqr((b2 + bd1)*(b3 + bd2) + a*a*(b2 + b3 + bd1 + bd2));
    return result;
    
}

double c11c10 ( double a , double bd1, double  xd1,double  bd2,double  xd2, double bd3,double xd3,double  b1,double x1){
    return (sqrt(bd1*bd2*bd3)*(4*(a*a + b1 + bd3)*((bd1 + bd2)*(b1 + bd3) +
                                                   a*a*(b1 + bd1 + bd2 + bd3))*(b1*(bd1 + bd2)*(x1 - xd3) +
                                                                                a*a*(bd1*xd1 + bd2*xd2 + b1*(x1 - xd3) - bd1*xd3 - bd2*xd3)) -
                               (1/bd1)*(2*a*a*((bd1 + bd2)*(b1 + bd3) + a*a*(b1 + bd1 + bd2 + bd3))*
                                        (2*bd1*bd2*(b1 + bd3)*(xd1 - xd2) - 2*a*a*bd1*(b1*(x1 - xd1) - bd2*xd1 -
                                                                                       bd3*xd1 + bd2*xd2 + bd3*xd3))) -
                               (1/bd2)*(2*a*a*((bd1 + bd2)*(b1 + bd3) + a*a*(b1 + bd1 + bd2 + bd3))*
                                        (-2*bd1*bd2*(b1 + bd3)*(xd1 - xd2) - 2*a*a*bd2*(bd1*xd1 + b1*(x1 - xd2) -
                                                                                        bd1*xd2 - bd3*xd2 + bd3*xd3))) - (1/(bd1*bd2))*
                               (2*(2*bd1*bd2*(b1 + bd3)*(xd1 - xd2) - 2*a*a*bd1*(b1*(x1 - xd1) -
                                                                                 bd2*xd1 - bd3*xd1 + bd2*xd2 + bd3*xd3))*
                                (-2*bd1*bd2*(b1 + bd3)*(xd1 - xd2) - 2*a*a*bd2*(bd1*xd1 + b1*(x1 - xd2) -
                                                                                bd1*xd2 - bd3*xd2 + bd3*xd3))*((-b1)*(bd1 + bd2)*(x1 - xd3) +
                                                                                                               a*a*((-bd1)*xd1 - bd2*xd2 + bd1*xd3 + bd2*xd3 + b1*(-x1 + xd3))))))/
    cube((bd1 + bd2)*(b1 + bd3) + a*a*(b1 + bd1 + bd2 + bd3));
    
}

double c11c11 (double a, double bd1,double xd1,  double bd2, double xd2, double bd3, double xd3,double bd4, double xd4){
    return (sqrt(bd1*bd2*bd3*bd4)*
            (8*a*a*a*a*sqr((bd1 + bd2)*(bd3 + bd4) + a*a*(bd1 + bd2 + bd3 + bd4)) +
             4*(a*a + bd1 + bd2)*(a*a + bd3 + bd4)*
             sqr((bd1 + bd2)*(bd3 + bd4) + a*a*(bd1 + bd2 + bd3 + bd4)) +
             (1/bd4)*(4*a*a*((bd1 + bd2)*(bd3 + bd4) + a*a*(bd1 + bd2 + bd3 + bd4))*
                      (-2*(bd1 + bd2)*bd3*bd4*(xd3 - xd4) - 2*a*a*bd4*(bd2*xd2 + bd3*xd3 +
                                                                       bd1*(xd1 - xd4) - bd2*xd4 - bd3*xd4))*(bd2*(bd3 + bd4)*(xd1 - xd2) +
                                                                                                              a*a*(bd3*xd1 + bd4*xd1 + bd2*(xd1 - xd2) - bd3*xd3 - bd4*xd4))) -
             8*(a*a + bd3 + bd4)*((bd1 + bd2)*(bd3 + bd4) + a*a*(bd1 + bd2 + bd3 + bd4))*
             ((bd1 + bd2)*bd3*(xd3 - xd4) + a*a*(bd2*xd2 + bd3*xd3 + bd1*(xd1 - xd4) -
                                                 bd2*xd4 - bd3*xd4))*((bd1 + bd2)*bd4*(xd3 - xd4) +
                                                                      a*a*((-bd2)*xd2 + bd2*xd3 + bd4*xd3 + bd1*(-xd1 + xd3) - bd4*xd4)) +
             8*a*a*((bd1 + bd2)*(bd3 + bd4) + a*a*(bd1 + bd2 + bd3 + bd4))*
             (bd2*(bd3 + bd4)*(xd1 - xd2) + a*a*(bd3*xd1 + bd4*xd1 + bd2*(xd1 - xd2) -
                                                 bd3*xd3 - bd4*xd4))*((bd1 + bd2)*bd4*(xd3 - xd4) +
                                                                      a*a*((-bd2)*xd2 + bd2*xd3 + bd4*xd3 + bd1*(-xd1 + xd3) - bd4*xd4)) +
             (1/(bd2*bd4))*(2*a*a*((bd1 + bd2)*(bd3 + bd4) +
                                   a*a*(bd1 + bd2 + bd3 + bd4))*(-2*(bd1 + bd2)*bd3*bd4*(xd3 - xd4) -
                                                                 2*a*a*bd4*(bd2*xd2 + bd3*xd3 + bd1*(xd1 - xd4) - bd2*xd4 - bd3*xd4))*
                            (-2*bd1*bd2*(bd3 + bd4)*(xd1 - xd2) - 2*a*a*bd2*(bd1*(xd1 - xd2) -
                                                                             bd3*xd2 - bd4*xd2 + bd3*xd3 + bd4*xd4))) +
             (1/bd2)*(4*(a*a + bd1 + bd2)*((bd1 + bd2)*(bd3 + bd4) +
                                           a*a*(bd1 + bd2 + bd3 + bd4))*(bd2*(bd3 + bd4)*(xd1 - xd2) +
                                                                         a*a*(bd3*xd1 + bd4*xd1 + bd2*(xd1 - xd2) - bd3*xd3 - bd4*xd4))*
                      (-2*bd1*bd2*(bd3 + bd4)*(xd1 - xd2) - 2*a*a*bd2*(bd1*(xd1 - xd2) -
                                                                       bd3*xd2 - bd4*xd2 + bd3*xd3 + bd4*xd4))) +
             (1/bd2)*(4*a*a*((bd1 + bd2)*(bd3 + bd4) + a*a*(bd1 + bd2 + bd3 + bd4))*
                      ((bd1 + bd2)*bd4*(xd3 - xd4) + a*a*((-bd2)*xd2 + bd2*xd3 + bd4*xd3 +
                                                          bd1*(-xd1 + xd3) - bd4*xd4))*(-2*bd1*bd2*(bd3 + bd4)*(xd1 - xd2) -
                                                                                        2*a*a*bd2*(bd1*(xd1 - xd2) - bd3*xd2 - bd4*xd2 + bd3*xd3 + bd4*xd4))) +
             (1/(bd2*bd4))*(4*(-2*(bd1 + bd2)*bd3*bd4*(xd3 - xd4) -
                               2*a*a*bd4*(bd2*xd2 + bd3*xd3 + bd1*(xd1 - xd4) - bd2*xd4 - bd3*xd4))*
                            (bd2*(bd3 + bd4)*(xd1 - xd2) + a*a*(bd3*xd1 + bd4*xd1 + bd2*(xd1 - xd2) -
                                                                bd3*xd3 - bd4*xd4))*((bd1 + bd2)*bd4*(xd3 - xd4) +
                                                                                     a*a*((-bd2)*xd2 + bd2*xd3 + bd4*xd3 + bd1*(-xd1 + xd3) - bd4*xd4))*
                            (-2*bd1*bd2*(bd3 + bd4)*(xd1 - xd2) - 2*a*a*bd2*(bd1*(xd1 - xd2) -
                                                                             bd3*xd2 - bd4*xd2 + bd3*xd3 + bd4*xd4)))))/
    sqr(sqr((bd1 + bd2)*(bd3 + bd4) + a*a*(bd1 + bd2 + bd3 + bd4)));
}


double c20c00 (double a, double b0, double x0, double b1, double x1, double b2, double x2, double b3, double x3){
    x1-= x0;
    x2-= x0;
    x3-= x0;

    return 1./b0*(2*a*b0*sqrt((b0 + b1)*(b2 + b3)*(1/a/a +
                                             4*(1/(b0 + b1) + 1/(b2 + b3))))*
            (-(4*a*a*(b2 + b3) + b1*(4*a*a + b2 + b3))*(4*a*a*(b2 + b3) + b1*(4*a*a + b2 + b3)) +
             b0*(2*b1*b1*(4*a*a + b2 + b3)*(4*a*a + b2 + b3)*x1*x1 - b1*(4*a*a + b2 + b3)*
                 (4*a*a + b2 + b3 + 16*a*a*b2*x1*x2 +
                  16*a*a*b3*x1*x3) +
                 4*a*a*(-4*a*a*b2 - 4*a*a*b3 + b2*b2*(-1 + 8*a*a*x2*x2) +
                        2*b2*b3*(-1 + 8*a*a*x2*x3) +
                        b3*b3*(-1 + 8*a*a*x3*x3)))))/
    pow(4*a*a*(b2 + b3) + b0*(4*a*a + b2 + b3) + b1*(4*a*a + b2 + b3),5./2.);
}


double c22c00 (double a, double b0, double x0, double b1, double x1, double b2, double x2, double b3, double x3){
    x1-= x0;
    x2-= x0;
    x3-= x0;


    return (4*b0*b1*sqrt(b0 + b1)*sqrt(b2 + b3)*
 sqrt((b0 + b1)*(b2 + b3)*(4*(1/(b0 + b1) + 1/(b2 + b3)) + 1/(a*a)))*a*
 (-4*(4*(b2 + b3)*a*a*a + b1*a*(b2 + b3 + 4*a*a))*(4*(b2 + b3)*a*a*a + b1*a*(b2 + b3 + 4*a*a))*(-4*(b2 + b3)*(b2 + b3)*a*a +
                                                                        b1*(-4*b2*a*a - 4*b3*a*a + b2*b2*(-1 + 8*x1*x1*a*a + 16*x1*x2*a*a +
                                                                                                                       8*x2*x2*a*a) + b3*b3*(-1 + 8*x1*x1*a*a + 16*x1*x3*a*a + 8*x3*x3*a*a) +
                                                                            2*b2*b3*(-1 + 8*x1*x1*a*a + 8*x2*x3*a*a + 8*x1*(x2 + x3)*a*a))) +
  b0*b0*b0*(b2 + b3 + 4*a*a)*(b2 + b3 + 4*a*a)*(4*b1*b1*b1*x1*x1*x1*x1*(b2 + b3 + 4*a*a)*(b2 + b3 + 4*a*a) -
                                   4*b1*b1*x1*x1*(b2 + b3 + 4*a*a)*(12*a*a + b2*(3 + 8*x1*x2*a*a) +
                                                                         b3*(3 + 8*x1*x3*a*a)) + 4*a*a*(b2*b2*(1 - 8*x2*x2*a*a) +
                                                                                                                      2*b2*(b3 + 2*a*a - 8*b3*x2*x3*a*a) + b3*(b3 + 4*a*a - 8*b3*x3*x3*a*a)) +
                                   b1*(48*a*a*a*a - 8*b3*a*a*(-3 + 4*x1*x1*a*a - 24*x1*x3*a*a) +
                                       b2*b2*(3 + 48*x1*x2*a*a + 8*x1*x1*a*a*(-1 + 8*x2*x2*a*a)) +
                                       b3*b3*(3 + 48*x1*x3*a*a + 8*x1*x1*a*a*(-1 + 8*x3*x3*a*a)) +
                                       2*b2*(4*a*a*(3 - 4*x1*x1*a*a + 24*x1*x2*a*a) +
                                             b3*(3 + 24*x1*(x2 + x3)*a*a + 8*x1*x1*a*a*(-1 + 8*x2*x3*a*a))))) +
  2*b0*b0*(b2 + b3 + 4*a*a)*(2*b1*b1*b1*x1*x1*(b2 + b3 + 4*a*a)*(b2 + b3 + 4*a*a)*
                                   (-12*a*a + b2*(-3 + 8*x1*x1*a*a + 8*x1*x2*a*a) +
                                    b3*(-3 + 8*x1*x1*a*a + 8*x1*x3*a*a)) - 8*(b2 + b3)*a*a*a*a*
                                   (b2*b2*(-3 + 16*x2*x2*a*a) + 2*b2*(-6*a*a + b3*(-3 + 16*x2*x3*a*a)) +
                                    b3*(-12*a*a + b3*(-3 + 16*x3*x3*a*a))) - b1*b1*(b2 + b3 + 4*a*a)*
                                   (-48*a*a*a*a + 8*b3*a*a*(-3 + 28*x1*x1*a*a) + b2*b2*(-3 + 128*x1*x1*x1*x2*a*a*a*a +
                                                                                                       8*x1*x1*a*a*(7 + 16*x2*x2*a*a)) + b3*b3*(-3 + 128*x1*x1*x1*x3*a*a*a*a +
                                                                                                                                                           8*x1*x1*a*a*(7 + 16*x3*x3*a*a)) + 2*b2*(4*a*a*(-3 + 28*x1*x1*a*a) +
                                                                                                                                                                                                               b3*(-3 + 64*x1*x1*x1*(x2 + x3)*a*a*a*a + 8*x1*x1*a*a*(7 + 16*x2*x3*a*a)))) +
                                   2*b1*a*a*(b2*b2*b2*(9 + 24*x2*x2*a*a + 16*x1*x2*a*a*(5 + 8*x2*x2*a*a) +
                                                          8*x1*x1*a*a*(-3 + 16*x2*x2*a*a)) +
                                                    b2*b2*(8*a*a*(9 - 12*x1*x1*a*a + 40*x1*x2*a*a + 12*x2*x2*a*a) +
                                                          b3*(3*(9 + 8*x2*x2*a*a + 16*x2*x3*a*a) + 8*x1*x1*a*a*(-9 + 16*x2*x2*a*a +
                                                                                                                                   32*x2*x3*a*a) + 16*x1*a*a*(10*x2 + 5*x3 + 24*x2*x2*x3*a*a))) +
                                                    b3*(144*a*a*a*a + 8*b3*a*a*(9 - 12*x1*x1*a*a + 40*x1*x3*a*a + 12*x3*x3*a*a) +
                                                        b3*b3*(9 + 24*x3*x3*a*a + 16*x1*x3*a*a*(5 + 8*x3*x3*a*a) +
                                                              8*x1*x1*a*a*(-3 + 16*x3*x3*a*a))) +
                                                    b2*(144*a*a*a*a - 16*b3*(-9*a*a + 4*(3*x1*x1 - 3*x2*x3 - 5*x1*(x2 + x3))*a*a*a*a) +
                                                        b3*b3*(3*(9 + 16*x2*x3*a*a + 8*x3*x3*a*a) + 8*x1*x1*a*a*(-9 + 32*x2*x3*a*a +
                                                                                                                                   16*x3*x3*a*a) + 16*x1*a*a*(5*x2 + 10*x3 + 24*x2*x3*x3*a*a))))) +
  b0*(-64*(b2 + b3)*(b2 + b3)*a*a*a*a*a*a*(b2*b2*(-3 + 8*x2*x2*a*a) +
                                  2*b2*(-6*a*a + b3*(-3 + 8*x2*x3*a*a)) +
                                  b3*(-12*a*a + b3*(-3 + 8*x3*x3*a*a))) + b1*b1*b1*(b2 + b3 + 4*a*a)*(b2 + b3 + 4*a*a)*
      (48*a*a*a*a - 8*b3*a*a*(-3 + 28*x1*x1*a*a + 24*x1*x3*a*a) +
       b2*b2*(3 - 48*x1*x2*a*a + 64*x1*x1*x1*x1*a*a*a*a + 128*x1*x1*x1*x2*a*a*a*a +
             8*x1*x1*a*a*(-7 + 8*x2*x2*a*a)) + b3*b3*(3 - 48*x1*x3*a*a + 64*x1*x1*x1*x1*a*a*a*a +
                                                                 128*x1*x1*x1*x3*a*a*a*a + 8*x1*x1*a*a*(-7 + 8*x3*x3*a*a)) +
       2*b2*(-4*a*a*(-3 + 28*x1*x1*a*a + 24*x1*x2*a*a) +
             b3*(3 - 24*x1*(x2 + x3)*a*a + 64*x1*x1*x1*x1*a*a*a*a + 64*x1*x1*x1*(x2 + x3)*a*a*a*a +
                 8*x1*x1*a*a*(-7 + 8*x2*x3*a*a)))) - 4*b1*b1*a*a*(b2 + b3 + 4*a*a)*
      (b2*b2*b2*(128*x1*x1*x1*x2*a*a*a*a - 3*(3 + 8*x2*x2*a*a) + 16*x1*x1*a*a*
             (5 + 16*x2*x2*a*a) + 32*x1*(x2*a*a + 4*x2*x2*x2*a*a*a*a)) +
       b2*b2*(8*a*a*(-9 + 40*x1*x1*a*a + 16*x1*x2*a*a - 12*x2*x2*a*a) +
             b3*(128*x1*x1*x1*(2*x2 + x3)*a*a*a*a - 3*(9 + 8*x2*x2*a*a + 16*x2*x3*a*a) +
                 16*x1*x1*a*a*(15 + 16*x2*x2*a*a + 32*x2*x3*a*a) +
                 32*x1*a*a*(2*x2 + x3 + 12*x2*x2*x3*a*a))) +
       b2*(-144*a*a*a*a + 16*b3*(-9*a*a + 4*(10*x1*x1 - 3*x2*x3 + 2*x1*(x2 + x3))*a*a*a*a) +
           b3*b3*(128*x1*x1*x1*(x2 + 2*x3)*a*a*a*a - 3*(9 + 16*x2*x3*a*a + 8*x3*x3*a*a) +
                 16*x1*x1*a*a*(15 + 32*x2*x3*a*a + 16*x3*x3*a*a) +
                 32*x1*a*a*(x2 + 2*x3 + 12*x2*x3*x3*a*a))) +
       b3*(-144*a*a*a*a + 8*b3*a*a*(-9 + 40*x1*x1*a*a + 16*x1*x3*a*a - 12*x3*x3*a*a) +
           b3*b3*(128*x1*x1*x1*x3*a*a*a*a - 3*(3 + 8*x3*x3*a*a) + 16*x1*x1*a*a*
                 (5 + 16*x3*x3*a*a) + 32*x1*(x3*a*a + 4*x3*x3*x3*a*a*a*a)))) +
      16*b1*a*a*a*a*(b2*b2*b2*b2*(9 + 16*x2*x2*a*a + 64*x2*x2*x2*x2*a*a*a*a + 8*x1*x1*a*a*
                              (-3 + 8*x2*x2*a*a) + 16*x1*(x2*a*a + 8*x2*x2*x2*a*a*a*a)) +
                        4*b2*b2*b2*(2*a*a*(9 - 12*x1*x1*a*a + 8*x1*x2*a*a + 8*x2*x2*a*a) +
                                b3*(9 + 8*x2*x2*a*a + 8*x2*x3*a*a + 64*x2*x2*x2*x3*a*a*a*a +
                                    8*x1*x1*a*a*(-3 + 4*x2*x2*a*a + 4*x2*x3*a*a) + 4*x1*a*a*
                                    (3*x2 + x3 + 8*x2*x2*x2*a*a + 24*x2*x2*x3*a*a))) +
                        4*b2*b3*(72*a*a*a*a + 2*b3*a*a*(27 - 36*x1*x1*a*a + 16*x2*x3*a*a + 8*x3*x3*a*a +
                                                                  8*x1*(x2 + 2*x3)*a*a) + b3*b3*(9 + 8*x2*x3*a*a + 8*x3*x3*a*a +
                                                                                                       64*x2*x3*x3*x3*a*a*a*a + 8*x1*x1*a*a*(-3 + 4*x2*x3*a*a + 4*x3*x3*a*a) +
                                                                                                       4*x1*a*a*(x2 + 3*x3 + 24*x2*x3*x3*a*a + 8*x3*x3*x3*a*a))) +
                        2*b2*b2*(72*a*a*a*a + 4*b3*a*a*(27 - 36*x1*x1*a*a + 8*x2*x2*a*a + 16*x2*x3*a*a +
                                                                 8*x1*(2*x2 + x3)*a*a) + b3*b3*(27 + 32*x2*x3*a*a + 8*x3*x3*a*a +
                                                                                                      8*x1*x1*a*a*(-9 + 4*x2*x2*a*a + 16*x2*x3*a*a + 4*x3*x3*a*a) +
                                                                                                      24*x1*a*a*(x2 + x3 + 8*x2*x2*x3*a*a + 8*x2*x3*x3*a*a) +
                                                                                                      8*x2*x2*(a*a + 24*x3*x3*a*a*a*a))) +
                        b3*b3*(144*a*a*a*a + 8*b3*a*a*(9 - 12*x1*x1*a*a + 8*x1*x3*a*a + 8*x3*x3*a*a) +
                              b3*b3*(9 + 16*x3*x3*a*a + 64*x3*x3*x3*x3*a*a*a*a + 8*x1*x1*a*a*(-3 + 8*x3*x3*a*a) +
                                    16*x1*(x3*a*a + 8*x3*x3*x3*a*a*a*a)))))))/
(sqrt(((b0 + b1)*(b2 + b3))/(4*(b2 + b3)*a*a + b0*(b2 + b3 + 4*a*a) +
                             b1*(b2 + b3 + 4*a*a)))*pow(4*(b2 + b3)*a*a + b0*(b2 + b3 + 4*a*a) +
                                                        b1*(b2 + b3 + 4*a*a),5));


}
double aGGCGG(double a , struct general_2index * pa){
    
    double b0 = pa->i[0].b0;
    double b1 = pa->i[0].b1;
    double x0 = pa->i[0].x0;
    double x1 = pa->i[0].x1;
    unsigned short l0 = pa->i[0].l0;
    unsigned short l1 = pa->i[0].l1;
    
    double b2 = pa->i[1].b0;
    double b3 = pa->i[1].b1;
    double x2 = pa->i[1].x0;
    double x3 = pa->i[1].x1;
    unsigned short l2 = pa->i[1].l0;
    unsigned short l3 = pa->i[1].l1;
    
    
    
    double f00 =  1./a* (4.*1.7724538509055159)/exp((b0*(b1*(b2 + b3)*sqr (x0 - x1) + b2*b3*sqr(x2 - x3)) + a*a*(b0*(b1*sqr(x0 - x1) + b2*sqr(x0 - x2) + b3*sqr(x0 - x3)) +   b1*(b2*sqr(x1 - x2) + b3*sqr(x1 - x3)) + b2*b3*sqr(x2 - x3)) +  b1*b2*b3*sqr(x2 - x3))/((b0 + b1)*(b2 + b3) + a*a*(b0 + b1 + b2 + b3)))/sqrt( sqrt( 1./b0/b1/b2/b3) * (b0+b1)*(b2+b3) *((b0 + b1)*(b2 + b3) + a*a*(b0 + b1 + b2 + b3))/(a*a*(b0 + b1)*(b2 + b3)));
    
    if ( l0 == 0 && l1 == 0 && l2 == 0 && l3 == 0 )
        return f00;
    /*signs added*/    /*now signed removed*/
    
    /*removed minus signs in this block*/
    if ( l0 == 1 && l1 == 0 && l2 == 0 && l3 == 0 )
        return c10c00(a,b0,x0,b1,x1,b2,x2,b3,x3)*f00;
    if ( l0 == 0 && l1 == 1 && l2 == 0 && l3 == 0 )
        return c10c00(a,b1,x1,b0,x0,b2,x2,b3,x3)*f00;
    if ( l0 == 0 && l1 == 0 && l2 == 1 && l3 == 0 )
        return c10c00(a,b2,x2,b3,x3,b0,x0,b1,x1)*f00;
    if ( l0 == 0 && l1 == 0 && l2 == 0 && l3 == 1 )
        return c10c00(a,b3,x3,b2,x2,b0,x0,b1,x1)*f00;
    
    
    
    if ( l0 == 1 && l1 == 1 && l2 == 0 && l3 == 0 )
        return c11c00(a,b0,x0,b1,x1, b2,x2,b3,x3 )*f00;
    if ( l0 == 0 && l1 == 0 && l2 == 1 && l3 == 1 )
        return c11c00(a, b2,x2,b3,x3,b0,x0,b1,x1 )*f00;
    
    if ( l0 == 1 && l1 == 0 && l2 == 1 && l3 == 0 )
        return c10c10(a, b0,x0,b1,x1,b2,x2,b3,x3 )*f00;
    if ( l0 == 1 && l1 == 0 && l2 == 0 && l3 == 1 )
        return c10c10(a, b0,x0,b1,x1,b3,x3,b2,x2 )*f00;
    if ( l0 == 0 && l1 == 1 && l2 == 1 && l3 == 0 )
        return c10c10(a, b1,x1,b0,x0,b2,x2,b3,x3 )*f00;
    if ( l0 == 0 && l1 == 1 && l2 == 0 && l3 == 1 )
        return c10c10(a, b1,x1,b0,x0,b3,x3,b2,x2 )*f00;
    
    
    /*removed minus signs in this block*/
    if(  l0 == 1 && l1 == 1 && l2 == 1 && l3 == 0 )
        return c11c10(a,b1,x1,b0,x0,b2,x2,b3,x3)*f00;
    if(  l0 == 1 && l1 == 1 && l2 == 0 && l3 == 1 )
        return c11c10(a,b1,x1,b0,x0,b3,x3,b2,x2)*f00;
    if(  l0 == 1 && l1 == 0 && l2 == 1 && l3 == 1 )
        return c11c10(a,b2,x2,b3,x3,b0,x0,b1,x1)*f00;
    if(  l0 == 0 && l1 == 1 && l2 == 1 && l3 == 1 )
        return c11c10(a,b2,x2,b3,x3,b1,x1,b0,x0)*f00;
    if ( l0 == 1 && l1 == 1 && l2 == 1 && l3 == 1 )
        return c11c11(a,b0,x0,b1,x1,b2,x2,b3,x3)*f00;
    
    
    if ( l0 == 2 && l1 == 0 && l2 == 0 && l3 == 0 )
        return c20c00(a,b0,x0,b1,x1,b2,x2,b3,x3)*f00;
    if ( l0 == 0 && l1 == 2 && l2 == 0 && l3 == 0 )
        return c20c00(a,b1,x1,b0,x0,b2,x2,b3,x3)*f00;
    if ( l0 == 0 && l1 == 0 && l2 == 2 && l3 == 0 )
        return c20c00(a,b2,x2,b3,x3,b0,x0,b1,x1)*f00;
    if ( l0 == 0 && l1 == 0 && l2 == 0 && l3 == 2 )
        return c20c00(a,b3,x3,b2,x2,b0,x0,b1,x1)*f00;

    
    if ( l0 == 2 && l1 == 2 && l2 == 0 && l3 == 0 )
        return c20c00(a,b0,x0,b1,x1,b2,x2,b3,x3)*f00;
    if ( l0 == 0 && l1 == 0 && l2 == 2 && l3 == 2 )
        return c20c00(a,b2,x2,b3,x3,b0,x0,b1,x1)*f00;

    return 0;
}


double collective( double beta ,struct general_2index * pa){
    
    if ( pa->gaussianAccelerationFlag ){
        return aGGCGG(beta,pa)/(2*pi);
    }

    
    pa->alpha = beta;
    double value= 0.,value2=0.;
    if ( pa->periodic == 1 ){
        
        if ( beta < 1e-9 ){
            return 0.;
        }
        double kSmall = (2*pi/pa->N1/max(pa->i[0].d,pa->i[1].d));
        INT_TYPE k;
        for ( k = - pa->N1 ; k <= pa->N1 ; k++){
            value += gaussianSinc(kSmall*k+pa->momentumShift,pa)*kSmall;
        }
    }else
        if ( pa->periodic == -1 ){
            value = gaussianSinc(0,pa);//like sheets of positive (external) charge...  yields a simple momentum integral.
        }else {

#ifdef APPLE
        quadrature_integrate_function g;
        g.fun = gaussianSincFunc;                               // Called to evaluate the function to integrate
        g.fun_arg = pa;                                  // Passed as first argument to the callback
        
        quadrature_integrate_options options;
        bzero(&options,sizeof(options));
        options.integrator = QUADRATURE_INTEGRATE_QAGS;    // Use QAGS: adaptive subdivision with convergence acceleration
        
        options.abs_tolerance = 1.0e-9;                    // Requested absolute tolerance on result
        options.max_intervals = 20;                        // Max number of intervals
        options.qag_points_per_interval = 25;
        quadrature_status status;
        double value,abs_error;
        
            if ( (pa->i[0].bra.basis == SincBasisElement && pa->i[0].ket.basis == SincBasisElement ) || (pa->i[1].bra.basis == SincBasisElement && pa->i[1].ket.basis == SincBasisElement )  || (pa->i[0].bra.basis == SincBasisElement && pa->i[0].ket.basis == nullFunction )|| (pa->i[0].bra.basis == nullFunction && pa->i[0].ket.basis == SincBasisElement ) || (pa->i[1].bra.basis == SincBasisElement && pa->i[1].ket.basis == nullFunction )|| (pa->i[1].bra.basis == nullFunction && pa->i[1].ket.basis == SincBasisElement ) ){
                INT_TYPE pt = imin(((pa->i[0].bra.basis == SincBasisElement) + (pa->i[0].ket.basis == SincBasisElement) ),((pa->i[1].bra.basis == SincBasisElement) + (pa->i[1].ket.basis == SincBasisElement) ));
                if ( pt == 0 )
                    pt = 1;

            value =  quadrature_integrate(&g, -pt*pi/max(pa->i[0].d,pa->i[1].d)+pa->momentumShift,pt*pi/max(pa->i[0].d,pa->i[1].d)+pa->momentumShift, &options, &status, &abs_error, 0, NULL);
            }
#else
        double abs_error;
        gsl_function F;
        F.function = &gaussianSinc;
        F.params = pa;
        
        gsl_integration_workspace * workspace= gsl_integration_workspace_alloc (1000);
        if ( (pa->i[0].bra.basis == SincBasisElement && pa->i[0].ket.basis == SincBasisElement ) || (pa->i[1].bra.basis == SincBasisElement && pa->i[1].ket.basis == SincBasisElement )  || (pa->i[0].bra.basis == SincBasisElement && pa->i[0].ket.basis == nullFunction )|| (pa->i[0].bra.basis == nullFunction && pa->i[0].ket.basis == SincBasisElement ) || (pa->i[1].bra.basis == SincBasisElement && pa->i[1].ket.basis == nullFunction )|| (pa->i[1].bra.basis == nullFunction && pa->i[1].ket.basis == SincBasisElement ) ){
                INT_TYPE pt = imin(((pa->i[0].bra.basis == SincBasisElement) + (pa->i[0].ket.basis == SincBasisElement) ),((pa->i[1].bra.basis == SincBasisElement) + (pa->i[1].ket.basis == SincBasisElement) ));
                if ( pt == 0 )
                    pt = 1;
                gsl_integration_qag (&F,  -pt*pi/max(pa->i[0].d,pa->i[1].d)+pa->momentumShift,  pt*pi/max(pa->i[0].d,pa->i[1].d)+pa->momentumShift, 1e-9, 1e-9,1000,6,workspace, &value, &abs_error);
            }
        else{
            gsl_integration_qagi (&F, 1e-9, 1e-9,1000,workspace, &value, &abs_error);
        }
        gsl_integration_workspace_free(workspace);
#endif
    }
    value2 = 0.0;
    if ( pa->fl->fn == Yukawa ){
        double m = pa->fl->param[2];
        value2  += exp(-sqr(m/2./beta));
    }
    else if ( pa->fl->fn == Morse ){
        double R = pa->fl->param[2];
        double a = pa->fl->param[3];
        value2  += - a * exp (     R * a - sqr(a/beta /2.))/sqr(beta);
        value2  +=   a * exp ( 2 * R * a - sqr(a/beta    ))/sqr(beta);
    }
    else if ( pa->fl->fn == Coulomb || pa->fl->fn == Pseudo || pa->fl->fn == nullFunction ){
        value2 += 1.;//
    } else {
    }
    
    
    value2 *= pa->fl->param[0];
    if ( value2 < 0 )
        value2 = -pow(-value2 , 1./SPACE);
    else
        value2 = pow(value2 , 1./SPACE);

    return value*value2;

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
    return collectives(beta,af)*collectives(beta,af+1)*collectives(beta,af+2);
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
    
    options.abs_tolerance = 1.0e-8;                    // Requested absolute tolerance on result
    options.max_intervals = 20;                        // Max number of intervals
    //options.qag_points_per_interval = 25;
    quadrature_status status;
    double  value,abs_error;
    value =  quadrature_integrate(&g, a, b, &options, &status, &abs_error, 0, NULL);
#else
    double value, abs_error;

    gsl_function F;
    F.function = &element;
    F.params = aAf;
    size_t neval;
    gsl_integration_workspace * workspace= gsl_integration_workspace_alloc (1000);

    if ( b < a)
        gsl_integration_qagiu(&F, a, 1e-8, 1e-8,1000,workspace, &value, &abs_error);
    else
        gsl_integration_qag (&F,  a,  b, 1e-8, 1e-8,1000,4,workspace, &value, &abs_error);
    gsl_integration_workspace_free(workspace);

#endif
    return value*(2*pi);

}






void mySeparateExactOne (struct field * f1, double scalar, enum division basis){
    
    
    INT_TYPE i,I1,I2,I3,I4,alpha,sp,space;
    INT_TYPE N1 = f1->sinc.N1;
    INT_TYPE N2 = N1*N1;
    INT_TYPE N12 = (N1-1)/2;
    double value,d = f1->sinc.d;
    if ( scalar < 0 )
        scalar = -pow(fabs(scalar), 1./SPACE);
    else
        scalar =  pow(fabs(scalar), 1./SPACE);
    
    f1->sinc.tulip[diagonalCube].header = Cube;
    
    
    tClear(f1, quadCube);
    f1->sinc.tulip[quadCube].Current[0] = 1;
    
    for ( alpha = 0 ; alpha < CanonicalRank(f1, interactionExchange, 0) ; alpha++){
        zero(f1, quadCube,0);
        
#pragma omp parallel for private (I1,I2,I3,I4,value,sp,space) schedule(dynamic,1)
        for ( I1 = 0 ;I1 < N1 ;I1++){// body 0 -1
            for ( I2 = 0 ;I2 < N1 ;I2++)// body 0 -2
                for ( I3= 0 ;I3 < N1 ;I3++)// body 1 -1
                    for ( I4 = 0 ;I4 < N1 ;I4++)// body 1 -2
                        for ( space = 0; space < SPACE ; space++){
                            value  = streams(f1, interactionExchange, 0,space)[ N2*N2 * alpha + (N2*(N1*I1+I3)+(N1*I2+I4))];
                            streams(f1, quadCube,0,space)[(N2*(N1*I1+I2)+(N1*I3+I4))]=value*scalar;
                        }
        }
        
        if ( basis ){
           // tTransformTo(f1, quadCube, quad,basis );
            tAddTwo(f1, interactionDirect, quad);
        } else {
            tAddTwo(f1, interactionDirect, quadCube);
        }
    }
}

void mySeparateExactTwo (struct field * f1, INT_TYPE periodic, double scalar,  enum division basis,INT_TYPE plus, INT_TYPE particle1,INT_TYPE particle2){
    //https://keisan.casio.com/exec/system/1329114617
    if ( periodic < 0 ){
        printf("error in periodic flag\n");
        return;
    }
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
    
    INT_TYPE I1,I2,I3,I4,beta,section;
    INT_TYPE N1 = f1->sinc.N1,space;
    INT_TYPE N12 = (f1->sinc.N1-1)/2;
    INT_TYPE N2 = f1->sinc.N1*f1->sinc.N1;
    INT_TYPE *dims1 = f1->sinc.dims[particle1];
    INT_TYPE *dims2 = f1->sinc.dims[particle2];

    
    
    enum functionType fn = f1->twoBody.func.fn;
    if ( fn == nullFunction )
        return ;
    
    struct general_2index g2;

    g2.realFlag = 1;
    double value,d = f1->sinc.d,g,x;
    double constant;
    INT_TYPE si,interval = f1->twoBody.func.interval;
    double * param   = f1->twoBody.func.param;
    getDescription(&f1->twoBody.func, scalar, stdout);
    
    double *gkX, *gkW;
    INT_TYPE ngk;
    f1->sinc.tulip[quadCube].header = Cube;
    
    tClear(f1, quadCube);
    f1->sinc.tulip[quadCube].Current[0] = 1;
    
    for ( section = 0 ; section < 2 ; section++){
        if ( imax(0,interval - section) == 0 ) {
                gkX = gk7X;
                gkW = gk7W;
                ngk = 7;
            }else if ( imax(0,interval - section) == 1 ){
                gkX = gk15X;
                gkW = gk15W;
                ngk = 15;
            }else if ( imax(0,interval - section) == 2 ){
                gkX = gk35X;
                gkW = gk35W;
                ngk = 35;
            }else if ( imax(0,interval - section) == 3 ){
                gkX = gk99X;
                gkW = gk99W;
                ngk = 99;
            }else {
                printf("oops\n");
                exit(0);
            }

            
            
        
            for ( beta = 0; beta < ngk ; beta++){
                
                if (1.){
                    
                    g = (gkX[beta]+1)/2.;
                    if ( section == 0 ){
                        x = g ;
                        constant = gkW[beta];
                    }
                    else {
                        x = ( g ) / (1. - g)+1 ;
                        constant = gkW[beta]/sqr(1.-g);
                    }
                    x *=  param[1];
                    constant *= 1./2.*param[1];
                }
                //divide 0->infinity
                zero(f1, quadCube,0);
                
                if ( fn == Pseudo){
                    if ( section == 1  )
                        continue;;
                }
                zero(f1, quadCube,0);
                
                constant *= scalar;//for purposes of scaling tw-body interactions by Ne...
                
                
                if ( constant < 0 )
                    constant = -pow(fabs(constant)*2.*pi, 1./SPACE);
                else
                    constant =  pow(fabs(constant)*2.*pi, 1./SPACE);;//*3.544907701811032;
                
              //  for ( space = 0 ; space < SPACE ; space++)//technically unnecessary , still in flux
                space = SPACE-1;

#ifdef OMP
#pragma omp parallel for private (si,I1,I2,I3,I4,g2,value) schedule(dynamic,1)
#endif
                for ( si = 0 ; si < dims1[space]*dims1[space]*dims2[space]*dims2[space] ; si++)
                
//
//                for ( I1 = 0 ;I1 < dims1[space] ;I1++){// body 0 -1
//                    for ( I2 = 0 ;I2 <  dims1[space] ;I2++)// body 0 -2
//                        for ( I3 = 0 ;I3 < dims2[space] ;I3++)// body 1 -1
//                            for ( I4 = 0 ;I4 < dims2[space] ;I4++)// body 1 -2
                            {
                                I1 = ( si ) % dims1[space];
                    
                                I2 = ( si / (dims1[space])) % dims1[space];
                                I3 = ( si /  (dims1[space]*dims1[space]) ) % dims2[space];
                                I4 = ( si / ( dims1[space]*dims1[space]*dims2[space]) ) % dims2[space];

                                g2.realFlag = 1;
                                g2.momentumShift = 0;
                                g2.gaussianAccelerationFlag = 0;
                                g2.point = 2;
                                g2.powSpace = 0;
                                g2.body = 2;
                                g2.N1 = N1;
                                g2.fl = & f1->twoBody.func;
                                
                                if ( space == 0 )
                                    g2.periodic = ( abs(periodic) ) % 2;
                                else if ( space == 1 )
                                    g2.periodic = ( abs(periodic) / 2 ) % 2;
                                else if (space == 2 )
                                    g2.periodic = ( abs(periodic) / 4) % 2;
                                g2.i[0].action = 0;
                                g2.i[1].action = 0;

#ifdef oldIndex
                                
                                g2.i[0].d = f1->sinc.d;
                                g2.i[1].d = f1->sinc.d;
                                if ( plus ) {
                                    //plus defines 1/(r+R)... and is essentially negative.
                                    
                                    g2.i[0].n = -(I1-N12);
                                    g2.i[0].m = -(I2-N12);
                                    g2.i[1].n = (I3-N12);
                                    g2.i[1].m = (I4-N12);
                                } else {
                                    g2.i[0].n = I1-N12;
                                    g2.i[0].m = I2-N12;
                                    g2.i[1].n = I3-N12;
                                    g2.i[1].m = I4-N12;
                                }
                                
                                
                                
#else
                                g2.i[0].d = grabBasis(f1, space, particle1, I1).length;
                                g2.i[1].d = grabBasis(f1, space, particle2, I3).length;

                                if ( plus ) {
                                    //plus defines 1/(r+R)... and is essentially negative.
                                    
                                    g2.i[0].bra = grabBasis ( f1, space, particle1, dims1[space]-I1-1);
                                    g2.i[0].ket = grabBasis ( f1, space, particle1, dims1[space]-I2-1);
                                    
                                    g2.i[1].bra = grabBasis ( f1, space, particle2, I3);
                                    g2.i[1].ket = grabBasis ( f1, space, particle2, I4);
                                } else {
                                    g2.i[0].bra = grabBasis ( f1, space, particle1, I1);
                                    g2.i[0].ket = grabBasis ( f1, space, particle1, I2);
                                    
                                    g2.i[1].bra = grabBasis ( f1, space, particle2, I3);
                                    g2.i[1].ket = grabBasis ( f1, space, particle2, I4);
                                }

#endif
                                
                                
                                value =  collectives(x,&g2);
                                value *= constant;
                                
                                streams(f1, quadCube,0,space)[dims1[space]*dims2[space]*(dims2[space]*I3+I1)+(dims2[space]*I4+I2)]=value;
                            }
                    
                
                printf("%d %d %15.36f \n",section,beta,cblas_dnrm2 ( dims1[space]*dims2[space]*dims1[space]*dims2[space] , streams(f1, quadCube, 0,space), 1 ));
                
                
//                if ( pow(1e-9,SPACE) >  pow(cblas_dnrm2 ( dims1[space]*dims2[space]*dims1[space]*dims2[space] , streams(f1, quadCube, 0,0), 1 ) ,SPACE)){
//                    printf("skip %lld\n",beta);
//                    continue;
//                }
                if ( plus )
                    tAddTw(f1, interactionExchangePlus, 0,quadCube,0);
                else
                    tAddTw(f1, interactionExchange, 0,quadCube,0);

            }
        }
//    printf("RMS %15.15f @ %d %d\n", tRMSDevRandom( f1, interactionExchange, periodic ,10000),10000,CanonicalRank(f1, interactionExchange, 0));
//    exit(0);
}

INT_TYPE separateExternal( struct calculation * c1,INT_TYPE periodic, INT_TYPE atom,double scalar, INT_TYPE dim, enum division basis , INT_TYPE particle1){
    //https://keisan.casio.com/exec/system/1329114617
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
        0.0,
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

    
    
    struct field * f1 = &c1->i.c;
    
    INT_TYPE section;
    INT_TYPE *dims1 = f1->sinc.dims[particle1];

    INT_TYPE gPoint = 2,i,beta,I1,I2,a,space,ii,spin;
    double constant,value,x,g;
    Stream_Type  *stream[3];
    struct function_label fl;
    struct general_2index g2;

    f1->sinc.tulip[diagonalCube].header = Cube;
    for ( i = 0; i < SPACE ; i++)
        stream[i] =  streams( f1, diagonalCube,0,i  );
    zero(f1, diagonalCube , 0);
    f1->sinc.tulip[diagonalCube].Current[0] = 1;
    
    
    /////IGNORE TRAIN COMMAND
    INT_TYPE *n1,*m1;
    if ( basis ){
        f1->sinc.tulip[linear].header = f1->sinc.tulip[basis].header-1;
        n1 = f1->sinc.Basis;
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
        n1 = f1->sinc.Basis;
        m1 = n1 ;
        f1->sinc.tulip[diagonalCube].header = Cube;
        f1->sinc.tulip[linear].header = Cube;
        
    }
    spin = 0;
    f1->sinc.tulip[linear].stop[0][0] = CanonicalRank(f1, linear, 0);
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
    getDescription(&f1->oneBody.func, scalar, stdout);

    double *gkX, *gkW,Z, gGX[10],gGW[10];
    INT_TYPE si,ngk,ai,powSpace[SPACE][10];
    
    f1->sinc.tulip[quadCube].header = Cube;
    struct basisElement boa;
    tClear(f1, quadCube);
    f1->sinc.tulip[quadCube].Current[0] = 1;
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
    for ( a = mA ; a <= xA; a++){
        ii = 0;
        printf("atom - %d--Z = %d--\n",  a,f1->atoms[a].label.Z);

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
            getDescription(&fl, 1., stdout);
            

            
        }
        else{
            fl.param[0] = param[0];
            fl.param[1] = param[1];
            fl.param[2] = param[2];
            fl.param[3] = param[3];
            fl.fn = fn;
            Z = f1->atoms[a].label.Z;
        }
        
        for ( section = 0 ; section < 9 ; section++){
            ngk = 0;
            if ( (( fn == LDA || fn == BLYP )&& section == 0) || ( ( fn != LDA && fn != BLYP ) && section < 2)  ){
                gPoint = 2;
                
                
                
                
                if ( imax(0,interval - section) == 0 ) {
                    gkX = gk7X;
                    gkW = gk7W;
                    ngk = 7;
                }else if ( imax(0,interval - section) == 1 ){
                    gkX = gk15X;
                    gkW = gk15W;
                    ngk = 15;
                }else if ( imax(0,interval - section) == 2 ){
                    gkX = gk35X;
                    gkW = gk35W;
                    ngk = 35;
                }else if ( imax(0,interval - section) == 3 ){
                    gkX = gk99X;
                    gkW = gk99W;
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
                        powSpace[0][0] = 0;
                        powSpace[1][0] = 0;
                        powSpace[2][0] = 0;
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
                    
                    if (1.){
                        
                        if ( section ==  0 ){
                            x = (gkX[beta]+1)/2. ;
                            constant = gkW[beta];
                            
                            x *=  fl.param[1];
                            constant *= 1./2.*fl.param[1];

                        }
                        else if ( section == 1 ){
                            g = (gkX[beta]+1)/2.;
                            x = ( g ) / (1. - g)+1 ;
                            constant = gkW[beta]/sqr(1.-g);
                            
                            x *=  fl.param[1];
                            constant *= 1./2.*fl.param[1];

                        }else {
                            x = gkX[beta];
                            constant = gkW[beta];
                        }
                    }
                    
                    if ( fn == Pseudo ){
                        if ( section == 1 )
                            continue;;
                    }
                    zero(f1, diagonalCube,0);
                    
                    constant *= -scalar;//attractive
                    
                    
                    if ( constant < 0 )
                        constant = -pow(fabs(constant)*2.*pi*fabs(Z), 1./SPACE);
                    else
                        constant =  pow(fabs(constant)*2.*pi*fabs(Z), 1./SPACE);

                for ( space = 0 ;space < SPACE ; space++)
#ifdef OMP
#pragma omp parallel for private (si,I1,I2,g2,value) schedule(dynamic,1)
#endif
                    for ( si = 0; si < dims1[space]*dims1[space]; si++)
//                        for ( I1 = 0; I1 < dims1[space] ; I1++)
//                            for ( I2 = 0; I2 < dims1[space] ; I2++)
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
                                g2.body = 1;
                                g2.fl = & fl;
                                
                                if ( space == 0 )
                                    g2.periodic = ( periodic ) % 2;
                                else if ( space == 1 )
                                    g2.periodic = ( periodic / 2 ) % 2;
                                else if (space == 2 )
                                    g2.periodic = ( periodic / 4) % 2;
                                if ( periodic < 0 )
                                    g2.periodic *= -1;//make it a simple periodic external field( plane-line-like).


#ifdef oldIndex

                                g2.i[0].n = I1-N12;
                                g2.i[0].m = I2-N12;
                                g2.i[1].x0 = f1->atoms[a].position[space+1];
                                g2.i[1].action = ( dim == space ) ;
                                
#else
                                g2.i[0].d = grabBasis(f1, space, particle1, I1).length;
                                g2.N1 = dims1[space];
                                
                                boa.basis = DiracDelta;
                                boa.origin = f1->atoms[a].position[space+1];
                                
                                g2.i[0].bra = grabBasis ( f1, space, particle1, I1);
                                g2.i[0].ket = grabBasis ( f1, space, particle1, I2);
                                
                                g2.i[1].bra = boa;
                                g2.i[1].ket = boa;
                                g2.i[1].ket.basis = nullBasisElement;
                                g2.i[1].action = ( dim == space ) ;

#endif
                                value = collectives(x, &g2);
                                value *= constant;
                                (stream[space])[I1*dims1[space]+I2] = value;
                            }
                    
                    va = 1;
                    for ( space = 0 ;space < SPACE ; space++)
                        va *= cblas_dnrm2 ( dims1[space]*dims1[space] , stream[space], 1 ) ;
                    if ( 1e-12 > va ){
                        printf("skip %d %d %d\n",a, section,beta);
                        continue;
                    }
                    tAddTw(f1, linear,0, diagonalCube,0);
                }
                
                f1->sinc.tulip[linear].stop[0][a] = CanonicalRank(f1, linear,0);
#if VERBOSE
            printf("stop %d %d\n", a, f1->sinc.tulip[linear].stop[0][a]);
#endif
            }
    }
    f1->sinc.tulip[linear].stop[0][f1->Na+1] = CanonicalRank(f1, linear,0);
    return 0;
}

INT_TYPE separateKinetic( struct field * f1, INT_TYPE periodic,enum division akinetic,  double amass, INT_TYPE particle1 ){
    INT_TYPE space,dim,I1,I2;
    INT_TYPE * dims1 = f1->sinc.dims[particle1];
    Stream_Type * stream;
    
    for ( space = 0 ;space < SPACE; space++){
        
        for ( dim = 0 ; dim < SPACE; dim++){
            stream =  streams( f1, akinetic, 0 , space ) + dims1[space]*dims1[space] * dim;//last one is COMPLEX:: cmpl!
            
            for ( I1 = 0 ; I1 < dims1[space] ; I1++)
                for( I2 = 0; I2 < dims1[space] ; I2++){
                    
                    if ( dim == space )
                        (stream )[dims1[space]*I1+I2] = - 0.5*Bd2B(grabBasis(f1, space, particle1, I1),grabBasis(f1, space, particle1, I2));
                    else
                        (stream )[dims1[space]*I1+I2] = BoB(grabBasis(f1, space, particle1, I1),grabBasis(f1, space, particle1, I2));

                    
                    //printf ("%d %d %d %f %f\n", particle1,I1,I2,sqr(f1->sinc.d)* Bd2B(grabBasis(f1, space, particle1, I1),grabBasis(f1, space, particle1, I2)), Sd2S(I1-I2));
                }
        }
    }
    
    f1->sinc.tulip[akinetic].Current[0]  = SPACE;

    tScale(f1, akinetic, 1./amass);
#if VERBOSE
    printf("kinetic %f\n", traceOne(f1, akinetic, 0));
#endif
    return 0;
}

INT_TYPE separateHarmonicExternal( struct calculation * c1,INT_TYPE periodic, double scalar,INT_TYPE dim, enum division basis, INT_TYPE particle1 ){
    
    struct field * f1 = &c1->i.c;
    
    INT_TYPE *dims1 = f1->sinc.dims[particle1];
    INT_TYPE i,alpha,I1,I2,space;
    Stream_Type  *stream[3];
    
    f1->sinc.tulip[diagonalCube].header = Cube;
    for ( i = 0; i < SPACE ; i++)
        stream[i] =  streams( f1, diagonalCube,0,i  );
    zero(f1, diagonalCube , 0);
    f1->sinc.tulip[diagonalCube].Current[0] = 1;
    
    
    /////IGNORE TRAIN COMMAND
    INT_TYPE *n1,*m1;
    if ( basis ){
        f1->sinc.tulip[harmonium].header = f1->sinc.tulip[basis].header-1;
        n1 = f1->sinc.Basis;
        m1 = n1;
        f1->sinc.tulip[diagonalCube].header = Cube;
    }
    //    else if ( header (f1, harmonium ) == SincSub  ){
    //        n1 = f1->sinc.n1;
    //        m1 = f1->sinc.n1;
    //        f1->sinc.tulip[diagonalCube].header = SincSub;
    //        f1->sinc.tulip[harmonium].header = SincSub;
    //    }
    else if ( header (f1, harmonium ) == Cube ){
        n1 = f1->sinc.Basis;
        m1 = n1 ;
        f1->sinc.tulip[diagonalCube].header = Cube;
    }
    else{
        printf("error with harmonium build\n");
        exit(0);
    }
    {
        for ( alpha = 0 ; alpha < SPACE ; alpha++){
            zero(f1, diagonalCube, 0);
            for ( space = 0 ;space < SPACE ; space++){
                for ( I1 = 0; I1 < dims1[space] ; I1++)
                    for ( I2 = 0; I2 < dims1[space] ; I2++)
                    {
                        if ( space == alpha ){
                            (stream[space])[I1*dims1[space]+I2] = 0;
                            //scalar * 0.500 * (c1->i.springConstant) *Bx2B(grabBasis(f1, space, particle1, I1),grabBasis(f1, space, particle1, I2));
;
                        } else {
                            (stream[space])[I1*dims1[space]+I2] =  BoB(grabBasis(f1, space, particle1, I1),grabBasis(f1, space, particle1, I2));
                        }
                    }
            }
            tAddTw(f1, harmonium,0, diagonalCube,0);
        }
    }
    if ( part(f1,harmonium ) < SPACE )
    {
        printf("harmonium\n");
        exit(0);
    }
    return 0;
}


double tTestTwoBody( struct field * f1, enum division mat,INT_TYPE periodic, INT_TYPE * p){
    INT_TYPE r,space,N1 = f1->sinc.N1;
    
    
    double sum = 0.,product;
    for ( r = 0; r < CanonicalRank(f1, mat, 0 );r++){
        product = 1;
        for ( space = 0 ; space < SPACE ; space++)
            product *= streams(f1, mat, 0,space)[p[4*space] + N1*p[4*space+1]+ N1*N1*p[4*space+2]+N1*N1*N1*p[4*space+3]+r*N1*N1*N1*N1];
        sum += product;
    }
    INT_TYPE i ;
    struct general_2index g3[3];
    for ( i = 0; i < 3; i++){
        g3[i].realFlag = 1;
        g3[i].i[0].d = f1->sinc.d;
        g3[i].i[1].d = f1->sinc.d;

        g3[i].i[0].n = p[i*4];
        g3[i].i[0].m = p[i*4+2];
        g3[i].i[1].n = p[i*4+1];
        g3[i].i[1].m = p[i*4+3];
        g3[i].point = 2;
        g3[i].fl = & f1->twoBody.func;
        g3[i].powSpace = 0;
        if ( i == 0 )
            g3[i].periodic = ( periodic ) % 2;
        else if ( i == 1 )
            g3[i].periodic = ( periodic / 2 ) % 2;
        else if (i == 2 )
            g3[i].periodic = ( periodic / 4) % 2;
        else
        {
            printf("here\n");
            exit(0);
        }
        g3[i].body = 2;
        g3[i].N1 = N1;
    }
   // printf("%1.15f \t %1.15f\n", sum, elementCal(1e-3,-1, g3));
    {
        double value;
        value = sum - elementCal(1e-3,-1, g3);
        if ( isnan(value)){
            printf(".");
            return 0.;
        } else
            return sqr(value);
    }
}

double tRMSDevRandom( struct field * f1, enum division mat, INT_TYPE periodic ,INT_TYPE Nc){
    INT_TYPE p[12],i,j;
    double sum = 0.;
    for ( i = 0; i < Nc ; i++){
        for ( j = 0 ; j < 12 ; j++)
            p[j] = rand()% f1->sinc.N1;
        sum += tTestTwoBody(f1, mat, periodic, p);
    
    }
    return sqrt(sum/Nc);
}

INT_TYPE buildElectronProtonInteraction ( struct field * f1, enum division mat){
    INT_TYPE space,r,i,j,n,m,N1 = f1->sinc.N1, N2 = f1->sinc.N1*f1->sinc.N1;
    double value;
    double coef = pow(f1->Ne,1./SPACE);
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
                    streams(f1, mat,0,space)[N1*i+j+r*N2] = -coef*value/N1;
            }
    }
    f1->sinc.tulip[mat].Current[0] = CanonicalRank(f1, interactionExchange, 0);

    if ( bodies(f1, eigenVectors)==two){
        tClear(f1, copy);
        tId(f1, copy,0);
        tScale(f1, copy,-3.*traceOne(f1, mat, 0)/(N1*N1*N1)/sqr(f1->Ne));
        //
        tAddTw(f1, mat, 0,copy,0);
    }
    else     if ( bodies(f1, eigenVectors)==three){
        tClear(f1, copy);
        tId(f1, copy,0);
        tScale(f1, copy,-6.*traceOne(f1, mat, 0)/(N1*N1*N1)/sqr(f1->Ne));
        //
        tAddTw(f1, mat,0, copy,0);
    }
    else     if ( bodies(f1, eigenVectors)==four){
        tClear(f1, copy);
        tId(f1, copy,0);
        tScale(f1, copy,-10.*traceOne(f1, mat, 0)/(N1*N1*N1)/sqr(f1->Ne));
        //
        tAddTw(f1, mat, 0,copy,0);
    }
    
    
    f1->sinc.tulip[mat].stop[0][0] = CanonicalRank(f1, mat,0);
    
    return 0;
}


void separateX ( struct field * f1,  double vectorDipole ){
    INT_TYPE space, i;
    INT_TYPE N1 = f1->sinc.N1;
    INT_TYPE N12 = (f1->sinc.N1-1)/2;
    INT_TYPE N2 = f1->sinc.N1*f1->sinc.N1;
    tClear(f1, X);
    zero(f1, X,0);
    if ( fabs(vectorDipole) < 1e-6 )
        return;
    tId(f1, X,0);
    for ( space = 0; space < 1 ; space++){
        for ( i = 0; i < N1 ; i++)
            streams(f1,X,0,space)[i*N1+i+space*N2] = f1->sinc.d*(i-N12)*vectorDipole;
    }
    return;
}//build X final
