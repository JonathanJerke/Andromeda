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

double Power ( double b, INT_TYPE n ){
    INT_TYPE i;
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

#if 1
DCOMPLEX cPower ( DCOMPLEX b, INT_TYPE n ){
    INT_TYPE i;
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

DCOMPLEX ovx (  DCOMPLEX Z, DCOMPLEX A,INT_TYPE bb, double kk , double dd, INT_TYPE LL, INT_TYPE b, double k,double d,INT_TYPE L ,INT_TYPE ll){
    INT_TYPE sum = 0;
    if ( (2*ll + LL) >= 0 ){
        if (   -k + kk + (1./d + (2*ll)/(dd*LL))*pi >= 0 )
            sum +=  1;

        if ( -k + kk + (-1./d + (2*ll)/(dd*LL))*pi >= 0 )
            sum -=  1;

    }
    if ( (2*ll - LL) >= 0 ){
        if (   -k + kk + (-1./d + (2*ll)/(dd*LL))*pi >= 0 )
            sum +=  1;

        if ( -k + kk + (1./d + (2*ll)/(dd*LL))*pi >= 0 )
            sum -=   1;

    }
    if ( sum ){
       return  sum*Z * cPower(A, ll);
    } else
        return 0.;
}



#else

DCOMPLEX ovx (  INT_TYPE bb, double kk , double dd, INT_TYPE LL, INT_TYPE b, double k,double d,INT_TYPE L ,INT_TYPE ll){
    INT_TYPE sum = 0;
    if ( (2*ll + LL) >= 0 ){
        if (   -k + kk + (1./d + (2*ll)/(dd*LL))*pi >= 0 )
            sum +=  1;

        if ( -k + kk + (-1./d + (2*ll)/(dd*LL))*pi >= 0 )
            sum -=  1;

    }
    if ( (2*ll - LL) >= 0 ){
        if (   -k + kk + (-1./d + (2*ll)/(dd*LL))*pi >= 0 )
            sum +=  1;

        if ( -k + kk + (1./d + (2*ll)/(dd*LL))*pi >= 0 )
            sum -=   1;

    }
    if ( sum ){
       return  sum*cexp(I*(b*d*k - bb*dd*kk - ((dd*kk*LL + 2*ll*pi)*(b*d - bb*dd))/(dd*LL)))/LL;
    } else
        return 0.;
}
#endif

DCOMPLEX dx ( INT_TYPE N,  INT_TYPE bb, double kk ,double dd,INT_TYPE LL,INT_TYPE ll){
    INT_TYPE nn;
    DCOMPLEX  dx = 1.;
    for (nn= 0; nn < N ; nn++)
        dx *= (dd*kk*LL+2*ll*pi)/dd/LL;
    return dx;
}

DCOMPLEX periodicBoost0 (INT_TYPE N, double P, INT_TYPE bb, double kk ,double dd, INT_TYPE LL,INT_TYPE b, double k ,double d,INT_TYPE L){
    DCOMPLEX  sum = 0.,Z,A;
    INT_TYPE ll = 0;
    Z = cexp(I*(b*d*k - bb*dd*kk - ((dd*kk*LL + 2*ll*pi)*(b*d - bb*dd))/(dd*LL)))/LL;
    ll = 1;
    A = (cexp(I*(b*d*k - bb*dd*kk - ((dd*kk*LL + 2*ll*pi)*(b*d - bb*dd))/(dd*LL)))/LL)/Z;

   // for ( ll = -(LL-1)/2; ll <= (LL-1)/2 ; ll++){
        for ( ll = -LL; ll <= LL ; ll++){
        sum += ovx (Z,A,bb,kk,dd,LL,b,k+P,d,L,ll)*dx(N,bb,kk,dd,LL,ll);
    }
    return sum;
}

DCOMPLEX periodicBoostBasisBasis( INT_TYPE N , double P, INT_TYPE bb, INT_TYPE kki ,double dd, INT_TYPE LL,  INT_TYPE b, INT_TYPE ki,double d,INT_TYPE L ){
    DCOMPLEX sum = 0.;
    double k = 2*pi/L/d * ki;
    double kk = 2*pi/LL/dd * kki;
    if ( ki != 0  && kki != 0 ){
        //b and bb shutdown!!!
        for ( b = 0 ; b < L ; b++)
            for ( bb = 0 ; bb < LL ; bb++){
                sum += ( sign(b)*sign(bb)* ei(-pi/L*1.*ki*b+pi/LL*1.*kki*bb) )  * periodicBoost0(N,P, bb,kk,dd,LL,b,k,d,L) /sqrt(L*LL) ;
            }
    }else if ( ki != 0 )   {
        //b shutdown!!!
        for ( b = 0 ; b < L ; b++)
            sum +=  ( sign(b)* ei(-pi/L*1.*ki*b) )  * periodicBoost0(N,P, bb,kk,dd,LL,b,k,d,L)/sqrt(L) ;

    }else if ( kki != 0 )   {
        //bb  shutdown!!!
        for ( bb = 0 ; bb < LL ; bb++){
            sum +=  periodicBoost0(N,P, bb,kk,dd,LL,b,k,d,L)* ( sign(bb)* ei(pi/LL*1.*kki*bb) /sqrt(LL)) ;
        }
    }else {
        sum += periodicBoost0 (N,P, bb,kk,dd,LL,b,k,d,L);
    }

    return sum;
}
#if 0
DCOMPLEX hyperGeometric (double gamma, INT_TYPE lambda, double delta){
    DCOMPLEX value = 0;
    double lambdad = lambda;
        
    
    double g,hg;
    
    if ( lambda % 2 == 0 ){
        g = gsl_sf_gamma((1.+lambdad)/2.);
        hg = gsl_sf_hyperg_1F1((1.+lambdad)/2.  ,0.5,-sqr(delta/gamma/2.));
    } else {
        g = gsl_sf_gamma(1.+lambdad/2.);
        hg = gsl_sf_hyperg_1F1( 1.+lambdad/2.    ,1.5,-sqr(delta/gamma/2.));
    }

    
    if ( lambda % 2 == 1 )
        value = I *  g*hg   *(delta/gamma);
    else
        value =      g*hg;
    
    value /= pow(gamma,lambdad)*gamma;
    return value;
}
#else
//Bill's addition

DCOMPLEX hyperGeometric (double gamma, INT_TYPE lambda, double delta){
    DCOMPLEX value = 0;
    double lambdad = lambda;
    

    switch (lambda) {
        case 0:
            value = 1;
            break;
        case 1:
            value = I * delta;
            break;
        case 2:
            value = -Power(delta,2) + 2*Power(gamma,2);
            break;
        case 3:
            value = -I * Power(delta,3) - 6*delta*Power(gamma,2);
            break;
        case 4:
            value = Power(delta,4) - 12*Power(delta,2)*Power(gamma,2) + 12*Power(gamma,4);
            break;
        case 5:
            value = I * Power(delta,5) - 20*Power(delta,3)*Power(gamma,2) + 60*delta*Power(gamma,4);
            break;
        case 6:
            value = -Power(delta,6) - 30*Power(delta,4)*Power(gamma,2) + 180*Power(delta,2)*Power(gamma,4) - 120*Power(gamma,6);
            break;
        case 7:
            value = -I*Power(delta,7) - 42*Power(delta,5)*Power(gamma,2) + 420*Power(delta,3)*Power(gamma,4) - 840*delta*Power(gamma,6);
            break;
        case 8:
            value = Power(delta,8) - 56*Power(delta,6)*Power(gamma,2) + 840*Power(delta,4)*Power(gamma,4) - 3360*Power(delta,2)*Power(gamma,6) + 1680*Power(gamma,8);
            break;
        case 9:
            value = I*Power(delta,9) - 72*Power(delta,7)*Power(gamma,2) +
            1512*Power(delta,5)*Power(gamma,4) -
            10080*Power(delta,3)*Power(gamma,6) + 15120*delta*Power(gamma,8);
            break;
        case 10:
            value = -Power(delta,10) - 90*Power(delta,8)*Power(gamma,2) +
            2520*Power(delta,6)*Power(gamma,4) -
            25200*Power(delta,4)*Power(gamma,6) +
            75600*Power(delta,2)*Power(gamma,8) - 30240*Power(gamma,10);
            break;
        case 11:
            value = -I*Power(delta,11) - 110*Power(delta,9)*Power(gamma,2) +
            3960*Power(delta,7)*Power(gamma,4) -
            55440*Power(delta,5)*Power(gamma,6) +
            277200*Power(delta,3)*Power(gamma,8) - 332640*delta*Power(gamma,10);
            break;
        case 12:
            value = Power(delta,12) - 132*Power(delta,10)*Power(gamma,2) +
            5940*Power(delta,8)*Power(gamma,4) -
            110880*Power(delta,6)*Power(gamma,6) +
            831600*Power(delta,4)*Power(gamma,8) -
            1995840*Power(delta,2)*Power(gamma,10) + 665280*Power(gamma,12);
            break;
        case 13:
            value = I*Power(delta,13) - 156*Power(delta,11)*Power(gamma,2) +
            8580*Power(delta,9)*Power(gamma,4) -
            205920*Power(delta,7)*Power(gamma,6) +
            2162160*Power(delta,5)*Power(gamma,8) -
            8648640*Power(delta,3)*Power(gamma,10) + 8648640*delta*Power(gamma,12);
            break;
        case 14:
            value = -Power(delta,14) - 182*Power(delta,12)*Power(gamma,2) +
            12012*Power(delta,10)*Power(gamma,4) -
            360360*Power(delta,8)*Power(gamma,6) +
            5045040*Power(delta,6)*Power(gamma,8) -
            30270240*Power(delta,4)*Power(gamma,10) +
            60540480*Power(delta,2)*Power(gamma,12) - 17297280*Power(gamma,14);
            break;
        case 15:
            value = -I*Power(delta,15) - 210*Power(delta,13)*Power(gamma,2) +
            16380*Power(delta,11)*Power(gamma,4) -
            600600*Power(delta,9)*Power(gamma,6) +
            10810800*Power(delta,7)*Power(gamma,8) -
            90810720*Power(delta,5)*Power(gamma,10) +
            302702400*Power(delta,3)*Power(gamma,12) -
            259459200*delta*Power(gamma,14) ;
            break;
        case 16:
            value = Power(delta,16) - 240*Power(delta,14)*Power(gamma,2) +
            21840*Power(delta,12)*Power(gamma,4) -
            960960*Power(delta,10)*Power(gamma,6) +
            21621600*Power(delta,8)*Power(gamma,8) -
            242161920*Power(delta,6)*Power(gamma,10) +
            1210809600*Power(delta,4)*Power(gamma,12) -
            2075673600*Power(delta,2)*Power(gamma,14) + 518918400*Power(gamma,16);
            break;
        case 17:
            value = I*Power(delta,17) - 272*Power(delta,15)*Power(gamma,2) +
            28560*Power(delta,13)*Power(gamma,4) -
            1485120*Power(delta,11)*Power(gamma,6) +
            40840800*Power(delta,9)*Power(gamma,8) -
            588107520*Power(delta,7)*Power(gamma,10) +
            4116752640*Power(delta,5)*Power(gamma,12) -
            11762150400*Power(delta,3)*Power(gamma,14) +
            8821612800*delta*Power(gamma,16) ;
            break;
        case 18:
            value = -Power(delta,18) - 306*Power(delta,16)*Power(gamma,2) +
            36720*Power(delta,14)*Power(gamma,4) -
            2227680*Power(delta,12)*Power(gamma,6) +
            73513440*Power(delta,10)*Power(gamma,8) -
            1323241920*Power(delta,8)*Power(gamma,10) +
            12350257920*Power(delta,6)*Power(gamma,12) -
            52929676800*Power(delta,4)*Power(gamma,14) +
            79394515200*Power(delta,2)*Power(gamma,16) -
            17643225600*Power(gamma,18);
            break;
        case 19:
            value = -I*Power(delta,19) - 342*Power(delta,17)*Power(gamma,2) +
            46512*Power(delta,15)*Power(gamma,4) -
            3255840*Power(delta,13)*Power(gamma,6) +
            126977760*Power(delta,11)*Power(gamma,8) -
            2793510720*Power(delta,9)*Power(gamma,10) +
            33522128640*Power(delta,7)*Power(gamma,12) -
            201132771840*Power(delta,5)*Power(gamma,14) +
            502831929600*Power(delta,3)*Power(gamma,16) -
            335221286400*delta*Power(gamma,18);
            break;
        case 20:
            value = Power(delta,20) - 380*Power(delta,18)*Power(gamma,2) +
            58140*Power(delta,16)*Power(gamma,4) -
            4651200*Power(delta,14)*Power(gamma,6) +
            211629600*Power(delta,12)*Power(gamma,8) -
            5587021440*Power(delta,10)*Power(gamma,10) +
            83805321600*Power(delta,8)*Power(gamma,12) -
            670442572800*Power(delta,6)*Power(gamma,14) +             2514159648000*Power(delta,4)*Power(gamma,16) -
            3352212864000*Power(delta,2)*Power(gamma,18) +
            670442572800*Power(gamma,20);
            break;
        case 21:
            value = I * Power(delta,21) - 420*Power(delta,19)*Power(gamma,2) +
            71820*Power(delta,17)*Power(gamma,4) -
            6511680*Power(delta,15)*Power(gamma,6) +
            341863200*Power(delta,13)*Power(gamma,8) -
            10666131840*Power(delta,11)*Power(gamma,10) +
            195545750400*Power(delta,9)*Power(gamma,12) -
            2011327718400*Power(delta,7)*Power(gamma,14) +
            10559470521600*Power(delta,5)*Power(gamma,16) -
            23465490048000*Power(delta,3)*Power(gamma,18) +
            14079294028800*delta*Power(gamma,20);
            break;
        case 22:
            value = -Power(delta,22) - 462*Power(delta,20)*Power(gamma,2) +
            87780*Power(delta,18)*Power(gamma,4) -
            8953560*Power(delta,16)*Power(gamma,6) +
            537213600*Power(delta,14)*Power(gamma,8) -
            19554575040*Power(delta,12)*Power(gamma,10) +
            430200650880*Power(delta,10)*Power(gamma,12) -
            5531151225600*Power(delta,8)*Power(gamma,14) +
            38718058579200*Power(delta,6)*Power(gamma,16) -             129060195264000*Power(delta,4)*Power(gamma,18) +
            154872234316800*Power(delta,2)*Power(gamma,20) -
            28158588057600*Power(gamma,22);
            break;
        case 23:
            value = -I*Power(delta,23) - 506*Power(delta,21)*Power(gamma,2) +
            106260*Power(delta,19)*Power(gamma,4) -
            12113640*Power(delta,17)*Power(gamma,6) +
            823727520*Power(delta,15)*Power(gamma,8) -
            34596555840*Power(delta,13)*Power(gamma,10) +
            899510451840*Power(delta,11)*Power(gamma,12) -
            14135164243200*Power(delta,9)*Power(gamma,14) +
            127216478188800*Power(delta,7)*Power(gamma,16) -
            593676898214400*Power(delta,5)*Power(gamma,18) +
            1187353796428800*Power(delta,3)*Power(gamma,20) -
            647647525324800*delta*Power(gamma,22);
            break;
        case 24:
            value = Power(delta,24) - 552*Power(delta,22)*Power(gamma,2) + 127512*Power(delta,20)*Power(gamma,4) -
            16151520*Power(delta,18)*Power(gamma,6) +
            1235591280*Power(delta,16)*Power(gamma,8) -
            59308381440*Power(delta,14)*Power(gamma,10) +
            1799020903680*Power(delta,12)*Power(gamma,12) -
            33924394183680*Power(delta,10)*Power(gamma,14) +
            381649434566400*Power(delta,8)*Power(gamma,16) -
            2374707592857600*Power(delta,6)*Power(gamma,18) +
            7124122778572800*Power(delta,4)*Power(gamma,20) -
            7771770303897600*Power(delta,2)*Power(gamma,22) +
            1295295050649600*Power(gamma,24);
            break;
    }

    value *= 1.7724538509055159/*sqrt(pi)*/*exp(-sqr( delta/2./gamma))/(pow(gamma,2*lambdad)*gamma)/pow(2.,lambdad);

    return value;
}



#endif


INT_TYPE maxZ( struct input * f1 ){
    INT_TYPE maxz = 0;
    INT_TYPE a,ai;
    for ( a= 1 ;a <= f1->Na ; a++)
        if ( maxz < f1->atoms[a].label.Z){
            maxz=f1->atoms[a].label.Z;
            ai = a;
        }
    return maxz;
}

INT_TYPE sumZ( struct input * f1 ){
    INT_TYPE a,ai=0;
    for ( a= 1 ;a <= f1->Na ; a++)
        ai +=  f1->atoms[a].label.Z;
    return ai;
}
