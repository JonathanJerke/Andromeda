/*
 *  eri.c
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


#include "eri.h"


int main (INT_TYPE argc , char * argv[]){


    struct function_label func1 ;
    
    func1.fn = Coulomb;//Yukawa, Erf
    func1.param[0] = 1.;//scalar out front
    func1.param[1] = 1.;//not used by elemCal...a quadrature parameter
    func1.param[2] = 1.;//parameter in Erf and Yukawa (scale/mass)
    func1.param[3] = 1.;//second a parameter, unused.
    func1.interval = 1;
    getDescription ( &func1 ,1.,stdout);
    
    INT_TYPE i ;
    INT_TYPE l;
    double x,b;
    struct general_2index g3[3];

#ifdef APPLE
    
    b = 1;
    x = 0;
    l = 2;
    
    
#else
    if ( argc != 4 ){
        printf("L,Distance, Gaussian\n");
        exit(0);
        
    }
     l = atoi(argv[1]);
    x = atof(argv[2]);
     b = atof(argv[3]);
    
    printf("ANGULAR %d \n %d-Body interaction", l , 2);
#endif
    
    {
        i = 0;
        g3[i].gaussianAccelerationFlag = 1;
        g3[i].point = 2;

//        defineGaussBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, ble[bl], 0.,             f1->rose[space].count1Basis,0)
    
        g3[i].i[0].bra = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, x, 0, l);
        g3[i].i[0].ket = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, x, 0, l);
        
        g3[i].i[1].bra = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, 0, 0, l);
        g3[i].i[1].ket = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, 0, 0, l);

        g3[i].fl = & func1;
        
    }
    {
        i = 1;
        g3[i].gaussianAccelerationFlag = 1;
        g3[i].point = 2;

        //        defineGaussBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, ble[bl], 0.,             f1->rose[space].count1Basis,0)
        
        g3[i].i[0].bra = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, x, 0, l);
        g3[i].i[0].ket = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, x, 0, l);
        
        g3[i].i[1].bra = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, 0, 0, l);
        g3[i].i[1].ket = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, 0, 0, l);
        
        g3[i].fl = & func1;
        
    }
    {
        i = 2;
        g3[i].gaussianAccelerationFlag = 1;
        g3[i].point = 2;

        //        defineGaussBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, ble[bl], 0.,             f1->rose[space].count1Basis,0)
        
        g3[i].i[0].bra = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, x, 0, l);
        g3[i].i[0].ket = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, x, 0, l);
        
        g3[i].i[1].bra = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, 0, 0, l);
        g3[i].i[1].ket = defineGaussBasis(nullNote, 1, GaussianBasisElement, b, 0, 0, l);
        
        g3[i].fl = & func1;
        
    }
    printf(" %1.15f\n", quadCal( g3));
    if(1){
        INT_TYPE nn;
    
    for ( nn = 0; nn < 10000; nn++)
        quadCal( g3);
    }
    else{
        INT_TYPE nn;
        
        INT_TYPE k = 0;
        for ( nn = 0; nn < 1000000000; nn++)
            k+= nn;
    }

}
