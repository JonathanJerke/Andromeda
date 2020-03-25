/*
 *  eri.c
 *
 *
 *  Copyright 2020 Jonathan Jerke and Bill Poirier.
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

#if 0
int main (int argc , char * argv[]){


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
#else

int main (int argc , char * argv[]){
    
    INT_TYPE fi,cmpl,lines = 0,num;
    size_t ms = MAXSTRING;
    char line0[MAXSTRING];
    char* line = line0;
    double pos;
    
    
    
    INT_TYPE n,i,j,N1;
    struct calculation c;
    struct field f;
    c = initCal();
    f = initField();

    struct runTime * rt = & c.rt;
    f.f.rt = rt;
    FILE * outME = fopen(argv[2],"w");
    FILE * outSINC = fopen(argv[3],"w");

    FILE * in = fopen(argv[1],"r");
    getline(&line, &ms, in);

    while ( ! feof(in) ) {
        sscanf(line, "%lf", &pos);
        c.i.Na = 1;
        c.i.atoms[1].label.Z = 1;
        c.i.atoms[1].position[1] = pos;
        printf("%f\n", pos);
        iModel(&c, &f);
        N1 = vectorLen(f.f, 0);
        double ar[N1*N1];

        for ( n = 0; n< N1*N1 ; n++)
            ar[n] = 0.;
        for ( n = 0; n < CanonicalRank(f.f, linear, 0) ; n++)
            cblas_daxpy(N1*N1, 1., streams(f.f,linear,0,0)+N1*N1*n, 1, ar, 1);
        
        
        for ( i = 0; i < N1 ; i++){
            for ( j = 0; j < N1 ; j++)
            {
                if ( i || j )
                    fprintf(outME,",");
                fprintf(outME,"%f", ar[i*N1+j]);
            }
            if ( i )
                fprintf(outSINC,",");
            fprintf(outSINC,"%f", Sinc(1, pos/f.i.d-(i-(N1-1)/2)));
        }
        fModel( &f.f);
        fprintf(outME, "\n");
        fprintf(outSINC, "\n");
        getline(&line, &ms, in);

    }
    fclose(outME);
    fclose(outSINC);
    fclose(in);
    
    return 0;
}
#endif
