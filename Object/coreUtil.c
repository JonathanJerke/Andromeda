/**
 *  coreUtil.c
 *
 *
 *  Copyright 2021 Jonathan Jerke and Bill Poirier.
 *  Ongoing support for this program is coordinated through quantumgalaxies.org.
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

#include "coreUtil.h"


division name (   sinc_label f1,   division label){
    return f1.name[f1.name[f1.name[f1.name[label].name].name].name].name;
}


inta part (   sinc_label f1 ,   division label){
    
    if ( label > f1.end || label < 0 ){
        printf("part, accessed name is outside the bounds\n");
    }
    return f1.name[name(f1,label)].Partition;
}

inta species (   sinc_label f1 ,   division label ){
      division u = name(f1,label);
    return f1.name[u].species;
}

bodyType Bodies (   sinc_label f1 ,   division label,inta space ){
    if ( species ( f1, label ) == vector )
        return f1.canon[space].body ;
    else
        return f1.name[name(f1,label)].space[space].body;
        
}

bodyType bodies (   sinc_label f1 ,   division label ){
    bodyType x=nada,u;
    inta space;
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.canon[space].body != nada)
    {
        u =Bodies(f1, label, space);
        if ( u > x )
            x = u ;
    }
    return x;
}


inta header (   sinc_label f1 ,   division label ){
    return 0;
}

inta length (   sinc_label f1 ,   division label, inta *lens ){
    
    inta space ;
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.canon[space].body != nada )
            lens[space] = alloc(f1,label ,space);
    return 0;
}

inta vectorLen(  sinc_label f1, inta space){
    if ( f1.canon[space].body == one )
        return f1.canon[space].count1Basis ;
    else if (f1.canon[space].body == two )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis  ;
    else if ( f1.canon[space].body == three )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis  ;
    else if ( f1.canon[space].body == four )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis  ;
    else
        return 0;
}

inta outerVectorLen(  sinc_label f1,   bodyType bd, inta space){
    if ( bd == one )
        return f1.canon[space].count1Basis ;
    else if ( bd == two )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis  ;
    else if ( bd == three )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis  ;
    else if ( bd == four )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis  ;
    else
        return 0;
}


inta vector1Len(  sinc_label f1, inta space){
    if ( f1.canon[space].body != nada )

        return f1.canon[space].count1Basis ;
    else
        return 0;
}

void length1(  sinc_label f1, inta *len){
    inta space ;
    for ( space = 0 ;space < SPACE;space++)
        if ( f1.canon[space].body != nada )
            len[space] = vector1Len(f1, space);
    
    return;
}

  division anotherLabel( sinc_label *f1,   inta particle,  bodyType body){
    switch(body){
        case nada:
            if ( f1->nullLabels.currLabel+1 > f1->nullLabels.maxLabel)
            {
                printf("add more 'names' %d\n", f1->nullLabels.maxLabel);
                fflush(stdout);
                exit(0);
            }
            else{
                division output = f1->nullLabels.head + f1->nullLabels.currLabel++;
                f1->name[output].linkNext = nullName;
                f1->name[output].chainNext = nullName;
                f1->name[output].loopNext = nullName;
                f1->name[output].multNext = nullName;

                f1->name[output].Current[0] = 0;
                f1->name[output].Current[1] = 0;
                f1->name[output].name = output;
                return output;
            }
        case one:
        case two:
            if ( f1->eikonLabels.currLabel+1 > f1->eikonLabels.maxLabel)
            {
                printf("add more 'eikons' %d\n",f1->eikonLabels.maxLabel);
                fflush(stdout);

                exit(0);
            }
            else{
                division output = f1->eikonLabels.head + f1->eikonLabels.currLabel++;
                f1->name[output].linkNext = nullName;
                f1->name[output].chainNext = nullName;
                f1->name[output].loopNext = nullName;
                f1->name[output].multNext = nullName;
                f1->name[output].Current[0] = 0;
                f1->name[output].Current[1] = 0;
                f1->name[output].name = output;
                return output;
            }
        default:
            break;

    }
    printf ("out of names\n");
    fflush(stdout);

    exit(0);
    return nullName;
};



inta matrixLen(  sinc_label f1,   bodyType body, inta space){
    if ( body == one )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis ;
    else if ( body == two )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis  ;
    else if ( body == three )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis  ;
    else if ( body == four )
        return f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis*f1.canon[space].count1Basis  ;
    else
        return 0;
}

inta topezOp(double origin, double lattice,  bodyType bd,inta act,   blockType cl,   blockType bl,  inta N1,floata * vector , inta pw, floata * vectorOut){
    ///COULD ALWAYS TRANSLATE WITH TRANSLATION OPERATOR
    inta n,m,m1,m2,m3;
    double sign = 1.,mult,sign1 =1.;
    if ( pw == 1 )
        sign1 = -1.;
    inta n1[7];
    inta perm[7];
    inta op[7];
    commandSA(bd,act,cl,bl,perm,op);
    
    switch(bd){
        case one:
            n1[0] = 1;
            for ( n = 0 ; n < N1 ; n++)
                vectorOut[n] = 0.;
        break;
        case two:
            n1[0] = 1;
            n1[1] = N1;
            for ( n = 0 ; n < N1*N1 ; n++)
                vectorOut[n] = 0.;
            break;
        case three:
            n1[0] = 1;
            n1[1] = N1;
            n1[2] = N1*N1;
            for ( n = 0 ; n < N1*N1*N1; n++)
                vectorOut[n] = 0.;
            break;
        case four:
            n1[0] = 1;
            n1[1] = N1;
            n1[2] = N1*N1;
            n1[3] = N1*N1*N1;
            for ( n = 0 ; n < N1*N1*N1*N1; n++)
                vectorOut[n] = 0.;
            break;

        default:
            break;

    }


    
    switch ( bd ){
        case one:
            //  vector -> vectorOut
            sign = 1.;
            mult = 0.;
            if ( pw == 0 ){
                cblas_dcopy(N1, vector, 1, vectorOut, 1);
            }else if ( pw == -1 ){
                for ( n = 0 ; n < N1 ; n++)
                    vectorOut[n] = (n*lattice-origin) * vector[n] ;
            }
            else if ( pw == -2 ){
                for ( n = 0 ; n < N1 ; n++)
                    vectorOut[n] = (n*lattice-origin) * (n*lattice-origin) * vector[n] ;
            }
            else {
                if (pw == 2 ){
                    cblas_dcopy(N1, vector, 1, vectorOut, 1);
                    cblas_dscal(N1, -pi*pi/3./lattice/lattice, vectorOut, 1);
                }

                for (n = 1 ; n < N1 ; n++){
                    if (pw == 1 ){
                        mult = 1./n/lattice;
                    }else if ( pw == 2){
                        mult = 2*1./(n*n)/lattice/lattice;
                    }
                    cblas_daxpy(N1-n, sign*mult, vector+n, 1, vectorOut, 1);
                    cblas_daxpy(N1-n, sign1*sign*mult, vector, 1, vectorOut+n,1);
                    sign *= -1;

                }
            }
            break;
        case two:
            if ( pw != 1 )
            for ( m = 0; m < N1 ; m++)
            {
                if ( pw == 0 ){
                    cblas_dcopy(N1, vector+n1[perm[op[1]]]*m, n1[perm[op[0]]], vectorOut+n1[op[1]]*m, n1[op[0]]);//must copy from perm
                }else if ( pw == -1 ){
                    for ( n = 0 ; n < N1 ; n++)
                        vectorOut[n*n1[op[0]]+n1[op[1]]*m] = (n*lattice-origin) *vector[n*n1[perm[op[0]]]+n1[perm[op[1]]]*m] ;
                }
                else if ( pw == -2 ){
                    for ( n = 0 ; n < N1 ; n++)
                        vectorOut[n*n1[op[0]]+n1[op[1]]*m] = (n*lattice-origin) *(n*lattice-origin) *vector[n*n1[perm[op[0]]]+n1[perm[op[1]]]*m] ;
                }
                else {
                        if (pw == 2) {
                            cblas_dcopy(N1, vector+n1[perm[op[1]]]*m, n1[perm[op[0]]], vectorOut+n1[op[1]]*m, n1[op[0]]);
                             cblas_dscal(N1, -pi*pi/3./lattice/lattice, vectorOut+n1[op[1]]*m, n1[op[0]]);
                        }
                }
            }
            
            
            if (pw == 2 || pw == 1 ) {
                for ( m = 0; m < N1 ; m++)
                {
                    sign = 1.;
                    mult = 0;
                    for (n = 1 ; n < N1 ; n++){
                        if (pw == 1 ){
                            mult = 1./n/lattice;
                        }else if ( pw == 2){
                            mult = 2*1./(n*n)/lattice/lattice;
                        }
                        cblas_daxpy(N1-n, sign*mult      , vector+n1[perm[op[0]]]*n+n1[perm[op[1]]]*m ,n1[perm[op[0]]], vectorOut+n1[op[1]]*m ,n1[op[0]]);
                        cblas_daxpy(N1-n, sign1*sign*mult, vector+n1[perm[op[1]]]*m , n1[perm[op[0]]], vectorOut+n1[op[0]]*n+n1[op[1]]*m ,n1[op[0]]);

                        sign *= -1;
                    }
                }
            }
            
            
            break;
        case three:
            if ( pw != 1 )

            for ( m = 0; m < N1 ; m++)
                for ( m2 = 0; m2 < N1 ; m2++)
                {
                    if ( pw == 0 ){
                        cblas_dcopy(N1, vector+n1[perm[op[1]]]*m+n1[perm[op[2]]]*m2, n1[perm[op[0]]], vectorOut+n1[op[1]]*m+n1[op[2]]*m2, n1[op[0]]);
                    }else if ( pw == -1 ){
                        for ( n = 0 ; n < N1 ; n++)
                            vectorOut[n*n1[op[0]]+n1[op[1]]*m+n1[op[2]]*m2] = (n*lattice-origin) * vector[n*n1[perm[op[0]]]+n1[perm[op[1]]]*m+n1[perm[op[2]]]*m2] ;
                    }
                    else if ( pw == -2 ){
                        for ( n = 0 ; n < N1 ; n++)
                            vectorOut[n*n1[op[0]]+n1[op[1]]*m+n1[op[2]]*m2] = (n*lattice-origin) * (n*lattice-origin) * vector[n*n1[perm[op[0]]]+n1[perm[op[1]]]*m+n1[perm[op[2]]]*m2] ;
                    }
                    else {
                        if (pw == 2) {
                             cblas_dcopy(N1, vector+n1[perm[op[1]]]*m+n1[perm[op[2]]]*m2, n1[perm[op[0]]], vectorOut+n1[op[1]]*m+n1[op[2]]*m2, n1[op[0]]);
                             cblas_dscal(N1, -pi*pi/3./lattice/lattice, vectorOut+n1[op[1]]*m+n1[op[2]]*m2, n1[op[0]]);
                         }
                     }
                }
            
            if (pw == 2 || pw == 1  ) {
                for ( m = 0; m < N1 ; m++)
                    for ( m2 = 0; m2 < N1 ; m2++)
                    {
                        sign = 1.;
                        mult = 0;

                            for (n = 1 ; n < N1 ; n++){
                                if (pw == 1 ){
                                    mult = 1./n/lattice;
                                }else if ( pw == 2){
                                    mult = 2*1./(n*n)/lattice/lattice;
                                }

                                       cblas_daxpy(N1-n, sign*mult, vector+n1[perm[op[0]]]*n+n1[perm[op[1]]]*m +n1[perm[op[2]]]*m2 ,n1[perm[op[0]]], vectorOut+n1[op[1]]*m +n1[op[2]]*m2,n1[op[0]]);
                                       cblas_daxpy(N1-n, sign1*sign*mult, vector+n1[perm[op[1]]]*m+n1[perm[op[2]]]*m2 ,n1[perm[op[0]]], vectorOut+n1[op[0]]*n +n1[op[1]]*m+n1[op[2]]*m2 ,n1[op[0]]);
                                    sign *= -1;

                                   }
                    }
            }
            break;
        case four:
            if ( pw != 1 )

            for ( m = 0; m < N1 ; m++)
                for ( m2 = 0; m2 < N1 ; m2++)
                    for ( m3 = 0; m3 < N1 ; m3++)

                {
                    if ( pw == 0 ){
                        cblas_dcopy(N1, vector+n1[perm[op[1]]]*m+n1[perm[op[2]]]*m2+n1[op[3]]*m3, n1[perm[op[0]]], vectorOut+n1[op[1]]*m+n1[op[2]]*m2+n1[op[3]]*m3, n1[op[0]]);
                    }else if ( pw == -1 ){
                        for ( n = 0 ; n < N1 ; n++)
                            vectorOut[n*n1[op[0]]+n1[op[1]]*m+n1[op[2]]*m2+n1[op[3]]*m3] = (n*lattice-origin) * vector[n*n1[perm[op[0]]]+n1[perm[op[1]]]*m+n1[perm[op[2]]]*m2+n1[op[3]]*m3] ;
                    }
                    else if ( pw == -2 ){
                        for ( n = 0 ; n < N1 ; n++)
                            vectorOut[n*n1[op[0]]+n1[op[1]]*m+n1[op[2]]*m2+n1[op[3]]*m3] = (n*lattice-origin) * (n*lattice-origin) * vector[n*n1[perm[op[0]]]+n1[perm[op[1]]]*m+n1[perm[op[2]]]*m2+n1[op[3]]*m3] ;
                    }
                    else {
                        if (pw == 2) {
                             cblas_dcopy(N1, vector+n1[perm[op[1]]]*m+n1[perm[op[2]]]*m2+n1[op[3]]*m3, n1[perm[op[0]]], vectorOut+n1[op[1]]*m+n1[op[2]]*m2+n1[op[3]]*m3, n1[op[0]]);
                             cblas_dscal(N1, -pi*pi/3./lattice/lattice, vectorOut+n1[op[1]]*m+n1[op[2]]*m2+n1[op[3]]*m3, n1[op[0]]);
                         }
                     }
                }
            
            if (pw == 2 || pw == 1 ) {
                for ( m1 = 0; m1 < N1 ; m1++)
                    for ( m2 = 0; m2 < N1 ; m2++)
                        for ( m3 = 0; m3 < N1 ; m3++)
                {
                        sign = 1.;
                        mult = 0;

                            for (n = 1 ; n < N1 ; n++){
                                if (pw == 1 ){
                                    mult = 1./n/lattice;
                                }else if ( pw == 2){
                                    mult = 2*1./(n*n)/lattice/lattice;
                                }

                                       cblas_daxpy(N1-n, sign*mult, vector+n1[perm[op[0]]]*n+n1[perm[op[1]]]*m1 +n1[perm[op[2]]]*m2+n1[op[3]]*m3 ,n1[perm[op[0]]], vectorOut+n1[op[1]]*m1 +n1[op[2]]*m2+n1[op[3]]*m3,n1[op[0]]);
                                       cblas_daxpy(N1-n, sign1*sign*mult, vector+n1[perm[op[1]]]*m1+n1[perm[op[2]]]*m2+n1[op[3]]*m3 ,n1[perm[op[0]]], vectorOut+n1[op[0]]*n +n1[op[1]]*m1+n1[op[2]]*m2 +n1[op[3]]*m3,n1[op[0]]);

                                    sign *= -1;

                                   }
                    }
            }
            break;
        default:
            break;

    }
        







    
//    if (bd == one ){
//
//
//        for (n = 1 ; n < N1 ; n++){
//            if (pw == 1 ){
//                mult = 1./n;
//            }else {
//                mult = 2*1./(n*n);
//            }
//            cblas_daxpy(N1-n, sign*mult, vector+n, 1, vectorOut, 1);
//            sign *= -1;
//            cblas_daxpy(N1-n, sign*mult, vector, 1, vectorOut+n,1);
//        }
//
//
//        //shift cblas_daxpy
//        return 0;
//        //with tolerances
//    }else if ( bd == two ){
//        for ( n = 0 ; n < N1*N1 ; n++)
//            vectorOut[n] = 0.;
//
//        switch( tv){
//            case tv1:
//                //shift cblas_daxpy
//                for ( m = 0 ; m < N1 ; m++){
//                    sign = 1.;
//                    for (n = 1 ; n < N1 ; n++){
//                        cblas_daxpy(N1-n, sign/n, vector+n+m*N1,1, vectorOut+m*N1, 1);
//                        sign *= -1;        //         m                   n
//                        cblas_daxpy(N1-n, sign/n, vector+m*N1, 1, vectorOut+n+m*N1,1);
//                    }
//                }
//                //with tolerances
//                return 0;
//            case tv2:
//                    //shift cblas_daxpy
//                    for (n = 1 ; n < N1 ; n++){
//                        cblas_daxpy(N1*N1-N1*n, sign/n, vector+n*N1, 1, vectorOut, 1);
//                        sign *= -1;
//                        cblas_daxpy(N1*N1-N1*n, sign/n,vector, 1,      vectorOut+n*N1,1);
//                    }
//                    //with tolerances
//                    return 0;
//
//
//        }
//
//
//    }
    return 0;
}

inta diagonalOp(  bodyType bd,  inta act,   blockType cl,   blockType bl, inta N1,floata * vector, floata * toep, floata* vectorOut){
    //COULD ALWAYS TRANSLATE WITH TRANSLATION OPERATOR
    inta perm[7],m,m2,m3,n;
    inta n1[7];
    inta op[7];
    switch (commandSA(bd,act,cl,bl,perm,op)){
        case one://OP
            switch(bd){
                case one:
                    n1[0] = 1;
                    cblas_dcopy(N1, vector, 1, vectorOut, 1);
                    cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,N1,0,toep,1, vectorOut,1);
                    return 0;
                case two:
                    n1[0] = 1;
                    n1[1] = N1;
                    for ( m = 0 ; m < N1 ; m++){
                        cblas_dcopy(N1, vector+n1[perm[op[1]]]*m, n1[perm[op[0]]], vectorOut+n1[op[1]]*m, n1[op[0]]);
                        cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,N1,0,toep,1, vectorOut+m*n1[op[1]],n1[op[0]]);
                    }
                    return 0;
                case three:
                    n1[0] = 1;
                    n1[1] = N1;
                    n1[2] = N1*N1;
                    for ( m = 0 ; m < N1 ; m++)
                        for ( m2 = 0 ; m2 < N1 ; m2++){
                            cblas_dcopy(N1, vector+m*n1[perm[op[1]]]+m2*n1[perm[op[2]]],n1[perm[op[0]]], vectorOut+n1[op[1]]*m+n1[op[2]]*m2, n1[op[0]]);
                            cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,N1,0,toep,1, vectorOut+n1[op[1]]*m+n1[op[2]]*m2, n1[op[0]]);
                        }
                    return 0;
                case four:
                    n1[0] = 1;
                    n1[1] = N1;
                    n1[2] = N1*N1;
                    n1[3] = N1*N1*N1;
                    for ( m = 0 ; m < N1 ; m++)
                        for ( m2 = 0 ; m2 < N1 ; m2++)
                            for ( m3 = 0 ; m3 < N1 ; m3++){
                            cblas_dcopy(N1, vector+m*n1[perm[op[1]]]+m2*n1[perm[op[2]]]+m3*n1[perm[op[3]]],n1[perm[op[0]]], vectorOut+n1[op[1]]*m+n1[op[2]]*m2+m3*n1[perm[op[3]]], n1[op[0]]);
                            cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,N1,0,toep,1, vectorOut+n1[op[1]]*m+n1[op[2]]*m2+m3*n1[perm[op[3]]], n1[op[0]]);
                        }
                return 0;
                    default:
                        break;

            }
        case two://OP
            switch(bd){
                case one:
                    printf("diagonalOp, not compatible with one!");
                    exit(0);
        //        case two:
//                    FOR EQUAL one == two...  but simply unncessary speed...
//                    n1[0] = 1;
//                    n1[1] = N1;
//
//                    for ( m = 0 ; m < N1 ; m++)
//                        cblas_dcopy(N1, vector+m*n1[perm[op[1]]],n1[perm[op[0]]], vectorOut+m*n1[op[1]],n1[op[0]]);
//
//                    //below is 2-body
//                    cblas_dscal(N1, toep[0], vectorOut, n1[op[0]]+n1[op[1]]);
//                    for (n = 1 ; n < N1 ; n++){
//                        {
//                            cblas_dscal(N1-n, toep[n], vectorOut+n1[op[0]]*n,  n1[op[0]]+n1[op[1]]);
//                            cblas_dscal(N1-n, toep[-n], vectorOut+n1[op[1]]*n,  n1[op[0]]+n1[op[1]]);
//                        }
//                    }
//                    return 0;
                    case two:
                        n1[0] = 1;
                        n1[1] = N1;
                        for ( n = 0 ; n < N1 ; n++){
                            cblas_dcopy(N1, vector+n1[perm[op[1]]]*n, n1[perm[op[0]]], vectorOut+n1[op[1]]*n, n1[op[0]]);
                                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,N1,0,toep+n*N1,1, vectorOut+n*n1[op[1]],n1[op[0]]);

                        }
                        return 0;
                
                case three:
                    n1[0] = 1;
                    n1[1] = N1;
                    n1[2] = N1*N1;
                    for ( m = 0 ; m < N1 ; m++)
                        for ( n = 0 ; n < N1 ; n++){
                            cblas_dcopy(N1, vector+n1[perm[op[1]]]*n+m*n1[perm[op[2]]], n1[perm[op[0]]], vectorOut+n1[op[1]]*n+m*n1[op[2]], n1[op[0]]);
                                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,N1,0,toep+n*N1,1, vectorOut+n*n1[op[1]]+m*n1[op[2]],n1[op[0]]);

                        }
                    return 0;
                    
                case four:
                    n1[0] = 1;
                    n1[1] = N1;
                    n1[2] = N1*N1;
                    n1[3] = N1*N1*N1;
                    for ( m3 = 0 ; m3 < N1 ; m3++)
                    for ( m = 0 ; m < N1 ; m++)
                        for ( n = 0 ; n < N1 ; n++){
                            cblas_dcopy(N1, vector+n1[perm[op[1]]]*n+m*n1[perm[op[2]]]+m3*n1[perm[op[3]]], n1[perm[op[0]]], vectorOut+n1[op[1]]*n+m*n1[op[2]]+m3*n1[perm[op[3]]], n1[op[0]]);
                                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,N1,0,toep+n*N1,1, vectorOut+n*n1[op[1]]+m*n1[op[2]]+m3*n1[perm[op[3]]],n1[op[0]]);

                        }
                    return 0;

                default:
                    break;

            }
        default:
            break;

    }
    return 0;
}


inta InvertOp(  bodyType bd, inta invert,inta N1,floata * vector, floata* vectorOut){
    ///INVERT PARTICLE 1
    inta m,m2,n;
    inta n1[7];
    switch(bd){
        case one:
            n1[0] = 1;
            for ( n = 0; n < N1 ;n++)
                vectorOut[N1-n-1] = vector[n];
            return 0;
        case two:
            switch( invert){
                case 1:
                    n1[0] = 1;
                    n1[1] = N1;
                    break;
                case 2:
                    n1[0] = N1;
                    n1[1] = 1;
                    break;
            }
            for ( m = 0 ; m < N1 ; m++){
                for ( n = 0; n < N1 ;n++)
                    vectorOut[(N1-n-1)*n1[0]+m*n1[1]] = vector[n*n1[0]+m*n1[1]];
            }
            return 0;
        case three:
            switch( invert){
                case 1:
                    n1[0] = 1;
                    n1[1] = N1;
                    n1[2] = N1*N1;
                    break;
                case 2:
                    n1[0] = N1;
                    n1[1] = 1;
                    n1[2] = N1*N1;
                    break;
                case 3:
                    n1[0] = N1*N1;
                    n1[1] = N1;
                    n1[2] = 1;
                    break;

            }
            for ( m = 0 ; m < N1 ; m++)
                for ( m2 = 0 ; m2 < N1 ; m2++){
                    for ( n = 0; n < N1 ;n++)
                        vectorOut[(N1-n-1)*n1[0]+m*n1[1]+m2*n1[2]] = vector[n*n1[0]+m*n1[1]+m2*n1[2]];
                }
        return 0;
    default:
            break;

    }
    return 0;
}


inta alloc (   sinc_label f1 ,   division label ,inta space){
    if  ( space <0 || space > SPACE ){
        printf("alloc, space outside bounds\n");
        exit(7);
    }
        if ( space == SPACE ){
        return 1;
        }else if (f1.canon[space].body == nada ){
            return 0;
        }else if ( species(f1, label ) == vector ){
            return vectorLen(f1, space);
        } else if ( species(f1, label ) == matrix ){
            return matrixLen(f1, f1.name[name(f1,label)].space[space].body, space);
        }else if ( species(f1, label ) == scalar ){
            return 1;
        }else if ( species(f1, label ) == outerVector ){
            return outerVectorLen(f1, f1.name[name(f1,label)].space[space].body, space);
        }else if ( species(f1, label ) >= eikon ){
//            if ( species(f1, label ) >= eikonKinetic ){
//                return 1;
//            }
            if (f1.name[name(f1,label)].space[space].body == one )
                return outerVectorLen(f1, two, space);
            else if (f1.name[name(f1,label)].space[space].body == two)
                return outerVectorLen(f1, two, space);
        }
    
        return 1;
}
inta pZero (   sinc_label * f1 ,   division label, inta spin ){
    return zero(*f1,label,spin);
}
inta zero (   sinc_label f1 ,   division label, inta spin ){
    //f1.name[label].Current = 0;
    inta i, space,M2[SPACE];
    length(f1, label, M2);
    
    
    for ( space = 0; space < SPACE ;space++)
        if ( f1.canon[space].body != nada){
        floata * pt = streams(f1,label, (spin) , space );

        for ( i = 0; i < M2[space]*part(f1, label) ; i++ ){
            *(pt+i) = 0.;
        }
    }
    return 0 ;
}

inta myZero (   sinc_label f1 ,   division label, inta spin ){
    //f1.name[label].Current = 0;
    inta i;

    for ( i = 0; i < part(f1, label) ; i++ ){
        myStreams(f1,label,spin)[i] = 0.;
    }
    return 0 ;
}
inta pClear (   sinc_label *f1 ,   division label ){
    return tClear(*f1,label);
}

inta tClear (   sinc_label f1 ,   division label ){
    inta spin ;
    for ( spin = 0; spin < MAX_CORE ; spin++){
        f1.name[label].Current[spin] = 0;
    }
    return 0 ;
}

//double volume (   input * f1 ){
//    return pow( vectorLen(f1->c.sinc,0)*f1->d,3. );
//}

inta CanonicalRank(   sinc_label f1 ,   division label , inta spin ){
    if ( label > f1.end ){
        printf("CanonicalRank, name outside bounds\n");
        fflush(stdout);
        exit(0);
    }
    
    if ( f1.name[label].name == label){
        if ( spin < spins(f1, label) ){
            if ( f1.name[f1.name[label].name].species >= eikon){
                return 1;
            } else {
                return f1.name[label].Current[spin];
            }
        }
            else
                return 0;
        }
    else {
//        printf("nameless %d %d %d\n", label, name(f1,label),f1.name[label].Current[spin]);
        return f1.name[label].Current[spin];
    }
}

inta CanonicalOperator( sinc_label f1, division label, inta spin ){
    inta rr = CanonicalRank(f1, name(f1,label), spin );
    division ll = f1.name[name(f1,label)].chainNext;
    while ( ll != nullName ){
        rr += CanonicalRank(f1, name(f1,ll), spin);
        ll = f1.name[name(f1,ll)].chainNext;
    }
    return rr;
}

inta Rank( sinc_label f1 ,   division label ){
    inta sp,ra=0;
    for ( sp = 0 ; sp < spins(f1, label);sp++)
        ra += CanonicalRank(f1, label, sp);
    return ra;
}

inta spins ( sinc_label f1 ,   division label ){
      spinType sp = f1.name[label].spinor;

#ifndef APPLE
#ifdef OMP
    if (sp == parallel)
            return  f1.rt->NLanes;
#endif
#endif
    
    
    if ( sp == real )
        return 1;
    else if ( sp == cmpl )
        return 2;
    else if (sp == parallel )
        return MAX_CORE;
    return 0;
}

void assignOneWithPointers(   sinc_label f1,   division oneMat , inta label){
    inta space;
      blockType tv ;
    f1.name[oneMat].species = matrix;
    {
        for ( space = 0; space < SPACE ; space++)
            if ((f1.canon[space].body != nada ) &&( f1.canon[space].label == label || label == 0)){
                f1.name[oneMat].space[space].body = one;
            }else {
                f1.name[oneMat].space[space].body = nada;
            }
    }

    for ( tv = tv1 ; tv <= tv4 ; tv++){
        f1.name[oneMat+tv].Current[0] = part(f1, oneMat);
        f1.name[oneMat+tv].spinor = spins(f1, oneMat);
        f1.name[oneMat+tv].species = matrix;
        f1.name[oneMat+tv].name = oneMat;
        for ( space = 0; space < SPACE ; space++)
            if ((f1.canon[space].body != nada ) &&( f1.canon[space].label == label || label == 0)){
                f1.name[oneMat+tv].space[space].block = tv;
            }else {
                f1.name[oneMat+tv].space[space].block = id0;
            }
    }
    
    
    
    return;
}

void assignParticle(  sinc_label f1,   division ma, inta label ,   bodyType ba ){
    inta space ;
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.canon[space].body != nada )
        if ( f1.canon[space].label == label || label == 0 )
            f1.name[ma].space[space].body = ba;
}

void assignTwoWithPointers(   sinc_label f1,   division twoMat ){
    inta space;
      blockType e ;
      blockType maxee = e12;
    if ( bodies(f1, twoMat) == three )
        maxee = e23;
    if ( bodies(f1, twoMat) == four )
        maxee = e34;
    for ( e = e12 ; e <= maxee ; e++){
        f1.name[twoMat+e].name = twoMat;
        f1.name[twoMat+e].species = matrix;
        f1.name[twoMat+e].spinor = spins(f1,twoMat);

        for ( space = 0; space < SPACE ; space++){
            f1.name[twoMat+e].space[space].block = e;
        }
    }
    return;
}

inta sizeofDivision(  sinc_label f1,   division head, inta space ){
    if ( space < SPACE && f1.canon[space].body == nada)
        return 0;
    if ( f1.name[head].Partition == 0 )
        return 0;
    if ( f1.name[head].memory == bufferAllocation && space < SPACE )
        return 0;
    if ( f1.name[head].memory == objectAllocation && space == SPACE )
        return 0;
    if ( f1.name[head].memory == noAllocation)
        return 0;
    if ( head == f1.end )
        return 0;
    return alloc( f1, head, space)*part(f1,head)*spins(f1,head);
}


void  fromBeginning(   sinc_label  f1 ,  division new,   division head ){
    if ( new > f1.end ){
        printf("fromBeginning, name outside bounds\n");
        fflush(stdout);
        exit(0);
    }
    inta space;
    
    if ( head == nullName ){
        
        for ( space = 0; space <= SPACE ; space++)
            f1.name[new].space[space].Address = 0;
        
    }else {
        for ( space = 0; space <= SPACE ; space++){
            f1.name[new].space[space].Address = f1.name[head].space[space].Address;
            if ( f1.name[head].space[space].Address == -1 )
            {
                printf("fromBeginning, something is wrong in Model.c!\n");
            }
        }
    }
    if ( f1.name[head].memory == objectAllocation){
        for ( space = 0; space < SPACE ; space++)
      //  if( f1.canon[space].body != nada)
        {
            f1.name[new].space[space].Address += sizeofDivision(f1,head,space);
        }
    } else if ( f1.name[head].memory == bufferAllocation){
        for ( space = SPACE; space <= SPACE ; space++)
        {
            f1.name[new].space[space].Address += sizeofDivision(f1,head,space);
        }
    } else {
        //nothing
    }
#if VERBOSE
    printf("||%d::", new);
    for ( space = 0; space <= SPACE ; space++)
       printf("%lld:", f1.name[new].space[space].Address);
    printf("\n\n");
#endif
    return;
}

floata*  myStreams (   sinc_label f1,   division label ,inta spin){
    inta space = SPACE;//buffer-space
    if ( spin < 0 || spin >= spins(f1, label)){
        printf("myStreams, has a spinor problem\n");
        exit(5);
    }

    inta leng = alloc(f1, name(f1,label),space);
    inta partit = part(f1, name(f1,label));
    if (f1.name[name(f1,label)].space[SPACE].Address == -1  ){
        
        exit(0);
    }

    if ( f1.name[label].memory == bufferAllocation){
        if ( name(f1,label) != label ){
            return f1.canon[space].stream+f1.name[name(f1,label)].space[space].Address + leng * partit * spin + leng*f1.name[label].Begin[spin];//HERE->BEGIN
        }
        else{
            return f1.canon[space].stream+f1.name[name(f1,label)].space[space].Address + leng * partit * spin  ;
        }
    }else {
        printf("not a buffer!\n %d",label);
        exit(0);
    }
    return NULL;
}

floata* pStreams (   sinc_label *f1,   division label ,inta spin, inta space ){
    return streams(*f1, label, spin, space);
}
floata* streams (   sinc_label f1,   division label ,inta spin, inta space ){
    floata * uu ;
    if ( spin < 0 || spin >= spins(f1, label)){
        printf("spins %d %d\n",label, spin);
        exit(5);
    }
    if ( space < 0 || space >= SPACE+1 ){
        printf("space\n");

        exit(4);
        
    }
    if ( label > f1.end ){
        printf("past end\n,%d",label);
        fflush(stdout);
        exit(0);
    }
    
    
    if (f1.name[name(f1,label)].space[space].Address == -1  ){
        printf("bad reference\t calling a division that is not allocated");
        printf("%d %d\n", label, name(f1,label));
        exit(0);
    }
    if ( f1.name[label].memory == objectAllocation){
        
        inta leng = alloc(f1, name(f1,label),space);
        inta partit = part(f1, name(f1,label));
        
        if ( name(f1,label) != label ){
            return f1.canon[space].stream+f1.name[name(f1,label)].space[space].Address + leng * partit * spin + leng*f1.name[label].Begin[spin] ;
        }
        else{
             uu =  f1.canon[space].stream+f1.name[name(f1,label)].space[space].Address + leng * partit * spin  ;
            return uu;
        }
    } else if ( f1.name[label].memory == bufferAllocation){
        myStreams(f1, label, spin);
    }
    
    return NULL;
}

void xsAdd ( double scalar , inta dim ,  sinc_label f1 ,   division targ ,inta tspin,  sinc_label  f2 ,   division orig,inta o,inta ospin ){
    inta M2 = alloc(f1, orig, dim);
    inta N2 = alloc(f1, targ, dim);
    inta flag = (N2 == M2),space=dim;

    if ( flag && f2.name[orig].memory == objectAllocation && f1.name[targ].memory == objectAllocation){
        cblas_dcopy(N2, streams(f2,orig,ospin,space)+N2*o,1,streams(f1,targ,tspin,space)+CanonicalRank(f1, targ, tspin)*N2,1);
        if ( scalar != 1. )
            cblas_dscal(N2, scalar, streams(f1,targ,tspin,space)+CanonicalRank(f1, targ, tspin)*N2,1);
    }
    else {
        printf("add\n");

        exit(0);
    }
}

void xsEqu ( double scalar , inta dim ,  sinc_label f1 ,   division targ ,inta t,inta tspin,inta dim2,  sinc_label  f2 ,   division orig,inta o,inta ospin ){
    inta M2 = alloc(f1, orig, dim);
    inta N2 = alloc(f1, targ, dim2);
    inta flag = (N2 == M2);
    
    if ( flag && f2.name[orig].memory == objectAllocation && f1.name[targ].memory == objectAllocation){
        cblas_dcopy(N2, streams(f2,orig,ospin,dim)+N2*o,1,streams(f1,targ,tspin,dim2)+t*N2,1);
        if ( scalar != 1. )
            cblas_dscal(N2, scalar, streams(f1,targ,tspin,dim2)+t*N2,1);
    }
    else {
        printf("add\n");
        
        exit(0);
    }
}


double xEqua ( sinc_label f1 , division targ ,inta tspin, sinc_label  f2 , division orig,inta ospin ){
    inta space;
    inta eb = CanonicalRank(f2,orig,ospin);
    inta M2[SPACE];
    length(f2, orig, M2);
    inta M1[SPACE];
    length(f1,targ,M1);
    
    
  //  if ( f2.name[orig].memory == objectAllocation && f1.name[targ].memory == objectAllocation)
    {
        
      //  if ( name(f1,targ) == targ  )
        {
            
            if ( (part(f1,targ) < eb) )
            {
                printf("tEqual.. memory %d %d %d %d\n", targ,part(f1,targ), orig,part(f1,orig) );
                printf("partition %d\n", eb);
                exit(0);
            }
            zero(f1,targ,tspin);

            for ( space = 0 ; space < SPACE ; space++)
                if ( f1.canon[space].body != nada){
                    if (M1[space] != M2[space]){
                        if ( Bodies(f2, orig,space) == one )
                            xOneBand(f2,space, orig,  ospin, f1, targ,tspin,0 );
                        else if ( Bodies(f2, orig,space) == two )
                            xTwoBand(f2, space,orig,  ospin, f1, targ,tspin,0 );
                        else if ( Bodies(f2, orig,space) == three )
                            xThreeBand(f2,space, orig,  ospin, f1, targ,tspin,0 );
                        else if (Bodies(f2, orig,space) == four )
                            xFourBand(f2,space, orig,  ospin, f1, targ,tspin,0 );
                        }
                    else {
                        {
                            cblas_dcopy(eb*M2[space], streams(f2,orig,ospin,space),1,streams(f1,targ,tspin,space),1);
                            f1.name[targ].Current[tspin] = f2.name[orig].Current[ospin];
                        }
                    }
                }
            //f1.name[name(f1,targ)].header = header(f2, name(f2,orig));
        }
    }
    return 0;
}

double traceOne(   sinc_label  f1 ,   division label , inta spin ){
    floata * base;
    double sum,sum2,product;
    inta l,i,space;
    inta N2[SPACE],N1[SPACE];
    length(f1, label,N2);
    for ( i = 0; i < SPACE ; i++)
        if ( f1.canon[i].body != nada)

        N1[i] = sqrt(N2[i]);

    
    if ( species(f1, name(f1,label)) != vector && bodies(f1,label) != two)
    if ( species(f1, name(f1,label)) != matrix ){
        printf("\nvector->matrix\n");
        return sqrt(pMatrixElement(f1, label ,spin,nullOverlap,0,label ,spin));
    }

    sum2 = 0.;
    for ( l = 0 ; l < CanonicalRank(f1,label,spin); l++)
    {
        product = 1.;
        for ( space = 0; space < SPACE ; space++)
            if ( f1.canon[space].body != nada )
            if ( N2[space] ){
                sum = 0.;
                base = streams(f1, label, spin, space )+l*N2[space];
                
                for ( i = 0; i < N1[space] ; i++ ){
                    sum += base[ i*N1[space]+i];
                }
                product *= sum;
                
            }
        sum2 += product;
    }
    
    
    return sum2;
}


inta ready (   sinc_label f1){
    inta readyMemory = 1;
    inta readyVector = 1;
    inta space;
    if ( ! f1.bootedMemory || f1.name == NULL )
        readyMemory = 0;
    
    if ( readyMemory )
        for ( space = 0 ; space <= SPACE ; space++)
            if ( f1.canon[space].stream == NULL )
            readyMemory = 0;
    
    
    if ( readyMemory )
        if ( CanonicalRank(f1, eigenVectors , 0 ) == 0 ){
            printf("passing over stage because vector is null\n");
            readyVector = 0;
        }
    
    return readyVector && readyMemory;
}

inta bootedQ (   sinc_label f1){
    inta readyMemory = 1;
    inta space;
    if ( ! f1.bootedMemory || f1.name == NULL )
        readyMemory = 0;
    
    if ( readyMemory )
        for ( space = 0 ; space <= SPACE ; space++)
            if ( f1.canon[space].stream == NULL )
                readyMemory = 0;
    
    
    return readyMemory;
}

inta balance (  sinc_label  f1,    division alloy, inta spin){
    inta l;
    assignCores(f1, 1);
    inta L1 = CanonicalRank(f1, alloy,spin);
    long double prod,norm[SPACE],mn,mx,snorm ;
    inta M2[SPACE],space;
    length(f1,alloy,M2);
    inta iOne = 1;
    {
        #ifdef OMP
        #pragma omp parallel for private (mn,mx,l,space,prod,snorm,norm) schedule(dynamic,1)
        #endif
        for ( l = 0; l < L1 ;l++){
            
            for ( space = 0; space < SPACE ; space++)
                if ( f1.canon[space].body != nada){
                    norm[space] = cblas_dnrm2(M2[space], streams(f1, alloy,spin,space)+l*M2[space],iOne);
                }
            mn = norm[1];
            mx = norm[1];
            for ( space = 1 ;space < SPACE ; space++)
                if ( f1.canon[space].body != nada){
                    if ( mn > norm[space])
                        mn = norm[space];
                    if ( mx < norm[space])
                        mx = norm[space];
                }

            
            if ( mx/mn > 10 ){
                prod = 1;
                for ( space = 0 ;space < SPACE ; space++)
                    if ( f1.canon[space].body != nada){
                        prod *= norm[space];
                    }
                for ( space = 0; space < SPACE ; space++)
                   if ( f1.canon[space].body != nada)
                   {
                       snorm = 1./norm[space] ;
                       if (! space )
                           snorm *= prod;
                       cblas_dscal(M2[space], snorm, streams(f1, alloy,spin,space)+l*M2[space],iOne);
                   }
            }
        }
    }
    return 0;
}
double tEqua (   sinc_label f1 ,   division targ ,inta tspin,   division orig,inta ospin ){
    return xEqua(f1, targ,tspin, f1, orig,ospin);
}

inta tEquals(   sinc_label f1 ,   division left ,   division right){
      spinType spl,spr;
    spl = f1.name[left].spinor;
    spr = f1.name[right].spinor;

    if ( spl > spr ){
        tClear(f1, left);
        tEqua(f1, left, 0, right , 0);
        return 1;
    }else if ( spl == spr ){
        tEqua(f1, left, 0, right , 0);
        if ( spl == cmpl )
            tEqua(f1, left, 1, right , 1);
        return 0;
    }else {
        tEqua(f1, left, 0, right , 0);
 
        
        printf("cmpl!\n");
      //  exit(0);
    }
    return 0;
    
}



inta tAddTwo(   sinc_label f1 ,   division left ,   division right){
      spinType spl,spr;
    spl = f1.name[left].spinor;
    spr = f1.name[right].spinor;
    
    if ( spl > spr ){
        tClear(f1, left);
        tAddTw(f1, left, 0, right , 0);
        return 1;
    }else if ( spl == spr ){
        tAddTw(f1, left, 0, right , 0);
        if ( spl == cmpl )
            tAddTw(f1, left, 1, right , 1);
        return 0;
    }else {
        tAddTw(f1, left, 0, right , 0);

        printf("cmpl!\n");
       // exit(0);
    }
    return 0;
}


  division defSpiralVector(   sinc_label *f1, inta spiralOp,   division ket){
      division opi = spiralOp;
    inta term = 0;
    double cweight,weight = 0.;
    while ( opi != nullName )
    {
        weight += 1.;
        opi = f1->name[opi].linkNext;
        term++;
    }
    
    
    opi = spiralOp;
      division buf,prev=0,spiral=0;
    inta t, curr = 0,sp,len,al[2];
    al[0] = f1->name[ket].Current[0];
    al[1] = f1->name[ket].Current[1];

    for ( t = 0 ; t < term ; t++){
            buf = anotherLabel(f1, all, nada);
            if ( ! prev ){
                spiral = buf;
            }else {
                f1->name[prev].linkNext = buf;
            }
        cweight = 1.;
        len = floor((part(*f1,ket)-term)*cweight/weight)+1;

        
        
        {
            for ( sp = 0  ; sp < f1->cmpl;sp++){
                f1->name[buf].name = ket;
                f1->name[buf].Partition = len;//not actually allocated,,,,it should never have its own name.
                f1->name[buf].Current[sp] = imin(len,al[sp]);
                f1->name[buf].Begin[sp]   = curr;
                al[sp] = imax(0,al[sp]-len);
            }
            curr += len;
        }
        prev = buf;
        opi = f1->name[opi].linkNext;
    }
    return spiral;
}


  division defRefVector(   sinc_label *f1, inta spiralOp,   division ket){
      division opi = spiralOp;
    inta term = 0;
    double cweight,weight = 0.;
    while ( opi != nullName )
    {
        weight += 1.;
        opi = f1->name[opi].linkNext;
        term++;
    }
        
    
    opi = spiralOp;
      division buf,prev=0,spiral=0;
    inta t, sp,len,al[2];
    al[0] = f1->name[ket].Current[0];
    al[1] = f1->name[ket].Current[1];

    for ( t = 0 ; t < term ; t++){
            buf = anotherLabel(f1, all, nada);
            if ( ! prev ){
                spiral = buf;
            }else {
                f1->name[prev].linkNext = buf;
            }
        cweight = 1.;
        len = floor((part(*f1,ket)-term)*cweight/weight)+1;

        //ALL POINTING AT SAME STRUCTURE
        {
            for ( sp = 0  ; sp < f1->cmpl;sp++){
                f1->name[buf].name = ket;
                f1->name[buf].Partition = len;//not actually allocated,,,,it should never have its own name.
                f1->name[buf].Current[sp] = al[sp];
                f1->name[buf].Begin[sp]   = 0;
                //al[sp] = imax(0,al[sp]-len);
            }
           // curr += len;
        }
        prev = buf;
        opi = f1->name[opi].linkNext;
    }
    return spiral;
}


division defSpiralMatrix( sinc_label *f1, division H){
    
    division pt = H,buf,prev=0,spiral=0;
    inta term = 0;
    do{
        if ( f1->name[pt].species == matrix ){
            buf = anotherLabel(f1, all, nada);
            if ( prev == 0 ){
                spiral = buf;
            }else {
                f1->name[prev].linkNext = buf;
            }
            f1->name[buf].name = pt;
            f1->name[buf].chainNext = pt;
            f1->name[buf].species = matrix;
            term++;
            prev = buf;
        }
        pt = f1->name[pt].linkNext;
    }while ( pt != nullName);
    buf = anotherLabel(f1, all, nada);
    f1->name[buf].species = scalar;
    f1->name[buf].name = nullName;

    return spiral;
}


division defSpiralGrid(   sinc_label *f1,   division bra, inta term, double diagonalPreference){
      division buf,prev=0,spiral=0;
    inta  t,tt;
      spinType sp;
    inta diagonal,offDiagonal,curr;
    
    
    for ( t = 0 ; t < term; t++){
        curr = 0;
        for ( tt = 0; tt < term ; tt++){
            buf=  anotherLabel(f1, all, nada);
            if ( prev == 0 ){
                spiral = buf;
            }else {
                f1->name[prev].chainNext = buf;
            }
            
            //Partition = diagonal + offDiagonal
            // diagonal = diagonalPreference * offDiagonal
//            if (0){
//                if ( term > 1 )
//                {
//                    diagonal = ceil((part(*f1, bra+t)-term)*diagonalPreference/(1+diagonalPreference))+1;
//                    offDiagonal = (part(*f1, bra+t)-diagonal)/(term-1);
//                }
//                else{
//                    diagonal = (part(*f1, bra+t)-term)+1;
//                    offDiagonal = 0;
//                }
//
//            }
          //  else
            {
                offDiagonal = (part(*f1, bra+t))/term;
                diagonal = part(*f1, bra+t) - (term-1) * offDiagonal;
                
                
            }
            for ( sp = 0  ; sp < f1->name[bra].spinor;sp++){
                f1->name[buf].name = bra+t;
                if ( t == tt ){
                    f1->name[buf].Current[sp] = 0;
                    f1->name[buf].Partition = diagonal;
                    f1->name[buf].Begin[sp] = curr;
                    curr += diagonal;
                }
                else{
                    f1->name[buf].Current[sp] = 0;
                    f1->name[buf].Partition = offDiagonal;
                    f1->name[buf].Begin[sp] = curr;
                    curr += offDiagonal;

                }
            }
            prev = buf;

        }
    }
    return spiral;
}


inta zeroSpiraly(   sinc_label f1,   division spiral){
    inta i,space,spin,flag = 0;
    floata * point;
      division spir = spiral;
    while ( spir != nullName ){
        for (spin = 0 ; spin < spins(f1,spiral);spin++){
            flag = 0;

            for (space = 0 ; space < SPACE ; space++)
                if ( f1.canon[space].body != nada){
                    point = streams(f1, spir, spin, space);
                    for ( i = f1.name[spir].Current[spin] ; i < f1.name[spir].Partition; i++)
                    {
                        flag = 1;
                        point[i] = 0.;
                    }
                }
            if ( flag )
                printf("zero %d->%d %d:%d-%d\n", spir,name(f1,spir), f1.name[spir].Begin[spin],f1.name[spir].Current[spin], f1.name[spir].Partition);
        }
        spir = f1.name[spir].chainNext;
    }
    return 0;
}


inta tScaleOne(   sinc_label f1,   division label,inta spin, double scalar ){
    
    if ( scalar == 1. )
        return 0;
    
    
    inta L1 = CanonicalRank(f1,label,spin);
    
    double prod;
    
    
    inta M2[SPACE];
    length(f1,label,M2);
    
    prod = 1.;
    
#if 0
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.canon[space].body != nada ){
            dimCount++;
        }
            
    
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.canon[space].body != nada ){
            cblas_dscal(L1*M2[space],pow(scale,1./dimCount), streams(f1,label,spin,space),1);
        }
    
    if ( scalar < 0 )
        cblas_dscal(L1*M2[0],-1., streams(f1,label,spin,0),1);

#else
    
        cblas_dscal(L1*M2[0],scalar, streams(f1,label,spin,0),1);

    
    
    
    
    
#endif
    
    
    
    
//
//    if ( SPACE == 3 ){
//        scale = pow ( scale , 1./3.);
//        if ( scalar < 0 )
//            scale *= -1;
//        cblas_dscal(L1*M2[0],scale, streams(f1,label,spin,0),1);
//        cblas_dscal(L1*M2[1],scale, streams(f1,label,spin,1),1);
//        cblas_dscal(L1*M2[2],scale, streams(f1,label,spin,2),1);
//    }else if ( SPACE == 2 ){
//        scale = pow ( scale , 1./2.);
//        cblas_dscal(L1*M2[1],scale, streams(f1,label,spin,1),1);
//        if ( scalar < 0 )
//            scale *= -1;
//        cblas_dscal(L1*M2[0],scale, streams(f1,label,spin,0),1);
//    }else if ( SPACE == 1 ){
//        cblas_dscal(L1*M2[0],scalar, streams(f1,label,spin,0),1);
//    }

    
    return 0;
}


inta tScale(   sinc_label f1,   division label, DCOMPLEX scalar ){
    if ( spins(f1, label ) == real ){
        
        if ( cimag(scalar) != 0. ){
            printf("errr in Scalar \n");
            printf("%f \n", cimag(scalar));
            exit(1);
        } else{
            tScaleOne(f1, label, 0, creal(scalar));
        }
    }else{
//        if ( fabs(cimag(scalar)) > f1.rt->TARGET && fabs(creal(scalar)) > f1.rt->TARGET){
//            tClear(f1, scalarTemp);
//            tAddTw(f1, scalarTemp, 0, label, 1);
//            tScaleOne(f1, scalarTemp, 0,  -cimag(scalar)/creal(scalar));
//            tAddTw(f1, scalarTemp, 0, label, 0);
//            tScaleOne(f1, scalarTemp, 0,  creal(scalar));
//            tCycleDecompostionListOneMP(-1, f1, scalarTemp, 0, NULL, label, 0, f1.rt->vCANON, part(f1,label), 1);
//
//            tClear(f1, scalarTemp);
//            tAddTw(f1, scalarTemp, 0, label, 1);
//            tScaleOne(f1, scalarTemp, 0, creal(scalar)/cimag(scalar));
//            tAddTw(f1, scalarTemp, 0, label, 0);
//            tScaleOne(f1, scalarTemp, 0,  cimag(scalar));
//            tCycleDecompostionListOneMP(-1, f1, scalarTemp, 0, NULL, label, 1, f1.rt->vCANON, part(f1,label), 1);
//
//        }else
        
        if ( fabs(cimag(scalar)) < f1.rt->THRESHOLD && fabs(creal(scalar)) > f1.rt->THRESHOLD)
        {
            tScaleOne(f1, label, 0, creal(scalar));
            tScaleOne(f1, label, 1, creal(scalar));
            
        }else if ( fabs(cimag(scalar)) > f1.rt->THRESHOLD && fabs(creal(scalar)) < f1.rt->THRESHOLD)
        {
            tScaleOne(f1, label, 0, cimag(scalar));
            tScaleOne(f1, label, 1, -cimag(scalar));
            
            tEqua(f1, scalarTemp, 0, label, 0);
            tEqua(f1, label, 0, label, 1);
            tEqua(f1, label, 1, scalarTemp, 0);
            
        }
    }
    
    return 0;
}

inta tAddTw(   sinc_label  f1 ,   division left, inta lspin,    division right , inta rspin){
    return xAddTw(f1,left, lspin, f1, right, rspin);
}

inta xAddTw(   sinc_label f1 ,   division left, inta lspin,  sinc_label f2 ,    division right , inta rspin){
    if ( CanonicalRank(f2,right,rspin) ){
        inta LL = f1.name[left].Current[lspin];
        inta LR = f2.name[right].Current[rspin];
        //printf("++ %lld %lld\n", LL,LR);
        inta MM = LL+LR;
        if ( MM > f1.name[left].Partition ){
            //tGetType allocation!
            printf("xAddTw, copy needs more allocation %d -> %d\n",right,left);
            exit(8);
        }
        inta M2[SPACE],space;
        length(f1, left, M2);
        if ( f1.name[right].memory == objectAllocation ){
            for ( space = 0; space < SPACE; space++)
                if (f1.canon[space].body != nada){
//                if ( (species(f1, right ) == vector) || (f1.name[right].space[space].body != nada && species(f1, right ) == matrix)){
                    cblas_dcopy(LR*M2[space], streams(f2,right,rspin,space),1,streams(f1,left,lspin,space)+LL*M2[space],1);
//                }
                
//                else if ((f1.name[right].space[space].body == nada && species(f1, right ) == matrix && f1.name[right].space[space].body == id0) ){
//
//                    inta i,l,n1 = outerVectorLen(f1,bodies(f1,name(f1,right)),space);//question!
//                    for ( i= 0 ; i < LR*M2[space] ; i++)
//                        (streams(f1,left,lspin,space)+LL*M2[space])[i] = 0;
//                    for ( l = 0; l < LR ; l++)
//                        for ( i = 0; i < n1 ; i++)
//                            (streams(f1,left,lspin,space)+(LL+l)*M2[space])[i*n1+i] = 1;
//                }
            
            
            
                }
        }else if ( f1.name[right].memory == bufferAllocation )
        {
            printf("xAddTw, buffers dont count\n");
            exit(8);
        }
            else if ( f1.name[right].memory == noAllocation ){
            printf("xAddTw, not allocated\n");
            exit(0);
        }
        f1.name[left].Current[lspin] += LR;
    }
    return 0;
}



//inta tPauli (   field * f1  ){
//    tClear(f1, PauliZ);
//    tClear(f1, PauliX);
//    if ( spins(f1,PauliZ) >= 2){
//        tId(f1, PauliZ,f1.arraySpin[f1.name[PauliZ].spinor][0][0]);
//        tScale(f1, PauliZ,-1);
//        tId(f1, PauliZ,f1.arraySpin[f1.name[PauliZ].spinor][1][1]);
//    }
//    if (spins(f1,PauliZ)  >= 3 ){
//        tId(f1, PauliX,f1.arraySpin[f1.name[PauliX].spinor][0][1] );
//    }
//    return 0;
//}

int double_cmp( const void * a ,const void * b ){
    const double* A = (const double* )a;
    const double* B = (const double* )b;
    if ( *A < * B )
        return 1;
    else if ( *A == *B )
        return 0;
    else return -1;
}


double levelDetermine ( inta M1 , double * array , double level){
    double *temp = array+M1*M1*M1,sum,psum;
    inta i;
    cblas_dcopy(M1*M1*M1, array, 1, temp, 1);
    qsort(temp, M1*M1*M1, sizeof(double), &double_cmp);
    sum = 0.;
    for ( i = 0 ;i < M1*M1*M1 ; i++)
        sum += temp[i];
    i= 0;
    psum = 0.;
    while (psum < level*sum ){
        psum += temp[i++];
    }
    
    return temp[i-1];
}


inta tId (   sinc_label f1 ,   division label,inta spin ){
    
    inta I1,I2,space;
    inta Current ;
    {
        
//        if ( f1.name[label].Current[spin] >= f1.name[label].Partition ){
//            printf("%d %d\n", label, spin);
//            printf("tryed to add to full array\n");
//            return 0;
//        }
        Current =  f1.name[label].Current[spin]++;
    }
    
    {
        if ( f1.name[label].species == vector || f1.name[label].species == outerVector){
            
            inta B1[SPACE];
            length(f1, label, B1);
            for ( space = 0; space < SPACE ; space++)
                if ( f1.canon[space].body != nada)
                {
                
                floata  * stream = streams(f1,label,spin,space)+Current*B1[space];
                for ( I2 = 0 ; I2 < B1[space] ; I2++){
                    stream[I2] = sign(I2);
                }
                    //stream[(B1[space]-1)/2]=0.;
            }
        }
        
        else if  ( f1.name[label].species == matrix ) {
            inta B1;
            

            
            for ( space = 0; space < SPACE ; space++)
                if ( f1.canon[space].body != nada)
                {
                B1 = outerVectorLen(f1,bodies(f1, name(f1,label)), space);
                floata * stream = streams(f1,label,spin,space)+Current*B1*B1;
                for ( I1 = 0 ; I1 < B1 ; I1++)
                    for ( I2 = 0 ; I2 < B1 ; I2++)
                        stream[I1*B1+I2] =0.;
                for ( I1 = 0 ; I1 < B1 ; I1++)
                    {
                        
                        stream[I1*B1+I1] = 1;
                        
                    }

                }
        }
        
        
    }
    return 0;
}

//inta tOv (   sinc_label f1 ,   division label,inta spin ){
//
//    inta I1,I2,space;
//    inta Current ;
//    {
//
//        if ( f1.name[label].Current[spin] >= f1.name[label].Partition ){
//            printf("%d %d\n", label, spin);
//            printf("tryed to add to full array\n");
//            return 0;
//        }
//        Current =  f1.name[label].Current[spin]++;
//    }
//
//    {
//        if ( f1.name[label].species == vector || f1.name[label].species == outerVector){
//
//            inta B1[SPACE];
//            length(f1, label, B1);
//            for ( space = 0; space < SPACE ; space++)
//                if ( f1.canon[space].body != nada)
//                {
//
//                    floata  * stream = streams(f1,label,spin,space)+Current*B1[space];
//                    for ( I2 = 0 ; I2 < B1[space] ; I2++){
//                        stream[I2] = sign(I2);
//                    }
//                }
//        }
//
//        else if  ( f1.name[label].species == matrix && bodies(f1, label) == one) {
//            inta B1;
//
//
//
//            for ( space = 0; space < SPACE ; space++)
//                if ( f1.canon[space].body != nada)
//                {
//                    B1 = outerVectorLen(f1,bodies(f1, name(f1,label)), space);
//                    floata * stream = streams(f1,label,spin,space)+Current*B1*B1;
//                    for ( I1 = 0 ; I1 < B1 ; I1++)
//                        for ( I2 = 0 ; I2 < B1 ; I2++)
//                            stream[I1*B1+I2] =creal(BoB(f1.canon[space].basisList[I1], f1.canon[space].basisList[I2]));
//
//                }
//        }
//
//
//    }
//    return 0;
//}


inta tReplace(   sinc_label f1 ,   division label,inta spin,inta space,inta l ){
    
    inta I1,I2;
    {
        
        if ( f1.name[label].Current[spin] < l )
            return 0;
       }
    {
        if ( f1.name[name(f1,label)].species == vector ){
            inta B1 = vectorLen(f1, space);
            {
                floata  * stream = streams(f1,label,spin,space)+l*B1;
                for ( I2 = 0 ; I2 < B1 ; I2++){
                    stream[I2] = sign(I2%2);
                }
            }
        }
        
        //        else if ( ( f1.name[label].species == matrix && bodies(f1,label)== two) ||  (f1.name[label].species == quartic && bodies(f1,label) == one)){
        //            for ( space = 0; space < SPACE ; space++){
        //                inta * B1;
        //                B1 = vectorLen(f1,label);
        //
        //
        //                floata * stream = streams(f1,label,spin,space)+Current*B1[space]*B1[space]*B1[space]*B1[space];
        //                for ( I1 = 0 ; I1 < B1[space]*B1[space] ; I1++)
        //                    for ( I2 = 0 ; I2 < B1[space]*B1[space] ; I2++)
        //                        stream[I1*B1[space]*B1[space]+I2] =0.;
        //
        //                for ( I1 = 0 ; I1 < B1[space]*B1[space] ; I1++)
        //                    for ( I2 = 0 ; I2 < B1[space]*B1[space] ; I2++)
        //                    {
        //
        //                        if ( I1==I2  )
        //                            stream[I1*B1[space]*B1[space]+I2] = 1;
        //
        //                    }
        //            }
        //        }
        else if  ( f1.name[name(f1,label)].species == matrix ) {
            inta B1;
            B1 = vectorLen(f1, space);
            {
                floata * stream = streams(f1,label,spin,space)+l*B1*B1;
                for ( I1 = 0 ; I1 < B1 ; I1++)
                    for ( I2 = 0 ; I2 < B1 ; I2++)
                        stream[I1*B1+I2] =0.;
                for ( I1 = 0 ; I1 < B1 ; I1++)
                    for ( I2 = 0 ; I2 < B1 ; I2++)
                    {
                        if ( I1==I2 )
                            stream[I1*B1+I2] = 1;
                    }
            }
        }else
        {
            printf("ask for upgrade with outer-vectors\n");
            exit(1);
        }
        
        
    }
    return 0;
}

inta tBoot (   sinc_label f1 ,   division label,inta spin,floata scale ){
    
    inta I1,I2,I3,I4,space;
    inta Current ;
    inta B1[SPACE];

    {
        
        if ( f1.name[label].Current[spin] >= f1.name[label].Partition )
            return 0;
        Current =  f1.name[label].Current[spin]++;
    }
    length1(f1, B1);

    for ( space = 0; space < SPACE ; space++)
        if ( f1.canon[space].body != nada)

        {
            if ( f1.canon[space].body == one){
            {
                            floata  * stream = streams(f1,label,spin,space)+Current*B1[space];
                                for ( I2 = 0 ; I2 < B1[space] ; I2++){
                                    stream[I2] = exp(-(I2-(B1[space]-1)/2)*(I2-(B1[space]-1)/2)*scale);
                                }
                        }
                    }else
            if (f1.canon[space].body == two){
                    {
                
                floata  * stream = streams(f1,label,spin,space)+Current*B1[space]*B1[space];
                for ( I1 = 0 ; I1< B1[space] ; I1++)
                    for ( I2 = 0 ; I2 < B1[space] ; I2++){
                        stream[I1*B1[space]+I2] = exp(-(I1-(B1[space]-1)/2)*(I1-(B1[space]-1)/2)*scale)*exp(-(I2-(B1[space]-1)/2)*(I2-(B1[space]-1)/2)*scale);
//              if ( I1 < I2 )
//                  stream[I1*B1[space]+I2]  *= -1;
//                if ( I1 == I2 )
//                    stream[I1*B1[space]+I2]  *= 0;
//
                    }
                
            }
        }
        else
            if ( f1.canon[space].body == three){
                    {
                
                floata  * stream = streams(f1,label,spin,space)+Current*B1[space]*B1[space]*B1[space];
                for ( I1 = 0 ; I1< B1[space] ; I1++)
                    for ( I2 = 0 ; I2 < B1[space] ; I2++)
                        for ( I3 = 0 ; I3 < B1[space] ; I3++){
                            stream[I3*B1[space]*B1[space]+I1*B1[space]+I2] = exp(-(I1-(B1[space]-1)/2)*(I1-(B1[space]-1)/2)*scale)*exp(-(I2-(B1[space]-1)/2)*(I2-(B1[space]-1)/2)*scale)*exp(-(I3-(B1[space]-1)/2)*(I3-(B1[space]-1)/2)*scale);
                }
            }
        }
            else
                if ( f1.canon[space].body == four){
                        {
                    
                    floata  * stream = streams(f1,label,spin,space)+Current*B1[space]*B1[space]*B1[space]*B1[space];
                    for ( I1 = 0 ; I1< B1[space] ; I1++)
                        for ( I2 = 0 ; I2 < B1[space] ; I2++)
                            for ( I3 = 0 ; I3 < B1[space] ; I3++)
                                for ( I4 = 0 ; I4 < B1[space] ; I4++)
                            {
                                stream[I4*B1[space]*B1[space]*B1[space]+I3*B1[space]*B1[space]+I1*B1[space]+I2] = exp(-(I1-(B1[space]-1)/2)*(I1-(B1[space]-1)/2)*scale)*exp(-(I2-(B1[space]-1)/2)*(I2-(B1[space]-1)/2)*scale)*exp(-(I3-(B1[space]-1)/2)*(I3-(B1[space]-1)/2)*scale)*exp(-(I4-(B1[space]-1)/2)*(I4-(B1[space]-1)/2));
                    }
                }
            }



      
        
    }
    return 0;
}

void SG( sinc_label f1, division vector ,inta spin,floata amplitude,  inta *gamma ){
    inta space,vc,vsp,msp,mss=0,vn1,v,m,n,Current,v0;
    floata variable ;
    bodyType body ;
    if ( f1.name[vector].Current[spin] >= f1.name[vector].Partition )
        return;
    Current =  f1.name[vector].Current[spin]++;
    msp = 0;
    for  ( space =0; space < SPACE ; space++){
        if ( f1.canon[space].body != nada){
            v0 = vectorLen(f1, space);
            for ( vc = 0; vc < v0 ; vc++){
                vsp = 1;
                variable = 1.0;
                vn1 = vector1Len(f1,space);
                mss = msp ;
                for ( body = one ; body <= f1.canon[space].body ; body++){
                    v = (vc/vsp)%vn1-(vn1-1)/2;
                    vsp *= vn1;
                    m = (gamma[mss]);
                    n = (gamma[mss+1]);
                    mss += 2;
                    ///n = 1 --> 1/2 which is the lowest level...
                    ///n = 3 -> 3/2
                    ///even n's are NOT to be used.
                    variable *= SymmetrizedGaussianInSinc(pi/f1.canon[space].particle[body].lattice,n,m,f1.canon[space].particle[body].lattice * v);
                }
                if ( ! space )
                    variable *= amplitude;
                streams(f1,vector,0,space)[vc+Current*v0] = variable;
            }
            msp = mss;
        }
    }
}

void GTO( sinc_label f1, division vector ,inta spin,floata amplitude, inta *gamma, floata *delta ){
    inta space,vc,vsp,msp,mss=0,vn1,v,n,Current,v0;
    floata variable,alpha,y ;
    bodyType body ;
    if ( f1.name[vector].Current[spin] >= f1.name[vector].Partition )
        return;
    Current =  f1.name[vector].Current[spin]++;

    msp = 0;
    for  ( space =0; space < SPACE ; space++){
        if ( f1.canon[space].body != nada){
            v0 = vectorLen(f1, space);
            for ( vc = 0; vc < v0 ; vc++){
                vsp = 1;
                variable = 1.0;
                vn1 = vector1Len(f1,space);
                mss = msp ;
                for ( body = one ; body <= f1.canon[space].body ; body++){
                    v = (vc/vsp)%vn1-(vn1-1)/2;
                    vsp *= vn1;
                    n = (gamma[mss]);
                    alpha = delta[2*mss];
                    y = delta[2*mss+1];
                    mss += 1;
                    variable *= GaussianInSinc(pi/f1.canon[space].particle[body].lattice,n,alpha,y,f1.canon[space].particle[body].lattice * v );
                }
                if ( ! space )
                    variable *= amplitude;

                streams(f1,vector,spin,space)[vc+Current*v0] = variable;
            }
            msp = mss;
        }
    }
}


void loopDetails(  sinc_label f1,   division loopHeader){
      division loopElement;
    
    loopElement = loopHeader;
    while ( loopElement != nullName ){
        printf("loop %d %d %d\n",loopHeader, loopElement, f1.name[loopElement].loopNext );
        loopElement = f1.name[loopElement].loopNext;
    }
}


void chainDetails(  sinc_label f1,   division chainHeader){
      division chainElement;
    
    chainElement = chainHeader;
    while ( chainElement != nullName ){
        printf("chain %d %d %d\n", chainHeader, chainElement, f1.name[chainElement].chainNext );
        loopDetails(f1, chainElement);
        chainElement = f1.name[chainElement].chainNext;
    }
}

void linkDetails(  sinc_label f1,   division linkHeader){
      division linkElement;
    
    linkElement = linkHeader;
    while ( linkElement != nullName ){
        printf("link %d %d %d\n",linkHeader, linkElement, f1.name[linkElement].linkNext );
        chainDetails(f1, linkElement);
        linkElement = f1.name[linkElement].linkNext;
    }
}



inta  countLinesFromFile(  calculation *c1,   field f1,inta location, inta *ir, inta *ix ){
    *ix = 0;
    inta fi,cmpl,lines = 0,num;
    int ms = MAXSTRING;
    char line0[MAXSTRING];
    char name[MAXSTRING];
    char name2[MAXSTRING];
    char line2[MAXSTRING];
    char title [MAXSTRING];

    char *line = line0;
    inta FIT,iva ;
    if ( location == 0 )
        FIT = f1.i.files ;
    else if ( location == 1 )
        FIT = f1.i.filesVectorOperator ;
    else
        exit(1);
    
    for ( fi =0 ; fi < FIT; fi++){
            if ( location == 0 )
                strcpy(name,f1.i.fileList[fi] );
            else if ( location == 1)
                strcpy(name,f1.i.fileVectorOperator[fi] );
            else
                exit(0);
            FILE * fp = fopen(name,"r");
            if ( fp == NULL ) {
                printf("cannot find file: %s\n", name);
             //   printf("perhaps system.h:MAX_FILE is too small\n currently %d\n", MAX_FILE);
                continue;
                //  exit(0);
            }
            fgets(line, ms, fp);
            while(!feof(fp))
            {
                if (! comment(line))
                {
                    if (  strstr(name,c1->name) != NULL){
                        printf(" cannot name inputs same as outputs\n");
                        printf("%s \t %s\n", name, c1->name);
                        exit(0);
                    }

                    for ( cmpl = 0 ; cmpl < f1.i.cmpl ; cmpl++){
                        strcpy(line2,line);
                        tFromReadToFilename(NULL, line2,  name2, f1.i.cmpl-1,cmpl,title,&num);
                        iva = inputFormat(f1.f, name2, nullName, 2);
                        *ir = imax(*ir,iva);
                        *ix += iva;
                    }
                    lines++;
                }
                fgets(line, ms, fp);
            }
            if ( fi > MAX_FILE)
            {
                printf("too many files, increase MAX_FILE\n");
                exit(0);
            }
            fclose(fp);
    }
    return lines;
}

#if 1

inta defineTerms (  calculation * c,   field *f,   division head, inta memory){
    
    
            sinc_label *f1 = &f->f;
            inta term=0,i,productIndex=0,index, intvType[MAX_PRODUCT];
            //tied to bra.
            division prevLink = head,newTerm[MAX_PRODUCT],walkLink = head ;
            
while ( f1->name[walkLink].chainNext != nullName)
                                                walkLink =f1->name[walkLink].chainNext;


            for ( i = 0; i < c->i.termNumber ; i++)
            {
                    if ( c->i.terms[i].headFlag == 1 ){
                            term++;
                            if ( memory== -1 || term == memory ){
                                if ( term > 1 ){
                                    for ( index = 0 ;index < productIndex ; index++){
                                        if ( intvType[index] == -1 ){
                                            ///mult
                                            if ( CanonicalRank(f->f, newTerm[index],0) > 1 ) {
                                                printf("oops, only coded multiply for 1 canon");
                                                exit(0);
                                            }
//printf("%d %d\n", walkLink, f1->name[walkLink].loopNext);
//        printf("%d %d %d\n", newTerm[index], f1->name[newTerm[index]].chainNext , f1->name[f1->name[newTerm[index]].chainNext].loopNext );
                    //    walkLink = f1->name[walkLink].chainNext;
                                            f1->name[f1->name[walkLink].loopNext].multNext = f1->name[f1->name[newTerm[index]].chainNext].loopNext;
                                        } else if ( intvType[index] >= 0 ){
                                            ///add
                                            while ( f1->name[walkLink].chainNext != nullName)
                                                walkLink =f1->name[walkLink].chainNext;
                                            f1->name[walkLink].chainNext = f1->name[newTerm[index]].chainNext;
while ( f1->name[walkLink].chainNext != nullName)
                                                walkLink =f1->name[walkLink].chainNext;
                                        }
                                
                                    }
                                
                            }
                            productIndex = 0;

                            f1->name[prevLink].linkNext = anotherLabel(f1, all, nada);
                            prevLink = f1->name[prevLink].linkNext;
                            f1->name[prevLink].species = matrix;
                            walkLink = prevLink;
while ( f1->name[walkLink].chainNext != nullName)
                                                walkLink =f1->name[walkLink].chainNext;
                            printf("newTerm\t%d\t%s\n", term,c->i.terms[i].desc);
                        }
                 }

                if ( memory== -1 || term == memory ){
                    newTerm[productIndex] = anotherLabel(&f->f, all, nada);
                    switch ( c->i.terms[i].type){
                        case 1:
                            buildConstant(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex], c->i.terms[i].label[0], 0, real);
                            break;
                        case 2:
                            buildLinear(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex],  c->i.terms[i].label[0], 0, real);
                            break;
                        case 3:
                            buildSpring(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex],  c->i.terms[i].label[0], 0, real);
                            break;
                        case 4:
                            buildDeriv(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex], c->i.terms[i].label[0],0, real);
                            break;
                        case 5:
                            buildKinetic(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex],  c->i.terms[i].label[0], 0, real);
                            break;
                        case 6:
                            break;
                        case 7:
                            buildElement(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex],  c->i.terms[i].label[0], 0, real,c->i.terms[i].bra,c->i.terms[i].ket);
                            break;
                        case 8:
                            buildExternalPotential(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex],  c->i.terms[i].label,c->i.terms[i].embed, 0, real,c->i.terms[i].mu,c->i.terms[i].atom);
                            break;
                        case 9:
                            buildPairWisePotential(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl,newTerm[productIndex], c->i.terms[i].label, c->i.terms[i].embed,0, real,c->i.terms[i].mu);
                            break;
                        case 10:
                            assignDiagonalMatrix(c,f,c->i.terms[i].filename,newTerm[productIndex]);
                            break;
                        
                        case 11:
                            buildNumber(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex],  c->i.terms[i].label[0], 0, real, c->i.terms[i].omega);
                            break;
                        case 12:
                            buildCreate(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex],  c->i.terms[i].label[0], 0, real, c->i.terms[i].omega);
                            break;
                        case 13:
                            buildDestroy(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, newTerm[productIndex],  c->i.terms[i].label[0], 0, real, c->i.terms[i].omega);
                            break;
                    }
                        intvType[productIndex] = c->i.terms[i].headFlag;
                        productIndex += 1;
                
                        if ( productIndex == MAX_PRODUCT ){
                            printf("Issue,  too many products(!)");
                            exit(0);
                        }

                    }

            }
    
    {
        for ( index = 0 ;index < productIndex ; index++){
            if ( intvType[index] == -1 ){
                ///mult
                if ( CanonicalRank(f->f, newTerm[index],0) > 1 ) {
                    printf("oops, only coded multiply for 1 canon");
                    exit(0);
                }
        walkLink = f1->name[walkLink].chainNext;
    ///    printf("%d %d\n", walkLink, f1->name[walkLink].loopNext);
    ///    printf("%d %d %d\n", newTerm[index], f1->name[newTerm[index]].chainNext , f1->name[f1->name[newTerm[index]].chainNext].loopNext );

                f1->name[f1->name[walkLink].loopNext].multNext =f1->name[f1->name[newTerm[index]].chainNext].loopNext;
            } else if ( intvType[index] >= 0 ){
                ///add
                while ( f1->name[walkLink].chainNext != nullName)
                    walkLink =f1->name[walkLink].chainNext;
                f1->name[walkLink].chainNext = f1->name[newTerm[index]].chainNext;

            }
        }
        
    }
    
    
    
            //analyzeChainElement(*f1,f->f.name[head].linkNext,0);
            return term;
        }

#else
inta defineTerms(  calculation * c,   field *f,   division head, inta memory){
    sinc_label *f1 = &f->f;
    inta term=0,i;
    //tied to bra.
    division prevLink = head;
    
    for ( i = 0; i < c->i.termNumber ; i++)
    {
            if ( c->i.terms[i].headFlag ){
                term++;
                
                if ( memory== -1 || term == memory ){
                    f1->name[prevLink].linkNext = anotherLabel(f1, all, nada);
                    prevLink = f1->name[prevLink].linkNext;
                    f1->name[prevLink].species = matrix;
                    printf("newTerm\n");
                }
         }
        if ( memory== -1 || term == memory ){
            switch ( c->i.terms[i].type){
                case 1:
                    buildConstant(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink, c->i.terms[i].label[0], 0, real);
                    break;
                case 2:
                    buildLinear(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label[0], 0, real);
                    break;
                case 3:
                    buildSpring(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label[0], 0, real);
                    break;
                case 4:
                    buildDeriv(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink, c->i.terms[i].label[0],0, real);
                    break;
                case 5:
                    buildKinetic(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label[0], 0, real);
                    break;
                case 6:
//                    buildClampKinetic(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label, 0, real);
                    break;
                case 7:
                    buildElement(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label[0], 0, real,c->i.terms[i].bra,c->i.terms[i].ket);
                    break;
                case 8:
                    buildExternalPotential(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label,c->i.terms[i].embed, 0, real,c->i.terms[i].mu,c->i.terms[i].atom);
                    break;
                case 9:
                    buildPairWisePotential(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl,prevLink, c->i.terms[i].label, c->i.terms[i].embed,0, real,c->i.terms[i].mu);
                    break;
                case 10:
                    assignDiagonalMatrix(c,f,c->i.terms[i].filename,prevLink);
                    break;
            }
            analyzeChainElement(*f1,prevLink,0);
        }
    }
    return term;
}
#endif

inta defineCores(  calculation *c,  field *f ){
    inta i;
    runTime * rt = & c->rt;
    
    f->f.rt = rt;
    
#ifdef APPLE
    rt->NLanes = MAX_CORE;
#endif
#ifdef OMP
    if ( c->i.omp > MAX_CORE ){
        printf("lanes > MAX_CORE\n");
        c->i.omp = MAX_CORE;
    }
    if ( c->i.omp == -1 ){
#ifdef MKL
        if ( c->i.mkl < 1 )
        {
            printf("set parallel");
            exit(0);
        }
        
        c->i.omp = MAX_CORE/c->i.mkl;
#else
        c->i.omp = MAX_CORE;
#endif
    }

    rt->NLanes = c->i.omp;
#pragma omp parallel for private (i)
    for ( i = 0; i < MAX_CORE ; i++){
        rt->NSlot = omp_get_num_threads();
    }
    if ( rt->NLanes > rt->NSlot ){
        printf("decrease lanes or increase your number of OMP cores\n");
        rt->NLanes = rt->NSlot;
    }
    
#ifdef MKL
    if ( rt->NSlot < c->i.omp*c->i.mkl )
    {
        printf("not enough slots for omp\n" );
        c->i.omp = rt->NSlot/c->i.mkl;
    }
    rt->NParallel = c->i.mkl;
    printf("parallel \t %d\n", rt->NParallel);

#endif
    printf("lanes \t %d\n",  rt->NLanes);

#endif
    return 0;
}

inta sumTo2(  sinc_label f1,double scalar, inta space,  blockType bl,   division mat,inta ms,   division sum,inta spin){
      blockType bol;
    inta I1,I2, I3, I4,r,N1;
    double value;
    for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
        for ( bol = tv1 ; bol <= tv2 ; bol++){
            {
                N1 = vector1Len(f1, space);
                floata * stream = streams(f1,sum,spin,space)+N1*N1*N1*N1*(CanonicalRank(f1, sum, spin)+r);
                for ( I1 = 0 ; I1 < N1 ; I1++)//body 0
                    for ( I2 = 0 ; I2 < N1 ; I2++)//body 0
                        for ( I3 = 0 ; I3 < N1 ; I3++)//body 1
                            for ( I4 = 0 ; I4 < N1 ; I4++)//body 1
                            if ( bol == bl)
                            {
                                value = 0.;
                                if ( bol == tv1 ){
                                    value  = streams(f1,mat,ms,space)[ I1*N1+I2 + r*N1*N1 ] * delta(I3-I4);
                                }else if (bol == tv2){
                                    value  = streams(f1,mat,ms,space)[ I3*N1+I4 + r*N1*N1 ] * delta(I1-I2);
                                }
                                if ( space == 0 )
                                    value *= scalar;
                                stream[ (I1+I3*N1)+ ( I2+I4*N1)*N1*N1 ] = value;
                            }
            }
        }
    return 0;
}

inta sumTo3(  sinc_label f1,double scalar,inta space,    blockType bl,   division mat,inta ms,   division sum,inta spin){

    inta n2[SPACE];
    length(f1, sum,n2);

    if ( bodies ( f1, sum ) == three){

        if ( bodies ( f1, mat ) == one ){
              blockType bol;
            inta I1,I2, I3, I4,I5,I6,r;
            inta n1[SPACE];
            length1(f1,n1);
            double value;

            for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
                for ( bol = tv1 ; bol <= tv3 ; bol++){
                    {
                        floata * stream = streams(f1,sum,spin,space)+n2[space]*(CanonicalRank(f1, sum, spin)+r);
                        if ( bol == bl )
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                            {
                                                value = 0;
                                                if ( bol == tv1 ){
                                                    value  = streams(f1,mat,ms,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6);
                                                }else if ( bol == tv2 ) {
                                                    value  = streams(f1,mat,ms,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6);
                                                }else if ( bol == tv3 ) {
                                                    value  = streams(f1,mat,ms,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4);
                                                }
                                                if ( space == 0 )
                                                value *= scalar;
                                                stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]] = value;
                                            }
                    }
                }
        }else if ( bodies ( f1, mat ) == two ){
              blockType pair ;
            inta I1,I2, I3, I4,I5,I6,r,ve;
            inta n1[SPACE];

            length1(f1, n1);

            double value;

            for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
                for ( pair = e12 ; pair <=e23  ; pair++){
                    {
                        floata * stream = streams(f1,sum,spin,space)+n2[space]*(CanonicalRank(f1, sum, spin)+r);
                        if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
                            printf("sumTo3, failed to allocate enough for Two body \n");
                            exit(0);
                        }
                        if ( pair == bl )
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)

                                            {
                                                value = 0;
                                                if ( pair == e12 ){
                                                    value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6);
                                                }else if ( pair == e13 ) {
                                                    value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4);
                                                }else if ( pair == e23 ) {
                                                    value  = streams(f1,mat,ms,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2);
                                                }
                                                ve = (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space];
                                                //                                                printf("%f %lld %lld\n", value,ve,n2[space]);
                                                //                                                fflush(stdout);
                                                if ( space == 0 )
                                                value *= scalar;
                                                stream[ ve ] = value;
                                                //                                                printf("x");
                                                //                                                fflush(stdout);

                                            }
                    }
                }
        }
        else {
            printf("Yo!");
            exit(0);
        }

    }
    return 0;
}








inta sumTo4(  sinc_label f1,double scalar,inta space,   blockType bl,   division mat,inta ms,   division sum,inta spin){

    inta n2[SPACE];
    length(f1, sum,n2);
    if (bodies(f1,sum) == four && f1.rt->calcType == electronicStuctureCalculation){

    if ( bodies ( f1, mat ) == one ){
          blockType bol ;
        inta I1,I2, I3, I4,I5,I6,I7,I8,r;
        inta n1[SPACE];
        length1(f1,n1);
        double value;

        for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
            for ( bol = tv1 ; bol <= tv4 ; bol++){
                {
                    floata * stream = streams(f1,sum,spin,space)+n2[space]*(CanonicalRank(f1, sum, spin)+r);
                    if ( bol == bl )
                    for ( I1 = 0 ; I1 < n1[space] ; I1++)
                        for ( I2 = 0 ; I2 < n1[space] ; I2++)
                            for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                    for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                        for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                            for ( I7 = 0 ; I7 < n1[space] ; I7++)
                                                for ( I8 = 0 ; I8 < n1[space] ; I8++)
                                                        
                                                {
                                                    value = 0;
                                                    if ( bol == tv1 ){
                                                        value  = streams(f1,mat,ms,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6)*delta(I7-I8);
                                                    }else if ( bol == tv2 ) {
                                                        value  = streams(f1,mat,ms,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6)*delta(I7-I8);
                                                    }else if ( bol == tv3 ) {
                                                        value  = streams(f1,mat,ms,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I7-I8);
                                                    }else if ( bol == tv4 ) {
                                                        value  = streams(f1,mat,ms,space)[ I7*n1[space]+I8 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I5-I6);
                                                    }
                                                    if ( space == 0 )
                                                    value *= scalar;
                                                    stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
                                                }
                }
            }
    }else if ( bodies ( f1, mat ) == two ){
          blockType pair;
        inta I1,I2, I3, I4,I5,I6,I7,I8,r;
        inta n1[SPACE];
        length1(f1, n1);

        double value;

        for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
            for ( pair = e12 ; pair <= e34 ; pair++){
                {
                    floata * stream = streams(f1,sum,spin,space)+n2[space]*(CanonicalRank(f1, sum, spin)+r);
                    if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
                        printf("sumTo4, failed to allocate enough for Two body \n");
                        exit(0);
                    }
                    if ( pair == bl )
                    for ( I1 = 0 ; I1 < n1[space] ; I1++)
                        for ( I2 = 0 ; I2 < n1[space] ; I2++)
                            for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                    for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                        for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                            for ( I7 = 0 ; I7 < n1[space] ; I7++)
                                                for ( I8 = 0 ; I8 < n1[space] ; I8++)

                                                {
                                                    //                                            0    e12,     1,3     2,4
                                                    //                                            1    e13,     1,5     2,6
                                                    //                                            2    e23,     3,5     4,6
                                                    //                                            3    e14,     1,7     2,8
                                                    //                                            4    e24,     3,7     4,8
                                                    //                                            5    e34      5,7     6,8


                                                    if ( pair == e12 ){
                                                        value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6)*delta(I7-I8);
                                                    }else if ( pair == e13 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I7-I8);
                                                    }else if ( pair == e23 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I7-I8);
                                                    }else if ( pair == e14 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I7) + (I2*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I5-I6);
                                                    }else if ( pair == e24 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I3*n1[space]+I7) + (I4*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I5-I6);
                                                    }else if ( pair == e34 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I5*n1[space]+I7) + (I6*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I3-I4);
                                                    }else {
                                                        printf ("rails!\n");
                                                        exit(0);
                                                    }
                                                    if ( space == 0 )
                                                    value *= scalar;
                                                    stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
                                                }
                }
            }
    }
    else {
        printf("sumTo4!");
        exit(0);
    }
    }
    return 0;
}


inta assignCores(  sinc_label f1, inta parallel ){
#ifdef OMP
#ifdef MKL
    inta nSlot = f1.rt->NSlot;
    inta nParallel = f1.rt->NParallel;
    inta nLanes = f1.rt->NLanes;
  
    inta omp;
    if ( parallel == 2 ){
        ///for serious serial MKL operations
        omp_set_num_threads(1);
        mkl_set_num_threads(nParallel*nLanes);
    }
    else if ( parallel == 1 ){
        ///mixed split between lanes and MLK parallelism
        omp_set_num_threads(nLanes);
        mkl_set_num_threads(nParallel);
    }
    else if ( parallel == 0 ){
        ///i suspect HDF5 behaves well here
        omp_set_num_threads(nLanes*nParallel);
        mkl_set_num_threads(1);
    }
#else
    ///not MKL
      inta nSlot = f1.rt->NSlot;
      inta nLanes = f1.rt->NLanes;
    
      inta omp;
      if ( parallel == 2 ){
          ///not recommended
          omp_set_num_threads(1);
      }
      else if ( parallel == 1 ){
          omp_set_num_threads(nLanes);
      }
      else if ( parallel == 0 ){
          omp_set_num_threads(nLanes);
      }
#endif
#else
    ///nothing to change
#endif
    return 0;

}


  basisElement_label grabBasis (  sinc_label f1, inta space, inta particle, inta elementIndex){
    ///new design decision,  drop higher SPACEs
      basisElement_label x;
    double length=0;
    length = f1.canon[space].particle[particle].lattice;
    //x.component = f1.canon[space].space +1;
    x.basis = f1.canon[space].basis;
    x.grid = f1.canon[space].count1Basis;
    if ( x.grid %2 == 1 )
        x.index = elementIndex ;
    else
        x.index = elementIndex ;
    x.index2= 0;
    x.length = length;
    ///left edge of grid
    x.origin = f1.canon[space].particle[particle].origin;
    return x;
}

double xOneBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label f2,   division out,inta s2,inta oldPeriodic){
    if ( ! allowQ(f1.rt, blockTransferBasisblock)){
        printf("blockTransferBasisblock Allow!\n");
        fflush(stdout);
        exit(0);
    }
    inta i,l,r,rank=0,p;
    inta n1[SPACE],N2;
    length1(f1,n1);
    inta n2[SPACE];
    length1(f2,n2);
    floata *band = myStreams(f1, bandBasis, rank);
    inta L1;
    f2.name[out].Current[s2] = 0;
    
    if ( f1.canon[space].body != nada){
            N2 = n2[space];
            L1 = n1[space];
            
            for ( l = 0; l < L1 ; l++){
                for ( i = 0 ; i < N2 ; i++)
                    for ( p = 1 ;p <= 1 ; p++)
                        {
                            ///backwards bc of dgemv
                            band[(p-1)*L1*N2+l*N2 + i] = BoB (grabBasis(f1, space, p,l),grabBasis(f2, space, p, i) );
                        }
            
            }
           ///add orthogonal directions...
            for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                ///band_ii'   Vector_i
                cblas_dgemv(CblasColMajor, CblasNoTrans, N2, L1, 1.,band,N2,streams( f1, vector1, s1,space )+r*L1,1, 0.,streams( f2, out, s2,space)+r*N2,1);
            }
        }
    f2.name[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    return 0.;
}

double xTwoBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label f2,   division out,inta s2,inta oldPeriodic){
    inta i,l,r,rank=0,p;
    inta n1[SPACE],N2;
    length1(f1,n1);
    inta n2[SPACE];
    length1(f2,n2);
    floata *band = myStreams(f1, bandBasis, rank);
    if ( ! allowQ(f1.rt, blockTransferBasisblock)){
        printf("blockTransferBasisblock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    
    floata *buffer = band+n1[0]*n2[0]*2;

    inta L1;
    f2.name[out].Current[s2] = 0;
    
    
        if ( f1.canon[space].body != nada){
            N2 = n2[space];
            L1 = n1[space];
            
            for ( l = 0; l < L1 ; l++){
                for ( i = 0 ; i < N2 ; i++)
                    for ( p = 1 ;p <= 2 ; p++)
                        {
                            band[(p-1)*L1*N2+i*L1 + l] = BoB (grabBasis(f1, space, p,l),grabBasis(f2, space, p, i) );
                        }
                }
            
            
            //add orthogonal directions...
            for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                //band_ii' Vector_ij = Vector_i'j
                cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans,N2,L1,L1,1.    ,band    ,L1,streams( f1, vector1, s1,space )+r*L1*L1,L1, 0.,  buffer  , N2   );
               
                //Vector_i'j band_jj' = Vector_i'j'
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans  ,N2,N2,L1,1.,buffer  ,N2,band+L1*N2  ,L1, 0.,  streams(f2, out, s2,space)+r*N2*N2 , N2   );
                
            }
        }
    f2.name[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    return 0.;
}


double xThreeBand (  sinc_label f1, inta space,  division vector1 ,inta s1,   sinc_label f2,   division out,inta s2,inta oldPeriodic){
    inta i,l,r,rank=0,p,k;
    inta n1[SPACE],N2;
    length1(f1,n1);
    inta n2[SPACE];
    length1(f2,n2);
    floata *band = myStreams(f1, bandBasis, rank);
    floata *buffer = band+n1[0]*n2[0]*3;
    floata *buffer2= buffer+n2[0]*n2[0]*n2[0];
    if ( ! allowQ(f1.rt, blockTransferBasisblock)){
        printf("blockTransferBasisblock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    inta L1;
    f2.name[out].Current[s2] = 0;
    
    
        if ( f1.canon[space].body != nada){
            N2 = n2[space];
            L1 = n1[space];
            
            for ( l = 0; l < L1 ; l++){
                for ( i = 0 ; i < N2 ; i++)
                    for ( p = 1 ;p <= 3 ; p++)
                        {
                            band[(p-1)*L1*N2+i*L1 + l] = BoB (grabBasis(f1, space, p,l),grabBasis(f2, space, p, i) );
                        }
            }
            
            
            //add orthogonal directions...
            for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                
                //particle 1
                //                //band_ii' Vector_ijk = Vector_i'jk

                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,N2,L1*L1,L1,1. ,band,L1,streams( f1, vector1, s1,space )+r*L1*L1*L1,L1, 0.,  buffer  , N2   );
                //factor third particle,
                // Vector_i'j--k band_jj'=> Vector_i'j'--k
                for ( k = 0 ; k < L1 ;k++)
                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans  ,N2,N2,L1,1.,buffer+N2*L1*k,N2,band+1*L1*N2 ,L1, 0.,  buffer2 + N2*N2*k, N2   );
                //Vector_i'j'k band_kk' = Vector_i'j'k'
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans  ,N2*N2,N2,L1,1. ,buffer2,N2*N2,band+2*L1*N2 ,L1, 0.,  streams( f2, out, s2,space)+r*N2*N2*N2 , N2*N2   );
                                
            }
        }
    f2.name[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    return 0.;
}


double xFourBand (  sinc_label f1,inta space,   division vector1 ,inta s1,   sinc_label f2,   division out,inta s2,inta oldPeriodic){
    inta i,l,r,rank=0,p,k;
    inta n1[SPACE],N2;
    length1(f1,n1);
    inta n2[SPACE];
    length1(f2,n2);
    floata *band = myStreams(f1, bandBasis, rank);
    floata *buffer = band+n2[0]*n2[0]*4;
    floata *buffer2= buffer+n2[0]*n2[0]*n2[0]*n2[0];
    if ( ! allowQ(f1.rt, blockTransferBasisblock)){
        printf("blockTransferBasisblock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    inta L1;
    f2.name[out].Current[s2] = 0;
    
    
        if ( f1.canon[space].body != nada){
            N2 = n2[space];
            L1 = n1[space];
            
            for ( l = 0; l < L1 ; l++){
                for ( i = 0 ; i < N2 ; i++)
                    for ( p = 1 ;p <= 4 ; p++)
                        {
                            band[(p-1)*L1*N2+i*L1 + l] = BoB (grabBasis(f1, space, p,l),grabBasis(f2, space, p, i) );
                        }
            }
            
            
            //add orthogonal directions...
            for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                
                //particle 1
                //                //band_i'i Vector_ijkl = Vector_i'jkl

                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,N2,L1*L1*L1,L1,1. ,band         ,L1,streams( f1, vector1, s1,space )+r*L1*L1*L1*L1,L1, 0.,  buffer  , N2   );
                //factor third-fourth particles,
                // Vector_i'j--kl band_jj'=> Vector_i'j'--kl
                for ( l = 0 ; l < L1 ;l++)
                    for ( k = 0 ; k < L1 ;k++)
                        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans  ,N2,N2,L1,1.,buffer+L1*N2*k+L1*L1*N2*l      ,N2,band+1*L1*N2 ,L1, 0.,  buffer2 + N2*N2*k+L1*N2*N2*l , N2   );
                
                //factor fourth particles,
                // Vector_i'j'k--l band_kk'=> Vector_i'j'k'--l
                for ( l = 0 ; l < L1 ;l++)
                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans  ,N2*N2,N2,L1,1.,buffer2+N2*N2*L1*l      ,N2*N2,band+2*L1*N2 ,L1, 0.,  buffer + N2*N2*N2*l , N2*N2   );
                
                //Vector_i'j'k'l band_ll' = Vector_i'j'k'l'
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans  ,N2*N2*N2,N2,L1,1. ,buffer                  ,N2*N2*N2,band+3*L1*N2 ,L1, 0.,  streams( f2, out, s2,space)+r*N2*N2*N2*N2 , N2*N2*N2   );
                                
            }
        }
    f2.name[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    return 0.;
}





/**
 *The printing output in D.kry
 *
 *@param c  The general parameters
 *@param f1 The container
 *@param Hamiltonian The all links will be printed
 *@param vector division of the input vector
 */
double printExpectationValues (  calculation *c,   sinc_label  f1 ,  division Hamiltonian  ,   division vector){
    mea totx,me,ov,edge;
    inta space,ed,i,spacer;
      bodyType body;
    if ( !allowQ(f1.rt, blockPrintStuffblock))
        return 0.;
    
    inta o;
    division op = defSpiralMatrix(&f1, Hamiltonian);
    division OpSpiral = defSpiralMatrix(&f1, Iterator);
    for (o = 0; f1.name[OpSpiral+o].species == matrix ; o++){
        printf("\nterm%d\n", o+1);
        analyzeChainElement(f1, OpSpiral+o ,0);
    }

    totx = 0.;
    printf("\n======Expectation========\n");
    if (CanonicalRank(f1, vector, 1))
        printf("\t%d (%d + i%d)\n", vector, CanonicalRank(f1, vector, 0), CanonicalRank(f1, vector, 1));
    else
        printf("\t(%d)\t%d\n",  CanonicalRank(f1, vector, 0),vector-eigenVectors);
    ov = pMatrixElement( f1, vector, 0, nullOverlap, 0, vector,0);
    printf("------Edges------\n\n");
    division header = anotherLabel(&f1, 0, nada);
    division mem    = anotherLabel(&f1, 0, one);
    char* desc [] = {"fourier","negative","positive"};
    for ( space = 0; space < SPACE ; space++){
        if ( f1.canon[space].basis != DiracDeltaElement)
        for (body = one ; body <=  f1.canon[space].body ; body++ )
            for ( ed = 0 ; ed < 3 ; ed++){
                f1.name[header].species = eikon;

                f1.name[header].loopNext = mem;
                f1.name[mem].species = eikonOuter;
                f1.name[mem].Current[0] = 1;
                f1.name[mem].space[space].body = one;

                for ( spacer = 0; spacer < SPACE ;spacer++)
                    f1.name[mem].space[spacer].block = id0;
                f1.name[mem].space[space].block = (blockType)body;
                zero(f1 , mem,0);
            switch ( ed ){
                case 0:
                    for ( i = 0; i < vector1Len(f1, space) ; i++)
                        streams(f1,mem,0,space)[i] = sign(i)*1./vector1Len(f1, space);
                    break;
                case 1:
                    streams(f1,mem,0,space)[0] = 1;
                    break;
                case 2:
                    streams(f1,mem,0,space)[vector1Len(f1, space)-1] = 1;
                    break;
            }
                edge = pMatrixElement(f1, vector, 0, header, 0, vector,0);
                printf("\t%d\t%d\t%s \t%1.12f\n", space+1,body,desc[ed], edge/ov);
            }
    }
        
    printf("------Terms------\n");
    inta sp,scr=0,cr =0;;
    inta terms = 0,oo=0;
    
    
    for (o = 0; f1.name[op+o].species == matrix ; o++)
    {
        
        while ( terms <= o && oo < c->i.termNumber ){
            if ( c->i.terms[oo].headFlag == 1 ) {
                terms++;
                printf("%s ", c->i.terms[oo].desc);
            }
            oo++;
        }
        me = 0.0;
        for ( sp = 0 ; sp < spins(f1,vector); sp++)
            me += pMatrixElement( f1, vector, sp, op+o, 0, vector,sp);
        cr = CanonicalOperator(f1, op+o, 0);
        scr += cr;
        printf("\t(%d)\t%6.12f\n", cr,creal(me/ov));
        fflush(stdout);
        totx += me/ov;
    }
    f1.name[vector].value.value = totx;
   {
#ifdef COMPLEXME
            printf("sum\t(%d)\t%6.12f\tI %6.6f\n",scr,creal(totx),cimag(totx));
#else
            printf("sum\t(%d)\t%6.12f\n",scr,(totx));
#endif
    }
        
    return totx;
}


/**
 *The matrix element calculator using multiply
 *
 *@param rank      OMP rank
 *@param f1          container
 *@param bra       division on left of operator
 *@param bspin  bra's spin
 *@param mat       A single link or Hamilonian Term
 *@param mspin mat's spin
 *@param ket        division on right of operator
 *@param kspin  ket's spin
 */
double tMatrixElements ( inta rank,  sinc_label  f1 , division bra, inta bspin,  division mat, inta mspin,  division ket, inta kspin ){
    mat = f1.name[mat].name;
    division holder ;
    inta holderRank, holderSpin;
    inta k,l,e,space;
    double prod,ME=0.;
    inta ca;
    ///outputFormat(f1, stdout, bra, 0);
    if ( rank )
        if ( ! allowQ(f1.rt, blockParallelMatrixElementblock)){
            printf("blockParallelMatrixElementblock Allow!\n");
            fflush(stdout);
            exit(0);
        }
    
    
    if ( mat == nullName || f1.name[mat].name == nullName)
        return 0.;
    
    if (mat == nullOverlap  ){
        ca = imax(1, f1.name[overlap].Current[mspin]);
    }
    else if ( f1.name[mat].species == matrix || f1.name[mat].species == vector || f1.name[mat].species >= eikon)
        ca = CanonicalOperator(f1,mat, mspin);
    else
        return 0.;
    
                for ( k = 0 ; k < CanonicalRank(f1, ket, kspin);k++){
                    for ( l = 0 ; l < ca;l++){
                        if ( mat == nullOverlap ){
                            if ( f1.name[overlap].Current[mspin] == 0 ){
                                holder = ket;
                                holderRank = k;
                                holderSpin = kspin;
                            }else {
                                tHX(rank, f1, overlap, l, mspin, 1., ket, k, kspin,canonicalmeVector, 0, rank);
                                
                                holder = canonicalmeVector;
                                holderRank = 0;
                                holderSpin = rank;
                            }
                        }
                        else {
                            tHX(rank, f1, mat, l, mspin, 1., ket, k, kspin,canonicalmeVector, 0, rank);
                            f1.name[canonicalmeVector].Current[0] = 1;
                            holder = canonicalmeVector;
                            holderRank = 0;
                            holderSpin = rank;
                            
                        }
                        for ( e = 0 ; e < CanonicalRank(f1, bra, bspin);e++){
                            prod = 1;
                            for ( space = 0 ; space < SPACE ; space++)
                                if ( f1.canon[space].body != nada){
                                    prod *= tDOT(rank, f1,space,CDT, bra, e, bspin,CDT, holder, holderRank, holderSpin);
                                }
                            ME += prod;
                        }
                    }
                }
        return ME;
}
            
        


/**
 *general matrix - vector multiply per dimension
 *
 *@param rank      OMP rank
 *@param f1          container
 *@param equals output division
 *@param espin  output spin
 *@param left  matrix
 *@param l  one of the matrices canonical ranks indexed
 *@param lspin matrix spin
 *@param right input division
 *@param r one of the input canonical ranks indexed
 *@param rspin input spin
 */
inta tGEMV (inta rank,    sinc_label  f1,   division equals, inta e, inta espin, floata scalar, division left,inta l,inta lspin,   division right,inta r, inta rspin ){
    if ( header(f1, left ) != header(f1, right ) ){
        printf("Two Head types GEMV\n %d %d %d %d %d %d",equals,header(f1, equals ) ,left,header(f1, left ) ,right,header(f1, right ) );
        exit(1);
    }
    if ( rank ) {
        if ( ! allowQ(f1.rt, blockParallelMultiplyblock)){
            printf("blockParallelMultiplyblock allow!\n");
            fflush(stdout);
            exit(0);
        }
    }
    
    division inT,outT,midT,laterT,initT,bufferT;
    inta space,inR,outR,inS,outS,midR,midS,laterR,laterS,initR,initS,bufferR,bufferS;
    f1.name[canonicalmvVector].Current[rank] = 0;
    f1.name[canonicalmv2Vector].Current[rank] = 1;
    f1.name[canonicalmv3Vector].Current[rank] = 0;
    if ( species(f1,left ) >= eikon && species(f1, right ) == vector ){
        char  in,out;
            in = 1;
            out = 1;
        
        initT = right;
        initR = r;
        initS = rspin;
        
        outT = equals;
        outR = e;
        outS = espin;
        
        midT = canonicalmv3Vector;
        midR = 0;
        midS = rank;
        
        laterT = canonicalmv3Vector;
        laterR = 1;
        laterS = rank;

        
        inT = canonicalmv3Vector;
        inR = 2;
        inS = rank;

        bufferT = canonicalmv3Vector;
        bufferR = 3;
        bufferS = rank;


            
                
                inta i;
        double * midP ;
        double * laterP ;
        double * initP  ;
        double * inP;
        double * outP ;
        double * bufferP;
#if VERBOSE
                    printf("p%d\n",su);
#endif
                    for ( space = 0; space < SPACE ; space++)
                    if ( f1.canon[space].body != nada ){
                        
                        division su = f1.name[left].loopNext;
                        division su0 = su;
                        inta timer = 0,xlxl=0;
                        ///PRODUCT!
                    
                        bodyType bd = Bodies(f1, right,space);
                        inta N1 = outerVectorLen(f1, one,space);
                        inta N2 = vectorLen(f1, space);

                            inta firstFlag =1 ;
                    
                    
                        
                        midP = streams(f1, midT, midS,space)+midR*N2;
                        laterP = streams(f1, laterT, laterS,space)+laterR*N2;
                        initP  = streams(f1, initT, initS,space)+initR*N2;
                        bufferP  = streams(f1, bufferT, bufferS,space)+bufferR*N2;
                        inP = NULL;
                        outP = streams(f1, outT, outS,space)+outR*N2;
                        while ( su != nullName ){

                        if ( firstFlag ){
                            inP = initP;
                            for ( i = 0 ; i < N2 ; i++){
                                outP[i] = 0.;
                                bufferP[i] = 0.;
                            }
                            firstFlag = 0;
                        }
#if 1
                        printf("in%d %f %d\n",su, cblas_dnrm2(N2, inP, 1),N1);
#endif
                    if ( f1.name[su].space[space].block == id0 )
                        xlxl = 1;
                    else
                        switch(species(f1,su)){
                            case eikonDiagonal:
                            case eikonKinetic:
                            case eikonConstant:
                            case eikonDeriv:
                            case eikonLinear:
                            case eikonSpring:
                            case eikonElement:
                            case eikonOuter:
                            case eikonSplit:
                            case eikonCreate:
                            case eikonDestroy:
                                xlxl = 1;
                                break;
                            case eikonSemiDiagonal:
                                switch(Bodies(f1, su,space)){
                                    case one:
                                        xlxl = 0;
                                        break;
                                    case two:
                                        xlxl = 4;
                                        break;
                                default:
                                    break;

                                }
                                break;
                            case eikonOffDiagonal:
                                switch(Bodies(f1, su,space)){
                                    case one:
                                        xlxl = 2;
                                        break;
                                    case two:
                                        xlxl = 4;
                                        break;
                                   default:
                                       break;

                                }
                                break;
                            default:
                                break;

                    }
                            for ( i = 0 ; i < N2 ; i++){
                                bufferP[i] = 0.;
                            }

                    for ( timer = 0 ; timer < xlxl ; timer++){
                        for ( i = 0 ; i < N2 ; i++){
                            midP[i] = 0.;
                            laterP[i] = 0.;
                        }
                        double flow = 1.;
                        if (space == 0 )
                            flow *= scalar;
                        
                         ///the notion is to buffer on mid and accumlate on out
                        inta N3 = alloc(f1, su, space);
                        double * suP  = streams(f1, su , lspin, space)+l*N3;
                        if ( f1.name[su].space[space].block == id0 ){
                            topezOp(0,0, bd,f1.name[su].space[space].act,tv1, tv1,N1,inP,0, laterP);
                        }else
                        if ( Bodies ( f1,su,space) == one ){
                             if ( species(f1,su) == eikonDeriv){
                                flow *= *suP;
                                topezOp(0,f1.canon[space].particle[f1.name[su].space[space].block].lattice, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,1, laterP);
                            } else
                        if ( species(f1,su) == eikonKinetic ){
                                flow *= *suP;
                                topezOp(0, f1.canon[space].particle[f1.name[su].space[space].block].lattice, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,2, laterP);
                            }
                        else if ( species(f1,su) == eikonConstant){
                            ///action can happen!
                            flow *= *suP;
                            topezOp(0,0, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,0, laterP);
                        }
                        else if ( species(f1,su) == eikonLinear){
                            flow *= *suP;
                            ///relative to grid only
                            floata center = 0.5*(f1.canon[space].count1Basis-1)*f1.canon[space].particle[f1.name[su].space[space].block].lattice;
                            topezOp(center, f1.canon[space].particle[f1.name[su].space[space].block].lattice, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,-1, laterP);
                        }
                        else if ( species(f1,su) == eikonCreate){
                            flow *= *suP*sqrt(OMEGA*0.500);
                            ///relative to grid only
                            floata center = 0.5*(f1.canon[space].count1Basis-1)*f1.canon[space].particle[f1.name[su].space[space].block].lattice;
                            topezOp(center, f1.canon[space].particle[f1.name[su].space[space].block].lattice, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,-1, laterP);
                            topezOp(0,f1.canon[space].particle[f1.name[su].space[space].block].lattice, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,1, midP);
                            cblas_daxpy(N2, -1./OMEGA, midP, 1, laterP, 1);

                        }
                        else if ( species(f1,su) == eikonDestroy){
                            flow *= *suP*sqrt(OMEGA*0.500);
                            ///relative to grid only
                            floata center = 0.5*(f1.canon[space].count1Basis-1)*f1.canon[space].particle[f1.name[su].space[space].block].lattice;
                            topezOp(center, f1.canon[space].particle[f1.name[su].space[space].block].lattice, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,-1, laterP);
                            topezOp(0,f1.canon[space].particle[f1.name[su].space[space].block].lattice, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,1, midP);
                            cblas_daxpy(N2, +1./OMEGA, midP, 1, laterP, 1);
                        }

                        else if ( species(f1,su) == eikonSpring){
                            flow *= *suP;
                            ///relative to grid only
                            floata center = 0.5* (f1.canon[space].count1Basis-1)*f1.canon[space].particle[f1.name[su].space[space].block].lattice;
                            topezOp(center, f1.canon[space].particle[f1.name[su].space[space].block].lattice, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,-2, laterP);

                        }
                        else if ( species(f1,su) == eikonElement){
                            flow *= *suP;
                            for ( i = 0 ;i < N2 ; i++)
                                laterP[i] = 0.;

                            ///not a spatial state, this is an electronic level... no permutations here...just dropout a vector.

                            laterP[f1.name[su].space[space].bra] = inP[f1.name[su].space[space].ket];

                        }
                        else if ( species(f1,su) == eikonOuter){
                            inta n1[MAXBODY] ,iii,ii,i;
                            double suX;
                            switch(f1.name[su].space[space].block){
                                case tv1:
                                    n1[0] = 1;
                                    n1[1] = N1;
                                    n1[2] = N1*N1;
                                    n1[3] = N1*N1*N1;

                                    break;
                                case tv2 :
                                    n1[1] = 1;
                                    n1[0] = N1;
                                    n1[2] = N1*N1;
                                    n1[3] = N1*N1*N1;

                                    break;
                                case tv3 :
                                    n1[1] = 1;
                                    n1[2] = N1;
                                    n1[0] = N1*N1;
                                    n1[3] = N1*N1*N1;

                                    break;
                                case tv4 :
                                    n1[3] = 1;
                                    n1[1] = N1;
                                    n1[2] = N1*N1;
                                    n1[0] = N1*N1*N1;
                                    break;

                                default:
                                    break;
                            }
                            for ( i = 0 ;i < N2 ; i++)
                                laterP[i] = 0.;
                            switch ( bd ){

                                case one:
                                    cblas_dcopy(N1, suP, 1, laterP, 1);
                                    flow *= cblas_ddot(N1, inP, 1, suP, 1);
                                    break;
                                case two:
                                    suX = 0.;
                                    for ( i = 0 ; i < N1 ; i++){
                                        cblas_dcopy(N1, suP, 1, laterP+i*n1[1], n1[0]);
                                        suX += cblas_ddot(N1, inP+i*n1[1], n1[0], suP, 1);
                                    }
                                    flow *= suX;
                                    break;
                                case three:
                                    suX = 0.;
                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( ii = 0 ; ii < N1 ; ii++){
                                            cblas_dcopy(N1, suP, 1, laterP+i*n1[1]+ii*n1[2], n1[0]);
                                            suX += cblas_ddot(N1, inP+i*n1[1]+ii*n1[2], n1[0], suP, 1);
                                        }
                                    flow *= suX;
                                    break;
                                case four:
                                    suX = 0.;
                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( ii = 0 ; ii < N1 ; ii++)
                                            for ( iii = 0 ; iii < N1 ; iii++){
                                                cblas_dcopy(N1, suP, 1, laterP+i*n1[1]+ii*n1[2]+iii*n1[3], n1[0]);
                                                suX += cblas_ddot(N1, inP+i*n1[1]+ii*n1[2]+iii*n1[3], n1[0], suP, 1);
                                        }
                                    flow *= suX;
                                    break;
                                default:
                                    break;

                            }
                        }

                        else
                        if ( species(f1,su) == eikonDiagonal ){
                            flow *= 1.;
                            diagonalOp(bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP, suP,laterP);
                        }else if ( species(f1,su) == eikonOffDiagonal ){
                            if ( timer == 0 ){
                                flow *= -1;
                                topezOp(0,1., bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP,1, midP);
                                diagonalOp(bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,midP, suP,laterP);
                            }
                            else if ( timer == 1 ){
                                flow *= 1;
                                diagonalOp(bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP, suP,midP);
                                topezOp(0,1., bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,midP, 1,laterP);
                            }
                    }else if ( species(f1,su) == eikonSplit ){
                        
                        if ( bd == one )
                            cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1, inP,1, 0.,laterP, 1  );
                        else if ( bd == two ){
                            inta i;
                            switch (f1.name[su].space[space].block){
                                case tv1:
                                    for ( i = 0 ; i < N1 ; i++)
                                        cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1, inP+N1*i,1, 0.,laterP+N1*i, 1  );
                                    break;
                                case tv2:
                                    for ( i = 0 ; i < N1 ; i++)
                                        cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1, inP+i,N1, 0.,laterP+i, N1  );
                                    break;
                                default:
                                    break;
                            }
                        }
                        else if ( bd == three ){
                            inta i,i2;
                            switch (f1.name[su].space[space].block){
                                case tv1:
                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( i2 = 0 ; i2 < N1 ; i2++)
                                            cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1, inP+N1*i+N1*N1*i2,1, 0.,laterP+N1*i+N1*N1*i2, 1  );
                                    break;
                                case tv2:
                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( i2 = 0 ; i2 < N1 ; i2++)
                                            cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1, inP+i+i2*N1*N1,N1, 0.,laterP+i+i2*N1*N1, N1  );
                                    break;
                                case tv3:
                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( i2 = 0 ; i2 < N1 ; i2++)
                                            cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1, inP+i+i2*N1,N1*N1, 0.,laterP+i+i2*N1,N1*N1  );
                                    break;
                                default:
                                    break;
                            }
                        }
                        else if ( bd == four ){
                            inta i,i2,i3;
                    
                            
                            
                            switch (f1.name[su].space[space].block){
                                case tv1:
                                  //  printf("su%d %f %d\n",su, cblas_dnrm2(N1*N1, suP, 1),N1);

                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( i2 = 0 ; i2 < N1 ; i2++)
                                    for ( i3 = 0 ; i3 < N1 ; i3++){
                                                cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1,
                                                            inP +i*N1 +i2*N1*N1 +i3*N1*N1*N1,1, 1.,laterP +i*N1 +i2*N1*N1 +i3*N1*N1*N1, 1  );
                                      //  printf("later%d %f %d\n",su, cblas_dnrm2(N1, laterP+i*N1 +i2*N1*N1 +i3*N1*N1*N1, 1),N1);
                                    }
                                 //   printf("in%d %f %d\n",su, cblas_dnrm2(N2, inP, 1),N1);

                                    break;
                                case tv2:
                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( i2 = 0 ; i2 < N1 ; i2++)
                                            for ( i3 = 0 ; i3 < N1 ; i3++)
                                                cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1,
                                                            inP +i +i2*N1*N1 +i3*N1*N1*N1,N1, 1.,laterP +i +i2*N1*N1 +i3*N1*N1*N1, N1  );
                                    break;
                                case tv3:
                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( i2 = 0 ; i2 < N1 ; i2++)
                                            for ( i3 = 0 ; i3 < N1 ; i3++)
                                                cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1,
                                                            inP +i +i2*N1 +N1*N1*N1*i3,N1*N1, 1.,laterP +i +i2*N1 +i3*N1*N1*N1, N1*N1  );
                                    break;
                                case tv4:
                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( i2 = 0 ; i2 < N1 ; i2++)
                                            for ( i3 = 0 ; i3 < N1 ; i3++)
                                                cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,suP, N1,
                                                            inP +i +i2*N1 +i3*N1*N1,N1*N1*N1, 1.,laterP +i+ i2*N1+ i3*N1*N1, N1*N1*N1  );
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                }
                   else  if ( Bodies ( f1,su,space) == two ){
                       if ( species(f1,su) == eikonDiagonal ){
                            flow *= 1;
                            {
                               ///intended for exchange multiply
                               diagonalOp(bd,f1.name[su].space[space].act,e12, f1.name[su].space[space].block,N1,inP, suP,laterP);
                            }
                        }else if ( species(f1,su) == eikonSemiDiagonal ){
                            
                                ///intended for exchange multiply

                                if ( timer == 0 ){
                                    flow *= -1;

                                    topezOp(0, 1, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP, 1,midP);
                                    diagonalOp(bd,f1.name[su].space[space].act,e12, f1.name[su].space[space].block,N1,midP, suP,laterP);

                               }

                                else if ( timer == 1 ){
                                    flow *= 1 ;
                                    topezOp(0,1, bd,f1.name[su].space[space].act,tv2, f1.name[su].space[space].block,N1,inP, 1,midP);
                                    diagonalOp(bd,f1.name[su].space[space].act,e12, f1.name[su].space[space].block,N1,midP, suP,laterP);
                                    
                                } else
                                 if ( timer == 2 ){
                                    flow *= 1;
                                     
                                     diagonalOp(bd,f1.name[su].space[space].act,e12, f1.name[su].space[space].block,N1,inP, suP,midP);
                                     topezOp(0,1, bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,midP, 1,laterP);

                                     
                                     
                                }
                            
                                 else if ( timer == 3 ){
                                     flow *= -1 ;
                                     
                                     diagonalOp(bd,f1.name[su].space[space].act,e12, f1.name[su].space[space].block,N1,inP, suP,midP);
                                     topezOp(0,1, bd,f1.name[su].space[space].act,tv2, f1.name[su].space[space].block,N1,midP, 1,laterP);

                                 }
                            
                        }else
                        if ( species(f1,su) == eikonOffDiagonal ) {
                            
                            if ( f1.name[su].space[space].block == d12 )
                            {
                                ///intended for direct multiply
                                ///bd = two; act = 1; block = d12  ##ALL DIRECTS## SAME way

                                if ( timer == 0 ){
                                flow *= 1;
                                //A
                                    cblas_dgemv(CblasColMajor, CblasNoTrans, N1, N1, 1.,suP,N1,inP,N1+1, 0.,laterP,N1+1);
                                    topezOp(0,1., bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,laterP, 1,midP);
                                    topezOp(0,1., bd,f1.name[su].space[space].act,tv2, f1.name[su].space[space].block,N1,midP, 1,laterP);

                                }
                            
                                else if ( timer == 1 ){
                                    flow *= 1 ;
                                    topezOp(0,1., bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP, 1,laterP);
                                    topezOp(0,1., bd,f1.name[su].space[space].act,tv2, f1.name[su].space[space].block,N1,laterP, 1,midP);
                                    cblas_dgemv(CblasColMajor, CblasNoTrans, N1, N1, 1.,suP,N1,midP,N1+1, 0.,laterP,N1+1);

                                }
                                else if ( timer == 2){
                                    flow *= -1 ;
                                    topezOp(0,1., bd,f1.name[su].space[space].act       ,tv1, f1.name[su].space[space].block,N1,inP, 1,laterP);
                                    cblas_dgemv(CblasColMajor, CblasNoTrans, N1, N1, 1.,suP,N1,laterP,N1+1, 0.,midP,N1+1);
                                    topezOp(0,1., bd         ,f1.name[su].space[space].act,tv2, f1.name[su].space[space].block,N1,midP, 1,laterP);

                                } else if (timer ==3 ){
                                    flow *= -1;
                                    topezOp(0,1., bd,f1.name[su].space[space].act       ,tv2, f1.name[su].space[space].block,N1                                  ,inP,   1,laterP);
                                    cblas_dgemv(CblasColMajor, CblasNoTrans, N1, N1, 1.,suP,N1,laterP,N1+1, 0.,midP,N1+1);
                                    topezOp(0,1., bd         ,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1  ,midP,  1,laterP);

                                }

                                
                                
                            }else
                                {
                                ///intended for exchange multiply

                                if ( timer == 0 ){
                                flow *= 1;
                                //A
                                    diagonalOp(bd,f1.name[su].space[space].act,e12, f1.name[su].space[space].block,N1,inP, suP,laterP);
                                    topezOp(0,1., bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,laterP, 1,midP);
                                    topezOp(0,1., bd,f1.name[su].space[space].act,tv2, f1.name[su].space[space].block,N1,midP, 1,laterP);

                                }
                            
                                else if ( timer == 1 ){
                                    flow *= 1 ;
                                    topezOp(0,1., bd,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1,inP, 1,laterP);
                                    topezOp(0,1., bd,f1.name[su].space[space].act,tv2, f1.name[su].space[space].block,N1,laterP, 1,midP);
                                    diagonalOp(bd,f1.name[su].space[space].act,e12, f1.name[su].space[space].block,N1,midP, suP,laterP);

                                }
                                else if ( timer == 2){
                                    flow *= -1 ;
                                    topezOp(0,1., bd,f1.name[su].space[space].act       ,tv1, f1.name[su].space[space].block,N1,inP, 1,laterP);
                                    diagonalOp(bd,f1.name[su].space[space].act    ,e12, f1.name[su].space[space].block,N1,laterP, suP,midP);
                                    topezOp(0,1., bd         ,f1.name[su].space[space].act,tv2, f1.name[su].space[space].block,N1,midP, 1,laterP);

                                } else if (timer ==3 ){
                                    flow *= -1;
                                    topezOp(0,1., bd,f1.name[su].space[space].act       ,tv2, f1.name[su].space[space].block,N1                                  ,inP,   1,laterP);
                                    diagonalOp(bd,f1.name[su].space[space].act    ,e12, f1.name[su].space[space].block,N1                                  ,laterP, suP,midP);
                                    topezOp(0,1., bd         ,f1.name[su].space[space].act,tv1, f1.name[su].space[space].block,N1  ,midP,  1,laterP);

                                    }
                                }
                            }
                        }
                        cblas_daxpy(N2, flow, laterP, 1, bufferP, 1);
                        ///sum collect to outP.
                    }
                            if ( f1.name[su].multNext != nullName){
                            printf("multi\n");
                su =f1.name[su].multNext;
                                inP = streams(f1, inT, inS,space)+inR*N2;
                                cblas_dcopy(N2, bufferP,1, inP,1);
                                ///link ot last multiply only
                            }else{
                                su = f1.name[su0].loopNext;//sum channel
                                su0 = su;
                                cblas_daxpy(N2, 1.0, bufferP, 1, outP, 1);
                                inP = initP;///link to begin
                            }
                        
                        }

                    }
                
            
        }
    else{ inta space;
        for ( space = 0 ; space < SPACE ; space++)
        if ( f1.canon[space].body != nada ){
                    if ( species(f1,left) == matrix && species(f1,right) == vector){
                        if ( bodies(f1,left) == bodies(f1,right))
                {
                    inta N1 = vectorLen(f1, space);
                    
                        cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,scalar,
                                streams(f1,left,lspin,space)+l*N1*N1, N1,
                                streams(f1,right,rspin,space)+r*N1,1, 0.,
                                streams(f1,equals,espin,space)+e*N1, 1  );
                    scalar = 1.0;

                }else if ( bodies(f1,left) < bodies(f1,right))
                {
                    inta N1 = outerVectorLen(f1,bodies(f1,left),space);
                    inta N2 = vectorLen(f1, space);
                        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,N1,N2/N1,N1,scalar,streams(f1, left, lspin,space)+l*N1*N1,N1,streams(f1,right,rspin,space)+r*N2,N1,0.,streams(f1, equals, espin,space)+e*N2,N1);
                    scalar = 1.0;

                }
            
            }
        
    else if ( species(f1, left) == diagonalMatrix && species(f1,right) == vector && species(f1,equals) == vector){
    ///for diagonal potential matrices, where left points at a vector
    
        inta N2 = vectorLen(f1, space);
        cblas_dcopy(N2, streams(f1, right,rspin,space)+r*N2, 1, streams( f1, equals, espin,space)+e*N2, 1);
        cblas_dscal(N2, scalar, streams( f1, equals, espin,space)+e*N2, 1);
        cblas_dtbmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, N2, 0,streams( f1, left, lspin,space )+l*N2, 1,streams( f1, equals, espin,space)+e*N2,1);
        scalar = 1.0;
    }
    else if ( species(f1, left) == vector && species(f1,right) == vector && species(f1,equals) == vector){
    ///specifically for geminals * transpose-geminals = geminal
    //if ( bodies(f1, left) == two && bodies(f1,in) == two && bodies(f1,out)== two )
    
        inta N1 = vector1Len(f1, space);
        inta N2 = N1*N1;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,N1,N1,N1,scalar,streams( f1, left, lspin,space )+l*N2,N1,streams(f1, right,rspin,space)+r*N2,N1, 0.,streams( f1, equals, espin,space)+e*N2,N1);
        scalar = 1.0;

    }
        }
    }
    return 0;
}

//inta tGEVV (inta rank,    sinc_label  f1,  inta space,  division equals,inta e,  inta espin,   division left,inta l,inta lspin,    division right,inta r, inta rspin ){
//    if ( header(f1, left ) != header(f1, right ) ){
//        printf("Two Head types GEVVo\n");
//        exit(0);
//    }
//    inta nl,nr,nm,space2 = space;
//      genusType sl = species(f1,left);
//      genusType sr = species(f1,right);
//    inta cr = CanonicalRank(f1, left, lspin);
//
//    if (( sl == outerVector ) && ( sr == vector ) ){
//          division inT,midT,outT;
//        inta inR,midR,outR,inS,midS,outS;
//        inT = right;
//        inR = r;
//        inS = rspin;
//
//        outT = equals;
//        outR = e;
//        outS = espin;
//        if ( in != 1 ){
//           tPermuteOne(rank, f1, space, in, right, r, rspin, canonicalvvVector,0, rank);
//           inT = canonicalvvVector;
//           inR = 0;
//           inS = rank;
//        }
//        midT = canonicalvv2Vector;
//        midR = 0;
//        midS = rank;
//
//
//        if (out != 1 ){
//            outT = canonicalvv3Vector;
//            outR = 0;
//            outS = rank;
//        }
//        if (sl == outerVector)
//            nl = outerVectorLen(f1, bodies(f1, left), space);
//        else
//            nl = vectorLen(f1, space);
//
//        if (sr == outerVector)
//            nr = outerVectorLen(f1, bodies(f1, right), space);
//        else
//            nr = vectorLen(f1, space);
//
//        if ( nl == nr )
//        {
//            cblas_dcopy(nl,streams(f1,inT,inS,space)+inR*nl ,1,streams(f1,outT,outS,space)+outR*nl ,1  );
//            cblas_dscal(nl,tDOT(rank, f1, space, 1, midT, midR,midS,1, outT, outR, outS),streams(f1,outT,outS,space)+outR*nl,1);
//        }
//            else if (nl < nr ) {
//                inta l1 = l%cr, l2 = (l/cr)% cr;
//                printf("%d**%d\n",l1,l2);
//                        //dgemv on right * left -> reduced dimensional intermediate
//                        nm = nr/nl;
//                cblas_dgemv( CblasColMajor, CblasNoTrans,  nm,  nl,1.,
//                        streams( f1, inT, inS,space )+inR*nr, nm,
//                        streams(f1, left, lspin,space)+l1*nl,1, 0.,
//                        streams( f1, midT, midS,space )+midR*nm, 1  );
//                cblas_dger(CblasColMajor, nl,nm, 1. , streams(f1,left,lspin,space) + l2 * nl,1 , streams(f1, midT,midS,space)+midR*nm,1, streams(f1, outT,outS,space)+outR*nr,nl);
//            }else {
//                printf("dimensional error\n");
//                exit(0);
//            }
//
//        if (out != 1 ){
//            tPermuteOne(rank, f1, space2, out, outT,outR,outS, equals, e,espin);
//        }
//
//    }
//    return 0;
//}

/**
 *general vector - vector inner product per dimension
 *
 *@param rank      OMP rank
 *@param f1          container
 *@param leftChar  left input group action  defalult 1
 *@param left  vector
 *@param l  one of the matrices canonical ranks indexed
 *@param lspin matrix spin
 *@param rightChar right input group action
 *@param right input division
 *@param r one of the input canonical ranks indexed
 *@param rspin input spin
*/
double tDOT (inta rank,    sinc_label  f1,inta dim,char leftChar,   division left,inta l,inta lspin, char rightChar,   division right,inta r, inta rspin ){
    inta space = dim,brab,bspin ,ketk, kspin;
    double prod = 0.;
      division bra,ket;
    f1.name[canonicaldotVector].Current[rank] = 0;
    f1.name[canonicaldot2Vector].Current[rank] = 0;
    if ( rank ){
        if ( rightChar != CDT || leftChar != CDT )
        if ( ! allowQ(f1.rt, blockParallelPermuteblock)){
            printf("blockParallelPermuteblock allow\n");
            fflush(stdout);
            exit(0);
        }
    }
    if ( rightChar != CDT || leftChar != CDT )
        if ( ! allowQ(f1.rt, blockPermutationsblock)){
            printf("blockPermutationsblock allow\n");
            fflush(stdout);
            exit(0);
        }

    
    if ( rightChar != CDT){
        tPermuteOne(rank, f1, space, rightChar, right, r, rspin, canonicaldotVector,0, rank);
        bra = canonicaldotVector;
        brab = 0;
        bspin = rank;
    }else {
        bra = right;
        brab = r;
        bspin = rspin;
    }
    if (leftChar != CDT ){
        tPermuteOne(rank, f1, space, leftChar, left, l, lspin, canonicaldot2Vector,0, rank);
        ket = canonicaldot2Vector;
        ketk = 0;
        kspin = rank;
        
    }else{
        ket = left;
        ketk = l;
        kspin = lspin;
    }
    inta al = alloc(f1, right, space );
    if ( alloc(f1, left, space ) == alloc(f1, right, space )){
        inta N1 = al;
        prod = cblas_ddot( N1 , streams(f1,bra,bspin,space)+brab*N1,1 , streams(f1,ket,kspin,space)+ketk*N1, 1);
    } else {
        printf("body count %d %d %d %d\n",left,alloc(f1, left, space ) , right,alloc(f1, right, space ) );
        exit(1);
    }
    return prod;
}


/**
 *Multiply by dimension
 *
 *@param rank      OMP rank
 *@param f1          container
 *@param left matrix
 *@param l matrix index, across all chained elements
 *@param im matrix spin
 *@param prod scalar multiply
 *@param[in] ket vector
 *@param k ket index
 *@param sp2 ket spin
 *@param[out] oket vector
 *@param o oket index
 *@param ospin spin
 */
inta tHX(  inta rank,   sinc_label f1 ,division left, inta l, inta im, double prod,   division ket , inta k, inta sp2,   division oket, inta o,inta ospin ){
    if ( left == nullName || species(f1, left ) == scalar){
        printf("*");
        return 0;
    }
    division in=nullName,out=nullName;
    inta inSp,outSp,inRank,outRank;
    
    inta lll;

        
        if ( rank ){
            ///check for parallel allocations
            if ( ! allowQ(f1.rt, blockParallelMultiplyblock)){
                printf("blockParallelMultiplyblock allow\n");
                fflush(stdout);
                exit(0);
            }
        }
        division ll = name(f1,left);
        inta mi = 0,xi=0;
            while ( ll != nullName){
                xi += CanonicalRank(f1, ll, im);
                lll =  l-mi;
                
                if ( mi <= l && l < xi ){
                        out =  oket;
                        outRank = o;
                        outSp = ospin;
                        
                        in = ket;
                        inRank = k;
                        inSp= sp2;
                        tGEMV(rank, f1,out,outRank,outSp,prod, ll, lll, im,in, inRank,inSp);
                        break;
                    }
                mi += CanonicalRank(f1, ll, im);
                ll = f1.name[name(f1,ll)].chainNext;
            }
    
    return 1;
}

/**
 *Multiply in total
 *
 *@param f1          container
 *@param[out] bra vector
 *@param left matrix
 *@param shiftFlag  for product * Hv +  sum v
 *@param[in] right vector
 *@param tolerance a number setting the absolute quality
 *@param relativeTolerance a number seting quality relative to magnitude of origin
 *@param condition Beylkin's condition (alpha)
 *@param threshold the smallest number
 *@param maxCycle the maxmium number of cycles in this routine
 *@param canon rank of output vector
*/
void tHXpY ( sinc_label f1 , division bra, division left,inta shiftFlag, division right , double tolerance , double relativeTolerance, double condition, double threshold, inta maxCycle, double maxCondition, inta canon, inta X1){
    double prod;
    inta rank0 = 0 ,rank;
    mea co2,coi;
    inta ilr,Ll,sp2,Rr,im,l , k,targSpin;
    
    if ( ! allowQ(f1.rt,blockTotalVectorBlock)){
        printf("blockTotalVectorBlock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    if (  right == totalVector){
        printf("you cannot feed totalVector into tHXpY\n");
        exit(0);
    }
    for ( targSpin = 0 ; targSpin < spins(f1, right ) ;targSpin++){
        zero(f1, totalVector, rank0);
        f1.name[totalVector].Current[rank0] =0;

        if (  bra != totalVector){
            if ( shiftFlag  == 1){
                tEqua(f1, totalVector, rank0, bra, targSpin);
            }
            else if ( shiftFlag  == -1){
                tEqua(f1, totalVector, rank0, right, targSpin);
                mea ME;
                ME = pMatrixElement(f1, right, targSpin, left, 0, right, targSpin);
                tScaleOne(f1, totalVector, rank0, -ME);
            }
            else
                f1.name[totalVector].Current[rank0] = 0;
        }
         {
            for ( im = 0; im < spins(f1, left ); im++)
                for ( sp2 = 0; sp2 < spins(f1,right); sp2++)
                {
                    if (sp2 == 1)
                        co2 = I;
                    else
                        co2 = 1;
                    if (im == 1 )
                        coi = I;
                    else
                        coi = 1;
                    if ( targSpin == 0 )
                        prod  = creal(co2 * coi );
                    else
                        prod = cimag(co2 * coi ) ;
                    if ( fabs(prod) > f1.rt->THRESHOLD ){
                        Rr = CanonicalRank ( f1, right,sp2 );
                        Ll = CanonicalOperator(f1, left, im);
                        inta su = f1.name[totalVector].Current[rank0];
                        
                        #ifdef OMP
                        #pragma omp parallel for private (ilr,l,k,rank) schedule(dynamic,1)
                        #endif
                            for ( ilr = 0; ilr <  Ll*Rr ; ilr++)
                            {

                            #ifdef OMP
                                    rank = omp_get_thread_num();
                            #else
                                    rank = 0;
                            #endif

                                l = ilr%Ll;
                                k = ilr/Ll;
                                tHX(rank, f1,left, l, im,prod, right, k, sp2,totalVector,ilr +su, rank0);
                            }
#if VERBOSE
                        printf("'lambda' -> %d %d %d\n", su , su+Ll*Rr,part(f1, totalVector ));
                        
#endif
                        f1.name[totalVector].Current[rank0]+= Ll*Rr;
                        if (f1.name[totalVector].Current[rank0] > part(f1, totalVector ) )
                        {
                            printf("'lambda' is too small\n");
                            fflush(stdout);
                            exit(1);
                        }
                    }
                }
        };
        if (  bra != totalVector){
            CanonicalRankDecomposition( f1,  NULL,totalVector, rank0, bra, targSpin, tolerance,relativeTolerance, condition,threshold,maxCycle,maxCondition, canon,0);
            if ( X1 > 0 && canon > 1 ){
                tEqua(f1, totalVector, rank0, bra, targSpin);
                CanonicalRankDecomposition( f1,  NULL,totalVector, rank0, bra, targSpin, tolerance,relativeTolerance, condition,threshold,maxCycle,maxCondition, canon,X1);
            }
        }
    }

    return;
}


/**
 *matrix element
 *upgraded in v9.3 for parallel operation
 *only one of these can run at a time.
 *
 *@param f1 container
 *@param alloy1 vector
 *@param spin1 vector's spin
 *@param op operator
 *@param ospin operator spin
 *@param alloy2 vector
 *@param spin2 vector's spin
 */
double pMatrixElement ( sinc_label  f1 ,   division alloy1 , inta spin1, division op, inta ospin, division alloy2 , inta spin2 ){
#ifndef CHERRY_PICKER
    return tMatrixElements(0, f1, alloy1, spin1, op, ospin, alloy2, spin2);
#else

    inta rank,i ,cl = CanonicalRank(f1, alloy1, spin1),cr = CanonicalRank(f1, alloy2, spin2);
    mea OV, ov[f1.rt->NLanes];
    
    
    
    
    division left=0;
    for ( rank = 0; rank < f1.rt->NLanes; rank++){
        ov[rank] = 0.;
        if ( ! rank )
            left = anotherLabel(&f1, 0, nada);
        else
            anotherLabel(&f1, 0, nada);

        f1.name[left+rank].Current[spin1] = 1;
        f1.name[left+rank].name = alloy1;
    }
    ///need to be in order to make each rank-array contiguous.
    division right=0 ;
    for ( rank = 0; rank < f1.rt->NLanes; rank++){
        if ( ! rank )
            right = anotherLabel(&f1, 0, nada);
        else
            anotherLabel(&f1, 0, nada);
        f1.name[right+rank].Current[spin2] = 1;
        f1.name[right+rank].name = alloy2;
    }

    
    #ifdef OMP
    #pragma omp parallel for private (i,rank) schedule(dynamic,1)
    #endif
                for ( i = 0 ;i < cl*cr; i++)
                {
    #ifdef OMP
                    rank = omp_get_thread_num();
    #else
                    rank = 0;
    #endif
                    f1.name[left+rank].Begin[spin1] = (i)%cl+f1.name[alloy1].Begin[spin1];
                    f1.name[right+rank].Begin[spin2] = (i/cl)+f1.name[alloy2].Begin[spin2];
                    ov[rank] += (tMatrixElements(rank, f1, left+rank, spin1, op, ospin, right+rank, spin2));
                }
    OV = 0.;
    for ( rank = 0; rank < f1.rt->NLanes; rank++)
        OV += ov[rank];

#ifdef COMPLEXME
    return (creal(OV));
#else
    return (OV);
#endif
#endif
}

/**
 *outer product of two vectors
 *
 *for building operators from vectors
 *
 *@param f1 container
 *@param space   in canon which has SPACE components
 *@param[in] vector any content
 *@param a vector1 spin
 *@param[in] vector2 any content
 *@param b vector's spin
 *@param[out] proj outer product matrix
 *@param c proj spin
*/
inta tOuterProductSuOne(   sinc_label  f1,inta space,  division vector , inta a,   division vector2,inta b,   division proj, inta c){
    inta ma = CanonicalRank(f1, vector,a), mb = CanonicalRank(f1, vector,b), zc = CanonicalRank(f1, proj, c),l,r,i;
    inta Np = alloc(f1, proj, space), n1 = alloc(f1, vector,space),n2 = alloc(f1, vector2,space);

    if ( species(f1, vector ) == matrix || species(f1, vector2) == matrix){
        printf("input arguments to Outer product need to be vectors\n");
        exit(0);
    }
    if ( zc + ma*mb > part(f1, proj) && proj < eigenVectors){
        printf("outerProductSu:  \n");
        exit(0);
    }
    if ( Np != n1 * n2 ){
        printf("suspect...");
        exit(1);
    }
    
    for ( l = 0 ; l < ma; l++)
        for ( r = 0; r < mb; r++)
                for ( i = 0; i < Np ; i++)
                    (streams(f1, proj,c,space)+(l*mb+r+zc)*Np)[i] = 0.;
    
    for ( l = 0 ; l < ma; l++)
        for ( r = 0; r < mb; r++)
                cblas_dger(CblasColMajor, n1,n2, 1. , streams(f1, vector,a,space)+l*n1,1, streams(f1, vector2,b,space)+r*n2,1, streams(f1, proj,c,space)+(l*mb+r+zc)*Np,n1);
    f1.name[proj].Current[c] += ma*mb;
    return 0;
}

//inta compressReplaceEikon(  sinc_label f1 ,   division eik ){
//    if ( species(f1, eik) != eikon){
//        return 1;
//    }
//    inta i = 0,ii;
//
//      division looper,chainer;
//      genus hidden;
//
//    for ( hidden = eikonDiagonal ; hidden <= eikonSemiDiagonal ; hidden++){
//        i = 0;
//        tClear(f1,copyVector);
//        tClear(f1,copyTwoVector);
//
//        for ( chainer = eik; chainer != nullName ; chainer= f1.name[chainer].chainNext){
//
//        for ( looper = chainer; looper != nullName ; looper= f1.name[looper].loopNext)
//            if ( f1.name[looper].species == hidden ){
//
//            tAddTw(f1, copyVector,  0, looper, 0);
//            i++;
//        }
//    }
//
//        for ( ii = 1 ; ii < 30 ; ii++){
//            tId(f1,copyTwoVector,0);
//            CanonicalRankDecomposition(0, f1, NULL, copyVector, 0, copyTwoVector, 0,  f1.rt->TOLERANCE,  f1.rt->relativeTOLERANCE,  f1.rt->ALPHA, f1.rt->THRESHOLD, f1.rt->MAX_CYCLE, <#inta canon#>)
//        }
//    }
//
//    return 0;
//}

void analyzeLoopElement( sinc_label f1, division loopElement, inta spin ){
#ifdef PRINT_LOOP_STRUCTURE
    division multer;
    for (multer = loopElement ; multer != nullName ; multer = f1.name[multer].multNext){
        printf("++ %d ", multer);
    }
    printf("\n");
#endif
    return;
}

void analyzeChainElement( sinc_label f1, division chainElement, inta spin ){
#ifdef PRINT_CHAIN_STRUCTURE
    if ( chainElement != nullName){

        printf("%d (%d)\n", chainElement, CanonicalRank(f1, chainElement , spin));
        division looper;
        for (looper = f1.name[chainElement].loopNext ; looper != nullName ; looper = f1.name[looper].loopNext){
            analyzeLoopElement( f1, looper,spin);
        }
        printf("\n");
        analyzeChainElement(f1, f1.name[chainElement].chainNext, spin);
    }
#endif
    return;

}
