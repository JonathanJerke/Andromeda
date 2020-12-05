/**
 *  coreUtil.c
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
                f1->name[output].Current[0] = 0;
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
                f1->name[output].Current[0] = 0;
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
    inta n,m,m2;
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
        default:
            break;

    }


    
    switch ( bd ){
        case one:
            //  vector -> vectorOut
            sign = 1.;
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
                        mult = 1./n/lattice/lattice;
                    }else {
                        mult = 2*1./(n*n)/lattice/lattice;
                    }
                    cblas_daxpy(N1-n, sign*mult, vector+n, 1, vectorOut, 1);

                    
                    cblas_daxpy(N1-n, sign1*sign*mult, vector, 1, vectorOut+n,1);
                    sign *= -1;

                }
            }
            break;
        case two:
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
            
            
            if (pw == 2 || pw == 1) {
                for ( m = 0; m < N1 ; m++)
                {
                    sign = 1.;
                    mult = 0;
                    for (n = 1 ; n < N1 ; n++){
                        if (pw == 1 ){
                            mult = 1./n/lattice/lattice;
                        }else if ( pw == 2 ){
                            mult = 2*1./(n*n)/lattice/lattice;
                        }
                        cblas_daxpy(N1-n, sign*mult, vector+n1[perm[op[0]]]*n+n1[perm[op[1]]]*m ,n1[perm[op[0]]], vectorOut+n1[op[1]]*m ,n1[op[0]]);
                        cblas_daxpy(N1-n, sign1*sign*mult, vector+n1[perm[op[1]]]*m , n1[perm[op[0]]], vectorOut+n1[op[0]]*n+n1[op[1]]*m ,n1[op[0]]);
                        sign *= -1;
                    }
                }
            }
            
            
            break;
        case three:
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
            
            if (pw == 2 || pw == 1) {
                for ( m = 0; m < N1 ; m++)
                    for ( m2 = 0; m2 < N1 ; m2++)
                    {
                        sign = 1.;
                        mult = 0;

                            for (n = 1 ; n < N1 ; n++){
                                       if (pw == 1 ){
                                           mult = 1./n/lattice/lattice;
                                       }else if ( pw == 2 ){
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
    inta perm[7],m,m2,n;
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
            if ( f1.name[label].species == eikon){
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
    inta found ;
      division ll = f1.name[name(f1,label)].chainNext,zz;
    while ( ll != nullName ){
        zz = label;
        found = 0;
        while ( zz != ll ){
            if (f1.name[zz].multId == f1.name[ll].multId)
                found = 1;
            zz = f1.name[name(f1,zz)].chainNext;
        }
        if ( ! found ){
            rr += CanonicalRank(f1, name(f1,ll), spin);//switch from product to addition!!!
        }
        
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

inta tBoot (   sinc_label f1 ,   division label,inta spin ){
    
    inta I1,I2,I3,space;
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
            if ( f1.name[label].species == vector && f1.canon[space].body == one){
            {
                            floata  * stream = streams(f1,label,spin,space)+Current*B1[space];
                                for ( I2 = 0 ; I2 < B1[space] ; I2++){
                                    stream[I2] = exp(-(I2-(B1[space]-1)/2)*(I2-(B1[space]-1)/2)*0.5);
                                }
                        }
                    }else
            if ( f1.name[label].species == vector && f1.canon[space].body == two){
                    {
                
                floata  * stream = streams(f1,label,spin,space)+Current*B1[space];
                for ( I1 = 0 ; I1< B1[space] ; I1++)
                    for ( I2 = 0 ; I2 < B1[space] ; I2++){
                        stream[I1*B1[space]+I2] = exp(-(I1-(B1[space]-1)/2)*(I1-(B1[space]-1)/2)*0.5)*exp(-(I2-(B1[space]-1)/2)*(I2-(B1[space]-1)/2)*0.5);
//              if ( I1 < I2 )
//                  stream[I1*B1[space]+I2]  *= -1;
//                if ( I1 == I2 )
//                    stream[I1*B1[space]+I2]  *= 0;
//
                    }
                
            }
        }
        else
            if ( f1.name[label].species == vector && f1.canon[space].body == three){
                    {
                
                floata  * stream = streams(f1,label,spin,space)+Current*B1[space];
                for ( I1 = 0 ; I1< B1[space] ; I1++)
                    for ( I2 = 0 ; I2 < B1[space] ; I2++)
                        for ( I3 = 0 ; I3 < B1[space] ; I3++){
                            stream[I3*B1[space]*B1[space]+I1*B1[space]+I2] = exp(-(I1-(B1[space]-1)/2)*(I1-(B1[space]-1)/2)*0.5)*exp(-(I2-(B1[space]-1)/2)*(I2-(B1[space]-1)/2)*0.5)*exp(-(I3-(B1[space]-1)/2)*(I3-(B1[space]-1)/2)*0.5);
                }
            }
        }


      
        
    }
    return 0;
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




int double_cmp( const void * a ,const void * b ){
    const double* A = (const double* )a;
    const double* B = (const double* )b;
    if ( *A < * B )
        return 1;
    else if ( *A == *B )
        return 0;
    else return -1;
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
                    buildConstant(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink, c->i.terms[i].label, 0, real);
                    break;
                case 2:
                    buildLinear(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label, 0, real);
                    break;
                case 3:
                    buildSpring(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label, 0, real);
                    break;
                case 4:
                    buildDeriv(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink, c->i.terms[i].label,0, real);
                    break;
                case 5:
                    buildKinetic(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label, 0, real);
                    break;
                case 6:
//                    buildClampKinetic(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label, 0, real);
                    break;
                case 7:
                    buildElement(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, prevLink,  c->i.terms[i].label, 0, real,c->i.terms[i].bra,c->i.terms[i].ket);
                    break;
                case 8:
                    buildExternalPotential(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl, c->i.terms[i].adjustOne,prevLink,  c->i.terms[i].label,c->i.terms[i].embed, 0, real,c->i.terms[i].mu,c->i.terms[i].atom);
                    break;
                case 9:
                    buildPairWisePotential(c, f1, c->i.terms[i].scalar,c->i.terms[i].invert,c->i.terms[i].act, c->i.terms[i].bl,c->i.terms[i].adjustOne,prevLink,  c->i.terms[i].label, c->i.terms[i].embed,0, real,c->i.terms[i].mu);
                    break;
                case 10:
                    assignDiagonalMatrix(c,f,c->i.terms[i].filename,prevLink);
                    break;
            }
        }
    }
    return term;
}

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
    x.component = f1.canon[space].space +1;
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
                //                //band_ii' Vector_ijkl = Vector_i'jkl

                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,N2,L1*L1*L1,L1,1. ,band         ,L1,streams( f1, vector1, s1,space )+r*L1*L1*L1*L1,L1, 0.,  buffer  , N2   );
                //factor third-fourth particles,
                // Vector_i'j--ml band_jj'=> Vector_i'j'--ml
                for ( l = 0 ; l < L1 ;l++)
                    for ( k = 0 ; k < L1 ;k++)
                        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans  ,N2,N2,L1,1.,buffer+L1*N2*k+L1*L1*N2*l      ,N2,band+1*L1*N2 ,L1, 0.,  buffer2 + N2*N2*k+L1*N2*N2*l , N2   );
                
                //factor fourth particles,
                // Vector_i'j'm--l band_k'k=> Vector_i'j'k'--l
                for ( l = 0 ; l < L1 ;l++)
                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans  ,N2*N2,N2,L1,1.,buffer2+N2*N2*L1*l      ,N2,band+2*L1*N2 ,L1, 0.,  buffer + N2*N2*N2*l , N2   );
                
                //Vector_i'j'k'l band_ll' = Vector_i'j'k'l'
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans  ,N2*N2*N2,N2,L1,1. ,buffer                  ,L1,band+3*L1*N2 ,L1, 0.,  streams( f2, out, s2,space)+r*N2*N2*N2*N2 , N2   );
                                
            }
        }
    f2.name[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    return 0.;
}


