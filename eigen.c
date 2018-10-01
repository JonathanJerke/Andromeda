/*
 *  eigen.c
 *
 *
 *  Copyright 2018 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University
 *  and the Robert A. Welch Foundation.
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


#include "eigen.h"


INT_TYPE tNBodyConstruction (struct calculation * c1, enum division build ,enum division eigen){
    struct field * f1 = &(c1->i.c);
    enum body bootBodies = c1->i.body;
    INT_TYPE cmpl,space,matrixNumber = c1->i.decomposeRankMatrix;
    //THRESING FLOOR
    tClear(f1, eigen);
    if ( c1->rt.printFlag  ){
        tScale(f1, density,-1);
//        tEqua(f1, squareTwo, 0, density, 0);
//        tPermute(0, f1, 'b', squareTwo, 0, density, 0);
        if (bodies ( f1,eigen ) == two)
            tCycleDecompostionOneMP(0, f1, density, 0, eigen, 0, f1->mem1->rt->CONVERGENCE, part(f1, eigen), -1);
    }
    
    else {
        
        for ( cmpl =  0; cmpl < 2 ; cmpl++ ) {
            tClear(f1, build);
            zero ( f1, eigen,cmpl);
            zero(f1,build,0);
            tClear(f1, copy);
            if ( cmpl == 0 ){
                tEqua(f1, copy ,0, kinetic,0);
                if ( c1->i.springFlag ){
                    tAddTw(f1,copy,0, harmonium ,0);
                }
                
                tSumMatrices(f1, build, copy);// B x (1+S) * 3
             //   printf("build %f\n", traceOne(f1, build, 0));
                tCycleDecompostionOneMP(0, f1, build, 0, eigen, 0, f1->mem1->rt->CONVERGENCE, part(f1, eigen), -1);
                
                if ( bootBodies > one ){
#if VERBOSE
                    printf("ie %d\n", CanonicalRank(f1, interactionExchange, 0));
#endif
                    
                    if ( CanonicalRank(f1, interactionExchange, 0)){
                        tEqua(f1, squareTwo,0, ocean(0,f1, interactionExchange,0,0),0);
                        tCycleDecompostionOneMP(0, f1, interactionExchange, 0,squareTwo, 0, f1->mem1->rt->CONVERGENCE, part(f1, eigen) , -1);
                        tSumMatrices(f1, build, squareTwo);//B2:1 || B3 : 3
                        tCycleDecompostionOneMP(0, f1, build, 0, eigen, 0, f1->mem1->rt->CONVERGENCE, part(f1, eigen), -1);
                    }
                }
                if ( f1->Na  ){
                    tClear(f1, copy);
                    tCycleDecompostionOneMP(0, f1, linear, 0,copy, 0, f1->mem1->rt->CONVERGENCE, matrixNumber, -1);
                    tSumMatrices(f1, build, copy);
                    tCycleDecompostionOneMP(0, f1, build, 0, eigen, 0, f1->mem1->rt->CONVERGENCE, part(f1, eigen), -1);
                }
            }else {
                if ( CanonicalRank(f1, kinetic,1 ) ){
                    tEqua(f1, copy ,0, kinetic,1);
                    tSumMatrices(f1, build, copy);// B x (1+S) * 3
                    tCycleDecompostionOneMP(0, f1, build, 0, eigen, 1, f1->mem1->rt->CONVERGENCE, part(f1, eigen), -1);
                }
            }
        }
    }
    INT_TYPE r,i;
    INT_TYPE n2[3];
    length(f1, eigen,n2);
    INT_TYPE *n1 = vectorLen(f1, eigen);
    fflush(stdout);
    double ze[SPACE] ,xe[SPACE];
    DCOMPLEX * hmat = (DCOMPLEX*)myStreams(f1, matrixHbuild,0);
    for ( r = 0; r < 1 ; r++){
        for ( space = 0; space < SPACE ; space++){
            for ( i = 0; i < n2[space]; i++)
                hmat[i] = streams(f1, eigen,0,space)[i] + I*streams(f1,eigen,1,space)[i];
            tzheev (0,f1,'V',n1[space],hmat,n1[space],streams(f1,foundationStructure,0,space)+r*n1[space]);
#if VERBOSE
            printf ("space%d eigen %f < %f\n", space,streams(f1,foundationStructure,0,space)[0] ,streams(f1,foundationStructure,0,space)[n1[space]-1] );
#endif
            for ( i = 0; i < n2[space]; i++){
                streams(f1, eigen,0,space)[i] = creal(hmat[i]);
                streams(f1, eigen,1,space)[i] = cimag(hmat[i]);
            }
        }
        
        
        
        ze[0] = 1e254;
        ze[1] = 1e254;
        ze[2] = 1e254;
        xe[0] = -1e254;
        xe[1] = -1e254;
        xe[2] = -1e254;
        
        
        for ( space = 0; space < SPACE ; space++){
            
            if ( ze[space] > streams(f1,foundationStructure,0,space)[r*n1[space]])
                ze[space] = streams(f1,foundationStructure,0,space)[r*n1[space]];
            if ( xe[space] < streams(f1,foundationStructure,0,space)[(r+1)*n1[space]-1] )
                xe[space] = streams(f1,foundationStructure,0,space)[(r+1)*n1[space]-1];
            
            for ( i = 0 ; i < n1[space] ; i++){
                if ( ze[space]*xe[space] > 0 && ze[space] < 0 )
                {
                    // printf("%lld %lld %lld %f \t",space,i,r*n1[space]+i,streams(f1,foundationStructure,0,space)[r*n1[space]+i]);
                    
                    streams(f1,foundationStructure,0,space)[r*n1[space]+i] = exp(-(streams(f1,foundationStructure,0,space)[r*n1[space]+i]-xe[space])/(ze[space]-xe[space]));
                    // printf("%f \n",streams(f1,foundationStructure,0,space)[r*n1[space]+i]);
                    
                    
                }else{
                    
                    streams(f1,foundationStructure,0,space)[r*n1[space]+i] = exp(-(streams(f1,foundationStructure,0,space)[r*n1[space]+i]-ze[space])/(xe[space]-ze[space]));
                    
                }
            }
        }
        
    }
    
    return 0;
}





INT_TYPE tBasisLevel( struct field * f1, enum division A ,INT_TYPE space, double lvl, INT_TYPE ops,enum division build,INT_TYPE xB){
    //GRID
    /// GRID
    ///// GRID
    ////////GRID
    INT_TYPE i,r;
    INT_TYPE classicalBasisSize;
    INT_TYPE n1[3];
    n1[0] = vectorLen(f1,A)[0];
    n1[1] = vectorLen(f1,A)[1];
    n1[2] = vectorLen(f1,A)[2];

    
    double value;
    classicalBasisSize = 0;
    for ( r = 0; r < 1 ; r++)
    {
        for ( i = 0 ; i < n1[space] ; i++)
                {
                    value =streams(f1,foundationStructure,0,space)[r*n1[space]+i];
                    if ( value > lvl ){
                        if( ops == 2 ){
                            
//                            {
//                                cblas_dcopy(n1[space],streams(f1, ocean(0,f1,A,r,0),0,space)+n1[space]*i,1,streams(f1, basis,0,space)+n1[space]*classicalBasisSize,1);
//                            }
                            //printf("%f\n", cblas_ddot(n1[space],streams(f1, basis,0,space)+n1[space]*classicalBasisSize,1,streams(f1, basis,0,space)+n1[space]*classicalBasisSize,1));

                        }
                        classicalBasisSize++;
                        if ( classicalBasisSize >= xB && ops== 2 ){
                            return xB;
                        }
                        
                    }
                }
    }
    if ( ops == -1 ){
        f1->sinc.tulip[build].header -= 1;
        if ( vectorLen(f1, build)[space]  != classicalBasisSize){//for now, this is just reading RDS
            printf("oop parity\n");
            exit(0);
        }
        f1->sinc.tulip[build].header += 1;
    }
    return classicalBasisSize;
}

INT_TYPE tFoundationLevel( struct field * f1, enum division A , double lvlm, double lvlx,INT_TYPE ops,enum division build,INT_TYPE xB, double lvl1, double lvl2, double lvl3,INT_TYPE *mmm){
    //GRID
    /// GRID
    ///// GRID
    ////////GRID
    INT_TYPE i,j,k,r1,r2,r3,space,ii,jj,kk,di;
    INT_TYPE sp, classicalBasisSize;
    INT_TYPE n1[3];
    n1[0] = vectorLen(f1,A)[0];
    n1[1] = vectorLen(f1,A)[1];
    n1[2] = vectorLen(f1,A)[2];
    double value;
    classicalBasisSize = 0;
   // printf("%d %f %f %f\n", ops,lvl,streams(f1,foundationStructure,0,1)[0],streams(f1,foundationStructure,0,1)[100]);

    ii = 0;
    for ( r1 = 0; r1 < 1 ; r1++)
        for ( i = 0 ; i < n1[0] ; i++)
            if ( streams(f1,foundationStructure,0,0)[i+r1*n1[0]] > lvl1 ){
                //printf("%lld\n",ii);
                jj = 0;
                for ( r2 = 0; r2 < 1 ; r2++)
                    for ( j = 0 ; j < n1[1] ; j++){
                        //  printf("%lld %lld %f %f\n", r2,j ,streams(f1,foundationStructure,0,1)[j+r2*n1[1]] , lvl2);
                        if ( streams(f1,foundationStructure,0,1)[j+r2*n1[1]] > lvl2 ){
                            //    printf("%lld %lld\n",ii,jj);
                            
                            kk = 0;
                            for ( r3 = 0; r3 < 1 ; r3++)
                                for ( k = 0 ; k < n1[2] ; k++)
                                    if ( streams(f1,foundationStructure,0,2)[k+r3*n1[2]] > lvl3 )
                                        
                                    {
                                        //             printf("%lld %lld %lld\n",ii,jj,kk);
                                        
                                        value = (streams(f1,foundationStructure,0,0)[i+r1*n1[0]] *
                                                 streams(f1,foundationStructure,0,1)[j+r2*n1[1]] *
                                                 streams(f1,foundationStructure,0,2)[k+r3*n1[2]]);
                                      // if ( ops < 0 )
                                          // printf("%f %d\n", value,ops);
                                        
                                        
                                        if ( lvlm < value && value <= lvlx ){
                                            if( ops == 2 ){
                                                tClear(f1,build+classicalBasisSize);
                                                for ( di = 0; di < part(f1,build+classicalBasisSize ); di++)
                                                {
                                                    f1->sinc.tulip[diagonalVectorA].header = Cube;
                                                    zero(f1, diagonalVectorA,0);
                                                    tClear(f1,diagonalVectorA );
                                                    f1->sinc.tulip[diagonalVectorA].Current[0] = 1;;

                                                    cblas_dcopy(n1[0], streams(f1, A, 0 , 0)+(i+di)*n1[0],1,streams(f1,diagonalVectorA,0,0),1);
                                                    cblas_dcopy(n1[1], streams(f1, A, 0 , 1)+(j+di)*n1[1],1,streams(f1,diagonalVectorA,0,1),1);
                                                    cblas_dcopy(n1[2], streams(f1, A, 0 , 2)+(k+di)*n1[2],1,streams(f1,diagonalVectorA,0,2),1);
                                                    tAddTw(f1, build+classicalBasisSize, 0,diagonalVectorA,0 );

                                                    
                                                    cblas_dcopy(n1[0], streams(f1, A, 1 , 0)+(i+di)*n1[0],1,streams(f1,diagonalVectorA,1,0),1);
                                                    cblas_dcopy(n1[1], streams(f1, A, 1 , 1)+(j+di)*n1[1],1,streams(f1,diagonalVectorA,1,1),1);
                                                    cblas_dcopy(n1[2], streams(f1, A, 1 , 2)+(k+di)*n1[2],1,streams(f1,diagonalVectorA,1,2),1);
                                                    f1->sinc.tulip[diagonalVectorA].Current[1] = 1;;
                                                    if ( inner(0, f1, diagonalVectorA, 1) < 1e-6)
                                                        f1->sinc.tulip[diagonalVectorA].Current[1] = 0;;

                                                    tAddTw(f1, build+classicalBasisSize, 1,diagonalVectorA,0 );
                                                    
                                                }
                                                tScale(f1, build+classicalBasisSize,1./magnitude(f1, build+classicalBasisSize));

//                                                if ( r1 == r2 && r2 == r3 )
//                                                    f1->sinc.tulip[build+classicalBasisSize].center = r1;
//                                                else
//                                                    f1->sinc.tulip[build+classicalBasisSize].center = -1;
                                                f1->sinc.tulip[build+classicalBasisSize].center = 0;
                                                
                                            }else if ( ops  == - (classicalBasisSize+1) ){
                                               // printf("-->%d\n",classicalBasisSize);

                                                mmm[0] = i;
                                                mmm[1] = j;
                                                mmm[2] = k;
                                                return 1;
                                            }
                                            classicalBasisSize++;

                                            if ( classicalBasisSize >= xB && ops== 2 ){
                                                return xB;
                                            }
                                            
                                        }
                                        kk++;
                                    }
                            jj++;
                        }
                    }
                ii++;
            }
    return 0;
    return classicalBasisSize;
}



INT_TYPE tFilter(struct field * f1, INT_TYPE Ve, INT_TYPE type, enum division usr){
    INT_TYPE ii,j,cmpl,rank,flag ;
    double value;
#ifdef OMP
#pragma omp parallel for private (ii,cmpl,rank,j) schedule(dynamic,1)
#endif
    for ( ii = 0; ii < Ve ; ii++)
    {
        
#ifdef OMP
        rank = omp_get_thread_num();
#else
        rank = 0;
#endif
        for ( cmpl = 0 ; cmpl < 2 ; cmpl++){
            tEqua(f1, copyTwoVector,rank, usr+ii,cmpl);
            f1->sinc.tulip[copyThreeVector].Current[rank] = 0;
            for ( j = 0 ; j < CanonicalRank(f1, copyTwoVector, rank);j++){
                f1->sinc.tulip[permutationVector].Current[rank] = 0;
                tBuildIrr(rank, f1, type, ocean(rank,f1,copyTwoVector,j,rank), rank, permutationVector, rank);
                if (fabs( inner(rank, f1, permutationVector, rank) ) > 1e-15 ){
                    f1->sinc.tulip[copyFourVector].Current[rank] = 0;

                    tCycleDecompostionOneMP(rank, f1, permutationVector, rank, copyFourVector, rank, f1->mem1->rt->vCONVERGENCE, 1, -1);
                    tAddTw(f1, copyThreeVector, rank, copyFourVector, rank);
                }
            }
            tCycleDecompostionOneMP(rank, f1, usr+ii, cmpl, copyThreeVector, rank, f1->mem1->rt->vCONVERGENCE, part(f1,usr+ii), -1);
        }
    }

	for ( ii = 0 ; ii < Ve ; ii++){
        value = magnitude(f1, usr+ii);

        if ( fabs(value) < 1e-15 ){
            printf("oopsy %f %d %d\n",value,ii,CanonicalRank(f1, usr+ii, 0));
            exit(0);
        }
        tScale(f1, usr+ii, 1./value);
    }
    return 0;
}

INT_TYPE tSelect(struct field * f1, INT_TYPE Ve, INT_TYPE type, enum division usr, enum division usa, INT_TYPE testFlag){
    INT_TYPE info,rank,maxEV = f1->mem1->rt->maxEV,stride = maxEV,n,m;;
    double value;
    
    tClear(f1, usr+Ve);
    tClear(f1, permutationVector);
    tBuildIrr(0, f1, type, usa, 0, permutationVector, 0);
    tBuildIrr(0, f1, type, usa, 1, permutationVector, 1);
    
    if (part(f1,usr+Ve) >= CanonicalRank(f1,permutationVector,0)){
        tEquals(f1, usr+Ve, permutationVector);
    } else {
        tCycleDecompostionOneMP(0, f1, permutationVector, 0, usr+Ve, 0, f1->mem1->rt->vCONVERGENCE, part(f1,usr+Ve), -1);
        tCycleDecompostionOneMP(0, f1, permutationVector, 1, usr+Ve, 1, f1->mem1->rt->vCONVERGENCE, part(f1,usr+Ve), -1);
    }
    value = magnitude(f1, usr+Ve);

    
    
    if ( value < 1e-6 || isnan(value) || isinf(value) )
        return 0;
    tScale(f1, usr+Ve, 1./value);

    DCOMPLEX * S = (DCOMPLEX*)(myStreams(f1, matrixSbuild, 0));
    DCOMPLEX * St = (DCOMPLEX*)(myStreams(f1, matrixSbuild, 0))+maxEV*maxEV;

    double * ov = myStreams(f1, twoBodyRitz, 0);
    assignCores(f1, 2);
#ifdef OMP
#pragma omp parallel for private (m,n,rank) schedule(dynamic,1)
#endif
    for ( n = 0; n < Ve+1 ; n++)
    {
        
#ifdef OMP
        rank = omp_get_thread_num();
#else
        rank = 0;
#endif
        
        for ( m = 0 ; m <= n   ; m++)    {
            S[n*stride+m] =
            (tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usr+n,0,'N', usr+m,0)
             +tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usr+n,1,'N', usr+m,1)
             +I*tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usr+n,1,'N', usr+m,0)
             -I*tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usr+n,0,'N', usr+m,1)
             );
            S[m*stride+n] = conj(S[n*stride+m]);
        }
        
    }
    
    assignCores(f1, 0);
    cblas_zcopy(maxEV*maxEV,S,1,St,1);
    tzheev(0, f1, 'N', Ve+1, St, stride, ov);

    if ( testFlag ){

        if ( ov[Ve]/ov[0] < 1e5 && ov[Ve]/ov[0] >0 ){

            if ( 1){
                printf(",%lld ",  Ve);
                fflush(stdout);

                return 1;
            }
            else{
                printf("selected vector failed to be in right class, %d %d\n", type,tClassify(0, f1, usr+Ve) );
                return 0;
            }
    }
    }
    
    return 0;
}


INT_TYPE tCollect (struct field * f1, INT_TYPE type,enum division usz, INT_TYPE target){
    INT_TYPE ii, Ve = 0;
    f1->sinc.tulip[diagonalVectorA].header = Cube;
    f1->sinc.tulip[diagonalVectorB].header = Cube;

    printf("Collect %d\n", target);
    char c,c0,c1;
    if ( f1->body == two || f1->body == one){
        c0 = 'N';
        c1 = 'N';
    }
    else if ( f1->body == three ){
        c0 = 'A';
        c1 = 'E';
    }
    else if ( f1->body == four ){
        c0 = 'a';
        c1 = 'w';
    }else {
        printf("bod)");
        exit(0);
    }
    if ( f1->body < four ){
        
        INT_TYPE mmm[3];
        INT_TYPE *n1 = vectorLen(f1, eigen);
        double value;
        for ( value = 1.0 ; value > 0 ; value-=0.01){
            ii = 0;
#if VERBOSE
            printf("%f\n", value);
#endif
            while (tFoundationLevel(f1, eigen,value-0.01,value ,-((ii++)+1),usz,target,0.,0.,0.,mmm) ){
                tClear(f1,diagonalVectorA);
                f1->sinc.tulip[diagonalVectorA].Current[0] = 1;;
                f1->sinc.tulip[diagonalVectorA].Current[1] = 1;;

                cblas_dcopy(n1[0], streams(f1, eigen, 0 , 0)+(mmm[0])*n1[0],1,streams(f1,diagonalVectorA,0,0),1);
                cblas_dcopy(n1[1], streams(f1, eigen, 0 , 1)+(mmm[1])*n1[1],1,streams(f1,diagonalVectorA,0,1),1);
                cblas_dcopy(n1[2], streams(f1, eigen, 0 , 2)+(mmm[2])*n1[2],1,streams(f1,diagonalVectorA,0,2),1);
                cblas_dcopy(n1[0], streams(f1, eigen, 1 , 0)+(mmm[0])*n1[0],1,streams(f1,diagonalVectorA,1,0),1);
                cblas_dcopy(n1[1], streams(f1, eigen, 1 , 1)+(mmm[1])*n1[1],1,streams(f1,diagonalVectorA,1,1),1);
                cblas_dcopy(n1[2], streams(f1, eigen, 1 , 2)+(mmm[2])*n1[2],1,streams(f1,diagonalVectorA,1,2),1);
                f1->sinc.tulip[diagonalVectorA].Current[1] = 1;;
                if ( inner(0, f1, diagonalVectorA, 1) < 1e-6)
                    f1->sinc.tulip[diagonalVectorA].Current[1] = 0;;

                //for (c = -1 ; c +c0 <= c1 ; c++)
                {
                    // tAddUpComponents(0, f1, diagonalVectorA, 0, diagonalVectorA, 0, up);
                    //printf("Y%f\n", up[type]);
                    //if ( fabs(up[type]) > 1e-3 )
                    {
                        Ve += tSelect(f1, Ve, type, usz, diagonalVectorA, 1);
                    }
                   // printf("***%d %d\n",Ve,target);
                    if ( Ve == target )
                        return Ve;
                }
            }
        }
    }else {
        INT_TYPE mmm[6],jj,ii;
        INT_TYPE *n1 = vectorLen(f1, eigen);
        double value,value2;
        for ( value = 1.0 ; value > 0 ; value-=0.01){
            
                ii = 0;
                while (tFoundationLevel(f1, eigen,value-0.01,value ,-((ii++)+1),usz,target,0.,0.,0.,mmm)){
                    for ( value2 = 1.0 ; value2 >= value ; value2-=0.01){
                    jj = 0;
                    while ( tFoundationLevel(f1, eigen,value2-0.01,value2 ,-((jj++)+1),usz,target,0.,0.,0.,mmm+3)  ){
                       // printf("%f %f %d %d\n", value,value2,ii,jj);
                        tClear(f1,diagonalVectorA);
                        f1->sinc.tulip[diagonalVectorA].Current[0] = 1;;
                        f1->sinc.tulip[diagonalVectorA].Current[1] = 1;;

                        tClear(f1,diagonalVectorB);
                        f1->sinc.tulip[diagonalVectorB].Current[0] = 1;;
                        f1->sinc.tulip[diagonalVectorB].Current[1] = 1;;

                        cblas_dcopy(n1[0], streams(f1, eigen, 0 , 0)+(mmm[0])*n1[0],1,streams(f1,diagonalVectorA,0,0),1);
                        cblas_dcopy(n1[1], streams(f1, eigen, 0 , 1)+(mmm[1])*n1[1],1,streams(f1,diagonalVectorA,0,1),1);
                        cblas_dcopy(n1[2], streams(f1, eigen, 0 , 2)+(mmm[2])*n1[2],1,streams(f1,diagonalVectorA,0,2),1);
                        cblas_dcopy(n1[0], streams(f1, eigen, 1 , 0)+(mmm[0])*n1[0],1,streams(f1,diagonalVectorA,1,0),1);
                        cblas_dcopy(n1[1], streams(f1, eigen, 1 , 1)+(mmm[1])*n1[1],1,streams(f1,diagonalVectorA,1,1),1);
                        cblas_dcopy(n1[2], streams(f1, eigen, 1 , 2)+(mmm[2])*n1[2],1,streams(f1,diagonalVectorA,1,2),1);
                        f1->sinc.tulip[diagonalVectorA].Current[1] = 1;;
                        if ( inner(0, f1, diagonalVectorA, 1) < 1e-6)
                            f1->sinc.tulip[diagonalVectorA].Current[1] = 0;;

                        
                        cblas_dcopy(n1[0], streams(f1, eigen, 0 , 0)+(mmm[0+3])*n1[0],1,streams(f1,diagonalVectorB,0,0),1);
                        cblas_dcopy(n1[1], streams(f1, eigen, 0 , 1)+(mmm[1+3])*n1[1],1,streams(f1,diagonalVectorB,0,1),1);
                        cblas_dcopy(n1[2], streams(f1, eigen, 0 , 2)+(mmm[2+3])*n1[2],1,streams(f1,diagonalVectorB,0,2),1);
                        cblas_dcopy(n1[0], streams(f1, eigen, 1 , 0)+(mmm[0+3])*n1[0],1,streams(f1,diagonalVectorB,1,0),1);
                        cblas_dcopy(n1[1], streams(f1, eigen, 1 , 1)+(mmm[1+3])*n1[1],1,streams(f1,diagonalVectorB,1,1),1);
                        cblas_dcopy(n1[2], streams(f1, eigen, 1 , 2)+(mmm[2+3])*n1[2],1,streams(f1,diagonalVectorB,1,2),1);
                        f1->sinc.tulip[diagonalVectorB].Current[1] = 1;;
                        if ( inner(0, f1, diagonalVectorB, 1) < 1e-6)
                            f1->sinc.tulip[diagonalVectorB].Current[1] = 0;;

                        tClear(f1,copyVector);
                        tOuterProductSu(f1, diagonalVectorA, 0, diagonalVectorB, 0, copyVector, 0);
                        tOuterProductSu(f1, diagonalVectorA, 1, diagonalVectorB, 1, copyVector, 1);
                        
                        for (c = -1 ; c +c0 <= c1 ; c++)
                        {
                            // tAddUpComponents(0, f1, diagonalVectorA, 0, diagonalVectorA, 0, up);
                            //printf("Y%f\n", up[type]);
                            //if ( fabs(up[type]) > 1e-3 )
                            {
                                Ve += tSelect(f1, Ve, type, usz, copyVector, 1);
                            }
                            // printf("***%d %d\n",Ve,target);
                            if ( Ve == target )
                                return Ve;
                        }
                    }
                }
            }
            
        }

    }
    return Ve;
}


INT_TYPE tSquareVectors(struct calculation * c1, INT_TYPE EV2, enum division usz,enum division usr ){
    INT_TYPE i,j;
    struct field *f1 = & c1->i.c;

    for ( i = 0; i < EV2 ; i++)
        for ( j = 0 ; j < EV2 ; j++)
        {
            tClear(f1,usz+i*EV2+j );
            tOuterProductSu(f1, usz+i, 0, usz+j, 0, usr+i*EV2+j, 0);
            tOuterProductSu(f1, usz+i, 1, usz+j, 1, usr+i*EV2+j, 1);
        }    
    return EV2*EV2;
}




INT_TYPE tOCSB (struct calculation * c1 , enum division usz )
{
    struct field *f1 = & c1->i.c;
    double lvl[5];
    INT_TYPE EV = 0,space;
    INT_TYPE target,count;
    double value,valueFloor;

    
        for ( space = 0; space < SPACE ; space++){
            INT_TYPE target= 0;
            double value;
            if ( f1->rds.flag == 1 )
                target = f1->rds.Basis[space];
            if ( f1->rds.flag == 2 )
                target = f1->rds.Basis2[space];
            if ( f1->rds.flag == 3 )
                target = f1->rds.Basis3[space];
            if ( f1->rds.flag == 4 )
                target = f1->rds.Basis4[space];
            INT_TYPE iter= 0;
            INT_TYPE count = target+1, prevCount,same = 0 ;
            double max = 1., min = 0.;
            do {
                prevCount = count;
                if ( iter ){
                    if (count < target  ){
                        max = 0.5*(max+min);
                    }else {
                        min = 0.5*(max+min);
                    }
                }
                value = 0.5*(max+min);
                count = tBasisLevel(f1, eigen, space, value, 0, basis, target);
                if ( count == prevCount )
                    same++;
                else
                    same = 0;
                iter++;
                if ( same > 100 )
                {
                    if ( count >= target ){
                        break;
                    } else {
                      //  printf("%lld -> %f %lld\n", space,value, count);
                        target = count;
                        break;
                        //exit(0);
                    }
                }
            }while ( 0 < abs( count - target )  || count < target  );
            lvl[space] = value;
            count = tBasisLevel(f1, eigen, space, value, 2, basis, target);
            printf("%lld -> %f %lld\n", space,value, count);
            fflush(stdout);
            
            if ( count != target ){
                printf("bailing\n");
                exit(0);
            }
        }
    
    //                        if(0){
    //                            // GRAM SCHMIDT
    //                            // GS
    //                            INT_TYPE * n1 = vectorLen(f1, usz);
    //                            INT_TYPE * m1 = vectorLen(f1,eigen);
    //                            INT_TYPE info;
    //                            for ( space = 0; space < SPACE ; space++){
    //
    //                                info = tdgeqr(0,f1,n1[space],m1[space],streams(f1, basis, 0 , space),m1[space]);
    //
    //                                for ( space = 0; space < SPACE ; space++)
    //                                    info = tdgeqr(0,f1,vectorLen(f1,usz)[space],EV,streams(f1,usz,0,space),vectorLen(f1, usz)[space] );
    //
    //                                if ( info != 0 ){
    //                                    printf("infos %lld %lld\n", info, space);
    //                                    exit(0);
    //                                }
    //                            }
    //                        }
    
    target = (c1->i.qFloor);
    INT_TYPE iter= 0;
    count = target;
    INT_TYPE prevCount,same = 0 ;
    double max = 1., min = 0.;
    do {
        prevCount = count;
        if ( iter    ){
            if (count < target+1  ){
                max = 0.5*(max+min);
            }else {
                min = 0.5*(max+min);
            }
        }
        value = 0.5*(max+min);
        
        if ( c1->i.sectors ){
            count = tFoundationLevel(f1, eigen,0.,value,0,usz,target,lvl[0],lvl[1],lvl[2],NULL);
        }
        
        if ( count == prevCount )
            same++;
        else
            same = 0;
        if ( same > 1000){
            if ( count >= target ){
                break;
            } else {
                printf("3 -> %f %lld\n", value, count);
                
                exit(0);
            }
        }
        iter++;
    } while ( (0.1*target < abs( count - target ) ) || count < target );

    c1->i.qFloor = target;
    valueFloor = value;
    EV += tFoundationLevel(f1, eigen,0.,valueFloor,2,usz ,target,lvl[0],lvl[1],lvl[2],NULL);
    printf("3 -> %f %lld\n", value, count);
    fflush(stdout);
    return EV;
}

INT_TYPE tGreatDivideIteration ( struct field * f1, enum division A , INT_TYPE I1, INT_TYPE I2, enum division usz, INT_TYPE foundation,INT_TYPE nMult, double shift){
    INT_TYPE expon,info;
    INT_TYPE rank ,i;
    //time_t start_t, lapse_t;
    double temp;
    //time(&start_t);
    INT_TYPE iii = 0;
    assignCores(f1, 1);

    for( expon = 1 ; foundation*expon < nMult  ; expon++){
        
#ifdef OMP
#pragma omp parallel for private (iii,rank,temp) schedule(dynamic,1)
#endif
        for ( iii = 0; iii < foundation ; iii++)
        {
            
#ifdef OMP
            rank = omp_get_thread_num();
#else
            rank = 0;
#endif
            
            
            {
                
                tEquals(f1,usz+iii+expon*foundation,usz+(expon-1)*foundation+iii );
                tHXpX(rank, f1, A, 0, 1.0, 0.0, usz+iii+expon*foundation, f1->mem1->rt->vCANON , part(f1,usz+(expon)*foundation+iii));
                
                temp = inner(rank, f1, usz+(expon)*foundation+iii, 0)+inner(rank, f1, usz+(expon)*foundation+iii, 1);
                
               // printf("%d(%d) -> %d(%d)::[%f]\n",(expon-1)*foundation+iii,part(f1,usz+(expon-1)*foundation+iii ), expon*foundation+iii,part(f1,usz+(expon)*foundation+iii ), temp- sqr( tMultiplyMP(rank, &info, f1, 1.0, -1, nullVector, 0, 'T', usz+(expon)*foundation+iii, 0, 'N', usz+(expon-1)*foundation+iii, 0) + tMultiplyMP(rank, &info, f1, 1.0, -1, nullVector, 0, 'T', usz+(expon)*foundation+iii, 1, 'N', usz+(expon-1)*foundation+iii, 1)) );
                tScale(f1, usz+iii+expon*foundation, 1./sqrt(temp));

//                if ( expon == 1 )
//                    ev[foundation] = tMultiplyMP(rank, &info,f1, 1., -1, nullVector, 0, 'T', usz+iii+expon*foundation, 0, 'N',usz+(expon-1)*foundation+iii, 0)+tMultiplyMP(rank,&info, f1, 1., -1, nullVector, 0, 'T', usz+iii+expon*foundation, 1, 'N',usz+(expon-1)*foundation+iii, 1);
//
                //fflush(stdout);
          //      ev[iii ]= magnitude(f1,usz+iii+expon*foundation);
//                printf(":: \t %f\n",ev[iii] );
//                fflush(stdout);
            }
            
        }
    }
    
    
    return nMult;
}

INT_TYPE tEdges(struct calculation *c1){
    struct field * f1 = &c1->i.c;
    INT_TYPE info;
    enum body bootBodies = f1->body;
    if ( 1 ){
        //EDGES ALT
        enum block b,bx;
        INT_TYPE iii,jjj=1,dim,type;
        double sum = 0;
        for ( type = 1 ;type <= 24 ; type++){
            jjj=1;
        for ( iii = 0 ; iii < c1->i.nStates ; iii++)
            if ( f1->sinc.tulip[eigenVectors+iii].symmetry  == type&& (! c1->i.type || type == c1->i.type))
            {
//            tTransformFrom(f1, eigenVectors+iii, squareVector, basis);
                bx = tv1;
             if ( bootBodies == two )
                bx = tv2;
            else if ( bootBodies == three )
                bx = tv3;
            else if ( bootBodies == four )
                bx = tv4;
                printf("Edge%d:\t%d\t",jjj++,type);
            for ( b = tv1 ; b <= bx; b++){
                printf("%d:\t",b+1);
                
                f1->sinc.tulip[edgeMatrix].blockType = b;
                sum = 0;
                for ( dim = 0 ; dim < 3 ;dim++){
                    
                    tMultiplyMP(0, &info,f1, 1.0, -1, productVector, 0, 'N', ocean(0, f1, edgeMatrix, dim, 0), 0, 'N', eigenVectors+iii, 0);
                    sum += sqr(tMultiplyMP(0, &info, f1, 1.0, -1, nullVector, 0, 'T', productVector, 0, 'N', eigenVectors+iii, 0));
                }
                printf("%1.9f\t", sqrt(sum));
                
                sum = 0;
                for ( dim = 3 ; dim < 6 ;dim++){
                    tMultiplyMP(0,  &info,f1, 1.0, -1, productVector, 0, 'N', ocean(0, f1, edgeMatrix, dim, 0), 0, 'N', eigenVectors+iii, 0);
                    sum += sqr(tMultiplyMP(0,  &info,f1, 1.0, -1, nullVector, 0, 'T', productVector, 0, 'N', eigenVectors+iii, 0));
                }
                printf("%1.9f\t", sqrt(sum));
            }
            printf("\n");
            }
        }
    }
    return 0;
}

//INT_TYPE tConvergeTest (struct calculation * c1, enum division input){
//
//
//    {
//        struct field * f1 = &c1->i.c;
//        double sum2,sum;
//        enum division Mat, leftP,ml[1000];
//        INT_TYPE Lane = MaxCore,iii,i,i2,j,j2,rk,s1=0,cat=0,rank,iiii= 0;
//        double Sa [c1->i.bRank  *c1->i.bRank];
//        double Saa [c1->i.bRank  ];
//
//        leftP = input;
//
//
//        do {
//
//            if ( linear == name(f1,leftP) ){
//                Mat = rivers(Lane++, f1, linear, cat );
//                f1->sinc.tulip[Mat].blockType = f1->sinc.tulip[leftP].blockType;
//            }else {
//                Mat = leftP;
//            }
//
//            if ( Lane -MaxCore >= MaxCore )
//            {
//                printf("neeed more \n");
//                exit(0);
//            }
//
//            if ( CanonicalRank(f1, Mat, 0)+CanonicalRank(f1, Mat, 1))
//                ml[iiii++] = Mat;
//
//            if (name(f1, leftP) == linear ){
//                cat++;
//
//                if ( cat > f1->Na ){
//                    leftP = f1->sinc.tulip[leftP].linkNext;
//                    cat = 0;
//                }
//
//            } else {
//                leftP = f1->sinc.tulip[leftP].linkNext;
//            }
//        } while ( leftP != nullName);
//
//
//
//        for ( iii = 0; iii < c1->i.nStates ; iii++){
//            sum = 0.;
//            sum2 = 0.;
//            rk = CanonicalRank(f1, eigenVectors+iii, 0);
//
//            for ( j = 0; j < iiii ; j++)
//            {
//
//#pragma omp parallel for private (i,i2,j2,rank) schedule(dynamic,1)
//                for ( i = 0; i < rk;i++)
//                {
//#if ARES
//
//                    rank = omp_get_thread_num();
//#else
//                    rank = 0;
//#endif
//
//                    //                                if ( Nb){
//                    //                                    tEqua(f1,squareVector, rank, ocean(rank,f1,eigenVectors+iii,i,0),0);
//                    //                                    tpProduct(rank, f1, ml[j], 0, squareVector, rank);
//                    //                                }
//                    //                                else
//                    tMultiplyMP(rank,f1,1., -1,squareVector, rank, 'N', ml[j] , 0,'N', ocean(rank,f1,eigenVectors+iii,i,0), s1);
//
//
//                    Saa[i] = tMultiplyMP(rank,f1,1., -1,nullVector, 0, CDT, squareVector, rank,'N', ocean(rank,f1,eigenVectors+iii,i,0),0);
//
//                    for ( j2 = 0 ; j2 < iiii ; j2++)
//                        for ( i2 = 0 ; i2 <= i   ; i2++){
//                            //                                        if ( Nb){
//                            //                                            tEqua(f1, totalVector, rank, ocean(rank,f1,eigenVectors+iii,i2,0),0);
//                            //                                            tpProduct(rank, f1, ml[j2] , 0, totalVector, rank);
//                            //                                        }
//                            //                                        else
//                            tMultiplyMP(rank,f1,1., -1,totalVector, rank, 'N',ml[j2] , 0,'N', ocean(rank,f1,eigenVectors+iii,i2,0), s1);
//
//                            Sa[i*rk+i2] = tMultiplyMP(rank,f1,1., -1,nullVector, 0, CDT, squareVector,rank,'N', totalVector, rank);
//
//                        }
//                }
//
//                for ( i = 0; i < rk ; i++){
//                    for ( i2 = 0; i2 < rk ; i2++){
//                        sum2 += Sa[i*rk+i2];
//                    }
//                    sum += Saa[i];
//                }
//            }
//
//            printf("Cv%lld \t %f\n", iii+1, sum2-sqr(sum));
//
//        }
//    }
//}



INT_TYPE tEigenLoad (struct field * f1, enum division A ,INT_TYPE type,  INT_TYPE Ne, enum division usz, INT_TYPE quantumBasisSize ,INT_TYPE foundation, INT_TYPE flag,  enum division outputSpace, enum division outputValues){
    INT_TYPE gvOut,prevBuild;
    time_t start_t, lapse_t;
    time(&start_t);
    INT_TYPE countLam = 0,countTot = 0,cl;
    enum division Mat;
    INT_TYPE cmpl,cmpl2,cmpl3,cat,iii = 0,maxEV = f1->mem1->rt->maxEV,rank;
    double * ritz = myStreams(f1, outputValues, 0);
    enum division el ;
    DCOMPLEX co ;
    DCOMPLEX *T  =  (DCOMPLEX *) myStreams(f1, matrixHbuild,0/*CORE RANK*/);
    double *H  =   myStreams(f1, vectorHbuild,0/*CORE RANK*/);
    double *vectors[2];
    vectors[0] =   myStreams(f1, matrixHbuild,0/*CORE RANK*/)+4*maxEV*maxEV;
    vectors[1] =   myStreams(f1, matrixHbuild,0/*CORE RANK*/)+5*maxEV*maxEV;
    DCOMPLEX *S  =  (DCOMPLEX *) myStreams(f1, matrixSbuild,0/*CORE RANK*/);
    INT_TYPE powerMat,stride;

    INT_TYPE info,n,m,s1,g,r,rr,rx,gx,a,a2;
    enum division leftP ;
    stride = maxEV;
    prevBuild = CanonicalRank(f1, matrixHbuild, 0);
    tClear  (f1, copyTwoVector);
    for ( n = prevBuild; n < quantumBasisSize ; n++)
    {
        tScale(f1, usz+n, 1./magnitude(f1, usz+n));
        H[n] = 1.;
    }
    
    assignCores(f1, 2);


#ifdef OMP
#pragma omp parallel for private (m,n,rank) schedule(dynamic,1)
#endif
    for ( n = prevBuild; n < quantumBasisSize ; n++)
    {

#ifdef OMP
        rank = omp_get_thread_num();
     #else
        rank = 0;
#endif

                for ( m = 0 ; m <= n   ; m++)    {
                    S[n*stride+m] =
                    (tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usz+m,0,'N', usz+n,0)
                    +tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usz+m,1,'N', usz+n,1)
                    +I*tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usz+m,1,'N', usz+n,0)
                    -I*tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usz+m,0,'N', usz+n,1)
                     )/H[n]/H[m];

                    S[m*stride+n] = conj(S[n*stride+m]);
                    {
                        T[n*stride+m] = 0.;
                        T[m*stride+n] = (T[n*stride+m]);
                    }
                }
            
        }
    if (0)
        for ( n = 0; n < quantumBasisSize ; n++){
            for ( m = 0; m < quantumBasisSize ; m++){
          //  if ( fabs(S[n*stride+m]-1.0)>1e-3 || isnan(S[n*stride+m] || isinf(S[n*stride+m))
            {
                printf("%lld :%f +i%f \n", n,creal(S[n*stride+m]), cimag(S[n*stride+m]) );
                fflush(stdout);
            }
        }
        }
    
            
    
            for ( s1 = 0 ; s1 < 1; s1++){
                leftP = A;
                cat = 0;

                do {

                    if ( linear == name(f1,leftP) ){
                        Mat = rivers(0, f1, linear, cat );
                        //printTrace(f1,Mat);
                        f1->sinc.tulip[Mat].blockType = f1->sinc.tulip[leftP].blockType;
                    }else {
                        Mat = leftP;
                    }
                    for ( cmpl = 0; cmpl < 2  ;cmpl++)
                        for ( cmpl2 = 0 ; cmpl2 < 2  ;cmpl2++)
                                if ( CanonicalRank(f1, name(f1, Mat), cmpl)){
                                    co = 1.;
                                    if ( cmpl )
                                        co *= I;
                                    if ( cmpl2 )
                                        co *= I;

                                    printf("%d-%d-%d\t%d\t%d\t%d \t %d\n",name(f1,Mat),cmpl,cmpl2, bodies(f1, Mat),f1->sinc.tulip[Mat].blockType,CanonicalRank(f1, name(f1,Mat), cmpl), cat );
                                    fflush(stdout);
    #ifdef OMP
    #pragma omp parallel for private (m,rank,n) schedule(dynamic,1)
    #endif
                                    for ( n = prevBuild; n < quantumBasisSize ; n++)
                                    {
                                        
#ifdef OMP
                                        rank = omp_get_thread_num();
#else
                                        rank = 0;
#endif
                                        for ( m = 0 ; m <=n   ; m++){
                                            (T+maxEV*maxEV)[n*stride+m] = matrixElements(rank, f1, usz+m, 0, 'N', Mat, cmpl, usz+n, cmpl2);
                                            (T+maxEV*maxEV)[n*stride+m] += -I*matrixElements(rank, f1, usz+m, 1, 'N', Mat, cmpl, usz+n, cmpl2);
                                            (T+maxEV*maxEV)[n*stride+m] *= co/H[m]/H[n];
                                        }
                                        
                                        
                                        
                                    }




                        for ( n = prevBuild; n < quantumBasisSize ; n++)
                            for ( m = 0; m <= n ; m++){
                                T[n*stride+m] += (T+maxEV*maxEV)[n*stride+m];
                                T[m*stride+n]  = conj(T[n*stride+m]);
                            }

                    }
                    for ( n = 0; n < quantumBasisSize ; n++)
                        for ( m = 0 ; m <=n   ; m++){
                            if ( cabs((T+maxEV*maxEV)[n*stride+m])  < 1e-6 ){
                                countLam++;
                            }
                            countTot++;
                        }


                    if (name(f1, leftP) == linear ){
                        cat++;

                        if ( cat > f1->Na ){
                            leftP = f1->sinc.tulip[leftP].linkNext;
                            cat = 0;
                        }

                    } else {
                        leftP = f1->sinc.tulip[leftP].linkNext;
                    }
                } while ( leftP != nullName);

            }
//    if(0)
//    for ( n = 0; n < quantumBasisSize-foundation; n++)
//        for ( m = 0 ; m <= n   ; m++){
//            if ( part(f1,usz+m) > minRank && part(f1, usz+n) > minRank && n % foundation == m % foundation){
//               // printf("%lld %lld %f %f\n", n/foundation,m/foundation,T[n*stride+m],T[n*stride+m+maxEV*maxEV*2]);
//                T[n*stride+m]= T[n*stride+m+maxEV*maxEV*2];
//                T[m*stride+n] = T[n*stride+m];
//            }
//        }
    //assignCores(f1, 1);
    
//    f1->mem1->rt->matrixElementsTotal += countTot;
//    f1->mem1->rt->matrixElementsZero += countLam;
//    time(&lapse_t);
//    f1->mem1->rt->buildTime += difftime(lapse_t, start_t);
//    time(&start_t);

    if (1){
        assignCores(f1, 0);

        char Job = 'N';
        if (flag)
            Job = 'V';

        cblas_zcopy(maxEV*maxEV , T , 1 , T+maxEV*maxEV , 1);
        cblas_zcopy(maxEV*maxEV , S , 1 , S+maxEV*maxEV , 1);

        tzheev(0, f1, 'N', quantumBasisSize, S+maxEV*maxEV, stride, ritz);
        
        if ( 1e-13 < ritz[quantumBasisSize-1]/ritz[0] && ritz[quantumBasisSize-1]/ritz[0] < 1e13 )
            printf("Condition Krylov \t %f \n",  ritz[quantumBasisSize-1]/ritz[0]);
        else {
            printf("Linear dependent!\n");
            return -1;
        }
        
        cblas_zcopy(maxEV*maxEV , S , 1 , S+maxEV*maxEV , 1);
        gvOut = tzhegv (0,f1,Job,quantumBasisSize,T+maxEV*maxEV,S+maxEV*maxEV,stride,ritz);
        assignCores(f1, 1);

        for ( iii = 0; iii < maxEV*maxEV ; iii++){
            vectors[0][iii] = creal((T+maxEV*maxEV)[iii]);
            vectors[1][iii] = cimag((T+maxEV*maxEV)[iii]);
        }
    }
    
    time(&lapse_t);

    //f1->mem1->rt->eigenTime += difftime(lapse_t, start_t);
    time(&start_t);
    
    {       //printf("\nHeader,Number,CEG,Linear,BODY,CLASS,WEIGHT\n");
        for ( iii = 0; iii < imin(quantumBasisSize,Ne) ; iii++)
        {
            f1->sinc.tulip[eigenVectors+iii].value = ritz[iii];
          //  printf("State%lld:%lld:,%lld ,%1.15f, %f, %d, %d , %f\n", iii+1, f1->sinc.N1,ritz[iii],iii+1,  sqr(cblas_dnrm2(quantumBasisSize, vectors[0]+(iii)*stride, 1))+sqr(cblas_dnrm2(quantumBasisSize, vectors[1]+(iii)*stride, 1)),bodies(f1,eigenVectors),type, deg(f1, type));
        }
        fflush(stdout);
    }
    if (flag > 1 ){
        double norms;
        if (1){
            f1->sinc.tulip[eigenList].name = usz;
            f1->sinc.tulip[eigenList].purpose = ptObject;
            
            f1->sinc.tulip[eigenList].ptRank[0] = 0;
            f1->sinc.tulip[eigenList].Current[0]= 0;
            for( g = 0; g < quantumBasisSize ; g++)
                for ( cmpl = 0 ;cmpl < spins(f1,usz+iii) ; cmpl++)
                    for ( r = 0; r < part(f1, usz+g); r++){
                        f1->sinc.tulip[eigenList].Current[0]++;//memory inline. need this.
                    }
            // printf("[[%d %d %lld %lld]]\n", eigenList, usz,f1->sinc.tulip[usz].Address, streams(f1, eigenList, 0, 0 )-streams(f1, usz, 0, 0 ));
            
            for ( cmpl = 0 ;cmpl < spins(f1,usz+iii) ; cmpl++)
                
                for ( powerMat = 1 ; powerMat <= 1 ; powerMat++){
                    
                    el = outputSpace + (powerMat-1)*Ne;
                    
#ifdef OMP
#pragma omp parallel for private (iii,rank,rr,g,r,cmpl2,norms) schedule(dynamic,1)
#endif
                    for ( iii = 0; iii < imin(quantumBasisSize,Ne) ; iii++)
                    {
                        
#ifdef OMP
                        rank = omp_get_thread_num();
#else
                        rank = 0;
#endif
                        
                        rr = 0;
                        norms = 0.;
                        for( g = 0; g < quantumBasisSize ; g++){
                            for ( cmpl2 = 0; cmpl2 < spins(f1, usz+g) ; cmpl2++)
                                for ( r = 0; r < part(f1, usz+g); r++){
                                    if (  r >= CanonicalRank(f1, usz+g, cmpl2)){
                                        myStreams(f1, canonicalBuffersC, rank)[rr++] = 0.0;
                                    }else{
                                        if( cmpl2 == cmpl ){
                                            norms += sqr((vectors[cmpl]+(iii)*stride)[g]);
                                            myStreams(f1, canonicalBuffersC, rank)[rr++] = (vectors[cmpl]+(iii)*stride)[g];
                                        }
                                        else
                                            myStreams(f1, canonicalBuffersC, rank)[rr++] = 0.0;//mask off other *SPIN*( or complex vector).
                                    }
                                //    printf("%d %d %d %f\n", usz+g,cmpl2,r,myStreams(f1, canonicalBuffersC, rank)[rr-1] );
                                }
                        }
                        if ( rr > part(f1, canonicalBuffersC))
                        {
                            printf("crap!\n %lld %lld",rr , part(f1, canonicalBuffersC));
                            exit(0);
                        }
                        f1->sinc.tulip[el+iii].Current[cmpl] = 0;
                     //   printf("\n%d %d %d %d\n", iii+1, cmpl, rr,CanonicalRank(f1, eigenList,0));
                        tCycleDecompostionListOneMP(rank,f1, eigenList, myStreams(f1, canonicalBuffersC,rank), el+iii, cmpl, f1->mem1->rt->vCANON, part(f1, el+iii), 1.);
                        
                       // printf("%d -> %f (%f)\n", iii+1, norms,H[iii]);
                    }
                }
        }
        time(&lapse_t);
        printf("\nSpectrum TIME \t %15.0f\n", difftime(lapse_t, start_t));
        
        
        
       if(1) {
            
#ifdef OMP
#pragma omp parallel for private (iii,rank) schedule(dynamic,1)
#endif
            for ( iii = 0; iii < imin(quantumBasisSize,Ne) ; iii++)
            {
                
#ifdef OMP
                rank = omp_get_thread_num();
#else
                rank = 0;
#endif
                tScale(f1, eigenVectors+iii, 1./sqrt( inner(rank, f1, eigenVectors+iii, 0)+inner(rank, f1, eigenVectors+iii, 1)));
                
                f1->sinc.tulip[eigenVectors+iii].symmetry = tClassify(rank,f1, eigenVectors+iii );
            }
            
            INT_TYPE kkk;
           
           
//           for ( kkk = 0; kkk < imin(quantumBasisSize,Ne) ; kkk++)
//               f1->sinc.tulip[eigenVectors+kkk].value2 =0.;
//#ifdef OMP
//#pragma omp parallel for private (kkk,rank,cmpl) schedule(dynamic,1)
//#endif
//
//            for ( kkk = 0; kkk < imin(quantumBasisSize,Ne) ; kkk++)
//                for ( cmpl = 0; cmpl < 2 ; cmpl++)
//            {
//#ifdef OMP
//                rank = omp_get_thread_num();
//#else
//                rank = 0;
//#endif
//
//                tMultiplyMP(rank, &info, f1, 1., -1, copyVector, rank, 'N', density, 0, 'N', eigenVectors+kkk, cmpl);
//
//                f1->sinc.tulip[eigenVectors+kkk].value2 +=  tMultiplyMP(0, &info, f1, 1., -1, nullVector, 0, 'T', copyVector, rank, 'N', eigenVectors+kkk, cmpl);
//            }
                
                
                
                
                
             //   printf("\nClass%lld:%lld:\t%lld \t %1.15f\t %d \t %1.1f \t", kkk+1, f1->sinc.N1,kkk+1, f1->sinc.tulip[eigenVectors+kkk].value ,f1->sinc.tulip[eigenVectors+kkk].symmetry,deg(f1,f1->sinc.tulip[eigenVectors+kkk].symmetry));
                
                
                
                //                if ( f1->sinc.tulip[eigenVectors+iii].symmetry2 != f1->sinc.tulip[eigenVectors+iii].symmetry ){
                //                    printf("\nBuild%lld:\t%lld \t %f\t %d \t %1.1f \t %d \t %1.1f ", iii+1+iNe,iii+1+iNe, ritz[iii+iNe],
                //                           f1->sinc.tulip[eigenVectors+iii].symmetry, deg(f1,f1->sinc.tulip[eigenVectors+iii].symmetry),
                //                           f1->sinc.tulip[eigenVectors+iii].symmetry2, deg(f1,f1->sinc.tulip[eigenVectors+iii].symmetry2));
                //                    f1->sinc.tulip[eigenVectors+iii].symmetry = f1->sinc.tulip[eigenVectors+iii].symmetry2;
                //                }
                //fflush(stdout);
            }
        
    }
//    time(&lapse_t);
//    f1->mem1->rt->lanczosTime += difftime(lapse_t, start_t);
//    f1->sinc.tulip[matrixHbuild].Current[0] = quantumBasisSize;
//    f1->sinc.tulip[matrixSbuild].Current[0] = quantumBasisSize;

    return 0;
}



