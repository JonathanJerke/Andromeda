/*
 *  eigen.c
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


#include "eigen.h"

double vale ( struct sortClass * f ){
    INT_TYPE * mm = f->mmm+f->i*6 ;
    return f->str[0][mm[0] + mm[3]*f->n1[0]] +f->str[1][mm[1] + mm[4]*f->n1[1]] +f->str[2][mm[2] + mm[5]*f->n1[2]];
}

double yale ( struct sortClass * f ){
    INT_TYPE * mm = f->mmm+f->i*6 ;
    return mm[3]+mm[4]*f->nG + mm[5]*f->nG*f->nG;
}


int sortComp (const void * elem1, const void * elem2)
{
    struct sortClass* f = ((struct sortClass*)elem1);
    struct sortClass* s = ((struct sortClass*)elem2);
    double valueF,valueS;
    valueF = vale(f);
    valueS = vale(s);
    if (valueF > valueS) return  1;
    if (valueF < valueS) return -1;
    return 0;
}

int sort2Comp (const void * elem1, const void * elem2)
{
    struct sortClass* f = ((struct sortClass*)elem1);
    struct sortClass* s = ((struct sortClass*)elem2);
    INT_TYPE valueF,valueS;
    valueF = yale(f);
    valueS = yale(s);
    if (valueF > valueS) return  1;
    if (valueF < valueS) return -1;
    return 0;
}

INT_TYPE t1BodyConstruction(struct calculation * c1, enum division eigen){
    assignCores(&c1->i.c,1);
    struct field * f1 = &(c1->i.c);
    enum body bootBodies = c1->i.body;
    INT_TYPE N1 = f1->sinc.N1,rank;
    INT_TYPE space,N2 = N1*N1,i,ii,iii,iv,v;
    double ar[N2],w[N1];
    

    for ( space = 0 ;space < SPACE ; space++){
        cblas_dcopy(N2, streams(f1, kinetic,0,space)+space*N2, 1, ar, 1);
        
        
        if ( c1->i.springFlag ){
            if ( space == 0 &&  ( c1->rt.runFlag )%2 == 0 ){
                cblas_daxpy(N2, 1.0,streams(f1, harmonium,0,0)+0*N2, 1, ar, 1);
            }
            if ( space == 1 &&  ( c1->rt.runFlag/2 )%2 == 0 ){
                cblas_daxpy(N2, 1.0,streams(f1, harmonium,0,1)+1*N2, 1, ar, 1);
            }
            if ( space == 2 &&  ( c1->rt.runFlag/4 )%2 == 0 ){
                cblas_daxpy(N2, 1.0,streams(f1, harmonium,0,2)+2*N2, 1, ar, 1);
            }
        }
        
        
        tdsyev(0, f1, 'V', N1, ar, N1, w);
        if ( bootBodies == one ){
            for ( i = 0  ; i < N1 ; i++){
                streams(f1,foundationStructure,0,space)[i] = w[i]-w[0];
                streams(f1,foundationStructure,1,space)[i] = 1;

            }
            cblas_dcopy(N2 , ar,1,streams(f1, eigen,0,space),1);
        }else if ( bootBodies == two ){
            
            
#ifdef OMP
#pragma omp parallel for private (v,rank,i,ii) schedule(dynamic,1)
#endif
            
            for ( v = 0 ; v < N2 ; v++){
#ifdef OMP
                rank = omp_get_thread_num();
#else
                rank = 0;
#endif
                
                i = (v)%N1;
                ii = (v/N1)%N1;
                
                // if ( l2*l2 > i*i+ii*ii )
                {
                    
                    cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,0), 1);
                    cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,1), 1);
                    cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,2), 1);
                    
                    f1->sinc.tulip[diagonal1VectorA].Current[rank] = 1;
                    cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,0), 1);
                    cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,1), 1);
                    cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,2), 1);
                    
                    f1->sinc.tulip[diagonal1VectorB].Current[rank] = 1;
                    f1->sinc.tulip[diagonalVectorA].Current[rank] = 0;
                    tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonalVectorA, rank);
                    
                    cblas_dcopy(N2, streams(f1, diagonalVectorA,rank,space), 1, streams(f1,eigen,0,space)+v*N2, 1);
                    
                    
                }
            }
            
            
            
            
            
        }else if ( bootBodies == three ){
#ifdef OMP
#pragma omp parallel for private (v,rank,i,ii,iii) schedule(dynamic,1)
#endif
            
            for ( v = 0 ; v < N2*N1 ; v++){
#ifdef OMP
                rank = omp_get_thread_num();
#else
                rank = 0;
#endif
                
                i = (v)%N1;
                ii = (v/N1)%N1;
                iii = (v/(N2))%N1;
                
                
                cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,0), 1);
                cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,1), 1);
                cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,2), 1);
                
                f1->sinc.tulip[diagonal1VectorA].Current[rank] = 1;
                cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,0), 1);
                cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,1), 1);
                cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,2), 1);
                
                f1->sinc.tulip[diagonal1VectorB].Current[rank] = 1;
                f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
                tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB,rank, diagonal2VectorA, rank);
                
                cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,0), 1);
                cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,1), 1);
                cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,2), 1);
                
                f1->sinc.tulip[diagonal1VectorC].Current[rank] = 1;
                
                f1->sinc.tulip[diagonalVectorA].Current[rank] = 0;
                tOuterProductSu(f1, diagonal2VectorA, rank, diagonal1VectorC, rank, diagonalVectorA, rank);
                
                cblas_dcopy(N1*N2, streams(f1, diagonalVectorA,rank,space), 1, streams(f1,eigen,0,space)+v*N1*N2, 1);
            }
        }else if ( bootBodies == four ){
            
            
#ifdef OMP
#pragma omp parallel for private (v,rank,i,ii,iii,iv) schedule(dynamic,1)
#endif
            
            for ( v = 0 ; v < N2*N2 ; v++){
#ifdef OMP
                rank = omp_get_thread_num();
#else
                rank = 0;
#endif
                i = (v)%N1;
                ii = (v/N1)%N1;
                iii = (v/(N2))%N1;
                iv = (v/(N2*N1))%N1;
                
                
                cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,0), 1);
                cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,1), 1);
                cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,2), 1);
                
                f1->sinc.tulip[diagonal1VectorA].Current[rank] = 1;
                cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,0), 1);
                cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,1), 1);
                cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,2), 1);
                
                f1->sinc.tulip[diagonal1VectorB].Current[rank] = 1;
                f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
                tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
                
                
                cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,0), 1);
                cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,1), 1);
                cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,2), 1);
                
                f1->sinc.tulip[diagonal1VectorC].Current[rank] = 1;
                
                cblas_dcopy(N1, ar+iv*N1, 1, streams(f1,diagonal1VectorD,rank,0), 1);
                cblas_dcopy(N1, ar+iv*N1, 1, streams(f1,diagonal1VectorD,rank,1), 1);
                cblas_dcopy(N1, ar+iv*N1, 1, streams(f1,diagonal1VectorD,rank,2), 1);
                
                f1->sinc.tulip[diagonal1VectorD].Current[rank] = 1;
                
                f1->sinc.tulip[diagonal2VectorB].Current[rank] = 0;
                f1->sinc.tulip[diagonalVectorA].Current[rank] = 0;
                
                tOuterProductSu(f1, diagonal1VectorC, rank, diagonal1VectorD, rank, diagonal2VectorB, rank);
                tOuterProductSu(f1, diagonal2VectorA, rank, diagonal2VectorB, rank, diagonalVectorA, rank);
                
                cblas_dcopy(N2*N2, streams(f1, diagonalVectorA,rank,space), 1, streams(f1,eigen,0,space)+v*N2*N2, 1);
                
                
            }
            
        }
    }
    f1->sinc.tulip[eigen].Current[0] = 1;
    
    
    return 0;
}



INT_TYPE tNBodyConstruction (struct calculation * c1, enum division build ,enum division eigen){
    struct field * f1 = &(c1->i.c);
    enum body bootBodies = c1->i.body;
    INT_TYPE cmpl,space,matrixNumber = c1->i.decomposeRankMatrix;
    //THRESING FLOOR
    tClear(f1, eigen);
//    if ( c1->rt.printFlag  ){
//        tScale(f1, density,-1);
////        tEqua(f1, squareTwo, 0, density, 0);
////        tPermute(0, f1, 'b', squareTwo, 0, density, 0);
//        if (bodies ( f1,eigen ) == two)
//            tCycleDecompostionOneMP(0, f1, density, 0, spam, 0, f1->mem1->rt->CONVERGENCE, 1, -1);
//    }
//
//    else
    {
    
        for ( cmpl =  0; cmpl < spins(f1, eigenVectors) ; cmpl++ ) {
            tClear(f1, build);
            zero ( f1,eigen,cmpl);
            zero(f1,build,0);
            tClear(f1, copy);
            if ( cmpl == 0 ){
                tEqua(f1, copy ,0, kinetic,0);
                if ( c1->i.springFlag ){
                    tAddTw(f1,copy,0, harmonium ,0);
                }
                
                tSumMatrices(f1, build, copy);// B x (1+S) * 3
                tCycleDecompostionOneMP(0, f1, build, 0, eigen, 0, f1->mem1->rt->CONVERGENCE, 1, -1);
                
                if ( bootBodies > one ){
                    
                    
                    if ( c1->i.bodyType > electron ){
                        if ( CanonicalRank(f1, interactionExchangePlus, 0)){
                            tEqua(f1, squareTwo,0, ocean(0,f1, interactionExchangePlus,0,0),0);
                            tCycleDecompostionOneMP(0, f1, interactionExchangePlus, 0,squareTwo, 0, f1->mem1->rt->CONVERGENCE, 1 , -1);
                            tSumMatrices(f1, build, squareTwo);//B2:1 || B3 : 3
                            tCycleDecompostionOneMP(0, f1, build, 0, eigen, 0, 1e-3, 1, -1);
                        }
                    }
                    
                    if ( CanonicalRank(f1, interactionExchange, 0)){
                        tEqua(f1, squareTwo,0, ocean(0,f1, interactionExchange,0,0),0);
                        tCycleDecompostionOneMP(0, f1, interactionExchange, 0,squareTwo, 0, f1->mem1->rt->CONVERGENCE, 1 , -1);
                        tSumMatrices(f1, build, squareTwo);//B2:1 || B3 : 3
                        tCycleDecompostionOneMP(0, f1, build, 0, eigen, 0, 1e-3, 1, -1);
                    }
                }
//                if ( f1->Na && matrixNumber  ){
//                    tClear(f1, copy);
//                    tCycleDecompostionOneMP(0, f1, linear, 0,copy, 0, f1->mem1->rt->CONVERGENCE, matrixNumber, -1);
//                    tSumMatrices(f1, build, copy);
//                    tCycleDecompostionOneMP(0, f1, build, 0, eigen, 0, f1->mem1->rt->CONVERGENCE, 1, -1);
//                }
            }else {
                if ( CanonicalRank(f1, kinetic,1 ) ){
                    tEqua(f1, copy ,0, kinetic,1);
                    tSumMatrices(f1, build, copy);// B x (1+S) * 3
                    tCycleDecompostionOneMP(0, f1, build, 0, eigen, 1, f1->mem1->rt->CONVERGENCE, 1, -1);
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
#if 1
            printf ("check: space%d eigen %f < %f\n", space,streams(f1,foundationStructure,0,space)[0] ,streams(f1,foundationStructure,0,space)[n1[space]-1] );
#endif
            for ( i = 0; i < n2[space]; i++){
                streams(f1, eigen,0,space)[i] = creal(hmat[i]);
                streams(f1, eigen,1,space)[i] = cimag(hmat[i]);
            }
        }
        
        
        
//        ze[0] = 1e254;
//        ze[1] = 1e254;
//        ze[2] = 1e254;
//        xe[0] = -1e254;
//        xe[1] = -1e254;
//        xe[2] = -1e254;
        
        
//        for ( space = 0; space < SPACE ; space++){
//
//            if ( ze[space] > streams(f1,foundationStructure,0,space)[r*n1[space]])
//                ze[space] = streams(f1,foundationStructure,0,space)[r*n1[space]];
//            if ( xe[space] < streams(f1,foundationStructure,0,space)[(r+1)*n1[space]-1] )
//                xe[space] = streams(f1,foundationStructure,0,space)[(r+1)*n1[space]-1];
//
//            for ( i = 0 ; i < n1[space] ; i++){
//                if ( ze[space]*xe[space] > 0 && ze[space] < 0 )
//                {
//                    // printf("%lld %lld %lld %f \t",space,i,r*n1[space]+i,streams(f1,foundationStructure,0,space)[r*n1[space]+i]);
//
//                    streams(f1,foundationStructure,0,space)[r*n1[space]+i] = exp(-(streams(f1,foundationStructure,0,space)[r*n1[space]+i]-xe[space])/(ze[space]-xe[space]));
//                    // printf("%f \n",streams(f1,foundationStructure,0,space)[r*n1[space]+i]);
//
//
//                }else{
//
//                    streams(f1,foundationStructure,0,space)[r*n1[space]+i] = exp(-(streams(f1,foundationStructure,0,space)[r*n1[space]+i]-ze[space])/(xe[space]-ze[space]));
//
//                }
//            }
//        }
        
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

INT_TYPE tFoundationLevel( struct field * f1, enum division A , double lvlm, double lvlx,INT_TYPE ops,enum division build,INT_TYPE xB, double lvl1, double lvl2, double lvl3,INT_TYPE *mmm, INT_TYPE irrep,double seekPower){
    //GRID
    /// GRID
    ///// GRID
    ////////GRID
    INT_TYPE v,lvl,i,j,k,r1,r2,r3,space,ii,jj,kk,di,xx[3];
    INT_TYPE vaMax=0,sp, classicalBasisSize,*mm;
    INT_TYPE n1[3];
    enum body bd = f1->body;
    INT_TYPE nG = tSize(bd);
    n1[0] = vectorLen(f1,A)[0];
    n1[1] = vectorLen(f1,A)[1];
    n1[2] = vectorLen(f1,A)[2];
    double value;
    classicalBasisSize = 0;
    //for( lvl = 0; lvl < imin(n1[0],imin(n1[1],n1[2])) ; lvl++)
    {
        
        for ( space = 0; space < 3 ;space++)
            for ( i = 0; i < n1[space] ; i++)
                if ( streams(f1,foundationStructure,1,space)[i] > 0. ){
                    xx[space] = i+1;
                }
        
        
        
        ii = 0;
        for ( r1 = 0; r1 < nG ; r1++)
            for ( i = 0 ; i < xx[0] ; i++)
                if ( streams(f1,foundationStructure,1,0)[i+r1*n1[0]] > 0.)
                    
                    
                {
                    
                    
                    
                    //printf("%lld\n",ii);
                    jj = 0;
                    for ( r2 = 0; r2 < nG ; r2++)
                        for ( j = 0 ; j < xx[1] ; j++){
                            //  printf("%lld %lld %f %f\n", r2,j ,streams(f1,foundationStructure,0,1)[j+r2*n1[1]] , lvl2);
                            if ( streams(f1,foundationStructure,1,1)[j+r2*n1[1]] > 0. )
                            {
                                
                                
                                //    printf("%lld %lld\n",ii,jj);
                                
                                kk = 0;
                                for ( r3 = 0; r3 < nG ; r3++)
                                    for ( k = 0 ; k < xx[2] ; k++)
                                        //if ( i == lvl || j == lvl || k == lvl )
                                        
                                        if ( streams(f1,foundationStructure,1,2)[k+r3*n1[2]] > 0.)
                                        {
                                            ii = r1;
                                            jj = r2;
                                            kk = r3;
                                            value = (streams(f1,foundationStructure,0,0)[i+ii*n1[0]] +
                                                     streams(f1,foundationStructure,0,1)[j+jj*n1[1]] +
                                                     streams(f1,foundationStructure,0,2)[k+kk*n1[2]]);
                                            
                                            if ( vaMax < value )
                                                vaMax = value;
                                            
                                            if ( lvlm <= value ){
                                                {
                                                    //                                                                printf("**%f\n", value);
                                                    
                                                    
                                                    if ( tIR ( bd,ii,jj,kk,irrep)){
                                                
                                                        
                                                        if ( ops  == 0 ){
#if VERBOSE
                                                            printf("%d-->%d : %d %d %d :: %d %d %d :: %f %f %f\n",irrep,classicalBasisSize,i,j,k,ii+1,jj+1,kk+1, streams(f1,foundationStructure,0,0)[i+ii*n1[0]] ,
                                                                   streams(f1,foundationStructure,0,1)[j+jj*n1[1]] ,
                                                                   streams(f1,foundationStructure,0,2)[k+kk*n1[2]]);
#endif
                                                            mm = mmm+classicalBasisSize*6;
                                                            
                                                            
                                                            mm[0] = i;
                                                            mm[1] = j;
                                                            mm[2] = k;
                                                            
                                                            mm[3] = ii;
                                                            mm[4] = jj;
                                                            mm[5] = kk;
                                                        }
                                                        classicalBasisSize++;
                                                        
                                                    }
                                                }
                                            }
                                        }
                            }
                        }
                }
    }
    
    
    
    if ( ops == 2 )
        return vaMax;
    else
        return classicalBasisSize;
}



INT_TYPE tFilter(struct field * f1, INT_TYPE Ve, INT_TYPE irrep, enum division usr){
    INT_TYPE ii,j,cmpl=0,rank,flag ,sp,nP = tPerms(f1->body),nG = tSize(f1->body);
    double value;
    assignCores(f1, 2);
#ifdef OMP
#pragma omp parallel for private (ii,sp,rank,cmpl) schedule(dynamic,1)
#endif
    for ( ii = 0; ii < Ve ; ii++)
    {
#ifdef OMP
        rank = omp_get_thread_num();
#else
        rank = 0;
#endif
        printf("rank %d %d %d\n", rank,ii+1,CanonicalRank(f1, usr+ii, 0) );
        fflush(stdout);

//#ifdef splitTag
//        if ( irrep ){
//            for ( sp = 0; sp < spins(f1, usr+ii);sp++)
//                for ( r = 0 ; r < CanonicalRank(f1, usr+ii, sp); r++){
//                    zero(f1, diagonalVector, rank);
//                    f1->sinc.tulip[diagonalVector].Current[rank] = 1;
//                    f1->sinc.tulip[permutationVector].Current[rank] = 0;
//                    tBuildIrr(rank, f1, -1, ocean(rank, f1, usr+ii, r, sp), sp, permutationVector, rank);
//                    for ( r2 = 0; r2 < nP ; r2++){
//                        cblas_daxpy(n1[0], get1(f1->body, (tPath(f1,usr+ii))%nG+1, r2), streams(f1, permutationVector,rank,0)+r2*n1[0], 1, streams(f1,diagonalVector, rank, 0),1);
//                        cblas_daxpy(n1[1], get1(f1->body, (tPath(f1,usr+ii)/nG)%nG+1, r2), streams(f1, permutationVector,rank,1)+r2*n1[1], 1, streams(f1,diagonalVector, rank, 1),1);
//                        cblas_daxpy(n1[2], get1(f1->body, (tPath(f1,usr+ii)/(nG*nG))%nG+1, r2), streams(f1, permutationVector,rank,2)+r2*n1[2], 1, streams(f1,diagonalVector, rank, 2),1);
//                        tEqua(f1, ocean(rank, f1, usr+ii, r, sp), sp, diagonalVector, rank);
//                    }
//                }
//            printf(" %d<---%d (%d) = %d\n",ii+usr, CanonicalRank(f1,usr+ii,cmpl),tPath(f1,usr+ii),f1->sinc.tulip[usr+ii].symmetry);
//        }
//#else
        if ( irrep && bodies(f1, usr+ii ) > one  ){
            for ( cmpl = 0 ; cmpl < spins(f1, usr+ii) ; cmpl++){
                f1->sinc.tulip[permutationVector].Current[rank] = 0;
                tBuildIrr(rank, f1, irrep+nP, usr+ii, cmpl, permutationVector, rank);
                tCycleDecompostionOneMP(rank, f1,permutationVector, rank, usr+ii,cmpl,f1->mem1->rt->vCONVERGENCE, part(f1,usr+ii), 1);
            }
        }
        
//#endif
        f1->sinc.tulip[usr+ii].symmetry = tClassify(rank, f1, usr+ii);
      //  printf(" %d->%d (%d) = %d\n",ii+usr, CanonicalRank(f1,usr+ii,cmpl),tPath(f1,usr+ii),f1->sinc.tulip[usr+ii].symmetry);
      //  fflush(stdout);
    }
    
    return 0;
}

INT_TYPE tSelect(struct field * f1, INT_TYPE Ve, INT_TYPE type, enum division usr, enum division usa, INT_TYPE testFlag){
    INT_TYPE sp,info,rank=0,maxEV = f1->sinc.maxEV,n,m;
    INT_TYPE stride = maxEV;
    double value;
    
	if ( ! CanonicalRank(f1, usa,0))
        return 0;


    if ( type ) {
        tClear(f1, usr+Ve);
        for ( sp = 0; sp < spins(f1, usr+Ve);sp++){
            f1->sinc.tulip[permutationVector].Current[rank] = 0;
            tBuildIrr(rank, f1, type, usa, sp, permutationVector, rank);
            tCycleDecompostionOneMP(rank, f1, permutationVector, rank, usr+Ve, sp, f1->mem1->rt->vCONVERGENCE, part(f1,usr+Ve), -1);
        }
        
        value = magnitude(f1, usr+Ve);
        
        
        if ( isnan(1./value) || isinf(1./value) )
           return 0;
        tScale(f1, usr+Ve, 1./value);

    } else {
	
        value = magnitude(f1, usa);
        
        
       if ( value < 1e-6 || isnan(value) || isinf(value) )
            return 0;
//        printf ( "tsel:%f\n", value);

        
        if ( testFlag >= 0 ){
            tEquals(f1, usr+Ve, usa);
            
            tScale(f1, usr+Ve, 1./value);
            
        }
    }
    if ( testFlag != 1 )
        return 1;

    
    DCOMPLEX * S = (DCOMPLEX*)(myStreams(f1, matrixSbuild, 0));

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
            tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usr+n,0,'N', usr+m,0);

             if ( spins(f1, usr+n) > 1 && (CanonicalRank(f1, usr+n, 1) ||CanonicalRank(f1, usr+m, 1))){
                 S[n*stride+m] +=
                 
                 +tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usr+n,1,'N', usr+m,1)
                 +I*tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usr+n,1,'N', usr+m,0)
                 -I*tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, CDT, usr+n,0,'N', usr+m,1)
                 ;
             }
            S[m*stride+n] = conj(S[n*stride+m]);
        }
        
    }
    if (0)
        for ( n = 0; n < Ve+1 ; n++){
            m = n ;
            //for ( m = 0; m < Ve+1 ; m++)
            {
                if ( cabs(S[n*stride+m]-1.0)>1e-3 || isnan(S[n*stride+m] || isinf(S[n*stride+m])))
                    //            {
                    printf("%lld %d :%f +i%f \n", n,m,creal(S[n*stride+m]), cimag(S[n*stride+m]) );
                fflush(stdout);
                //            }
            }
        }

    
    assignCores(f1, 0);
  //  cblas_zcopy(maxEV*stride,S,1,S+maxEV*stride,1);
    tzheev(0, f1, 'N', Ve+1,S/*+maxEV*stride*/, stride, ov);
    if ( testFlag ){

        if ( ov[Ve]/ov[0] < f1->mem1->rt->TOL && ov[Ve]/ov[0] >0 && ov[0] > 0 ){
                printf("(%f,%f,%d)\n",  ov[Ve],ov[0],Ve+1 );
                fflush(stdout);

                return 1;
        } else {
            //printf("* %1.15f\n", ov[Ve]/ov[0]);
        }
    }
    
    return 0;
}


INT_TYPE tCollect (struct field * f1, INT_TYPE irrep,enum division usz, INT_TYPE target,double seekPower){
    INT_TYPE  Ve = 0;
    f1->sinc.tulip[diagonalVectorA].header = Cube;
    f1->sinc.tulip[diagonalVectorB].header = Cube;
    INT_TYPE nG = tSize(f1->body);
    printf("Collect %d\nN1 %d\nBody %d\n Target %d\n\n", target, f1->sinc.N1, f1->body,target);
    if ( 1 ){
        INT_TYPE flag = 0,ct = 0;
        double min0= 0,max0= tFoundationLevel(f1, build,0,0,2,0,target,1e9,1e9,1e9,NULL,irrep,seekPower) ;
        double min = min0, max =  max0,vx= 1000.,va = 1.2;
        
        while(1){
            if ( flag ){
                if ( ct < target*va )
                    min = 0.5*(max+min);
                else
                    max = 0.5*(max+min);
            }
            ct = tFoundationLevel(f1, build,0,0.5*(min+max),1,usz,target,1e9,1e9,1e9,NULL,irrep,seekPower);
            printf("va %f %d %d %d %f\n", 0.5*(min+max),ct ,flag++, target,va);
            if ( max-min < 1e-9  ){
                printf("conv");
                exit(0);
            }
            if ( ct >= va*target && ct <= target*va*vx)
            {
            
                Ve = 0;
                INT_TYPE mmm[6*ct],i,*mm;
                ct = tFoundationLevel(f1, build,0,0.5*(min+max),0,usz,target,1e9,1e9,1e9,mmm,irrep,seekPower);
                INT_TYPE *n1 = vectorLen(f1, build);
                struct sortClass sc [ct];
                for ( i = 0; i < ct ; i++){
                    sc[i].n1[0] = n1[0];
                    sc[i].n1[1] = n1[1];
                    sc[i].n1[2] = n1[2];
                    sc[i].i = i;
                    sc[i].nG = nG;
                    sc[i].mmm = mmm;
                    sc[i].str[0] = streams(f1,foundationStructure,0,0);
                    sc[i].str[1] = streams(f1,foundationStructure,0,1);
                    sc[i].str[2] = streams(f1,foundationStructure,0,2);
                }
                
                qsort(sc, ct, sizeof(struct sortClass), &sortComp);
                
                
                for ( i = 0; i < ct ; i++){
                    mm = mmm+sc[i].i*6;
                    tClear(f1,diagonalVectorA);
                    f1->sinc.tulip[diagonalVectorA].Current[0] = 1;;
                    
                    cblas_dcopy(n1[0], streams(f1, build, 0 , 0)+(mm[0])*n1[0]+(mm[3])*n1[0]*n1[0],1,streams(f1,diagonalVectorA,0,0),1);
                    cblas_dcopy(n1[1], streams(f1, build, 0 , 1)+(mm[1])*n1[1]+(mm[4])*n1[1]*n1[1],1,streams(f1,diagonalVectorA,0,1),1);
                    cblas_dcopy(n1[2], streams(f1, build, 0 , 2)+(mm[2])*n1[2]+(mm[5])*n1[2]*n1[2],1,streams(f1,diagonalVectorA,0,2),1);
             
#if VERBOSE
                    printf("%d:: %d %d %d :: %d %d %d %f %f\n", i,mm[0],mm[1],mm[2],mm[3]+1,mm[4]+1,mm[5]+1,vale(&sc[i]),magnitude(f1, diagonalVectorA) );
#endif
                
                    Ve =  tSASplit(f1, irrep, Ve,target, usz, diagonalVectorA);
                    if ( Ve == target )
                        break;
                }
#ifdef splitTag
                if ( Ve == target ){
                    Ve = 0;
                    ct = i+3;
                    qsort(sc, ct, sizeof(struct sortClass), &sort2Comp);
                    
                    for ( i = 0; i < ct ; i++){
                        mm = mmm+sc[i].i*6;
                        tClear(f1,diagonalVectorA);
                        f1->sinc.tulip[diagonalVectorA].Current[0] = 1;;
                        
                        cblas_dcopy(n1[0], streams(f1, build, 0 , 0)+(mm[0])*n1[0]+(mm[3])*n1[0]*n1[0],1,streams(f1,diagonalVectorA,0,0),1);
                        cblas_dcopy(n1[1], streams(f1, build, 0 , 1)+(mm[1])*n1[1]+(mm[4])*n1[1]*n1[1],1,streams(f1,diagonalVectorA,0,1),1);
                        cblas_dcopy(n1[2], streams(f1, build, 0 , 2)+(mm[2])*n1[2]+(mm[5])*n1[2]*n1[2],1,streams(f1,diagonalVectorA,0,2),1);
                        
                        Ve2 =  tSASplit(f1, irrep, Ve,target, usz, diagonalVectorA);
                        for ( j = Ve ; j < Ve2 ; j++){
                            printf("%d : p%d\n", usz+j, mm[5]*nG*nG + mm[4]*nG + mm[3]);
                            f1->sinc.tulip[usz+j].path = mm[5]*nG*nG + mm[4]*nG + mm[3];
                        }
                        Ve = Ve2;
                        if ( Ve == target )
                            break;

                    }
                    
                }
#endif
                if ( Ve == target )
                    break;
            }            
        }
    if ( target != Ve ){
        printf("ack no !\n");
        exit(0);
    }
    }
    
    return Ve;
}


INT_TYPE tSASplit ( struct field * f1, INT_TYPE irrep , INT_TYPE Ve , INT_TYPE target,enum division usz, enum division vector){
    INT_TYPE map[24],nDeg=0,ii;
    
    if ( f1->body == one || irrep == 0 ){
        nDeg = 1;
        map[1] = 0;
    }else
        
        if ( f1->body == two ){
            if ( irrep == 1 ){
                map[1] = 1;
                nDeg = 1;
            }else
                if ( irrep == 2 ){
                    nDeg = 1;
                    map[1] = 2;
                }
        }else
            if ( f1->body== three ){
                if ( irrep == 1 ){
                    map[1] = 1;
                    nDeg = 1;
                }else
                    if ( irrep == 2 ){
                        map[1] = 2;
                        nDeg = 1;
                    } else {
                        nDeg = 4;
                        map[1] = 3;
                        map[2] = 4;
                        map[3] = 5;
                        map[4] = 6;
                    }
            }
            else if ( f1->body == four ){
                //
                if ( irrep == 1 ){
                    nDeg = 1;
                    map [1] = 1;
                }else if ( irrep == 2 ){
                    map[1] = 2;
                    nDeg = 1;
                } else if ( irrep == 3 ){
                    map[1] = 3;
                    map[2] = 4;
                    map[3] = 5;
                    map[4] = 6;
                    nDeg = 4;
                } else if ( irrep == 4){
                    map[1] = 7;
                    map[2] = 8;
                    map[3] = 9;
                    map[4] = 10;
                    map[5] = 11;
                    map[6] = 12;
                    map[7] = 13;
                    map[8] = 14;
                    map[9] = 15;
                    nDeg = 9;
                }else if ( irrep == 5 ){
                    map[1] = 16;
                    map[2] = 17;
                    map[3] = 18;
                    map[4] = 19;
                    map[5] = 20;
                    map[6] = 21;
                    map[7] = 22;
                    map[8] = 23;
                    map[9] = 24;
                    nDeg = 9;
                }
            }
    
    
    for ( ii = 1 ; ii <= nDeg ; ii++)
        if ( tSelect(f1, Ve, map[ii] , usz, vector, 1)){
            Ve++;
            if ( Ve == target )
                return Ve;
        
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


INT_TYPE tGreatDivideIteration ( struct field * f1, enum division A , INT_TYPE I1, INT_TYPE I2, enum division usz, INT_TYPE foundation,INT_TYPE nMult, INT_TYPE shift){
    INT_TYPE expon,info;
    INT_TYPE rank ,i;
    //time_t start_t, lapse_t;
    double temp,temp2, sum = 0;
    //time(&start_t);
    INT_TYPE iii = 0;
    assignCores(f1, 1);

    for( expon = 1 ; foundation*expon < nMult  ; expon++){
        
#ifdef OMP
#pragma omp parallel for private (iii,rank,temp,temp2) schedule(dynamic,1) reduction(+:sum)
#endif
        for ( iii = 0; iii < foundation ; iii++)
        {
            
#ifdef OMP
            rank = omp_get_thread_num();
#else
            rank = 0;
#endif
            
            
            {
                //printf("%d: %d %d %d %1.15f\n", rank,iii+1,part(f1,usz+(expon)*foundation+iii),part(f1,usz+(expon-1)*foundation+iii) ,f1->mem1->rt->vCANON);
                tEquals(f1,usz+iii+expon*foundation,usz+(expon-1)*foundation+iii );
                tHXpX(rank, f1, A, 0, 1.0, 0.0, usz+iii+expon*foundation, f1->mem1->rt->vCANON , part(f1,usz+(expon)*foundation+iii));
           //	printf("usr %d : %d %d\n", usz+iii+expon*foundation, CanonicalRank(f1, usz+expon*foundation+iii,0),part(f1, usz+foundation*expon+iii));     
                temp = inner(rank, f1, usz+(expon)*foundation+iii, 0)+inner(rank, f1, usz+(expon)*foundation+iii, 1);
                temp2 = temp- sqr( tMultiplyMP(rank, &info, f1, 1.0, -1, nullVector, 0, 'T', usz+(expon)*foundation+iii, 0, 'N', usz+(expon-1)*foundation+iii, 0) + tMultiplyMP(rank, &info, f1, 1.0, -1, nullVector, 0, 'T', usz+(expon)*foundation+iii, 1, 'N', usz+(expon-1)*foundation+iii, 1));
#if VERBOSE
                printf("%d(%d) -> %d(%d)::[%f]\n",(expon-1)*foundation+iii,part(f1,usz+(expon-1)*foundation+iii ), expon*foundation+iii,part(f1,usz+(expon)*foundation+iii ),temp2  );
#endif
                sum += temp2;
                tScale(f1, usz+iii+expon*foundation, 1./sqrt(temp));

//                if ( expon == 1 )
//                    ev[foundation] = tMultiplyMP(rank, &info,f1, 1., -1, nullVector, 0, 'T', usz+iii+expon*foundation, 0, 'N',usz+(expon-1)*foundation+iii, 0)+tMultiplyMP(rank,&info, f1, 1., -1, nullVector, 0, 'T', usz+iii+expon*foundation, 1, 'N',usz+(expon-1)*foundation+iii, 1);
//
                //fflush(stdout);
          //      ev[iii ]= magnitude(f1,usz+iii+expon*foundation);
#if 1
                printf("%d :: \t %f\t %d\n",iii,temp2, CanonicalRank(f1, usz+(expon)*foundation+iii, 0)+CanonicalRank(f1, usz+(expon)*foundation+iii, 1));
#endif
                f1->sinc.tulip[usz+(expon)*foundation+iii].value3 = temp2;
                fflush(stdout);
            }
            
        }
    }
    
    
    
    
    
    INT_TYPE sum2 = 0;
    expon = 1;
    for ( iii = 0; iii < foundation ; iii++)
        sum2 += CanonicalRank(f1, usz+(expon)*foundation+iii, 0)+CanonicalRank(f1, usz+(expon)*foundation+iii, 1);
    
    printf("GREATER:\t %f \t %f\n", sum2*1. / foundation , sum );
    
    
    
    return nMult;
}

INT_TYPE tMinorDivideIteration ( struct field * f1, enum division A , INT_TYPE I1, INT_TYPE I2, enum division usz, INT_TYPE foundation,INT_TYPE nMult, double shift){
    INT_TYPE expon,info;
    INT_TYPE rank ,i;
    INT_TYPE iii = 0;
    assignCores(f1, 1);
    
    for( expon = 1 ; foundation*expon < nMult  ; expon++){
        
#ifdef OMP
#pragma omp parallel for private (iii,rank) schedule(dynamic,1)
#endif
        for ( iii = 0; iii < foundation ; iii++)
        {
            
#ifdef OMP
            rank = omp_get_thread_num();
#else
            rank = 0;
#endif
            
            
            {
                if ( CanonicalRank(f1, usz+(expon-1)*foundation+iii, 0))
                    tEqua(f1,usz+iii+expon*foundation,0,ocean(rank, f1, usz+(expon-1)*foundation+iii, rand()%CanonicalRank(f1, usz+(expon-1)*foundation+iii, 0) , 0),0 );
                if ( CanonicalRank(f1, usz+(expon-1)*foundation+iii, 1))
                    tEqua(f1,usz+iii+expon*foundation,1,ocean(rank, f1, usz+(expon-1)*foundation+iii, rand()%CanonicalRank(f1, usz+(expon-1)*foundation+iii, 1) , 1),1 );

                tHXpX(rank, f1, A, 0, 1.,0., usz+iii+expon*foundation, f1->mem1->rt->vCANON , part(f1,usz+iii+expon*foundation)/shift);
            }
            
        }
    }

    return nMult;
}



INT_TYPE tSAboot(struct calculation *c1){
    struct field *f1 = &c1->i.c;
    INT_TYPE irrep;
    {
        
        assignCores(f1, 2);
        
        
        
//        if (c1->i.body != one ){
//            tNBodyConstruction ( c1, build,  eigen);
//        }
//        else
        {
            t1BodyConstruction ( c1, eigen);
        }
        
        INT_TYPE *n1 =vectorLen(f1, eigen);
        INT_TYPE   xyz1,i1,j1,k1,p, xyz = imin(n1[0],pow(c1->i.side,bodies(f1,eigenVectors)));
        INT_TYPE i,j,k,nG = tSize(f1->body), v,rr,r2,info,r,rank=0,space,irrepm,LN2[3];
        INT_TYPE hits[nG];
        for ( i = 0 ; i < nG ; i++)
            hits[i] = 0;
        length(f1, eigen, LN2);
        double vo,va,vm;
        for ( irrepm = 0 ; irrepm <  nG ; irrepm++){
            
            for ( r = 0 ; r < n1[0]; r++){
                streams(f1,foundationStructure,0,0)[r+(irrepm)*n1[0]] = 0;
                streams(f1,foundationStructure,0,1)[r+(irrepm)*n1[1]] = 0;
                streams(f1,foundationStructure,0,2)[r+(irrepm)*n1[2]] = 0;

                streams(f1,foundationStructure,1,0)[r+(irrepm)*n1[0]] = 0;
                streams(f1,foundationStructure,1,1)[r+(irrepm)*n1[1]] = 0;
                streams(f1,foundationStructure,1,2)[r+(irrepm)*n1[2]] = 0;
            }
        }
        
        
            
        for ( xyz1 = 1 ; xyz1 <= xyz ; xyz1++){
            
#ifdef OMP
#pragma omp parallel for private (r,r2,p,i1,j1,k1,rank,space,va,vo,i,j,k,irrepm) schedule(static,1)
#endif
            for ( r = 0 ; r < xyz1*xyz1; r++)
                
            {
#ifdef OMP
                rank = omp_get_thread_num();
#else
                rank = 0;
#endif
                
                i1 = (r)%xyz1;
                j1 = (r/xyz1)%xyz1;
                k1 = xyz1-1;
                
                if ( i1 <= j1 && j1 <= k1 )
                    if ( i1*i1 + j1*j1 + k1*k1 < c1->i.l2*c1->i.l2 )
                        for ( p = 0 ; p < 6 ; p++){
                            //                            if ( ! p )
                            //                                printf ("--->%d %d:: %d %d %d\n", xyz,r,i1,j1,k1);
                            
                            i = i1;
                            j = j1;
                            k = k1;
                            
                            if ( p == 0 ){
                                i = i1;
                                j = j1;
                                k = k1;
                            } else
                                if ( i1 != j1 || j1 != k1 || i1 != k1 ){
                                    
                                    if ( p == 1 ){
                                        i = j1;
                                        j = k1;
                                        k = i1;
                                    }else if ( p == 2 ){
                                        i = k1;
                                        j = i1;
                                        k = j1;
                                    }else if ( i1 != j1 && j1 != k1  && i1 != k1 ){
                                        if ( p == 3 ){
                                            i = j1;
                                            j = i1;
                                            k = k1;
                                        }else if ( p == 4 ){
                                            i = k1;
                                            j = j1;
                                            k = i1;
                                        }else if ( p == 5 ){
                                            i = i1;
                                            j = k1;
                                            k = j1;
                                        }
                                    }else {
                                        continue;
                                    }
                                }else {
                                    continue;
                                }
                           // printf("%d %d %d\n", i,j,k);
                            f1->sinc.tulip[diagonalVectorA].Current[rank] = 1;
                            space = 0;
                            cblas_dcopy(n1[space], streams(f1,eigen,0,space)+i*n1[space], 1, streams(f1, diagonalVectorA,rank,space), 1);
                            space = 1;
                            cblas_dcopy(n1[space], streams(f1,eigen,0,space)+j*n1[space], 1, streams(f1, diagonalVectorA,rank,space), 1);
                            space = 2;
                            cblas_dcopy(n1[space], streams(f1,eigen,0,space)+k*n1[space], 1, streams(f1, diagonalVectorA,rank,space), 1);
                            
                            f1->sinc.tulip[permutationVector].Current[rank]= 0;
                            
                            tBuildIrr(rank, f1, -1, diagonalVectorA, rank, permutationVector, rank);
                            
                            
                            zero(f1, diagonalVectorA, rank);
                            for ( irrepm = 0; irrepm < nG ; irrepm++){
                                for ( r2 = 0; r2 < CanonicalRank(f1, permutationVector, rank); r2++){
                                    cblas_daxpy(n1[0], get1(f1->body, irrepm+1, r2), streams(f1, permutationVector,rank,0)+r2*n1[0], 1, streams(f1,diagonalVectorA, rank, 0),1);
                                    cblas_daxpy(n1[1], get1(f1->body, irrepm+1, r2), streams(f1, permutationVector,rank,1)+r2*n1[1], 1, streams(f1,diagonalVectorA, rank, 1),1);
                                    cblas_daxpy(n1[2], get1(f1->body, irrepm+1, r2), streams(f1, permutationVector,rank,2)+r2*n1[2], 1, streams(f1,diagonalVectorA, rank, 2),1);
                                    //YES< THIS IS CORRECT<  the components are separate!!!
                                    //zero?
                                }
                                
                                vo = inner(rank,f1, diagonalVectorA,rank);
                           //     printf( "%d %d %f\n", irrepm+1,CanonicalRank(f1, diagonalVectorA,rank),vo);
                                
                                if ( vo > 0.01 && hits[irrepm] < xyz){
                                    tScaleOne(f1, diagonalVectorA,rank, 1./sqrt(vo));
                                    
                                    cblas_dcopy(n1[0], streams(f1, diagonalVectorA,rank,0), 1, streams(f1,build,0,0)+i*n1[0]+(irrepm)*LN2[0], 1);
                                    cblas_dcopy(n1[1], streams(f1, diagonalVectorA,rank,1), 1, streams(f1,build,0,1)+j*n1[1]+(irrepm)*LN2[1], 1);
                                    cblas_dcopy(n1[2], streams(f1, diagonalVectorA,rank,2), 1, streams(f1,build,0,2)+k*n1[2]+(irrepm)*LN2[2], 1);
                                    tEqua(f1, diagonalVector,rank, diagonalVectorA,rank);
                                    tHXpX(rank, f1, Ha, 0, 1.0, 0, diagonalVector, 1e-6, 2);
                                    va = tMultiplyMP(rank, &info, f1, 1., -1, nullVector, 0, 'T', diagonalVectorA, rank, 'N', diagonalVector, rank);
                                    space = 0;
                                    streams(f1,foundationStructure,0,space)[i+(irrepm)*n1[space]] += va;
                                    streams(f1,foundationStructure,1,space)[i+(irrepm)*n1[space]] += 1.;
                                    
                                    space = 1;
                                    streams(f1,foundationStructure,0,space)[j+(irrepm)*n1[space]] += va;
                                    streams(f1,foundationStructure,1,space)[j+(irrepm)*n1[space]] += 1.;
                                    
                                    space = 2;
                                    streams(f1,foundationStructure,0,space)[k+(irrepm)*n1[space]] += va;
                                    streams(f1,foundationStructure,1,space)[k+(irrepm)*n1[space]] += 1.;
//                                    printf( "%d : %d : %d %d %d %f\n", hits[irrepm],irrepm+1, i,j,k,va);
                                    hits[irrepm]++;
                                }
                            }
                        }
            }
            
        }
            
            for ( i = 0 ; i < xyz ; i++)
                for ( irrepm = 0 ; irrepm < nG ; irrepm++)
                for ( space = 0; space < SPACE ; space++)
                    if ( streams(f1,foundationStructure,1,space)[i+(irrepm)*n1[space]] > 0. )
                        streams(f1,foundationStructure,0,space)[i+(irrepm)*n1[space]] /= streams(f1,foundationStructure,1,space)[i+(irrepm)*n1[space]];
            
        
//        for ( space = 0; space < SPACE ; space++)
//            for ( irrepm = 0 ; irrepm < nG ; irrepm++){
//                rr = 0;
//                for ( r = 0 ; r < n1[space]; r++){
//                    if ( streams(f1,foundationStructure,1,space)[r+(irrepm)*n1[space]] > 0. ){
//                        if( r != rr ){
//                            printf ( "s%d t%d   : %d ->  %d\n",space, (irrepm),r,rr);
//                            cblas_dcopy(n1[space], streams(f1,build,0,space)+r*n1[space]+(irrepm)*LN2[space], 1, streams(f1,build,0,space)+rr*n1[space]+(irrepm)*LN2[space],1);
//                            streams(f1,foundationStructure,0,space)[rr+(irrepm)*n1[space]] = streams(f1,foundationStructure,0,space)[r+(irrepm)*n1[space]];
//                            streams(f1,foundationStructure,1,space)[rr+(irrepm)*n1[space]] = streams(f1,foundationStructure,1,space)[r+(irrepm)*n1[space]];
//
//                        }
//                        rr++;
//                    }
//                }
//                for ( r = rr ; r < n1[space]; r++){
//                    streams(f1,foundationStructure,1,space)[r+(irrepm)*n1[space]] = 0;
//                }
//
//            }

        {
            
            for ( space = 0; space < SPACE ; space++){
                va = INFINITY;
                for ( irrepm = 0 ; irrepm < nG ; irrepm++)
                    for ( v = 0 ; v < n1[space] ; v++){
                        if ( va > streams(f1,foundationStructure,0,space)[v+(irrepm)*n1[space]]&&
                            streams(f1,foundationStructure,1,space)[v+(irrepm)*n1[space]]> 0.)
                            va = streams(f1,foundationStructure,0,space)[v+(irrepm)*n1[space]];
                    }
                printf("-min->%f\n",va);
                for ( irrepm = 0; irrepm < nG ; irrepm++)
                    for ( v = 0 ; v < n1[space] ; v++)
                        if (streams(f1,foundationStructure,1,space)[v+(irrepm)*n1[space]]> 0.)
                        {
                            printf("s%d t%d %d %f\n",space, irrepm+1,v, streams(f1,foundationStructure,0,space)[v+(irrepm)*n1[space]]);
                            streams(f1,foundationStructure,0,space)[v+(irrepm)*n1[space]] -= va;

                        }
            }
            
        }
    }
    return 0;
}



INT_TYPE tEdges(struct calculation *c1){
    struct field * f1 = &c1->i.c;
    INT_TYPE info;
    enum body bootBodies = f1->body;
    if ( 1 ){
        //EDGES ALT
        enum block b,bx;
        INT_TYPE iii,jjj=1,dim,irrep;
        double sum = 0;
        for ( irrep = 1 ;irrep <= 5 ; irrep++){
            jjj=1;
            for ( iii = 0 ; iii < c1->i.nStates ; iii++)
                if ( f1->sinc.tulip[eigenVectors+iii].symmetry  == irrep&& (! c1->i.irrep || irrep == c1->i.irrep))
                {
                    //            tTransformFrom(f1, eigenVectors+iii, squareVector, basis);
                    bx = tv1;
                    if ( bootBodies == two )
                        bx = tv2;
                    else if ( bootBodies == three )
                        bx = tv3;
                    else if ( bootBodies == four )
                        bx = tv4;
                    printf("Edge%d:\t%d:%d:\t",jjj++,irrep,f1->sinc.N1);
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



INT_TYPE tEigenCycle (struct field * f1, enum division A ,char permutation,  INT_TYPE Ne, enum division usz, INT_TYPE quantumBasisSize ,INT_TYPE iterations, INT_TYPE foundation, INT_TYPE irrep,INT_TYPE flag,  enum division outputSpace, enum division outputValues){
    INT_TYPE gvOut,prevBuild;
    time_t start_t, lapse_t;
    
    time(&start_t);
    INT_TYPE countLam = 0,countTot = 0,cl;
    enum division Mat;
    INT_TYPE cmpl,cmpl2,cmpl3,cat,iii = 0,maxEV = f1->sinc.maxEV,rank;
    INT_TYPE stride = maxEV;
    double * ritz = myStreams(f1, outputValues, 0);
    double * overlap = myStreams(f1, conditionOverlapNumbers, 0);
    enum division el ;
    DCOMPLEX co ;
    DCOMPLEX *T  =  (DCOMPLEX *) myStreams(f1, matrixHbuild,0/*CORE RANK*/);
    double *H  =   myStreams(f1, vectorHbuild,0/*CORE RANK*/);
    double *vectors[2];
    vectors[0] =   myStreams(f1, matrixHbuild,0/*CORE RANK*/)+4*maxEV*stride;
    vectors[1] =   myStreams(f1, matrixHbuild,0/*CORE RANK*/)+5*maxEV*stride;
    DCOMPLEX *S  =  (DCOMPLEX *) myStreams(f1, matrixSbuild,0/*CORE RANK*/);
    INT_TYPE powerMat;

    INT_TYPE info,n,m,s1,g,r,rr,rx,gx,a,a2;
    enum division leftP ;
    prevBuild = 0* CanonicalRank(f1, matrixHbuild, 0);
    tClear  (f1, copyTwoVector);
    for ( n = prevBuild; n < quantumBasisSize ; n++)
    {
#if VERBOSE
        printf("m%d %f\n", usz+n,magnitude(f1, usz+n) );
#endif
    //    tScale(f1, usz+n, 1./magnitude(f1, usz+n));
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
                    (tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, permutation, usz+n,0,'N', usz+m,0))/H[n]/H[m];;
                    if ( spins(f1, usz+n) > 1 && (CanonicalRank(f1, usz+n, 1) ||CanonicalRank(f1, usz+m, 1)) ){
                        S[n*stride+m] +=
                        (tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, permutation, usz+n,1,'N', usz+m,1)
                         +I*tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, permutation, usz+n,1,'N', usz+m,0)
                         -I*tMultiplyMP(rank,&info,f1,1., -1,nullVector, 0, permutation, usz+n,0,'N', usz+m,1)
                         )/H[n]/H[m];;
                    }
                    S[m*stride+n] = conj(S[n*stride+m]);
                    {
                        T[n*stride+m] = 0.;
                        T[m*stride+n] = (T[n*stride+m]);
                    }
                }
            
        }
    if (0)
        for ( n = 0; n < quantumBasisSize ; n++){
            m = n ;
           // for ( m = 0; m < quantumBasisSize ; m++)
            {
                if ( cabs(S[n*stride+m]-1.0)>1e-3 || isnan(S[n*stride+m] || isinf(S[n*stride+m])))
                    //            {
                    printf("%lld %d :%f +i%f \n", n,m,creal(S[n*stride+m]), cimag(S[n*stride+m]) );
                fflush(stdout);
                //            }
            }
        }
    
            
    
            for ( s1 = 0 ; s1 < 1; s1++){
                leftP = A;
                cat = 0;

                do {

                    if ( linear == name(f1,leftP) ){
                        Mat = rivers(0, f1, linear, cat );
                        f1->sinc.tulip[Mat].blockType = f1->sinc.tulip[leftP].blockType;
                    }else {
                        Mat = leftP;
                    }
                    for ( cmpl = 0; cmpl < spins(f1,Mat)  ;cmpl++)
                        if ( CanonicalRank(f1, name(f1,Mat), cmpl)){

                        for ( cmpl2 = 0 ; cmpl2 < spins(f1,usz)  ;cmpl2++)
                                if ( CanonicalRank(f1,  Mat, cmpl)){
                                    co = 1.;
                                    if ( cmpl )
                                        co *= I;
                                    if ( cmpl2 )
                                        co *= -I;
#if 1
                                    printf("**%d-%d-%d-%d\t%d\t%d\t%d-%d\t %d %f\n",Mat,name(f1,Mat),cmpl,cmpl2, bodies(f1, Mat),f1->sinc.tulip[Mat].blockType,CanonicalRank(f1, Mat, cmpl),f1->sinc.tulip[Mat].ptRank[cmpl], cat,traceOne(f1, Mat, cmpl) );
                                    fflush(stdout);
#endif
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
                                            (T+stride*maxEV)[n*stride+m] = co/H[m]/H[n]*matrixElements(rank, f1, permutation,usz+n, 'N', Mat, cmpl, usz+m, cmpl2);
                                        }
                                        
                                        
                                        
                                    }

                        for ( n = prevBuild; n < quantumBasisSize ; n++)
                            for ( m = 0; m <= n ; m++){
                                T[n*stride+m] += (T+stride*maxEV)[n*stride+m];
                                T[m*stride+n]  = conj(T[n*stride+m]);
                            }

                    }
                    for ( n = 0; n < quantumBasisSize ; n++)
                        for ( m = 0 ; m <=n   ; m++){
                            if ( cabs((T+stride*maxEV)[n*stride+m])  < 1e-6 ){
                                countLam++;
                            }
                            countTot++;
                        }

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

        cblas_zcopy(maxEV*stride , T , 1 , T+maxEV*stride , 1);
        cblas_zcopy(maxEV*stride , S , 1 , S+maxEV*stride , 1);
        tzheev(0, f1, 'N', quantumBasisSize, S+maxEV*stride, stride, overlap);
        //printf("%f\n", f1->mem1->rt->TOL);
        if ( ! flag )
            printf("GREATER: \t %d \t %f \n",quantumBasisSize,overlap[quantumBasisSize-1]/overlap[0] );

        if (  (overlap[0] > 0.) && overlap[quantumBasisSize-1]/overlap[0] < f1->mem1->rt->TOL && flag <= 2)
            printf(" Krylov-2\t %f \n",  overlap[quantumBasisSize-1]/overlap[0]);
        else         if (  (overlap[0] > 0.) && overlap[quantumBasisSize-1]/overlap[0] < f1->mem1->rt->TOL && flag == 3)
            printf(" Krylov-3 \t %f \n",  overlap[quantumBasisSize-1]/overlap[0]);
        else         if (  (overlap[0] > 0.) && overlap[quantumBasisSize-1]/overlap[0] < f1->mem1->rt->TOL && flag == 4)
            printf(" Krylov-4 \t %f \n",  overlap[quantumBasisSize-1]/overlap[0]);
        else {
            printf("Linear dependent! %f\t%1.16f\n",overlap[quantumBasisSize-1],overlap[0]);
            return quantumBasisSize;
        }
        
        cblas_zcopy(maxEV*stride , S , 1 , S+maxEV*stride , 1);
        gvOut = tzhegv (0,f1,Job,quantumBasisSize,T+maxEV*stride,S+maxEV*stride,stride,ritz);
        printf("eigenSolved\n");
        fflush(stdout);
        assignCores(f1, 1);
        for ( iii = 0; iii < maxEV*stride ; iii++){
            vectors[0][iii] = creal((T+maxEV*stride)[iii]);
            vectors[1][iii] = cimag((T+maxEV*stride)[iii]);
        }
    }
    
    time(&lapse_t);

    //f1->mem1->rt->eigenTime += difftime(lapse_t, start_t);
    time(&start_t);
    
    {       //printf("\nHeader,Number,CEG,Linear,BODY,CLASS,WEIGHT\n");
        for ( iii = 0; iii < imin(quantumBasisSize,Ne) ; iii++)
        {
            f1->sinc.tulip[eigenVectors+iii].value = ritz[iii];
            printf("Press%lld:%lld:,%lld ,%1.15f, %f, %d, %d , %f\n", iii+1, f1->sinc.N1,iii+1,  ritz[iii],sqr(cblas_dnrm2(quantumBasisSize, vectors[0]+(iii)*stride, 1))+sqr(cblas_dnrm2(quantumBasisSize, vectors[1]+(iii)*stride, 1)),bodies(f1,eigenVectors),irrep, deg(f1, irrep));
        }
    //    fflush(stdout);
    }


    
    
    if (flag == 3 ){
        if ( Ne > quantumBasisSize ){
            printf ("warning basis too small\n");
            exit(0);
        }
        el = outputSpace;
        
        
        INT_TYPE u;
        printf(" foundation %d < quan %d\n", foundation, quantumBasisSize );
        for ( u = 0; u < quantumBasisSize ;u++){
            if ( u < foundation )
                printf("*");
            printf ("USZ %d ---%d)\n", u+usz, tPath(f1, usz+u));
            
            
        }

        for ( u = 0; u < Ne ; u++)
            f1->sinc.tulip[el+u].path = -1;
        
        INT_TYPE sp,nG= tSize(f1->body),gi,gf,path,fpath,ipath,i,h,xx,six,ii2,jj2,kk2,fx,nm,v;
        double mag2,norms,*pointers[MaxCore];
        
        f1->sinc.tulip[eigenList].ptRank[0] = 0;
        f1->sinc.tulip[eigenList].purpose = ptObject;
        {
            ipath = 0;
            fpath = 0;
            for ( v = 0 ; v < nG*nG*nG ; v++ )
            {
                ii2 = v % nG;
                jj2 = (v/nG)%nG;
                kk2 = (v/(nG*nG))%nG;
                        if (tIR(f1->body,ii2, jj2 , kk2, irrep))
                        {
                            path = v;
                            fx = 0;
                            six = 0;//should be clustered..
                            for ( xx = 0; xx < quantumBasisSize ; xx++)
                                if ( tPath(f1, usz+xx) == path ){
                                    if  ( xx < foundation )
                                        if (! (fx++) ){
                                            f1->sinc.tulip[eigenList].name = usz+xx;
                                    //        printf("begin ");
                                        }
                                    six+= spins(f1, usz+xx)*part(f1, usz+xx);
                                   // printf ("%d,%d %d\n",  usz+xx, spins(f1, usz+xx),part(f1, usz+xx));
                                }
                           // printf(" fx %d\n", fx);
                            f1->sinc.tulip[eigenList].Current[0] = six;
                            
                            
                            if ( CanonicalRank(f1, eigenList, 0) )
                                for ( cmpl = 0; cmpl < spins(f1, usz) ; cmpl++)
                                {
                                    printf("%d:+:%d \n", path, CanonicalRank(f1, eigenList, 0));
                                    
#ifdef OMP
#pragma omp parallel for private (iii,rank,rr,h,i,nm,cmpl2,r,sp,mag2) schedule(dynamic,1)
#endif
                                    for ( iii = 0; iii < Ne  ; iii++)
                                    {
#ifdef OMP
                                        rank = omp_get_thread_num();
#else
                                        rank = 0;
#endif
                                        pointers[rank] = myStreams(f1, canonicalBuffersC, rank);
                                        rr = 0;
                                        //loop over actual memory order...
                                        for ( h = 0; h < fx ;h++)
                                            for ( i = 0; i < quantumBasisSize ; i+=foundation)
                                            {
                                                nm = h+fpath + i ;
                                                if ( nm >= quantumBasisSize ){
                                                    printf("**%d %d %d %d **", nm,h,fpath, i );
                                                    exit(0);
                                                }
                                                    
                                              //      printf ( "%d:: %d <-",iii+1, nm );
                                               // fflush(stdout);
                                                for ( cmpl2 = 0; cmpl2 < spins(f1,usz + nm ) ; cmpl2++)
                                                    for ( r = 0; r < part(f1,usz + nm ); r++){
                                                        if ( r < CanonicalRank(f1,usz + nm,cmpl) &&  cmpl2 == cmpl && tPath (f1,usz + nm) == path){
                                                            pointers[rank][rr++] = (vectors[cmpl]+(iii)*stride)[nm];
                                                        }else{
                                                            pointers[rank][rr++] = 0.0;
                                                        }
                                                    //    printf ( "->%d\n", rr );
                                                    }
                                                //fflush(stdout);

                                            }
                                     //   printf("%d %d ==%f\n", path, iii+1,cblas_dnrm2 ( rr, pointers[rank],1) );
                                        if ( rr != CanonicalRank(f1, eigenList, 0) ){
                                            printf("%d %d exit\n", rr, CanonicalRank(f1, eigenList, 0) );
                                            exit(0);
                                        }
                                        if ( cblas_dnrm2 ( rr, pointers[rank],1) > 1e-6 ){
                                            if ( rr > part(f1, canonicalBuffersC))
                                            {
                                                printf("crap!\n %lld %lld",rr , part(f1, canonicalBuffersC));
                                                exit(0);
                                            }
                                            nm = ipath * Ne + iii;
                                      //      fflush(stdout);

                                            tCycleDecompostionListOneMP(rank,f1, eigenList, pointers[rank],el+nm, cmpl, f1->mem1->rt->vCANON, part(f1, el+nm), 1.);
                                            f1->sinc.tulip[el+nm].path = path;
                                            mag2 = 0;
                                            for ( sp = 0; sp < spins(f1, el+nm) ; sp++)
                                                mag2 += inner(rank, f1, el+nm, sp);
                                          //  printf("mag %f\n", sqrt(mag2));
                                            
                                            tScale(f1, el+nm, 1./sqrt(mag2));
                               //             printf("assign %d <= %d\n", nm , path);
                                        }else {

                                        }
                                    }
                                }
                            
                            ipath += 1;
                            fpath += fx;
                        }
        }
        }
        if ( fpath != foundation ){
            printf ("wanring !!! fpath %d foundation %d", fpath, foundation);
            fflush(stdout);
            exit(0);
        }

        time(&lapse_t);
        printf("\nSpectrum TIME \t %15.0f\n", difftime(lapse_t, start_t));
        fflush(stdout);

    }else if ( flag == 4 )
        {
            double norms,*pointers[MaxCore];
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
                        //                    printf("begin\n");
                        //                    fflush(stdout);
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
                            pointers[rank] = myStreams(f1, canonicalBuffersC, rank);
                            
                            //                        printf( "%d %d\n",iii+1, cblas_idamax(stride, (vectors[cmpl]+(iii)*stride), 1));
                            
                            rr = 0;
                            norms = 0.;
                            for( g = 0; g < quantumBasisSize ; g++){
                                for ( cmpl2 = 0; cmpl2 < spins(f1, usz+g) ; cmpl2++)
                                    for ( r = 0; r < part(f1, usz+g); r++){
                                        if (  r >= CanonicalRank(f1, usz+g, cmpl2)){
                                            pointers[rank][rr++] = 0.0;
                                        }else{
                                            if( cmpl2 == cmpl ){
                                                norms += sqr((vectors[cmpl]+(iii)*stride)[g]);
                                                pointers[rank][rr++] = (vectors[cmpl]+(iii)*stride)[g];
                                            }
                                            else
                                                pointers[rank][rr++] = 0.0;//mask off other *SPIN*( or complex vector).
                                        }
                                        //    printf("%d %d %d %f\n", usz+g,cmpl2,r,myStreams(f1, canonicalBuffersC, rank)[rr-1] );
                                    }
                            }
                            if ( rr > part(f1, canonicalBuffersC))
                            {
                                printf("crap!\n %lld %lld",rr , part(f1, canonicalBuffersC));
                                fflush(stdout);

                                exit(0);
                            }
                            tCycleDecompostionListOneMP(rank,f1, eigenList, pointers[rank],el+iii   , cmpl, f1->mem1->rt->vCANON, part(f1, el+iii), 1.);
                        }
                    }
            }
            time(&lapse_t);
            printf("\nFinal Spectrum TIME \t %15.0f\n", difftime(lapse_t, start_t));
            fflush(stdout);
        }
    //f1->sinc.tulip[matrixHbuild].Current[0] = quantumBasisSize;
    //f1->sinc.tulip[matrixSbuild].Current[0] = quantumBasisSize;

    return 0;
}



