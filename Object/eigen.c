/**
 *  eigen.c
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


#include "eigen.h"


/**
 *Project vector into an irrep
 */
inta tFilter(  sinc_label f1, inta Ve, inta irrep,   division usr){
    inta ii,cmpl=0,rank;
    rank = 0;
    if ( ! allowQ(f1.rt,blockTotalVectorBlock)){
        printf("blockTotalVectorBlock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    if ( irrep ){
        
        printf("Symmetry Adaption -> %d-irrep\n", irrep );
        for( ii= 0 ; ii < Ve ; ii++){
            for ( cmpl = 0 ; cmpl < spins(f1, usr+ii) ; cmpl++){
                
                f1.name[totalVector].Current[rank] = 0;
                tBuildIrr(rank, f1, irrep, usr+ii, cmpl, totalVector, rank);
                CanonicalRankDecomposition( f1, NULL,totalVector, rank,  usr+ii, cmpl, f1.rt->TOLERANCE,f1.rt->relativeTOLERANCE, f1.rt->ALPHA,f1.rt->THRESHOLD, f1.rt->MAX_CYCLE,f1.rt->XCONDITION, part(f1,usr+ii),0);
            }
        }
    }
    
    for ( ii = 0; ii < Ve ; ii++)
    {
        f1.name[usr+ii].value.symmetry = tClassify( f1, usr+ii);
    }
    return 0;
}

/**
 *Filter and cap condition number, 'maxCondition', of vectors read into program
*/
inta tSelect(  sinc_label  f1, inta Ve, inta irrep,   division usr, inta testFlag){
    inta sp,rank=0,maxEV = f1.maxEV,n,m;
    inta stride = maxEV;
    if ( ! CanonicalRank(f1, usr+Ve,0)){
        printf("tSelect, null Origin\n");
        return 0;
    }

    if ( ! allowQ(f1.rt,blockMatrixElementsblock)){
        printf("blockMatrixElementsblock Allow!\n");
        fflush(stdout);
        exit(0);
    }
    if ( ! allowQ(f1.rt,blockTotalVectorBlock)){
        printf("blockTotalVectorBlock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    

    if ( irrep ) {
        for ( sp = 0; sp < spins(f1, usr+Ve);sp++){
            f1.name[totalVector].Current[rank] = 0;
            tBuildIrr(rank, f1, irrep, usr+Ve, sp, totalVector, rank);
            CanonicalRankDecomposition( f1, NULL,totalVector, rank,usr+Ve,sp, f1.rt->TOLERANCE,f1.rt->relativeTOLERANCE, f1.rt->ALPHA,f1.rt->THRESHOLD, f1.rt->MAX_CYCLE,f1.rt->XCONDITION, part(f1,name(f1,usr+Ve)),0);
        }

    }
    if ( testFlag != 1 )
        return 1;

    
    mea * S = (mea*)(myStreams(f1, matrixSbuild, 0));
    myZero(f1, matrixSbuild,0);
    double * ov = myStreams(f1, twoBodyRitz, 0);
    assignCores(f1, 1);
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
            S[n*stride+m] = tMatrixElements(rank, f1, usr+n, 0, nullOverlap, 0, usr+m, 0);
            if ( m < n ){
                S[m*stride+n] = (S[n*stride+m]);
            }
        }
        
    }
    assignCores(f1, 2);
    inta complete;
#ifdef COMPLEXME
    complete = tzheev(0, f1, 'N', Ve+1,S, stride, ov);
#else
    complete = tdsyev(0,f1,'N',Ve+1,S,stride,ov);
#endif
    if ( complete ){
        printf("eigensolve of overlap failed\n");
    } else 

    if ( testFlag ){
        if ( ov[Ve]/ov[0] < f1.rt->XCONDITION && ov[Ve]/ov[0] >0 && ov[0] > 0 ){
                printf("\t**%f**\t%d\n",  ov[Ve]/ov[0],Ve+1 );
                fflush(stdout);

                return 1;
        } else {
        }
    }
    
    return 0;
}

/**
 *Build Ritz Matrix Elements
*/
inta tBuildMatrix (inta minusFlag,   sinc_label  f1,   division A ,    division usz, inta quantumBasisSize){
    inta in,prevBuild;
    inta rank;
    inta i,stride = f1.maxEV;
    if ( ! allowQ(f1.rt, blockMatrixElementsblock)){
        printf("blockMatrixElementsblock Allow!\n");
        fflush(stdout);
        exit(0);
    }
    
    
    mea *T  =  (mea *) myStreams(f1, matrixHbuild,0/*CORE RANK*/);
    mea *S  =  (mea *) myStreams(f1, matrixSbuild,0/*CORE RANK*/);
    assignCores(f1, 1);

    if ( part(f1, matrixHbuild) * sizeof(double) < stride*stride*sizeof(mea) ){
        printf("somehow there is not enough allocation for matrixHbuild\n");
        exit(0);
    }
    if ( part(f1, matrixSbuild) * sizeof(double) < stride*stride*sizeof(mea) ){
        printf("somehow there is not enough allocation for matrixSbuild\n");
        exit(0);
    }
    
    for ( i = 0 ; i < stride*stride ; i++){
        S[i] = 0.;
        T[i] = 0.;
    }

    inta n,m;
      division leftP ;
    prevBuild = 0;
//    for ( n = prevBuild; n < quantumBasisSize ; n++)
//    {
//        if ( magnitude(f1, usz+n,0) == 0. )
//        {
//            printf("The vector %d is zero!\n",n);
//            exit(0);
//        }
//    }

        
#ifdef OMP
#pragma omp parallel for private (in,m,rank,n) schedule(dynamic,1)
#endif
            for ( in = 0 ;in < quantumBasisSize * quantumBasisSize; in++)
            {
                
#ifdef OMP
                rank = omp_get_thread_num();
#else
                rank = 0;
#endif
                m = in % quantumBasisSize;
                n = (in/quantumBasisSize) % quantumBasisSize;
                if ( m<=n ){
                    S[n*stride+m] = tMatrixElements(rank, f1, usz+n,0, nullOverlap, 0, usz+m, 0);
                }
            }
          
    

    leftP = A;

    do {
        if ( f1.name[leftP].species == matrix ){
#ifdef CHERRY_PICKER
    #ifdef OMP
    #pragma omp parallel for private (in,m,rank,n) schedule(dynamic,1)
    #endif
                for ( in = 0 ;in < quantumBasisSize * quantumBasisSize; in++)
                {
                    
    #ifdef OMP
                    rank = omp_get_thread_num();
    #else
                    rank = 0;
    #endif
                    m = in % quantumBasisSize;
                    n = (in/quantumBasisSize) % quantumBasisSize;
                    if ( m<=n ){
                        T[n*stride+m] += tMatrixElements(rank, f1, usz+n,0, leftP, 0, usz+m, 0);
                    }
                }
#else
                    for ( m = 0 ;m < quantumBasisSize; m++)
                    {
                        
                        tHXpY(f1, totalVector, leftP, 0, usz+m, 0, 0, 0, 0, 0,1e6, CanonicalRank(f1,usz+m,0)*CanonicalRank(f1, leftP, 0), CanonicalRank(f1,usz+m,0)*CanonicalRank(f1, leftP, 0));
#ifdef OMP
#pragma omp parallel for private (rank,n) schedule(dynamic,1)
#endif
                        for ( n = 0 ;n < quantumBasisSize; n++){
                            #ifdef OMP
                                            rank = omp_get_thread_num();
                            #else
                                            rank = 0;
                            #endif

                            if ( m<=n ){
                                T[n*stride+m] += tMatrixElements(rank, f1, usz+n,0, nullOverlap, 0, totalVector, 0);
                            }
                        }
                    }
    
#endif
                }
              leftP = f1.name[leftP].linkNext;
          } while ( leftP != nullName);
    
    
    
    if(1)
    for ( in = 0 ;in < quantumBasisSize * quantumBasisSize; in++){
                    m = in % quantumBasisSize;
                    n = (in/quantumBasisSize) % quantumBasisSize;
                    if ( m < n ){
            #ifdef COMPLEXME
                        T[m*stride+n] = conj(T[n*stride+m] ) ;
                        S[m*stride+n] = conj(S[n*stride+m ]) ;
                #else
                        T[m*stride+n] = (T[n*stride+m] ) ;
                        S[m*stride+n] =(S[n*stride+m ]) ;
                #endif
                    }
                    
    }
    
#ifdef COMPLEXME
    if (minusFlag){
        DCOMPLEX one = -1.;
        cblas_zscal(quantumBasisSize*quantumBasisSize,&one,T,1);
    }
#else
    if (minusFlag){
        mea one = -1.;
        cblas_dscal(quantumBasisSize*quantumBasisSize,one,T,1);
    }
#endif
    
    return stride;
    
}
    

/**
 *Ritz eigensolve
 */
inta tSolveMatrix (inta typer,   sinc_label  f1,inta Ne,  division usz, inta quantumBasisSize,   division outputValues){

    if ( ! allowQ(f1.rt, blockMatrixElementsblock)){
        printf("blockMatrixElementsblock Allow!\n");
        fflush(stdout);
        exit(0);
    }

            inta qs;
            inta iii = 0,maxEV = f1.maxEV;
            inta stride = maxEV;
            double * ritz = myStreams(f1, outputValues, 0);
            mea *T  =  (mea *) myStreams(f1, matrixHbuild,0/*CORE RANK*/);
            mea *S  =  (mea *) myStreams(f1, matrixSbuild,0/*CORE RANK*/);
    
            qs = quantumBasisSize;
            
            if (1){
                assignCores(f1, 2);
                inta complete;
                char Job = 'V';
#ifdef COMPLEXME
                complete = tzhegv (0,f1,Job,qs,T,S,stride,ritz);
              //  cblas_zcopy(stride*stride, t, 1, T, 1);
#else
                complete = tdsygv (0,f1,Job,qs,T,S,stride,ritz);
               // cblas_dcopy(stride*stride, t, 1, T, 1);
#endif
                if ( complete ){
                    printf("eigensolve -Ritz- failed\t%d\n",complete);
                    exit(0);
                    return 0;
                }
                
    #if VERBOSE
                printf("eigenSolved\n");
    #endif
                fflush(stdout);
            }
            {
                for ( iii = 0; iii < imin(qs,Ne) ; iii++)
                {
                    f1.name[eigenVectors+iii].value.value = ritz[iii];
                    printf("\nPress%d:,%1.15f, %f", iii+1, f1.name[eigenVectors+iii].value.value,cblas_dznrm2(qs,T+iii*stride, 1));
                    fflush(stdout);
                }
            }
        
        //}

    return 0;
}
