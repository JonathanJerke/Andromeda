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
                if ( f1.rt->dynamic > 0 ){
                    tEqua(f1, totalVector, rank, usr+ii, cmpl);
                    CanonicalRankDecomposition( f1, NULL, totalVector, rank, usr+ii, cmpl, f1.rt->TOLERANCE,f1.rt->relativeTOLERANCE, f1.rt->ALPHA,f1.rt->THRESHOLD, f1.rt->MAX_CYCLE,f1.rt->XCONDITION, part(f1,usr+ii),f1.rt->dynamic);
                }
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
    complete = tdsyev(0,'N',Ve+1,S,stride,ov);
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
inta tBuildMatrix (inta minusFlag,   sinc_label  f1,   division A ,   division usz, inta quantumBasisSize){
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
                    return 1;
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

/**
 *Solve one coordinate of an Optimal Separable Component Basis Hamiltonian
 *
 *
 *See Bill's Papers for definition:  Minimize all but degree of freedom interested in.
 *@param f1 container
 *@param hamiltonian  ordered by spiral
 *@param space to consider
 *@param index of the basis element to project in space
 */
floata tComponent( sinc_label f1, division hamiltonian, inta space, inta index){
    if ( ! allowQ(f1.rt, blockComponentblock)){
        printf("blockComponentblock Allow!\n");
        fflush(stdout);
        exit(0);
    }
    division A = defSpiralMatrix(&f1, hamiltonian);
    floata value1,value2,curr ,prev,target;
    inta dim ,N1,i,op;
    f1.name[component].Current[0] = 1;
    ///Initialize
    for ( dim = 0 ; dim < SPACE ; dim++ )
        if ( f1.canon[dim].body != nada )
        {
            floata *pt = streams(f1, component, 0, dim);
            N1 = vectorLen(f1, dim);
            if ( dim == space){
                for ( i = 0 ; i < N1; i++ )
                    pt[i] = 0;
                pt[index] = 1;
            } else {
                for ( i = 0 ; i < N1; i++ )
                    pt[i] = 1;
                
            }
    }
    curr = 100.;
    do{
        ///cycle...
        f1.name[componentTotal].Current[0] = 0;
        {
            floata norm = sqrt(pMatrixElement(f1, component ,0,nullOverlap,0,component ,0));
            tScaleOne(f1, component, 0, 1./norm);
        }
        for (op = 0; f1.name[A+op].species == matrix ; op++)
            tHXpY(f1, componentTotal, f1.name[A+op].name, op, component, f1.rt->TOLERANCE, f1.rt->relativeTOLERANCE, f1.rt->ALPHA, f1.rt->THRESHOLD, f1.rt->MAX_CYCLE, f1.rt->XCONDITION, part(f1,componentTotal), 0);
        prev = curr;
        curr = pMatrixElement(f1, componentTotal, 0, nullOverlap, 0, component, 0);
        tEqua(f1, component, 0, componentTotal, 0);

        for ( dim = 0 ; dim < SPACE ; dim++ )
            if ( f1.canon[dim].body != nada )
            {
                floata *pt = streams(f1, component, 0, dim);
                N1 = vectorLen(f1, dim);
                if ( dim == space){
                    for ( i = 0 ; i < N1; i++ )
                        pt[i] = 0;
                    pt[index] = 1;
                }
            }

        target = max(f1.rt->relativeTOLERANCE*curr,f1.rt->TOLERANCE);
    } while( fabs(curr-prev) > target);
    value1 = curr;
    value2 = 100;
    if ( fabs(value1) > f1.rt->THRESHOLD ){
        curr = 100.;
        do{
            ///cycle...
            {
                floata norm = sqrt(pMatrixElement(f1, component ,0,nullOverlap,0,component ,0));
                tScaleOne(f1, component, 0, 1./norm);
                
                tEqua(f1, componentTotal, 0, component, 0);
                tScaleOne(f1, componentTotal, 0, -value1);

            }
            for (op = 0; f1.name[A+op].species == matrix ; op++)
                tHXpY(f1, componentTotal, f1.name[A+op].name, 1, component, f1.rt->TOLERANCE, f1.rt->relativeTOLERANCE, f1.rt->ALPHA, f1.rt->THRESHOLD, f1.rt->MAX_CYCLE, f1.rt->XCONDITION, part(f1,componentTotal), f1.rt->dynamic);
            
            prev = curr;
            curr = pMatrixElement(f1, componentTotal, 0, nullOverlap, 0, component, 0)+value1;
            tEqua(f1, component, 0, componentTotal, 0);
            floata *pt = streams(f1, component, 0, space);
            N1 = vectorLen(f1, space);
                for ( i = 0 ; i < N1; i++ )
                    pt[i] = 0;
            pt[index] = 1;
            target = max(f1.rt->relativeTOLERANCE*curr,f1.rt->TOLERANCE);
        } while( fabs(curr-prev) > target);
    }
    value2 = curr;
    return min(value1,value2);
        
}


/**
 *Report diagonal matrix element value of Hamiltonian
 *
 *Bill said, 'do the simple stuff first'
 *All other dimensions will be occupied by an uncorrelated gaussian at origin
 *@param f1 container
 *@param hamiltonian  ordered by spiral
 *@param space to consider
 *@param index of the basis element to project in space
 */
floata tComponentPoint( calculation * c,sinc_label f1, division hamiltonian, inta space, inta index){
    if ( ! allowQ(f1.rt, blockComponentblock)){
        printf("blockComponentblock Allow!\n");
        fflush(stdout);
        exit(0);
    }
    floata value;
    inta dim ,N1,i;
    ///Initialize
    tClear(f1, component);
    tBoot(f1, component, 0, 1);
    for ( dim = 0 ; dim < SPACE ; dim++ )
        if ( f1.canon[dim].body != nada )
        {
            floata *pt = streams(f1, component, 0, dim);
            N1 = vectorLen(f1, dim);
            if ( dim == space )
            {
                for ( i = 0 ; i < N1; i++ )
                    pt[i] = 0;
                
                pt[index] = 1;
            }
        }
    value = printExpectationValues(c, f1, hamiltonian, component);
    return value;
        
}
