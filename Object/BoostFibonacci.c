/**
*  BoostFibonacci.c
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

#include "Fibonacci.h"

/**
 * Controls to initiate and run canonicalRankDecomposition
 *
 *@param f0          container
 *@param coeff allows for rescaling the origin by each canonical rank
 *@param[in] origin the division with more canonical ranks
 *@param os         spin of origin to consider
 *@param[out] alloy  the division with less canonical ranks, to be trained
 *@param spin      spin of alloy to consider
 *@param tolerance a number setting the absolute quality
 *@param relativeTolerance a number seting quality relative to magnitude of origin
 *@param condition Beylkin's condition (alpha)
 *@param threshold the smallest number
 *@param maxCycle the maxmium number of cycles in this routine
*/
double CanonicalRankDecomposition (  sinc_label  f0 ,double * coeff,   division origin,inta os,  division alloy,inta spin, double tolerance, double relativeTolerance, double condition,double threshold, inta maxCycle ,double maxCondition, inta canon ){
    inta rank;
    division G = nullName;
    inta ii,n,m,c,g,G1 = CanonicalRank(f0, origin, os);
    
    inta L1 = canon,GG1;
        
    inta *iiii[2][2];
    inta *iii[2][2];
    
    if ( ! G1 ){
        printf("CanonicalRankDecomposition, Origin is empty\n");
        return 0;
    }
    if ( G1 <= L1 ){
        tEqua(f0, alloy, spin, origin, os);
        return 0;
    }
    
    ///ADDITION BEGIN
    inta space;
    field F1;
    F1.f = f0;
    for ( space = 0; space < SPACE ; space++)
        if ( F1.f.canon[space].body != nada ){
            F1.f.canon[space].count1Basis = imin( G1,vectorLen(f0, space));
            F1.f.canon[space].body = one;
        }
    
    
    F1.f.bootedMemory = 0;
    for ( space = 0; space < SPACE ; space++)
        F1.f.canon[space].stream = NULL;
    F1.f.name = NULL;

    F1.f.boot = noMatrices;
    F1.i.body = one;
    F1.i.canonRank = G1;
    F1.i.cmpl = real;
    F1.i.collect = 0;
    F1.i.files = 0;
    F1.i.filesVectorOperator = 0;
    F1.i.filter = 0;
    F1.i.flex = 0;
    F1.i.iRank = 0;
    F1.i.irrep = 0;
    F1.i.Iterations = 1;
    F1.i.matrices = 0;
    F1.i.nOperator = 0;
    F1.i.nStates = 1;
    F1.i.OpIndex = 0;
    F1.i.qFloor = 0;
    F1.i.xRank = 0;

    calculation c1;
    c1.rt = *f0.rt;
    resetA(&c1.rt);
    allowQ(&c1.rt, blockCopyBlock);
    allowQ(&c1.rt, blockTransferBasisblock);
    allowQ(&c1.rt, blockMatrixElementsblock);
    allowQ(&c1.rt, blockPermutationsblock);
    allowQ(&c1.rt, blockParallelMultiplyblock);
    allowQ(&c1.rt, blockParallelMatrixElementblock);
    allowQ(&c1.rt, blockParallelPermuteblock);
    allowQ(&c1.rt, blockPrintStuffblock);
    allowQ(&c1.rt, blockTotalVectorParallelBlock);

    F1.f.eikonLabels.currLabel = 0;
    F1.f.nullLabels.currLabel = 0;

    c1.rt.NLanes = f0.rt->NLanes;
    c1.i.Angstroms = 0;
    c1.i.lambda = 1;
    c1.i.minIterationPrint =0;
    c1.i.Na = 0;
    c1.i.numNames = G1+1000;
    c1.i.numVectors = 0;
    c1.i.RAMmax = 10000;
    c1.i.shiftFlag = 0;
    c1.i.termNumber = 0;
    iModel(&c1, &F1);
    
    
    zero(F1.f, totalVector, 0);
    
    #ifdef OMP
    #pragma omp parallel for private (space,rank) schedule(dynamic,1)
    #endif
    for ( space = 0; space < SPACE ; space++)
        if ( f0.canon[space].body != nada){
            #ifdef OMP
                    rank = omp_get_thread_num();
            #else
                    rank = 0;
            #endif

            if ( G1 < vectorLen(f0, space) )
                tdgeqr(rank, f0, G1, vectorLen(f0, space), streams(f0,origin,os,space), vectorLen(f0, space), myStreams(f0,CanonicalBuffers,rank), streams(F1.f,totalVector,0,space), G1);
            else
                cblas_dcopy(G1*vectorLen(f0, space), streams(f0,origin,os,space), 1, streams(F1.f,totalVector,0,space), 1);
        }
    
    ///ADDITION END


    {
        floata * me;
        me = myStreams(F1.f, CanonicalBuffers, 0);

           ///RENAME
           division origin = totalVector;
            
           division alloy = eigenVectors;
           inta spin = 0,os = 0;
           zero(F1.f,alloy,spin);

           ///RENAME END
        inta M2[SPACE];
        length(F1.f, origin, M2);

    for ( g = 0 ; g < G1 ; g++){
        if ( G== nullName )
            G = anotherLabel(&F1.f, 0, nada);
        else
            anotherLabel(&F1.f, 0, nada);///will be in place
        F1.f.name[G+g].name = origin;
        F1.f.name[G+g].Partition = 1;
        for ( c = 0 ; c < MAX_CORE; c++){
            F1.f.name[G+g].Begin[c] = 0;
            F1.f.name[G+g].Current[c] = 0;
        }
        F1.f.name[G+g].Begin[os] = g;
        F1.f.name[G+g].Current[os] = 1;
        if ( g + 1 < G1 )
            F1.f.name[G+g].chainNext = G+g+1;
        else
            F1.f.name[G+g].chainNext =nullName;

    }
       

       for ( ii= 0; ii < G1 ; ii++){
                   
                   rank = 0;
           
           me[G1*ii+ii]  = tMatrixElements(rank, F1.f, G+ii, 0, nullOverlap, 0, G+ii, 0);
       }
    
    inta iv = 0;
    for (ii= 0; ii < G1 ; ii++)
        if ( me[G1*ii+ii] > f0.rt->THRESHOLD ){
            if ( ii > iv ){
                for ( space = 0; space < SPACE ; space++)
                    if ( F1.f.canon[space].body != nada){
                        cblas_dcopy(M2[space], streams(F1.f,origin,0,space)+ii*M2[space], 1, streams(F1.f,origin,0,space)+iv*M2[space], 1);
                    }
            }
            iv++;
        }
        GG1 = G1;
        G1 = iv;

    for ( ii= 0; ii < G1*G1 ; ii++){

        m = ii%G1;
        n = (ii/G1)%G1;
                
        rank = 0;
        
        if ( m <= n ){
            me[G1*n+m]  = tMatrixElements(rank, F1.f, G+n, 0, nullOverlap, 0, G+m, 0);
            me[G1*m+n]  = me[G1*n+m];
        }
    }

    
    
    iii[0][0] = malloc(L1 * sizeof( inta ));
    iii[0][1] = malloc(L1 * sizeof( inta ));
    iii[1][0] = malloc(L1 * sizeof( inta ));
    iii[1][1] = malloc(L1 * sizeof( inta ));

    iiii[0][0] = malloc(L1 * sizeof( inta ));
    iiii[0][1] = malloc(L1 * sizeof( inta ));
    iiii[1][0] = malloc(L1 * sizeof( inta ));
    iiii[1][1] = malloc(L1 * sizeof( inta ));
    
    
    for ( ii = 0; ii < L1 ; ii++)
    {
        iii[0][0][ii] = ii;
        iii[0][1][ii] = ii+1;
        iii[1][0][ii] = ii*( G1/L1 );
        iii[1][1][ii] = (ii+1)*(G1/L1);
    }
    iii[1][1][L1-1] = G1;

        F1.f.name[origin].Current[os] = G1;
        F1.f.name[alloy].Current[spin] = L1;
        inta nRun = L1;
        
        inta run,r2;
        
        while ( nRun > 0 ){
            {
                for ( r2 = 0; r2 < nRun; r2++ ){
                    canonicalRankDecomposition(F1.f, coeff, G1,me, origin,iii[1][0][r2],iii[1][1][r2], os,nRun == L1,alloy, iii[0][0][r2],iii[0][1][r2],spin,tolerance, relativeTolerance,condition,maxCondition,maxCycle);
                }
            }
            //merge
            r2 = 0;
            if ( nRun % 2 == 0 ){
                for ( run = 0; run < nRun ; run+=2){
                    iiii[1][0][r2] = iii[1][0][run];
                    iiii[1][1][r2] = iii[1][1][run+1];
                    iiii[0][0][r2] = iii[0][0][run];
                    iiii[0][1][r2] = iii[0][1][run+1];
                    r2++;
                }
            }
            else
            {
                for ( run = 0; run < nRun-1 ; run+=2){
                    iiii[1][0][r2] = iii[1][0][run];
                    iiii[1][1][r2] = iii[1][1][run+1];
                    iiii[0][0][r2] = iii[0][0][run];
                    iiii[0][1][r2] = iii[0][1][run+1];
                    r2++;
                }
                if ( r2 ) {
                    iiii[1][1][r2-1] = iii[1][1][nRun-1];
                    iiii[0][1][r2-1] = iii[0][1][nRun-1];
                }
                
            }
            if ( nRun == 1 )
                nRun= 0;
            else
                nRun = r2;
            
            for ( run = 0; run < nRun ; run++){
                iii[1][0][run] = iiii[1][0][run];
                iii[1][1][run] = iiii[1][1][run];
                iii[0][0][run] = iiii[0][0][run];
                iii[0][1][run] = iiii[0][1][run];
            }
            
        }
    }
    free(iiii[0][0]);
    free(iiii[0][1]);
    free(iiii[1][0]);
    free(iiii[1][1]);
    free(iii[0][0]);
    free(iii[0][1]);
    free(iii[1][0]);
    free(iii[1][1]);

    ///ADDITION BEGIN
    #ifdef OMP
    #pragma omp parallel for private (space) schedule(dynamic,1)
    #endif
    for ( space = 0; space < SPACE ; space++)
        if ( f0.canon[space].body != nada){
            if ( GG1 < vectorLen(f0, space) ){
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, vectorLen(f0, space), L1, G1, 1., streams(f0,origin,os,space), vectorLen(f0, space), streams(F1.f,eigenVectors,0,space), G1, 0., streams(f0,alloy,spin,space),vectorLen(f0, space) );
            } else {
                cblas_dcopy(L1*vectorLen(f0, space),streams(F1.f,eigenVectors,0,space),1,streams(f0,alloy,spin,space),1);
            }
        }
    f0.name[alloy].Current[spin] = L1;
    fModel( &F1.f);
    ///ADDITION END
    
    return 0.;
}
