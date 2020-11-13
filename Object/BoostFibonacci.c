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
double CanonicalRankDecomposition (  sinc_label  f0 ,double * coeff,   division origin,inta os,  division alloy,inta spin, double tolerance, double relativeTolerance, double condition,double threshold, inta maxCycle ,double maxCondition, inta canon ,inta X1){
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
    if ( G1 <= L1-X1 ){
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
    blockA(&c1.rt, blockCopyBlock);
    blockA(&c1.rt, blockTransferBasisblock);
    blockA(&c1.rt, blockMatrixElementsblock);
    blockA(&c1.rt, blockPrintStuffblock);
    blockA(&c1.rt, blockTotalVectorParallelBlock);

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
    
    if ( ! allowQ(F1.f.rt,blockTotalVectorBlock)){
        printf("blockTotalVectorBlock Allow!\n");
        fflush(stdout);
        exit(0);
    }
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
            if ( part(f0,CanonicalBuffers) < G1 ){
                printf("CanonicalBuffers0 \n");
                exit(0);
            }

            if ( G1 < vectorLen(f0, space) )
                tdgeqr(rank, f0, G1, vectorLen(f0, space), streams(f0,origin,os,space), vectorLen(f0, space), myStreams(f0,CanonicalBuffers,rank), streams(F1.f,totalVector,0,space), G1);
            else
                cblas_dcopy(G1*vectorLen(f0, space), streams(f0,origin,os,space), 1, streams(F1.f,totalVector,0,space), 1);
        }
    
    ///ADDITION END


    {
        floata * me;
        me = myStreams(F1.f, CanonicalBuffers, 0);
        if ( part(F1.f,CanonicalBuffers) < G1*G1 ){
            printf("CanonicalBuffers1 \n");
            exit(0);
        }
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
           me[ii]  = pMatrixElement( F1.f, G+ii, 0, nullOverlap, 0, G+ii, 0);
       }
    
        
        ///first determine a matrix-rank.
    inta iv = 0;
    for (ii= 0; ii < G1 ; ii++)
        if ( me[ii] > f0.rt->THRESHOLD ){
            if ( ii > iv ){
                for ( space = 0; space < SPACE ; space++)
                    if ( F1.f.canon[space].body != nada){
                        cblas_dcopy(M2[space], streams(F1.f,origin,os,space)+ii*M2[space], 1, streams(F1.f,origin,os,space)+iv*M2[space], 1);
                    }
            }
            iv++;
        }
        GG1 = G1;
        G1 = iv;

        
        ///now determine a overlap of terms
    for ( ii= 0; ii < G1*G1 ; ii++){

        m = ii%G1;
        n = (ii/G1)%G1;
                        
        if ( m <= n ){
            me[G1*n+m]  = pMatrixElement(F1.f, G+n, 0, nullOverlap, 0, G+m, 0);
            me[G1*m+n]  = me[G1*n+m];
        }
    }

        if ( X1 == 0 ){
        
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
                
            free(iiii[0][0]);
            free(iiii[0][1]);
            free(iiii[1][0]);
            free(iiii[1][1]);
            free(iii[0][0]);
            free(iii[0][1]);
            free(iii[1][0]);
            free(iii[1][1]);
        }else {
            ///assume you have already worked out decent solutions, now try to decompose them further...
            F1.f.name[origin].Current[os] = G1;
            F1.f.name[alloy].Current[spin] = canon;

            tEqua(F1.f, alloy, spin, origin, os);

            inta i,x1;
            floata target,Glen = 0.,curr;
            for ( i = 0; i < G1*G1 ; i++)
                Glen += me[i];
            
            for ( x1 = 1 ; x1 <= X1 ; x1++){
                F1.f.name[alloy].Current[spin] = canon-x1;
                printf("attempting to decompose to %d \n", canon-x1);

                curr = canonicalRankDecomposition(F1.f, coeff, G1, me, origin, 0, canon, os, 0, alloy, 0, canon-x1, spin, tolerance, relativeTolerance, condition, maxCondition, maxCycle);
                target = max( tolerance , Glen*relativeTolerance );
                if ( curr > target){
                    L1 -=x1+1;
                    printf("stop decomposing at %d above %f\n", canon-x1+1, target);
                    F1.f.name[alloy].Current[spin] = canon-x1+1;

                    canonicalRankDecomposition(F1.f, coeff, G1, me, origin, 0, canon, os, 0, alloy, 0, canon-x1+1, spin, tolerance, relativeTolerance, condition, maxCondition, maxCycle);
                    break;
                }
                printf("decomposed to %d \t %f \t %f\n", canon-x1,curr,Glen);
            }
        }
    }
    ///ADDITION BEGIN
    #ifdef OMP
    #pragma omp parallel for private (space) schedule(dynamic,1)
    #endif
    for ( space = 0; space < SPACE ; space++)
        if ( f0.canon[space].body != nada){
            if ( GG1 < vectorLen(f0, space) ){
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, vectorLen(f0, space), L1, vectorLen(F1.f, space), 1., streams(f0,origin,os,space), vectorLen(f0, space), streams(F1.f,eigenVectors,0,space), vectorLen(F1.f,space), 0., streams(f0,alloy,spin,space),vectorLen(f0, space) );
            } else {
                cblas_dcopy(L1*vectorLen(f0, space),streams(F1.f,eigenVectors,0,space),1,streams(f0,alloy,spin,space),1);
            }
        }
    f0.name[alloy].Current[spin] = L1;
    fModel( &F1.f);
    ///ADDITION END
    
    return 0.;
}
