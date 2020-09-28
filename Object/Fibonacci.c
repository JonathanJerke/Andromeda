/**
*  Fibonacci.c
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
 *@param rank      OMP rank
 *@param f1          container
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
double CanonicalRankDecomposition ( inta rank,  sinc_label  f1 ,double * coeff,   division origin,inta os,  division alloy,inta spin, double tolerance, double relativeTolerance, double condition,double threshold, inta maxCycle ,double maxCondition, inta canon ){
    
    inta L1 = canon;
        
    inta *iiii[2][2];
    inta *iii[2][2];
    floata  *me = myStreams(f1, CanonicalBuffers, 0);
    
    division G = nullName;
    inta ii,n,m,c,g,G1 = CanonicalRank(f1, origin, os);
    
    if ( ! G1 ){
        printf("CanonicalRankDecomposition, Origin is empty\n");
        return 0;
    }
    if ( G1 <= L1 ){
        tEqua(f1, alloy, spin, origin, os);
        return 0;
    }

    for ( g = 0 ; g < G1 ; g++){
        if ( G== nullName )
            G = anotherLabel(&f1, 0, nada);
        else
            anotherLabel(&f1, 0, nada);//will be in place
        f1.name[G+g].name = origin;
        f1.name[G+g].Partition = 1;
        for ( c = 0 ; c < MAX_CORE; c++){
            f1.name[G+g].Begin[c] = 0;
            f1.name[G+g].Current[c] = 0;
        }
        f1.name[G+g].Begin[os] = g;
        f1.name[G+g].Current[os] = 1;
        if ( g + 1 < G1 )
            f1.name[G+g].chainNext = G+g+1;
        else
            f1.name[G+g].chainNext =nullName;

    }
    assignCores(f1, 1);
       
       #ifdef OMP
       #pragma omp parallel for private (ii,n,m,rank) schedule(dynamic,1)
       #endif

       for ( ii= 0; ii < G1*G1 ; ii++){

           m = ii%G1;
           n = (ii/G1)%G1;
                   
           #ifdef OMP
                   rank = omp_get_thread_num();
           #else
                   rank = 0;
           #endif
           
           if ( m <= n ){
               me[G1*n+m]  = tMatrixElements(rank, f1, G+n, 0, nullOverlap, 0, G+m, 0);
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
    
    f1.name[alloy].Current[spin] = L1;
    inta nRun = L1;
    
    inta run,r2;
    
    while ( nRun > 0 ){
        {
            for ( r2 = 0; r2 < nRun; r2++ ){
                canonicalRankDecomposition(f1, coeff, G1,me, origin,iii[1][0][r2],iii[1][1][r2], os,nRun == L1,alloy, iii[0][0][r2],iii[0][1][r2],spin,tolerance, relativeTolerance,condition,maxCondition,maxCycle);
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

    return 0.;
}
