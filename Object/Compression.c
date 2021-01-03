/**
*  Compression.c
*
*
*  Copyright 2021 Jonathan Jerke and Bill Poirier.
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


#include "Compression.h"


/**
 *Canonical Rank Compression or Alternating Least Squares (ALS)
 *
 *@param[in] spatial Defining how origin breaks down to alloy
 *@param[in] cofact allows for rescaling the origin by each canonical rank
 *@param f1 origin container
 *@param[in] GG          internal from ASTER, will speedup processes
 *@param[in] origin the division with more canonical ranks
 *@param l1         Allowance for starting somewhere else but the first index
 *@param l2         last index of origin
 *@param os         spin of origin to consider
 *@param f2 alloy container
 *@param neo    flag for replacing last canonRank
 *@param[in,out] alloy  the division with less canonical ranks, to be trained
 *@param l3         alloy's first index
 *@param l4         alloy's last index
 *@param spin      spin of alloy to consider
 *@param tolerance a number setting the absolute quality
 *@param relativeTolerance a number seting quality relative to magnitude of origin
 *@param condition Beylkin's condition (alpha)
 *@param maxCondition the smallest number
 *@param maxCycle the maxmium number of cycles in this routine
 */
floata canonicalRankCompression( inta  spatial[SPACE][SPACE], floata * cofact,sinc_label  f1 ,inta G,floata *GG, division origin,inta l1,inta l2,inta os, sinc_label  f2 ,inta neo,division alloy ,inta l3 , inta l4, inta spin ,double tolerance,double relativeTolerance, double condition,double maxCondition, inta maxCycle){
    if (l2 < l1 || l4 < l3 )
    {
        printf("indices out of order!\n");
        printf("%d %d %d %d\n", l1,l2,l3,l4);
        printf("%d %d\n", origin, alloy);
        exit(0);
    }
    if ( ! allowQ(f1.rt,blockTrainVectorsblock)){
        printf("blockTrainVectorsblock allow 1!\n");
        fflush(stdout);
        exit(0);
    }
    if ( ! allowQ(f2.rt,blockTrainVectorsblock)){
        printf("blockTrainVectorsblock allow 2!\n");
        fflush(stdout);
        exit(0);
    }

    inta xM = 0;
    {inta space;
    for ( space = 0; space < SPACE ; space++)
        if( f1.canon[space].body != nada )
            if ( xM < vectorLen(f1, space))
                xM = vectorLen(f1,space);
    }
    
    
    inta rank = 0;
    double xprod = 0.;
    inta xOriginIndex=0;
    double sum2 ,iGG=0,iGF=0,iFF=0;
    inta flagAllDim,spaces = 0,spaces2 = 0;
    inta g,l,count = 1,space2;
    double target;
    floata * array[SPACE],curr=1e6,prev;
    floata * array2[SPACE],*norm[SPACE];
    floata * guide, *track,*tracker;
    inta space,spaceX=0,space0=0,flipSignFlag=0;
    inta L1 = l4-l3;
    inta M2[SPACE];
    inta M1[SPACE];
    inta G1 = l2-l1;
    for ( space = 0; space < SPACE ; space++)
        if ( f1.canon[space].body != nada )
            spaces++;
    for ( space2 = 0; space2 < SPACE ; space2++)
        if ( f2.canon[space2].body != nada )
            spaces2++;

    if ( G1 == 0 ){
        printf("CD: empty origin %d _%d -> %d _%d\n", origin, os , alloy, spin);
        return 0;
    }
    
    length(f1,alloy,M1);
    length(f2,alloy,M2);
    inta info;
    division canonicalStore;
    inta LS1 = L1;
    {
        for ( space2 = 0; space2 < SPACE ; space2++)
            if( f2.canon[space2].body != nada )
            {
                ///MAY WANT TO SEPARATE OWNERSHIP OF BUFFERS BETWEEN f1 and f2
                array[space2] =  streams(f2, canonicalBuffers, rank , space2);
                norm[space2] = array[space2] + L1*L1;
            }else{
                array[space2] = NULL;
                norm[space2] = NULL;
            }
            for ( space = 0; space < SPACE ; space++){
                if( f1.canon[space].body != nada )
                    array2[space] =  streams(f1, canonicalBuffers, rank , space);
                else
                    array2[space] = NULL;
            }
        guide =  myStreams(f1, guideBuffer, rank );
        track =  myStreams(f1, trackBuffer, rank );
        tracker =  myStreams(f1, trackBuffer, rank )+L1*L1*2;

        
        if (   L1*L1+L1 >  part(f2,canonicalBuffers) || L1*G1 > part(f1,canonicalBuffers)|| G1*L1 > part(f1,guideBuffer) || L1*L1*2+L1 > part(f1,trackBuffer)){
#if 1
            printf("canonicalRankDecomposition, mem prob with canonicalBuffers, guideBuffer, or trackBuffer\n %d %d \n", L1, G1);
            printf("canonicalBuffers %d %d \n guideBuffer %d\n trackBuffer %d\n",part(f1,canonicalBuffers),part(f2,canonicalBuffers),part(f1,guideBuffer),part(f1,trackBuffer) );
            fflush(stdout);
#endif
            exit(0);
        }
        
        
        canonicalStore = canonicalBuffersD;
        
        if ( part(f1, canonicalStore) < 3 * xM )
        {
            printf("canonicalRankDecomposition, mem prob with canonicalBuffersB\n");
            fflush(stdout);

            exit(0);
        }
        
        
    }
    floata *pt[f1.rt->NLanes] ;
    floata *ot[f1.rt->NLanes] ;
    floata *qt[f1.rt->NLanes] ;

    {
        inta r;
        for ( r = 0 ; r < f1.rt->NLanes ; r++)
        {///CHANGED FROM LENGTH L1 to G1
            pt[r] = myStreams(f1,canonicalStore , r);
            ot[r] = myStreams(f1,canonicalStore , r)+xM;
            qt[r] = myStreams(f1,canonicalStore , r)+xM*2;
        }
    }
    
    inta dim0 = 0;
    inta dim[SPACE];

    floata *** originStream = malloc(SPACE * sizeof(floata**));
    floata *** alloyStream = malloc(SPACE * sizeof(floata**));
    inta * originIndex = malloc( G1* sizeof(inta));
    zero(f1, canonicalVector, rank);

    for ( space = 0; space < SPACE;space++)
        if ( f1.canon[space].body != nada )
    {
        originStream[space] = malloc((G1+1)*sizeof(floata*));
        originStream[space][G1] = streams(f1,canonicalVector,rank,space);//

        division nIter;
        inta n,ni,nc;
        nIter = origin;
        ni = 0;
        nc = f1.name[nIter].Current[os];
        for ( n = 0; n < G1+l1 ; n++){
            while( n >= ni+nc ){
                ni += nc;

                   nIter = f1.name[nIter].chainNext;
                nc = f1.name[nIter].Current[os];

                   if ( nIter == nullName){
                       printf("somehow the pointers in Origin lost track, you probably have too many canonical ranks\n");
                       exit(0);
                       }
               }
            if ( n >= l1 ){
                originStream[space][n-l1] = streams(f1, nIter, os, space) + (n - ni)*M1[space] ;
                originIndex[n-l1] = f1.name[nIter].Begin[os]+(n - ni);
//                if  ( ! space)
//                    printf("%d%d<%d %d %d>\n",origin,alloy, n-l1, n, originIndex[n-l1]);
                if ( cofact != NULL && ! space )
                    cblas_daxpy(M1[space], cofact[originIndex[n-l1]] , originStream[space][n-l1], 1, originStream[space][G1], 1);
                else
                    cblas_daxpy(M1[space],1. , originStream[space][n-l1], 1, originStream[space][G1], 1);

#if VERBOSE
                printf("%d>%d<%ld\n", n-l1,space, (originStream[space][n-l1]-streams(f1,origin,os,space)));
                fflush(stdout);
#endif
            }

        }
    }
        for ( space = 0; space < SPACE;space++)
            if ( f2.canon[space].body != nada )
            {
                alloyStream[space] = malloc(L1*sizeof(floata*));

                inta n,ni,nc;
                division nIter;
                nIter = alloy;
                nc = f2.name[nIter].Current[spin];
                ni = 0;
                for ( n = 0; n < L1+l3 ; n++){
                    
                    while ( n >= ni+nc ){
                        ni += nc;
                            nIter = f2.name[nIter].chainNext;
                        nc = f2.name[nIter].Current[spin];

                            if ( nIter == nullName){
                                printf("somehow the pointers in Alloy lost track, you probably have too many canonical ranks\n");
                                exit(0);
                                }
                    }
                    if ( n >= l3){
                        alloyStream[space][n-l3] = streams(f2, nIter, spin, space) +(n - ni)*M2[space] ;
#if VERBOSE
                        printf("%d<%d>%ld\n", n-l3,space, (alloyStream[space][n-l3]-streams(f2, alloy, spin, space)));
                        fflush(stdout);
#endif
                    }

                }
            }
    
    
        
    
        
        if ( GG == NULL ){
            ///intended for single continuous origin, even if pointers are all messy.  Tied to ASTER PROCEDURE.
            ///to avoid redundancy of computation.
            floata product,sum=0. ;
            inta l,ll,space;
            for ( l = 0 ; l < G1 ; l++)
                for ( ll = 0; ll< G1 ; ll++){
                    product = 1.;
                    for (space = 0; space < SPACE ; space++)
                        if ( f1.canon[space].body != nada )
                        product *= cblas_ddot(M1[space], originStream[space][l], 1, originStream[space][ll], 1);
                    if ( cofact != NULL ){
                        product *= (cofact)[originIndex[l]]*(cofact)[originIndex[ll]];
                    }
                    if ( product > xprod){
                        xprod = product;
                        xOriginIndex = ll;
                    }
                    sum += product;
                }
            iGG = sum;
            
            
            
            ///here one could argue for iGG begin too small.
            
            } else {
                floata product,sum=0. ;
                inta l,ll;
                for ( l = 0 ; l < G1 ; l++)
                    for ( ll = 0; ll< G1 ; ll++){
                            product = GG[originIndex[l]*G+originIndex[ll]];
                            if ( cofact != NULL ){
                                product *= (cofact)[originIndex[l]]*(cofact)[originIndex[ll]];
                            }
                        if ( product > xprod){
                            xprod = product;
                            xOriginIndex = ll;
                        }

                            sum += product;
                    }
                iGG = sum;
            }
#if VERBOSE
            printf("iGG %f\n",iGG);
            fflush(stdout);
#endif

            target = max( iGG*relativeTolerance , tolerance );
    
            if(neo)
            { inta space2,space,ii,i,stride;
                for ( space = 0; space < SPACE;space++)
                    if ( f1.canon[space].body != nada ){
                        stride = 1;
                        ///insert simplistic component separation here
                        for ( space2 = 0; space2 < SPACE;space2++)
                            if ( f2.canon[space2].body != nada )
                            
                            if ( spatial[space][space2])
                            { floata *pt =alloyStream[space2][0];
                                inta mm,kk = 1;
                                for ( mm = 0 ; mm < M2[space2] ; mm++)
                                    pt[mm] = cos(kk*2*pi*mm*1./M2[space2]);
                                kk++;
                            }
                    }
                ///separate output ..which means first rank is done...
            }
    
        
            
        ///orders
        { inta space,dim2;
            dim0=0;
            for ( dim2 = 0; dim2 < SPACE ; dim2++)
                if ( f2.canon[dim2].body != nada)
                    dim[(dim2+spaceX)%(spaces2)] = dim0++;
        
            for ( space = 0; space < SPACE ; space++)
                if ( f1.canon[space].body != nada)
                    if ( spatial[space][dim[0]] ){
                        space0 = space;
                        break;
                    }
        }///orders

        flagAllDim = 1;
    
        while ( 1 ){
                        
            ///get norms
            { inta m,space2;
                for ( space2 = 0; space2 < dim0 ; space2++)
                    if ( f2.canon[space2].body != nada)
                        if ( flagAllDim == 1 || space2 == dim[0])
                        {
                        inta kk=1;
                        for ( m = 0; m < L1; m++){
                            norm[space2][ m ] = cblas_dnrm2(M2[space2], alloyStream[space2][m],1);
#if VERBOSE
                            printf("norm %f\n", norm[space2][m]);
#endif
                            if ( norm[space2][ m ] == 0. ){
                                floata *pt =alloyStream[space2][m];
                                inta mm;
                                for ( mm = 0 ; mm < M2[space2] ; mm++)
                                    pt[mm] = cos(kk*2*pi*mm*1./M2[space2]);
                                kk++;
                            }
                            norm[space2][ m ] = cblas_dnrm2(M2[space2], alloyStream[space2][m],1);
                        }
                        }
                    
            }///end norms
        
            ///get arrays
            { inta m,n,space2,bufferDim;
                for ( space2 = 0; space2 < dim0 ; space2++)
                    if ( f2.canon[space2].body != nada)
                        if ( flagAllDim == 1 || space2 == dim[0])
                        {
                        for ( m = 0; m < L1; m++){
                            if ( norm[space2][m] == 0. )
                                printf("oops %d %d",space2,m);
                            cblas_dscal(M2[space2], 1./(norm[space2][m]),alloyStream[space2][m], 1);
                            array[space2][ m*LS1 + m ]  = 1.;
                            for ( n = 0; n < m ; n++){
                                array[space2][ n*LS1 + m ] = cblas_ddot(M2[space2], alloyStream[space2][n],1,alloyStream[space2][m],1);
                                array[space2][ m*LS1 + n ] = array[space2][ n*LS1 + m ];
                            }
                        }
                }///end all array
                    
                
                ///all  array2
                { inta rank = 0;
                    for ( space = 0; space < SPACE ; space++)
                        if ( f1.canon[space].body != nada)
                            if ( flagAllDim == 1 || space == space0)
                            {
                            for ( m = 0; m < L1; m++)
                                for ( n = 0; n < G1 ; n++){
                                    floata* bufferPointer = originStream[space][n];
                                    floata* bufferResource = pt[rank];
                                    bufferDim = M1[space];

                                    for ( space2 = 0 ; space2 < SPACE ; space2++)
                                        if ( f2.canon[space2].body != nada)
                                            if ( spatial[space][space2] )
                                                {
                                                if ( bufferDim == M2[space2] ){
                                                    array2[space][ n*LS1 + m ] = cblas_ddot(M2[space2],bufferPointer,1,alloyStream[space2][m],1);
                                                } else {
                                                    bufferDim /= M2[space2];
                                                    cblas_dgemv(CblasColMajor, CblasNoTrans, bufferDim, M2[space2], 1.,bufferPointer,bufferDim,alloyStream[space2][m],1,0.,bufferResource,1);
                                                    bufferPointer = bufferResource;
                                                    if ( bufferPointer == pt[rank])
                                                        bufferResource = ot[rank];
                                                    else
                                                        bufferResource = pt[rank];
                                                    }
                                            }
                                }
                        }
                }///end array2
            }
            ///end arrays
        
            for ( space2 = 0 ; space2 < SPACE ; space2++)
                if ( f2.canon[space2].body != nada)
                    if ( space2 != dim[0]){
                        inta m;
                        for ( m = 0; m < L1; m++){
                            norm[space2][m] = 1.;
                        }
                    }
            ///orders
            { inta space,dim2;
                dim0=0;
                for ( dim2 = 0; dim2 < SPACE ; dim2++)
                    if ( f2.canon[dim2].body != nada)
                        dim[(dim2+spaceX)%(spaces2)] = dim0++;
            
                for ( space = 0; space < SPACE ; space++)
                    if ( f1.canon[space].body != nada)
                        if ( spatial[space][dim[0]] ){
                            space0 = space;
                            break;
                        }
            }///orders

            flagAllDim = 0;
            
            
            ///APPLY ARRAYS BELOW...TO DEFINE OPERATIVE TERMS
            
            ///all track
            { inta dim2;
                for ( l =0  ; l < L1 ; l++)
                    cblas_dcopy(L1, array[dim[1]]+l*LS1,1,track+l*LS1,1);
                for ( dim2 = 2; dim2 < dim0 ; dim2++)
                    if ( f2.canon[dim[dim2]].body != nada)
                        for ( l =0  ; l < L1 ; l++)
                            cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1, 0,array[dim[dim2]]+l*LS1,1, track+l*LS1,1 );
                cblas_dcopy(LS1*LS1, track, 1, track+LS1*LS1, 1);
                for ( l = 0 ; l < L1; l++){
                    track[ l*LS1 + l ] += condition ;
                }
             }///end all track
            
            ///all guide
            { inta n ;
                for ( n = 0; n < G1*L1 ; n++)
                    guide[n] = 1;
                
                if ( cofact != NULL )
                    for ( g = 0; g < G1 ; g++)
                        for ( l = 0; l < L1 ; l++)
                            guide[g*L1+l] *= (cofact)[originIndex[g]];
            
                for ( n = 0; n < G1 ; n++){
                    for ( space = 0; space < SPACE ; space++)
                        if ( f1.canon[space].body != nada){
                            if ( space != space0 ){
                                    cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1, 0,array2[space]+n*LS1,1, guide+n*L1,1 );
                            }
                        }
                }
            }///all guide
            
            
            
            
            
            
            /// Vectors  L1 x G1
            /// list...  L1 x M2 ==   ( cross * gstream**T )
            {
                info = tdpotrf(L1, track,LS1);
                if ( info != 0 ){
#if 1
//printf("Linear Dependence failure %f\n",1./iCondition);
                    fflush(stdout);

                    printf("potrf-Linear Dependence failure %d\n",info);
                    fflush(stdout);
#endif
                    for ( space = 0; space < SPACE;space++){
                        if ( f1.canon[space].body != nada ){
                            free(originStream[space]);
                        }
                        if ( f2.canon[space].body != nada ){
                            free(alloyStream[space]);
                        }
                    }
                    free(alloyStream);
                    free(originStream);
                    free(originIndex);

                    return -1;
                }

                inta space2 ;

                inta i, ii,lll,stride,n,m;
                inta bufferDim ;
                floata* bufferPointer;
                floata* bufferResource;
//#ifdef OMP
//#pragma omp parallel for private (ii,m,n,rank)
//#endif
                
                for ( ii = 0; ii < M2[dim[0]] ; ii++)
                {
#ifdef OMP
                        rank = 0;//omp_get_thread_num();
#else
                        rank = 0;
#endif
                    ///guide * ot ::  like   (L1,G1) * (G1) => L1
                                  
                    for ( m = 0; m < L1 ; m++)
                       tracker[m] = 0.;
                        
                    for ( n = 0; n < G1 ; n++)
                         {///recipe for awesomeness...
                            ///!buffer

                            stride = 1;
                            for ( space2 = 0 ; space2 < SPACE ; space2++)
                                if ( f2.canon[space2].body != nada )
                                    if ( spatial[space0][space2]){
                                        if ( space2 == dim[0] ){
                                            break;
                                        }else {
                                            stride *= M2[space2];
                                        }
                                    }
                            ///needs to know its only occuring ii..
                             ///M1[space0] = N^3
                             ///----
                             ///first dimension, stride = 1
                             ///ii + N * a + N * N * b  =>  a + N * b
                             ///
                             ///second dimension, stride = N
                             ///a + N * ii + N * N * b  =>  a + N * b
                             ///
                             ///third dimension, stride = N*N
                             ///a + N * b + N * N * ii  =>  a + N * b

                            for ( i = 0 ; i < M1[space0]; i += stride*M2[dim[0]] ){
                                cblas_dcopy( stride , originStream[space0][n]+ii*stride+i, stride*M2[dim[0]], qt[rank]+i/M2[dim[0]], 1 );
                            }

                             bufferPointer = qt[rank];
                             bufferResource = pt[rank];
                             for ( m = 0; m < L1 ; m++){
                                bufferDim = M1[space0]/M2[dim[0]];
                                for ( space2 = 0; space2 < SPACE ; space2++)
                                        if ( f2.canon[space2].body != nada )
                                            if ( spatial[space0][space2] ){
                                                if ( dim[0] != space2 ){
                                                    if ( bufferDim == M2[space2] ){
                                                        tracker[m] += cblas_ddot(M2[space2], bufferPointer, 1, alloyStream[space2][m], 1)*(guide)[m+L1*n];
#if VERBOSE
                                                        printf("%d %d %f\n",ii,space2,tracker[m]);
#endif
                                                    } else {
                                                        bufferDim /= M2[space2] ;
                                                        cblas_dgemv(CblasColMajor, CblasNoTrans, bufferDim, M2[space2],     1.,bufferPointer,bufferDim,alloyStream[space2][m],1,0.,bufferResource,1) ;
                                                    }
                                                    bufferPointer = bufferResource ;
                                                    if ( bufferPointer == pt[rank] )
                                                        bufferResource = ot[rank] ;
                                                    else
                                                        bufferResource = pt[rank] ;
                                                }
                                            }
                            }
                         }
                    if ( ! info )
                        if (tdpotrs(L1,  1, track,LS1,  tracker ,LS1 )){
                            info = 1;
                        }
                    
                    for ( m = 0 ;m < L1 ; m++){
                        alloyStream[dim[0]][m][ii] = tracker[m];
                    }
                }///end ii
            }
                
                if ( info != 0 ){
    #if 1
                    printf("potrs-Linear Dependence failure %d\n",info);
                    fflush(stdout);
    #endif
                    for ( space = 0; space < SPACE;space++){
                        if ( f1.canon[space].body != nada ){
                            free(originStream[space]);
                        }
                        if ( f2.canon[space].body != nada ){
                            free(alloyStream[space]);
                        }
                    }
                    free(alloyStream);
                    free(originStream);
                    free(originIndex);

                    return -1;
                }

            
            
            ///at end of loop,,,
            if ( dim[0] == spaces2-1 ){
                
                ///get norms
                { inta m,space2;
                    for ( space2 = 0; space2 < dim0 ; space2++)
                        if ( f2.canon[space2].body != nada)
                            if ( flagAllDim == 1 || space2 == dim[0])
                            {
                                for ( m = 0; m < L1; m++){
                                    norm[space2][ m ] = cblas_dnrm2(M2[space2], alloyStream[space2][m],1);
                                }
                            }
                        
                }///end norms
                
                
                
                
                
                ///get sum-sqr
                { inta n,m;
                    sum2 = 0.;
                    for ( m = 0; m < L1 ; m++){
                        sum2 += norm[0][m] * norm[0][m];
                    }///end sum-sqr
                    
                    
                    ///all inner track
                    for ( m =0  ; m < L1 ; m++)
                        cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1, 0,array[dim[0]]+m*LS1,1, track+LS1*LS1+m*LS1,1 );

                    ///all inner guide
                    for ( n = 0; n < G1 ; n++ ){
                        cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1, 0,array2[space0]+n*LS1,1, guide+n*L1,1 );
                    }
                
                ///all inner guide
                ///FF
                    { inta l,ll;
                        iFF = 0.;
                        for ( l = 0; l < L1 ; l++)
                            for ( ll = 0 ; ll < L1 ; ll++)
                            {
                                iFF += norm[0][l]*(track+LS1*LS1)[ l*LS1+ll ]*norm[0][ll];
                            }
                    
#if VERBOSE
                    printf("iFF %1.15f\n",iFF);
                    fflush(stdout);
#endif
                    }
                  
                ///GF
                    { inta n,m;
                        iGF = 0.;
                        for ( n = 0 ; n < G1 ; n++ )
                            for ( m = 0; m < L1 ; m++ ){
                                iGF += (guide)[ n*L1+m ]*norm[0][m];
                            }
    #if VERBOSE
                    printf("iGF %1.15f\n",iGF);
                    fflush(stdout);
    #endif
                    }
                    
                    prev =curr;
                    curr = fabs(iGG + iFF - 2 * iGF);
                    printf("%f + %f + %f = %f\n", iGG,iFF,-2*iGF,curr);

            }///get inners
                
                
            {
                
#if VERBOSE
            printf("%d %f %d %d\n", count, curr,G1,L1);
            fflush(stdout);
#endif
            if ( fabs(curr) < target ){
#ifdef VERBOSE_ALS
                printf("Beylkin's condition number %f\n", sum2/iGG);

                printf("cycles %d\n distance-2 %1.15f \t magnitude-2 (%f) \t %d->%d\n", count,(curr),(iGG),G1,L1);
#endif
                if ( flipSignFlag )
                    tScaleOne(f2, alloy,spin, -1);
                
                for ( space = 0; space < SPACE;space++){
                    if ( f1.canon[space].body != nada ){
                        free(originStream[space]);
                    }
                    if ( f2.canon[space].body != nada ){
                        free(alloyStream[space]);
                    }
                }
                  free(alloyStream);
                  free(originStream);
                  free(originIndex);
                return (curr);
            }
                
            if ( count > maxCycle )
                {
#ifdef VERBOSE_ALS
                    printf("Beylkin's condition number %f\n", sum2/iGG);

                    printf("Failed in maxCycles %d\n distance-2 %1.15f \t magnitude-2 (%f) \t %d->%d\n", count,(curr),(iGG),G1,L1);
#endif
                    if ( flipSignFlag )
                        tScaleOne(f2, alloy,spin, -1);

                    
                    for ( space = 0; space < SPACE;space++){
                        if ( f1.canon[space].body != nada ){
                            free(originStream[space]);
                        }
                        if ( f2.canon[space].body != nada ){
                            free(alloyStream[space]);
                        }
                    }
                      free(alloyStream);
                      free(originStream);
                      free(originIndex);
                    return (curr);
                }
            }
            }//end of loop

            count++;
            flipSignFlag = 0;
            
            spaceX++;
        }
    return 0;
}


/**
 * Controls to initiate and run canonicalRankCompression
 *
 *@param spatial Defining how origin breaks down to alloy
 *@param f1          container
 *@param[in] origin the division with more canonical ranks
 *@param os         spin of origin to consider
 *@param f2         container
 *@param[out] alloy  the division with less canonical ranks, to be trained
 *@param spin      spin of alloy to consider
 *@param tolerance a number setting the absolute quality
 *@param relativeTolerance a number seting quality relative to magnitude of origin
 *@param condition Beylkin's condition (alpha)
 *@param threshold the smallest number
 *@param maxCycle the maxmium number of cycles in this routine
*/
double CanonicalRankCompression ( inta spatial[SPACE][SPACE],  double * cofact, sinc_label  f1 , division origin,inta os, sinc_label  f2 , division alloy,inta spin,  double tolerance ,  double relativeTolerance, double condition,double threshold, inta maxCycle ,double maxCondition, inta canon,inta X1 ){
    inta rank=0;
    inta L1 = canon;
        
    inta *iiii[2][2];
    inta *iii[2][2];
    floata  *me = myStreams(f1, CanonicalBuffers, 0);
    
    division G = nullName;
    inta ii,n,m,c,g,G1 = CanonicalRank(f1, origin, os);
    
    if ( ! G1 ){
    //    printf("CanonicalRankDecomposition, Origin is empty\n");
        f1.name[alloy].Current[spin] = 0;
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
                canonicalRankCompression(spatial, cofact,f1, G1,me, origin,iii[1][0][r2],iii[1][1][r2], os,f2,nRun == L1,alloy, iii[0][0][r2],iii[0][1][r2],spin,tolerance, relativeTolerance,condition,maxCondition,maxCycle);
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
