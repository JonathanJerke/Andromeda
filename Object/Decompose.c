/**
*  Decompose.c
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

#include "Decompose.h"



/**
 *Canonical Rank Decomposition or Alternating Least Squares (ALS)
 *
 *
 *
 *@param f1          container
 *@param[in] cofact allows for rescaling the origin by each canonical rank
 *@param[in] GG          internal from ASTER, will speedup processes
 *@param[in] origin the division with more canonical ranks
 *@param l1         Allowance for starting somewhere else but the first index
 *@param l2         last index of origin
 *@param os         spin of origin to consider
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
floata canonicalRankDecomposition( sinc_label  f1 , floata * cofact,inta G,floata *GG,   division origin,inta l1,inta l2,inta os, inta neo,division alloy ,inta l3 , inta l4,  inta spin ,double tolerance,double relativeTolerance, double condition,double maxCondition, inta maxCycle){
    if (l2 < l1 || l4 < l3 )
    {
        printf("indices out of order!\n");
        printf("%d %d %d %d\n", l1,l2,l3,l4);
        printf("%d %d\n", origin, alloy);
        exit(0);
    }
    if ( ! allowQ(f1.rt,blockTrainVectorsblock)){
        printf("blockTrainVectorsblock allow!\n");
        fflush(stdout);
        exit(0);
    }
    inta rank = 0;
    double xprod = 0.;
    inta xOriginIndex=0;
    assignCores(f1, 1);
    double sum2 ,iGG=0,iGF=0,iFF=0;
    inta flagAllDim,spaces = 0;
    inta g,l,count = 1;
    double target;
    floata * array[SPACE],curr=1e6,prev;
    floata * array2[SPACE],*norm[SPACE];
    floata * guide, *track;
    inta space,space0,flipSignFlag=0;
    inta L1 = l4-l3;
    inta M2[SPACE];
    inta G1 = l2-l1;
    
    if ( G1 == 0 ){
        printf("CD: empty origin %d _%d -> %d _%d\n", origin, os , alloy, spin);
        return 0;
    }
    for ( space = 0; space < SPACE ; space++)
        if ( f1.canon[space].body != nada )
            spaces++;
    
        
    length(f1,alloy,M2);
#if VERBOSE
    printf("m2 %d %d %d\n", M2[0],M2[1],M2[2]);
    printf("rank %d\n", rank);
    fflush(stdout);
#endif
    inta info;
      division canonicalStore;
    inta LS1 = L1;
    {
        for ( space = 0; space < SPACE ; space++)
        if( f1.canon[space].body != nada )
        {
            array[space] =  streams(f1, canonicalBuffers, rank , space);
            array2[space] =  array[space] + L1*L1;
            norm[space] = array2[space] + G1*L1;
        }else {
            array[space] = NULL;
            array2[space] = NULL;
            norm[space] =NULL;
        }
        guide =  myStreams(f1, guideBuffer, rank );
        track =  myStreams(f1, trackBuffer, rank );
        
        
        if (  L1*G1 + L1*L1+L1 >  part(f1,canonicalBuffers)|| G1*L1 > part(f1,guideBuffer) || L1*L1*2 > part(f1,trackBuffer)){
#if 1
            printf("canonicalRankDecomposition, mem prob with canonicalBuffers, guideBuffer, or trackBuffer\n %d %d \n", L1, G1);
            fflush(stdout);
#endif
            exit(0);
        }
        
        canonicalStore = canonicalBuffersB;
        
        if ( part(f1, canonicalStore) < L1 + G1)
        {
            printf("canonicalRankDecomposition, mem prob with canonicalBuffersB\n");
            fflush(stdout);

            exit(0);
        }
        
        
    }
    floata *pt[f1.rt->NLanes] ;
    floata *ot[f1.rt->NLanes] ;
    {
        inta r;
        for ( r = 0 ; r < f1.rt->NLanes ; r++)
        {
            pt[r] = myStreams(f1,canonicalStore , r);
            ot[r] = myStreams(f1,canonicalStore , r)+L1;
        }
    }
    
    inta dim0 = 0;
    
    floata *** originStream = malloc(SPACE * sizeof(floata**));
    floata *** alloyStream = malloc(SPACE * sizeof(floata**));
    inta * originIndex = malloc( G1* sizeof(inta));
    zero(f1, canonicalVector, rank);

    for ( space = 0; space < SPACE;space++)
        if ( f1.canon[space].body != nada ){
    {
        originStream[space] = malloc((G1+1)*sizeof(floata*));
        originStream[space][G1] = streams(f1,canonicalVector,rank,space);//

        division nIter;
        inta n,ni,nc ;
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
                originStream[space][n-l1] = streams(f1, nIter, os, space) + (n - ni)*M2[space] ;
                originIndex[n-l1] = f1.name[nIter].Begin[os]+(n - ni);
//                if  ( ! space)
//                    printf("%d%d<%d %d %d>\n",origin,alloy, n-l1, n, originIndex[n-l1]);
                if ( cofact != NULL && ! space )
                    cblas_daxpy(M2[space], cofact[originIndex[n-l1]] , originStream[space][n-l1], 1, originStream[space][G1], 1);
                else
                    cblas_daxpy(M2[space],1. , originStream[space][n-l1], 1, originStream[space][G1], 1);

#if VERBOSE
                printf("%d><%ld\n", n-l1, (originStream[space][n-l1]-streams(f1,origin,os,space)));
                fflush(stdout);
#endif
            }

        }
    
        
        {
            alloyStream[space] = malloc(L1*sizeof(floata*));

            inta n,ni,nc;
            division nIter;
            nIter = alloy;
            nc = f1.name[nIter].Current[spin];
            ni = 0;
            for ( n = 0; n < L1+l3 ; n++){
                
                while ( n >= ni+nc ){
                    ni += nc;
                        nIter = f1.name[nIter].chainNext;
                    nc = f1.name[nIter].Current[spin];

                        if ( nIter == nullName){
                            printf("somehow the pointers in Alloy lost track, you probably have too many canonical ranks\n");
                            exit(0);
                            }
                }
                if ( n >= l3){
                    alloyStream[space][n-l3] = streams(f1, nIter, spin, space) +(n - ni)*M2[space] ;
#if VERBOSE
                    printf("%d<>%ld\n", n-l3, (alloyStream[space][n-l3]-streams(f1, alloy, spin, space)));
                    fflush(stdout);
#endif
                }

            }
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
                    product *= cblas_ddot(M2[space], originStream[space][l], 1, originStream[space][ll], 1);
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
    target = max( iGG*relativeTolerance , tolerance );
    
    
            if(1){
                inta ct = 0;
                for ( space = 0; space < SPACE;space++)
                    if ( f1.canon[space].body != nada ){
                        ct++;
                        if ( neo ){
                            cblas_dcopy(M2[space], originStream[space][G1], 1, alloyStream[space][L1-1], 1);
                        }
                    }
                if ( ct == 1 ){
                    //NO SOP!
                    for ( space = 0; space < SPACE;space++)
                        if ( f1.canon[space].body != nada ){
                            free(alloyStream[space]);
                            free(originStream[space]);
                        }
                    free(alloyStream);
                    free(originStream);
                    free(originIndex);
                    ///automatically did the explicit summation.
                    
                    return 0.;
                    //return 1;
                    
                }
                ///separate output ..which means first rank is done...
            }
    {
        
        if ( L1 >= G1 ){
            for ( l = 0; l < G1 ; l++)
                for ( space = 0 ; space < SPACE ; space++)
                   if ( f1.canon[space].body != nada )
                    cblas_dcopy(M2[space], originStream[space][l],1,alloyStream[space][l],1);
#if VERBOSE
            printf("short");
#endif
            for ( space = 0; space < SPACE;space++)
                if ( f1.canon[space].body != nada ){
                    free(alloyStream[space]);
                    free(originStream[space]);
                }
            free(alloyStream);
            free(originStream);
            free(originIndex);
            return 0.;
            //return G1;
        }
//        if ( iGG < tolerance ){
//#if VERBOSE
//            printf("small\n");
//            fflush(stdout);
//#endif
//            for ( space = 0; space < SPACE;space++)
//                if ( f1.canon[space].body != nada ){
//                    free(alloyStream[space]);
//                    free(originStream[space]);
//                }
//            free(alloyStream);
//            free(originStream);
//            free(originIndex);
//
//            return 0;
//        }
    }
    
    inta dim[SPACE],ll;

    {
        
        space0 = 0;
        
        dim0=0;
        for ( space = 0; space < SPACE ; space++)
            if ( f1.canon[space].body != nada)
                dim[(space0+space)%(spaces)] = dim0++;
        
        
        {
            inta m,space;

        
        ///get norms
        for ( space = 0; space < dim0 ; space++)
            if ( f1.canon[space].body != nada){
                inta kk=1;
                for ( m = 0; m < L1; m++){
                    norm[space][ m ] = cblas_dnrm2(M2[space], alloyStream[space][m],1);
                    if ( norm[space][m] == 0. ){
                        floata *pt =alloyStream[space][m];
                        inta mm;
                        for ( mm = 0 ; mm < M2[space] ; mm++)
                            pt[mm] = cos(kk*2*pi*mm*1./M2[space]);
                        kk++;
                    }
                    norm[space][ m ] = cblas_dnrm2(M2[space], alloyStream[space][m],1);
                }
            }
            
        }
        
        ///get inners
        for ( space = 0; space < dim0 ; space++)
            if ( f1.canon[space].body != nada){
            
                inta m,n;
                for ( m = 0; m < L1; m++){
#if VERBOSE
                    for ( l= 0 ; l < M2[space] ; l++)
                        if (isnan (alloyStream[space][m][l] ) || isinf(alloyStream[space][m][l]))
                            printf("alloy error\n");
#endif
                    if ( norm[space][m] == 0. )
                        printf("oops %d %d",space,m);
                    cblas_dscal(M2[space], 1./(norm[space][m]),alloyStream[space][m], 1);
                    array[space][ m*LS1 + m ]  = 1.;
                    for ( n = 0; n < m ; n++){
                        array[space][ n*LS1 + m ] = cblas_ddot(M2[space], alloyStream[space][n],1,alloyStream[space][m],1);
                        array[space][ m*LS1 + n ] = array[space][ n*LS1 + m ];
                    }
                    for ( n = 0; n < G1 ; n++){
#if VERBOSE
                            for ( l= 0 ; l < M2[space] ; l++)
                                if (isnan (originStream[space][n][l] ) || isinf(originStream[space][n][l]))
                                    printf("origin error\n");
#endif
                        array2[space][ n*LS1 + m ] = cblas_ddot(M2[space],originStream[space][n],1,alloyStream[space][m],1);
                        
                    }
                }
            }
        }
            
        space0 = 1;
        while (1 ){
            ///all computed inner products
            ///            printf("dim[0]=%d\n", dim[0]);
            dim0 = 0;
            for ( space = 0; space < SPACE ; space++)
                if ( f1.canon[space].body != nada)
                    dim[(space0+space)%(spaces)] = dim0++;
                       
                       
                   
            for ( l =0  ; l < L1 ; l++)
                cblas_dcopy(L1, array[dim[1]]+l*LS1,1,track+l*LS1,1);
            for ( space = 2; space < dim0 ; space++)
                if ( f1.canon[dim[space]].body != nada)
                    for ( l =0  ; l < L1 ; l++)
                        cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1, 0,array[dim[space]]+l*LS1,1, track+l*LS1,1 );
           
            cblas_dcopy(LS1*LS1, track, 1, track+LS1*LS1, 1);
            for ( l = 0 ; l < L1; l++)
                    track[ l*LS1 + l ] += condition ;

//            if(0){
//                double ev[LS1];
//
//                printf("\n{");
//                for ( l = 0; l < L1*L1 ; l++)
//                {
//                        printf("%f,",track[l]);
//                }
//                printf("}\n");
//
//
//                tdsyev(0, f1, 'N', L1, track, LS1, ev);
//
//
//
////                cblas_dcopy(LS1*LS1, track+LS1*LS1, 1, track, 1);
////                for ( l = 0 ; l < L1; l++)
////                        track[ l*LS1 + l ] += condition ;
////
////                iCondition = tdpocon(0, f1, L1, track, LS1);
////                printf("c %f %f\n", ev[L1-1]/ev[0],iCondition);
//                iCondition = 1;
//                if ( 1./maxCondition > ev[0]/ev[L1-1] || isnan(ev[0]/ev[L1-1])||isinf(ev[0]/ev[L1-1])  ){
//                    #if 1
//                    printf("Linear Dependence failure %f\n",1./iCondition);
//                                        fflush(stdout);
//                    #endif
//                                        for ( space = 0 ; space < SPACE ; space++)
//                                            if ( f1.canon[space].body != nada ){
//                                                free(alloyStream[space]);
//                                                free(originStream[space]);
//                                            }
//                                        free(alloyStream);
//                                        free(originStream);
//                                        free(originIndex);
//
//                                        return -1;
//                }
//                cblas_dcopy(LS1*LS1, track+LS1*LS1, 1, track, 1);
//                for ( l = 0 ; l < L1; l++ )
//                    track[ l*LS1 + l ] += condition ;
//            }
            
            
            
            for ( l = 0; l < G1 ; l++)
                cblas_dcopy(L1, array2[dim[1]]+l*LS1,1,guide+l*L1,1);
            for ( space = 2; space < dim0 ; space++)
                if ( f1.canon[dim[space]].body != nada){
                    for ( l = 0; l < G1 ; l++)
                        cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1, 0,array2[dim[space]]+l*LS1,1, guide+l*L1,1 );
                    
                }
            

            if ( cofact != NULL )
                for ( g = 0; g < G1 ; g++)
                    for ( l = 0; l < L1 ; l++)
                        guide[g*L1+l] *= (cofact)[originIndex[g]];
            
            /// Vectors  L1 x G1
            // list...  L1 x M2 ==   ( cross * gstream**T )
            
            {
                info = tdpotrf(L1, track,LS1);
                if ( info != 0 ){
#if 1
//printf("Linear Dependence failure %f\n",1./iCondition);
                    fflush(stdout);

                    printf("potrf-Linear Dependence failure %d\n",info);
                    fflush(stdout);
#endif
                    for ( space = 0 ; space < SPACE ; space++)
                        if ( f1.canon[space].body != nada ){
                            free(alloyStream[space]);
                            free(originStream[space]);
                        }
                    free(alloyStream);
                    free(originStream);
                    free(originIndex);

                    return -1;
                }

                inta lll;
#ifdef OMP
#pragma omp parallel for private (ll,lll,rank)
#endif
                for ( ll = 0; ll < M2[dim[0]] ; ll++){
                    #ifdef OMP
                            rank = omp_get_thread_num();
                    #else
                            rank = 0;
                    #endif
                    for ( lll = 0 ;lll < G1 ; lll++){
                        ot[rank][lll] = originStream[dim[0]][lll][ll];
                    }
                    ///guide * ot ::  like   (L1,G1) * (G1) => L1
                                        
                    cblas_dgemv(CblasColMajor, CblasNoTrans,L1,G1,1.,guide,L1,ot[rank],1, 0., pt[rank],1);
                    if ( ! info )
                    if (tdpotrs(L1,  1, track,LS1,  pt[rank],LS1 )){
                        info = 1;
                    }
                    for ( lll = 0 ;lll < L1 ; lll++){
                        alloyStream[dim[0]][lll][ll] = pt[rank][lll];
                    }
                }
                
                if ( info != 0 ){
    #if 1
             //       printf("Linear Dependence failure %f\n",1./iCondition);
             //       fflush(stdout);

                    printf("potrs-Linear Dependence failure %d\n",info);
                    fflush(stdout);

    #endif
                    for ( space = 0 ; space < SPACE ; space++)
                        if ( f1.canon[space].body != nada ){
                            free(alloyStream[space]);
                            free(originStream[space]);
                        }
                    free(alloyStream);
                    free(originStream);
                    free(originIndex);

                    return -1;
                }

            }
            

            if ( dim[0] == spaces-1 ){
                
                        
                
                
                double prod;
                                  inta ll;
                                  sum2 = 0.;
                                  for ( ll = 0; ll < LS1 ; ll++){
                                      prod = 1.;
                                      for ( space = 0 ; space < SPACE ; space++)
                                          if ( f1.canon[space].body != nada )
                                              prod *= norm[space][ll] * norm[space][ll];
                                      sum2 += prod;
                                  }
                
                
                for ( l = 0; l < G1 ; l++){
                    cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1, 0,array2[dim[0]]+l*LS1,1, guide+l*L1,1 );
                }
                for ( l =0  ; l < L1 ; l++)
                    cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1, 0,array[dim[0]]+l*LS1,1, track+LS1*LS1+l*LS1,1 );

                iFF = 0.; iGF = 0.;
                
                for ( l = 0; l < L1 ; l++)
                    for ( ll = 0 ; ll < L1 ; ll++)
                            iFF += (norm[0][l]*(track+LS1*LS1)[ l*LS1+ll]*norm[0][ll]);



#if VERBOSE
                printf("iFF %f\n",iFF);
                fflush(stdout);
#endif
                for ( l = 0 ; l < G1 ; l++)
                    for ( ll = 0; ll < L1 ; ll++){
                        prod = 1.;
                        if ( cofact != NULL )
                            prod *= (cofact)[originIndex[l]];

                        iGF += prod * (guide)[ l*L1+ll] * norm[0][ll] ;

        
                    }
#if VERBOSE
                printf("iGF %f\n",iGF);
                fflush(stdout);
#endif
                prev =curr;
                if (fabs(iGG+iFF - 2 * iGF) < fabs(iGG+iFF + 2 * iGF) ){
                    curr = fabs(iGG+iFF - 2 * iGF);
                    flipSignFlag = 0;
                    }else
                    {
                    curr = fabs(iGG+iFF + 2 * iGF);
                    flipSignFlag = 1;
                }

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
                    tScaleOne(f1, alloy,spin, -1);
                
                for ( space = 0; space < SPACE;space++)
                    if ( f1.canon[space].body != nada ){
                        free(alloyStream[space]);
                        free(originStream[space]);
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
                        tScaleOne(f1, alloy,spin, -1);

                    
                    for ( space = 0; space < SPACE;space++)
                        if ( f1.canon[space].body != nada ){
                            free(alloyStream[space]);
                            free(originStream[space]);
                      }
                    free(alloyStream);
                    free(originStream);
                    free(originIndex);
                    return (curr);
                }
            }
            }
                        
            flagAllDim = 0;
            count++;
        
            
            {
                {
                    inta m,space;
                    
                    for ( space = 0; space < dim0 ; space++)
                        if ( f1.canon[space].body != nada){
                            for ( m = 0; m < L1; m++){
                                norm[space][ m ] = cblas_dnrm2(M2[space], alloyStream[space][m],1);
                            }
                        }
                }
                
                
                            
                inta space ;
                for ( space = 0 ; space < SPACE ; space++)
                    if( f1.canon[space].body != nada ){

                    if ( flagAllDim == 1 || space == dim[0])
                {
                    
                
                    inta m,n;
                    
                    
    #ifdef OMP
        #pragma omp parallel for private (m)
                    #endif
                    for ( m = 0; m < L1; m++){
                        if ( flipSignFlag && space == dim[0] )
                            cblas_dscal(M2[space], -1./(norm[space][m]),alloyStream[space][m], 1);
                        else
                            cblas_dscal(M2[space], 1./(norm[space][m]),alloyStream[space][m], 1);
                    }
                    
#ifdef OMP
#pragma omp parallel for private (m,n)
#endif
                for ( m = 0; m < L1; m++){
                        array[space][ m*LS1 + m ]  = 1.;
                        for ( n = 0; n < m ; n++){
                            array[space][ n*LS1 + m ] = cblas_ddot(M2[space], alloyStream[space][n],1,alloyStream[space][m],1);
                            array[space][ m*LS1 + n ] = array[space][ n*LS1 + m ];
                        }
                        for ( n = 0; n < G1 ; n++){
                            array2[space][ n*LS1 + m ] = cblas_ddot(M2[space],originStream[space][n],1,alloyStream[space][m],1);
                       }
                    }
                }else {
                    inta m;
                    for ( m = 0; m < L1; m++){
                        norm[space][m] = 1.;
                    }
                }
                    }
                flipSignFlag = 0;

            }
            
            
            

             space0++;
            
        }
    return 0;
}


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
    
    inta L1;
    if ( X1 == 0 )
      L1 = canon;
    else
      L1  = CanonicalRank(f0, alloy, spin);
    
    inta GG1;
    if ( L1 - X1 <= 0 ){
        X1 = L1-1;
    }
        
    inta *iiii[2][2];
    inta *iii[2][2];
    
    if ( ! G1 ){
        //    printf("CanonicalRankDecomposition, Origin is empty\n");
            f0.name[alloy].Current[spin] = 0;
        return 0;
    }
    if ( (G1 <= L1) && !X1  ){
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
    c1.i.Lambda = 1;
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
    
//    #ifdef OMP
//    #pragma omp parallel for private (space,rank) schedule(dynamic,1)
//    #endif
    for ( space = 0; space < SPACE ; space++)
        if ( f0.canon[space].body != nada){
//            #ifdef OMP
//                    rank = omp_get_thread_num();
//            #else
            rank = 0;
//            #endif
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

    floata * me;
    me = myStreams(F1.f, CanonicalBuffers, 0);
    inta M2[SPACE];
    length(F1.f, totalVector, M2);

    {

        if ( part(F1.f,CanonicalBuffers) < G1*G1 ){
            printf("CanonicalBuffers1 \n");
            exit(0);
        }
           ///RENAME
           inta spin = 0,os = 0;
           zero(F1.f,eigenVectors,spin);

           ///RENAME END

    for ( g = 0 ; g < G1 ; g++){
        if ( G== nullName )
            G = anotherLabel(&F1.f, 0, nada);
        else
            anotherLabel(&F1.f, 0, nada);///will be in place
        F1.f.name[G+g].name = totalVector;
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
    }
    {
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
                        cblas_dcopy(M2[space], streams(F1.f,totalVector,os,space)+ii*M2[space], 1, streams(F1.f,totalVector,os,space)+iv*M2[space], 1);
                    }
            }
            iv++;
        }
    GG1 = G1;
    G1 = iv;
    
    F1.f.name[totalVector].Current[os] = G1;
    if ( (G1 <= L1) && !X1 ){
        L1 = G1;
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
        return 0;
    }else if ( (G1 < L1) && X1 ){
        ///this alternative should never happen ,because origin is an output from the same procedure (w/ X1=0), now its just conditiionoing
    }

        
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

                    F1.f.name[totalVector].Current[os] = G1;
                    F1.f.name[eigenVectors].Current[spin] = L1;
                    inta nRun = L1;
                    
                    inta run,r2;
                    
                    while ( nRun > 0 ){
                        {
                            for ( r2 = 0; r2 < nRun; r2++ ){
                                canonicalRankDecomposition(F1.f, coeff, G1,me, totalVector,iii[1][0][r2],iii[1][1][r2], os,nRun == L1,eigenVectors, iii[0][0][r2],iii[0][1][r2],spin,tolerance, relativeTolerance,condition,maxCondition,maxCycle);
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
            F1.f.name[totalVector].Current[os] = G1;
            F1.f.name[eigenVectors].Current[spin] = L1;

            tEqua(F1.f, eigenVectors, spin, totalVector, os);

            inta i,x1;
            floata target,Glen = 0.,curr;
            for ( i = 0; i < G1*G1 ; i++)
                Glen += me[i];
            
            for ( x1 = 1 ; x1 <= X1 ; x1++){
#ifdef VERBOSE_ALS
                printf("attempting to decompose to %d \n", L1-1);
#endif
                curr = canonicalRankDecomposition(F1.f, coeff, G1, me, totalVector, 0, G1, os, 0, eigenVectors, 0, L1-1, spin, tolerance, relativeTolerance, condition, maxCondition, maxCycle);
                target = max( tolerance , Glen*relativeTolerance );
                if ( curr > target){
#ifdef VERBOSE_ALS
                    printf("stop decomposing at %d above %f\n", L1, target);
#endif
                    canonicalRankDecomposition(F1.f, coeff, G1, me, totalVector, 0, G1, os, 0, eigenVectors, 0, L1, spin, tolerance, relativeTolerance, condition, maxCondition, maxCycle);
                    break;
                }
                L1--;
                F1.f.name[eigenVectors].Current[spin] = L1;
#ifdef VERBOSE_ALS
                printf("decomposed to %d \t %f \t %f\n", L1,curr,Glen);
#endif
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
