/**
 *  mAls.c
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

#include "mAls.h"


/**
 *Canonical Rank Decomposition or Alternating Least Squares (ALS)
 *
 *
 *
 *Mean to be run for fixed input lengths from ASTER
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
inta canonicalRankDecomposition( sinc_label  f1 , floata * cofact,inta G,floata *GG,   division origin,inta l1,inta l2,inta os, inta neo,division alloy ,inta l3 , inta l4,  inta spin ,double tolerance,double relativeTolerance, double condition,double maxCondition, inta maxCycle){
    if (l2 < l1 || l4 < l3 )
    {
        printf("indices out of order!\n");
        printf("%d %d %d %d\n", l1,l2,l3,l4);
        printf("%d %d\n", origin, alloy);
        exit(0);
    }
    double iCondition;
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
#if VERBOSE
    printf("G %d L %d\n", G1,L1 );
    fflush(stdout);
#endif
    
    
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
                inta ct = 0,ll;
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
                    return 1;
                    
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
            return G1;
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
                for ( m = 0; m < L1; m++){
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

            if(0){
                double ev[LS1];

                printf("\n{");
                for ( l = 0; l < L1*L1 ; l++)
                {
                        printf("%f,",track[l]);
                }
                printf("}\n");

                
                tdsyev(0, f1, 'N', L1, track, LS1, ev);
                
                
                
//                cblas_dcopy(LS1*LS1, track+LS1*LS1, 1, track, 1);
//                for ( l = 0 ; l < L1; l++)
//                        track[ l*LS1 + l ] += condition ;
//
//                iCondition = tdpocon(0, f1, L1, track, LS1);
//                printf("c %f %f\n", ev[L1-1]/ev[0],iCondition);
                iCondition = 1;
                if ( 1./maxCondition > ev[0]/ev[L1-1] || isnan(ev[0]/ev[L1-1])||isinf(ev[0]/ev[L1-1])  ){
                    #if 1
                    printf("Linear Dependence failure %f\n",1./iCondition);
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
                cblas_dcopy(LS1*LS1, track+LS1*LS1, 1, track, 1);
                for ( l = 0 ; l < L1; l++ )
                    track[ l*LS1 + l ] += condition ;
            }
            
            
            
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
                    printf("Linear Dependence failure %f\n",1./iCondition);
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
                    printf("Linear Dependence failure %f\n",1./iCondition);
                    fflush(stdout);

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
#if 1
                printf("Beylkin's condition number %f\n", sum2/iGG);

                printf("cycles %d\n distance %1.15f \t magnitude (%f) \t %d->%d\n", count,sqrt(curr),sqrt(iGG),G1,L1);
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
                
                
                
                
                
                
                return L1;
            }
                
            if ( count > maxCycle )
                {
#if 1
                    printf("Beylkin's condition number %f\n", sum2/iGG);

                    printf("Failed in maxCycles %d\n distance %1.15f \t magnitude (%f) \t %d->%d\n", count,sqrt(curr),sqrt(iGG),G1,L1);
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
                    return L1;
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
 *The printing output in D.kry
 *
 *@param c  The general parameters
 *@param f1 The container
 *@param Hamiltonian The all links will be printed
 *@param vector division of the input vector
 */
double printExpectationValues (  calculation *c,   sinc_label  f1 ,  division Hamiltonian  ,   division vector){
    mea totx,me,ov,edge;
    inta space,ed,i,spacer;
      bodyType body;
    if ( !allowQ(f1.rt, blockPrintStuffblock))
        return 0.;
    

    totx = 0.;
    printf("\n======Expectation========\n");
    if (CanonicalRank(f1, vector, 1))
        printf("\t%d (%d + i%d)\n", vector, CanonicalRank(f1, vector, 0), CanonicalRank(f1, vector, 1));
    else
        printf("\t(%d)\t%d\n",  CanonicalRank(f1, vector, 0),vector-eigenVectors);
    ov = pMatrixElement( f1, vector, 0, nullOverlap, 0, vector,0);
    printf("------Edges------\n\n");
    division header = anotherLabel(&f1, 0, nada);
    division eik    = anotherLabel(&f1, 0, nada);
    division mem    = anotherLabel(&f1, 0, one);
    char* desc [] = {"fourier","negative","positive"};
    for ( space = 0; space < SPACE ; space++){
        for (body = one ; body <=  f1.canon[space].body ; body++ )
            for ( ed = 0 ; ed < 3 ; ed++){
                f1.name[header].species = matrix;
                f1.name[header].chainNext = eik;
                f1.name[eik].loopNext = mem;
                f1.name[eik].species = eikon;
                f1.name[mem].species = eikonOuter;
                f1.name[mem].Current[0] = 1;
                f1.name[eik].space[space].body = one;

                for ( spacer = 0; spacer < SPACE ;spacer++)
                    f1.name[eik].space[spacer].block = id0;
                f1.name[eik].space[space].block = (blockType)body;
                zero(f1 , mem,0);
            switch ( ed ){
                case 0:
                    for ( i = 0; i < vector1Len(f1, space) ; i++)
                        streams(f1,mem,0,space)[i] = sign(i)*1./vector1Len(f1, space);
                    break;
                case 1:
                    streams(f1,mem,0,space)[0] = 1;
                    break;
                case 2:
                    streams(f1,mem,0,space)[vector1Len(f1, space)-1] = 1;
                    break;
            }
                edge = pMatrixElement(f1, vector, 0, header, 0, vector,0);
                printf("\t%d\t%d\t%s \t%1.15f\n", space+1,body,desc[ed], edge/ov);
            }
    }
        
    printf("------Terms------\n");

    division op = defSpiralMatrix(&f1, Hamiltonian);
    inta o,scr=0,cr =0;;
    inta terms = 0,oo=0;
    
    
    for (o = 0; f1.name[op+o].species == matrix ; o++)
    {
        
        while ( terms <= o && oo < c->i.termNumber ){
            if ( c->i.terms[oo].headFlag ) {
                terms++;
            }
            printf("%s ", c->i.terms[oo].desc);
            oo++;
        }
        
        me = pMatrixElement( f1, vector, 0, op+o, 0, vector,0);
        cr = CanonicalOperator(f1, f1.name[op+o].name, 0);
        scr += cr;
        printf("\t(%d)\t%6.6f\n", cr,creal(me/ov));
        fflush(stdout);
        totx += me/ov;
    }
    f1.name[vector].value.value = totx;
   {
#ifdef COMPLEXME
            printf("sum\t(%d)\t%6.6f\tI %6.6f\n",scr,creal(totx),cimag(totx));
#else
            printf("sum\t(%d)\t%6.6f\n",scr,(totx));
#endif
    }
        
    return totx;
}


/**
 *The matrix element calculator using multiply
 *
 *@param rank      OMP rank
 *@param f1          container
 *@param bra       division on left of operator
 *@param bspin  bra's spin
 *@param mat       A single link or Hamilonian Term
 *@param mspin mat's spin
 *@param ket        division on right of operator
 *@param kspin  ket's spin
 */
double tMatrixElements ( inta rank,  sinc_label  f1 ,   division bra, inta bspin,  division mat, inta mspin,  division ket, inta kspin ){
    
      division holder ;
    inta holderRank, holderSpin;
    inta k,l,e,space;
    double prod,ME=0;
    inta ca;
    
    if ( rank )
        if ( ! allowQ(f1.rt, blockParallelMatrixElementblock)){
            printf("blockParallelMatrixElementblock Allow!\n");
            fflush(stdout);
            exit(0);
        }
    
    
    if ( mat == nullName || f1.name[mat].name == nullName)
        return 0.;
    
    if (mat == nullOverlap  )
        ca = 1;
    else if ( f1.name[mat].species == matrix )
        ca = CanonicalOperator(f1,f1.name[mat].name, mspin);
    else
        return 0.;
    
                for ( k = 0 ; k < CanonicalRank(f1, ket, kspin);k++){
                    for ( l = 0 ; l < ca;l++){
                        if ( mat == nullOverlap ){
                            holder = ket;
                            holderRank = k;
                            holderSpin = kspin;
                        }
                        else {
                            tHX(rank, f1, f1.name[mat].name, l, mspin, 1., ket, k, kspin,canonicalmeVector, 0, rank);
                            
                            holder = canonicalmeVector;
                            holderRank = 0;
                            holderSpin = rank;
                            
                        }
                        for ( e = 0 ; e < CanonicalRank(f1, bra, bspin);e++){
                            prod = 1;
                            for ( space = 0 ; space < SPACE ; space++)
                                if ( f1.canon[space].body != nada){
                                    prod *= tDOT(rank, f1,space,CDT, bra, e, bspin,CDT, holder, holderRank, holderSpin);
                                }
                            ME += prod;
                        }
                    }
                }
        return ME;
}
            
        


/**
 *general matrix - vector multiply per dimension
 *
 *@param rank      OMP rank
 *@param f1          container
 *@param space   in canon which has SPACE components
 *@param equals output division
 *@param espin  output spin
 *@param left  matrix
 *@param l  one of the matrices canonical ranks indexed
 *@param lspin matrix spin
 *@param right input division
 *@param r one of the input canonical ranks indexed
 *@param rspin input spin
 */
inta tGEMV (inta rank,    sinc_label  f1, inta space,   division equals, inta e, inta espin,  division left,inta l,inta lspin,   division right,inta r, inta rspin ){
    if ( header(f1, left ) != header(f1, right ) ){
        printf("Two Head types GEMV\n %d %d %d %d %d %d",equals,header(f1, equals ) ,left,header(f1, left ) ,right,header(f1, right ) );
        exit(1);
    }
    if ( rank ) {
        if ( ! allowQ(f1.rt, blockParallelMultiplyblock)){
            printf("blockParallelMultiplyblock allow!\n");
            fflush(stdout);
            exit(0);
        }
    }
    
    bodyType bd = Bodies(f1, right,space);
    division inT,outT;
    inta inR,outR,inS,outS;
    f1.name[canonicalmvVector].Current[rank] = 0;
    f1.name[canonicalmv2Vector].Current[rank] = 1;
    f1.name[canonicalmv3Vector].Current[rank] = 0;
    if ( species(f1,left )== eikon && species(f1, right ) == vector ){
        char  in,out;
            in = 1;
            out = 1;
        
        inT = right;
        inR = r;
        inS = rspin;
        
        outT = equals;
        outR = e;
        outS = espin;
//        if ( in != 1 ){
//           tPermuteOne(rank, f1, space, in, right, r, rspin, canonicalmvVector,0, rank);
//           inT = canonicalmvVector;
//           inR = 0;
//           inS = rank;
//        }
//        if (out != 1 ){
//            outT = canonicalmv2Vector;
//            outR = 0;
//            outS = rank;
//        }
        
          division midT = canonicalmv3Vector;
        inta midR = 0;
        inta midS = rank;
        
          division laterT = canonicalmv3Vector;
        inta laterR = 1;
        inta laterS = rank;

        if ( species(f1,left) == eikon ){
                //NOT id0
                //only eikons show multiple SOP structures per output SOP canonical rank
            
            
                inta N1 = outerVectorLen(f1, one,space);
                inta N2 = vectorLen(f1, space);
                
                inta i;
                double * midP = streams(f1, midT, midS,space)+midR*N2;
                double * laterP = streams(f1, laterT, laterS,space)+laterR*N2;
                double * inP  = streams(f1, inT, inS,space)+inR*N2;

                double * outP = streams(f1, outT, outS,space)+outR*N2;
            for ( i = 0 ; i < N2 ; i++)
                outP[i] = 0.;

#if VERBOSE
            printf("in %f %d\n", cblas_dnrm2(N2, inP, 1),N1);
#endif
                  division su = f1.name[left].loopNext;//direct summation per component!
                inta timer = 0,xlxl=0;
                while ( su != nullName ){
                    for ( i = 0 ; i < N2 ; i++){
                        laterP[i] = 0.;
                        midP[i] = 0.;
                    }
                    switch(species(f1,su)){
                        case eikonDiagonal:
                        case eikonKinetic:
                        case eikonConstant:
                        case eikonDeriv:
                        case eikonLinear:
                        case eikonSpring:
                        case eikonElement:
                        case eikonOuter:
                            xlxl = 1;
                            break;
                        case eikonSemiDiagonal:
                            switch(Bodies(f1, left,space)){
                                case one:
                                    xlxl = 0;
                                    break;
                                case two:
                                    xlxl = 4;
                                    break;
                            default:
                                break;

                            }
                            break;
                        case eikonOffDiagonal:
                            switch(Bodies(f1, left,space)){
                                case one:
                                    xlxl = 2;
                                    break;
                                case two:
                                    xlxl = 4;
                                    break;
                               default:
                                   break;

                            }
                            break;
                        default:
                            break;

                    }
                    for ( timer = 0 ; timer < xlxl ; timer++){
                        double flow = 1.;
                         ///the notion is to buffer on mid and accumlate on out
                        inta N3 = alloc(f1, su, space);
                        double * suP  = streams(f1, su , lspin, space)+l*N3;
                    if ( Bodies ( f1,left,space) == one ){
                         if ( species(f1,su) == eikonDeriv){
                            flow *= *suP;
                            topezOp(0,f1.canon[space].particle[f1.name[left].space[space].block].lattice, bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP,1, laterP);
                        } else
                        if ( species(f1,su) == eikonKinetic ){
                            flow *= *suP;
                            topezOp(0, f1.canon[space].particle[f1.name[left].space[space].block].lattice, bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP,2, laterP);
                        }
                        else if ( species(f1,su) == eikonConstant){
                            ///action can happen!
                            flow *= *suP;
                            topezOp(0,0, bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP,0, laterP);
                        }
                        else if ( species(f1,su) == eikonLinear){
                            flow *= *suP;
                            ///relative to grid only
                            floata center = (f1.canon[space].count1Basis-1)*f1.canon[space].particle[f1.name[left].space[space].block].lattice;
                            topezOp(center, f1.canon[space].particle[f1.name[left].space[space].block].lattice, bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP,-1, laterP);
                        }
                        else if ( species(f1,su) == eikonSpring){
                            flow *= *suP;
                            ///relative to grid only
                            floata center = 0.5* (f1.canon[space].count1Basis-1)*f1.canon[space].particle[f1.name[left].space[space].block].lattice;
                            topezOp(center, f1.canon[space].particle[f1.name[left].space[space].block].lattice, bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP,-2, laterP);

                        }
                        else if ( species(f1,su) == eikonElement){
                            flow *= *suP;
                            for ( i = 0 ;i < N2 ; i++)
                                laterP[i] = 0.;

                            ///not a spatial state, this is an electronic level... no permutations here...just dropout a vector.

                            laterP[f1.name[left].space[space].bra] = inP[f1.name[left].space[space].ket];

                        }
                        else if ( species(f1,su) == eikonOuter){
                            inta n1[MAXBODY] ,ii,i;
                            double su;
                            switch(f1.name[left].space[space].block){
                                case tv1:
                                    n1[0] = 1;
                                    n1[1] = N1;
                                    n1[2] = N1*N1;
                                    break;
                                case tv2 :
                                    n1[1] = 1;
                                    n1[0] = N1;
                                    n1[2] = N1*N1;
                                    break;
                                case tv3 :
                                    n1[1] = 1;
                                    n1[2] = N1;
                                    n1[0] = N1*N1;
                                    break;
                                default:
                                    break;
                            }
                            for ( i = 0 ;i < N2 ; i++)
                                laterP[i] = 0.;
                            switch ( bd ){
                                case one:
                                    cblas_dcopy(N1, suP, 1, laterP, 1);
                                    flow *= cblas_ddot(N1, inP, 1, suP, 1);
                                    break;
                                case two:
                                    su = 0.;
                                    for ( i = 0 ; i < N1 ; i++){
                                        cblas_dcopy(N1, suP, 1, laterP+i*n1[1], n1[0]);
                                        su += cblas_ddot(N1, inP+i*n1[1], n1[0], suP, 1);
                                    }
                                    flow *= su;
                                    break;
                                case three:
                                    su = 0.;
                                    for ( i = 0 ; i < N1 ; i++)
                                        for ( ii = 0 ; ii < N1 ; ii++){
                                            cblas_dcopy(N1, suP, 1, laterP+i*n1[1]+ii*n1[2], n1[0]);
                                            su += cblas_ddot(N1, inP+i*n1[1]+ii*n1[2], n1[0], suP, 1);
                                        }
                                    flow *= su;
                                    break;
                                default:
                                    break;

                            }
                        }

                        else
                        if ( species(f1,su) == eikonDiagonal ){
                            flow *= 1.;
                            diagonalOp(bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP, suP,laterP);
                        }else if ( species(f1,su) == eikonOffDiagonal ){
                            if ( timer == 0 ){
                                flow *= -1;
                                topezOp(0,1., bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP,1, midP);
                                diagonalOp(bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,midP, suP,laterP);
                            }
                            else if ( timer == 1 ){
                                flow *= 1;
                                diagonalOp(bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP, suP,midP);
                                topezOp(0,1., bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,midP, 1,laterP);
                            }
                    }
                }
                   else  if ( Bodies ( f1,left,space) == two ){
                       if ( species(f1,su) == eikonDiagonal ){
                            flow *= 1;
                             diagonalOp(bd,f1.name[left].space[space].act,e12, f1.name[left].space[space].block,N1,inP, suP,laterP);
                        }else if ( species(f1,su) == eikonSemiDiagonal ){

                            if ( timer == 0 ){
                                flow *= -1;

                                topezOp(0, 1, bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP, 1,midP);
                                diagonalOp(bd,f1.name[left].space[space].act,e12, f1.name[left].space[space].block,N1,midP, suP,laterP);

                           }

                            else if ( timer == 1 ){
                                flow *= 1 ;
                                topezOp(0,1, bd,f1.name[left].space[space].act,tv2, f1.name[left].space[space].block,N1,inP, 1,midP);
                                diagonalOp(bd,f1.name[left].space[space].act,e12, f1.name[left].space[space].block,N1,midP, suP,laterP);
                                
                            } else
                             if ( timer == 2 ){
                                flow *= 1;
                                 
                                 diagonalOp(bd,f1.name[left].space[space].act,e12, f1.name[left].space[space].block,N1,inP, suP,midP);
                                 topezOp(0,1, bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,midP, 1,laterP);

                                 
                                 
                            }
                        
                             else if ( timer == 3 ){
                                 flow *= -1 ;
                                 
                                 diagonalOp(bd,f1.name[left].space[space].act,e12, f1.name[left].space[space].block,N1,inP, suP,midP);
                                 topezOp(0,1, bd,f1.name[left].space[space].act,tv2, f1.name[left].space[space].block,N1,midP, 1,laterP);

                             }
                        }else
                        if ( species(f1,su) == eikonOffDiagonal ) {
                            if ( timer == 0 ){
                            flow *= 1;
                            //A
                                diagonalOp(bd,f1.name[left].space[space].act,e12, f1.name[left].space[space].block,N1,inP, suP,laterP);
                                topezOp(0,1., bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,laterP, 1,midP);
                                topezOp(0,1., bd,f1.name[left].space[space].act,tv2, f1.name[left].space[space].block,N1,midP, 1,laterP);

                            }
                        
                            else if ( timer == 1 ){
                                flow *= 1 ;
                                topezOp(0,1., bd,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1,inP, 1,laterP);
                                topezOp(0,1., bd,f1.name[left].space[space].act,tv2, f1.name[left].space[space].block,N1,laterP, 1,midP);
                                diagonalOp(bd,f1.name[left].space[space].act,e12, f1.name[left].space[space].block,N1,midP, suP,laterP);

                            }
                            else if ( timer == 2){
                                flow *= -1 ;
                                topezOp(0,1., bd,f1.name[left].space[space].act       ,tv1, f1.name[left].space[space].block,N1,inP, 1,laterP);
                                diagonalOp(bd,f1.name[left].space[space].act    ,e12, f1.name[left].space[space].block,N1,laterP, suP,midP);
                                topezOp(0,1., bd         ,f1.name[left].space[space].act,tv2, f1.name[left].space[space].block,N1,midP, 1,laterP);

                            } else if (timer ==3 ){
                                flow *= -1;
                                topezOp(0,1., bd,f1.name[left].space[space].act       ,tv2, f1.name[left].space[space].block,N1                                  ,inP,   1,laterP);
                                diagonalOp(bd,f1.name[left].space[space].act    ,e12, f1.name[left].space[space].block,N1                                  ,laterP, suP,midP);
                                topezOp(0,1., bd         ,f1.name[left].space[space].act,tv1, f1.name[left].space[space].block,N1  ,midP,  1,laterP);

                            }
                            }
                         }
                        if ( f1.name[left].space[space].act < 0 ){
#if VERBOSE
                            printf("invert\n");
#endif
                            InvertOp(bd,-f1.name[left].space[space].act, N1, laterP, midP);
                            cblas_daxpy(N2, flow, midP , 1, outP, 1);
                        }else {
                            cblas_daxpy(N2, flow, laterP, 1, outP, 1);
                       }
                    }
                    su = f1.name[su].loopNext;//sum channel
                }
                
        }
//        else{
//
//        if ( bodies(f1,left) == bodies(f1,right))
//        {
//            inta N1 = vectorLen(f1, space);
//
//                cblas_dgemv( CblasColMajor, CblasNoTrans,  N1, N1,1.,
//                        streams( f1, left, lspin,space )+l*N1*N1, N1,
//                        streams(f1, inT, inS,space)+inR*N1,1, 0.,
//                        streams( f1, outT, outS,space2 )+outR*N1, 1  );
//
//        }else if ( bodies(f1,left) < bodies(f1,right))
//        {
//            inta N1 = outerVectorLen(f1,bodies(f1,left),space);
//            inta N2 = vectorLen(f1, space);
//                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,N1,N2/N1,N1,1.,streams( f1, left, lspin,space )+l*N1*N1,N1,streams(f1, inT, inS,space)+inR*N2,N1, 0.,     streams( f1, outT, outS,space2)+outR*N2,N1);
//        }
//        }
//        if (out != 1 ){
//            tPermuteOne(rank, f1, space2, out, outT,outR,outS, equals, e,espin);
//        }
    }
    return 0;
}

//inta tGEVV (inta rank,    sinc_label  f1,  inta space,  division equals,inta e,  inta espin,   division left,inta l,inta lspin,    division right,inta r, inta rspin ){
//    if ( header(f1, left ) != header(f1, right ) ){
//        printf("Two Head types GEVVo\n");
//        exit(0);
//    }
//    inta nl,nr,nm,space2 = space;
//      genusType sl = species(f1,left);
//      genusType sr = species(f1,right);
//    inta cr = CanonicalRank(f1, left, lspin);
//
//    if (( sl == outerVector ) && ( sr == vector ) ){
//          division inT,midT,outT;
//        inta inR,midR,outR,inS,midS,outS;
//        inT = right;
//        inR = r;
//        inS = rspin;
//        
//        outT = equals;
//        outR = e;
//        outS = espin;
//        if ( in != 1 ){
//           tPermuteOne(rank, f1, space, in, right, r, rspin, canonicalvvVector,0, rank);
//           inT = canonicalvvVector;
//           inR = 0;
//           inS = rank;
//        }
//        midT = canonicalvv2Vector;
//        midR = 0;
//        midS = rank;
//
//        
//        if (out != 1 ){
//            outT = canonicalvv3Vector;
//            outR = 0;
//            outS = rank;
//        }
//        if (sl == outerVector)
//            nl = outerVectorLen(f1, bodies(f1, left), space);
//        else
//            nl = vectorLen(f1, space);
//        
//        if (sr == outerVector)
//            nr = outerVectorLen(f1, bodies(f1, right), space);
//        else
//            nr = vectorLen(f1, space);
//
//        if ( nl == nr )
//        {
//            cblas_dcopy(nl,streams(f1,inT,inS,space)+inR*nl ,1,streams(f1,outT,outS,space)+outR*nl ,1  );
//            cblas_dscal(nl,tDOT(rank, f1, space, 1, midT, midR,midS,1, outT, outR, outS),streams(f1,outT,outS,space)+outR*nl,1);
//        }
//            else if (nl < nr ) {
//                inta l1 = l%cr, l2 = (l/cr)% cr;
//                printf("%d**%d\n",l1,l2);
//                        //dgemv on right * left -> reduced dimensional intermediate
//                        nm = nr/nl;
//                cblas_dgemv( CblasColMajor, CblasNoTrans,  nm,  nl,1.,
//                        streams( f1, inT, inS,space )+inR*nr, nm,
//                        streams(f1, left, lspin,space)+l1*nl,1, 0.,
//                        streams( f1, midT, midS,space )+midR*nm, 1  );
//                cblas_dger(CblasColMajor, nl,nm, 1. , streams(f1,left,lspin,space) + l2 * nl,1 , streams(f1, midT,midS,space)+midR*nm,1, streams(f1, outT,outS,space)+outR*nr,nl);
//            }else {
//                printf("dimensional error\n");
//                exit(0);
//            }
//                     
//        if (out != 1 ){
//            tPermuteOne(rank, f1, space2, out, outT,outR,outS, equals, e,espin);
//        }
//
//    }
//    return 0;
//}

/**
 *general vector - vector inner product per dimension
 *
 *@param rank      OMP rank
 *@param f1          container
 *@param leftChar  left input group action  defalult 1
 *@param left  vector
 *@param l  one of the matrices canonical ranks indexed
 *@param lspin matrix spin
 *@param rightChar right input group action
 *@param right input division
 *@param r one of the input canonical ranks indexed
 *@param rspin input spin
*/
double tDOT (inta rank,    sinc_label  f1,inta dim,char leftChar,   division left,inta l,inta lspin, char rightChar,   division right,inta r, inta rspin ){
    inta space = dim,brab,bspin ,ketk, kspin;
    double prod = 0.;
      division bra,ket;
    f1.name[canonicaldotVector].Current[rank] = 0;
    f1.name[canonicaldot2Vector].Current[rank] = 0;
    if ( rank ){
        ///check for parallel allocations
        if ( ! allowQ(f1.rt, blockParallelPermuteblock)){
            printf("blockParallelPermuteblock allow\n");
            fflush(stdout);
            exit(0);
        }
    }
    
    
    if ( rightChar != CDT){
        tPermuteOne(rank, f1, space, rightChar, right, r, rspin, canonicaldotVector,0, rank);
        bra = canonicaldotVector;
        brab = 0;
        bspin = rank;
    }else {
        bra = right;
        brab = r;
        bspin = rspin;
    }
    if (leftChar != CDT ){
        tPermuteOne(rank, f1, space, leftChar, left, l, lspin, canonicaldot2Vector,0, rank);
        ket = canonicaldot2Vector;
        ketk = 0;
        kspin = rank;
        
    }else{
        ket = left;
        ketk = l;
        kspin = lspin;
    }
    inta al = alloc(f1, right, space );
    if ( alloc(f1, left, space ) == alloc(f1, right, space )){
        inta N1 = al;
        prod = cblas_ddot( N1 , streams(f1,bra,bspin,space)+brab*N1,1 , streams(f1,ket,kspin,space)+ketk*N1, 1);
    } else {
        printf("body count %d %d %d %d\n",left,alloc(f1, left, space ) , right,alloc(f1, right, space ) );
        exit(1);
    }
    return prod;
}


/**
 *Multiply by dimension
 *
 *@param rank      OMP rank
 *@param f1          container
 *@param left matrix
 *@param l matrix index, across all chained elements
 *@param im matrix spin
 *@param prod scalar multiply
 *@param[in] ket vector
 *@param k ket index
 *@param sp2 ket spin
 *@param[out] oket vector
 *@param o oket index
 *@param ospin spin
 */
inta tHX(  inta rank,   sinc_label f1 ,division left, inta l, inta im, double prod,   division ket , inta k, inta sp2,   division oket, inta o,inta ospin ){
    if ( left == nullName || species(f1, left ) == scalar){
        printf("*");
        return 0;
    }
    
    division in=nullName,out=nullName;
    inta inSp,outSp,inRank,outRank;
    
    inta N2,flag,lll,found;
    inta dim,iter=0;
    if ( species(f1, left) == matrix && species(f1,ket) == vector && species(f1,oket) == vector){
        
        if ( rank ){
            ///check for parallel allocations
            if ( ! allowQ(f1.rt, blockParallelMultiplyblock)){
                printf("blockParallelMultiplyblock allow\n");
                fflush(stdout);
                exit(0);
            }
        }

        
        
          division ll = f1.name[left].chainNext,xx,zz;
                    inta mi = 0,xi=0;
                    while ( ll != nullName){
                        //chain will tie various separate terms into a single entity
                        zz = f1.name[left].chainNext;
                        found = 0;
                        while ( zz != ll ){
                            if (f1.name[zz].multId == f1.name[ll].multId)
                                found = 1;
                            zz = f1.name[zz].chainNext;
                        }
                        if ( ! found ){

                        
                        xi += CanonicalRank(f1, ll, im);
                        lll =  l-mi;
                        if ( mi <= l && l < xi )
                            {
                                
                                xx = ll;
                                while ( xx != nullName ){
                                    if ( f1.name[ll].multId == f1.name[xx].multId){
                                        iter++;
                                        if ( iter == 1 ){
                                            out =  oket;
                                            outRank = o;
                                            outSp = ospin;
                                            
                                            in = ket;
                                            inRank = k;
                                            inSp= sp2;
                                            
                                        }else if ( (iter%2) == 0 ){
                                            //SWAPPED!!@!
                                            out = multiplyVector;
                                            outRank = 0;
                                            outSp = rank;
                                            
                                            in = oket;
                                            inRank = o;
                                            inSp = ospin;
                                            
                                        }else {
                                            
                                            out = oket;
                                            outRank = o;
                                            outSp = ospin;
                                            
                                            in = multiplyVector;
                                            inRank = 0;
                                            inSp = rank;

                                        }
                                        
                                        
                                        flag = 1;
                                        for ( dim = 0 ; dim < SPACE ; dim++)
                                            if ( f1.canon[dim].body != nada){
                                                N2 = alloc(f1, out, dim);
                                                if ( f1.name[xx].space[dim].block != id0 )
                                                {
                                                    tGEMV(rank, f1, dim,out,outRank,outSp, xx, lll, im,in, inRank,inSp);
                                                }else{
                                                    topezOp(0,0, Bodies(f1,in,dim),f1.name[xx].space[dim].act,tv1, tv1,outerVectorLen(f1, one,dim),streams(f1,in,inSp, dim)+inRank*N2,0, streams(f1,out,outSp, dim)+outRank*N2);
                                                }
                                                if ( flag ){
                                                    cblas_dscal(N2, prod, streams(f1,out,outSp, dim)+outRank*N2, 1);
                                                    flag = 0;
                                                }

                                            }
                                    
                                        
                                        
//                                    else if( species(f1, left) == outerVector){
//
//                                        flag = 1;
//                                        for ( dim = 0 ; dim < SPACE ; dim++)
//                                            if ( f1.canon[dim].body != nada){
//                                                N2 = alloc(f1, ket, dim);
//                                                {
//                                                    tGEVV(rank, f1, dim,oket, o,ospin, xx, lll, im, ket, k, sp2);
//                                                    if ( flag ){
//                                                        cblas_dscal(N2, prod, streams(f1,oket,ospin, dim)+o*N2, 1);
//                                                        flag = 0;
//                                                        }
//                                                    }
//                                                }
//                                            }
                                        }
                                    xx = f1.name[xx].chainNext;
                                }
                            }
                            mi += CanonicalRank(f1, ll, im);

                        }
                        ll = f1.name[ll].chainNext;

                }
        
            if ( iter > 1 && (iter%2) == 0 ){
                for ( dim = 0 ; dim < SPACE ; dim++)
                    if ( f1.canon[dim].body != nada){
                        N2 = alloc(f1, ket, dim);
                        cblas_dcopy(N2, streams(f1,multiplyVector, rank, dim), 1, streams(f1,oket,ospin,dim)+o*N2, 1);
                    }
            }

    }else if ( species(f1, left) == vector && species(f1,ket) == vector && species(f1,oket) == vector){
        ///specifically for geminals * transpose-geminals = geminal
        if ( bodies(f1, left) == two && bodies(f1,in) == two && bodies(f1,out)== two ){
            inta space;
            for ( space = 0; space < SPACE ; space++)
                if ( f1.canon[space].body != nada ){
                    inta N1 = vector1Len(f1, space);
                    inta N2 = N1*N1;
                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,N1,N1,N1,1.,streams( f1, left, im,space )+l*N2,N1,streams(f1, ket,sp2,space)+k*N2,N1, 0.,streams( f1, oket, ospin,space)+o*N2,N1);
            }
        }
    }
    
    return 1;
}

/**
 *Multiply in total
 *
 *@param f1          container
 *@param[out] bra vector
 *@param left matrix
 *@param shiftFlag  for product * Hv +  sum v
 *@param[in] right vector
 *@param tolerance a number setting the absolute quality
 *@param relativeTolerance a number seting quality relative to magnitude of origin
 *@param condition Beylkin's condition (alpha)
 *@param threshold the smallest number
 *@param maxCycle the maxmium number of cycles in this routine
 *@param canon rank of output vector
*/
void tHXpY ( sinc_label f1 , division bra, division left,inta shiftFlag, division right , double tolerance , double relativeTolerance, double condition, double threshold, inta maxCycle, double maxCondition, inta canon, inta X1){
    double prod;
    inta rank0 = 0 ,rank;
    mea co2,coi;
    inta ilr,Ll,sp2,Rr,im,l , k,targSpin;
    division pt,Mat;
    
    if ( ! allowQ(f1.rt,blockTotalVectorBlock)){
        printf("blockTotalVectorBlock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    if (  right == totalVector){
        printf("you cannot feed totalVector into tHXpY\n");
        exit(0);
    }
    for ( targSpin = 0 ; targSpin < spins(f1, right ) ;targSpin++){
        pt = left;
        zero(f1, totalVector, rank0);
        f1.name[totalVector].Current[rank0] =0;

        if (  bra != totalVector){
            if ( shiftFlag  == 1){
                tEqua(f1, totalVector, rank0, bra, targSpin);
            }
            else if ( shiftFlag  == -1){
                tEqua(f1, totalVector, rank0, right, targSpin);
                mea ME;
                ME = tMatrixElements(rank0, f1, right, targSpin, left, 0, right, targSpin);
                tScaleOne(f1, totalVector, rank0, -ME);
            }
            else
                f1.name[totalVector].Current[rank0] = 0;
        }
         {
            Mat = pt;
            for ( im = 0; im < spins(f1, Mat ); im++)
                for ( sp2 = 0; sp2 < spins(f1,right); sp2++)
                {
                    if (sp2 == 1)
                        co2 = I;
                    else
                        co2 = 1;
                    if (im == 1 )
                        coi = I;
                    else
                        coi = 1;
                    if ( targSpin == 0 )
                        prod  = creal(co2 * coi );
                    else
                        prod = cimag(co2 * coi ) ;
                    if ( fabs(prod) > f1.rt->THRESHOLD ){
                        Rr = CanonicalRank ( f1, right,sp2 );
                        Ll = CanonicalOperator(f1, Mat, im);
                        inta su = f1.name[totalVector].Current[rank0];
                        
                        #ifdef OMP
                        #pragma omp parallel for private (ilr,l,k,rank) schedule(dynamic,1)
                        #endif
                            for ( ilr = 0; ilr <  Ll*Rr ; ilr++)
                            {

                            #ifdef OMP
                                    rank = omp_get_thread_num();
                            #else
                                    rank = 0;
                            #endif

                                l = ilr%Ll;
                                k = ilr/Ll;
                                tHX(rank, f1,Mat, l, im,prod, right, k, sp2,totalVector,ilr +su, rank0);
                            }
#if VERBOSE
                        printf("'lambda' -> %d %d %d\n", su , su+Ll*Rr,part(f1, totalVector ));
                        
#endif
                        f1.name[totalVector].Current[rank0]+= Ll*Rr;
                        if (f1.name[totalVector].Current[rank0] > part(f1, totalVector ) )
                        {
                            printf("'lambda' is too small\n");
                            fflush(stdout);
                            exit(1);
                        }
                    }
                }
        };
        
        if (  bra != totalVector){
            CanonicalRankDecomposition( f1,  NULL,totalVector, rank0, bra, targSpin, tolerance,relativeTolerance, condition,threshold,maxCycle,maxCondition, canon);
        }
    }
    return;
}


/**
 *matrix element
 *upgraded in v9.3 for parallel operation
 *only one of these can run at a time.
 *
 *@param f1 container
 *@param alloy1 vector
 *@param spin1 vector's spin
 *@param op operator
 *@param ospin operator spin
 *@param alloy2 vector
 *@param spin2 vector's spin

 */
double pMatrixElement ( sinc_label  f1 ,   division alloy1 , inta spin1, division op, inta ospin, division alloy2 , inta spin2 ){
    inta rank,i ,cl = CanonicalRank(f1, alloy1, spin1),cr = CanonicalRank(f1, alloy2, spin2);
    mea OV, ov[f1.rt->NLanes];
    
    division left=0;
    for ( rank = 0; rank < f1.rt->NLanes; rank++){
        ov[rank] = 0.;
        if ( ! rank )
            left = anotherLabel(&f1, 0, nada);
        else
            anotherLabel(&f1, 0, nada);

        f1.name[left+rank].Current[spin1] = 1;
        f1.name[left+rank].name = alloy1;
    }
    ///need to be in order to make each rank-array contiguous.
    division right=0 ;
    for ( rank = 0; rank < f1.rt->NLanes; rank++){
        if ( ! rank )
            right = anotherLabel(&f1, 0, nada);
        else
            anotherLabel(&f1, 0, nada);
        f1.name[right+rank].Current[spin2] = 1;
        f1.name[right+rank].name = alloy2;
    }

    
    #ifdef OMP
    #pragma omp parallel for private (i,rank) schedule(dynamic,1)
    #endif
                for ( i = 0 ;i < cl*cr; i++)
                {
    #ifdef OMP
                    rank = omp_get_thread_num();
    #else
                    rank = 0;
    #endif
                    f1.name[left+rank].Begin[spin1] = (i)%cl+f1.name[alloy1].Begin[spin1];
                    f1.name[right+rank].Begin[spin2] = (i/cl)+f1.name[alloy2].Begin[spin2];
                    ov[rank] += (tMatrixElements(rank, f1, left+rank, spin1, op, ospin, right+rank, spin2));
                }
    OV = 0.;
    for ( rank = 0; rank < f1.rt->NLanes; rank++)
        OV += ov[rank];

#ifdef COMPLEXME
    return (creal(OV));
#else
    return (OV);
#endif
}

/**
 *outer product of two vectors
 *
 *for building operators from vectors
 *
 *@param f1 container
 *@param space   in canon which has SPACE components
 *@param[in] vector any content
 *@param a vector1 spin
 *@param[in] vector2 any content
 *@param b vector's spin
 *@param[out] proj outer product matrix
 *@param c proj spin
*/
inta tOuterProductSuOne(   sinc_label  f1,inta space,  division vector , inta a,   division vector2,inta b,   division proj, inta c){
    inta ma = CanonicalRank(f1, vector,a), mb = CanonicalRank(f1, vector,b), zc = CanonicalRank(f1, proj, c),l,r,i;
    inta Np = alloc(f1, proj, space), n1 = alloc(f1, vector,space),n2 = alloc(f1, vector2,space);

    if ( species(f1, vector ) == matrix || species(f1, vector2) == matrix){
        printf("input arguments to Outer product need to be vectors\n");
        exit(0);
    }
    if ( zc + ma*mb > part(f1, proj) && proj < eigenVectors){
        printf("outerProductSu:  \n");
        exit(0);
    }
    if ( Np != n1 * n2 ){
        printf("suspect...");
        exit(1);
    }
    
    for ( l = 0 ; l < ma; l++)
        for ( r = 0; r < mb; r++)
                for ( i = 0; i < Np ; i++)
                    (streams(f1, proj,c,space)+(l*mb+r+zc)*Np)[i] = 0.;
    
    for ( l = 0 ; l < ma; l++)
        for ( r = 0; r < mb; r++)
                cblas_dger(CblasColMajor, n1,n2, 1. , streams(f1, vector,a,space)+l*n1,1, streams(f1, vector2,b,space)+r*n2,1, streams(f1, proj,c,space)+(l*mb+r+zc)*Np,n1);
    f1.name[proj].Current[c] += ma*mb;
    return 0;
}

//inta compressReplaceEikon(  sinc_label f1 ,   division eik ){
//    if ( species(f1, eik) != eikon){
//        return 1;
//    }
//    inta i = 0,ii;
//
//      division looper,chainer;
//      genus hidden;
//
//    for ( hidden = eikonDiagonal ; hidden <= eikonSemiDiagonal ; hidden++){
//        i = 0;
//        tClear(f1,copyVector);
//        tClear(f1,copyTwoVector);
//
//        for ( chainer = eik; chainer != nullName ; chainer= f1.name[chainer].chainNext){
//
//        for ( looper = chainer; looper != nullName ; looper= f1.name[looper].loopNext)
//            if ( f1.name[looper].species == hidden ){
//
//            tAddTw(f1, copyVector,  0, looper, 0);
//            i++;
//        }
//    }
//
//        for ( ii = 1 ; ii < 30 ; ii++){
//            tId(f1,copyTwoVector,0);
//            CanonicalRankDecomposition(0, f1, NULL, copyVector, 0, copyTwoVector, 0,  f1.rt->TOLERANCE,  f1.rt->relativeTOLERANCE,  f1.rt->ALPHA, f1.rt->THRESHOLD, f1.rt->MAX_CYCLE, <#inta canon#>)
//        }
//    }
//
//    return 0;
//}

