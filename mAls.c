/*
 *  mAls.c
 *
 *
 *  Copyright 2019 Jonathan Jerke and Bill Poirier.
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

#include "coreUtil.h"
#include "mAls.h"

INT_TYPE normalize (struct field * f1,  enum division alloy, INT_TYPE spin, INT_TYPE space){
    INT_TYPE l,flag=0;
    double norm;
    INT_TYPE M2[SPACE];
    length(f1,alloy,M2);
    struct name_label u = f1->sinc.tulip[alloy];
    
    INT_TYPE iOne = 1;
    {
        for ( l = 0; l < CanonicalRank(f1, alloy,spin) ;l++){
            norm = cblas_dnrm2(M2[space], streams(f1, alloy,spin,space)+l*M2[space],iOne);
            if ( fabs(norm) > 0.  ){
                norm = 1./norm ;
                cblas_dscal(M2[space], norm, streams(f1, alloy,spin,space)+l*M2[space],iOne);
            }else {
#if VERBOSE
                printf("ACK! %d %lld %lld\n", alloy,spin,l);
#endif
                
                if (1 ){
                    tReplace(f1, alloy, spin , space,l );
                    norm = cblas_dnrm2(M2[space], streams(f1, alloy,spin,space)+l*M2[space],iOne);
                    if ( fabs(norm) > 0.  ){
                        norm = 1./norm ;
                        cblas_dscal(M2[space], norm, streams(f1, alloy,spin,space)+l*M2[space],iOne);
                    }else {
                        printf("WTF\n");
                        exit(0);
                    }
                } else {
                    flag = 1;
                    if (f1->sinc.tulip[alloy].Current[spin] > 2 )
                        f1->sinc.tulip[alloy].Current[spin]--;
                    else{
                        tReplace(f1, alloy, spin, space, l);
                    }
                }
            }
        }
    }
    return flag;
    
}

INT_TYPE spread (struct field * f1, enum division origin, INT_TYPE os, enum division alloy, INT_TYPE spin, INT_TYPE space, Stream_Type * output,Stream_Type * output2){
    INT_TYPE m,n;
    if (normalize(f1,alloy,spin,space) )
        return 1;
    INT_TYPE L1 = CanonicalRank(f1, alloy,spin);
    INT_TYPE G1 = CanonicalRank(f1, origin,os);
    INT_TYPE M2[SPACE];
    length(f1,alloy,M2);
    double *Alloy =streams(f1,alloy,spin,space);
    double *Origin =streams(f1,origin,os,space);
//    printf("spread\n");
//    fflush(stdout);
    {
        
        for ( m = 0; m < L1; m++)
            for ( n = 0; n <=m ; n++){
                output[ m*L1 + n ] = cblas_ddot(M2[space], Alloy+m*M2[space] ,1,Alloy+n*M2[space] ,1);
                output[ n*L1 + m ] = output[ m*L1 + n ];
            }
        if ( output2 != NULL )
        for ( m = 0; m < L1; m++)
            for ( n = 0; n < G1 ; n++)
                output2[ n*L1 + m ] = cblas_ddot(M2[space], Alloy+m*M2[space] ,1,Origin+n*M2[space] ,1);
    }
    
    //printf("%f__\n",output[0]);
    return 0;
}

INT_TYPE pSpread (struct field * f1, enum division origin, INT_TYPE os, enum division alloy, INT_TYPE spin, INT_TYPE space, Stream_Type * output,Stream_Type * output2){
    INT_TYPE m,n;
    if (normalize(f1,alloy,spin,space) )
        return 1;
    INT_TYPE L1 = CanonicalRank(f1, alloy,spin);
    INT_TYPE nm,G1 = CanonicalRank(f1, origin,os);
    INT_TYPE M2[SPACE];
    length(f1,alloy,M2);
    double *Alloy =streams(f1,alloy,spin,space);
    double *Origin =streams(f1,origin,os,space);
//    printf("Pspread\n");
//    fflush(stdout);

    {
#pragma omp parallel for private (m,n,nm)
        for ( nm = 0; nm < L1*L1; nm++){
            m = (nm)%L1;
            n = (nm/L1)%L1;
            if ( n <= m ){
                output[ nm ] = cblas_ddot(M2[space], Alloy+m*M2[space] ,1,Alloy+n*M2[space] ,1);
                output[ m*L1 + n ] = output[ nm ];
            }
        }
        if ( output2 != NULL )
#pragma omp parallel for private (m,n,nm)
            for ( nm = 0; nm < L1*G1; nm++){
                m = (nm)%L1;
                n = (nm/L1)%G1;
                    output2[ nm ] = cblas_ddot(M2[space], Alloy+m*M2[space] ,1,Origin+n*M2[space] ,1);
            }
    }
    return 0;
}


INT_TYPE balance (struct field * f1,  enum division alloy, INT_TYPE spin){
    INT_TYPE l;
    INT_TYPE L1 = CanonicalRank(f1, alloy,spin),sign[SPACE],signs;
    double sum,snorm,value;
    long double prod,factor,norm[SPACE] ;
    INT_TYPE n1,ii,M2[SPACE],space,spaces=0;
    length(f1,alloy,M2);
    INT_TYPE iOne = 1;
    //#pragma omp parallel for private(l,norm,flag,ii,i,trace,space)
    if ( species(f1, alloy ) == vector ){
      //IGNORE
        return 0;
    }
    {
        for ( l = 0; l < L1 ;l++){
            
            for ( space = 0; space < SPACE ; space++)
                if ( f1->sinc.rose[space].body != nada){
                    spaces++;
                    norm[space] = cblas_dnrm2(M2[space], streams(f1, alloy,spin,space)+l*M2[space],iOne);
                    sum = 0.;
                    n1 = outerVectorLen(f1, species(f1,alloy), space);
                    for ( ii = 0; ii < n1; ii++)
                        sum += ( streams(f1, alloy,spin,space)+l*M2[space])[ii*n1+ii];
                    if ( sum >= 0 )
                        sign[space] = 1;
                    else
                        sign[space] = -1;
                }
            
            
            prod = 1;
           // printf("%d ::", alloy);
            for ( space = 0 ;space < SPACE ; space++)
                if ( f1->sinc.rose[space].body != nada){
                    prod *= sign[space]*norm[space];
              //      printf("%f ", sign[space]*norm[space]);
                    
                }
          //  printf("\n");
            if ( prod >= 0 )
                signs = 1;
            else
                signs = -1;
            factor = powl( fabs(prod),1./spaces);
            for ( space = 0; space < SPACE ; space++)
               if ( f1->sinc.rose[space].body != nada)
               {
                snorm = factor/norm[space]*signs/sign[space] ;
                cblas_dscal(M2[space], snorm, streams(f1, alloy,spin,space)+l*M2[space],iOne);
               }
            
        }
    }
    return 0;
}

double canonicalListDecompositionMP( INT_TYPE rank,struct field * f1 , Stream_Type * cofact, enum division origin,INT_TYPE os,   enum division alloy ,  INT_TYPE spin ,double tolerance,double magn, INT_TYPE preferred){
    if ( tolerance < 0 ){
        tolerance = -(tolerance);
    }
    INT_TYPE yes = 0,singleFlag = 0;
    if ( rank < 0){
        singleFlag = 1;
        rank = 0;
    }
    double value,value2 =1,cross,other2,value3;
    //printf("%d (%d)\n", alloy,CanonicalRank(f1, alloy, 0));
    //fflush(stdout);
    INT_TYPE g,l,count = 1 ,spaces=0;
    
    INT_TYPE ns = iNS;
    double ALPHA = f1->mem1->rt->ALPHA;
//    if ( cofact != NULL )
//    for ( l = 0 ; l < CanonicalRank(f1, origin, os); l++)
//        printf("%d %f\n", l, cofact[l]);
//    printf("\n");
    Stream_Type * array[SPACE];
    Stream_Type * array2[SPACE];
    Stream_Type * guide, *track;
    INT_TYPE space,space0;
    INT_TYPE L1 = CanonicalRank(f1,alloy,spin);;
    INT_TYPE G1 =CanonicalRank(f1,origin,os);
    if ( G1 == 0 ){
        f1->sinc.tulip[alloy].Current[spin] = 0;
        printf("CD: empty origin %d _%lld -> %d _%lld\n", origin, os , alloy, spin);
        return 0;
    }
    for ( space = 0; space < SPACE ; space++)
        if ( f1->sinc.rose[space].body != nada )
            spaces++;

    if ( L1 == 0 ){
        // printf("CD: zero length %lld %lld %lld\n", origin,alloy,spin);
        tId(f1, alloy, spin);
        L1 = 1;
    }
    
    INT_TYPE M2[SPACE];
    length(f1,alloy,M2);
    INT_TYPE info;
    double subValue2 = 1,subValue3;
    double rcond;
    enum division canonicalStore;
    INT_TYPE T1 = G1*L1;
    {//REF CORE NUMBER
        //GET CORE... ALLOCATE TO CORE-LANE.
        
        
        for ( space = 0; space < SPACE ; space++){
            array[space] =  streams(f1, canonicalBuffers, rank , space);
            array2[space] =  streams(f1, canonicalBuffers, rank , space) + L1*L1;
            
        }
        guide =  myStreams(f1, guideBuffer, rank );
        track =  myStreams(f1, trackBuffer, rank );

        
        if (  T1 + L1*L1 >  part(f1,canonicalBuffers)){
#if 1
            printf("mem prob\n %d %d \n", L1, G1);
#endif
            exit(0);
        }
        
        if ( species(f1, origin ) == matrix )
            canonicalStore = canonicalBuffersBM;
            else
        canonicalStore = canonicalBuffersB;
        
        if ( part(f1, canonicalStore) < L1  )
        {
            printf("list alloc\n");
            exit(0);
        }
        
        
    }
    INT_TYPE dim0 = 0;

    INT_TYPE dim[SPACE],ll;
   
    {

        space0 = 0;
        
        dim0=0;
        if ( preferred == -1 ){
            for ( space = 0; space < SPACE ; space++)
                if ( f1->sinc.rose[space].body != nada)
                    dim[(space0+space)%(spaces)] = dim0++;
            
        } else {
            for ( space = 0; space < SPACE ; space++)
                if ( f1->sinc.rose[space].body != nada)
                    
                    dim[dim0++] = (preferred+space)%(spaces);
            

        }
        for ( space = 0; space < dim0 ; space++)
            if ( f1->sinc.rose[dim[space]].body != nada){
                
                
                if ( !singleFlag){
                 if ( spread(f1,origin,os,alloy,spin,dim[space],array[dim[space]],array2[dim[space]]) )
                    return -1;
                }else
                    if ( pSpread(f1,origin,os,alloy,spin,dim[space],array[dim[space]],array2[dim[space]]) )
                        return -1;

        }
//        {
//            INT_TYPE spa;
//            for ( spa = 0 ; spa < SPACE ; spa++)
//                if ( f1->sinc.rose[spa].body != nada )
//                    printf("*%1.3f*", cblas_dnrm2(L1*M2[spa], streams(f1, alloy, spin, spa), 1)/sqrt(L1));
//            printf("\n");
//        }
        
        space0 = 0;
        while (1 ){
            //all computed inner products
//            printf("dim[0]=%d\n", dim[0]);
            dim0 = 0;
            if ( preferred == -1 ){
                for ( space = 0; space < SPACE ; space++)
                    if ( f1->sinc.rose[space].body != nada)

                        dim[(space0+space)%(spaces)] = dim0++;
                
            }
            cblas_dcopy(L1*L1, array[dim[1]],1,track,1);
            for ( space = 2; space < dim0 ; space++)
                if ( f1->sinc.rose[dim[space]].body != nada)
                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1*L1, 0,array[dim[space]],1, track,1 );
            if ( L1 >1 )
                for ( l = 0 ; l < L1; l++)
                    track[ l*L1 + l ] += ALPHA;

            if (space0 == 0 && L1 > 1){
                cblas_dcopy(L1*L1, track, 1, track+L1*L1, 1);
                tdsyev(rank, f1, 'V', L1, track+L1*L1, L1, track+2*L1*L1);
                
                rcond = (track+2*L1*L1)[L1-1]/(track+2*L1*L1)[0];
                if (  isnan(rcond) || isinf(rcond) || rcond > 1e7){
#if VERBOSE
                    printf("Warning: condition of Beylkin\n");
                    fflush(stdout);
#endif
                     return -1;
                }
                //maybe one at a time???

            }
            cblas_dcopy(T1, array2[dim[1]],1,guide,1);
            for ( space = 2; space < dim0 ; space++)
                if ( f1->sinc.rose[dim[space]].body != nada)
                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,T1, 0,array2[dim[space]],1, guide,1 );
            if ( cofact != NULL )
                for ( g = 0; g < G1 ; g++)
                    for ( l = 0; l < L1 ; l++)
                        guide[g*L1+l] *= cofact[g];

            // Vectors  L1 x G1
            // list...  L1 x M2 ==   ( cross * gstream**T )
            
            
            
            double *pt = myStreams(f1,canonicalStore , rank);
            
            if ( 0 ) {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,L1,M2[dim[0]],G1,1,guide,L1,streams(f1,origin,os,dim[0]),M2[dim[0]], 0, pt, L1 );
                
               // guide //  L1 * G1
                // origin // M2 * G1
                info = tdpotrf(L1, track);
                if ( info != 0 ){
#if VERBOSE
                    printf("info failure \n");
#endif
                    return -1;
                }

                info = tdpotrs(L1,  M2[dim[0]], track,  pt,L1 );
                if ( info != 0 ){
#if VERBOSE
                    printf("info failture\n");
#endif
                    return -1;
                }
                transpose ( L1, M2[dim[0]] , pt , streams(f1,alloy,spin,dim[0]));
            } else {
                info = tdpotrf(L1, track);
                if ( info != 0 ){
#if VERBOSE
                    printf("info failure \n");
#endif
                    return -1;
                }

                for ( ll = 0; ll < M2[dim[0]] ; ll++){
                    cblas_dgemv(CblasColMajor, CblasNoTrans,L1,G1,1.,guide,L1,streams(f1,origin,os,dim[0])+ll,M2[dim[0]], 0., pt,1);
                    if ( 0 ){
                        cblas_dcopy(L1*L1, track, 1, track+L1*L1, 1);
                        info = tdgels(rank, f1, L1, 1, track+L1*L1,  pt, L1);
                        if ( info != 0 ){
#if VERBOSE
                            printf("info failure \n");
#endif
                            return -1;
                        }

                    } else {
                        
                        info = tdpotrs(L1,  1, track,  pt,L1 );
                        if ( info != 0 ){
#if VERBOSE
                            printf("info failture\n");
#endif
                            return -1;
                        }

                    }
                    cblas_dcopy(L1, pt, 1, streams(f1,alloy,spin,dim[0])+ll, M2[dim[0]]);

                }

            }
            
            
            subValue3 = subValue2;
            
//            {
//                INT_TYPE spa;
//                for ( spa = 0 ; spa < SPACE ; spa++)
//                    if ( f1->sinc.rose[spa].body != nada )
//                        printf("%1.3f ", cblas_dnrm2(L1*M2[spa], streams(f1, alloy, spin, spa), 1)/sqrt(L1));
//                printf("\n");
//            }
            
            subValue2 =  cblas_dnrm2(L1*M2[dim[0]], streams(f1, alloy, spin, dim[0]), 1)/magn/sqrt(L1);
            
            if ( fabs(subValue2 - subValue3) < tolerance ){
                yes++;
            }else {
                yes = max(yes-2,0);
            }

            if ( yes > 100 ){
              //  printf("ct %d\n", count);
                return 0;
            }
            
            ///printf("(%3.3f %3.3f)\n", subValue2, subValue3);
            count++;
            if ( space0 > 10000 )
            {
               // printf("ct(fail) %d\n", count);
                return 0;
            }
            if ( !singleFlag){
                if ( spread(f1,origin,os,alloy,spin,dim[0],array[dim[0]],array2[dim[0]]) )
                    return -1;
            }else
                if ( pSpread(f1,origin,os,alloy,spin,dim[0],array[dim[0]],array2[dim[0]]) )
                    return -1;

            space0++;

        }
    }
    return -1;
}



double tCycleDecompostionListOneMP ( INT_TYPE rank, struct field * f1 , enum division origin,INT_TYPE os, double * coeff, enum division alloy,INT_TYPE spin,  double tolerance , INT_TYPE maxRun , double power  ){
    INT_TYPE ct,ft;
    double toleranceAdjust = 1.,past = 1e9,prev=1e9,value,value2,cross,other2 ;//= cblas_dnrm2 ( part(f1, origin), coeff, 1);
    if ( coeff == NULL )
        value2 = (inner(imax(0,rank),f1, origin,os));
    else
        value2 = tInnerListMP(imax(0,rank), f1, origin, coeff);
    if (! CanonicalRank(f1, origin, os ) )
        return 0;
//    if ( value2 < f1->mem1->rt->TARGET ){
//        f1->sinc.tulip[alloy].Current[spin]=0;
//        return 0;
//    }
    do{
        ct = canonicalListDecompositionMP(rank, f1, coeff, origin, os,alloy, spin, toleranceAdjust*tolerance,value2,-1);
    
        if ( coeff == NULL ){
            value = distanceOne(imax(0,rank), f1, origin, os, alloy, spin);
        }else {
            cross = tInnerVectorListMP(imax(0,rank), f1, origin, coeff, alloy, spin) ;
            other2 = inner(imax(0,rank), f1, alloy, spin);
            value = fabs(value2 - 2. * cross + other2 );
        }
        if ( value >= prev ){
            ft = 1;
            printf("%d (%d) %f -- /%f/\n",alloy,CanonicalRank(f1, alloy, spin),value,value2);
            fflush(stdout);
            if ( fabs(value ) < f1->mem1->rt->TARGET*fabs(value2)  ){
                return 0;
            }
        }else {
            ft = 0;
        }
        if (ct == -1 ){
#if 1
            printf("List bailed \n");
#endif
            f1->sinc.tulip[alloy].Current[spin]--;
            canonicalListDecompositionMP(rank, f1, coeff, origin, os,alloy, spin, toleranceAdjust*tolerance,value2,-1);
            return 0;
        }
        else if ( ct == 1 || ! ft ){
            toleranceAdjust /= 1.1;

        }
        else if ( ct == 0 ){
            if ( CanonicalRank(f1, alloy, spin ) +1 <=  maxRun ){
               
                if ( past < value )
                    return 0;
                
                tId(f1, alloy, spin);
                
                
                past = value;
            }else {
                //may wnat to cut back teh failure here...
                return 0;
            }
        }
        
        prev = value;
        
    }while(1);
    return 0.;
}

double tInnerVectorListMP( INT_TYPE rank, struct field * f1 , enum division origin, double * coeff, enum division vector,INT_TYPE spin ){
    
    double sum = 0.,prod;
    INT_TYPE i = 0,space,vv;
    for ( i = 0; i < CanonicalRank(f1, origin,0 ) ; i++ )
        for ( vv = 0; vv < CanonicalRank(f1, vector,spin  ) ; vv++ ) {
            {
                prod = 1;
                for ( space = 0 ; space < SPACE ; space++)
                    if ( f1->sinc.rose[space].body != nada )
                    prod *= tMultiplyOne(rank, f1, space, nullName, 0, origin, i, 0, vector,vv, spin);
                
                sum += coeff[i]*prod;
            }
        }
    return sum;
}

double tInnerListMP( INT_TYPE rank, struct field * f1 , enum division origin, double * coeff ){
    
    double sum = 0.,prod;
    INT_TYPE j,i = 0,space;
    for ( i = 0; i < CanonicalRank(f1, origin,0 ) ; i++ )
        for ( j = 0 ;j < CanonicalRank(f1, origin, 0); j++){
            prod = 1;
            for ( space = 0 ; space < SPACE ; space++)
                if ( f1->sinc.rose[space].body != nada )

                prod *= tMultiplyOne(rank, f1, space, nullName, 0, origin, i, 0, origin, j, 0);
            
            sum += coeff[i]*coeff[j]*prod ;
    }
    
    return sum;
}

INT_TYPE printExpectationValues (struct field * f1 , enum division ha  , enum division vector){
    DCOMPLEX co,expat,totx,me,ov,exov;
    totx = 0.;
    enum division leftP = ha,Mat;
    if (CanonicalRank(f1, vector, 1))
        printf("%d (%d + i%d)\n", vector, CanonicalRank(f1, vector, 0), CanonicalRank(f1, vector, 1));
    else
        printf("%d (%d)\n", vector, CanonicalRank(f1, vector, 0));

    do {
        Mat = leftP;
        me = 0.;
        ov = 0.;
        //outputFormat(f1, stdout, proton1, 0);
        matrixElements(0, f1, vector,  Mat,  vector,&me,&ov);
        if ( Rank(f1, name(f1,Mat))){
            if (CanonicalRank(f1, name(f1,Mat), 1))
                printf("%d\t ::\t %d \t| (%d + i%d)\t\t%f %f\n",bodies(f1, Mat),name(f1,Mat), CanonicalRank(f1, name(f1,Mat), 0),CanonicalRank(f1, name(f1,Mat), 1),creal(me/ov), cimag(me/ov));
                else
                    printf("%d\t ::\t %d \t| (%d) \t\t%f\n",bodies(f1, Mat),name(f1,Mat), CanonicalRank(f1, name(f1,Mat), 0),creal(me/ov));
            }
        fflush(stdout);
        totx += me/ov;
        leftP = f1->sinc.tulip[leftP].linkNext;
        
    } while ( leftP != nullName);
    f1->sinc.tulip[vector].value.value = totx;

    if (fabs(cimag(totx)) > 0. )
        printf("\t\t\t%f\t%f\n",creal(totx),cimag(totx));
    else
        printf("\t\t\t%f\n",creal(totx));

    return 0;
}


void matrixElements ( INT_TYPE rank,struct field * f1 , enum division bra, enum division mat, enum division ket, DCOMPLEX *ME,DCOMPLEX *OV ){
    INT_TYPE r,l,e,dim,sp,im,sp2,im2,l2,spm,spx,sp2m,sp2x;
    double prod,bracket[SPACE];
    double va,va2,va3;
    DCOMPLEX co,co2,coi,coi2;
    if ( f1->sinc.tulip[bra].spinor == parallel || f1->sinc.tulip[ket].spinor == parallel){
        spm = rank;
        spx = rank+1;
        sp2m = rank;
        sp2x= rank+1;
    }
    else {
        spm = 0;
        spx = spins(f1, bra);
        sp2m = 0;
        sp2x = spins(f1,ket);
    }
    *OV = 0.;
    if ( species(f1, mat ) != outerVector ){
        co = 1;
        for ( sp = spm ; sp < spx;sp++){
            if (sp == 1 && f1->sinc.tulip[bra].spinor == cmpl)
                co *= -I;
            co2 = 1.;
            for ( e = 0 ; e < CanonicalRank(f1, bra, sp);e++)

            for ( sp2 = sp2m ; sp2 < sp2x;sp2++){
                if (sp2 == 1&& f1->sinc.tulip[ket].spinor == cmpl)
                    co2 *= I;
                for ( r = 0 ; r < CanonicalRank(f1, ket, sp2);r++){
                    prod = 1.;
                    for ( dim = 0 ; dim < SPACE ; dim++)
                        if ( f1->sinc.rose[dim].body != nada){
                            bracket[dim] = tDOT(rank, f1,dim,'T' , bra, e,sp,'N' ,  ket, r, sp2);
                            prod *= bracket[dim];
                        }
                    *OV += co*co2*prod;
                    if ( mat == nullName )
                        continue;
                    coi=1.;
                    for ( im = 0 ; im < spins(f1, mat);im++){
                        if (im == 1 && f1->sinc.tulip[mat].spinor == cmpl)
                            coi*= I;
                        
                        for ( l = 0 ; l < CanonicalRank(f1, name(f1,mat), im);l++){
                            
                            prod = 1;
                            for ( dim = 0 ; dim < SPACE ; dim++)
                                if ( f1->sinc.rose[dim].body != nada){
                                    
                                    ////bra-ket
                                    if ( f1->sinc.tulip[mat].space[dim].block == id0 ){
                                        if ( im == 0 )
                                            prod *= bracket[dim];
                                        else
                                            prod = 0;
                                    }
                                    else {
                                        f1->sinc.tulip[canonicalmeVector].Current[rank] =0;
                                        tGEMV(rank, f1, dim,canonicalmeVector, rank, mat, l, im, ket, r, sp2);
//                                        if ( mat == proton1 ){
//                                            f1->sinc.tulip[canonicalmeVector].Current[rank] =1;
//                                            outputFormat(f1, stdout, canonicalmeVector, rank);
//                                            f1->sinc.tulip[canonicalmeVector].Current[rank] =0;
//
//                                            outputFormat(f1, stdout, ket, sp2);
//
//                                        }
                                        prod *= tDOT(rank, f1,dim,'T', bra, e, sp,'N', canonicalmeVector, 0, rank);
//                                         va = tDOT(rank, f1,dim,'T', bra, e, sp,'N',ket,r,sp2);
//                                         va2 = tDOT(rank, f1,dim,'T', bra, e, sp,'N',bra,e,sp);
//                                        va3 = tDOT(rank, f1,dim,'T', ket,r, sp2,'N',ket,r,sp2);
                                    }
                                }
                            *ME += co2*co *coi* prod;
                        }
                    }
                }
            }
        }
    }
    else {
            co = 1;
            for ( sp = 0 ; sp < spins(f1, bra);sp++){
                if (sp == 1&& f1->sinc.tulip[bra].spinor == cmpl)
                    co *= -I;
                for ( e = 0 ; e < CanonicalRank(f1, bra, sp);e++){

                co2 = 1.;
                for ( sp2 = 0 ; sp2 < spins(f1, ket);sp2++){
                    if (sp2 == 1&& f1->sinc.tulip[ket].spinor == cmpl)
                        co2 *= I;
                        for ( r = 0 ; r < CanonicalRank(f1, ket, sp2);r++){
                            prod = 1.;
                            for ( dim = 0 ; dim < SPACE ; dim++)
                            if ( f1->sinc.rose[dim].body != nada){
                                bracket[dim] = tDOT(rank, f1,dim,'T' , bra, e,sp,'N' ,  ket, r, sp2);
                                prod *= bracket[dim];
                            }
                            *OV += co*co2*prod;
                            coi=1.;
                            for ( im = 0 ; im < spins(f1, mat);im++){
                                if (im == 1 )
                                    coi*= I;
                                for ( l = 0 ; l < CanonicalRank(f1, name(f1,mat), im);l++){
                                    coi2 = 1.;
                                    for ( im2 = 0 ; im2 < spins(f1, mat);im2++){
                                        if (im2 == 1 )
                                            coi2*= -I;
                                        for ( l2 = 0 ; l2 < CanonicalRank(f1, name(f1,mat), im2);l2++){
                                            prod = 1;
                                            for ( dim = 0 ; dim < SPACE ; dim++)
                                            if ( f1->sinc.rose[dim].body != nada){
                                                enum division c1 = complement + (f1->sinc.rose[dim].body - f1->sinc.tulip[mat].space[dim].body);
                                                enum division c2 = complementTwo + (f1->sinc.rose[dim].body - f1->sinc.tulip[mat].space[dim].body);

                                                if ( c1 < complement ){
                                                    printf("wrong dimensionality   %d ", c1 - complement);
                                                    exit(0);
                                                }
                                                
                                                if ( c1 != complement ){
                                                    f1->sinc.tulip[c1].Current[rank] =0;
                                                    f1->sinc.tulip[c2].Current[rank] =0;
                                                    
                                                    tGEVV(rank, f1, dim, c1, rank, 'T', mat, l, im, 'N', ket, r, sp2);
                                                    tGEVV(rank, f1, dim, c2, rank, 'T', mat, l2, im2, 'N', bra, e, sp);
                                                    
                                                    prod *= tDOT(rank, f1, dim, 'T', c2, 0, im2, 'N', c1, 0,im);
                                                }else {
                                                    prod *= tDOT(rank, f1, dim, 'T', mat, l, im, 'N', ket, r, sp2)*
                                                    tDOT(rank, f1, dim, 'T', mat, l2, im2, 'N', bra, e, sp);
                                                }
                                        }
                                    }
                                    *ME += -co2*co *coi* coi2*prod;
                                    }
                                }
                            }
                        }
                }
            }
        }
    }
    return ;
}



//double canonicalMultiplyMP( INT_TYPE rank,struct field * f1 , INT_TYPE begin, INT_TYPE end,char mc,enum division mat,INT_TYPE ms, enum division vec, INT_TYPE vs,   enum division alloy ,  INT_TYPE spin ,double tolerance){
//    if ( tolerance < 0 ){
//        tolerance = -(tolerance);
//    }
//    if ( vec == productVector || alloy == productVector ){
//        printf("canonicalMulMP\n");
//        exit(0);
//    }
//    INT_TYPE ns = iNS;
//    Stream_Type * array[SPACE];
//    Stream_Type * array2[SPACE];
//    INT_TYPE space,space0,space1,space2;
//    INT_TYPE l,count = 1 ;
//    INT_TYPE L1 = CanonicalRank(f1,alloy,spin);;
//    INT_TYPE R1 = CanonicalRank(f1, mat,ms);
//    INT_TYPE V1 = CanonicalRank(f1, vec, vs);
//    if ( R1*V1 == 0 ){
//     //   printf("%d %d %d\n", mat, f1->sinc.tulip[mat].name, name(f1,mat));
//     //   printf("%d %d %d\n", CanonicalRank(f1,mat,0), CanonicalRank(f1, name(f1,mat),0));
//
//        f1->sinc.tulip[alloy].Current[spin] = 0;
//        return 1;
//    }
//    double ALPHA = f1->mem1->rt->ALPHA;
//
//    if ( L1 == 0 ){
//        //   printf("CD: zero length %lld %lld %lld\n", origin,alloy,spin);
//        tId(f1, alloy, spin);
//        L1 = 1;
//        if (! CanonicalRank(f1, alloy, spin)){
//            printf("hurm... %d\n",alloy);
//            exit(0);
//        }
//    }
//    INT_TYPE M2[SPACE];
//    length(f1,alloy,M2);
//    double value2= 1;
//    INT_TYPE info;
//    double rcond;
//
//    {//REF CORE NUMBER
//        //GET CORE... ALLOCATE TO CORE-LANE.
//
//        for ( space = 0; space < SPACE ; space++){
//            array[space] =  streams(f1, canonicalBuffers, rank , space);
//            array2[space] =  streams(f1, canonicalBuffers, rank , space) + L1*L1;
//
//        }
//
//        if (  L1*L1 + L1*imax(V1,M2[0])  >  part(f1,canonicalBuffers)){
//#if 1
//            printf("mem prob\n %lld+ %lld  > %lld \n",R1*V1,L1*L1,part(f1,canonicalBuffers) );
//#endif
//            exit(0);
//        }
//
////            canonicalStore = canonicalBuffersB;
////            canonicalStore = canonicalBuffersBM;
//
//    }
//
//    enum division alloyBak;
//    if ( species (f1, alloy ) == vector && bodies(f1,alloy) == one){
//        alloyBak = trainVector;
//    }
//    else     if ( species (f1, alloy ) == vector && bodies(f1,alloy) == two){
//        alloyBak = trainVector2;
//    }
//    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == one){
//
//        f1->sinc.tulip[trainMatrix].header = header(f1,alloy);
//
//        alloyBak = trainMatrix;
//    }else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == two){
//
//        f1->sinc.tulip[trainMatrix2].header = header(f1,alloy);
//
//        alloyBak = trainMatrix2;
//    }
//    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == three){
//
//        f1->sinc.tulip[trainVector3].header = header(f1,alloy);
//
//        alloyBak = trainVector3;
//    }
//    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == four){
//
//        f1->sinc.tulip[trainVector4].header = header(f1,alloy);
//
//        alloyBak = trainVector4;
//    }
//    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == three){
//
//        f1->sinc.tulip[trainMatrix3].header = header(f1,alloy);
//
//        alloyBak = trainMatrix3;
//    }
//    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == four){
//
//        f1->sinc.tulip[trainMatrix4].header = header(f1,alloy);
//
//        alloyBak = trainMatrix4;
//    }
////    else if ( species(f1, alloy ) == quartic ){
////        f1->sinc.tulip[trainQuartic].header = header(f1,alloy);
////
////        alloyBak = trainQuartic;
//    else{
//        printf("cd: tObject species %d \n", alloy);
//
//        exit(0);
//
//    }
//    // printf(":: %d %d -- %d %d \n", origin , os, alloy ,spin );
//    space2 = 1;
////    if ( normalize(f1,alloy,spin,space2) )
////        return -1;
////
////    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, L1, M2[space2], 1.0, streams(f1, alloy,spin,space2), M2[space2], streams(f1,alloy,spin,space2), M2[space2], 0.0, array[space2], L1);
//
//
//    do{
//        if ( count % ns == 0 )
//            tEqua(f1, alloyBak,rank, alloy, spin );
//
//        for ( space0 = 0; space0 < SPACE ;space0++){
//
//            if ( SPACE == 3 ){
//                if ( space0 == 0 ){
//                    space0 = 0;
//                    space1 = 1;
//                    space2 = 2;
//                } else if ( space0 == 1 ){
//                    space0 = 1;
//                    space1 = 2;
//                    space2 = 0;
//                }else if ( space0 == 2 ){
//                    space0 = 2;
//                    space1 = 0;
//                    space2 = 1;
//                } else {
//                    printf("cd space");
//                    exit(0);
//                }
//            }
//            else if ( SPACE == 2 ){
//                if ( space0 == 0 ){
//                    space0 = 0;
//                    space1 = 1;
//                } else if ( space0 == 1 ){
//                    space0 = 1;
//                    space1 = 0;
//                }else {
//                    printf("cd space");
//                    exit(0);
//                }
//            }
//
//
//
//
//            if ( SPACE == 3 ){
//                //normal 3d calculator
//
//                if ( normalize(f1,alloy,spin,space2) )
//                    return -1;
//
//                if ( normalize(f1,alloy,spin,space1) )
//                    return -1;
//
//
//                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, L1, M2[space2], 1.0, streams(f1, alloy,spin,space2), M2[space2], streams(f1,alloy,spin,space2), M2[space2], 0.0, array[space2], L1);
//
//                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, L1, M2[space1], 1.0, streams(f1, alloy,spin,space1), M2[space1], streams(f1,alloy,spin,space1), M2[space1], 0.0, array[space1], L1);
//
//
//                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1*L1, 0,array[space2],1, array[space1],1 );
//                //replace Matrix
//
//                if ( L1 > 1 )
//                    for ( l = 0 ; l < L1; l++){
//                        if ( fabs(array[space1][ l*L1 + l ]-1.0) > 1e-6 ||  isnan ( array[space1][l*L1+l] )  ||  isinf ( array[space1][l*L1+l] ) ){
//#if VERBOSE
//                            printf("%lld:%f (%lld)-\n", l,array[space1][ l*L1 + l ],rank);
//                            fflush(stdout);
//#endif
//                        }
//                        array[space1][ l*L1 + l ] += ALPHA;
//                    }
//
//            }else if ( SPACE == 2 ){
//                //normal 3d calculator
//
//
//                if ( normalize(f1,alloy,spin,space1) )
//                    return -1;
//
//
//                if(0 ){
//                    spread(f1, nullName,0, alloy, spin, space1, array[space1], NULL);
//                } else {
//                    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, L1, M2[space1], 1.0, streams(f1, alloy,spin,space1), M2[space1], streams(f1,alloy,spin,space1), M2[space1], 0.0, array[space1], L1);
//                }
//
//                if ( L1 > 1 )
//                    for ( l = 0 ; l < L1; l++){
//                        if ( fabs(array[space1][ l*L1 + l ]-1.0) > 1e-6 ||  isnan ( array[space1][l*L1+l] )  ||  isinf ( array[space1][l*L1+l] ) ){
//#if VERBOSE
//                            printf("%lld:%f (%lld)-\n", l,array[space1][ l*L1 + l ],rank);
//                            fflush(stdout);
//#endif
//                        }
//                        array[space1][ l*L1 + l ] += ALPHA;
//                    }
//
//            }
//
//            info = tdpotrf(L1, array[space1]);
//
//            if ( info != 0 ){
//#if 1
//                printf("info failure %lld %lld %d %lld %d %lld\n", info, rank,mat,vec, alloy, spin);
//#endif
//                return 1;
//            }
//            rcond = tdpocon(rank, f1, L1, array[space1] );
////            if ( f1->mem1->rt->targetCondition > 0 )
//////
////            ALPHA /= min(1e1,max(1e-1, 1 - log(rcond*f1->mem1->rt->targetCondition) ));
//
//
//            if (  isnan(rcond) || isinf(rcond) || rcond < 1e-12){
//#if VERBOSE
//                printf("Warning: condition of Beylkin %d->%d is bad %16.16f\n", mat, alloy, rcond);
//                fflush(stdout);
//#endif
//                return 1;
//            }
//
//            INT_TYPE r,r2,ll,dimmer;
//            for ( r = 0; r < L1*M2[space0] ; r++)
//                array2[space0][r] = 0.;
////            assignView(4, f1, mat, 1);
//            if (name(f1, mat)< f1->sinc.vectorOperator){
//                for ( r = begin ; r < imin(R1,end) ; r++){
//                    f1->sinc.tulip[productVector].Current[rank] = 0;
//                    for (ll = 0 ;ll < CanonicalRank(f1, vec, vs);ll++ ){
//                            for ( dimmer = 0 ;dimmer < SPACE ; dimmer++)
//                                tMultiplyOne(rank, f1, dimmer, productVector, rank, mat, r, ms, vec, ll, vs);
//                        f1->sinc.tulip[productVector].Current[rank]++;
//                    }
////                    tMultiplyMP(rank, &info, f1, 1., -1, productVector, rank, mc,ocean(4,f1,r,ms),ms, 'N', vec, vs);
//                    {
//                        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, V1, M2[space1], 1.0,  streams(f1,alloy,spin,space1), M2[space1], streams(f1, productVector,rank,space1), M2[space1],0.0, array2[space1], L1);
//
//                        if ( SPACE > 2 ){
//                            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, V1, M2[space2], 1.0,  streams(f1,alloy,spin,space2), M2[space2], streams(f1, productVector,rank,space2), M2[space2],0.0, array2[space2], L1);
//
//                            cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,V1*L1, 0,array2[space2],1, array2[space1],1 );
//                        }
//                    }//form <g,f>
//
//                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,L1,M2[space0],V1,1.0,array2[space1],L1,streams(f1,productVector,rank,space0),M2[space0], 1.00000, array2[space0], L1 );
//                }
//            }else {
//                assignView(0, f1, mat, 1);
//                    for ( r = 0 ; r < R1 ; r++)
//                        for ( r2 = 0 ; r2 < R1 ; r2++){
//
//                            tMultiplyMP(rank, &info, f1, 1., -1, complement, rank, 'N',ocean(0,f1,r,ms),ms, 'N', vec, vs);
//                            tMultiplyMP(rank, &info, f1, 1., -1, productVector, rank, 'N',ocean(0,f1,r2,ms),ms, 'N', complement, rank);
//                            tScaleOne(f1, productVector,rank, -1);
//                            {
//                                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, V1, M2[space1], 1.0,  streams(f1,alloy,spin,space1), M2[space1], streams(f1, productVector,rank,space1), M2[space1],0.0, array2[space1], L1);
//
//                                if ( SPACE > 2 ){
//                                    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, V1, M2[space2], 1.0,  streams(f1,alloy,spin,space2), M2[space2], streams(f1, productVector,rank,space2), M2[space2],0.0, array2[space2], L1);
//
//                                    cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,V1*L1, 0,array2[space2],1, array2[space1],1 );
//                                }
//                            }//form <g,f>
//
//                            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,L1,M2[space0],V1,1.0,array2[space1],L1,streams(f1,productVector,rank,space0),M2[space0], 1.00000, array2[space0], L1 );
//                        }
//                }
//            info = tdpotrs(L1,  M2[space0], array[space1],  array2[space0] );
//            if ( info != 0 ){
//#if VERBOSE
//                printf("L1 %lld \n", L1);
//                printf("M2 %lld \n", M2[0]);
//                printf("R1 %lld \n", V1*L1);
//#endif
//                printf("info failture %lld\n", info);
//
//                exit(0);
//                return 1;
//            }
//
//            transpose ( L1, M2[space0] , array2[space0] , streams(f1,alloy,spin,space0));
//
//        }
//        balance(f1, alloy , spin );
//
//        if ( count % ns == 0 ){
//            value2 = distanceOne(rank, f1, alloy, spin,alloyBak, rank);
//        }
//        count++;
//    } while (  value2 >tolerance ) ;
//    return 0;
//
//}


//double tCycleMultiplyMP ( INT_TYPE rank,struct field * f1 , INT_TYPE begin, INT_TYPE end,char mc,enum division mat,INT_TYPE ms, enum division vec, INT_TYPE vs,   enum division alloy ,  INT_TYPE spin ,double tolerance, INT_TYPE maxRun, double power){
//    if ( ! CanonicalRank(f1, mat, ms ) )//likely has zer
//        return 0;
//
//    do{
//        if (canonicalMultiplyMP(rank, f1,begin,end, mc,mat,ms, vec, vs,alloy, spin, tolerance) ){
//#if VERBOSE
//            printf("Multiply  bailed %d %d %d %d -- %d\n",mat,ms,alloy,spin , CanonicalRank(f1, alloy, spin));
//            fflush(stdout);
//#endif
//            f1->sinc.tulip[alloy].Current[spin]--;
//            canonicalMultiplyMP(rank, f1,begin,end, mc,mat,ms, vec, vs,alloy, spin, tolerance);
//            return 0;
//        }
//        if ( CanonicalRank(f1, alloy, spin ) < maxRun )
//            tId(f1, alloy, spin);
//        else
//            return 0;
//    }while(1);
//    return 0.;
//}



INT_TYPE tOuterProductSu( struct field * f1,enum division vector , INT_TYPE a, enum division vector2,INT_TYPE b, enum division proj, INT_TYPE c){
    INT_TYPE ma = CanonicalRank(f1, vector,a), mb = CanonicalRank(f1, vector,b), zc = CanonicalRank(f1, proj, c),l,r,space,i;
    INT_TYPE n[SPACE];
    length(f1,vector,n);
    INT_TYPE N2[SPACE];
    INT_TYPE n2[SPACE];
    length(f1, a, n2);
    length(f1, proj, N2);
    if ( species(f1, vector ) == matrix || species(f1, vector2) == matrix){
        printf("outer\n");
        exit(0);
    }
    if ( zc + ma*mb > part(f1, proj) && proj < eigenVectors){
        printf("outerProductSu:  \n");
        exit(0);
    }
    
    for ( l = 0 ; l < ma; l++)
        for ( r = 0; r < mb; r++)
            for ( space = 0; space < SPACE ; space++)
                for ( i = 0; i < N2[space] ; i++)
                    (streams(f1, proj,c,space)+(l*mb+r+zc)*N2[space])[i] = 0.;
    
    for ( l = 0 ; l < ma; l++)
        for ( r = 0; r < mb; r++)
            for ( space = 0; space < SPACE; space++)
                cblas_dger(CblasColMajor, n[space],n[space], 1. , streams(f1, vector,a,space)+l*n2[space],1, streams(f1, vector2,b,space)+r*n2[space],1, streams(f1, proj,c,space)+(l*mb+r+zc)*N2[space],n[space]);
    f1->sinc.tulip[proj].Current[c] += ma*mb;
    return 0;
}

INT_TYPE tOuterProductSuOne( struct field * f1,INT_TYPE space,enum division vector , INT_TYPE a, enum division vector2,INT_TYPE b, enum division proj, INT_TYPE c){
    INT_TYPE ma = CanonicalRank(f1, vector,a), mb = CanonicalRank(f1, vector,b), zc = CanonicalRank(f1, proj, c),l,r,i;
    INT_TYPE N2 = alloc(f1, proj, space), n2 = alloc(f1, vector,space), n[SPACE];
  length1(f1,n);


    if ( species(f1, vector ) == matrix || species(f1, vector2) == matrix){
        printf("outer\n");
        exit(0);
    }
    if ( zc + ma*mb > part(f1, proj) && proj < eigenVectors){
        printf("outerProductSu:  \n");
        exit(0);
    }
    
    for ( l = 0 ; l < ma; l++)
        for ( r = 0; r < mb; r++)
                for ( i = 0; i < N2 ; i++)
                    (streams(f1, proj,c,space)+(l*mb+r+zc)*N2)[i] = 0.;
    
    for ( l = 0 ; l < ma; l++)
        for ( r = 0; r < mb; r++)
                cblas_dger(CblasColMajor, n[space],n[space], 1. , streams(f1, vector,a,space)+l*n2,1, streams(f1, vector2,b,space)+r*n2,1, streams(f1, proj,c,space)+(l*mb+r+zc)*N2,n[space]);
//    f1->sinc.tulip[proj].Current[c] += ma*mb;
    return 0;
}

INT_TYPE tGEMV (INT_TYPE rank,  struct field * f1, INT_TYPE space, enum division equals, INT_TYPE espin,enum division left,INT_TYPE l,INT_TYPE lspin, enum division right,INT_TYPE r, INT_TYPE rspin ){
    if ( header(f1, left ) != header(f1, right ) ){
        printf("Two Head types\n");
        exit(1);
    }
    f1->sinc.tulip[canonicalmvVector].Current[rank] = 0;
    f1->sinc.tulip[canonicalmv2Vector].Current[rank] = 0;
    f1->sinc.tulip[canonicalmv3Vector].Current[rank] = 0;
    if ( species(f1,left) == matrix&& species(f1, right ) == vector && f1->sinc.tulip[left].space[space].block != id0 ){
       char  in =  matrixAction( bodies(f1, right), f1->sinc.tulip[left].space[space].block,1);
       char  out = matrixAction( bodies(f1, right), f1->sinc.tulip[left].space[space].block,-1);
       // printf("(%d %d), %c %c\n",l,r ,in, out);

        if ( iv1 <= f1->sinc.tulip[left].space[space].block && f1->sinc.tulip[left].space[space].block <= iv4 ){
            tPermuteOne(rank, f1, space, in, right, r, rspin, canonicalmvVector, rank);
            tMultiplyOne(rank, f1, space, canonicalmv2Vector, rank, inversion, 0, 0, canonicalmvVector, 0, rank);
            tMultiplyOne(rank,  f1, space,canonicalmv3Vector, rank, left,l,lspin, canonicalmv2Vector,0, rank);
            tMultiplyOne(rank, f1, space, canonicalmvVector, rank, inversion, 0, 0, canonicalmv3Vector, 0, rank);
            tPermuteOne(rank, f1, space, out, canonicalmvVector, 0, rank, equals, espin);
        }else {
            tPermuteOne(rank, f1, space, in, right, r, rspin, canonicalmvVector, rank);
            tMultiplyOne(rank,  f1, space,canonicalmv2Vector, rank, left,l,lspin, canonicalmvVector,0, rank);
            tPermuteOne(rank, f1, space, out, canonicalmv2Vector, 0, rank, equals, espin);
        }
    }else {
        printf("%d\n", species(f1, right ) );
        printf("%d %d %d",species(f1,left) == matrix, species(f1, right ) == vector , f1->sinc.tulip[left].space[space].block != id0);
        printf("odd\n");
        exit(1);
    }
    return 0;
}

INT_TYPE tGEMM (INT_TYPE rank,  struct field * f1, INT_TYPE space, enum division equals, INT_TYPE espin,enum division left,INT_TYPE l,INT_TYPE lspin, enum division right, INT_TYPE r,INT_TYPE rspin ){
    
    if ( header(f1, left ) != header(f1, right ) ){
        printf("Two Head types\n");
        exit(1);
    }
    if ( species(f1,left) == matrix && species(f1, right ) == matrix){
        tMultiplyOne(rank,  f1, space,equals, espin, left,l,lspin, right,r, rspin);
    }
    return 0;
}

INT_TYPE tGEVV (INT_TYPE rank,  struct field * f1,  INT_TYPE space,enum division equals, INT_TYPE espin,INT_TYPE leftChar, enum division left,INT_TYPE l,INT_TYPE lspin, INT_TYPE rightChar, enum division right,INT_TYPE r, INT_TYPE rspin ){
    if ( header(f1, left ) != header(f1, right ) ){
        printf("Two Head types\n");
        exit(0);
    }
    f1->sinc.tulip[canonicalvvVector].Current[rank] = 0;
    f1->sinc.tulip[canonicalvv2Vector].Current[rank] = 0;
    f1->sinc.tulip[canonicalvv3Vector].Current[rank] = 0;
    printf("head\n");
    if ( species(f1,left) == outerVector&& species(f1, right ) == vector){
        tPermuteOne(rank, f1, space, rightChar, right, r, rspin, canonicalvvVector, rank);
        if ( leftChar != 'T' ){
            printf("!T");
            
            tPermuteOne(rank, f1, space, leftChar, left, l, lspin, canonicalvv2Vector, rank);
            tMultiplyOne(rank,  f1, space,equals, espin, canonicalvv2Vector,0,rank, canonicalvvVector,0, rank);
        } else {
            printf("T");
            tMultiplyOne(rank,  f1, space,equals, espin, left,l,lspin, canonicalvvVector,0, rank);
        }
    }
    return 0;
}


//could parallelize.
double tDOT (INT_TYPE rank,  struct field * f1,INT_TYPE dim,char leftChar, enum division left,INT_TYPE l,INT_TYPE lspin, char rightChar, enum division right,INT_TYPE r, INT_TYPE rspin ){
    INT_TYPE space = dim;
    double prod = 0.;
    f1->sinc.tulip[canonicaldotVector].Current[rank] = 0;
    f1->sinc.tulip[canonicaldot2Vector].Current[rank] = 0;
    f1->sinc.tulip[canonicaldot3Vector].Current[rank] = 0;

        if ( leftChar != 'T' && rightChar != 'N'){
            tPermuteOne(rank, f1, space, rightChar, right, r, rspin, canonicaldotVector, rank);
            tPermuteOne(rank, f1, space, leftChar, left, l, lspin, canonicaldot2Vector, rank);
            prod = tMultiplyOne(rank,  f1, space,nullScalar, 0, canonicaldot2Vector,0,rank, canonicaldotVector,0, rank);
            
        } else if ( leftChar == 'T' && rightChar != 'N'){
            tPermuteOne(rank, f1, space, rightChar, right, r, rspin, canonicaldotVector, rank);
            prod = tMultiplyOne(rank,  f1, space,nullScalar, 0, left,l,lspin, canonicaldotVector,0, rank);
        }else if ( leftChar != 'T' && rightChar == 'N'){
                tPermuteOne(rank, f1, space, leftChar, left, l, lspin, canonicaldot2Vector, rank);
                prod = tMultiplyOne(rank,  f1, space,nullScalar, 0, canonicaldot2Vector,0,rank, right,r, rspin);
        }else if ( leftChar == 'T' && rightChar == 'N'){
                prod = tMultiplyOne(rank,  f1, space,nullScalar, 0, left,l,lspin, right,r, rspin);
        }
    
    return prod;
}


double tMultiplyMP (INT_TYPE rank, INT_TYPE * info, struct field * f1,double number, INT_TYPE beta,  enum division equals, INT_TYPE espin ,char leftChar, enum division left,INT_TYPE lspin, char rightChar,enum division right, INT_TYPE rspin ){
    INT_TYPE l,r,dim;
    double sum,prod;
    if ( beta == -1 )
        f1->sinc.tulip[equals].Current[espin] = 0;
    *info = 0;
    enum genus a = species(f1,left), b = species(f1,right);
    if ( number != 1.0 )
        exit(10);

    sum = 0.;
    for ( l = 0 ; l < CanonicalRank(f1, left, lspin) ; l++)
        for ( r = 0 ; r < CanonicalRank(f1, right, rspin) ; r++){
            prod = 1.;
            for ( dim = 0; dim < SPACE ; dim++)
            if ( f1->sinc.rose[dim].body != nada )
            {
                if ( species(f1,equals)== scalar)
                    prod *= tDOT(rank, f1, dim,leftChar, left,l, lspin, rightChar, right,r, rspin);
                else  if ( species(f1,left) == vector&& species(f1, right ) == vector){
                    tGEVV(rank, f1, dim, equals, espin,leftChar, left, l,lspin, rightChar, right,r, rspin);
                }
                else if ( species(f1,left) == matrix && species(f1, right ) == matrix){
                    tGEMM(rank, f1, dim, equals, espin, left, l,lspin, right,r, rspin);
                    
                } else     if ( species(f1,left) == matrix&& species(f1, right ) == vector){
                    tGEMV(rank, f1, dim, equals, espin, left,l, lspin, right,r, rspin);
                }
            }
            sum += prod ;
            f1->sinc.tulip[equals].Current[espin]++;

        }
    if ( species(f1,equals)== scalar ){
        return number * sum ;
    }
    return 0.;
}

double tMultiplyOne (INT_TYPE rank, struct field * f1,INT_TYPE space,  enum division equals,INT_TYPE espin , enum division left,INT_TYPE l,INT_TYPE lspin, enum division right,INT_TYPE r, INT_TYPE rspin){
    INT_TYPE LN2[SPACE],RN2[SPACE],EN2[SPACE] ;
    if ( f1->sinc.rose[space].body == nada)
        return 0.;
    INT_TYPE N1[SPACE];
    length(f1,left,LN2 );
    length(f1,right,RN2 );
    length(f1,equals,EN2 );
    INT_TYPE cur = CanonicalRank(f1, equals, espin);
    //0: scalar -->   dot
    //leading dimensions contract with leading dimensions...(comma repected)
    //1:                                        i,j | j > = i   (gemv)
    //2: vector from either matrix * vector :    i,j |jj'> = ij' (gemm)
    //3: vector form vector * vector :   (insist) ik | i> = k  (gemv)
    //4: matrix from matrix * matrx :            i,j * j,k = i,k (gemm)
    

    if (  (species(f1,equals) == scalar && (species(f1, left) == species(f1, right)|| ((species(f1, left) + species(f1, right)) == vector +outerVector)) && bodies(f1,left) == bodies(f1,right))){
      
//        if ( left == eigenVectors && right == canonicalmeVector ){
//        printf("((((");
//        outputFormat(f1, stdout, left, lspin);
//        outputFormat(f1, stdout, right, rspin);
//        printf("))))");
//        }
        
        double va = cblas_ddot( LN2[space] , streams(f1,left,lspin,space)+l*LN2[space],1 , streams(f1,right,rspin,space)+r*RN2[space], 1);
        
        return va;
    }else
                if (  (species(f1, equals) == vector && species(f1, left) == matrix && species(f1, right) == vector && bodies(f1,left) == bodies(f1,right)))
                {

                    cblas_dgemv( CblasColMajor, CblasNoTrans,  EN2[space], EN2[space],1.,
                                streams( f1, left, lspin,space )+l*LN2[space], EN2[space],
                                streams(f1, right, rspin,space)+r*RN2[space],1, 0.,
                                streams( f1, equals, espin,space )+cur*EN2[space], 1  );
                    
                }else if (  (species(f1, equals) == vector && species(f1, left) == matrix && species(f1, right) == vector && bodies(f1,left) < bodies(f1,right)))
                {
                    length(f1,left,N1);
                    {
                        INT_TYPE i;
                    for ( i = 0 ; i < SPACE ; i++)
                        N1[i] = sqrt(N1[i]);
                    }
                    if ( RN2[space]%N1[space] != 0 ){
                        printf("watch mods R\n");
                        exit(1);
                    }
                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,N1[space],RN2[space]/N1[space],N1[space],1.,streams( f1, left, lspin,space )+l*LN2[space],N1[space],streams(f1, right, rspin,space)+r*RN2[space],N1[space], 0., streams( f1, equals, espin,space)+cur*EN2[space],N1[space]);

                    
                }else if (  (species(f1, equals) == vector && species(f1, left) == matrix && species(f1, right) == vector && bodies(f1,name(f1,left)) > bodies(f1,name(f1,right))))
                {
                    struct name_label l = f1->sinc.tulip[left];
                    struct name_label r = f1->sinc.tulip[right];
                    struct name_label e = f1->sinc.tulip[equals];

                    printf("body count\n");
                    exit(0);
                }
                else if ( species(f1, equals) == outerVector && species(f1, left)  ==outerVector && vector == species(f1, right)){
                    length(f1, right,N1);
                    printf("here\n");
                    if ( bodies(f1, left ) < bodies(f1,right)){
                        printf("complecationx\n");
                        exit(1);
                    }
                  //      return tMultiplyOne(rank, f1, space, equals,  espin, right, r, rspin, left, l, lspin);
                    
                    //( bodies(f1, left ) >= bodies(f1,right))
                    if ( LN2[space]%N1[space] != 0 ){
                        printf("watch mods L\n");
                        exit(1);
                    }
                    cblas_dgemv( CblasColMajor, CblasTrans,  LN2[space]/N1[space], N1[space],1.,
                                streams( f1, left, lspin,space )+l*LN2[space], N1[space],
                                streams(f1, right, rspin,space)+r*RN2[space],1, 0.,
                                streams( f1, equals, espin,space )+cur*EN2[space], 1  );

                }
                else if (  (species(f1, equals) == matrix && species(f1, left)  ==matrix && matrix == species(f1, right)&& bodies(f1,left) == bodies(f1,right))){
                    length(f1, eigenVectors,N1);

                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,N1[space],N1[space],N1[space],1.,streams( f1, left, lspin,space )+l*LN2[space],N1[space],streams(f1, right, rspin,space)+r*RN2[space],N1[space], 0., streams( f1, equals, espin,space)+cur*EN2[space],N1[space]);
                }
    else   {
        struct name_label n = f1->sinc.tulip[equals];

        printf("tM: wrong type %d=  %d * %d \n", equals, left, right);
        exit(0);
    }
    return 0.;
}



void tHYpY(  INT_TYPE rank, struct field * f1 ,INT_TYPE targSpin, enum division left, INT_TYPE l, double product, double productCmpl, enum division ket , INT_TYPE k, enum division oket, INT_TYPE ospin ){
    
    if ( left == nullName || species(f1, left ) == scalar)
        return;
    INT_TYPE N2;
    double prod ;
    INT_TYPE dim,im,sp2,im2,ll1,ll2;
    DCOMPLEX co2,coi,coi2,pro = product + I * productCmpl;
    if ( species(f1, left ) != outerVector ){
        co2 = 1.;
        for ( sp2 = 0 ; sp2 < spins(f1, ket);sp2++){
            if (sp2 == 1)
                co2 *= I;
            if ( k >= CanonicalRank(f1, ket, sp2 ))
                continue;
            coi=1.;
            for ( im = 0 ; im < spins(f1, left);im++){
                if (im == 1 )
                    coi*= I;
                if ( l >= CanonicalRank(f1, name(f1,left), im)  )
                    continue;
                if (( fabs(creal(co2 * coi * pro)) > f1->mem1->rt->TARGET && targSpin == 0 )|| ( fabs(cimag(co2 * coi * pro)) > f1->mem1->rt->TARGET && targSpin == 1 )){
                    for ( dim = 0 ; dim < SPACE ; dim++)
                        if ( f1->sinc.rose[dim].body != nada){
                            N2 = alloc(f1, ket, dim);
                            if ( f1->sinc.tulip[left].space[dim].block == id0 ){
                                xsAdd(1.,dim, f1, oket, ospin, f1, ket, k, sp2);
                            }
                            else {
                                tGEMV(rank, f1, dim,oket, ospin, left, l, im, ket, k, sp2);
                            }//both procedues will add to end of oket no incremen
                            if ( targSpin == 0 )
                                prod = creal(co2 * coi*pro );
                            else
                                prod = cimag(co2 * coi*pro );
                            cblas_dscal(N2, prod, streams(f1,oket,ospin, dim)+CanonicalRank(f1, oket, ospin)*N2, 1);
                        }
                    f1->sinc.tulip[oket].Current[ospin] += 1;//potential total addition of 4
                }
                
                
            }
            
        }
    }
    else {
        
        ll1 = l % part(f1, left);
        ll2 = (l/part(f1,left)) % part(f1, left);

        co2 = 1.;
        for ( sp2 = 0 ; sp2 < spins(f1, ket);sp2++){
            if (sp2 == 1)
                co2 *= I;
            if ( k >= CanonicalRank(f1, ket, sp2 ))
                continue;

            coi=1.;
            for ( im = 0 ; im < spins(f1, left);im++){
                if (im == 1 )
                    coi*= I;
                coi2 = 1.;
                if ( ll1 >= CanonicalRank(f1, left, im ))
                    continue;
                for ( im2 = 0 ; im2 < spins(f1, left);im2++){
                    if (im2 == 1 )
                        coi2*= -I;
                    if ( ll2 >= CanonicalRank(f1, left, im2 ))
                        continue;
                    
                    if (( fabs(creal(co2 * coi*coi2*pro)) > f1->mem1->rt->TARGET && targSpin == 0 )|| ( fabs(cimag(co2 * coi*coi2*pro )) > f1->mem1->rt->TARGET && targSpin == 1 )){
                        for ( dim = 0 ; dim < SPACE ; dim++)
                            if ( f1->sinc.rose[dim].body != nada){
                                N2 = alloc(f1, ket, dim);

//                                if ( f1->sinc.tulip[left].space[dim].block == id0 ){
//                                    xsAdd(1.,dim, f1, oket, ospin, f1, ket, k, sp2);
//                                }
//                                else
                                {
                                
                                    enum division c1 = complement + (f1->sinc.rose[dim].body - f1->sinc.tulip[left].space[dim].body);
                                    if ( c1 < complement ){
                                        printf("wrong dimensionality %d ", c1 - complement);
                                        exit(0);
                                    }

                                    if ( c1 != complement ){
                                        f1->sinc.tulip[c1].Current[rank] = 0.;
                                        tGEVV(rank, f1, dim,c1, rank, 'N',left, ll1, im,'N', ket, k, sp2);
                                        tGEVV(rank, f1, dim,oket, ospin, 'T',left, ll2, im2,'N', c1, 0, rank);
                                    }else {
                                        xsAdd(tDOT(rank, f1, dim, 'T', left, ll1, im, 'N', ket, k, sp2), dim, f1, oket, ospin, f1, left, ll2, im2);
                                    }
                                }//both procedues will add to end of oket...but do not increment
                                
                                if ( targSpin == 0 )
                                    prod = creal(co2 * coi*coi2*pro );
                                else
                                    prod = cimag(co2 * coi*coi2*pro );
                                
                                cblas_dscal(N2, prod, streams(f1,oket,ospin, dim)+CanonicalRank(f1, oket, ospin)*N2, 1);
                            }
                        f1->sinc.tulip[oket].Current[ospin] += 1;//potential total addition of 8

                    }
                    
                    
                }
                
            }
        }
    }
    return ;
}
    
void tHXpX (  INT_TYPE rank, struct field * f1 , enum division left,INT_TYPE shiftFlag, double product, double productCmpl,  enum division right ,  double tolerance , INT_TYPE maxRun  ){
    INT_TYPE l , k,targSpin;
    enum division pt,Mat;
    struct name_label *fz = &f1->sinc.tulip[totalFuzzyVector];
    struct name_label *rg = &f1->sinc.tulip[right];

    struct name_label nm ;
    f1->sinc.tulip[totalVector].Current[rank] = 0;
    f1->sinc.tulip[totalFuzzyVector].Current[rank] = 0;

    if ( right == productVector || right == totalVector){
        printf("oops\n");
        exit(0);
    }
    //cpy in if shifted.
    for ( targSpin = 0 ; targSpin < spins(f1, right ) ;targSpin++){
        if ( shiftFlag && targSpin == 0 )
            tEqua(f1, totalFuzzyVector, rank, right, targSpin);
        pt = left;
        //cycle over all links in left.
        do {
            nm = f1->sinc.tulip[pt];
#if VERBOSE
            Mat = pt;
            printf("%d-%d\t%d\t%d\t %f || \t",Mat,name(f1,Mat),bodies(f1, Mat),CanonicalRank(f1, name(f1,Mat), cmpl),traceOne(f1, name(f1,Mat), cmpl));
            INT_TYPE space;
            printf(":");
            for ( space = 0; space < SPACE ; space++)
                printf("%d:", f1->sinc.tulip[Mat].space[space].block );
            printf("\n");
            fflush(stdout);
#endif
            for ( k = 0 ; k < part ( f1, right ) ; k++)
                for ( l = 0; l < part(f1, name(f1,pt) ) ; l++)

                {
            //        printf("%d %d %d\n", pt,CanonicalRank(f1, totalFuzzyVector, 0),CanonicalRank(f1, totalVector, 0));

                tHYpY(rank, f1, targSpin,pt, l, product, productCmpl, right, k, totalVector, rank);
                struct name_label t = f1->sinc.tulip[totalVector];
                {
                    if ( CanonicalRank(f1, totalVector, rank) >= part (f1, totalVector )-8){
                        tEqua(f1, productVector, rank, right, targSpin);
                        tCycleDecompostionListOneMP(rank, f1, totalVector, rank,NULL, productVector, rank, tolerance, maxRun, 1.);
                        tAddTw(f1, totalFuzzyVector, rank, productVector, rank);
                        f1->sinc.tulip[totalVector].Current[rank] = 0;
                    }
                }
            }
            pt = f1->sinc.tulip[pt].linkNext;
            
        }while ( pt != nullName );
        
        tEqua(f1, productVector, rank, right, targSpin);
        tCycleDecompostionListOneMP(rank, f1, totalVector, rank, NULL,productVector, rank, tolerance, maxRun, 1.);
        tAddTw(f1, totalFuzzyVector, rank, productVector, rank);

        tCycleDecompostionListOneMP(rank, f1, totalFuzzyVector, rank,NULL, right, targSpin, tolerance, maxRun, 1.);
        
    }
    struct name_label r = f1->sinc.tulip[right];

    return;
}

double distanceOne(INT_TYPE rank,struct field * f1 , enum division alloy , INT_TYPE spin , enum division alloyBak, INT_TYPE spin2){
    double value,value2;
    INT_TYPE info;
    value = 0.5*( tMultiplyMP(rank,&info,f1,1., -1,nullName, 0 ,'T', alloy, spin,'N', alloy ,spin));
    value += 0.5*(tMultiplyMP(rank,&info,f1,1., -1,nullName, 0 ,'T', alloyBak, spin2,'N', alloyBak ,spin2));
    value2 = 2*fabs(value-tMultiplyMP(rank,&info,f1,1., -1,nullName, 0 , 'T', alloyBak, spin2,'N', alloy ,spin) );
    return (fabs(value2));
    
}

double inner(INT_TYPE rank,struct field * f1 , enum division alloy , INT_TYPE spin ){
    INT_TYPE info;
    return tMultiplyMP(rank,&info,f1,1., -1,nullName, 0 ,'T', alloy, spin,'N', alloy ,spin);;
}

double magnitude ( struct field * f1 , enum division alloy ){
   // printf ("%f \t %f\n", inner ( 0 , f1, alloy , 0 ),inner ( 0 , f1, alloy ,1 ));
    double sum = 0.;
    INT_TYPE i;
    for ( i = 0; i < spins(f1, alloy); i++)
        sum +=  inner ( 0 , f1, alloy , i );
    return sqrt(sum);
}



INT_TYPE ready ( struct calculation * c ){
    INT_TYPE readyMemory = 1;
    INT_TYPE readyVector = 1;
    INT_TYPE space;
    if ( ! c->mem.bootedMemory || c->i.c.sinc.tulip == NULL )
        readyMemory = 0;
    
    if ( readyMemory )
        for ( space = 0 ; space <= SPACE ; space++)
            if ( c->i.c.sinc.rose[space].stream == NULL )
            readyMemory = 0;
    
    
    if ( readyMemory )
        if ( CanonicalRank(&c->i.c, eigenVectors , 0 ) == 0 ){
            printf("passing over stage because vector is null\n");
            readyVector = 0;
        }
    
    return readyVector && readyMemory;
}




INT_TYPE xConstructFoundation (struct calculation * calc , enum division usr, INT_TYPE UR, struct calculation * calc2, enum division usz, INT_TYPE UZ ,INT_TYPE mx){
    INT_TYPE *f = (INT_TYPE*)myStreams(&calc2->i.c, dsyBuffers, 0);
    INT_TYPE final,mdi,ii,iii,iv,rank = 0,cmpl,i;
    
    ii = 0;
    iii= 0;
    iv = 0;
    for ( i = 0; i < UZ ; i++)
        if ( CanonicalRank(&calc2->i.c, usz+i, 0)+CanonicalRank(&calc2->i.c, usz+i, 1))
        //if ( calc2->i.c.sinc.tulip[usz+i].value.symmetry == calc2->i.irrep|| ! calc2->i.irrep)
        {
//            fflush(stdout);
            if ( mx <= (CanonicalRank(&calc2->i.c,usz+i,0)+CanonicalRank(&calc2->i.c,usz+i,1)) ){
                f[ii++]  = i ;

                iii += mx ;
                iv = imax(iv , (CanonicalRank(&calc2->i.c,usz+i,0)+CanonicalRank(&calc2->i.c,usz+i,1))/mx + !(!((CanonicalRank(&calc2->i.c,usz+i,0)+CanonicalRank(&calc2->i.c,usz+i,1)%mx ))) ) ;
           //    printf("%d -%d - %d - %d\n",iv, iii,ii,i);

            }
        }
    
    if ( ! UR )
        return iii;
    if ( UR == -1 )
        return iv;

    if ( iii != UR ){
    }

    for ( i = 0; i < ii ; i++){
        for ( cmpl = 0; cmpl < spins(&calc2->i.c,usz) ; cmpl++)
        {
            xEqua(&calc->i.c, usr+i, cmpl, &calc2->i.c, usz+f[i], cmpl);
   
#if VERBOSE
            INT_TYPE info;
for ( iii = 0 ; iii <= i ; iii++)
    printf("%d %d %1.15f\n", i,iii, tMultiplyMP(0, &info, &calc->i.c, 1., -1, nullName, 0, 'T', usr+i, 0, 'N', usr+iii, 0));
#endif
        }
        
    }
    final = 0;
    for ( i = 0; i < ii*mx ; i++){
        double value;
        
        value = magnitude(&calc->i.c, usr+i);
        if ( value > 0. ){
            tScale(&calc->i.c , usr+i,1./value);
            final++;
        }
        else
            tClear(&calc->i.c,usr+i);
    }
#if VERBOSE
    printf("complete transfer %d\n", ii*mx);
#endif
    fflush(stdout);
    return final;
}
