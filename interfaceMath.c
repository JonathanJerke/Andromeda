/*
 *  interfaceMath.c
 *
 *
 *  Copyright 2020 Jonathan Jerke and Bill Poirier.
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

#include "interfaceMath.h"



double Sinc( double d , double x){
    double dimless , arg = pi*(x/d);
    if(  fabs(arg)< 0.001 ){
        dimless = 1 - 1./6. * arg * arg + 1. /120. * arg * arg * arg * arg ;
    } else {
        dimless = sin( arg )/arg;
    }
    return dimless;
}

DCOMPLEX periodicSinc ( double d , double x, double momentum, INT_TYPE N1 ){
    double dimless =1., arg = pi*(x/d);
    
    
    if ( arg != 0.0 )
        dimless = sin( arg )/2./N1 *(1./tan(arg/2./N1)+tan(arg/2./N1));
//    if (  isnan(dimless) ||  isinf(dimless))
//        return 1.;
    
    
    return dimless * ei( momentum * d * x );
}

double SS( double d1 , double x , double d2, double y )    {
    if ( d1 > d2 ){
        if ( x == y )
            return sqrt(d2/d1);
        else
            return Sinc(d1,x-y) * sqrt(d2/d1);
    }
    if ( x == y )
        return sqrt(d1/d2);
    else
        return Sinc(d2,x-y) * sqrt(d1/d2);
    return 0;
}

//double periodicSS( double d1 , double x ,INT_TYPE N1, double d2, double y , INT_TYPE N2)    {
//    if ( d1 > d2 ){
//        if ( x == y )
//            return sqrt(d2/d1);
//        else
//            return periodicSinc(d1,x-y,N1) * sqrt(d2/d1);
//    }
//    if ( x == y )
//        return sqrt(d1/d2);
//    else
//        return periodicSinc(d2,x-y,N2) * sqrt(d1/d2);
//    return 0;
//}

//double OVERLAP ( struct basisElement *b1,  struct basisElement *b2){
//    if ( b1->basis == SincBasisElement && b2->basis == SincBasisElement)
//        return SS(b1->length, b1->length*b1->index - b1->origin, b2->length, b2->length*b2->index - b2->origin);
//    else if ( b1->basis == periodicSincBasisElement && b2->basis == periodicSincBasisElement)
//            return periodicSS(b1->length, b1->length*b1->index - b1->origin, b1->auxIndex,b2->length, b2->length*b2->index - b2->origin,b2->auxIndex);
//    else
//        return 0;
//}
//
//INT_TYPE dimBasisElement  (struct field * f1, INT_TYPE body ){
//    INT_TYPE ibe,i,ct;
//    if ( body > f1->body )
//        return 0;
//    INT_TYPE count[SPACE] ;
//    for ( i = 0;i < SPACE ; i++)
//        count[i] = 0;
//    
//    for ( ibe = 0 ; ibe < f1->boaL   ; ibe++)
//    {
//        if ( f1->boa[ibe].body == body )
//            count[f1->boa[ibe].dim]++;
//        
//    }
//    ct = 0;
//    for ( i = 0;i < SPACE ; i++)
//        if ( count[i] )
//            ct++;
//
//    return ct;
//    
//}
//
//INT_TYPE assocBasisElement (struct field * f1 ){
//    INT_TYPE count[MAXASSOC];
//    INT_TYPE ibe,i,ct;
//    for ( i = 0;i < MAXASSOC ; i++)
//        count[i] = 0;
//    
//    for ( ibe = 0 ; ibe < f1->boaL   ; ibe++)
//    {
//        {
//            if (f1->boa[ibe].association < 0 || f1->boa[ibe].association >= MAXASSOC )
//                exit(0);
//            else
//                count[f1->boa[ibe].association]++;
//        }
//    }
//    ct = 0;
//    for ( i = 0;i < MAXASSOC ; i++)
//        if ( count[i] )
//            ct++;
//    return ct;
//    
//}
//INT_TYPE lenBasisElement   (struct field * f1, INT_TYPE body ,INT_TYPE dim, INT_TYPE assoc){
//    INT_TYPE count=0;
//    INT_TYPE ibe;
//    if ( body > f1->body )
//        return 0;
//    for ( ibe = 0; ibe < f1->boaL; ibe++)
//    {
//        if ( f1->boa[ibe].body == body && f1->boa[ibe].association == assoc && f1->boa[ibe].dim == dim ){
//            count++;
//        }
//    }
//    return count;
//}
//
//struct basisElement * getBasisElement (struct field * f1, INT_TYPE body ,INT_TYPE dim , INT_TYPE assoc, INT_TYPE l ){
//    INT_TYPE count=0;
//    INT_TYPE ibe;
//    if ( body > f1->body )
//        return 0;
//    for ( ibe = 0; ibe < f1->boaL; ibe++)
//    {
//        if ( f1->boa[ibe].body == body && f1->boa[ibe].association == assoc && f1->boa[ibe].dim == dim ){
//            if ( l == count ){
//                return &f1->boa[ibe];
//            }
//            count++;
//        }
//    }
//    return NULL;
//}
//
//void arrayBasis(struct field * f1,INT_TYPE *B ){
//    INT_TYPE ibe,i,space,ct;
//    for ( space = 0; space < SPACE ; space++)
//        B[space] = 0;
//    for ( ibe = 0 ; ibe < f1->boaL   ; ibe++)
//    {
//        B[f1->boa[ibe].dim]++;
//    }
//    return;
//}


void transpose(INT_TYPE N, INT_TYPE M, Stream_Type * orig, Stream_Type* targ){
#if VERBOSE
    printf("transpose\n");
    fflush(stdout);
#endif

#ifndef MKL
    INT_TYPE i ,j;
    for ( i = 0 ; i < N ; i++)
        for ( j = 0; j < M ; j++){
            targ[i*M+j] = orig[j*N+i];
        }
#else
    mkl_domatcopy ('C', 'T',N, M,1.,orig,N,targ,M);
#endif
#if VERBOSE
    printf("transpose-end\n");
    fflush(stdout);
#endif

}

INT_TYPE tdsyev( INT_TYPE rank, struct sinc_label f1, char job , INT_TYPE n, double * ar, INT_TYPE ns , double * w ){
#ifdef APPLE

    char charU = 'U';
    INT_TYPE info,lbuffer = part(f1, dsyBuffers);
    dsyev_ ( &job, &charU, &n , ar , &ns , w , myStreams(f1, dsyBuffers,rank ), &lbuffer, &info );
    return info;
#else
    return LAPACKE_dsyev(LAPACK_COL_MAJOR, job, 'U', n, ar, ns , w);
#endif
}

INT_TYPE tzheev( INT_TYPE rank, struct sinc_label f1, char job , INT_TYPE n, DCOMPLEX * ar, INT_TYPE ns , double * w ){
    INT_TYPE i;

#ifdef APPLE
    char charU = 'U';
    INT_TYPE info,lbuffer = part(f1, dsyBuffers);
  zheev_ (&job, &charU, &n , (__CLPK_doublecomplex *)ar , &ns , w ,(__CLPK_doublecomplex *) (myStreams(f1, dsyBuffers,rank )+3*n), &lbuffer,  myStreams(f1, dsyBuffers,rank ),&info );
    return info;
#else
//    for ( i = 0 ; i < n*ns ;i++)
//        if ( isinf ( ar[i]) || isnan(ar[i]) ){
//            printf ("error with inptu to heev\n");
//            exit(0);
//        }
    
    if ( n == 0 ||ns == 0 )
    {
        printf ( "empty lists \n");
        exit(0);
    }
    return LAPACKE_zheev(LAPACK_COL_MAJOR, job, 'U', n, (DCOMPLEX_PRIME *)ar, ns , w );
#endif
}


INT_TYPE tdsygv( INT_TYPE rank, struct sinc_label f1, char job , INT_TYPE n, double * sr, double * ar, INT_TYPE ns , double * w ){
    INT_TYPE info;
    INT_TYPE type = 1;
    char charU = 'U';
#ifdef APPLE
    INT_TYPE lbuffer = part(f1, dsyBuffers);
    dsygv_(&type, &job, &charU,&n,sr,&ns,ar,&ns, w,myStreams(f1, dsyBuffers,rank ),&lbuffer, &info );
#else
    info =  LAPACKE_dsygv(LAPACK_COL_MAJOR,1, job, 'U', n,sr,ns, ar, ns , w );
#endif
    return info;

}

                     
INT_TYPE tzhegv( INT_TYPE rank, struct sinc_label f1, char job , INT_TYPE n,DCOMPLEX * sr, DCOMPLEX * ar, INT_TYPE ns , double * w ){
    INT_TYPE info;
//    INT_TYPE liwork = 3+5*n;
//    INT_TYPE lrwork = 1+5*n+2*n*n;
//    INT_TYPE lwork = 2*n+n*n;
//    INT_TYPE iwork[liwork];
//    double rwork[lrwork];
//    DCOMPLEX work[lwork];
    INT_TYPE type = 1;
    char charU = 'U';
#ifdef APPLE
    INT_TYPE i;
//    for ( i = 0 ; i  < ns * n ; i++){
//        if (( isnan(creal(sr[i])) || isnan(cimag(sr[i])))|| ( isinf(creal(sr[i])) || isinf(cimag(sr[i]))))
//          //  printf("%f %f\n", sr[i]);
//    }
//    for ( i = 0 ; i  <  n ; i++){
//        if( ( isnan(creal(sr[i])) || isnan(cimag(sr[i])))||( isinf(creal(sr[i])) || isinf(cimag(sr[i]))))
//
//       // printf("%f\n", cimag(ar[i*ns+i]));
//    }
//    if (0){
    INT_TYPE lbuffer = (part(f1,dsyBuffers)-10*n)/2;
    zhegv_(&type, &job, &charU,&n,(__CLPK_doublecomplex *)sr,&ns,(__CLPK_doublecomplex *)ar,&ns, w,(__CLPK_doublecomplex *) (myStreams(f1, dsyBuffers,rank ))+10*n, &lbuffer,  myStreams(f1, dsyBuffers,rank ),&info );
//    } else {
//        zhegvd_(&type, &job , &charU, &n,(__CLPK_doublecomplex *)sr,&ns,(__CLPK_doublecomplex *)ar,&ns, w, (__CLPK_doublecomplex *)(work), &lwork, rwork, &lrwork, iwork, &liwork, &info);
//    }
#else
    info =  LAPACKE_zhegv(LAPACK_COL_MAJOR,1, job, 'U', n,(DCOMPLEX_PRIME *)sr,ns,(DCOMPLEX_PRIME *) ar, ns , w );
#endif
    
    return info;
}
                               

                              
INT_TYPE tdgeqr( INT_TYPE rank, struct sinc_label f1,INT_TYPE len, INT_TYPE n, double * ar, INT_TYPE ns ,double *w){
    INT_TYPE info,info2;
    if ( len > n ){
        printf("%d -> %d\n", len, n);
        len = n;
    }else if ( len <= 0 ){
        printf("cancel tdgeqr\n");
        return -1;
    }
    //len is number of elements to consider, if len > n, then you may as well remove some of len...its over-linear dependent.
#if VERBOSE
    printf("eqr\n");
    fflush(stdout);
#endif

    printf("tdgeqr %d %d %d %d\n", rank, len, n, ns );
#ifdef APPLE
    INT_TYPE lbuffer = part(f1, dsyBuffers);
    dgeqrf_(&n,& len, ar, &ns, w, myStreams(f1, dsyBuffers,rank ), &lbuffer, &info);
    dorgqr_(&n, &len, &len, ar, &ns, w, myStreams(f1, dsyBuffers,rank ), &lbuffer, &info2);
#else
    info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, n,len, ar, ns, w);
    info2 = LAPACKE_dorgqr(LAPACK_COL_MAJOR,n,len,len,ar,ns,w);
#endif
    return info * 100000 + info2;

}


INT_TYPE tdpotrf ( INT_TYPE L1, double * array ) {
    INT_TYPE info;
    char charU = 'U';
    
#if VERBOSE
    printf("rtf\n");
    fflush(stdout);
#endif

#ifdef APPLE
        dpotrf_( &charU , &L1, array, &L1, &info ) ;

#else
info =  LAPACKE_dpotrf(LAPACK_COL_MAJOR,'U',L1,  array, L1 );
#endif
    if ( info != 0 ){
        printf("POTRF %d\n", info);
        INT_TYPE i ;
        for ( i = 0 ; i < L1 *L1 ; i++)
            printf("potrf %f\n", array[i]);
        exit(1);
    }
    return info;
}

INT_TYPE tdgels ( INT_TYPE rank,struct sinc_label  f1 , INT_TYPE L1, INT_TYPE M2, double * array, double * arrayo ,INT_TYPE inc){
    char charT = 'N';
    INT_TYPE lbuffer = part(f1, dsyBuffers),info=0;
    dgels_(&charT, &L1, &L1,& M2, array,& L1, arrayo, &inc, myStreams(f1, dsyBuffers,rank ), &lbuffer, &info);
    return info;
}


INT_TYPE tdpotrs ( INT_TYPE L1, INT_TYPE M2, double * array, double * arrayo ,INT_TYPE inc){
    INT_TYPE info;
    char charU = 'U';
#if VERBOSE
    printf("pot\n");
    fflush(stdout);
#endif

#ifdef APPLE
dpotrs_(&charU,&L1,  &M2, array, &L1, arrayo, &inc, &info );
#else
    INT_TYPE i;
    for ( i = 0; i < L1*M2 ; i++){
        if ( isnan(arrayo[i]) || isinf(arrayo[i])){
            printf("warning inf, try reducing the rank of vectors\n");
            fflush(stdout);
            exit(0);
        }
    }
    info =  LAPACKE_dpotrs(LAPACK_COL_MAJOR,'U',L1,  M2, array, L1, arrayo, inc );
#endif
#if VERBOSE
    printf("pot-end\n");
    fflush(stdout);
#endif
    return info;

}

double tdpocon (INT_TYPE rank,struct sinc_label  f1,  INT_TYPE L1 , double * Matrix ){
    INT_TYPE info;
    char charU = 'U';
    double rcond;
  // return 1.0;
#if VERBOSE
    printf("pop\n");
    fflush(stdout);
#endif
    double norm1 = cblas_dasum(L1*L1, Matrix, 1);
    cblas_dcopy(L1*L1,Matrix,1,myStreams(f1, dsyBuffers,rank ),1);
    #ifdef APPLE
    INT_TYPE lbuffer = part(f1, dsyBuffers)-L1*L1;
    dpocon_( &charU, &L1, myStreams(f1, dsyBuffers,rank ), &L1, &norm1, &rcond, myStreams(f1, dsyBuffers,rank )+L1*L1,&lbuffer ,&info );
#else
    INT_TYPE lbuffer = part(f1, dsyBuffers)-L1*L1;
    info =  LAPACKE_dpocon(LAPACK_COL_MAJOR,charU,L1,  myStreams(f1, dsyBuffers,rank ), L1, norm1, &rcond);
#endif
//    printf("--%f\n", rcond);

    if ( info != 0 ){
        printf ("pop corn ! %d\n",info);
        exit(0);
    }
#if VERBOSE
    printf("end pop %1.15f\n",rcond);
    fflush(stdout);
#endif
    return rcond;
}


//tdgesvd(0, f1, N1*N1,N2*N2, streams(f1, quadCube,0,space), streams(f1, interactionExchangePlus,0,space)+CanonicalRank(f1, interactionExchangePlus, 0)*N1*N1, streams(f1, interactionExchangePlus,0,space2)+CanonicalRank(f1, interactionExchangePlus, 0)*N2*N2)
#if 1
INT_TYPE tdgesvd ( INT_TYPE rank, struct sinc_label f1 ,  INT_TYPE M1, INT_TYPE M2, Stream_Type * ge, Stream_Type * m1, Stream_Type* m2 ){
    INT_TYPE info,Mm =imin(M1,M2),i;
    char Job = 'S';
    INT_TYPE lbuffer = part(f1, dsyBuffers)-M1-M2;
    Stream_Type * sg = myStreams(f1, dsyBuffers,rank );
    Stream_Type * buf = myStreams(f1, canonicalBuffersBM, rank);
    Stream_Type * work = myStreams(f1, dsyBuffers,rank )+M1+M2;
    
//    {
//        double sum = 0.;
//        sum += myStreams(f1, oneByOneBuffer,0)[0];
//        printf("%f ",sum);
//    }
    
    
    dgesvd_(&Job, &Job, &M1, &M2, ge, &M1, sg, m1, &M1, buf, &Mm, work, &lbuffer, &info);
    if (info ){
        printf("svd %d\n",info);
        lbuffer = -1;
        
      //  dgesvd_(&Job, &Job, &M1, &M2, ge, &M1, sg, m1, &M1, buf, &Mm, work, &lbuffer, &info);
        printf("%f %d\n", work[0],part(f1, dsyBuffers)-M1-M2);
//
//        Job = 'N';
//        dgesvd_(&Job, &Job, &M1, &M2, ge, &M1, sg, m1, &M1, buf, &Mm, work, &lbuffer, &info);
        for ( i = 1; i < Mm+1 ; i++)
            printf("%f\t", work[i]);
        exit(9);
    }
    
//    {
//        INT_TYPE ii;
////    for ( ii = 0 ; ii < min(M1,M2); ii++)
//        printf("%d-%f\n", ii, sg[ii]);
//    }
    
    transpose(Mm, M2, buf, m2);
    
    
    for ( i = 0; i < Mm; i++){
          cblas_dscal(M1, sg[i], m1+i*M1, 1);

//             cblas_dscal(M1, (sg[i]), m1+i*M1, 1);
    }
    
//    {
//        INT_TYPE l ;
//        double sum = 0;
//        for ( l = 0; l < Mm ; l++)
//            sum += m1[l*M1]*m2[l*M2];
//
//        printf("%f \n",sum);
//
//    }
    return Mm;
}
#else
INT_TYPE tdgesvd ( INT_TYPE rank, struct sinc_label f1 ,  INT_TYPE M1, INT_TYPE M2, Stream_Type * ge, Stream_Type * m1, Stream_Type* m2 ){
    INT_TYPE info,Mm =imin(M1,M2),i;
    char Job = 'S';
    INT_TYPE lbuffer = part(f1, dsyBuffers)-M1-M2;
    Stream_Type * sg = myStreams(f1, dsyBuffers,rank );
    Stream_Type * buf = myStreams(f1, canonicalBuffersBM, rank);
    Stream_Type * work = myStreams(f1, dsyBuffers,rank )+M1+M2;
    // printf("%lu %lu %lu \n",sg-ge, buf-sg,work-buf);
    dgesvd_(&Job, &Job, &M1, &M2, ge, &M1, sg, m1, &M1, buf, &Mm, work, &lbuffer, &info);
    // printf("alloc %d : %d %d %d\n", Mm*M2,alloc(f1, interaction1Plus, 0), part(f1,interaction1Plus),alloc(f1, interaction1Plus, 0)* part(f1,interaction1Plus));
    if (info ){
        printf("svd %d\n",info);
        exit(9);
    }
    INT_TYPE Mt = 0;
    // cblas_dcopy(M2, myStreams(f1, canonicalBuffersBM, rank), 1, m2, 1);
    transpose(Mm, M2, myStreams(f1, canonicalBuffersBM, rank), m2);
    for ( i = 0; i < Mm; i++){
        //   printf("%d %1.15f\n", i, sg[i] );
        cblas_dscal(M1, sg[i], m1+i, M1);
    }
    
    
    return Mm;
}

#endif



INT_TYPE tInverse( struct sinc_label f1, INT_TYPE n, double * ar){
    INT_TYPE info=0,info2=0;
    
#ifndef APPLE
    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,ar,n,(INT_TYPE*)myStreams(f1,dsyBuffers,0));
    info2 = LAPACKE_dgetri(LAPACK_COL_MAJOR,n,ar,n,(INT_TYPE*)myStreams(f1,dsyBuffers,0));
#endif
    return info+1000*info2;
}

INT_TYPE tzInverse( struct sinc_label f1, INT_TYPE n, DCOMPLEX * ar){
    INT_TYPE info=0,info2=0;
    
#ifndef APPLE
    info = LAPACKE_zgetrf(LAPACK_COL_MAJOR,n,n,ar,n,(INT_TYPE*)myStreams(f1,dsyBuffers,0));
    info2 = LAPACKE_zgetri(LAPACK_COL_MAJOR,n,ar,n,(INT_TYPE*)myStreams(f1,dsyBuffers,0));
#endif
    return info+1000*info2;
}
