/**
 *  interfaceMath.c
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

double pSinc( double d , double x, inta N1){
    inta i,N12 = (N1-1)/2;
    double su=0.;
    for ( i = -N12; i <= N12 ; i++)
        su += cos( 2.*pi*i * x/d /N1);
    return 1./(N1)*su;
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

double pSS( double d1 , double x , inta N1, double d2, double y, inta N2 )    {
    if ( d1 > d2 ){
        if ( x == y )
            return sqrt(d2/d1);
        else
            return pSinc(d1,x-y,N1) * sqrt(d2/d1);
    }
    if ( x == y )
        return sqrt(d1/d2);
    else
        return pSinc(d2,x-y,N2) * sqrt(d1/d2);
    return 0;
}



void transpose(inta N, inta M, floata * orig, floata* targ){

#ifndef MKL
    inta i ,j;
    for ( i = 0 ; i < N ; i++)
        for ( j = 0; j < M ; j++){
            targ[i*M+j] = orig[j*N+i];
        }
#else
    mkl_domatcopy ('C', 'T',N, M,1.,orig,N,targ,M);
#endif

}

inta tdsyev( inta rank, char job , inta n, double * ar, inta ns , double * w ){
    
    if ( n == 0 ||ns == 0 )
    {
        printf ( "tdsyev, either there are no vectors or no stride\n");
        exit(0);
    }


#ifdef APPLE

    char charU = 'U';
    inta lwork = 5 * n;
    floata work [lwork];
    
    inta info;
    dsyev_ ( &job, &charU, &n , ar , &ns , w , &work[0], &lwork, &info );
    return info;
#else
    return LAPACKE_dsyev(LAPACK_COL_MAJOR, job, 'U', n, ar, ns , w);
#endif
}

inta tzheev( inta rank,   sinc_label f1, char job , inta n, DCOMPLEX * ar, inta ns , double * w ){

#ifdef APPLE
    char charU = 'U';
    inta info,lbuffer = part(f1, dsyBuffers)-3*n;
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
        printf ( "tzheev, empty \n");
        exit(0);
    }
    return LAPACKE_zheev(LAPACK_COL_MAJOR, job, 'U', n, (DCOMPLEX_PRIME *)ar, ns , w );
#endif
}


inta tdsygv( inta rank,   sinc_label f1, char job , inta n, double * sr, double * ar, inta ns , double * w ){
    inta info;
    inta type = 1;
    char charU = 'U';
#ifdef APPLE
    inta lbuffer = n*5;
    double buffer [lbuffer];
    dsygv_(&type, &job, &charU,&n,sr,&ns,ar,&ns, w,&buffer[0],&lbuffer, &info );
#else
    info =  LAPACKE_dsygv(LAPACK_COL_MAJOR,1, job, 'U', n,sr,ns, ar, ns , w );
#endif
    return info;

}

                     
inta tzhegv( inta rank,   sinc_label f1, char job , inta n,DCOMPLEX * sr, DCOMPLEX * ar, inta ns , double * w ){
    inta info;
//    inta liwork = 3+5*n;
//    inta lrwork = 1+5*n+2*n*n;
//    inta lwork = 2*n+n*n;
//    inta iwork[liwork];
//    double rwork[lrwork];
//    DCOMPLEX work[lwork];
    inta type = 1;
    char charU = 'U';
#ifdef APPLE
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
    inta lbuffer = (part(f1,dsyBuffers)-10*n)/2;
    zhegv_(&type, &job, &charU,&n,(__CLPK_doublecomplex *)sr,&ns,(__CLPK_doublecomplex *)ar,&ns, w,(__CLPK_doublecomplex *) (myStreams(f1, dsyBuffers,rank ))+10*n, &lbuffer,  myStreams(f1, dsyBuffers,rank ),&info );
//    } else {
//        zhegvd_(&type, &job , &charU, &n,(__CLPK_doublecomplex *)sr,&ns,(__CLPK_doublecomplex *)ar,&ns, w, (__CLPK_doublecomplex *)(work), &lwork, rwork, &lrwork, iwork, &liwork, &info);
//    }
#else
    info =  LAPACKE_zhegv(LAPACK_COL_MAJOR,1, job, 'U', n,(DCOMPLEX_PRIME *)sr,ns,(DCOMPLEX_PRIME *) ar, ns , w );
#endif
    
    return info;
}
                               

                              
inta tdgeqr( inta rank,   sinc_label f1,inta len, inta n, double * ar, inta ns ,double *w, double *xr , inta xs ){
    inta info,info2,i;
    if ( len > n ){
        printf("Gram Schmidt: %d -> %d\n", len, n);
        len = n;
    }else if ( len <= 0 ){
        printf("cancel tdgeqr\n");
        return -1;
    }
    ///len is number of elements to consider, if len > n, then you may as well remove some of len...its over-linear dependent.

#ifdef APPLE
    inta lbuffer = part(f1, dsyBuffers);
    dgeqrf_(&n,& len, ar, &ns, w, myStreams(f1, dsyBuffers,rank ), &lbuffer, &info);
    for ( i = 0 ; i < len ; i++){
        cblas_dcopy(i+1, ar+ns*i, 1, xr+xs*i, 1);
    }
    dorgqr_(&n, &len, &len, ar, &ns, w, myStreams(f1, dsyBuffers,rank ), &lbuffer, &info2);
#else
    info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, n,len, ar, ns, w);
    for ( i = 0 ; i < len ; i++){
        cblas_dcopy(i+1, ar+ns*i, 1, xr+xs*i, 1);
    }
    info2 = LAPACKE_dorgqr(LAPACK_COL_MAJOR,n,len,len,ar,ns,w);
    
#if 0
    ///mapping back...
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, len, len, 1., ar, ns, xr, xs, 0., out, ns);
#endif
#endif
    return info * 100000 + info2;

}


inta tdpotrf ( inta L1, double * array ,inta LS1) {
    inta info;
    char charU = 'U';
    
#ifdef APPLE
    dpotrf_( &charU , &L1, array, &LS1, &info ) ;
#else
    info =  LAPACKE_dpotrf(LAPACK_COL_MAJOR,'U',L1,  array, LS1 );
#endif
    return info;
}

inta tdgels ( inta rank,  sinc_label  f1 , inta L1, inta M2, double * array, double * arrayo ,inta inc){
    char charT = 'N';
    inta lbuffer = part(f1, dsyBuffers),info=0;
    dgels_(&charT, &L1, &L1,& M2, array,& L1, arrayo, &inc, myStreams(f1, dsyBuffers,rank ), &lbuffer, &info);
    return info;
}


inta tdpotrs ( inta L1, inta M2, double * array,inta LS1, double * arrayo ,inta inc){
    inta info;
    char charU = 'U';

#ifdef APPLE
    dpotrs_(&charU,&L1,  &M2, array, &LS1, arrayo, &LS1, &info );
#else
    info =  LAPACKE_dpotrs(LAPACK_COL_MAJOR,'U',L1,  M2, array, LS1, arrayo, LS1 );
#endif
    return info;

}

double tdpocon (inta rank,  sinc_label  f1,  inta L1 , double * Matrix , inta Stride){
    inta info;
    char charU = 'U';
    double rcond;
    inta l;
    ///not used, but may want to use stride,
    double norm1=0.;
    for ( l = 0 ; l < L1 ; l++)
        norm1 += cblas_dasum(L1, Matrix+l*Stride, 1);
    #ifdef APPLE
    inta lbuffer = part(f1, dsyBuffers);
    dpocon_( &charU, &L1, Matrix, &L1, &norm1, &rcond, myStreams(f1, dsyBuffers,rank ),&lbuffer ,&info );
#else
    info =  LAPACKE_dpocon(LAPACK_COL_MAJOR,charU,L1,  Matrix, Stride, norm1, &rcond);
#endif

    if ( info != 0 ){
        printf ("pop corn ! %d\n",info);
        exit(0);
    }
    return rcond;
}

inta tInverse(   sinc_label f1, inta n, double * ar){
    inta info=0,info2=0;
    
#ifndef APPLE
    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,ar,n,(inta*)myStreams(f1,dsyBuffers,0));
    info2 = LAPACKE_dgetri(LAPACK_COL_MAJOR,n,ar,n,(inta*)myStreams(f1,dsyBuffers,0));
#else
    inta z[n];
    inta lwork = part(f1,dsyBuffers);
    dgetrf_(&n, &n, ar, &n, z, &info);
    dgetri_(&n, ar, &n, z, (floata*)myStreams(f1,dsyBuffers,0), &lwork, &info2);
#endif
    
    
    return info+1000*info2;
}

inta tzInverse(   sinc_label f1, inta n, DCOMPLEX * ar){
    inta info=0,info2=0;
    
#ifndef APPLE
    info = LAPACKE_zgetrf(LAPACK_COL_MAJOR,n,n,(DCOMPLEX_PRIME *)ar,n,(inta*)myStreams(f1,dsyBuffers,0));
    info2 = LAPACKE_zgetri(LAPACK_COL_MAJOR,n,(DCOMPLEX_PRIME *)ar,n,(inta*)myStreams(f1,dsyBuffers,0));
#endif
    return info+1000*info2;
}

inta tLowdin( inta n , double *ar, double *lowdinVec, double * lowdinMatrix ){
    inta i;
    double w[n];
    tdsyev(0, 'V', n, ar, n, &w[0]);
    if ( w[0] <= 0. ){
        printf("non invertable overlap\n");
        exit(0);
    }
    for ( i = 0 ; i < n*n; i++){
        lowdinVec[i] = 0.;
        lowdinMatrix[i] = 0.;
    }
    for ( i = 0 ; i < n ; i++){
        cblas_dger(CblasColMajor, n,n, sqrt(w[i]) , ar+i*n,1, ar+i*n,1, lowdinVec,n);
        cblas_dger(CblasColMajor, n,n, 1./sqrt(w[i]) , ar+i*n,1, ar+i*n,1, lowdinMatrix,n);
    }
    return 0;
}
