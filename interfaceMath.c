/*
 *  interfaceMath.c
 *
 *
 *  Copyright 2018 Jonathan Jerke and Bill Poirier.
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

double periodicSinc ( double d , double x, INT_TYPE N1 ){
    double dimless =1., arg = pi*(x/d);
    
    
    if ( arg != 0.0 )
        dimless = sin( arg )/2./N1 *(1./tan(arg/2./N1)+tan(arg/2./N1));
    if (  isnan(dimless) ||  isinf(dimless))
        return 1.;
    
    
    return dimless;;
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

double periodicSS( double d1 , double x ,INT_TYPE N1, double d2, double y , INT_TYPE N2)    {
    if ( d1 > d2 ){
        if ( x == y )
            return sqrt(d2/d1);
        else
            return periodicSinc(d1,x-y,N1) * sqrt(d2/d1);
    }
    if ( x == y )
        return sqrt(d1/d2);
    else
        return periodicSinc(d2,x-y,N2) * sqrt(d1/d2);
    return 0;
}


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

INT_TYPE tdsyev( INT_TYPE rank, struct field * f1, char job , INT_TYPE n, double * ar, INT_TYPE ns , double * w ){
#ifdef APPLE

    char charU = 'U';
    INT_TYPE info,lbuffer = part(f1, dsyBuffers);
    dsyev_ ( &job, &charU, &n , ar , &ns , w , myStreams(f1, dsyBuffers,rank ), &lbuffer, &info );
    return info;
#else
    return LAPACKE_dsyev(LAPACK_COL_MAJOR, job, 'U', n, ar, ns , w);
#endif
}

INT_TYPE tzheev( INT_TYPE rank, struct field * f1, char job , INT_TYPE n, DCOMPLEX * ar, INT_TYPE ns , double * w ){
#ifdef APPLE
    
    char charU = 'U';
    INT_TYPE info,lbuffer = part(f1, dsyBuffers);
  zheev_ (&job, &charU, &n , (__CLPK_doublecomplex *)ar , &ns , w ,(__CLPK_doublecomplex *) (myStreams(f1, dsyBuffers,rank )+3*n), &lbuffer,  myStreams(f1, dsyBuffers,rank ),&info );
    return info;
#else
    return LAPACKE_zheev(LAPACK_COL_MAJOR, job, 'U', n, (DCOMPLEX_PRIME *)ar, ns , w );
#endif
}


INT_TYPE tdsygv( INT_TYPE rank, struct field * f1, char job , INT_TYPE n, double * sr, double * ar, INT_TYPE ns , double * w ){
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

                     
INT_TYPE tzhegv( INT_TYPE rank, struct field * f1, char job , INT_TYPE n,DCOMPLEX * sr, DCOMPLEX * ar, INT_TYPE ns , double * w ){
    INT_TYPE info;
    INT_TYPE type = 1;
    char charU = 'U';
#ifdef APPLE
    INT_TYPE lbuffer = part(f1, dsyBuffers);
    zhegv_(&type, &job, &charU,&n,(__CLPK_doublecomplex *)sr,&ns,(__CLPK_doublecomplex *)ar,&ns, w,(__CLPK_doublecomplex *) (myStreams(f1, dsyBuffers,rank )+6*n), &lbuffer,  myStreams(f1, dsyBuffers,rank ),&info );
#else
    info =  LAPACKE_zhegv(LAPACK_COL_MAJOR,1, job, 'U', n,(DCOMPLEX_PRIME *)sr,ns,(DCOMPLEX_PRIME *) ar, ns , w );
#endif
    
    return info;
}
                               

                              
INT_TYPE tdgeqr( INT_TYPE rank, struct field * f1,INT_TYPE len, INT_TYPE n, double * ar, INT_TYPE ns ,double *w){
    INT_TYPE info,info2;
    if ( len > n ){
        printf("%d -> %d\n", len, n);
        len = n;
    }
    //len is number of elements to consider, if len > n, then you may as well remove some of len...its over-linear dependent.
#if VERBOSE
    printf("eqr\n");
    fflush(stdout);
#endif


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
    
    return info;
}



INT_TYPE tdpotrs ( INT_TYPE L1, INT_TYPE M2, double * array, double * arrayo ){
    INT_TYPE info;
    char charU = 'U';
#if VERBOSE
    printf("pot\n");
    fflush(stdout);
#endif

#ifdef APPLE
dpotrs_(&charU,&L1,  &M2, array, &L1, arrayo, &L1, &info );
#else
    INT_TYPE i;
    for ( i = 0; i < L1*M2 ; i++){
        if ( isnan(arrayo[i]) || isinf(arrayo[i])){
            printf("warning inf, try reducing the rank of vectors\n");
            fflush(stdout);
            exit(0);
        }
    }
    info =  LAPACKE_dpotrs(LAPACK_COL_MAJOR,'U',L1,  M2, array, L1, arrayo, L1 );
#endif
#if VERBOSE
    printf("pot-end\n");
    fflush(stdout);
#endif
    return info;

}

double tdpocon (INT_TYPE rank,struct field * f1,  INT_TYPE L1 , double * Matrix ){
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
    if ( info != 0 ){
        printf ("pop corn ! %lld\n",info);
        exit(0);
    }
#if VERBOSE
    printf("end pop %1.15f\n",rcond);
    fflush(stdout);
#endif
    return rcond;
}


