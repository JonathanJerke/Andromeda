/**
 *  interfaceMath.c
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

//DCOMPLEX periodicSinc ( double d , double x, double momentum, inta N1 ){
//    double dimless =1., arg = pi*(x/d);
//    
//    
//    if ( arg != 0.0 )
//        dimless = sin( arg )/2./N1 *(1./tan(arg/2./N1)+tan(arg/2./N1));
////    if (  isnan(dimless) ||  isinf(dimless))
////        return 1.;
//    
//    
//    return dimless * ei( momentum * d * x );
//}

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

inta tdsyev( inta rank,   sinc_label f1, char job , inta n, double * ar, inta ns , double * w ){
    
    if ( n == 0 ||ns == 0 )
    {
        printf ( "tdsyev, either there are no vectors or no stride\n");
        exit(0);
    }


#ifdef APPLE

    char charU = 'U';
    inta info,lbuffer = part(f1, dsyBuffers);
    dsyev_ ( &job, &charU, &n , ar , &ns , w , myStreams(f1, dsyBuffers,rank ), &lbuffer, &info );
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
    inta lbuffer = part(f1, dsyBuffers);
    dsygv_(&type, &job, &charU,&n,sr,&ns,ar,&ns, w,myStreams(f1, dsyBuffers,rank ),&lbuffer, &info );
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
                               

                              
inta tdgeqr( inta rank,   sinc_label f1,inta len, inta n, double * ar, inta ns ,double *w){
    inta info,info2;
    if ( len > n ){
        printf("%d -> %d\n", len, n);
        len = n;
    }else if ( len <= 0 ){
        printf("cancel tdgeqr\n");
        return -1;
    }
    ///len is number of elements to consider, if len > n, then you may as well remove some of len...its over-linear dependent.

    printf("tdgeqr %d %d %d %d\n", rank, len, n, ns );
#ifdef APPLE
    inta lbuffer = part(f1, dsyBuffers);
    dgeqrf_(&n,& len, ar, &ns, w, myStreams(f1, dsyBuffers,rank ), &lbuffer, &info);
    dorgqr_(&n, &len, &len, ar, &ns, w, myStreams(f1, dsyBuffers,rank ), &lbuffer, &info2);
#else
    info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, n,len, ar, ns, w);
    info2 = LAPACKE_dorgqr(LAPACK_COL_MAJOR,n,len,len,ar,ns,w);
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
    
    ///not used, but may want to use stride,
    double norm1 = cblas_dasum(L1*L1, Matrix, 1);
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
