/*
 *  mAls.c
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

#include "coreUtil.h"
#include "mAls.h"

INT_TYPE normalize (struct field * f1,  enum division alloy, INT_TYPE spin, INT_TYPE space){
    INT_TYPE l,flag=0;
    double norm;
    INT_TYPE M2[3];
    length(f1,alloy,M2);
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
    INT_TYPE M2[3];
    length(f1,alloy,M2);
    double *Alloy =streams(f1,alloy,spin,space);
    double *Origin =streams(f1,origin,os,space);
    {
        
        for ( m = 0; m < L1; m++)
            for ( n = 0; n <=m ; n++){
                
                output[ m*L1 + n ] = cblas_ddot(M2[space], Alloy+m*M2[space] ,1,Alloy+n*M2[space] ,1);
                output[ n*L1 + m ] = output[ m*L1 + n ];
            }
        for ( m = 0; m < L1; m++)
            for ( n = 0; n < G1 ; n++)
                output2[ n*L1 + m ] = cblas_ddot(M2[space], Alloy+m*M2[space] ,1,Origin+n*M2[space] ,1);
    }
    //printf("%f__\n",output[0]);
    return 0;
}



INT_TYPE balance (struct field * f1,  enum division alloy, INT_TYPE spin){
    INT_TYPE l;
    INT_TYPE L1 = CanonicalRank(f1, alloy,spin),sign[3],signs;
    double snorm,value;
    long double factor,norm[3] ;
    INT_TYPE M2[3],space;
    length(f1,alloy,M2);
    INT_TYPE iOne = 1;
    //#pragma omp parallel for private(l,norm,flag,ii,i,trace,space)
    
    
    if ( SPACE == 3 ){
        for ( l = 0; l < L1 ;l++){
            
            for ( space = 0; space < SPACE ; space++){
                norm[space] = cblas_dnrm2(M2[space], streams(f1, alloy,spin,space)+l*M2[space],iOne);
                value = cblas_idamax(M2[space], streams(f1, alloy,spin,space)+l*M2[space], 1);
                if ( value == 0. )
                    sign[space] = 0;
                if ( value > 0. )
                    sign[space] = 1;
                else
                    sign[space] = -1;
            }
            
            
            if ( sign[0] * sign[1] * sign[2] == 0 ){
                // very unlikely...
                //if so , all set to zero...
            }
            
            factor = pow( norm[0]*norm[1]*norm[2],0.333333333333333333333333333333333333333333333333333);
            
            for ( space = 0; space < SPACE ; space++){
                
                if ( space == 0 )
                    signs = sign[1] * sign[2];
                else if ( space == 1 )
                    signs = sign[0] * sign[2];
                else
                    signs = sign[0] * sign[1];
                snorm = signs*factor/norm[space] ;
                cblas_dscal(M2[space], snorm, streams(f1, alloy,spin,space)+l*M2[space],iOne);
            }
            
        }
    }else if ( SPACE == 2){
        for ( l = 0; l < L1 ;l++){
            
            for ( space = 0; space < SPACE ; space++){
                norm[space] = cblas_dnrm2(M2[space], streams(f1, alloy,spin,space)+l*M2[space],iOne);
                value = cblas_idamax(M2[space], streams(f1, alloy,spin,space)+l*M2[space], 1);
                if ( value == 0. )
                    sign[space] = 0;
                if ( value > 0. )
                    sign[space] = 1;
                else
                    sign[space] = -1;
            }
            
            
            if ( sign[0] * sign[1]  == 0 ){
                // very unlikely...
                //if so , all set to zero...
            }
            
            factor = pow( norm[0]*norm[1],1./SPACE);
            
            for ( space = 0; space < SPACE ; space++){
                
                if ( space == 0 )
                    signs = sign[1] ;
                else
                    signs = 1;
                snorm = signs*factor/norm[space] ;
                cblas_dscal(M2[space], snorm, streams(f1, alloy,spin,space)+l*M2[space],iOne);
            }
            
        }
    }
    return 0;
}

double canonicalDecompositionMP( INT_TYPE rank,struct field * f1 , enum division origin,INT_TYPE os,   enum division alloy ,  INT_TYPE spin ,double tolerance){
    if ( tolerance < 0 ){
        //printf("%d->%d\n", origin, alloy);
        tolerance = -(tolerance);
        
        //   exit(0);
    }
    
    INT_TYPE ns = 9;
    Stream_Type * array[3];
    Stream_Type * array2[3];
    INT_TYPE space,space0,space1,space2;
    //printf("%lld %lld %lld\n", spin,  f1->sinc.tulip[origin].Current[spin], f1->sinc.tulip[alloy].Current[spin]);
    //fflush(stdout);
    INT_TYPE l,count = 1 ;
    INT_TYPE L1 = CanonicalRank(f1,alloy,spin);;
    INT_TYPE G1 =CanonicalRank(f1,origin,os);
    if ( G1 == 0 ){
        f1->sinc.tulip[alloy].Current[spin] = 0;
#if VERBOSE
        printf("CD: empty origin %d _%lld -> %d _%lld\n", origin, os , alloy, spin);
#endif
        return 0;
    }
    
    if ( L1 == 0 ){
     //   printf("CD: zero length %lld %lld %lld\n", origin,alloy,spin);
        tId(f1, alloy, spin);
        L1 = 1;
        if (! CanonicalRank(f1, alloy, spin)){
            printf("uhm... %d\n",alloy);
            exit(0);
        }
    }
    double value2=0.,ALPHA = f1->mem1->rt->ALPHA;
    INT_TYPE M2[3];
    length(f1,alloy,M2);
    INT_TYPE info;
    double rcond;
    enum division canonicalStore;
    INT_TYPE T1 = G1*L1;
    
    {//REF CORE NUMBER
        //GET CORE... ALLOCATE TO CORE-LANE.
        
        for ( space = 0; space < SPACE ; space++){
            array[space] =  streams(f1, canonicalBuffers, rank , space);
            array2[space] =  streams(f1, canonicalBuffers, rank , space) + L1*L1;
        }
        
        if (  T1 + L1*L1 >  part(f1,canonicalBuffers)){
#if 1
            printf("mem prob\n %lld+ %lld  > %lld \n",T1,L1*L1,part(f1,canonicalBuffers) );
#endif
            exit(0);
        }
        
        if( (species(f1, alloy ) == vector && bodies(f1, alloy ) == bodies(f1,canonicalBuffersB )) ||(species(f1, alloy ) == matrix && bodies(f1, alloy ) == one && f1->body > one)){
            canonicalStore = canonicalBuffersB;
        }else if( species(f1, alloy ) == vector && bodies(f1, alloy ) == bodies(f1,canonicalBuffersBX ) ){
            canonicalStore = canonicalBuffersBX;
        }
        else {
            canonicalStore = canonicalBuffersBM;
        }
        
        
        if ( part(f1, canonicalStore) < L1 || alloc(f1, canonicalStore) < alloc(f1, alloy) )
        {
#if 1
            printf("species-%d\nbodies%d\n", species(f1, canonicalStore), bodies(f1, canonicalStore));
            printf("%d -> %d (%lld) <= %d \n", origin, alloy,L1,part(f1, canonicalStore));
            printf("memory PROB\n");
#endif
            exit(0);
        }
        
        
    }
    
    
    char side = CDT;
    
    
    enum division alloyBak;
    if ( species (f1, alloy ) == vector && bodies(f1,alloy) == one){
        alloyBak = trainVector;
    }
    else     if ( species (f1, alloy ) == vector && bodies(f1,alloy) == two){
        alloyBak = trainVector2;
    }
    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == one){
        
        f1->sinc.tulip[trainMatrix].header = header(f1,alloy);
        
        alloyBak = trainMatrix;
    }else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == two){
        
        f1->sinc.tulip[trainMatrix2].header = header(f1,alloy);
        
        alloyBak = trainMatrix2;
    }
    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == three){
        
        f1->sinc.tulip[trainVector3].header = header(f1,alloy);
        
        alloyBak = trainVector3;
    }
    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == four){
        
        f1->sinc.tulip[trainVector4].header = header(f1,alloy);
        
        alloyBak = trainVector4;
    }
    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == three){
        
        f1->sinc.tulip[trainMatrix3].header = header(f1,alloy);
        
        alloyBak = trainMatrix3;
    }
    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == four){
        
        f1->sinc.tulip[trainMatrix4].header = header(f1,alloy);
        
        alloyBak = trainMatrix4;
    }
    else if ( species(f1, alloy ) == quartic ){
        f1->sinc.tulip[trainQuartic].header = header(f1,alloy);
        
        alloyBak = trainQuartic;
    } else{
        printf("cd: tObject species %d \n", alloy);
        
        exit(0);
        
    }
    if (  ( purpose(f1, alloy ) < tObject   )){
        printf("cd: tObject status %d \n", alloy);
        exit(0);
        
    }
   // printf(":: %d %d -- %d %d \n", origin , os, alloy ,spin );

    
    if ( spread(f1,origin,os,alloy,spin,1,array[1],array2[1]) )
        return -1;
    
    do{
        if ( count % ns == 0 )
            tEqua(f1, alloyBak,rank, alloy, spin );
        
        for ( space0 = 0; space0 < SPACE ;space0++){
            
            if ( SPACE == 3 ){
                if ( space0 == 0 ){
                    space0 = 0;
                    space1 = 1;
                    space2 = 2;
                } else if ( space0 == 1 ){
                    space0 = 1;
                    space1 = 2;
                    space2 = 0;
                }else if ( space0 == 2 ){
                    space0 = 2;
                    space1 = 0;
                    space2 = 1;
                } else {
                    printf("cd space");
                    exit(0);
                }
            }
            else if ( SPACE == 2 ){
                if ( space0 == 0 ){
                    space0 = 0;
                    space1 = 1;
                } else if ( space0 == 1 ){
                    space0 = 1;
                    space1 = 0;
                }else {
                    printf("cd space");
                    exit(0);
                }
            }
            
            
//
//            if ( normalize(f1,alloy,spin,space1) )
//                return -1;
            
            
            
            if ( SPACE == 3 ){
                //normal 3d calculator
                
                if ( spread(f1,origin,os,alloy,spin,space1,array[space1],array2[space1]) )
                    return -1;

                if ( spread(f1,origin,os,alloy,spin,space2,array[space2],array2[space2]) )
                    return -1;
                
                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1*L1, 0,array[space2],1, array[space1],1 );
                //replace Matrix
                
                if ( L1 >1 )
                    for ( l = 0 ; l < L1; l++){
                        if ( fabs(array[space1][ l*L1 + l ]-1.0) > 1e-6 ||  isnan ( array[space1][l*L1+l] )  ||  isinf ( array[space1][l*L1+l] ) ){
#if VERBOSE
                            printf("%lld:%f (%lld)-\n", l,array[space1][ l*L1 + l ],rank);
                            fflush(stdout);
#endif
                        }
                        array[space1][ l*L1 + l ] += ALPHA;
                    }
                
                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,T1   , 0,array2[space2],1, array2[space1],1 );
                
                //replace Vectors
            }else if (SPACE == 2){
                //plane-calculator ..  third dimenison is absent, allowing for simplication.
                
                
                if ( spread(f1,origin,os,alloy,spin,space1,array[space1],array2[space1]) )
                    return -1;

                if ( L1 >1 )
                    for ( l = 0 ; l < L1; l++){
#if VERBOSE

                        if ( fabs(array[space1][ l*L1 + l ]-1.0) > 1e-6 ||  isnan ( array[space1][l*L1+l] )  ||  isinf ( array[space1][l*L1+l] ) ){
                            printf("%lld:%f (%lld)-\n", l,array[space1][ l*L1 + l ],rank);
                            fflush(stdout);
                        }
#endif

                        array[space1][ l*L1 + l ] += ALPHA;
                    }

            }

            info = tdpotrf(L1, array[space1]);

            if ( info != 0 ){
#if VERBOSE
                printf("info failure %lld %lld %d %lld %d %lld\n", info, rank,origin,os, alloy, spin);
#endif
                return 1;
            }
            rcond = tdpocon(rank, f1, L1, array[space1] );
//            if ( f1->mem1->rt->targetCondition > 0 )
//                ALPHA /= min(1e1,max(1e-1, 1 - log(rcond*f1->mem1->rt->targetCondition) ));

            if (  isnan(rcond) || isinf(rcond) || rcond < 1e-12){
#if VERBOSE
                printf("Warning: condition of Beylkin %d->%d is bad %16.16f\n", origin, alloy, rcond);
                fflush(stdout);
#endif
                return 1;
            }
            
            if(1){
                INT_TYPE i;
                
                if(VERBOSE){
                    for ( i = 0; i < G1*M2[space0] ; i++)
                        if ( isnan(streams(f1,origin,os,space0)[i])||isinf(streams(f1,origin,os,space0)[i]))
                        {
                            printf("inf\n");
                            fflush(stdout);
                            exit(0);
                        }
                    for ( i = 0; i < T1 ; i++)
                        if ( isnan( array2[space1][i])||isinf( array2[space1][i]))
                        {
                            printf("inf array\n");
                            fflush(stdout);
                            
                            exit(0);
                        }
                }
                
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,L1,M2[space0],G1,1,array2[space1],L1,streams(f1,origin,os,space0),M2[space0], 0., myStreams(f1, canonicalStore, rank), L1 );
                if(VERBOSE){
                    for ( i = 0; i < L1*M2[space0] ; i++)
                        if ( isnan( myStreams(f1, canonicalStore, rank)[i])||isinf( myStreams(f1, canonicalStore, rank)[i]))
                        {
                            printf("inf can\n");
                            fflush(stdout);
                            
                            exit(0);
                        }
                }
            }
            
            info = tdpotrs(L1,  M2[space0], array[space1],  myStreams(f1, canonicalStore, rank) );
            if ( info != 0 ){
#if VERBOSE
                printf("L1 %lld \n", L1);
                printf("M2 %lld \n", M2[0]);
                printf("R1 %lld \n", CanonicalRank(f1, origin, os));
                printf("info failture %lld\n", info);
#endif
                return 1;
            }
            
            
            
            
            transpose ( L1, M2[space0] , myStreams(f1, canonicalStore, rank) , streams(f1,alloy,spin,space0));
            
            
            
        }
        balance(f1, alloy , spin );

        if ( count % ns == 0 ){
            ns++;
            value2 = distanceOne(rank, f1, alloy, spin,alloyBak, rank);
        }
        count++;

    } while (  value2 >tolerance && ns < 16) ;
    return 0;
}

double tCycleDecompostionOneMP ( INT_TYPE rank, struct field * f1 , enum division origin, INT_TYPE os, enum division alloy,INT_TYPE spin,  double tolerance , INT_TYPE maxRun , double power  ){
    double value , value2 = sqrt(inner(rank, f1, origin, os));
    if ( maxRun >= CanonicalRank(f1, origin, os)){
        tEqua(f1, alloy,spin, origin, os );
     //   printf("echo %d\n", CanonicalRank(f1, origin, os));
		return 0.;
    }
    
   while (1){
        if (canonicalDecompositionMP(rank, f1, origin, os,alloy, spin, tolerance*value2) ){
#if 1
            printf("bailed %d %d %d %d -- %d\n",origin,0,alloy,spin , CanonicalRank(f1, alloy, spin));
            fflush(stdout);
#endif
            f1->sinc.tulip[alloy].Current[spin]--;
            return 1;
        }
        value = distanceOne(rank, f1, alloy, spin, origin, os);
    //	printf("%1.15f \t %1.15f : %d <%d %d\n", value, f1->mem1->rt->TARGET*value2,CanonicalRank(f1, alloy, spin ),part(f1, alloy),maxRun);
	    if ( fabs(value ) > f1->mem1->rt->TARGET*value2 && CanonicalRank(f1, alloy, spin ) < maxRun )
            tId(f1, alloy, spin);
        else
            return 0;
    }
    
    return 0;
}

double canonicalListDecompositionMP( INT_TYPE rank,struct field * f1 , Stream_Type * cofact, enum division origin,INT_TYPE os,   enum division alloy ,  INT_TYPE spin ,double tolerance,double magn){
    os = 0;
    if ( tolerance < 0 ){
        tolerance = -(tolerance);
        
        //   exit(0);
    }
    //printf("%d (%d)\n", alloy,CanonicalRank(f1, alloy, 0));
    //fflush(stdout);
    INT_TYPE g,l,l2,count = 0 ;
    
    INT_TYPE ns = 9;
    double ALPHA = f1->mem1->rt->ALPHA;

    Stream_Type * array[SPACE];
    Stream_Type * array2[SPACE];
    Stream_Type * array3;
    INT_TYPE space,space0,space1,space2;
    //    printf("%lld %lld %lld\n", spin,  f1->sinc.tulip[origin].Current[spin], f1->sinc.tulip[alloy].Current[spin]);
    INT_TYPE L1 = CanonicalRank(f1,alloy,spin);;
    INT_TYPE G1 =CanonicalRank(f1,origin,os);
    if ( G1 == 0 ){
        f1->sinc.tulip[alloy].Current[spin] = 0;
        printf("CD: empty origin %d _%lld -> %d _%lld\n", origin, os , alloy, spin);
        return 0;
    }
    
    if ( L1 == 0 ){
        // printf("CD: zero length %lld %lld %lld\n", origin,alloy,spin);
        tId(f1, alloy, spin);
        L1 = 1;
    }
    
    INT_TYPE M2[3];
    length(f1,alloy,M2);
    double value2 = 1;
    INT_TYPE info;
    double rcond;
    enum division canonicalStore;
    INT_TYPE T1 = G1*L1;
    
    {//REF CORE NUMBER
        //GET CORE... ALLOCATE TO CORE-LANE.
        
        
        for ( space = 0; space < SPACE ; space++){
            array[space] =  streams(f1, canonicalBuffers, rank , space);
            array2[space] =  streams(f1, canonicalBuffers, rank , space) + L1*L1;
        }
        
        
        if (  T1 + L1*L1 >  part(f1,canonicalBuffers)){
#if 1
            printf("mem prob\n %lld+ %lld  > %lld \n",T1,L1*L1,part(f1,canonicalBuffers) );
#endif
            exit(0);
        }
        if (species(f1, alloy ) == vector && bodies(f1, alloy ) == bodies(f1,canonicalBuffersB ) ) {
            canonicalStore = canonicalBuffersB;
        }else {
            printf("dow!\n");
            printf("%d %d %d %d\n", origin, alloy,species(f1, alloy ),bodies(f1, alloy ));
            exit(0);
        }
        if ( part(f1, canonicalStore) < L1 || alloc(f1, canonicalStore) < alloc(f1, alloy) )
        {
#if 1
            printf("%d -> %d (%lld: %lld) \n", origin, alloy,L1,part(f1, canonicalStore));
            printf("memory PROB %lld %lld \n",alloc(f1, canonicalStore) , alloc(f1, alloy));
#endif
            exit(0);
        }
        
        
    }
    
    

    
    enum division alloyBak;
    if ( species (f1, alloy ) == vector && bodies(f1,alloy) == one){
        alloyBak = trainVector;
    }
    else     if ( species (f1, alloy ) == vector && bodies(f1,alloy) == two){
        alloyBak = trainVector2;
    }
    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == one){
        
        f1->sinc.tulip[trainMatrix].header = header(f1,alloy);
        
        alloyBak = trainMatrix;
    }else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == two){
        
        f1->sinc.tulip[trainMatrix2].header = header(f1,alloy);
        
        alloyBak = trainMatrix2;
    }
    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == three){
        
        f1->sinc.tulip[trainVector3].header = header(f1,alloy);
        
        alloyBak = trainVector3;
    }
    
    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == four){
        
        f1->sinc.tulip[trainVector4].header = header(f1,alloy);
        
        alloyBak = trainVector4;
    }
    
    else if ( species(f1, alloy ) == quartic ){
        f1->sinc.tulip[trainQuartic].header = header(f1,alloy);
        
        alloyBak = trainQuartic;
    } else{
        printf("cd: tObject species %d \n", alloy);
        
        exit(0);
        
    }
    
    
    if ( spread(f1,origin,os,alloy,spin,1,array[1],array2[1]) )
        return -1;
    do{
        if ( count % ns == 0 )
            tEqua(f1, alloyBak,rank, alloy, spin );
        
        for ( space0 = 0; space0 < SPACE ;space0++){
            
            if ( SPACE == 3 ){
                if ( space0 == 0 ){
                    space0 = 0;
                    space1 = 1;
                    space2 = 2;
                } else if ( space0 == 1 ){
                    space0 = 1;
                    space1 = 2;
                    space2 = 0;
                }else if ( space0 == 2 ){
                    space0 = 2;
                    space1 = 0;
                    space2 = 1;
                } else {
                    printf("cd space");
                    exit(0);
                }
            }
            else if ( SPACE == 2 ){
                if ( space0 == 0 ){
                    space0 = 0;
                    space1 = 1;
                } else if ( space0 == 1 ){
                    space0 = 1;
                    space1 = 0;
                }else {
                    printf("cd space");
                    exit(0);
                }
            }
            
            
            
            
            
            
            //1
//            if ( normalize(f1,alloy,spin,space1) )
//                return -1;
            //            if ( mpSpread(f1,origin,os,alloy,spin,space1,array[space1],array2[space1]) )
            //                return -1;
            
            if ( SPACE == 3 ){
                //normal 3d calculator
                if ( spread(f1,origin,os,alloy,spin,space1,array[space1],array2[space1]) )
                    return -1;

                
                if ( spread(f1,origin,os,alloy,spin,space2,array[space2],array2[space2]) )
                    return -1;
                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1*L1, 0,array[space2],1, array[space1],1 );
                //replace Matrix
                if ( L1 >1 )
                    for ( l = 0 ; l < L1; l++){
#if VERBOSE
                        if ( fabs(array[space1][ l*L1 + l ]-1.0) > 1e-6 ||  isnan ( array[space1][l*L1+l] )  ||  isinf ( array[space1][l*L1+l] ) ){
                            printf("%lld:%f (%lld)-\n", l,array[space1][ l*L1 + l ],rank);
                            fflush(stdout);
                        }
#endif
                        array[space1][ l*L1 + l ] += ALPHA;
                    }
                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,T1   , 0,array2[space2],1, array2[space1],1 );
                //replace Vectors
            }else if (SPACE == 2){
                //plane-calculator ..  third dimenison is absent, allowing for simplication.
                
                if ( normalize(f1,alloy,spin,space1) )
                    return -1;

                if ( spread(f1,origin,os,alloy,spin,space1,array[space1],array2[space1]) )
                    return -1;
                
                if ( L1 >1 )
                    for ( l = 0 ; l < L1; l++){
                        if ( fabs(array[space1][ l*L1 + l ]-1.0) > 1e-6 ||  isnan ( array[space1][l*L1+l] )  ||  isinf ( array[space1][l*L1+l] ) ){
#if VERBOSE
                            printf("%lld:%f (%lld)-\n", l,array[space1][ l*L1 + l ],rank);
                            fflush(stdout);
#endif
                        }
                        array[space1][ l*L1 + l ] += ALPHA;
                    }

            }
            
            for ( g = 0; g < G1 ; g++)
                for ( l = 0; l < L1 ; l++){
                    array2[space1][g*L1+l] *= cofact[g];
                }
            
            
            
            // Vectors  L1 x G1
            
            // list...  L1 x M2 ==   ( cross * gstream**T )
            
            
            
            if ( 1 )
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,L1,M2[space0],G1,1,array2[space1],L1,streams(f1,origin,os,space0),M2[space0], 0, myStreams(f1, canonicalStore, rank), L1 );
            
            
            info = tdpotrf(L1, array[space1]);
            
            if ( info != 0 ){
#if 1
                printf("info failure %lld %lld %d %lld %d %lld %f\n", info, rank,origin,os, alloy, spin, ALPHA);
#endif
                return 1;
            }
            rcond = tdpocon(rank, f1, L1, array[space1] );
           // printf("%f %d %f\n", rcond,L1,ALPHA);
//
//            if ( f1->mem1->rt->targetCondition > 0 && L1 > 1 )
//                ALPHA /= min(1e1,max(1e-1, 1 - log(rcond*f1->mem1->rt->targetCondition) ));

            if (  isnan(rcond) || isinf(rcond) || rcond < 1e-12){
#if 1
                printf("Warning: condition of Beylkin %d->%d is bad %16.16f\n", origin, alloy, rcond);
                fflush(stdout);
#endif
                return 1;
            }

            info = tdpotrs(L1,  M2[space0], array[space1],  myStreams(f1, canonicalStore, rank) );
            if ( info != 0 ){
#if 1
                printf("L1 %lld \n", L1);
                printf("M2 %lld \n", M2[0]);
                
                printf("info failture %lld\n", info);
#endif
                return 1;
            }
            transpose ( L1, M2[space0] , myStreams(f1, canonicalStore, rank) , streams(f1,alloy,spin,space0));
        }
        balance(f1, alloy , spin );

        if ( count % ns == 0 ){
            ns++;
            value2 = distanceOne(rank, f1, alloy, spin,alloyBak, rank);
        }
        count++;
    } while (  value2 > tolerance && ns < 16) ;
    
    return 0;
}

double tCycleDecompostionListOneMP ( INT_TYPE rank, struct field * f1 , enum division origin, double * coeff, enum division alloy,INT_TYPE spin,  double tolerance , INT_TYPE maxRun , double power  ){
    double value,value2 = cblas_dnrm2 ( CanonicalRank(f1, origin, 0 ), coeff, 1);
    if ( value2 < 1e-6 )
        return 0;
    value2 = tInnerListMP(rank, f1, origin, coeff);
    //printf("%f\n", value2);
    do{
        if (canonicalListDecompositionMP(rank, f1, coeff, origin, 0,alloy, spin, tolerance,value2) ){
#if 1
            printf("bailed %d %d %d %d -- %d\n",origin,0,alloy,spin , CanonicalRank(f1, alloy, spin));
            fflush(stdout);
#endif
            f1->sinc.tulip[alloy].Current[spin]--;
            return 1;
        }
        value = fabs(value2 - 2. * tInnerVectorListMP(rank, f1, origin, coeff, alloy, spin) + inner(rank, f1, alloy, spin) );
     //  printf("%d %d %f %f %f %d\n", alloy, spin,value ,tInnerVectorListMP(rank, f1, origin, coeff, alloy, spin),inner(rank, f1, alloy, spin) , CanonicalRank(f1, alloy, spin));
        if ( fabs(value ) > f1->mem1->rt->TARGET*value2 && CanonicalRank(f1, alloy, spin ) < maxRun )
                tId(f1, alloy, spin);
            else
                return 0;
        //printf(".");
    }while(1);
    return 0.;
}

double tInnerVectorListMP( INT_TYPE rank, struct field * f1 , enum division origin, double * coeff, enum division vector,INT_TYPE spin ){
    
    double sum = 0.;
    INT_TYPE i = 0,info;
    
    for ( i = 0; i < CanonicalRank(f1, origin,0 ) ; i++ ) {
        tEqua(f1, copyVector, rank, ocean(rank, f1, origin, i, 0), 0);

        sum += coeff[i]*tMultiplyMP(rank, &info, f1, -1, 1.0, nullVector, 0, 'T', copyVector,rank, 'N', vector, spin);
    }
    
    return sum;
}

double tInnerListMP( INT_TYPE rank, struct field * f1 , enum division origin, double * coeff ){
    
    double sum = 0.;
    INT_TYPE j,i = 0,info;
    
    for ( i = 0; i < CanonicalRank(f1, origin,0 ) ; i++ )
        for ( j = 0 ;j < CanonicalRank(f1, origin, 0); j++){
            tEqua(f1, copyVector, rank, ocean(rank, f1, origin, i, 0), 0);
            tEqua(f1, copyTwoVector, rank, ocean(rank, f1, origin, j, 0), 0);

            sum += coeff[i]*coeff[j]*tMultiplyMP(rank, &info, f1, -1, 1.0, nullVector, 0, 'T', copyVector, rank, 'N', copyTwoVector, rank);
    }
    
    return sum;
}

INT_TYPE printExpectationValues (struct field * f1 , enum division ha  , enum division vector){
    INT_TYPE cat = 0,cmpl,cmpl2,totcan[2];
    DCOMPLEX co,expat,totx=0.;
    enum division leftP = ha,Mat;
    totcan[0] = 0;
    totcan[1] = 0;
    double norm2 = inner(0, f1, vector, 0)+inner(0, f1, vector, 1);
    do {
        
        if ( linear == name(f1,leftP) ){
            Mat = rivers(0, f1, linear, cat );
            f1->sinc.tulip[Mat].blockType = f1->sinc.tulip[leftP].blockType;
        }else {
            Mat = leftP;
        }
        for ( cmpl = 0; cmpl < spins(f1,Mat)  ;cmpl++)
            if ( CanonicalRank(f1, name(f1,Mat), cmpl)){
                expat = 0.;
                for ( cmpl2 = 0 ; cmpl2 < spins(f1,vector)  ;cmpl2++)
                    if ( CanonicalRank(f1,  Mat, cmpl)){
                        co = 1.;
                        if ( cmpl )
                            co *= I;
                        if ( cmpl2 )
                            co *= -I;
                        totcan[cmpl] += CanonicalRank(f1, Mat, cmpl);
                        expat += co* matrixElements(0, f1, 'T', vector, 'N', Mat, cmpl, vector, cmpl2);
                        
                    }
                totx += expat/(DCOMPLEX)norm2;
                printf("**%d-%d-%d-%d\t%d\t%d\t%d-%d\t%d\t%f\t%f\n",Mat,name(f1,Mat),cmpl,cmpl2, bodies(f1, Mat),f1->sinc.tulip[Mat].blockType,CanonicalRank(f1, Mat, cmpl),f1->sinc.tulip[Mat].ptRank[cmpl], cat,creal(expat/norm2),cimag(expat/norm2));
            }
        
        if (name(f1, leftP) == linear ){
            cat++;
            
            if ( cat > f1->Na ){
                leftP = f1->sinc.tulip[leftP].linkNext;
                cat = 0;
            }
            
        } else {
            leftP = f1->sinc.tulip[leftP].linkNext;
        }
    } while ( leftP != nullName);
    printf("**%d-%d-%d-%d\t%d\t%d\t%d-%d\t%d\t%f\t%f\n",0,0,0,0, bodies(f1, vector),0,totcan[0],totcan[1], 0,creal(totx),cimag(totx));

    return 0;
}

DCOMPLEX matrixElements ( INT_TYPE rank,struct field * f1 , char perm, enum division uec, char mc, enum division mat,INT_TYPE ms, enum division vec,INT_TYPE vs){
    INT_TYPE info,r;
    double value;
    if ( CanonicalRank(f1, vec, vs) * (CanonicalRank(f1, uec, 0 ) +CanonicalRank(f1, uec,1)  )== 0 )
        return 0.;
    if ( vec == productVector || uec == productVector){
        printf("me\n");
        exit(0);
    }
    DCOMPLEX sum = 0.;
    
    for ( r = 0 ; r < CanonicalRank(f1, mat, ms);r++){
        tMultiplyMP(rank, &info, f1, 1., -1, productVector, rank, mc, ocean(rank, f1, mat, r, ms) ,ms, 'N', vec, vs);
        value = tMultiplyMP(rank, &info, f1, 1., -1, nullVector, 0, 'T', uec ,0, 'N', productVector,rank );
        if ( CanonicalRank(f1, uec,1 ))
            value += I*tMultiplyMP(rank, &info, f1, 1., -1, nullVector, 0, 'T', uec ,1, 'N', productVector,rank );
        sum += value;
    }
    return sum;
}



double canonicalMultiplyMP( INT_TYPE rank,struct field * f1 , char mc,enum division mat,INT_TYPE ms, enum division vec, INT_TYPE vs,   enum division alloy ,  INT_TYPE spin ,double tolerance){
    if ( tolerance < 0 ){
        //printf("%d->%d\n", origin, alloy);
        tolerance = -(tolerance);
        
        //   exit(0);
    }
    if ( vec == productVector || alloy == productVector ){
        printf("canonicalMulMP\n");
        exit(0);
    }
    INT_TYPE ns = 9;
    Stream_Type * array[3];
    Stream_Type * array2[3];
    INT_TYPE space,space0,space1,space2;
    INT_TYPE l,count = 1 ;
    INT_TYPE L1 = CanonicalRank(f1,alloy,spin);;
    INT_TYPE R1 = CanonicalRank(f1, mat,ms);
    INT_TYPE V1 = CanonicalRank(f1, vec, vs);
    if ( R1*V1 == 0 ){
     //   printf("%d %d %d\n", mat, f1->sinc.tulip[mat].name, name(f1,mat));
     //   printf("%d %d %d\n", CanonicalRank(f1,mat,0), CanonicalRank(f1, name(f1,mat),0));

        f1->sinc.tulip[alloy].Current[spin] = 0;
        return 1;
    }
    double ALPHA = f1->mem1->rt->ALPHA;

    if ( L1 == 0 ){
        //   printf("CD: zero length %lld %lld %lld\n", origin,alloy,spin);
        tId(f1, alloy, spin);
        L1 = 1;
        if (! CanonicalRank(f1, alloy, spin)){
            printf("hurm... %d\n",alloy);
            exit(0);
        }
    }
    INT_TYPE M2[3];
    length(f1,alloy,M2);
    double value2= 1;
    INT_TYPE info;
    double rcond;
    
    {//REF CORE NUMBER
        //GET CORE... ALLOCATE TO CORE-LANE.
        
        for ( space = 0; space < SPACE ; space++){
            array[space] =  streams(f1, canonicalBuffers, rank , space);
            array2[space] =  streams(f1, canonicalBuffers, rank , space) + L1*L1;

        }
        
        if (  L1*L1 + L1*imax(V1,M2[0])  >  part(f1,canonicalBuffers)){
#if 1
            printf("mem prob\n %lld+ %lld  > %lld \n",R1*V1,L1*L1,part(f1,canonicalBuffers) );
#endif
            exit(0);
        }
        
//            canonicalStore = canonicalBuffersB;
//            canonicalStore = canonicalBuffersBX;
//            canonicalStore = canonicalBuffersBM;
        
    }
    
    enum division alloyBak;
    if ( species (f1, alloy ) == vector && bodies(f1,alloy) == one){
        alloyBak = trainVector;
    }
    else     if ( species (f1, alloy ) == vector && bodies(f1,alloy) == two){
        alloyBak = trainVector2;
    }
    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == one){
        
        f1->sinc.tulip[trainMatrix].header = header(f1,alloy);
        
        alloyBak = trainMatrix;
    }else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == two){
        
        f1->sinc.tulip[trainMatrix2].header = header(f1,alloy);
        
        alloyBak = trainMatrix2;
    }
    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == three){
        
        f1->sinc.tulip[trainVector3].header = header(f1,alloy);
        
        alloyBak = trainVector3;
    }
    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == four){
        
        f1->sinc.tulip[trainVector4].header = header(f1,alloy);
        
        alloyBak = trainVector4;
    }
    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == three){
        
        f1->sinc.tulip[trainMatrix3].header = header(f1,alloy);
        
        alloyBak = trainMatrix3;
    }
    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == four){
        
        f1->sinc.tulip[trainMatrix4].header = header(f1,alloy);
        
        alloyBak = trainMatrix4;
    }
    else if ( species(f1, alloy ) == quartic ){
        f1->sinc.tulip[trainQuartic].header = header(f1,alloy);
        
        alloyBak = trainQuartic;
    } else{
        printf("cd: tObject species %d \n", alloy);
        
        exit(0);
        
    }
    if (  ( purpose(f1, alloy ) < tObject   )){
        printf("cd: tObject status %d \n", alloy);
        exit(0);
        
    }
    // printf(":: %d %d -- %d %d \n", origin , os, alloy ,spin );
    space2 = 1;
//    if ( normalize(f1,alloy,spin,space2) )
//        return -1;
//
//    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, L1, M2[space2], 1.0, streams(f1, alloy,spin,space2), M2[space2], streams(f1,alloy,spin,space2), M2[space2], 0.0, array[space2], L1);

    
    do{
        if ( count % ns == 0 )
            tEqua(f1, alloyBak,rank, alloy, spin );
        
        for ( space0 = 0; space0 < SPACE ;space0++){
            
            if ( SPACE == 3 ){
                if ( space0 == 0 ){
                    space0 = 0;
                    space1 = 1;
                    space2 = 2;
                } else if ( space0 == 1 ){
                    space0 = 1;
                    space1 = 2;
                    space2 = 0;
                }else if ( space0 == 2 ){
                    space0 = 2;
                    space1 = 0;
                    space2 = 1;
                } else {
                    printf("cd space");
                    exit(0);
                }
            }
            else if ( SPACE == 2 ){
                if ( space0 == 0 ){
                    space0 = 0;
                    space1 = 1;
                } else if ( space0 == 1 ){
                    space0 = 1;
                    space1 = 0;
                }else {
                    printf("cd space");
                    exit(0);
                }
            }
            
            
            
//            if ( normalize(f1,alloy,spin,space1) )
//                return -1;
            //?? MAY NOT NEED
            
            
            
            if ( SPACE == 3 ){
                //normal 3d calculator
                
                if ( normalize(f1,alloy,spin,space2) )
                    return -1;
                
                if ( normalize(f1,alloy,spin,space1) )
                    return -1;

                
                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, L1, M2[space2], 1.0, streams(f1, alloy,spin,space2), M2[space2], streams(f1,alloy,spin,space2), M2[space2], 0.0, array[space2], L1);
                
                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, L1, M2[space1], 1.0, streams(f1, alloy,spin,space1), M2[space1], streams(f1,alloy,spin,space1), M2[space1], 0.0, array[space1], L1);

                
                cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,L1*L1, 0,array[space2],1, array[space1],1 );
                //replace Matrix
                
                if ( L1 > 1 )
                    for ( l = 0 ; l < L1; l++){
                        if ( fabs(array[space1][ l*L1 + l ]-1.0) > 1e-6 ||  isnan ( array[space1][l*L1+l] )  ||  isinf ( array[space1][l*L1+l] ) ){
#if VERBOSE
                            printf("%lld:%f (%lld)-\n", l,array[space1][ l*L1 + l ],rank);
                            fflush(stdout);
#endif
                        }
                        array[space1][ l*L1 + l ] += ALPHA;
                    }
                
            }else if ( SPACE == 2 ){
                //normal 3d calculator
                
                
                if ( normalize(f1,alloy,spin,space1) )
                    return -1;
                
                
                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, L1, M2[space1], 1.0, streams(f1, alloy,spin,space1), M2[space1], streams(f1,alloy,spin,space1), M2[space1], 0.0, array[space1], L1);
                
                
                if ( L1 > 1 )
                    for ( l = 0 ; l < L1; l++){
                        if ( fabs(array[space1][ l*L1 + l ]-1.0) > 1e-6 ||  isnan ( array[space1][l*L1+l] )  ||  isinf ( array[space1][l*L1+l] ) ){
#if VERBOSE
                            printf("%lld:%f (%lld)-\n", l,array[space1][ l*L1 + l ],rank);
                            fflush(stdout);
#endif
                        }
                        array[space1][ l*L1 + l ] += ALPHA;
                    }
                
            }
            
            info = tdpotrf(L1, array[space1]);
            
            if ( info != 0 ){
#if 1
                printf("info failure %lld %lld %d %lld %d %lld\n", info, rank,mat,vec, alloy, spin);
#endif
                return 1;
            }
            rcond = tdpocon(rank, f1, L1, array[space1] );
//            if ( f1->mem1->rt->targetCondition > 0 )
////
//            ALPHA /= min(1e1,max(1e-1, 1 - log(rcond*f1->mem1->rt->targetCondition) ));
            
            
            if (  isnan(rcond) || isinf(rcond) || rcond < 1e-12){
#if 1
                printf("Warning: condition of Beylkin %d->%d is bad %16.16f\n", mat, alloy, rcond);
                fflush(stdout);
#endif
                return 1;
            }
            
            INT_TYPE r;
            for ( r = 0; r < L1*M2[space0] ; r++)
                array2[space0][r] = 0.;
            for ( r = 0 ; r < R1 ; r++){
            //    printf("%d %d\n", name(f1, mat), name(f1,ocean(rank,f1,mat,r,ms)));
                //printf("(%d %d %d) \n", mat , purpose(f1,namat), name(f1,mat ));
                //fflush(stdout);
//                if (  purpose(f1,mat) == ptObject ){
//                    if ( name(f1,mat ) == linear ){
//                        tEqua(f1 , copy, rank, ocean(rank,f1,mat,r,ms),ms);
//                        tMultiplyMP(rank, &info, f1, 1., -1, productVector, rank, mc,copy,rank, 'N', vec, vs);
//                    }
//                    else
//                        exit(0);
//                }else {
                
                    tMultiplyMP(rank, &info, f1, 1., -1, productVector, rank, mc,ocean(rank,f1,mat,r,ms),ms, 'N', vec, vs);
//                }
                
                {
                    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, V1, M2[space1], 1.0,  streams(f1,alloy,spin,space1), M2[space1], streams(f1, productVector,rank,space1), M2[space1],0.0, array2[space1], L1);
                    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, L1, V1, M2[space2], 1.0,  streams(f1,alloy,spin,space2), M2[space2], streams(f1, productVector,rank,space2), M2[space2],0.0, array2[space2], L1);

                    cblas_dtbmv(CblasColMajor, CblasUpper,CblasNoTrans,CblasNonUnit,V1*L1, 0,array2[space2],1, array2[space1],1 );
                }//form <g,f>
                
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,L1,M2[space0],V1,1.0,array2[space1],L1,streams(f1,productVector,rank,space0),M2[space0], 1.00000, array2[space0], L1 );
                
            }
            info = tdpotrs(L1,  M2[space0], array[space1],  array2[space0] );
            if ( info != 0 ){
#if 1
                printf("L1 %lld \n", L1);
                printf("M2 %lld \n", M2[0]);
                printf("R1 %lld \n", V1*L1);
                printf("info failture %lld\n", info);
#endif
                exit(0);
                return 1;
            }
            
            transpose ( L1, M2[space0] , array2[space0] , streams(f1,alloy,spin,space0));
        
        }
        balance(f1, alloy , spin );
        
        if ( count % ns == 0 ){
            value2 = distanceOne(rank, f1, alloy, spin,alloyBak, rank);
            ns++;
        }
        count++;
    } while (  value2 >tolerance && ns < 16) ;
    return 0;
    
}


double tCycleMultiplyMP ( INT_TYPE rank,struct field * f1 , char mc,enum division mat,INT_TYPE ms, enum division vec, INT_TYPE vs,   enum division alloy ,  INT_TYPE spin ,double tolerance, INT_TYPE maxRun, double power){
    double value,value2 ;
    if ( ! CanonicalRank(f1, mat, ms ) )//likely has zer
        return 0;

    do{
        if (canonicalMultiplyMP(rank, f1, mc,mat,ms, vec, vs,alloy, spin, tolerance) ){
#if 1
            printf("bailed %d %d %d %d -- %d\n",mat,ms,alloy,spin , CanonicalRank(f1, alloy, spin));
            fflush(stdout);
#endif
            return 1;
        }
        if ( CanonicalRank(f1, alloy, spin ) < maxRun )
            tId(f1, alloy, spin);
        else
            return 0;
    }while(1);
    return 0.;
}



INT_TYPE tOuterProductSu( struct field * f1,enum division vector , INT_TYPE a, enum division vector2,INT_TYPE b, enum division proj, INT_TYPE c){
    INT_TYPE ma = CanonicalRank(f1, vector,a), mb = CanonicalRank(f1, vector,b), zc = CanonicalRank(f1, proj, c),l,r,space,i;
    INT_TYPE *n;
    if ( bodies(f1, vector ) == one )
        n = f1->sinc.Basis;
    else if ( bodies(f1, vector ) == two )
        n = f1->sinc.Basis2;
    else if ( bodies(f1, vector ) == three )
        n = f1->sinc.Basis3;
    else if ( bodies(f1, vector ) == four )
        n = f1->sinc.Basis4;
    else
        exit(0);
    INT_TYPE N2[3];
    INT_TYPE n2[3];
    length(f1, a, n2);
    length(f1, proj, N2);
    if ( species(f1, vector ) == matrix || species(f1, vector2) == matrix){
        printf("outer\n");
        exit(0);
    }
    if ( zc + ma*mb > part(f1, proj) && proj < eigenVectors){
        printf("outerProductSu:  %lld > %lld\n",  zc + ma*mb , part(f1, proj) );
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

double tMultiplyMP (INT_TYPE rank, INT_TYPE * info, struct field * f1,double number, INT_TYPE beta,  enum division equals, INT_TYPE espin ,char leftChar, enum division left,INT_TYPE lspin, char rightChar,enum division right, INT_TYPE rspin){
    INT_TYPE LN2[SPACE],RN2[SPACE],EN2[SPACE],space,l,r,ii ,type ;
    INT_TYPE flag;
    number = 1.;
    *info = 0;
    enum CBLAS_TRANSPOSE leftSymbol,rightSymbol,skewSymbol;
    
    if (leftChar == 'N')
        leftSymbol = CblasNoTrans;
    else if (leftChar == 'T')
        leftSymbol = CblasTrans;
    if (rightChar == 'N')
        rightSymbol = CblasNoTrans;
    else if (rightChar == 'T')
        rightSymbol = CblasTrans;
    
    
    
    if ( espin < 0 || lspin < 0 || rspin < 0 ){
        printf("poor spins!\n");
        exit(0);
    }
    INT_TYPE LL = CanonicalRank(f1,left,lspin);
    INT_TYPE LR = CanonicalRank(f1,right,rspin);
    
    
    INT_TYPE *N1,NX[SPACE],A[SPACE],B[SPACE],AA[SPACE],BB[SPACE],flagTranspose,flagTranspose2,flagTranspose3,flagTranspose4;
    INT_TYPE *n1 = f1->sinc.Basis;
    length(f1,left,LN2 );
    length(f1,right,RN2 );
    length(f1,equals,EN2 );
    
    INT_TYPE MM = LL*LR;
    INT_TYPE reduce,gamma;
    
    double sum =0.,product;
    if ( beta < 0 ){
        reduce = 1;
        gamma = -(beta+1);
        beta = 0;
        
    } else{
        reduce = 0;
        gamma = 0;
    }
    
    if ( header(f1, left ) != header(f1, right )){
        printf("Two Head types\n");
        
        printf("%d -%lld (%lld) %d (%lld)\n", left, name(f1,left),header(f1,left), right,header(f1,right));
        exit(0);
        return 0.;
        
    }
    
    
    Stream_Type * array[SPACE][5];
    
    if (MM+beta*MM +gamma*CanonicalRank(f1,equals,espin) > part(f1, equals) && species(f1, equals) != scalar && purpose(f1, equals)!= ptObject){
        //tGetType allocation!
        printf("tMultiply name %d %d %d\n", equals, left, right);
        printf("%lld >  %lld: %d %d %u\n",MM+beta*MM+gamma*CanonicalRank(f1,equals,espin) , part(f1, equals),LL,LR, streams(f1, left,0,0) );
        printf("%lld %lld %lld %lld\n",part(f1, equals), LL, LR,species(f1, equals));
        printf("%lld %lld \n",part(f1, left),part(f1, right) );
        printf("%c %d %d\n", leftChar, lspin, rspin);
       // exit(0);
       // return 1;
    }
    
    type = -1;
    if (  (species(f1,equals) == scalar && species(f1, left) == species(f1, right) && bodies(f1,left) == bodies(f1,right))){
        type = 0;
        N1 = vectorLen(f1, left);
        
        if ( rightChar != 'N' )
        {
            printf("fix your character! %d %d %d\n",equals, left , right);
            exit(0);
        }
        flag = (leftChar == 'N' && species(f1, left ) == matrix&& bodies(f1,left)  == one );
        if  ( ! flag )
            flag = 2 *(leftChar != 'T' && species(f1, left ) == vector && bodies(f1,left)  == two );
        if ( ! flag )
            flag = 3 *(leftChar != 'T'&&  species(f1, left ) == vector && bodies(f1,left)  == three );
        if ( ! flag )
            flag = 4 *( leftChar != 'T' && species(f1, left ) == vector && bodies(f1,left)  == four );
        
        if ( flag == 1 )
        {
            
            if ( species(f1, tensorBuffers ) >= vector && bodies(f1,tensorBuffers)  >= one && part(f1, tensorBuffers)>=1){
                for ( space =0 ; space < SPACE ; space++)
                    array[space][0] = streams( f1, tensorBuffers, rank,space );
            }else {
                printf("oops\n");
                exit(0);
            }
            
        }
        if ( flag == 2)
        {
            N1 = n1;

            //            if ( header (f1,left ) != Cube ){
            //                printf("op");
            //                exit(0);
            //            }
            if ( ( species(f1, tensorBuffers ) >= vector && bodies(f1,tensorBuffers)  >= two && part(f1, tensorBuffers)>=1) ){
//                printf("tens %d %d %d %d\n", name(f1, tensorBuffers), bodies(f1, tensorBuffers),species(f1, tensorBuffers),part(f1,tensorBuffers));
//                printf("left %d %d %d %d\n", name(f1, left), bodies(f1, left),species(f1, left),part(f1,left));
//                printf("right %d %d %d %d\n", name(f1, right), bodies(f1, right),species(f1, right),part(f1,right));

                for ( space =0 ; space < SPACE ; space++)
                    array[space][0] = streams( f1, tensorBuffers, rank,space );
            }else {
                printf("oopss\n");
                exit(0);
            }
            
        }
        if ( flag == 3)
        {
            INT_TYPE r;
            N1 = n1;
            
            
            
            //            if ( header (f1,left ) != Cube ){
            //                printf("op");
            //                exit(0);
            //            }
            
            
            //
            //            train[0] = 'T';//(1)            123
            //            train[1] = 'A';//(123)          231
            //            train[2] = 'B';//(123).(123)    312
            //            train[3] = 'C';//(12)           213
            //            train[4] = 'D';//(13)           321
            //            train[5] = 'E';//(23)           132
            
            flagTranspose3 = 0;
            for ( space = 0; space < SPACE ; space++){
                A[space] = N1[space]*N1[space];
                B[space] = N1[space];
            }
            
            if ( leftChar == 'T' ){
                flagTranspose = 0;
                flagTranspose2 = 0;
            }
            
            else  if ( leftChar == 'A' ){
                flagTranspose = 1;// a | b c
                
                flagTranspose2 = 0;
                for ( space = 0; space < SPACE ; space++){
                    A[space] = N1[space];
                    B[space] = N1[space]*N1[space];
                }
            }else  if ( leftChar == 'B' ){
                flagTranspose = 1; // a b | c
                flagTranspose2 = 0;
                for ( space = 0; space < SPACE ; space++){
                    A[space] = N1[space]*N1[space];
                    B[space] = N1[space];
                }
            } else if ( leftChar == 'C' ){
                flagTranspose = 0;
                flagTranspose2 = 1;
            } else  if ( leftChar == 'D' ){
                flagTranspose = 1; // a b | c
                flagTranspose2 = 1;
                for ( space = 0; space < SPACE ; space++){
                    A[space] = N1[space]*N1[space];
                    B[space] = N1[space];
                }
            }else  if ( leftChar == 'E' ){
                flagTranspose = 1; // a b | c
                flagTranspose2 = 0;
                flagTranspose3 = 1;
                for ( space = 0; space < SPACE ; space++){
                    A[space] = N1[space]*N1[space];
                    B[space] = N1[space];
                }
            }
            
            else {
                printf("unknown flag %c\n",leftChar);
                exit(0);
            }
            
            
            
            if ( ( species(f1, tensorBuffers ) >= vector && bodies(f1,tensorBuffers)  >= three && part(f1, tensorBuffers)>=2) ){
                for ( space =0 ; space < SPACE ; space++)
                    for ( r = 0; r < 2 ; r++)
                        array[space][r] = streams( f1, tensorBuffers, rank,space)+LN2[space]*r;
            }else {
                printf("oopss\n");
                exit(0);
            }
            
        }
        else if (flag == 4 ){
            {
                INT_TYPE r;
                N1 = n1;
                flagTranspose = 0;
                flagTranspose2 = 0;
                flagTranspose3 = 0;
                flagTranspose4 = 0;
                
                for ( space = 0; space < SPACE ; space++){
                    A[space] = N1[space]*N1[space];
                    B[space] = N1[space]*N1[space];
                }
                for ( space = 0; space < SPACE ; space++){
                    AA[space] = N1[space];
                    BB[space] = N1[space]*N1[space]*N1[space];
                }
                //i <+> l
                //j <+> k
                
                
                
                if ( 1 ){
                    
                    // 0,0,0,0 : i,j,k,l   :: T 0
                    // 0,2,1,2 : i,j,l,k   :: 'a'24
                    // 1,3,1,2 : i,k,j,l   :: 'b'20 i,j,k,l->j,i,k,l->l,j,i,k->j,l,i,k->i,k,j,l
                    // 1,1,0,0 : i,k,l,j   :: 'c'4
                    // 0,3,1,0 : i,l,j,k   :: 'd'12
                    // 1,2,1,3 : i,l,k,j   :: 'e'22
                    // 1,0,0,0 : j,i,k,l   :: 'f'2
                    // 1,2,1,2 : j,i,l,k   :: 'g' 23
                    // 0,3,1,2 : j,k,i,l   :: 'h'19
                    // 0,1,0,0 : j,k,l,i   :: 'i'3
                    // 1,3,1,0 : j,l,i,k   :: 'j'13
                    // 0,2,1,3 : j,l,k,i   :: 'k'21
                    // 0,2,1,1 : k,i,j,l   :: 'l'15
                    // 1,1,1,0 : k,i,l,j   :: 'm'10
                    // 1,2,1,1 : k,j,i,l   :: 'n'16
                    // 0,1,1,0 : k,j,l,i   :: 'o'9
                    // 0,2,0,0 : k,l,i,j   :: 'p'5
                    // 1,1,0,1 : k,l,j,i   :: 'q'14
                    //c 0,3,0,0 : l,i,j,k   :: 'r'7
                    //c 1,3,1,1 : l,i,k,j   :: 's'18 i,j,k,l ->  j,i,k,l -> l,j,i,k -> j,l,i,k -> l,i,k,j
                    //c 1,3,0,0 : l,j,i,k   :: 't'8
                    //c 0,3,1,1 : l,j,k,i   :: 'u'17 i,j,k,l ->  l,i,j,k -> i,l,j,k -> l,j,k,i
                    //c 0,2,1,0 : l,k,i,j   :: 'v'6
                    //c 1,2,1,0 : l,k,j,i   :: 'w'11
                    
                    
                    
                    if ( leftChar == 'T' ){//(i,j,k,l)
                        
                    }
                    // 0,2,1,2 : i,j,l,k   :: 'a'24
                    else if ( leftChar == 'a' ){//(i,j,l,k)
                        flagTranspose2 = 1;
                        flagTranspose3 = 1;
                        flagTranspose4 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                        
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space]*N1[space];
                            BB[space] = N1[space]*N1[space];
                        }
                        
                    }
                    // 1,3,1,2 : i,k,j,l   :: 'b'20 i,j,k,l->j,i,k,l->l,j,i,k->j,l,i,k->i,k,j,l
                    else if ( leftChar == 'b' ){//(i,k,j,l)
                        flagTranspose = 1;
                        flagTranspose2 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                        flagTranspose3 = 1;
                        flagTranspose4 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space]*N1[space];
                            BB[space] = N1[space]*N1[space];
                        }
                        
                        
                    }
                    // 1,1,0,0 : i,k,l,j   :: 'c'4
                    else  if ( leftChar == 'c' ){//(i,k,l,j)
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space]*N1[space];
                        }
                    }
                    // 0,3,1,0 : i,l,j,k   :: 'd'12
                    else  if ( leftChar == 'd' ){//(i,l,j,k)
                        flagTranspose2 = 1;
                        flagTranspose3 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                    }
                    // 1,2,1,3 : i,l,k,j   :: 'e'22
                    else  if ( leftChar == 'e' ){//(i,l,k,j) 6
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                        flagTranspose3 = 1;
                        flagTranspose4 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space]*N1[space]*N1[space];
                            BB[space] = N1[space];
                        }
                        
                    }
                    // 1,0,0,0 : j,i,k,l   :: 'f'2
                    else  if ( leftChar == 'f' ){//(j,i,k,l)
                        flagTranspose = 1; // a b | c
                    }
                    // 1,2,1,2 : j,i,l,k   :: 'g' 23
                    else  if ( leftChar == 'g' ){//(j,i,l,k)
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                        flagTranspose3 = 1;
                        flagTranspose4 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space]*N1[space];
                            BB[space] = N1[space]*N1[space];
                        }
                        
                    }
                    // 0,3,1,2 : j,k,i,l   :: 'h'19
                    else  if ( leftChar == 'h' ){//(j,k,i,l)
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                        
                        flagTranspose3 = 1; // a b | c
                        flagTranspose4 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space]*N1[space];
                            BB[space] = N1[space]*N1[space];
                        }
                        
                    }
                    // 0,1,0,0 : j,k,l,i   :: 'i'3
                    else  if ( leftChar == 'i' ){//(j,k,l,i)
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space]*N1[space];
                        }
                    }
                    // 1,3,1,0 : j,l,i,k   :: 'j'13
                    else  if ( leftChar == 'j' ){//(j,l,i,k)
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        flagTranspose3 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                    }
                    // 0,2,1,3 : j,l,k,i   :: 'k'21
                    else  if ( leftChar == 'k' ){//(j,l,k,i)12
                        flagTranspose2 = 2;
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                        flagTranspose3 = 1;
                        flagTranspose4 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space]*N1[space]*N1[space];
                            BB[space] = N1[space];
                        }
                        
                    }
                    // 0,2,1,1 : k,i,j,l   :: 'l'15
                    else  if ( leftChar == 'l' ){//(k,i,j,l)
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                        
                        flagTranspose3 = 1;
                        flagTranspose4 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space];
                            BB[space] = N1[space]*N1[space]*N1[space];
                        }
                        
                    }
                    // 1,1,1,0 : k,i,l,j   :: 'm'10
                    else  if ( leftChar == 'm' ){//(k,i,l,j)
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        flagTranspose3 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space]*N1[space];
                        }
                    }
                    // 1,2,1,1 : k,j,i,l   :: 'n'16
                    else  if ( leftChar == 'n' ){//(k,j,i,l)
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                        
                        flagTranspose3 = 1;
                        flagTranspose4 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space];
                            BB[space] = N1[space]*N1[space]*N1[space];
                        }
                        
                        
                    }
                    // 0,1,1,0 : k,j,l,i   :: 'o'9
                    else  if ( leftChar == 'o' ){//(k,j,l,i)
                        flagTranspose2 = 1;
                        flagTranspose3 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space]*N1[space];
                        }
                    }
                    // 0,2,0,0 : k,l,i,j   :: 'p'5
                    else  if ( leftChar == 'p' ){//(k,l,i,j)
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                    }
                    // 1,1,0,1 : k,l,j,i   :: 'q'14
                    else  if ( leftChar == 'q' ){//(k,l,j,i)//18
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space]*N1[space];
                        }
                        flagTranspose4 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space];
                            BB[space] = N1[space]*N1[space]*N1[space];
                        }
                        
                    }
                    // 0,3,0,0 : l,i,j,k   :: 'r'7
                    else  if ( leftChar == 'r' ){//(l,i,j,k)
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                    }
                    // 1,3,1,1 : l,i,k,j   :: 's'18
                    else  if ( leftChar == 's' ){//(l,i,k,j)
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                        flagTranspose3 = 1; // a b | c
                        flagTranspose4 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space];
                            BB[space] = N1[space]*N1[space]*N1[space];
                        }
                        
                    }
                    // 1,3,0,0 : l,j,i,k   :: 't'8
                    else  if ( leftChar == 't' ){//(l,j,i,k)
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                    }
                    // 0,3,1,1 : l,j,k,i   :: 'u'17
                    else  if ( leftChar == 'u' ){//(l,j,k,i)
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                        
                        flagTranspose3 = 1;
                        flagTranspose4 = 1;
                        for ( space = 0; space < SPACE ; space++){
                            AA[space] = N1[space];
                            BB[space] = N1[space]*N1[space]*N1[space];
                        }
                    }
                    // 0,2,1,0 : l,k,i,j   :: 'v'6
                    else  if ( leftChar == 'v' ){//(l,k,i,j)
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                        flagTranspose3 = 1;
                    }
                    // 1,2,1,0 : l,k,j,i   :: 'w'11
                    else  if ( leftChar == 'w' ){//(l,k,j,i)//END
                        flagTranspose = 1; // a b | c
                        flagTranspose2 = 1;
                        flagTranspose3 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                    }
                    
                } else {
                    
                    
                    // 0,0,0,0 : i,j,k,l   :: T 0           l,k,j,i
                    // 0,2,1,2 : i,j,l,k   :: 'a'24         k,l,j,i  1,0,0,0    'f'
                    // 1,3,1,2 : i,k,j,l   :: 'b'           l,j,k,i  **         'b'
                    // 1,1,0,0 : i,k,l,j   :: 'c'4          j,l,k,i  1,3,0,0    'l'
                    // 0,3,1,0 : i,l,j,k   :: 'd'12         k,j,l,i             'h'
                    // 1,2,1,3 : i,l,k,j   :: 'e'22         j,k,l,i             'n'
                    // 1,0,0,0 : j,i,k,l   :: 'f'2          l,k,i,j             'a'
                    // 1,2,1,2 : j,i,l,k   :: 'g' 23        k,l,i,j             'g'
                    // 0,3,1,2 : j,k,i,l   :: 'h'19         l,i,k,j             'd'
                    // 0,1,0,0 : j,k,l,i   :: 'i'3          i,l,k,j             'r'
                    // 1,3,1,0 : j,l,i,k   :: 'j'13         k,i,l,j             'j'
                    // 0,2,1,3 : j,l,k,i   :: 'k'21         i,k,l,j             't'
                    // 0,2,1,1 : k,i,j,l   :: 'l'15         l,j,i,k             'c'
                    // 1,1,1,0 : k,i,l,j   :: 'm'10         j,l,i,k             'm'
                    // 1,2,1,1 : k,j,i,l   :: 'n'16         l,i,j,k             'e'
                    // 0,1,1,0 : k,j,l,i   :: 'o'9          i,l,j,k             's'
                    // 0,2,0,0 : k,l,i,j   :: 'p'5          j,i,l,k             'p'
                    // 1,1,0,1 : k,l,j,i   :: 'q'14         i,j,l,k             'v'
                    //c 0,3,0,0 : l,i,j,k   :: 'r'7         k,j,i,l             'i'
                    //c 1,3,1,1 : l,i,k,j   :: 's'          j,k,i,l             'o'
                    //c 1,3,0,0 : l,j,i,k   :: 't'8         k,i,j,l             'k'
                    //c 0,3,1,1 : l,j,k,i   :: 'u'          i,k,j,l             'u'
                    //c 0,2,1,0 : l,k,i,j   :: 'v'6         j,i,k,l             'q'
                    //c 1,2,1,0 : l,k,j,i   :: 'w'11        i,j,k,l             'w'
//                    {
//
//
//                        if ( leftChar == 'T' ){//(i,j,k,l)
//
//                        }
//                        // 0,2,1,2 : i,j,l,k   :: 'a'24
//                        else if ( leftChar == 'f' ){//(i,j,l,k)
//                            flagTranspose2 = 1;
//                            flagTranspose3 = 1;
//                            flagTranspose4 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space];
//                                B[space] = N1[space]*N1[space];
//                            }
//
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space]*N1[space];
//                                BB[space] = N1[space]*N1[space];
//                            }
//
//                        }
//                        // 1,3,1,2 : i,k,j,l   :: 'b'20 i,j,k,l->j,i,k,l->l,j,i,k->j,l,i,k->i,k,j,l
//                        else if ( leftChar == 'b' ){//(i,k,j,l)
//                            flagTranspose = 1;
//                            flagTranspose2 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space]*N1[space];
//                                B[space] = N1[space];
//                            }
//                            flagTranspose3 = 1;
//                            flagTranspose4 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space]*N1[space];
//                                BB[space] = N1[space]*N1[space];
//                            }
//
//
//                        }
//                        // 1,1,0,0 : i,k,l,j   :: 'c'4
//                        else  if ( leftChar == 'l' ){//(i,k,l,j)
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space];
//                                B[space] = N1[space]*N1[space]*N1[space];
//                            }
//                        }
//                        // 0,3,1,0 : i,l,j,k   :: 'd'12
//                        else  if ( leftChar == 'h' ){//(i,l,j,k)
//                            flagTranspose2 = 1;
//                            flagTranspose3 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space]*N1[space];
//                                B[space] = N1[space];
//                            }
//                        }
//                        // 1,2,1,3 : i,l,k,j   :: 'e'22
//                        else  if ( leftChar == 'n' ){//(i,l,k,j) 6
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space];
//                                B[space] = N1[space]*N1[space];
//                            }
//                            flagTranspose3 = 1;
//                            flagTranspose4 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space]*N1[space]*N1[space];
//                                BB[space] = N1[space];
//                            }
//
//                        }
//                        // 1,0,0,0 : j,i,k,l   :: 'f'2
//                        else  if ( leftChar == 'a' ){//(j,i,k,l)
//                            flagTranspose = 1; // a b | c
//                        }
//                        // 1,2,1,2 : j,i,l,k   :: 'g' 23
//                        else  if ( leftChar == 'g' ){//(j,i,l,k)
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space];
//                                B[space] = N1[space]*N1[space];
//                            }
//                            flagTranspose3 = 1;
//                            flagTranspose4 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space]*N1[space];
//                                BB[space] = N1[space]*N1[space];
//                            }
//
//                        }
//                        // 0,3,1,2 : j,k,i,l   :: 'h'19
//                        else  if ( leftChar == 'd' ){//(j,k,i,l)
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space]*N1[space];
//                                B[space] = N1[space];
//                            }
//
//                            flagTranspose3 = 1; // a b | c
//                            flagTranspose4 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space]*N1[space];
//                                BB[space] = N1[space]*N1[space];
//                            }
//
//                        }
//                        // 0,1,0,0 : j,k,l,i   :: 'i'3
//                        else  if ( leftChar == 'r' ){//(j,k,l,i)
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space];
//                                B[space] = N1[space]*N1[space]*N1[space];
//                            }
//                        }
//                        // 1,3,1,0 : j,l,i,k   :: 'j'13
//                        else  if ( leftChar == 'j' ){//(j,l,i,k)
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//                            flagTranspose3 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space]*N1[space];
//                                B[space] = N1[space];
//                            }
//                        }
//                        // 0,2,1,3 : j,l,k,i   :: 'k'21
//                        else  if ( leftChar == 't' ){//(j,l,k,i)12
//                            flagTranspose2 = 2;
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space];
//                                B[space] = N1[space]*N1[space];
//                            }
//                            flagTranspose3 = 1;
//                            flagTranspose4 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space]*N1[space]*N1[space];
//                                BB[space] = N1[space];
//                            }
//
//                        }
//                        // 0,2,1,1 : k,i,j,l   :: 'l'15
//                        else  if ( leftChar == 'c' ){//(k,i,j,l)
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space];
//                                B[space] = N1[space]*N1[space];
//                            }
//
//                            flagTranspose3 = 1;
//                            flagTranspose4 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space];
//                                BB[space] = N1[space]*N1[space]*N1[space];
//                            }
//
//                        }
//                        // 1,1,1,0 : k,i,l,j   :: 'm'10
//                        else  if ( leftChar == 'm' ){//(k,i,l,j)
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//                            flagTranspose3 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space];
//                                B[space] = N1[space]*N1[space]*N1[space];
//                            }
//                        }
//                        // 1,2,1,1 : k,j,i,l   :: 'n'16
//                        else  if ( leftChar == 'e' ){//(k,j,i,l)
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space];
//                                B[space] = N1[space]*N1[space];
//                            }
//
//                            flagTranspose3 = 1;
//                            flagTranspose4 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space];
//                                BB[space] = N1[space]*N1[space]*N1[space];
//                            }
//
//
//                        }
//                        // 0,1,1,0 : k,j,l,i   :: 'o'9
//                        else  if ( leftChar == 's' ){//(k,j,l,i)
//                            flagTranspose2 = 1;
//                            flagTranspose3 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space];
//                                B[space] = N1[space]*N1[space]*N1[space];
//                            }
//                        }
//                        // 0,2,0,0 : k,l,i,j   :: 'p'5
//                        else  if ( leftChar == 'p' ){//(k,l,i,j)
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space];
//                                B[space] = N1[space]*N1[space];
//                            }
//                        }
//                        // 1,1,0,1 : k,l,j,i   :: 'q'14
//                        else  if ( leftChar == 'v' ){//(k,l,j,i)//18
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space];
//                                B[space] = N1[space]*N1[space]*N1[space];
//                            }
//                            flagTranspose4 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space];
//                                BB[space] = N1[space]*N1[space]*N1[space];
//                            }
//
//                        }
//                        // 0,3,0,0 : l,i,j,k   :: 'r'7
//                        else  if ( leftChar == 'i' ){//(l,i,j,k)
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space]*N1[space];
//                                B[space] = N1[space];
//                            }
//                        }
//                        // 1,3,1,1 : l,i,k,j   :: 's'18
//                        else  if ( leftChar == 'o' ){//(l,i,k,j)
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space]*N1[space];
//                                B[space] = N1[space];
//                            }
//                            flagTranspose3 = 1; // a b | c
//                            flagTranspose4 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space];
//                                BB[space] = N1[space]*N1[space]*N1[space];
//                            }
//
//                        }
//                        // 1,3,0,0 : l,j,i,k   :: 't'8
//                        else  if ( leftChar == 'k' ){//(l,j,i,k)
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space]*N1[space];
//                                B[space] = N1[space];
//                            }
//                        }
//                        // 0,3,1,1 : l,j,k,i   :: 'u'17
//                        else  if ( leftChar == 'u' ){//(l,j,k,i)
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space]*N1[space];
//                                B[space] = N1[space];
//                            }
//
//                            flagTranspose3 = 1;
//                            flagTranspose4 = 1;
//                            for ( space = 0; space < SPACE ; space++){
//                                AA[space] = N1[space];
//                                BB[space] = N1[space]*N1[space]*N1[space];
//                            }
//                        }
//                        // 0,2,1,0 : l,k,i,j   :: 'v'6
//                        else  if ( leftChar == 'q' ){//(l,k,i,j)
//                            flagTranspose2 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space];
//                                B[space] = N1[space]*N1[space];
//                            }
//                            flagTranspose3 = 1;
//                        }
//                        // 1,2,1,0 : l,k,j,i   :: 'w'11
//                        else  if ( leftChar == 'w' ){//(l,k,j,i)//END
//                            flagTranspose = 1; // a b | c
//                            flagTranspose2 = 1;
//                            flagTranspose3 = 1;
//
//                            for ( space = 0; space < SPACE ; space++){
//                                A[space] = N1[space]*N1[space];
//                                B[space] = N1[space]*N1[space];
//                            }
//                        }
//                        else {
//                            printf("unknown flag %c\n",leftChar);
//                            exit(0);
//                        }
//
//                    }
                }
                
                
                
                if ( ( species(f1, tensorBuffers ) >= vector && bodies(f1,tensorBuffers)  >= four && part(f1, tensorBuffers)>=2) ){
                    for ( space =0 ; space < SPACE ; space++)
                        for ( r = 0; r < 2 ; r++)
                            array[space][r] = streams( f1, tensorBuffers, rank,space)+LN2[space]*r;
                }else {
                    printf("oopss\n");
                    exit(0);
                }
                
            }
        }
    }else
        if (  (species(f1, equals) == matrix && species(f1, left) == vector && species(f1, right) == vector && bodies(f1,left) == bodies(f1,right)&& bodies(f1, equals) == two))
        {//forms two-body density matrix.
            N1 = vectorLen(f1, equals);
            if ( bodies (f1,left ) == three ){
                NX[0] =f1->sinc.N1;
                if ( SPACE > 1 )
                NX[1] = f1->sinc.N1;
                if ( SPACE > 2 )
                NX[2] = f1->sinc.N1;
            } else if ( bodies(f1, left ) == four ){
                NX[0] = f1->sinc.N1*f1->sinc.N1;
                if ( SPACE > 1 )
                NX[1] = f1->sinc.N1*f1->sinc.N1;
                if ( SPACE > 2 )
                NX[2] = f1->sinc.N1*f1->sinc.N1;
            } else {
                if ( beta == -1 )
                    f1->sinc.tulip[equals].Current[espin] = 0;
                tOuterProductSu(f1, left, lspin, right, rspin, equals, espin);
                return 0.0;

            }
            type = 20;
            if ( rightChar != 'N' )
            {
                printf("fix your character! %d %d %d\n",equals, left , right);
                exit(0);
            }
            
        }else
            if (  (species(f1, equals) == matrix && species(f1, left) == vector && species(f1, right) == vector && bodies(f1,left) == bodies(f1,right)&& bodies(f1, equals) == one) )
            {//forms one-body density matrix.
                N1 = vectorLen(f1, equals);
                if ( bodies(f1, left ) == two ){
                    NX[0] = N1[0];
                    if ( SPACE > 1 )
                    NX[1] = N1[1];
                    if ( SPACE > 2 )

                    NX[2] = N1[2];
                }else if ( bodies (f1,left ) == three ){
                    NX[0] = N1[0]*N1[0];
                    if ( SPACE > 1 )
                    NX[1] = N1[1]*N1[1];
                    if ( SPACE > 2 )
                    NX[2] = N1[2]*N1[2];
                } else if ( bodies(f1, left ) == four ){
                    NX[0] = N1[0]*N1[0]*N1[0];
                    if ( SPACE > 1 )
                    NX[1] = N1[0]*N1[1]*N1[1];
                    if ( SPACE > 2 )
                    NX[2] = N1[0]*N1[2]*N1[2];
                } else  {
                    if ( beta == -1 )
                        f1->sinc.tulip[equals].Current[espin] = 0;
                    tOuterProductSu(f1, left, lspin, right, rspin, equals, espin);
                    return 0.0;
                }
                type = 20;
                if ( rightChar != 'N' )
                {
                    printf("fix your character! %d %d %d\n",equals, left , right);
                    exit(0);
                }
                
            }
    
    
            else
                if (  (species(f1, equals) == vector && species(f1, left) == matrix && species(f1, right) == vector && bodies(f1,left) == bodies(f1,right)))
                {
                    N1 = vectorLen(f1, left);
                    type = 1;
                    if ( rightChar != 'N' )
                    {
                        printf("fix your character! %d %d %d\n",equals, left , right);
                        exit(0);
                    }
                    
                } else if (  (species(f1, equals) == vector && species(f1, left)  ==vector && vector == species(f1, right)&& bodies(f1,left) ==two && one == bodies(f1,right)&& one == bodies(f1,equals))){
                    N1 = vectorLen(f1, right);
                    type = 1;
                }
                else if (  (species(f1, equals) == matrix && species(f1, left)  ==matrix && matrix == species(f1, right)&& bodies(f1,left) == bodies(f1,right))){
                    N1 = vectorLen(f1, left);
                    
                    type = 2;
                }
                else if (  (species(f1, equals) == matrix && species(f1, left)  ==quartic && matrix ==species(f1, right) && bodies(f1,left) == bodies(f1,right))){
                    N1 = vectorLen(f1, left);
                    if ( rightChar != 'N' )
                    {
                        printf("fix your character! %d %d %d\n",equals, left , right);
                        exit(0);
                    }
                    
                    type = 3;
                }else if (  (species(f1, equals) == quartic && species(f1, left)  ==quartic && quartic ==species(f1, right) && bodies(f1,left) == bodies(f1,right))){
                    N1 = vectorLen(f1, left);
                    
                    type = 4;
                } else if (  (species(f1, equals) == vector && species(f1, left)  ==matrix && vector == species(f1, right) && bodies(f1,left)  == one&&bodies(f1,right)  == two)){
                    N1 = vectorLen(f1, left);
                    //TWO BODIES
                    type = 5;
                    if ( rightChar != 'N' )
                    {
                        printf("fix your character! %d %d %d\n",equals, left , right);
                        exit(0);
                    }
                    
                    if ( f1->sinc.tulip[left].blockType == tv1 ){
                        skewSymbol = CblasNoTrans;
                        leftSymbol  = CblasNoTrans;
                    }else if ( f1->sinc.tulip[left].blockType == tv2 ){
                        skewSymbol = CblasTrans;
                        if ( leftChar == 'T')
                            leftSymbol = CblasNoTrans;
                        else if (leftChar == 'N')
                            leftSymbol = CblasTrans;
                    } else {
                        printf("type2 failure %d %d\n", left, f1->sinc.tulip[left].blockType);
                        exit(0);
                    }
                }
                else if (  (species(f1, equals) == vector && species(f1, left)  ==matrix && vector == species(f1, right) &&bodies(f1,left)  == one&&bodies(f1,right)  == four)){
                    N1 = n1;
                    if ( rightChar != 'N' )
                    {
                        printf("fix your character! %d %d %d\n",equals, left , right);
                        exit(0);
                    }
                    
                    if ( f1->sinc.tulip[left].blockType == tv1 ){
                        flagTranspose = 0;
                    }else  if ( f1->sinc.tulip[left].blockType == tv2 ){
                        flagTranspose = 1;
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space]*N1[space];
                        }
                    }else  if ( f1->sinc.tulip[left].blockType == tv3 ){
                        flagTranspose = 1;
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                    }else  if ( f1->sinc.tulip[left].blockType == tv4 ){
                        flagTranspose = 1;
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                    }else {
                        printf("type4 failure %d %d\n", left, f1->sinc.tulip[left].blockType);
                        exit(0);
                    }
                    
                    if ( bodies(f1, vectorCubeBuffers) != four)
                    {
                        printf("allocate more to vectorCubeBuffers\n");
                    }
                    
                    type = 6;
                }
                else if (  (species(f1, equals) == vector && species(f1, left)  ==matrix && vector == species(f1, right) &&bodies(f1,left)  == two&&bodies(f1,right)  == four)){
                    N1 = n1;
                    skewSymbol = CblasNoTrans;
                    if ( rightChar != 'N' )
                    {
                        printf("fix your character! %d %d %d\n",equals, left , right);
                        exit(0);
                    }
                    
                    if ( f1->sinc.tulip[left].blockType == e12 ){
                        flagTranspose = 0;
                        flagTranspose2 = 0;
                        
                    }else  if ( f1->sinc.tulip[left].blockType == e23 ){
                        flagTranspose = 1;
                        flagTranspose2 = 0;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space]*N1[space];
                        }
                    }else  if ( f1->sinc.tulip[left].blockType == e34 ){
                        flagTranspose = 1;
                        flagTranspose2 = 0;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                    }else  if ( f1->sinc.tulip[left].blockType == e14 ){
                        flagTranspose = 1;
                        flagTranspose2 = 0;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                    }else  if ( f1->sinc.tulip[left].blockType == e13 ){
                        flagTranspose = 1;
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space]*N1[space];
                        }
                    }else  if ( f1->sinc.tulip[left].blockType == e24 ){
                        flagTranspose = 1;
                        flagTranspose2 = 1;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                        skewSymbol = CblasTrans;
                        if ( leftChar == 'T')
                            leftSymbol = CblasNoTrans;
                        else if (leftChar == 'N')
                            leftSymbol = CblasTrans;
                        
                    }else {
                        printf("type6 failure %d %d\n", left, f1->sinc.tulip[left].blockType);
                        exit(0);
                    }
                    
                    if ( bodies(f1, vectorCubeBuffers) != four)
                    {
                        printf("allocate more to vectorCubeBuffers\n");
                        exit(0);
                    }
                    
                    type = 7;
                }
                else if (  (species(f1, equals) == vector &&
                            species(f1, left)  ==matrix &&
                            vector == species(f1, right) &&
                            bodies(f1,left)  == two  &&
                            bodies(f1,right)  == three)){
                    N1 = n1;
                    skewSymbol = CblasNoTrans;//not used.
                    leftSymbol = CblasNoTrans;
                    
                    if ( rightChar != 'N' )
                    {
                        printf("fix your character! %d %d %d\n",equals, left , right);
                        exit(0);
                    }
                    
                    if ( f1->sinc.tulip[left].blockType == e12 ){
                        flagTranspose = 0;
                        flagTranspose2 = 0;
                        
                    }else  if ( f1->sinc.tulip[left].blockType == e23 ){// 1| 23
                        flagTranspose = 1;
                        flagTranspose2 = 0;
                        
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                    }else  if ( f1->sinc.tulip[left].blockType == e13 ){// 12 | 3
                        if (0){
                            flagTranspose = 1;
                            flagTranspose2 = 1;
                            
                            for ( space = 0; space < SPACE ; space++){
                                A[space] = N1[space];
                                B[space] = N1[space]*N1[space];
                            }
                        }else {
                            if ( leftChar == 'T')
                                leftSymbol = CblasNoTrans;
                            else if (leftChar == 'N')
                                leftSymbol = CblasTrans;
                            
                            
                            flagTranspose = 1;
                            flagTranspose2 = 0;
                            for ( space = 0; space < SPACE ; space++){
                                A[space] = N1[space]*N1[space];
                                B[space] = N1[space];
                            }
                            
                        }
                    } else {
                        printf("type3 failure %d %d\n", left, f1->sinc.tulip[left].blockType);
                        exit(0);
                    }
                    
                    if ( bodies(f1, vectorCubeBuffers) != three)
                    {
                        printf("allocate more to vectorCubeBuffers\n");
                        exit(0);
                    }
                    
                    type = 8;
                }
                else if (  (species(f1, equals) == vector && species(f1, left)  ==matrix && vector == species(f1, right) &&bodies(f1,left)  == one&&bodies(f1,right)  == three)){
                    N1 = n1;
                    flagTranspose2 = 0;
                    leftSymbol = CblasNoTrans;
                    if ( rightChar != 'N' )
                    {
                        printf("fix your character! %d %d %d\n",equals, left , right);
                        exit(0);
                    }
                    
                    if ( f1->sinc.tulip[left].blockType == tv1 ){
                        flagTranspose = 0;
                    }else  if ( f1->sinc.tulip[left].blockType == tv2 ){
                        flagTranspose = 1;// a | b c
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space];
                            B[space] = N1[space]*N1[space];
                        }
                    }else  if ( f1->sinc.tulip[left].blockType == tv3 ){
                        flagTranspose = 1; // a b | c
                        for ( space = 0; space < SPACE ; space++){
                            A[space] = N1[space]*N1[space];
                            B[space] = N1[space];
                        }
                    }else {
                        printf("type3-3 failure %d %d\n", left, f1->sinc.tulip[left].blockType);
                        exit(0);
                    }
                    
                    if ( bodies(f1, vectorCubeBuffers) != three)
                    {
                        printf("allocate more to vectorCubeBuffers\n");
                    }
                    
                    type = 9;
                }else if (species(f1, equals) == matrix && bodies(f1, equals) == one && species(f1, left)  == vector && vector == species(f1, right)&& bodies(f1, left) == bodies(f1,right)){
                    N1 = n1;
                    if ( bodies(f1, left ) == one ){
                        tClear(f1, equals);
                        tOuterProductSu(f1, left, lspin, right, rspin, equals, espin);
                        return 0.;
                    }
                    if ( bodies(f1, left) == two)
                    {
                        for ( space = 0; space < SPACE ; space++)
                            A[space] = N1[space];
                    }
                    if ( bodies(f1, left) == three)
                    {
                        for ( space = 0; space < SPACE ; space++)
                            A[space] = N1[space]*N1[space];
                    }
                    if ( bodies(f1, left) == four)
                    {
                        for ( space = 0; space < SPACE ; space++)
                            A[space] = N1[space]*N1[space]*N1[space];
                    }
                    
                    type = 10;
                }
                else if (species(f1, equals) == matrix && bodies(f1, equals) == two && species(f1, left)  == vector && vector == species(f1, right)&& bodies(f1, left) == bodies(f1,right)){
                    N1 = vectorLen(f1,equals);
                    if ( bodies(f1, left )== one ){
                        printf("OU!!!\n");
                        exit(0);
                    }
                    if ( bodies(f1, left ) == two ){
                        tClear(f1, equals);
                        tOuterProductSu(f1, left, lspin, right, rspin, equals, espin);
                        return 0.;
                    }
                    if ( bodies(f1, left) == three)
                    {
                        for ( space = 0; space < SPACE ; space++)
                            A[space] = N1[space];
                        
                    }
                    if ( bodies(f1, left) == four)
                    {
                        for ( space = 0; space < SPACE ; space++)
                            A[space] = N1[space];//*N1[space];
                    }
                    
                    type = 10;
                }
    
    
    
                else if (  (species(f1, equals) == vector && species(f1, left)  ==vector && vector == species(f1, right)&& bodies(f1,left) == two&& two == bodies(f1,right)&& two == bodies(f1,equals))){
                    N1 = f1->sinc.Basis;
                    
                    type = 2;
                }
    
                else if (  (species(f1, equals) == vector && species(f1, left)  ==vector && vector == species(f1, right) && bodies(f1,left) == bodies(f1,right) + bodies(f1, equals))){
                    
                    //useful for forming higher dimensional inner products.
                    N1 = f1->sinc.Basis;
                    
                    if ( bodies(f1, right ) == one ){
                        NX[0] = N1[0];
                        if ( SPACE > 1 )

                        NX[1] = N1[1];
                        if ( SPACE > 2 )

                        NX[2] = N1[2];
                    } else
                        if ( bodies(f1, right ) == two ){
                            NX[0] = N1[0]*N1[0];
                            if ( SPACE > 1 )

                            NX[1] = N1[1]*N1[1];
                            if ( SPACE > 2 )

                            NX[2] = N1[2]*N1[2];
                        } else if ( bodies (f1, right ) == three ){
                            NX[0] = N1[0]*N1[0]*N1[0];
                            if ( SPACE > 1 )

                            NX[1] = N1[1]*N1[1]*N1[1];
                            if ( SPACE > 2 )

                            NX[2] = N1[2]*N1[2]*N1[2];
                        }else {
                            printf("high body count!\n");
                            exit(0);
                        }
                    N1 = vectorLen(f1,equals);
                    
                    type = 30;
                }
    
    if ( type == -1 )
    {
        
        printf("tM: wrong type %d=  %d * %d \n", equals, left, right);
        printf("Tm: %lld %lld %lld\n", species(f1, equals), species(f1, left), species(f1,right));
        printf("Tm: %d %d %d\n", bodies(f1, equals), bodies(f1, left), bodies(f1,right));
        
        exit(0);
    }
    if ( type == 0 ){
        double *pleft[SPACE],*pright[SPACE];
        
        sum = 0.;
        if ( flag == 1){
            
            {
                for ( space =0 ; space < spaces(f1,left) ; space++){
                    pleft[space] = streams( f1, left, lspin,space );
                    pright[space] = streams( f1, right, rspin,space );
                }
                for ( l = 0 ; l < LL ; l++){
                    for ( space = 0; space < spaces(f1,left) ;space++){//check they have same lengths...
                        {
                            transpose(N1[space],N1[space],pleft[space]+l*LN2[space],array[space][0]);
                           // mkl_domatcopy ('C', 'T',N1[space], N1[space],1.,pleft[space]+l*LN2[space],N1[space],array[space][0],N1[space]);
                            
                        }
                    }
                    for ( r = 0 ; r < LR ; r++){
                        product = 1.;
                        for ( space = 0; space < spaces(f1,left) ;space++){//check they have same lengths...
                            product *= cblas_ddot( LN2[space] ,array[space][0],1 , pright[space]+r*RN2[space], 1);
                        }
                        sum += product;//*co*coeff(f1,right,rspin,r);
                    }
                }
            }
            
            
            
            
        }else         if ( flag == 2){
            
            {
                for ( l = 0 ; l < LL ; l++){
                    for ( space =0 ; space < spaces(f1,left) ; space++){
                        pleft[space] = streams( f1, left, lspin,space )+l*LN2[space];
                        pright[space] = streams( f1, right, rspin,space );
                    }
                    //printf("%d %d %d\n", N1[0],LN2[0],RN2[0]);
                    for ( space = 0; space < spaces(f1,left) ;space++){//check they have same lengths...
                        transpose(N1[space],N1[space],pleft[space],array[space][0]);
                     //   printf("%lld mem %lld\n", pleft[space]-streams(f1,foundationStructure,0,space),LN2[space]);
                    }
                    //#pragma omp parallel for private (r, product, space) reduction(+:sum)
                    for ( r = 0 ; r < LR ; r++){
                        product = 1.;
                        for ( space = 0; space < spaces(f1,left) ;space++){//check they have same lengths...
                            //   printf("%lld %lld %lld %f - %f %f\n",LN2[space], pleft[space]-(streams( f1, left, lspin,space )+l*LN2[space]),streams( f1, right, rspin,space )-streams( f1, left, lspin,space ),cblas_ddot( LN2[space] ,pleft[space],1 , streams( f1, right, rspin,space )+r*RN2[space], 1),cblas_dnrm2( LN2[space] ,pleft[space],1),cblas_dnrm2( LN2[space] ,streams( f1, right, rspin,space )+r*RN2[space],1));
                            
                            product *= cblas_ddot( LN2[space] ,array[space][0],1 , pright[space]+r*RN2[space], 1);
                        }
                        sum += product;
                    }
                }
            }
            
            
            
            
        }else if (flag == 3){
            INT_TYPE bs,o;
            for ( l = 0 ; l < LL ; l++){
                
                for ( space = 0; space < spaces(f1,left) ;space++){
                    //
                    bs = 0;
                    pleft[space] = streams( f1, left, lspin,space )+l*LN2[space];
                    
                    
                    if (flagTranspose2){
                        for ( o = 0 ; o < N1[space] ;o++)
                            transpose(N1[space], N1[space],pleft[space]+N1[space]*N1[space]*o ,
                                      array[space][bs]+N1[space]*N1[space]*o);
                        
                          //  mkl_domatcopy ('C', 'T',N1[space], N1[space],
//                                           1.,
//                                           pleft[space]+N1[space]*N1[space]*o ,N1[space],
//                                           array[space][bs]+N1[space]*N1[space]*o,N1[space]);
                        pleft[space] = array[space][bs++];
                        bs = bs % 2;
                    }
                    if ( flagTranspose  ){
                        transpose(A[space], B[space],pleft[space] ,array[space][bs]);

                       // mkl_domatcopy ('C', 'T',A[space], B[space],
//                                       1.,
//                                       pleft[space] ,A[space],
//                                       array[space][bs],B[space]);
                        pleft[space] = array[space][bs++];
                        bs = bs % 2;
                        
                    }
                    
                    if ( flagTranspose3){
                        for ( o = 0 ; o < N1[space] ;o++)
                            transpose(N1[space], N1[space],pleft[space]+N1[space]*N1[space]*o ,
                                      array[space][bs]+N1[space]*N1[space]*o);

//                            mkl_domatcopy ('C', 'T',N1[space], N1[space],
//                                           1.,
//                                           pleft[space]+N1[space]*N1[space]*o ,N1[space],
//                                           array[space][bs]+N1[space]*N1[space]*o,N1[space]);
                        pleft[space] = array[space][bs++];
                        bs = bs % 2;
                        
                    }
                }
                
                
                
                for ( r = 0 ; r < LR ; r++){
                    product = 1.;
                    for ( space = 0; space < spaces(f1,left) ;space++){
                        product *= cblas_ddot( LN2[space] ,pleft[space],1 , streams( f1, right, rspin,space )+r*RN2[space], 1);
                    }
                    sum += product;
                    
                }
                
            }
            
            
            
        } else if (flag == 4){
            INT_TYPE bs,o;
            for ( l = 0 ; l < LL ; l++){
                
                for ( space = 0; space < spaces(f1,left) ;space++){
                    bs = 0;
                    //
                    pleft[space] = streams( f1, left, lspin,space )+l*LN2[space];
                    pright[space] = streams( f1, right, rspin,space );
                    
                    
                    if (flagTranspose){
                        for ( o = 0 ; o < N1[space]*N1[space] ;o++)
                            transpose(N1[space], N1[space],pleft[space]+N1[space]*N1[space]*o ,
                                      array[space][bs]+N1[space]*N1[space]*o);
                        pleft[space] = array[space][bs++];
                        bs = bs % 2;
                    }
                    
                    if ( flagTranspose2  ){
                        transpose(A[space], B[space],pleft[space] ,array[space][bs]);
                        pleft[space] = array[space][bs++];
                        bs = bs % 2;
                        
                        
                    }
                    
                    if ( flagTranspose3){
                        for ( o = 0 ; o < N1[space]*N1[space] ;o++)
                            transpose(N1[space], N1[space],pleft[space]+N1[space]*N1[space]*o ,
                                      array[space][bs]+N1[space]*N1[space]*o);

//                            mkl_domatcopy ('C', 'T',N1[space], N1[space],
//                                           1.,
//                                           pleft[space]+N1[space]*N1[space]*o ,N1[space],
//                                           array[space][bs]+N1[space]*N1[space]*o,N1[space]);
                        pleft[space] = array[space][bs++];
                        bs = bs % 2;
                        
                        
                    }
                    
                    if ( flagTranspose4  ){
                       
                        transpose(AA[space], BB[space],pleft[space] ,array[space][bs]);

                        pleft[space] = array[space][bs++];
                        bs = bs % 2;
                        
                        
                    }
                    
                }
                
                
                for ( r = 0 ; r < LR ; r++){
                    product = 1.;
                    
                    for ( space = 0; space < spaces(f1,left) ;space++){
                        //printf("%lld %lld %lld %f - %f %f\n",LN2[space], pleft[space]-(streams( f1, left, lspin,space )+l*LN2[space]),streams( f1, right, rspin,space )-streams( f1, left, lspin,space ),cblas_ddot( LN2[space] ,pleft[space],1 , streams( f1, right, rspin,space )+r*RN2[space], 1),cblas_dnrm2( LN2[space] ,pleft[space],1),cblas_dnrm2( LN2[space] ,streams( f1, right, rspin,space )+r*RN2[space],1));
                        
                        product *= cblas_ddot( LN2[space] ,pleft[space],1 , pright[space]+r*RN2[space], 1);
                    }
                    sum += product;
                    
                }
                
            }
            
        }
        else
            if ( 1 )
            {
                for ( space =0 ; space < spaces(f1,left) ; space++){
                    pleft[space] = streams( f1, left, lspin,space );
                    pright[space] = streams( f1, right, rspin,space );
                }
                //#pragma omp parallel for private (l,r, product, space) reduction(+:sum)
                for ( l = 0 ; l < LL ; l++){
                    for ( r = 0 ; r < LR ; r++){
                        product = 1.;
                        for ( space = 0; space < spaces(f1,left) ;space++){//check they have same lengths...
                            product *= cblas_ddot( LN2[space] , pleft[space]+l*LN2[space],1 , pright[space]+r*RN2[space], 1);
//                            INT_TYPE i ;
//                            if ( left ==lanesc)
//                            {
//                            for ( i = 0; i < LN2[space];i++)
//                                printf("%1.3f:",pleft[space][i]);
//                            printf("\n");
//                            }
                            }
                        sum += product;
//                        if ( left ==lanesc)
//                        printf("%d %d %f %d %u\n", l,r,product,LN2[0],pleft[0]);
                    }
                }
            }
    }
    else {
        
        for ( l = 0 ; l < LL ; l++)
            for ( r = 0 ; r < LR ; r++){
                ii = r*LL+l+beta*MM+gamma*CanonicalRank(f1,equals,espin);
                if ( type == 1 ){
                    //  N1
                    //=
                    //  N1 x N1
                    //  N1
                    {
                        for ( space = 0; space < spaces(f1,left) ;space++){
                            cblas_dgemv( CblasColMajor, leftSymbol,  N1[space], N1[space],1.,
                                        streams( f1, left, lspin,space )+l*LN2[space], N1[space],
                                        streams(f1, right, rspin,space)+r*RN2[space],1, 0.,
                                        streams( f1, equals, espin,space )+ii*EN2[space], 1  );
                            //  printf("%f\n", cblas_dnrm2(EN2[space], streams( f1, equals, espin,space )+ii*EN2[space], 1));
                            
                        }
                    }
                    
                }
                
                
                
                
                else if ( type == 2 ){
                    //  N1 x N1
                    //=
                    //  N1 x N1
                    //  N1 x N1
                    {
                        for ( space = 0; space < spaces(f1,left) ;space++){
                            cblas_dgemm(CblasColMajor, leftSymbol, rightSymbol,N1[space],N1[space],N1[space],1.,streams( f1, left, lspin,space )+l*LN2[space],N1[space],streams(f1, right, rspin,space)+r*RN2[space],N1[space], 0., streams( f1, equals, espin,space )+ii*EN2[space], N1[space]);
                        }
                    }
                }
                else if ( type == 20 ){
                    //  N1 x N1
                    //=
                    //  N1 x NX
                    //  NX x N1
                    {
                        for ( space = 0; space < spaces(f1,left) ;space++){
                            cblas_dgemm(CblasColMajor, leftSymbol, rightSymbol,N1[space],N1[space],NX[space],1.,streams( f1, left, lspin,space )+l*LN2[space],NX[space],streams(f1, right, rspin,space)+r*RN2[space],NX[space], 0., streams( f1, equals, espin,space )+ii*EN2[space], N1[space]);
                        }
                    }
                }
                
                else if ( type == 30 ){
                    //  N1 x N1
                    //=
                    //  N1 x NX
                    //  NX x N1
                    {
                        for ( space = 0; space < spaces(f1,left) ;space++){
                            if ( leftSymbol == CblasNoTrans )
                                cblas_dgemv(CblasColMajor, leftSymbol,N1[space],NX[space],1.,streams( f1, left, lspin,space )+l*LN2[space],N1[space],streams(f1, right, rspin,space)+r*RN2[space],1, 0., streams( f1, equals, espin,space )+ii*EN2[space], 1);
                            else
                                cblas_dgemv(CblasColMajor, leftSymbol,NX[space],N1[space],1.,streams( f1, left, lspin,space )+l*LN2[space],NX[space],streams(f1, right, rspin,space)+r*RN2[space],1, 0., streams( f1, equals, espin,space )+ii*EN2[space], 1);
                            
                        }
                    }
                }
                
                else if ( type == 3 ){
                    //  N2
                    //=
                    //  N2 x N2
                    //  N2
                    {
                        for ( space = 0; space < spaces(f1,left) ;space++){
                            cblas_dgemv( CblasColMajor, CblasNoTrans, N1[space]*N1[space], N1[space]*N1[space],1., streams( f1, left, lspin,space )+ l*LN2[space] , N1[space]*N1[space], streams(f1, right, rspin,space)+r*RN2[space],1, 0., streams( f1, equals, espin,space )+ii*EN2[space], 1  );
                        }
                    }
                }
                else if ( type == 4 ){
                    //  N2 x N2
                    //=
                    //  N2 x N2
                    //  N2 x N2
                    {
                        for ( space = 0; space < spaces(f1,left) ;space++){
                            cblas_dgemm(CblasColMajor, leftSymbol, rightSymbol,N1[space]*N1[space],N1[space]*N1[space],N1[space]*N1[space],1.,streams( f1, left, lspin,space )+l*LN2[space],N1[space]*N1[space],streams(f1, right, rspin,space)+r*RN2[space],N1[space]*N1[space], 0., streams( f1, equals, espin,space )+ii*EN2[space], N1[space]*N1[space]);
                        }
                    }
                }
                else if ( type == 5 ){
                    //TWO_BODY VECTOR , ONE-BODY MATRIX
                    //  N1 x N1
                    //=
                    //  N1 x N1
                    //  N1 x N1
                    //  1 body vector operations
                    {
                        for ( space = 0; space < spaces(f1,left) ;space++){
                            if (skewSymbol == CblasNoTrans )
                                cblas_dgemm(CblasColMajor,leftSymbol, CblasNoTrans,N1[space],N1[space],N1[space],1.,streams( f1, left, lspin,space )+l*LN2[space],N1[space],streams(f1, right, rspin,space)+r*RN2[space],N1[space], 0., streams( f1, equals, espin,space )+ii*EN2[space], N1[space]);
                            else
                                //special case,  general solution is in type == 6
                                cblas_dgemm(CblasColMajor, CblasNoTrans, leftSymbol,N1[space],N1[space],N1[space],1.,streams(f1, right, rspin,space)+r*RN2[space],N1[space],streams( f1, left, lspin,space )+l*LN2[space],N1[space], 0., streams( f1, equals, espin,space )+ii*EN2[space], N1[space]);
                            
                        }
                    }
                }
                else if ( type == 6 ){
                    //FOUR BODY vector ,  ONE-BODY matrix
                    
                    for ( space = 0; space < spaces(f1, left) ;space++){
                        
                        if ( flagTranspose ){
                            transpose(A[space], B[space],streams(f1, right, rspin,space)+r*RN2[space] ,streams(f1,vectorCubeBuffers,rank,0 ));

                            
                            cblas_dgemm(CblasColMajor, leftSymbol,CblasNoTrans,
                                        N1[space],N1[space]*N1[space]*N1[space],N1[space],
                                        1.,
                                        streams(f1, left, lspin,space )+l*LN2[space],N1[space],
                                        streams(f1,vectorCubeBuffers,rank,0 ),N1[space],
                                        0.,
                                        streams(f1,vectorCubeBuffers,rank,1 ), N1[space]);
                            
                            transpose(B[space], A[space],streams(f1,vectorCubeBuffers,rank,1 ),streams(f1, equals, espin,space )+ii*RN2[space]);
                        } else {
                            
                            cblas_dgemm(CblasColMajor, leftSymbol, CblasNoTrans,N1[space],N1[space]*N1[space]*N1[space],N1[space],
                                        1.,
                                        streams( f1, left, lspin,space)+l*LN2[space],N1[space],
                                        streams(f1, right, rspin,space)+r*RN2[space],N1[space],
                                        0.,
                                        streams( f1, equals, espin,space )+ii*EN2[space], N1[space]);
                            
                        }
                        
                    }
                }
                else if ( type == 7 ){
                    //FOUR BODY vector ,  TWO-BODY matrix
                    INT_TYPE o;
                    for ( space = 0; space < spaces(f1, left) ;space++){
                        
                        if ( flagTranspose && flagTranspose2 ){
                            for ( o = 0; o < N1[space]*N1[space] ; o++)
                                transpose(N1[space], N1[space],streams(f1, right, rspin,space)+r*RN2[space]+o*N1[space]*N1[space],
                                          streams(f1,vectorCubeBuffers,rank,0 )+o*N1[space]*N1[space]);
                            
                            transpose(A[space], B[space], streams(f1,vectorCubeBuffers,rank,0 ), streams(f1,vectorCubeBuffers,rank,1 ));
                            if ( skewSymbol == CblasNoTrans )
                                cblas_dgemm(CblasColMajor, leftSymbol, CblasNoTrans,N1[space]*N1[space],N1[space]*N1[space],N1[space]*N1[space],
                                            1.,
                                            streams( f1, left, lspin,space )+l*LN2[space],N1[space]*N1[space],
                                            streams(f1,vectorCubeBuffers,rank,1 ),N1[space]*N1[space],
                                            0.,
                                            streams(f1,vectorCubeBuffers,rank,2 ),N1[space]*N1[space]);
                            else
                                cblas_dgemm(CblasColMajor,CblasNoTrans, leftSymbol, N1[space]*N1[space],N1[space]*N1[space],N1[space]*N1[space],
                                            1.,
                                            streams(f1,vectorCubeBuffers,rank,1 ),N1[space]*N1[space],
                                            streams( f1, left, lspin,space )+l*LN2[space],N1[space]*N1[space],
                                            0.,
                                            streams(f1,vectorCubeBuffers,rank,2 ),N1[space]*N1[space]);
                            
                            
                            
                            
                            transpose(B[space], A[space],streams(f1,vectorCubeBuffers,rank,2 ),streams(f1,vectorCubeBuffers,rank,1 ));
//                            mkl_domatcopy ('C', 'T',B[space], A[space],
//                                           1.,
//                                           streams(f1,vectorCubeBuffers,rank,2 ),B[space],
//                                           streams(f1,vectorCubeBuffers,rank,1 ),A[space]);
                            
                            for ( o = 0; o < N1[space]*N1[space] ; o++)
                                transpose(N1[space], N1[space],streams(f1,vectorCubeBuffers,rank,1 )+o*N1[space]*N1[space],streams(f1, equals, espin,space )+ii*EN2[space]+o*N1[space]*N1[space]);
                                
//                                mkl_domatcopy ('C', 'T',N1[space], N1[space],
//                                               1.,
//                                               streams(f1,vectorCubeBuffers,rank,1 )+o*N1[space]*N1[space],N1[space],
//                                               streams(f1, equals, espin,space )+ii*EN2[space]+o*N1[space]*N1[space],N1[space]);
                            
                        } else
                            if ( flagTranspose && ! flagTranspose2 ){
                                transpose(A[space], B[space],streams(f1, right, rspin,space)+r*RN2[space] ,streams(f1,vectorCubeBuffers,rank,0 ));
                                
//                                mkl_domatcopy ('C', 'T',A[space], B[space],
//                                               1.,
//                                               streams(f1, right, rspin,space)+r*RN2[space] ,A[space],
//                                               streams(f1,vectorCubeBuffers,rank,0 ),B[space]);
                                
                                cblas_dgemm(CblasColMajor, leftSymbol, CblasNoTrans,N1[space]*N1[space],N1[space]*N1[space],N1[space]*N1[space],
                                            1.,
                                            streams( f1, left, lspin,space )+l*LN2[space],N1[space]*N1[space],
                                            streams(f1,vectorCubeBuffers,rank,0 ),N1[space]*N1[space],
                                            0.,
                                            streams(f1,vectorCubeBuffers,rank,1 ),N1[space]*N1[space]);
                                
                                transpose(B[space], A[space],streams(f1,vectorCubeBuffers,rank,1 ),streams(f1, equals, espin,space )+ii*EN2[space]);
//                                mkl_domatcopy ('C', 'T',B[space], A[space],
//                                               1.,
//                                               streams(f1,vectorCubeBuffers,rank,1 ),B[space],
//                                               streams(f1, equals, espin,space )+ii*EN2[space],A[space]);
                            } else if ( (! flagTranspose) && (! flagTranspose2 ) ){
                                
                                cblas_dgemm(CblasColMajor, leftSymbol, CblasNoTrans,N1[space]*N1[space],N1[space]*N1[space],N1[space]*N1[space],
                                            1.,
                                            streams( f1, left, lspin,space )+l*LN2[space],N1[space]*N1[space],
                                            streams(f1, right, rspin,space)+r*RN2[space],N1[space]*N1[space],
                                            0.,
                                            streams( f1, equals, espin,space )+ii*EN2[space], N1[space]*N1[space]);
                                
                            } else {
                                printf("oops\n");
                                exit(0);
                            }
                        
                    }
                }
                else if ( type == 8 ){
                    //THREE BODY vector ,  TWO-BODY matrix
                    INT_TYPE o;
                    for ( space = 0; space < spaces(f1, left) ;space++){
                        
                        //(13)
                        if ( flagTranspose && flagTranspose2 ){
                            
                            
                            for ( o = 0 ; o < N1[space] ;o++)
                                transpose(N1[space], N1[space],streams(f1, right, rspin,space)+r*RN2[space]+N1[space]*N1[space]*o,streams(f1,vectorCubeBuffers,rank,0 )+N1[space]*N1[space]*o);
//                                mkl_domatcopy ('C', 'T',N1[space], N1[space],
//                                               1.,
//                                               streams(f1, right, rspin,space)+r*RN2[space]+N1[space]*N1[space]*o ,N1[space],
//                                               streams(f1,vectorCubeBuffers,rank,0 )+N1[space]*N1[space]*o,N1[space]);
                            //213
                            
                            transpose(A[space], B[space],streams(f1,vectorCubeBuffers,rank,0 ),streams(f1,vectorCubeBuffers,rank,1 ));
//                            mkl_domatcopy ('C', 'T',A[space], B[space],
//                                           1.,
//                                           streams(f1,vectorCubeBuffers,rank,0 ),A[space],
//                                           streams(f1,vectorCubeBuffers,rank,1 ),B[space]);
                            //
                            
                            //132
                            
                            
                            //                            if ( skewSymbol == CblasNoTrans )
                            cblas_dgemm(CblasColMajor, leftSymbol, CblasNoTrans,N1[space]*N1[space],N1[space],N1[space]*N1[space],
                                        1.,
                                        streams( f1, left, lspin,space )+l*LN2[space],N1[space]*N1[space],
                                        streams(f1,vectorCubeBuffers,rank,1 ),N1[space]*N1[space],
                                        0.,
                                        streams(f1,vectorCubeBuffers,rank,2 ),N1[space]*N1[space]);
                            // (13)2
                            //                            else
                            //                                cblas_dgemm(CblasColMajor,CblasNoTrans, leftSymbol, N1[space]*N1[space],N1[space],N1[space]*N1[space],
                            //                                            1.,
                            //                                            streams(f1,vectorCubeBuffers,rank,1 ),N1[space]*N1[space],
                            //                                            streams( f1, left, lspin,space )+l*LN2[space],N1[space]*N1[space],
                            //                                            0.,
                            //                                            streams(f1,vectorCubeBuffers,rank,2 ),N1[space]*N1[space]);
                            
                            transpose(B[space], A[space],streams(f1,vectorCubeBuffers,rank,2 ),streams(f1,vectorCubeBuffers,rank,1 ));
                            // 2(13)
                            
                            
                            for ( o = 0 ; o < N1[space] ;o++)
                                transpose(N1[space], N1[space],streams(f1,vectorCubeBuffers,rank,1 )+N1[space]*N1[space]*o,streams(f1, equals, espin,space )+ii*EN2[space]+N1[space]*N1[space]*o);
                            // (1)2(3)
                        } else
                            
                            //(23)
                            if ( flagTranspose && ! flagTranspose2 ){
                                transpose(A[space], B[space], streams(f1, right, rspin,space)+r*RN2[space], streams(f1,vectorCubeBuffers,rank,0 ));
                                
                                cblas_dgemm(CblasColMajor, leftSymbol, CblasNoTrans,N1[space]*N1[space],N1[space],N1[space]*N1[space],
                                            1.,
                                            streams( f1, left, lspin,space )+l*LN2[space],N1[space]*N1[space],
                                            streams(f1,vectorCubeBuffers,rank,0 ),N1[space]*N1[space],
                                            0.,
                                            streams(f1,vectorCubeBuffers,rank,1 ), N1[space]*N1[space]);
                                
                                transpose(B[space], A[space], streams(f1,vectorCubeBuffers,rank,1 ), streams(f1, equals, espin,space )+ii*EN2[space]);
                            } else
                                //(13)
                                if ( (! flagTranspose) && (! flagTranspose2 ) ){
                                    
                                    cblas_dgemm(CblasColMajor, leftSymbol, CblasNoTrans,N1[space]*N1[space],N1[space],N1[space]*N1[space],
                                                1.,
                                                streams( f1, left, lspin,space )+l*LN2[space],N1[space]*N1[space],
                                                streams(f1, right, rspin,space)+r*RN2[space],N1[space]*N1[space],
                                                0.,
                                                streams( f1, equals, espin,space )+ii*EN2[space], N1[space]*N1[space]);
                                    
                                } else {
                                    printf("oops\n");
                                    exit(0);
                                }
                        
                    }
                }
                else if ( type == 9 ){
                    //THREE BODY vector ,  ONE-BODY matrix
                    
                    for ( space = 0; space < spaces(f1, left) ;space++){
                        
                        if ( flagTranspose ){
                            //if A is 3D -> (2)
                            //if A is 6D -> (3)
                            transpose( A[space], B[space], streams(f1, right, rspin,space)+r*RN2[space], streams(f1,vectorCubeBuffers,rank,0 ));
//                            mkl_domatcopy ('C', 'T',A[space], B[space],
//                                           1.,
//                                           streams(f1, right, rspin,space)+r*RN2[space] ,A[space],
//                                           streams(f1,vectorCubeBuffers,rank,0 ),B[space]);
                            
                            cblas_dgemm(CblasColMajor, leftSymbol,CblasNoTrans,
                                        N1[space],N1[space]*N1[space],N1[space],
                                        1.,
                                        streams(f1, left, lspin,space )+l*LN2[space],N1[space],
                                        streams(f1,vectorCubeBuffers,rank,0 ),N1[space],
                                        0.,
                                        streams(f1,vectorCubeBuffers,rank,1 ), N1[space]);
                            
                            transpose(B[space], A[space], streams(f1,vectorCubeBuffers,rank,1 ), streams(f1, equals, espin,space )+ii*EN2[space]);
                        } else {
                            //(1)
                            cblas_dgemm(CblasColMajor, leftSymbol, CblasNoTrans,N1[space],N1[space]*N1[space],N1[space],
                                        1.,
                                        streams( f1, left, lspin,space)+l*LN2[space],N1[space],
                                        streams(f1, right, rspin,space)+r*RN2[space],N1[space],
                                        0.,
                                        streams( f1, equals, espin,space )+ii*EN2[space], N1[space]);
                            
                        }
                        
                    }
                }
                else if ( type == 10 ){
                    //  N1 x N1
                    //=
                    //  N1 x NA
                    //  N1 x NA ^T
                    {
                        for ( space = 0; space < spaces(f1,left) ;space++){
                            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,N1[space],N1[space],A[space],1.,streams( f1, left, lspin,space )+l*LN2[space],N1[space],streams(f1, right, rspin,space)+r*RN2[space],N1[space], 0., streams( f1, equals, espin,space )+ii*EN2[space], N1[space]);
                        }
                    }
                }
            }
    }
    if ( species(f1,equals) != scalar  ){
        f1->sinc.tulip[equals].Current[espin] = MM+beta*MM+gamma*CanonicalRank(f1,equals,espin);
    }
    if ( isnan( sum ) || isinf(sum)){
#if VERBOSE
        printf("no!%lld %lld %lld %f %lld %lld %lld %lld %d %d %d %lld %lld %lld %c %c\n",reduce, MM, gamma,number, LL,LR, LN2[0], RN2[0], equals,left, right,rank,lspin,rspin,leftChar,rightChar);
#endif
        *info = 1;
        return 0.;
    }
    
    return sum;
}

// equals = product * left . right + shift * right
void tHXpX (  INT_TYPE rank, struct field * f1 , enum division left,INT_TYPE shiftFlag, double product, double productCmpl,  enum division right ,  double tolerance , INT_TYPE maxRun  ){
    
    if ( right == copyVector || right == copyTwoVector|| right == copyThreeVector|| right == copyFourVector || right == totalVector  ){
        printf("conflict\n");
        exit(0);
    }
    if ( right == diagonalVector ){
        tEqua(f1, copyThreeVector,rank, right,rank );
        f1->sinc.tulip[copyFourVector].Current[rank] = 0;
    }else {
        tEqua(f1, copyThreeVector,rank, right,0 );
        tEqua(f1, copyFourVector,rank, right,1 );
    }
    INT_TYPE info,spin,cmpl,cmpl2,skipFlag,alist=0;
    time_t start_t;
    enum division Mat,Vec;
    if ( maxRun < 1 )
        return;
    if ( !CanonicalRank(f1, copyThreeVector,rank ) &&! CanonicalRank(f1, copyFourVector,rank ) )
        return;
    
    for ( spin = 0 ; spin < 1 + (right != diagonalVector) ; spin++){
        
        f1->sinc.tulip[totalVector].header = header(f1, right);
        f1->sinc.tulip[totalVector].Current[rank] = 0;
        time(&start_t);
        
        if ( shiftFlag ){
            if ( spin ==0 ){ //Y
                tAddTw(f1, totalVector, rank,copyThreeVector,rank);
            }else
                if ( spin == 1 ){
                    tAddTw(f1, totalVector, rank,copyFourVector,rank);
                }
        }
        for ( cmpl = 0; cmpl < 2 ; cmpl++)//HX
            for ( cmpl2 = 0; cmpl2 < 2 ; cmpl2++)
                if (
                    (((!spin && cmpl == cmpl2) || (spin && cmpl != cmpl2 )) && (fabs(product) > 1e-6)  )
                    ||
                    (((!spin && cmpl != cmpl2) || (spin && cmpl == cmpl2 )) && (fabs(productCmpl) > 1e-6) )
                    )
                {
                    enum division leftP = left;
                    alist = 0;
                    skipFlag = 0;
                    do {
                        if ( linear == name(f1,leftP) ){
                            Mat = rivers(rank, f1, linear, alist );
                            f1->sinc.tulip[Mat].blockType = f1->sinc.tulip[leftP].blockType;
                        }else {
                            Mat = leftP;
                        }

                    
                        if ( cmpl2 == 0 ){
                            Vec = copyThreeVector;
                        }else if ( cmpl2 == 1 ){
                            Vec =  copyFourVector;
                        }else {
                            printf("cmpl2 %d\n", cmpl2);
                            exit(0);
                        }
#if VERBOSE
                        printf("%d-%d-%d-%d\t%d\t%d\t%d-%d-%d\t %d %f\n",Mat,name(f1,Mat),cmpl,cmpl2, bodies(f1, Mat),f1->sinc.tulip[Mat].blockType,CanonicalRank(f1, name(f1,Mat), cmpl),CanonicalRank(f1, Mat, cmpl),f1->sinc.tulip[Mat].ptRank[cmpl], alist,traceOne(f1, Mat, cmpl) );
                        fflush(stdout);
#endif

                        if ( CanonicalRank(f1,Mat, cmpl)&& CanonicalRank(f1, name(f1,Mat), cmpl)&&CanonicalRank(f1, Vec, rank) ){
                            tEqua(f1, copyTwoVector, rank, Vec, rank);
                            tCycleMultiplyMP(rank, f1, 'N',Mat, cmpl,Vec,rank, copyTwoVector, rank, f1->mem1->rt->vCANON, maxRun, -1.);
                            {
                                if(((!spin && cmpl == cmpl2) || (spin && cmpl != cmpl2 )) && (fabs(product) > 1e-6)  ){
                                    if ( cmpl == 1 && cmpl2 == 1 && spin == 0 ){
                                        tEqua(f1, copyVector, rank,copyTwoVector,rank);
                                        tScaleOne(f1, copyVector,rank, -product);
                                        tAddTw(f1, totalVector,rank, copyVector,rank);
                                    } else {
                                        tEqua(f1, copyVector, rank,copyTwoVector,rank);
                                        tScaleOne(f1, copyVector,rank, product);
                                        tAddTw(f1, totalVector,rank, copyVector,rank);

                                    }
                                }
                                if(((!spin && cmpl != cmpl2) || (spin && cmpl == cmpl2 )) && (fabs(productCmpl) > 1e-6) )
                                {
                                    if (( cmpl == 0 && cmpl2 == 0 && spin == 1 )){
                                        tEqua(f1, copyVector, rank,copyTwoVector,rank);
                                        tScaleOne(f1, copyVector,rank, productCmpl);
                                        tAddTw(f1, totalVector,rank, copyVector,rank);
                                        
                                    } else {
                                        tEqua(f1, copyVector, rank,copyTwoVector,rank);
                                        tScaleOne(f1, copyVector,rank, -productCmpl);
                                        tAddTw(f1, totalVector,rank, copyVector,rank);
                                    }
                                }
#if VERBOSE
                                printf(" : %d %d: %d %d = o%d %f %d (%d)/n",Mat,CanonicalRank(f1, name(f1,Mat),cmpl),cmpl,cmpl2,spin,product,CanonicalRank(f1,right, cmpl2), alist);
                                fflush(stdout);
#endif

                            }
                        }
                        if (name(f1, leftP) == linear ){
                            alist++;
                            
                            if ( alist > f1->Na ){
                                leftP = f1->sinc.tulip[leftP].linkNext;
                                alist = 0;
                            }
                            
                        } else {
                            leftP = f1->sinc.tulip[leftP].linkNext;
                        }

#if VERBOSE
                        if (CanonicalRank(f1, totalVector, rank) ){
                        printf("[%f:%d]\n", inner(rank, f1, totalVector, rank),CanonicalRank(f1, totalVector, rank));
                        }
#endif

                    } while ( leftP != nullName);
                    
                    //PRODUCT
                    //spin//cmpl//cmpl2//
                    //0     0       0                           1
                    //0     1       1                          -1
                    //1     1       0                           1
                    //1     0       1                           1
                    //PRODUCTCMPL
                    //(iP) (cmpl)( cmpl2)
                    
                    //spin//cmpl//cmpl2//
                    //0     1       0                           -1
                    //0     0       1                           -1
                    //1     0       0                            1
                    //1     1       1                           -1
                    
                    //SHIFT
                    //spin////cmpl2//
                    //0            0                               1
                    //1            1                               1
                    //SHIFTCMPL
                    //spin///cmpl2//
                    //0            1                              -1
                    //1            0                               1
                    
                
            }
        if ( CanonicalRank(f1, totalVector, rank)){
            
#if VERBOSE
           printf("%f____", tMultiplyMP(rank, &info, f1, 1.0, -1, nullVector, 0, 'T', totalVector, rank, 'N', totalVector, rank));
#endif
            tMultiplyMP(rank, &info, f1, 1.0, -1, nullVector, 0, 'T', totalVector, rank, 'N', totalVector, rank);
            if ( !info ){
                if ( right == diagonalVector )
                    tCycleDecompostionOneMP(rank,f1, totalVector, rank, right,rank, f1->mem1->rt->vCONVERGENCE, maxRun, -1.);
                else
                    tCycleDecompostionOneMP(rank,f1, totalVector, rank, right,spin, f1->mem1->rt->vCONVERGENCE, maxRun, -1.);
            }
        }
    }
    
#if 0
    time(&lapse_t);
    if ( maxRun > 0 )
        f1->mem1->rt->timeGEMV[rank][maxRun-1] += difftime(lapse_t, start_t);
    time(&start_t);
    if ( maxRun > 0 )
        f1->mem1->rt->distBeylkin[rank][maxRun-1] += distanceOne(rank, f1, totalVector, rank, equals, spin);
    if ( maxRun > 0 )
        f1->mem1->rt->count[rank][maxRun-1]++;
    time(&lapse_t);
    if ( maxRun > 0 )
        f1->mem1->rt->timeBeylkin[rank][maxRun-1] += difftime(lapse_t, start_t);
#endif
    return ;
}

double distanceOne(INT_TYPE rank,struct field * f1 , enum division alloy , INT_TYPE spin , enum division alloyBak, INT_TYPE spin2){
    double value,value2;
    INT_TYPE info;
    value = 0.5*( tMultiplyMP(rank,&info,f1,1., -1,nullMatrix, 0 ,'T', alloy, spin,'N', alloy ,spin)+tMultiplyMP(rank,&info,f1,1., -1,nullMatrix, 0 ,'T', alloyBak, spin2,'N', alloyBak ,spin2));
    value2 = 2*fabs(value-tMultiplyMP(rank,&info,f1,1., -1,nullMatrix, 0 , 'T', alloyBak, spin2,'N', alloy ,spin) );
    return value2;
    
}

double inner(INT_TYPE rank,struct field * f1 , enum division alloy , INT_TYPE spin ){
    INT_TYPE info;
    return tMultiplyMP(rank,&info,f1,1., -1,nullMatrix, 0 ,'T', alloy, spin,'N', alloy ,spin);;
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


INT_TYPE tConstructDensity(struct calculation * calc , INT_TYPE ct ){
    INT_TYPE i,info,rank=0;
    if ( part(&calc->i.c,squareTwo) < ct*part(&calc->i.c,copyThree))
    {
        printf("consider more ranks on squareTwo\n");
        exit(0);
    }
    for ( i = 0 ; i < ct ; i++){
        calc->i.c.sinc.tulip[copyTwo].Current[0] = 0;
        tMultiplyMP(rank, &info,&calc->i.c, 1.0, -2, copyTwo, 0, 'T', eigenVectors+i, 1, 'N' ,eigenVectors+i , 1);
        tMultiplyMP(rank, &info,&calc->i.c, 1.0, -2, copyTwo, 0, 'T', eigenVectors+i, 0, 'N' ,eigenVectors+i , 0);
        tCycleDecompostionOneMP(rank, &calc->i.c, copyTwo, 0, copyThree, 0,calc->rt.CANON, part(&calc->i.c, copyThree), -1);
        printf("trace-%d %f\n",i+1,traceOne(&calc->i.c,copyTwo,0));
        fflush(stdout);
        tAddTw(&calc->i.c,squareTwo, 0,copyThree,0);
    }
    tCycleDecompostionOneMP(rank, &calc->i.c, squareTwo, 0, density, 0,calc->rt.CANON, part(&calc->i.c, density), -1);
    tScale(&calc->i.c,density,1./traceOne(&calc->i.c,density,0));
    
    return ct;
}


double tLanczosConvergence (struct field * f1 ,enum division Ha, enum division Vec, enum division usr){
    INT_TYPE sumR = 0,info;
    enum division Mat = Ha;
    do {
        sumR += CanonicalRank(f1, Mat ,0) + CanonicalRank(f1, Mat ,1);
        Mat = f1->sinc.tulip[Mat].linkNext;
    } while ( Mat != nullName );
    
    if ( sumR * CanonicalRank(f1, Vec,0) > part(f1, canonicalBuffersC)){
        printf("consider making a bigger pot!\n");
        return -1.;
    }
    
    enum division buffer = ocean(0, f1, eigenVectors, 0, 0);
    f1->sinc.tulip[buffer].name = usr;
    tClear(f1,buffer);
    double sum = 0.,sum2=0.;
    tScale(f1, Vec, 1./magnitude(f1, Vec));
    
   // for ( cmpl = 0; cmpl < 2 ; cmpl++)
    {
        Mat = Ha;
        do {
//            for ( cmpl2 = 0; cmpl2 < 2 ; cmpl2++){
//                if ( cmpl == 0 )
//                    cmpl3 = cmpl2;
//                else
//                    cmpl3 = ! cmpl2;
                tMultiplyMP(0, &info, f1, 1., -2, buffer, 0, 'N', Mat, 0, 'N', Vec,0 );
            
          //  printf("buffers %d %d \n",Mat, CanonicalRank(f1, buffer, 0));
            Mat = f1->sinc.tulip[Mat].linkNext;
        } while ( Mat != nullName );
        sum2 += inner(0, f1, buffer, 0);
        sum += sqr(tMultiplyMP(0, &info, f1, 1., -1, nullVector, 0, 'T', buffer, 0, 'N', Vec,0));
    }
   // printf("%f\t %f\n", sum2 , sum);
    return sum2 - sum;
}



INT_TYPE xConstructFoundation (struct calculation * calc , enum division usr, INT_TYPE UR, struct calculation * calc2, enum division usz, INT_TYPE UZ ,INT_TYPE mx){
    
   // printf("%d %d %d\n", UR,UZ,mx);
    INT_TYPE mdi,ii,iii,iv,rank = 0,cmpl,i;
    
    
    //    for ( i = 0 ; i < UR ; i++)
    //        sumR +=  part(&calc->i.c, usr+i);
    //
    //    for ( i = 0 ; i < UZ ; i++)
    //        sumEr +=  CanonicalRank(&calc2->i.c,usz+i, 0);
    //    for ( i = 0 ; i < UZ ; i++)
    //        sumEc +=  CanonicalRank(&calc2->i.c,usz+i, 1);
    //
    //    for ( i = 0; i < imin(UZ,UR);i++)
    //        if ( part(&calc->i.c, usr+i ) < CanonicalRank(&calc2->i.c, usz+i, 0) || CanonicalRank(&calc->i.c, usz+i, 0) > mx)
    //            flag++;
    
    enum division f[UZ];
    ii = 0;
    iii= 0;
    iv = 0;
    for ( i = 0; i < UZ ; i++)
        if ( CanonicalRank(&calc2->i.c, usz+i, 0)+CanonicalRank(&calc2->i.c, usz+i, 1))
        if ( calc2->i.c.sinc.tulip[usz+i].symmetry == calc2->i.irrep|| ! calc2->i.irrep){
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
       // printf ("oopsy %d %d\n",iii,UR);
       // exit(0);
    }
    
//#ifdef OMP
//#pragma omp parallel for private (rank,i,mdi,cmpl) schedule(dynamic,1)
//#endif
    for ( i = 0; i < ii ; i++){
//#ifdef OMP
//        rank = omp_get_thread_num();
//#else
//        rank = 0;
//#endif
        cmpl = 0;
        ///            printf("%d %d \n", cmpl,i);
        //            fflush(stdout);
       // for ( cmpl = 0; cmpl < spins(f1,usz) ; cmpl++)
        {
            zero(&calc->i.c, copyVector,rank);
            calc->i.c.sinc.tulip[copyVector].Current[rank] = 0;
            if ( bodies(&calc2->i.c, usz+f[i]) == two )
                xTwoBand(&calc2->i.c, usz+f[i], cmpl, &calc->i.c, copyVector,rank,calc->rt.runFlag );
            else if ( bodies(&calc2->i.c, usz+f[i]) == three )
                xThreeBand(&calc2->i.c, usz+f[i], cmpl, &calc->i.c, copyVector,rank,calc->rt.runFlag );
            else if ( bodies(&calc2->i.c, usz+f[i]) == four )
                xFourBand(&calc2->i.c, usz+f[i], cmpl, &calc->i.c, copyVector,rank,calc->rt.runFlag );
            else if ( bodies(&calc2->i.c, usz+f[i]) == one )
                xOneBand(&calc2->i.c, usz+f[i], cmpl, &calc->i.c, copyVector,rank,calc->rt.runFlag );

//            for ( mdi = 0; mdi < mx ; mdi++)
//                calc->i.c.sinc.tulip[usr+i*mx+mdi].Current[cmpl] = 0;
            mdi = 0;
//            if  ( mx != 1 ){
//                for ( j = 0; j < CanonicalRank(&calc->i.c,copyVector,rank);j++){
//                    printf("%d %d %d %d %d ::%f::\n",i,f[i],j,((mdi)%mx),CanonicalRank(&calc->i.c,copyVector,rank),magnitude(&calc->i.c,ocean(rank, &calc->i.c, copyVector, j, rank)));
//
//                    tEqua(&calc->i.c, copyTwoVector,rank, ocean(rank, &calc->i.c, copyVector, j, rank),rank);
//                    //                fflush(stdout);
//                    tAddTw(&calc->i.c, usr+i*mx+((mdi)%mx), cmpl, copyTwoVector, rank);
//                    mdi++;
//                }
//            } else {
//printf("%d : %d _ %d\n", usr+ i , tPath(&calc->i.c, usr+i),tPath(&calc2->i.c,usz+f[i]));
                tEqua(&calc->i.c, usr+i*mx+((mdi)%mx), 0, copyVector, rank);
//            }
            //LOST VECTORS
            //            if( CanonicalRank(&calc->i.c, copyVector, rank) > part(&calc->i.c, usr+i)){
            //                la = CanonicalRank(&calc->i.c, copyVector, rank);
            //                calc->i.c.sinc.tulip[copyVector].Current[rank] = part(&calc->i.c,usr+i);
            //                tEqua(&calc->i.c, usr+i,cmpl, copyVector ,rank);
            //                calc->i.c.sinc.tulip[copyVector].Current[rank] = la;
            //                canonicalDecompositionMP(rank, &calc->i.c, copyVector, rank, usr+i, cmpl, calc->rt.vCONVERGENCE);
            //            } else {
            //                tEqua(&calc->i.c, usr+i,cmpl, copyVector ,rank);
            //            }
        }
        
    }
    for ( i = 0; i < ii*mx ; i++){
        tScale(&calc->i.c , usr+i,1./magnitude(&calc->i.c, usr+i));
    }
//    for ( k = 0; k < ii*mx ; k++){
//        for ( k2 = 0; k2 < ii*mx ; k2++)
//            printf("%f ..", tMultiplyMP(rank, &info, &calc->i.c, 1., -1, nullVector, 0, 'T', usr+k, 0, 'N', usr+k2, 0));
//        printf("\n");
//    }
//
//    for ( k = 0; k < 20 ; k++){
//        for ( k2 = 0; k2 < 20 ; k2++)
//            printf("%f ..", tMultiplyMP(rank, &info, &calc->i.c, 1., -1, nullMatrix, 0, 'T', ocean(rank,&calc->i.c,interactionExchange,k,0), 0, 'N', ocean(rank+1,&calc->i.c,interactionExchange,k2,0), 0));
//        printf("\n");
//    }

    return ii*mx;
}
