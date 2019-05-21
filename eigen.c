/*
 *  eigen.c
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


#include "eigen.h"

double vale ( struct sortClass * f ){
    INT_TYPE * mm = f->mmm+f->i*SPACE*2 ;
    double value=0.;
    INT_TYPE space;
    for ( space = 0; space < SPACE ; space++)
        if ( f->n1[space])
            value +=  f->str[space][mm[2*space] + mm[2*space+1]*f->n1[space]];
    return value;
}

double yale ( struct sortClass * f ){
    INT_TYPE * mm = f->mmm+f->i*6 ;
    return mm[3]+mm[4]*f->nG + mm[5]*f->nG*f->nG;
}


int sortComp (const void * elem1, const void * elem2)
{
    struct sortClass* f = ((struct sortClass*)elem1);
    struct sortClass* s = ((struct sortClass*)elem2);
    double valueF,valueS;
    valueF = vale(f);
    valueS = vale(s);
    if (valueF > valueS) return  1;
    if (valueF < valueS) return -1;
    return 0;
}

int sort2Comp (const void * elem1, const void * elem2)
{
    struct sortClass* f = ((struct sortClass*)elem1);
    struct sortClass* s = ((struct sortClass*)elem2);
    INT_TYPE valueF,valueS;
    valueF = yale(f);
    valueS = yale(s);
    if (valueF > valueS) return  1;
    if (valueF < valueS) return -1;
    return 0;
}

INT_TYPE tBoot1Construction(struct calculation * c1, enum division eigen){
    assignCores(&c1->i.c,1);
    struct field * f1 = &(c1->i.c);
    enum bodyType bootBodies = c1->rt.body;
    INT_TYPE n1[SPACE];
    length1(f1,n1);
    INT_TYPE N1,rank;
    INT_TYPE space,N2,i,ii,iii,iv,v;
    double * ar,*w;
    for ( space = 0 ;space < SPACE ; space++)
        if( f1->sinc.rose[space].body != nada){
            
            N1 = n1[space];
            N2 = N1*N1;
            zero(f1,copy,0);
            ar = myStreams(f1, canonicalBuffersBM, 0);
            w = ar + N2;
//            enum division pt = f1->sinc.tulip[Ha].linkNext;
//            do{
//                if ( f1->sinc.tulip[pt].space[space].body == one && f1->sinc.tulip[pt].space[space].block != id0 ){
//                    for ( i = 0; i < CanonicalRank(f1, pt, 0);i++)
//                        cblas_daxpy(N2, 1., streams(f1,pt,0,space), 1,ar,1 );
//                }
//
//                pt = f1->sinc.tulip[pt].linkNext;
//            }while ( pt != nullName);
//
            if ( c1->rt.calcType == electronicStuctureCalculation )
                cblas_dcopy(N2, streams(f1, kinetic,0,space)+space*N2, 1, ar, 1);
            else if ( c1->rt.calcType == clampProtonElectronCalculation ){
                if ( space < COMPONENT )
                    cblas_dcopy(N2, streams(f1, kinetic,0,space)+space*N2, 1, ar, 1);
                else {
                    zero(f1,kineticMass,0);
                    cblas_dcopy(N2,streams(f1, kineticMass,0,space), 1, ar, 1);
                    for ( i = 0; i < CanonicalRank(f1, protonRepulsion, 0);i++)
                        cblas_daxpy(N2, 1., streams(f1,protonRepulsion,0,space)+i*N2, 1,ar,1 );
                    // N12 >= 3
//                    exit(0);
                }
            }
//            if ( c1->i.springFlag ){
//                if ( space == 0 &&  ( c1->rt.runFlag )%2 == 0 ){
//                    cblas_daxpy(N2, 1.0,streams(f1, harmonium,0,0)+0*N2, 1, ar, 1);
//                }
//                if ( space == 1 &&  ( c1->rt.runFlag/2 )%2 == 0 ){
//                    cblas_daxpy(N2, 1.0,streams(f1, harmonium,0,1)+1*N2, 1, ar, 1);
//                }
//                if ( space == 2 &&  ( c1->rt.runFlag/4 )%2 == 0 ){
//                    cblas_daxpy(N2, 1.0,streams(f1, harmonium,0,2)+2*N2, 1, ar, 1);
//                }
//            }
            
            
            tdsyev(0, f1, 'V', N1, ar, N1, w);
            


            if ( bootBodies == one ){
                for ( i = 0  ; i < N1 ; i++)
                    if ( w[i] - w[0] < c1->i.level ){
                        streams(f1,foundationStructure,0,space)[i] = w[i]-w[0];
                        
                        streams(f1,foundationStructure,1,space)[i] = 1;
                    }
                cblas_dcopy(N2 , ar,1,myStreams(f1, bill1+space,0),1);
            }else if ( bootBodies == two ){
                
                
#ifdef OMP
#pragma omp parallel for private (v,rank,i,ii) schedule(dynamic,1)
#endif
                
                for ( v = 0 ; v < N2 ; v++){
#ifdef OMP
                    rank = omp_get_thread_num();
#else
                    rank = 0;
#endif
                    
                    i = (v)%N1;
                    ii = (v/N1)%N1;
                    if ( w[i]-w[0] + w[ii] - w[0] < c1->i.level ){

                    // if ( l2*l2 > i*i+ii*ii )
                    
                     //   printf("%d %d %f %f\n", i,ii, c1->i.level, w[i]-w[0] + w[ii] - w[0]);
                        
                        cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,space), 1);
                        f1->sinc.tulip[diagonal1VectorA].Current[rank] = 1;
                        cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,space), 1);
                        f1->sinc.tulip[diagonal1VectorB].Current[rank] = 1;
                        f1->sinc.tulip[diagonalVectorA].Current[rank] = 0;
                        tOuterProductSuOne(f1,space, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonalVectorA, rank);
                        cblas_dcopy(N2, streams(f1, diagonalVectorA,rank,space), 1, myStreams(f1,bill1+space,0)+v*N2, 1);
                        
                        
                        streams(f1,foundationStructure,0,space)[v] = w[i]-w[0] + w[ii] - w[0];
                        //
                        streams(f1,foundationStructure,1,space)[v] = 1;
                        //
                    }
                }
                
                
                
                
                
            }else if ( bootBodies == three ){
#ifdef OMP
#pragma omp parallel for private (v,rank,i,ii,iii) schedule(dynamic,1)
#endif
                
                for ( v = 0 ; v < N2*N1 ; v++){
#ifdef OMP
                    rank = omp_get_thread_num();
#else
                    rank = 0;
#endif
                    
                    i = (v)%N1;
                    ii = (v/N1)%N1;
                    iii = (v/(N2))%N1;
                    
                    if ( w[i]-w[0] + w[ii] - w[0]+w[iii]-w[0] < c1->i.level )
                    {
                    cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,space), 1);
                    double va = cblas_dnrm2(N1, streams(f1,diagonal1VectorA,rank,space), 1);

                    f1->sinc.tulip[diagonal1VectorA].Current[rank] = 1;
                    cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,space), 1);
                    va = cblas_dnrm2(N1, streams(f1,diagonal1VectorB,rank,space), 1);

                    f1->sinc.tulip[diagonal1VectorB].Current[rank] = 1;
                    f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
                    
                    struct name_label d2 = f1->sinc.tulip[diagonal2VectorA];
                    
                    
                    tOuterProductSuOne(f1, space,diagonal1VectorA, rank, diagonal1VectorB,rank, diagonal2VectorA, rank);
                    va = cblas_dnrm2(N1*N1, streams(f1,diagonal2VectorA,rank,space), 1);
                    cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,space), 1);
                    
                    f1->sinc.tulip[diagonal1VectorC].Current[rank] = 1;
                    f1->sinc.tulip[diagonal2VectorA].Current[rank] = 1;

                    f1->sinc.tulip[diagonalVectorA].Current[rank] = 0;
                    tOuterProductSuOne(f1,space, diagonal2VectorA, rank, diagonal1VectorC, rank, diagonalVectorA, rank);
                    
                    cblas_dcopy(N1*N2, streams(f1, diagonalVectorA,rank,space), 1, myStreams(f1,bill1+space,0)+v*N1*N2, 1);
                    
                    streams(f1,foundationStructure,0,space)[v] = w[i]-w[0] + w[ii] - w[0]+w[iii]-w[0];
                    streams(f1,foundationStructure,1,space)[v] = 1;
                    }
                }
            }else if ( bootBodies == four ){
                
                
#ifdef OMP
#pragma omp parallel for private (v,rank,i,ii,iii,iv) schedule(dynamic,1)
#endif
                
                for ( v = 0 ; v < N2*N2 ; v++){
#ifdef OMP
                    rank = omp_get_thread_num();
#else
                    rank = 0;
#endif
                    i = (v)%N1;
                    ii = (v/N1)%N1;
                    iii = (v/(N2))%N1;
                    iv = (v/(N2*N1))%N1;
                    
                    if( w[i]-w[0] + w[ii] - w[0]+w[iii]-w[0] + w[iv]-w[0]< c1->i.level ){
                    cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,space), 1);
                    
                    f1->sinc.tulip[diagonal1VectorA].Current[rank] = 1;
                    cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,space), 1);
                    
                    f1->sinc.tulip[diagonal1VectorB].Current[rank] = 1;
                    f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
                    tOuterProductSuOne(f1, space,diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
                    
                    cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,space), 1);
                    
                    f1->sinc.tulip[diagonal1VectorC].Current[rank] = 1;
                    
                    cblas_dcopy(N1, ar+iv*N1, 1, streams(f1,diagonal1VectorD,rank,space), 1);
                    
                    f1->sinc.tulip[diagonal1VectorD].Current[rank] = 1;
                    
                    f1->sinc.tulip[diagonal2VectorB].Current[rank] = 0;
                    f1->sinc.tulip[diagonalVectorA].Current[rank] = 0;
                    
                    tOuterProductSuOne(f1, space,diagonal1VectorC, rank, diagonal1VectorD, rank, diagonal2VectorB, rank);
                    tOuterProductSuOne(f1,space, diagonal2VectorA, rank, diagonal2VectorB, rank, diagonalVectorA, rank);
                    
                    cblas_dcopy(N2*N2, streams(f1, diagonalVectorA,rank,space), 1, myStreams(f1,bill1+space,0)+v*N2*N2, 1);
                    
                    streams(f1,foundationStructure,0,space)[v] = w[i]-w[0] + w[ii] - w[0]+w[iii]-w[0] + w[iv]-w[0];
                    streams(f1,foundationStructure,1,space)[v] = 1;
                    }
                }
                
            }
        }
    return 0;
}

INT_TYPE tMap (struct calculation * c1 ){
    struct field * f1 = & c1->i.c;
    size_t ms = MAXSTRING;
    char line0[MAXSTRING];
    char name[MAXSTRING];

    char *line = line0;
    INT_TYPE rank = 0;
    FILE  * list = NULL;

    
    sprintf(name, "%s.body", c1->name);
    list = fopen(name, "r");
    if ( list == NULL ){
        printf("nop");
        exit(0);
    }
    DCOMPLEX one = 1.;
    int lines= 0,n[6],r[6],space,N1;

    printVector(c1,c1->name,c1->name,-1,0, &one);

    getline(&line, &ms, list);
    while(!feof(list)){
        enum bodyType outBody = three;
        if ( outBody == three ){
            sscanf(line, "(%d,%d)(%d,%d)(%d,%d)", &n[0],&r[0], &n[1],&r[1],&n[2],&r[2]);
            
            for ( space = 0; space < SPACE ; space++)
                if ( f1->sinc.rose[space].body != nada )
                {
                    N1 = f1->sinc.rose[space].count1Basis;
                    cblas_dcopy(N1, streams(f1,f1->sinc.user + n[0],0,space)+N1*r[0], 1, streams(f1,diagonal1VectorA,rank,space), 1);
                    cblas_dcopy(N1, streams(f1,f1->sinc.user + n[1],0,space)+N1*r[1], 1, streams(f1,diagonal1VectorB,rank,space), 1);
                    cblas_dcopy(N1, streams(f1,f1->sinc.user + n[2],0,space)+N1*r[2], 1, streams(f1,diagonal1VectorC,rank,space), 1);

                }
            f1->sinc.tulip[diagonal1VectorA].Current[rank] = 1;
            f1->sinc.tulip[diagonal1VectorB].Current[rank] = 1;
            f1->sinc.tulip[diagonal1VectorC].Current[rank] = 1;

            
            f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
            f1->sinc.tulip[diagonal3VectorA].Current[rank] = 0;

            tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
            tOuterProductSu(f1, diagonal2VectorA, rank, diagonal1VectorC, rank, diagonal3VectorA, rank);
            if ( magnitude(f1, diagonal3VectorA) > 0.01 ){
                printf("%d %d %d | %d %d %d ", n[0],n[1],n[2],r[0],r[1],r[2]);
                printf("mag %f\n",magnitude(f1, diagonal3VectorA) );
                tFilename(c1->name, lines+1, outBody, 0, 0, name);
                FILE * outVector = fopen(name, "w");
                outputFormat(f1, outVector, diagonal3VectorA, rank);
                fclose(outVector);
                printVector(c1,c1->name,c1->name,lines,0, &one);
            }
        }
    
        lines++;
        getline(&line, &ms, list);

    }
    
    fclose(list);
    return 0;
}

INT_TYPE tSlam (struct field * f1,INT_TYPE allc, enum division vl, double fmax2){
    INT_TYPE tot =0,space,t,n1[SPACE];
    tot =  tFoundationLevel(f1, nullName, 0, fmax2, 1, nullName, 0, 0, 0, 0, NULL, 0, 0);
    
    for ( space = 0 ; space < SPACE ; space++)
        n1[space] = vectorLen(f1, space);
    if ( allc < tot ){
        printf("increase foundation %d or decrease levelLevel %f", tot, fmax2);
        exit(0);
    }
    INT_TYPE * mmm = malloc(sizeof(INT_TYPE ) * tot * SPACE *2),*mm;
    tot =  tFoundationLevel(f1, nullName, 0, fmax2, 0, nullName, 0, 0, 0, 0,mmm, 0, 0);

    for ( t = 0; t < tot ; t++){
        mm = mmm + 2*SPACE * t ;
        for ( space = 0; space < SPACE ; space++)
            if ( f1->sinc.rose[space].body != nada)
                cblas_dcopy(n1[space], myStreams(f1, bill1+space, 0)+(mm[2*space])*n1[space]+(mm[2*space+1])*n1[space]*n1[space],1,streams(f1,vl+t,0,space),1);
        f1->sinc.tulip[vl+t].Current[0] = 1;
    }
    free(mmm);
    printf("loaded %d alloc %d\n", tot, allc);
    return tot;
}


INT_TYPE tBootManyConstruction (struct calculation * c1){
    struct field * f1 = &(c1->i.c);
    enum bodyType bootBodies = c1->rt.body;
    INT_TYPE sp,cmpl,space,matrixNumber = c1->i.decomposeRankMatrix;
    //THRESING FLOOR
    tClear(f1, eigen);
    INT_TYPE r,i,im;
    INT_TYPE n2[SPACE];
    length(f1, eigen,n2);
    double mx= 1e99;
    INT_TYPE n1[SPACE];
    length(f1,eigenVectors,n1);
    DCOMPLEX * hmat = (DCOMPLEX*)myStreams(f1, matrixHbuild,0), sum ,minus = -1.;
    cmpl = 0;
        //for ( cmpl = 0; cmpl < spins(f1,eigen);cmpl++){
            tClear(f1, build);
            tClear(f1,eigen);
            zero(f1,eigen,cmpl);

            if (cmpl == 0 ){
                tClear(f1, copy);
                zero(f1,copy,0);
                tClear(f1,copy);
                tCycleDecompostionListOneMP(-1, f1, linear, 0,NULL, copy, 0, f1->mem1->rt->CANON, c1->i.decomposeRankMatrix, -1);
                tAddTw(f1, copy ,0, kinetic,0);
                tAddTw(f1,copy,0, vectorMomentum ,0);
                tSumMatrices(f1, build,0, copy);// B x (1+S) * 3
                if ( c1->rt.calcType == clampProtonElectronCalculation ){
                    zero(f1,copy,0);
                    tClear(f1,copy);
                    tId(f1,copy,0);
                    INT_TYPE N2 = vectorLen(f1, 3)*vectorLen(f1,3);
                        for ( i = 0; i < CanonicalRank(f1, protonRepulsion, 0);i++)
                            cblas_daxpy(N2, 1., streams(f1,protonRepulsion,0,3)+i*N2, 1,streams(f1,copy,0,3),1 );
                    tSumMatrices(f1, build,0, copy);// B x (1+S) * 3

                    tClear(f1, copy);
                    tCycleDecompostionListOneMP(-1, f1, interactionExchangePlus, 0,NULL,copy, 0, f1->mem1->rt->CANON, 1 , -1);
                    tSumMatrices(f1, build,0, copy);//B2:1 || B3 : 3
                }
                
                
                if ( bootBodies > one ){
                    if ( c1->rt.runFlag == 0 ){
                        if ( CanonicalRank(f1, interactionExchange, 0)){
                            tClear(f1,squareTwo);
                            tCycleDecompostionListOneMP(-1, f1, interactionExchange, 0,NULL,squareTwo, 0, f1->mem1->rt->CANON, 1 , -1);
                            tSumMatrices(f1, build,0, squareTwo);//B2:1 || B3 : 3
                        }
                    }else {
                        if ( CanonicalRank(f1, interactionEwald, 0)){
                            tClear(f1,squareTwo);
                            tCycleDecompostionListOneMP(-1, f1, interactionEwald, 0,NULL,squareTwo, 0, f1->mem1->rt->CANON, 1 , -1);
                            tSumMatrices(f1, build,0, squareTwo);//B2:1 || B3 : 3
                        }
                    }
                }
                
            }
//            else {
//                if ( CanonicalRank(f1, vectorMomentum,1 ) ){
//                    tEqua(f1, copy ,0, vectorMomentum,1);
//                    tSumMatrices(f1, build,0, copy);// B x (1+S) * 3
//                }
//            }
            //balance(f1,build,0);

            for ( space =  0; space < SPACE ; space++ )
                if ( f1->sinc.rose[space].body != nada){
                    tClear(f1,eigen);
                    zero(f1,eigen,0);
                    while ( CanonicalRank(f1, eigen, cmpl) < part(f1,eigen)){
                        tId(f1 ,eigen,cmpl);
                        canonicalListDecompositionMP(0, f1, NULL, build, 0, eigen, cmpl, f1->mem1->rt->CANON, magnitude(f1, build), space);
                     //   printf("%f-- \n", traceOne(f1,build,0));
                    }
                    for ( r = 0; r <part(f1,eigen) ; r++){
                        for ( i = 0; i < n2[space]; i++){
                            hmat[i] =0.;
                            for ( sp = 0 ; sp < 1;sp++)
                                if ( sp == 0 )
                                    hmat[i] += streams(f1, eigen,sp,space)[i+r*n2[space]] ;
                                else
                                    hmat[i] += I*streams(f1, eigen,sp,space)[i+r*n2[space]] ;
                            
                        }
//                        sum = 0.;
//                        for (i = 0; i < n1[space];i++){
//                            sum += hmat[i*n1[space]+i] ;
//                          //  printf("%d %d %f \n", r,i,creal(hmat[i*n1[space]+i]));
//                        }
                        //   if ( creal(sum ) < 0 )
                        //     cblas_zscal(n2[space], &minus, hmat, 1);
                        
                        tzheev (0,f1,'V',n1[space],hmat,n1[space],streams(f1,foundationStructure,0,space)+r*n1[space]);
                        
                        for ( i = 0; i < n2[space]; i++){
                            myStreams(f1, bill1+space,0)[i+r*n2[space]] = creal(hmat[i]);
                            if ( spins (f1, eigen) > 1 )
                                myStreams(f1, bill1+space,1)[i+r*n2[space]] = cimag(hmat[i]);
                        }
                        for (i = 0; i < n1[space];i++)
                            if (((streams(f1,foundationStructure,0,space)+r*n1[space])[i])-((streams(f1,foundationStructure,0,space)+r*n1[space])[0]) < c1->i.level  ){
                                (streams(f1,foundationStructure,1,space)+r*n1[space])[i] = 1;
                                if ( mx >(streams(f1,foundationStructure,0,space)+r*n1[space])[i]  )
                                    mx = (streams(f1,foundationStructure,0,space)+r*n1[space])[i]  ;
                            }
                        
                    }
                    for ( r = 0; r <part(f1,eigen) ; r++){
                        for (i = 0; i < n1[space];i++)
                            (streams(f1,foundationStructure,0,space)+r*n1[space])[i] -= mx;
                    }
                }
    return 0;
}



INT_TYPE tFoundationLevel( struct field * f1, enum division A , double lvlm, double lvlx,INT_TYPE ops,enum division build,INT_TYPE xB, double lvl1, double lvl2, double lvl3,INT_TYPE *mmm, INT_TYPE irrep,double seekPower){
    //GRID
    /// GRID
    ///// GRID
    ////////GRID
    INT_TYPE i,j,k,r1,r2,r3,space,ii,jj,kk,xx[SPACE];
    INT_TYPE vaMax=0, classicalBasisSize,*mm;
    enum bodyType bd = f1->mem1->rt->body;
    INT_TYPE nG = part(f1,eigen);//tSize(bd);
    INT_TYPE n1[SPACE];
    INT_TYPE r[SPACE];
    INT_TYPE GGG,g[SPACE],iii,flag;
    for ( space = 0 ; space < SPACE ; space++)
        n1[space] = vectorLen(f1, space);
    double value;
    classicalBasisSize = 0;
    //for( lvl = 0; lvl < imin(n1[0],imin(n1[1],n1[2])) ; lvl++)
        
    for ( space = 0; space < SPACE ;space++){
            if ( n1[space] ){
                xx[space] = 1;
                for ( i = 0; i < n1[space] ; i++)
                if ( streams(f1,foundationStructure,1,space)[i] > 0.5 ){
                    xx[space] = i+1;
                }
            }else
                xx[space] = 1;
       // printf("xx %d %d\n", xx[space],n1[space]);
    }
        INT_TYPE mx = 1;
        for ( space = 0 ; space < SPACE ; space++)
                mx *= nG*xx[space];
        
        for ( iii = 0; iii < mx ; iii++){
            
            GGG = 1;
            for ( space = 0; space < SPACE ; space++)
                {
                    g[space] = (iii/GGG) % nG;
                    GGG *= nG;
                    r[space] = (iii/GGG) % xx[space];
                    GGG *= xx[space];
                }

            
            
            flag = 1;
            value = 0;
            for ( space = 0; space < SPACE ; space++)
                if ( n1[space] != 0 ) {
                if (  streams(f1,foundationStructure,1,space)[r[space]+g[space]*n1[space]] > 0.5    ){
                    value += (streams(f1,foundationStructure,0,space)[r[space]+g[space]*n1[space]]);
                } else {
                    flag = 0;
                }
                }
            
            if ( flag )
            if ( vaMax < value )
                vaMax = value;
            
            
            if (flag &&  lvlm <= value && value < lvlx){
                {
                    //                                                                printf("**%f\n", value);
                    
                    
                    if ( 1 ){
                        //if ( SPACE <= 3 && ! tIR ( bd,g[0],g[1],g[2],irrep) ){
                        //    continue;
                       // }
                        if ( ops  == 0 ){
#if VERBOSE
                            printf("**%d-->%d : %d %d %d :: %d %d %d :: %f\n",irrep,classicalBasisSize,r[0],r[1],r[2],g[0],g[1],g[2],value);

#endif
                            mm = mmm+classicalBasisSize*SPACE*2;
                            
                            for ( space = 0 ; space < SPACE ; space++){
                                mm[2*space] = r[space];
                                mm[2*space+1] = g[space];
                            }
                        }
                       // printf("%d-->%d : %d %d %d :: %d %d %d :: %f\n",irrep,classicalBasisSize,r[0],r[1],r[2],g[0],g[1],g[2],value);

                        classicalBasisSize++;
                    }
                }
        }
        }
    if ( ops == 2 )
        return vaMax;
    else
        return classicalBasisSize;
}



INT_TYPE tFilter(struct field * f1, INT_TYPE Ve, INT_TYPE irrep, enum division usr){
    INT_TYPE ii,cmpl=0,rank,nP = tPerms(f1->mem1->rt->body);
    assignCores(f1, 2);
#ifdef OMP
#pragma omp parallel for private (ii,rank,cmpl) schedule(dynamic,1)
#endif
    for ( ii = 0; ii < Ve ; ii++)
    {
#ifdef OMP
        rank = omp_get_thread_num();
#else
        rank = 0;
#endif

        if ( irrep && bodies(f1, usr+ii ) > one  ){
            for ( cmpl = 0 ; cmpl < spins(f1, usr+ii) ; cmpl++){
                f1->sinc.tulip[permutationVector].Current[rank] = 0;
                tBuildIrr(rank, f1, irrep+nP, usr+ii, cmpl, permutationVector, rank);
        //        tCycleDecompostionGridOneMP(rank, f1, permutationVector, rank, NULL,usr+ii , cmpl, f1->mem1->rt->TARGET, part(f1,usr+ii), 1.);
                tCycleDecompostionListOneMP(rank, f1, permutationVector, rank, NULL,  usr+ii,cmpl,f1->mem1->rt->TARGET, part(f1,usr+ii), 1);
            }
        }
        
//#endif
        f1->sinc.tulip[usr+ii].value.symmetry = tClassify(rank, f1, usr+ii);
      //  printf(" %d->%d (%d) = %d\n",ii+usr, CanonicalRank(f1,usr+ii,cmpl),tPath(f1,usr+ii),f1->sinc.tulip[usr+ii].symmetry);
      //  fflush(stdout);
    }
    
    return 0;
}

INT_TYPE tSelect(struct field * f1, INT_TYPE Ve, INT_TYPE type, enum division usr, enum division usa, INT_TYPE testFlag){
    INT_TYPE sp,info,rank=0,maxEV = f1->sinc.maxEV,n,m;
    INT_TYPE stride = maxEV;
    double value;
    
	if ( ! CanonicalRank(f1, usa,0))
        return 0;


    if ( type ) {
        tClear(f1, usr+Ve);
        for ( sp = 0; sp < spins(f1, usr+Ve);sp++){
            f1->sinc.tulip[permutationVector].Current[rank] = 0;
            tBuildIrr(rank, f1, type, usa, sp, permutationVector, rank);
            tCycleDecompostionListOneMP(rank, f1, permutationVector, rank, NULL,usr+Ve, sp, f1->mem1->rt->vCANON, part(f1,usr+Ve), -1);
        }
        
        value = magnitude(f1, usr+Ve);
        
        
        if ( isnan(1./value) || isinf(1./value) )
           return 0;
        tScale(f1, usr+Ve, 1./value);

    } else {
	
        value = magnitude(f1, usa);
        
        
       if ( value < 1e-6 || isnan(value) || isinf(value) )
            return 0;
//        printf ( "tsel:%f\n", value);

        
        if ( testFlag >= 0 ){
            for ( sp = 0; sp < 1;sp++){

                tEqua(f1, usr+Ve,sp, usa,sp);
            }
            tScale(f1, usr+Ve, 1./value);
            
        }
    }
    if ( testFlag != 1 )
        return 1;

    
    DCOMPLEX HH, * S = (DCOMPLEX*)(myStreams(f1, matrixSbuild, 0));
    myZero(f1, matrixSbuild,0);
    double * ov = myStreams(f1, twoBodyRitz, 0);
    assignCores(f1, 2);
#ifdef OMP
#pragma omp parallel for private (m,n,rank,HH) schedule(dynamic,1)
#endif
    for ( n = 0; n < Ve+1 ; n++)
    {
        
#ifdef OMP
        rank = omp_get_thread_num();
#else
        rank = 0;
#endif
        
        for ( m = 0 ; m <= n   ; m++)    {
            matrixElements(rank, f1, usr+n, nullName, usr+m, &HH, S+n*stride+m);
            S[m*stride+n] = conj(S[n*stride+m]);
        }
        
    }
    
    assignCores(f1, 0);
    tzheev(0, f1, 'N', Ve+1,S, stride, ov);
    if ( testFlag ){

        if ( ov[Ve]/ov[0] < f1->mem1->rt->TOL && ov[Ve]/ov[0] >0 && ov[0] > 0 ){
                printf("\t%f\t%d\n",  ov[Ve]/ov[0],Ve+1 );
                fflush(stdout);

                return 1;
        } else {
        }
    }
    
    return 0;
}


INT_TYPE tCollect (struct field * f1, INT_TYPE irrep,enum division usz, INT_TYPE target,double seekPower){
    INT_TYPE  Ve = 0;
    f1->sinc.tulip[diagonalVectorA].header = Cube;
    f1->sinc.tulip[diagonalVectorB].header = Cube;
    struct name_label bd = f1->sinc.tulip[build];
    struct name_label eg = f1->sinc.tulip[eigen];

    INT_TYPE space, nG = part(f1,eigen);
    INT_TYPE flag = 1,ct = 0;
    double min0= 0,max0= 1+tFoundationLevel(f1, build,0,0,2,0,target,1e9,1e9,1e9,NULL,irrep,seekPower) ;
    double min = min0, max =  max0,vx= 100.,va = 1.;
    printf("\n\n\t| Seek %d Body vectors in %d components \t|\n", bodies(f1, eigenVectors), SPACE);
    printf("\t| Greater than Target \t: %6d\t\t|\n",target);
    
    ct = tFoundationLevel(f1, build,0,max0,1,0,target,1e9,1e9,1e9,NULL,irrep,seekPower);
    while(1){
        if ( flag ){
            if ( ct < target*va )
                min = 0.5*(max+min);
            else
                max = 0.5*(max+min);
        }
        ct = tFoundationLevel(f1, build,0,0.5*(min+max),1,usz,target,1e9,1e9,1e9,NULL,irrep,seekPower);
        printf("\t| Current \t: %6d \t %f\t\t|\n",ct,0.5*(min+max) );
        if ( max-min < 1e-9  ){
            printf("conv");
            exit(0);
        }
        if ( ct >= va*target && ct <= target*va*vx)
        {
            
            Ve = 0;
            INT_TYPE *mmm = malloc(sizeof(INT_TYPE ) *2*SPACE*ct),i,*mm;
            ct = tFoundationLevel(f1, build,0,0.5*(min+max),0,usz,target,1e9,1e9,1e9,mmm,irrep,seekPower);
            INT_TYPE n1[SPACE];
            {
                INT_TYPE space ;
                for ( space = 0 ; space < SPACE ; space++)
                    n1[space]= vectorLen(f1, space);
            }
            struct sortClass * sc = malloc ( sizeof( struct sortClass )*ct); ;
            for ( i = 0; i < ct ; i++){
                for ( space = 0; space < SPACE ; space++){
                    sc[i].n1[space] = n1[space];
                    sc[i].str[space] = streams(f1,foundationStructure,0,space);
                }
                sc[i].i = i;
                sc[i].nG = nG;
                sc[i].mmm = mmm;
            }
            
            qsort(sc, ct, sizeof(struct sortClass), &sortComp);
            
            
            for ( i = 0; i < ct ; i++){
                mm = mmm+sc[i].i*2*SPACE;
                tClear(f1,diagonalVectorA);
                f1->sinc.tulip[diagonalVectorA].Current[0] = 1;;
                
                for ( space = 0; space < SPACE ; space++)
                    if ( f1->sinc.rose[space].body != nada)
                    cblas_dcopy(n1[space], myStreams(f1, bill1+space, 0)+(mm[2*space])*n1[space]+(mm[2*space+1])*n1[space]*n1[space],1,streams(f1,diagonalVectorA,0,space),1);
                
                Ve =  tSASplit(f1, irrep, Ve,target, usz, diagonalVectorA);
                if ( Ve == target )
                    break;
            }
            
            if ( Ve == target )
                break;
            else
            {
                printf("instead try floorFlag %d\n", Ve);
                va *= 1.5;
                max = max0;
                printf("\t| increasing buffer region to %1.3f\t|\n", va);

                if ( va > 2 )
                    exit(2);
            }
            free(mmm);
            free ( sc);
        }
        
    }
    if ( target != Ve ){
        printf("ack no !\n");
        exit(0);
    }
    
    return Ve;
}


INT_TYPE tSASplit ( struct field * f1, INT_TYPE irrep , INT_TYPE Ve , INT_TYPE target,enum division usz, enum division vector){
    INT_TYPE map[24],nDeg=0,ii;
    
    if ( f1->mem1->rt->body == one || irrep == 0 ){
        nDeg = 1;
        map[1] = 0;
    }else
        
        if ( f1->mem1->rt->body == two ){
            if ( irrep == 1 ){
                map[1] = 1;
                nDeg = 1;
            }else
                if ( irrep == 2 ){
                    nDeg = 1;
                    map[1] = 2;
                }
        }else
            if ( f1->mem1->rt->body== three ){
                if ( irrep == 1 ){
                    map[1] = 1;
                    nDeg = 1;
                }else
                    if ( irrep == 2 ){
                        map[1] = 2;
                        nDeg = 1;
                    } else {
                        nDeg = 4;
                        map[1] = 3;
                        map[2] = 4;
                        map[3] = 5;
                        map[4] = 6;
                    }
            }
            else if ( f1->mem1->rt->body == four ){
                //
                if ( irrep == 1 ){
                    nDeg = 1;
                    map [1] = 1;
                }else if ( irrep == 2 ){
                    map[1] = 2;
                    nDeg = 1;
                } else if ( irrep == 3 ){
                    map[1] = 3;
                    map[2] = 4;
                    map[3] = 5;
                    map[4] = 6;
                    nDeg = 4;
                } else if ( irrep == 4){
                    map[1] = 7;
                    map[2] = 8;
                    map[3] = 9;
                    map[4] = 10;
                    map[5] = 11;
                    map[6] = 12;
                    map[7] = 13;
                    map[8] = 14;
                    map[9] = 15;
                    nDeg = 9;
                }else if ( irrep == 5 ){
                    map[1] = 16;
                    map[2] = 17;
                    map[3] = 18;
                    map[4] = 19;
                    map[5] = 20;
                    map[6] = 21;
                    map[7] = 22;
                    map[8] = 23;
                    map[9] = 24;
                    nDeg = 9;
                }
            }
    
    
    for ( ii = 1 ; ii <= nDeg ; ii++)
        
        if ( tSelect(f1, Ve, map[ii] , usz, vector, 1)){
            Ve++;
            if ( Ve == target )
                return Ve;
        
        }
    return Ve;
}


INT_TYPE tSquareVectors(struct calculation * c1, INT_TYPE EV2, enum division usz,enum division usr ){
    INT_TYPE i,j;
    struct field *f1 = & c1->i.c;

    for ( i = 0; i < EV2 ; i++)
        for ( j = 0 ; j < EV2 ; j++)
        {
            tClear(f1,usz+i*EV2+j );
            tOuterProductSu(f1, usz+i, 0, usz+j, 0, usr+i*EV2+j, 0);
            tOuterProductSu(f1, usz+i, 1, usz+j, 1, usr+i*EV2+j, 1);
        }    
    return EV2*EV2;
}


INT_TYPE tGreatDivideIteration (INT_TYPE translateFlag , double realPart, struct field * f1, enum division A , INT_TYPE I1, INT_TYPE I2, enum division usz, INT_TYPE foundation,INT_TYPE nMult, INT_TYPE shift){
    INT_TYPE expon,info;
    INT_TYPE rank ;
    //time_t start_t, lapse_t;
    DCOMPLEX temp,temp2, sum = 0,vhhv,vhv;
    //time(&start_t);
    INT_TYPE iii = 0;
    assignCores(f1, 1);
    
    for( expon = 1 ; foundation*expon < nMult  ; expon++){
        
//#ifdef OMP
//#pragma omp parallel for private (iii,rank,vhhv,vhv) schedule(dynamic,1)
//#endif
        for ( iii = 0; iii < foundation ; iii++)
        {
            
//#ifdef OMP
//            rank = omp_get_thread_num();
//#else
            rank = 0;
//#endif
            
            
            {
                tEquals(f1,usz+iii+expon*foundation,usz+(expon-1)*foundation+iii );//initialize
//                matrixElements(rank, f1, usz+iii+expon*foundation, nullName, usz+iii+expon*foundation, NULL, &vhhv);
//                tScale(f1, usz+iii+expon*foundation, 1./sqrt(cabs(vhhv)));

                tHXpX(rank, f1, A, translateFlag, translateFlag*realPart+ (!translateFlag)*1., 0.0, usz+iii+expon*foundation, f1->mem1->rt->TARGET , part(f1,usz+(expon)*foundation+iii),1 == 1);
                
                pMatrixElements( f1, usz+iii+expon*foundation, nullName, usz+iii+expon*foundation, NULL, &vhhv);
//                matrixElements(rank, f1, usz+iii+expon*foundation, nullName, usz+iii+(expon-1)*foundation, NULL, &vhv);

              //  printf("%d\t<\t%f\t | \t %f \t>\n", iii+1, creal(vhv), creal(vhhv) - sqr(creal(vhv)));
                if ( cabs(vhhv) > 0. )
                    tScale(f1, usz+iii+expon*foundation, 1./sqrt(cabs(vhhv)));
                else {
                    printf("oops norms \n");
                    exit(0);
                }
            }
            
        }
    }
    
    
    
    
    
    INT_TYPE sum2 = 0;
    expon = 1;
    for ( iii = 0; iii < foundation ; iii++)
        sum2 += CanonicalRank(f1, usz+(expon)*foundation+iii, 0)+CanonicalRank(f1, usz+(expon)*foundation+iii, 1);
    
    printf("-------\n\tAve-Canonical Rank \t| %2.1f \n", sum2*1. / foundation );    
    
    return nMult; 
}

INT_TYPE tLesserDivideIteration ( struct field * f1, enum division A , INT_TYPE I1, INT_TYPE I2, enum division usz, INT_TYPE foundation,INT_TYPE nMult, INT_TYPE shift);

INT_TYPE tEdges(struct calculation *c1, enum division vector){
    struct field * f1 = &c1->i.c;
    INT_TYPE info,spatial;
    DCOMPLEX ov,me;
    enum bodyType bootBodies = f1->body;
    if ( 1 ){
        //EDGES ALT
        enum block b,bx;
        INT_TYPE iii,jjj=1,dim,irrep;
        double sum = 0;
        //for ( irrep = 0 ;irrep <= 5 ; irrep++)
        {
              //  if ((! c1->i.irrep || f1->sinc.tulip[vector].value.symmetry  == irrep) && irrep == c1->i.irrep)
                {
                    bx = tv1;
                    if ( bootBodies == two )
                        bx = tv2;
                    else if ( bootBodies == three )
                        bx = tv3;
                    else if ( bootBodies == four )
                        bx = tv4;
                    for ( b = tv1 ; b <= bx; b++){
                        
                        for ( spatial = 0 ; spatial < 2 ; spatial++){
                            printf("electron %d:%d\t",b,spatial);

                            sum = 0;
                            for ( dim = 0; dim < COMPONENT ; dim++){
                                tClear(f1,edgeElectronMatrix );
                                if ( spatial == 0 )
                                    tEnd(f1, edgeElectronMatrix, 0, dim);
                                if ( spatial == 1 )
                                    tAlt(f1, edgeElectronMatrix, 0, dim);
                                me = 0.;
                                enum division u = edgeElectronMatrix+b;
                                struct name_label uu = f1->sinc.tulip[u];
                                pMatrixElements( f1, vector, edgeElectronMatrix+b, vector, &me, &ov);
                                printf("%1.8f ", (creal(me/ov)));
                                sum += (creal(me/ov));
                            }
                            printf(": %1.8f\n", sum);
                        }

                        if ( c1->rt.calcType == clampProtonElectronCalculation)
                            for ( spatial = 0 ; spatial < 2 ; spatial++){
                                printf("proton %d:%d\t",b,spatial);

                                sum = 0;
                                for ( dim = COMPONENT; dim < 2*COMPONENT ; dim++)
                                    if ( f1->sinc.rose[dim].body != nada ){
                                        
                                    tClear(f1,edgeProtonMatrix );
                                    if ( spatial == 0 )
                                        tEnd(f1, edgeProtonMatrix, 0, dim);
                                    if ( spatial == 1 )
                                        tAlt(f1, edgeProtonMatrix, 0, dim);
                                    me = 0.;
                                    pMatrixElements( f1, vector, edgeProtonMatrix+b, vector, &me, &ov);
                                    printf("%1.8f ", creal(me/ov));
                                    sum += (creal(me/ov));
                                }
                                printf(": %1.8f\n", sum);
                            }
                        
                    }
                    printf("\n\n");
                }
        }
    }
    return 0;
}

//INT_TYPE tConvergeTest (struct calculation * c1, enum division input){
//
//
//    {
//        struct field * f1 = &c1->i.c;
//        double sum2,sum;
//        enum division Mat, leftP,ml[1000];
//        INT_TYPE Lane = MaxCore,iii,i,i2,j,j2,rk,s1=0,cat=0,rank,iiii= 0;
//        double Sa [c1->i.bRank  *c1->i.bRank];
//        double Saa [c1->i.bRank  ];
//
//        leftP = input;
//
//
//        do {
//
//            if ( linear == name(f1,leftP) ){
//                Mat = rivers(Lane++, f1, linear, cat );
//                f1->sinc.tulip[Mat].block = f1->sinc.tulip[leftP].;
//            }else {
//                Mat = leftP;
//            }
//
//            if ( Lane -MaxCore >= MaxCore )
//            {
//                printf("neeed more \n");
//                exit(0);
//            }
//
//            if ( CanonicalRank(f1, Mat, 0)+CanonicalRank(f1, Mat, 1))
//                ml[iiii++] = Mat;
//
//            if (name(f1, leftP) == linear ){
//                cat++;
//
//                if ( cat > f1->Na ){
//                    leftP = f1->sinc.tulip[leftP].linkNext;
//                    cat = 0;
//                }
//
//            } else {
//                leftP = f1->sinc.tulip[leftP].linkNext;
//            }
//        } while ( leftP != nullName);
//
//
//
//        for ( iii = 0; iii < c1->i.nStates ; iii++){
//            sum = 0.;
//            sum2 = 0.;
//            rk = CanonicalRank(f1, eigenVectors+iii, 0);
//
//            for ( j = 0; j < iiii ; j++)
//            {
//
//#pragma omp parallel for private (i,i2,j2,rank) schedule(dynamic,1)
//                for ( i = 0; i < rk;i++)
//                {
//#if ARES
//
//                    rank = omp_get_thread_num();
//#else
//                    rank = 0;
//#endif
//
//                    //                                if ( Nb){
//                    //                                    tEqua(f1,squareVector, rank, ocean(rank,f1,eigenVectors+iii,i,0),0);
//                    //                                    tpProduct(rank, f1, ml[j], 0, squareVector, rank);
//                    //                                }
//                    //                                else
//                    tMultiplyMP(rank,f1,1., -1,squareVector, rank, 'N', ml[j] , 0,'N', ocean(rank,f1,eigenVectors+iii,i,0), s1);
//
//
//                    Saa[i] = tMultiplyMP(rank,f1,1., -1,nullVector, 0, CDT, squareVector, rank,'N', ocean(rank,f1,eigenVectors+iii,i,0),0);
//
//                    for ( j2 = 0 ; j2 < iiii ; j2++)
//                        for ( i2 = 0 ; i2 <= i   ; i2++){
//                            //                                        if ( Nb){
//                            //                                            tEqua(f1, totalVector, rank, ocean(rank,f1,eigenVectors+iii,i2,0),0);
//                            //                                            tpProduct(rank, f1, ml[j2] , 0, totalVector, rank);
//                            //                                        }
//                            //                                        else
//                            tMultiplyMP(rank,f1,1., -1,totalVector, rank, 'N',ml[j2] , 0,'N', ocean(rank,f1,eigenVectors+iii,i2,0), s1);
//
//                            Sa[i*rk+i2] = tMultiplyMP(rank,f1,1., -1,nullVector, 0, CDT, squareVector,rank,'N', totalVector, rank);
//
//                        }
//                }
//
//                for ( i = 0; i < rk ; i++){
//                    for ( i2 = 0; i2 < rk ; i2++){
//                        sum2 += Sa[i*rk+i2];
//                    }
//                    sum += Saa[i];
//                }
//            }
//
//            printf("Cv%lld \t %f\n", iii+1, sum2-sqr(sum));
//
//        }
//    }
//}



INT_TYPE tEigenCycle (struct field * f1, enum division A ,char permutation,  INT_TYPE Ne, enum division usz, INT_TYPE quantumBasisSize ,INT_TYPE iterations, INT_TYPE foundation, INT_TYPE irrep,INT_TYPE flag,  enum division outputSpace, enum division outputValues){
    INT_TYPE in,gvOut,prevBuild;
    time_t start_t, lapse_t;
    myZero(f1, matrixHbuild, 0);
    myZero(f1, matrixSbuild, 0);
    time(&start_t);
    enum division Mat;
    INT_TYPE cmpl,cmpl2,cmpl3,cat,iii = 0,maxEV = f1->sinc.maxEV,rank;
    INT_TYPE stride = maxEV;
    double * ritz = myStreams(f1, outputValues, 0);
    double * overlap = myStreams(f1, conditionOverlapNumbers, 0);
    enum division el ;
    DCOMPLEX *T  =  (DCOMPLEX *) myStreams(f1, matrixHbuild,0/*CORE RANK*/);
    DCOMPLEX *S  =  (DCOMPLEX *) myStreams(f1, matrixSbuild,0/*CORE RANK*/);
    INT_TYPE powerMat;

    INT_TYPE info,n,m,s1,g,r,rr,rx,gx,a,a2;
    enum division leftP ;
    prevBuild = 0;
    for ( n = prevBuild; n < quantumBasisSize ; n++)
    {
        if ( magnitude(f1, usz+n) == 0. )
        {
            printf("normaly, I would punch you!\n");
            exit(0);
        }
 //       tScale(f1, usz+n, 1./magnitude(f1, usz+n));
    }
    
    assignCores(f1, 2);

        leftP = A;
        
        do {
            Mat = leftP;
        
#if 1
            if ( Rank(f1,name(f1, Mat))){
                INT_TYPE space;
                for ( space = 0 ;space < SPACE ; space++)
                    printf("%2d", f1->sinc.tulip[name(f1,Mat)].space[space].body);
                if (CanonicalRank(f1, name(f1,Mat), 1))
                    printf("\t::\t %d \t| (%d+i%d)\t",name(f1,Mat), CanonicalRank(f1, name(f1,Mat), 0),CanonicalRank(f1, name(f1,Mat), 1));
                else
                    printf("\t::\t %d \t| (%d)\t",name(f1,Mat), CanonicalRank(f1, name(f1,Mat), 0));
                for ( space = 0 ;space < SPACE ; space++)
                    printf("%2d", f1->sinc.tulip[Mat].space[space].block);
                printf("\t%f\n", traceOne(f1,name(f1, Mat),0));
                printf("\n");
            }
            fflush(stdout);
#endif
#ifdef OMP
#pragma omp parallel for private (in,m,rank,n) schedule(dynamic,1)
#endif
            for ( in = 0 ;in < quantumBasisSize * quantumBasisSize; in++)
            {
                
#ifdef OMP
                rank = omp_get_thread_num();
#else
                rank = 0;
#endif
                m = in % quantumBasisSize;
                n = (in/quantumBasisSize) % quantumBasisSize;
                if ( m <= n ){
                    matrixElements(rank, f1, usz+n, Mat, usz+m, T+n*stride+m, S+n*stride+m);
                    T[m*stride+n]= conj(T[n*stride+m] ) ;
                    S[m*stride+n] = conj(S[n*stride+m ]) ;
                }
            }
            
            leftP = f1->sinc.tulip[leftP].linkNext;
        } while ( leftP != nullName);

  
    if(0){
        INT_TYPE ii;
        for ( ii = 0 ; ii < quantumBasisSize ; ii++){
            printf("%d %f %f d \n", ii+1, creal(T[ii*stride+ii]), creal(S[ii*stride+ii]));
            printf("%d %f %f x \n", ii+1, creal(T[ii*stride]), creal(S[ii*stride]));
        }
    }
    
    if ( 0){
        printf("\n\nH={");
        for ( m = 0; m < quantumBasisSize ; m++){
            if ( m )
                printf("},{");
            else
                printf("{");
            for ( n =0 ; n < quantumBasisSize ; n++){
                if ( n )
                    printf(",");
                printf("%f+I*%f", creal(T[n*stride+m]),cimag(T[n*stride+m]));
            
            }
        }
        printf("}};");

        
    }
    if ( 0){
        printf("S={");
        for ( m = 0; m < quantumBasisSize ; m++){
            if ( m )
                printf("},{");
            else
                printf("{");
            for ( n =0 ; n < quantumBasisSize ; n++){
                if ( n )
                    printf(",");
                printf("%f+I*%f", creal(S[n*stride+m]),cimag(S[n*stride+m]));
                
            }
        }
        printf("}};\n\n");
        
        
    }


    //    if(0)
//    for ( n = 0; n < quantumBasisSize-foundation; n++)
//        for ( m = 0 ; m <= n   ; m++){
//            if ( part(f1,usz+m) > minRank && part(f1, usz+n) > minRank && n % foundation == m % foundation){
//               // printf("%lld %lld %f %f\n", n/foundation,m/foundation,T[n*stride+m],T[n*stride+m+maxEV*maxEV*2]);
//                T[n*stride+m]= T[n*stride+m+maxEV*maxEV*2];
//                T[m*stride+n] = T[n*stride+m];
//            }
//        }
    //assignCores(f1, 1);
    
//    f1->mem1->rt->matrixElementsTotal += countTot;
//    f1->mem1->rt->matrixElementsZero += countLam;
//    time(&lapse_t);
//    f1->mem1->rt->buildTime += difftime(lapse_t, start_t);
//    time(&start_t);]
      if (1){
        assignCores(f1, 0);

        char Job = 'N';
//        if (flag)
            Job = 'V';

       // cblas_zcopy(maxEV*stride , T , 1 , g+maxEV*stride , 1);
          if ( 0 ){
              cblas_zcopy(maxEV*stride , S , 1 , S+maxEV*stride , 1);
              tzheev(0, f1, 'N', quantumBasisSize, S+maxEV*stride, stride, overlap);
              if (  (overlap[0] > 0.) && overlap[quantumBasisSize-1]/overlap[0] < f1->mem1->rt->TOL && flag <= 2)
                  printf(" Krylov-2\t %f \n",  overlap[quantumBasisSize-1]/overlap[0]);
              else         if (  (overlap[0] > 0.) && overlap[quantumBasisSize-1]/overlap[0] < f1->mem1->rt->TOL && flag == 3)
                  printf(" Krylov-3 \t %f \n",  overlap[quantumBasisSize-1]/overlap[0]);
              else         if (  (overlap[0] > 0.) && overlap[quantumBasisSize-1]/overlap[0] < f1->mem1->rt->TOL && flag == 4)
                  printf(" Krylov-4 \t %f \n",  overlap[quantumBasisSize-1]/overlap[0]);
              else {
                  printf("Linear dependent! %f\t%1.16f\n",overlap[quantumBasisSize-1],overlap[0]);
                  return quantumBasisSize;
              }
              
              cblas_zcopy(maxEV*stride , S , 1 , S+maxEV*stride , 1);
          }
        gvOut = tzhegv (0,f1,Job,quantumBasisSize,T,S,stride,ritz);
#if VERBOSE
          printf("eigenSolved\n");
#endif
        fflush(stdout);
        assignCores(f1, 1);
    }
    time(&lapse_t);

    //f1->mem1->rt->eigenTime += difftime(lapse_t, start_t);
    time(&start_t);
    
    {       //printf("\nHeader,Number,CEG,Linear,BODY,CLASS,WEIGHT\n");
        for ( iii = 0; iii < imin(quantumBasisSize,Ne) ; iii++)
        {
            f1->sinc.tulip[eigenVectors+iii].value.value = ritz[iii];
            printf("Press%d:,%1.15f, %f\n", iii+1, ritz[iii],cblas_dznrm2(quantumBasisSize,T+iii*stride, 1));
        }
    //    fflush(stdout);
    }



    
//    if (flag == 3 ){
//        if ( Ne > quantumBasisSize ){
//            printf ("warning basis too small\n");
//            exit(0);
//        }
//        el = outputSpace;
//
//
//        INT_TYPE u;
//        printf(" foundation %d < quan %d\n", foundation, quantumBasisSize );
//        for ( u = 0; u < quantumBasisSize ;u++){
//            if ( u < foundation )
//                printf("*");
//            printf ("USZ %d ---%d)\n", u+usz, tPath(f1, usz+u));
//
//
//        }
//
////        for ( u = 0; u < Ne ; u++)
////            f1->sinc.tulip[el+u].path = -1;
//
//        INT_TYPE sp,nG= tSize(f1->body),path,fpath,ipath,i,h,xx,six,ii2,jj2,kk2,fx,nm,v;
//        double mag2,*pointers[MaxCore];
//
//        f1->sinc.tulip[eigenList].ptRank[0] = 0;
//        f1->sinc.tulip[eigenList].purpose = ptObject;
//        {
//            ipath = 0;
//            fpath = 0;
//            for ( v = 0 ; v < nG*nG*nG ; v++ )
//            {
//                ii2 = v % nG;
//                jj2 = (v/nG)%nG;
//                kk2 = (v/(nG*nG))%nG;
//                        if (tIR(f1->body,ii2, jj2 , kk2, irrep))
//                        {
//                            path = v;
//                            fx = 0;
//                            six = 0;//should be clustered..
//                            for ( xx = 0; xx < quantumBasisSize ; xx++)
//                                if ( tPath(f1, usz+xx) == path ){
//                                    if  ( xx < foundation )
//                                        if (! (fx++) ){
//                                            f1->sinc.tulip[eigenList].name = usz+xx;
//                                    //        printf("begin ");
//                                        }
//                                    six+= spins(f1, usz+xx)*part(f1, usz+xx);
//                                   // printf ("%d,%d %d\n",  usz+xx, spins(f1, usz+xx),part(f1, usz+xx));
//                                }
//                           // printf(" fx %d\n", fx);
//                            f1->sinc.tulip[eigenList].Current[0] = six;
//
//
//                            if ( CanonicalRank(f1, eigenList, 0) )
//                                for ( cmpl = 0; cmpl < spins(f1, usz) ; cmpl++)
//                                {
//                                    printf("%d:+:%d \n", path, CanonicalRank(f1, eigenList, 0));
//
//#ifdef OMP
//#pragma omp parallel for private (iii,rank,rr,h,i,nm,cmpl2,r,sp,mag2) schedule(dynamic,1)
//#endif
//                                    for ( iii = 0; iii < Ne  ; iii++)
//                                    {
//#ifdef OMP
//                                        rank = omp_get_thread_num();
//#else
//                                        rank = 0;
//#endif
//                                        pointers[rank] = myStreams(f1, canonicalBuffersC, rank);
//                                        rr = 0;
//                                        //loop over actual memory order...
//                                        for ( h = 0; h < fx ;h++)
//                                            for ( i = 0; i < quantumBasisSize ; i+=foundation)
//                                            {
//                                                nm = h+fpath + i ;
//                                                if ( nm >= quantumBasisSize ){
//                                                    printf("**%d %d %d %d **", nm,h,fpath, i );
//                                                    exit(0);
//                                                }
//
//                                              //      printf ( "%d:: %d <-",iii+1, nm );
//                                               // fflush(stdout);
//                                                for ( cmpl2 = 0; cmpl2 < spins(f1,usz + nm ) ; cmpl2++)
//                                                    for ( r = 0; r < part(f1,usz + nm ); r++){
//                                                        if ( r < CanonicalRank(f1,usz + nm,cmpl) &&  cmpl2 == cmpl && tPath (f1,usz + nm) == path){
//                                                            pointers[rank][rr++] = (vectors[cmpl]+(iii)*stride)[nm];
//                                                        }else{
//                                                            pointers[rank][rr++] = 0.0;
//                                                        }
//                                                    //    printf ( "->%d\n", rr );
//                                                    }
//                                                //fflush(stdout);
//
//                                            }
//                                     //   printf("%d %d ==%f\n", path, iii+1,cblas_dnrm2 ( rr, pointers[rank],1) );
//                                        if ( rr != CanonicalRank(f1, eigenList, 0) ){
//                                            printf("%d %d exit\n", rr, CanonicalRank(f1, eigenList, 0) );
//                                            exit(0);
//                                        }
//                                        if ( cblas_dnrm2 ( rr, pointers[rank],1) > 1e-6 ){
//                                            if ( rr > part(f1, canonicalBuffersC))
//                                            {
//                                                printf("crap!\n %lld %lld",rr , part(f1, canonicalBuffersC));
//                                                exit(0);
//                                            }
//                                            nm = ipath * Ne + iii;
//                                      //      fflush(stdout);
//
//                                            tCycleDecompostionListOneMP(rank,f1, eigenList, pointers[rank],el+nm, cmpl, f1->mem1->rt->vCANON, part(f1, el+nm), 1.);
////                                            f1->sinc.tulip[el+nm].path = path;
//                                            mag2 = 0;
//                                            for ( sp = 0; sp < spins(f1, el+nm) ; sp++)
//                                                mag2 += inner(rank, f1, el+nm, sp);
//                                          //  printf("mag %f\n", sqrt(mag2));
//
//                                            tScale(f1, el+nm, 1./sqrt(mag2));
//                               //             printf("assign %d <= %d\n", nm , path);
//                                        }else {
//
//                                        }
//                                    }
//                                }
//
//                            ipath += 1;
//                            fpath += fx;
//                        }
//        }
//        }
//        if ( fpath != foundation ){
//            printf ("wanring !!! fpath %d foundation %d", fpath, foundation);
//            fflush(stdout);
//            exit(0);
//        }
//
//        time(&lapse_t);
//        printf("\nSpectrum TIME \t %15.0f\n", difftime(lapse_t, start_t));
//        fflush(stdout);
//
//    }else
        if ( flag == 4 )
        {
            double norms,*pointers[MaxCore];
            if (1){
                f1->sinc.tulip[eigenList].name = usz;
                f1->sinc.tulip[eigenList].species = vector;
                f1->sinc.tulip[eigenList].Partition= 0;
                for( g = 0; g < quantumBasisSize ; g++)
                    for ( cmpl = 0 ;cmpl < spins(f1,usz+g) ; cmpl++)
                        for ( r = 0; r < part(f1, usz+g); r++){
                            f1->sinc.tulip[eigenList].Partition++;//memory inline. need this.
                        }
                // printf("[[%d %d %lld %lld]]\n", eigenList, usz,f1->sinc.tulip[usz].Address, streams(f1, eigenList, 0, 0 )-streams(f1, usz, 0, 0 ));
                for ( cmpl = 0 ;cmpl < spins(f1,eigenVectors) ; cmpl++)
                    
                    for ( powerMat = 1 ; powerMat <= 1 ; powerMat++){
                        
                        el = outputSpace + (powerMat-1)*Ne;
                        //                    printf("begin\n");
                        //                    fflush(stdout);
//#ifdef OMP
//#pragma omp parallel for private (iii,rank,rr,g,r,cmpl2) schedule(dynamic,1)
//#endif
                        for ( iii = 0; iii < imin(quantumBasisSize,Ne) ; iii++)
                        {
                            
//#ifdef OMP
//                            rank = omp_get_thread_num();
//#else
                            rank = 0;
//#endif
                            pointers[rank] = myStreams(f1, canonicalBuffersC, rank);
                            
                            //                        printf( "%d %d\n",iii+1, cblas_idamax(stride, (vectors[cmpl]+(iii)*stride), 1));
                            
                            rr = 0;
                            norms = 0.;
                            for( g = 0; g < quantumBasisSize ; g++){
                                for ( cmpl2 = 0; cmpl2 < spins(f1, usz+g) ; cmpl2++)
                                    for ( r = 0; r < part(f1, usz+g); r++){
                                        if (  r >= CanonicalRank(f1, usz+g, cmpl2)){
                                            pointers[rank][rr++] = 0.0;
                                        }else{
                                            if( cmpl2 == cmpl ){
                                                if ( cmpl == 0 )
                                                    pointers[rank][rr++] = creal((T + iii*stride)[g]);
                                                else if ( cmpl == 1 )
                                                    pointers[rank][rr++] = cimag((T + iii*stride)[g]);
                                            }
                                            else
                                                pointers[rank][rr++] = 0.0;//mask off other *SPIN*( or complex vector).
                                        }
                                          //  printf("%d %d %d %f\n", usz+g,cmpl2,r,myStreams(f1, canonicalBuffersC, rank)[rr-1] );
                                    }
                            }
                            
//                            for ( g = 0 ; g < 35 ; g++)
//                                printf("%1.2f,",pointers[rank][g]);
//                            printf("\n\n");
                            if ( rr > part(f1, canonicalBuffersC))
                            {
                                printf("crap!\n %lld %lld",rr , part(f1, canonicalBuffersC));
                                fflush(stdout);

                                exit(0);
                            }
                            tCycleDecompostionGridOneMP(-2,f1, eigenList,0, pointers[rank],el+iii   , cmpl, f1->mem1->rt->vCANON, part(f1, el+iii), f1->mem1->rt->powDecompose);
                        
//                            INT_TYPE jj;
  //                          for ( jj = 0 ;jj <= iii ; jj++)
    //                            printf("_%f_", tMultiplyMP(rank, &info, f1, 1., -1, nullName, 0, 'T', el+jj, 0,'N',el+iii, 0));
                        
                        
                        }
                    }
            }
            time(&lapse_t);
            printf("\nFinal Spectrum TIME \t %15.0f\n", difftime(lapse_t, start_t));
            fflush(stdout);
        }
    
    
    //f1->sinc.tulip[matrixHbuild].Current[0] = quantumBasisSize;
    //f1->sinc.tulip[matrixSbuild].Current[0] = quantumBasisSize;

    return 0;
}



