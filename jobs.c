/*
 *  jobs.c
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

#include "jobs.h"
#undef __cplusplus
#include "Faddeeva.cc"


INT_TYPE countHam ( struct calculation *c1 , struct field f1 ){
    INT_TYPE spin,sum = 0;
    enum division di;
    
    for ( spin = 0 ; spin < f1.f.cmpl ; spin++){
        tClear(f1.f,hamiltonian);

        for ( di = Ha ; di!= nullName; di = f1.f.tulip[di].linkNext){
            if ( CanonicalRank(f1.f,f1.f.tulip[di].name, spin)){
                printf("%d %d %d %f %d\n", di, f1.f.tulip[di].name, spin, traceOne(f1.f, di, spin),CanonicalRank(f1.f, di, spin));
                sum++;
            }
        
        }
    }
    return sum;
}

INT_TYPE foundation1(struct calculation *c1, struct field f1){
    INT_TYPE EV;
    f1.i.Iterations = 1;
    if ( 1 ){
        iModel(c1,&f1);
       // separateKinetic(f1.f, 0,kinetic, 1,electron);
        tBoot1Construction(c1,f1.f ,build);
        if ( ! tSortBoot(c1,f1.f,build) )
            EV =   tSlam(f1.f,f1.i.qFloor,f1.f.user,c1->i.level);
        else
            EV = 0;
        
        
#ifdef OVERFLAG
            tInnerTest(f1.f, kinetic, copy);

            INT_TYPE i;
            for ( i = 0; i < EV ; i++){
                printf("%f\n", magnitude(f1.f, f1.f.user+i));
            }   //
////testSAAgain(f1.f, f1.f.user+i);
//            }
#endif
//            tEigenCycle(1,f1.f,overlap1,CDT, f1.i.nStates, f1.f.user,EV,0, EV,0,1,eigenVectors,twoBodyRitz);

    }else {
        exit(1);
    }
    
    INT_TYPE ii,flag = 1;;
#ifndef OVERFLAG
        print(c1,f1,1,0,EV , f1.f.user);
#endif
    fModel(&f1.f);

    return EV;
}

INT_TYPE foundationM(struct calculation *c1, struct field f1){
    INT_TYPE EV = 0,space;
    f1.i.Iterations = 1;
    

    if ( allowQ(f1.f.rt, blockAllocateHamiltonianBlock) ){
        iModel(c1,&f1);
        
        INT_TYPE spin ;
        enum division onem;
        if ( bodies(f1.f,hamiltonian ) == one)
            onem = diagonalCube;
        else
            onem = quadCube;
        tClear(f1.f, onem);
        
        {
            for ( spin = 0 ; spin < f1.f.cmpl ; spin++){
                tClear(f1.f, hamiltonian);
                tId(f1.f , onem,spin);
                enum division ha1 = Ha;
                while ( ha1 != nullName ){
                    printf("%f--", traceOne(f1.f, ha1, 0));
                    if ( bodies(f1.f, name(f1.f,ha1) ) == one ){
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)
                                sumTo2(f1.f, 1.,space,f1.f.tulip[ha1].space[space].block, name(f1.f,ha1), spin, hamiltonian, 0);
                        f1.f.tulip[hamiltonian].Current[0] += CanonicalRank(f1.f, ha1, 0);
                    }
                    else if ( bodies (f1.f,name(f1.f,ha1 ) ) == two )
                        tAddTw(f1.f, hamiltonian, 0, name(f1.f,ha1), spin);
                    ha1 = f1.f.tulip[ha1].linkNext;
                }
                
                tCycleDecompostionGridOneMP(-1, f1.f, hamiltonian, 0, NULL,onem  , spin, c1->rt.CANON, 1, 0);
                switch(bodies(f1.f,eigen)){
                    case one:
                        tEqua(f1.f, eigen, spin, onem, spin);
                        break;
                    case two:
                        tEqua(f1.f, eigen, spin, onem, spin);
                        break;
                    case three:
                        tClear(f1.f, build);
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)
                                sumTo3(f1.f,1.0, space, e12,onem, spin,build, 0);
                            
                        f1.f.tulip[build].Current[0] += 1;
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)
                                sumTo3(f1.f, 1.0,space, e13,onem, spin,build, 0);
                        f1.f.tulip[build].Current[0] += 1;
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)
                                sumTo3(f1.f,1.0, space, e23,onem, spin,build, 0);
                        f1.f.tulip[build].Current[0] += 1;
        
                        tId(f1.f , eigen,spin);
                        tCycleDecompostionGridOneMP(-2, f1.f, build, 0, NULL,eigen , spin, c1->rt.CANON, 1, 0);
                        break;
                    case four:
                        tClear(f1.f, build);
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)
                                sumTo4(f1.f,1.0, space, e12,onem, spin,build, 0);
                                f1.f.tulip[build].Current[0] += 1;
                                
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)

                                sumTo4(f1.f,1.0, space, e13,onem, spin,build, 0);
                                f1.f.tulip[build].Current[0] += 1;
                                
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)

                                sumTo4(f1.f,1.0, space, e23,onem, spin,build, 0);
                                f1.f.tulip[build].Current[0] += 1;
                                
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)

                                sumTo4(f1.f,1.0, space, e14,onem, spin,build, 0);
                                f1.f.tulip[build].Current[0] += 1;
                                
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)

                                sumTo4(f1.f, 1.0,space, e24,onem, spin,build, 0);
                                f1.f.tulip[build].Current[0] += 1;
                        for (space = 0; space < SPACE ;space++)
                            if ( f1.f.rose[space].body != nada)

                                sumTo4(f1.f,1.0, space, e34,onem, spin,build, 0);
                                f1.f.tulip[build].Current[0] += 1;
                                
                            
                        tId(f1.f , eigen,spin);
                        tCycleDecompostionGridOneMP(-2, f1.f, build, 0, NULL,eigen , spin, c1->rt.CANON, 1, 0);
                        break;
                }
            }
            tBootManyConstruction(c1,f1.f ,eigen);
            EV =   tSlam(f1.f,f1.i.qFloor,f1.f.user,c1->i.level);
            print(c1,f1,1,0,EV , f1.f.user);
            
        }
        fModel(&f1.f);
        
    }

    return EV;
}



INT_TYPE krylov ( struct calculation *c1, struct field f1){
    INT_TYPE EV = 0,i,fi,next;

#ifndef APPLE
    f1.i.qFloor = countLinesFromFile(c1,f1,0,&f1.i.iRank, &f1.i.xRank);
    //count canonical-rank...
#endif
    
    f1.i.nStates =f1.i.Iterations  ;
   
 iModel(c1,&f1);
    countHam(c1,f1);

    for ( fi =0 ; fi < f1.i.files ; fi++){
        tLoadEigenWeights (c1,f1, f1.i.fileList[fi], &EV,f1.f.user , f1.i.collect);
    }
        if (EV == 0 ){
#ifdef APPLE
            EV =1 ;
            tId(f1.f, f1.f.user, 0);
            
            
#else
        print(c1,f1,1,0,0,eigenVectors);
        return 1;
#endif
        }
    INT_TYPE RdsSize = EV,iterator=0;

    if(1){
        printf ("Step \t%d\n", iterator);
        INT_TYPE cmpl,g ;
        //        for ( iii = 0; iii < EV ; iii++){
        //            printf ( "\n Vector \t%d \n", iii+1);
        //            for ( sp = 0 ; sp < spins(f1.f, f1.f.user); sp++)
        //                tCycleDecompostionGridOneMP(-1, f1.f, f1.f.user+iii, sp, NULL, eigenVectors+iii, sp, c1->rt.vCANON, part(f1.f,eigenVectors+iii), -1);
        
        for ( cmpl = 0; cmpl < spins(f1.f, f1.f.user) ; cmpl++){
            tClear(f1.f,totalVector);
            if ( f1.f.cat && f1.i.filter  && f1.i.irrep )
            {
                printf("cat filter\n");
                fflush(stdout);
                for( g = 0; g < EV ; g++)
                    tBuild3Irr(0, f1.f, f1.i.irrep, f1.f.user+g, cmpl, totalVector, 0);
            }
//            else if ( f1.i.filter  && f1.i.irrep ){
//                for( g = 0; g < EV ; g++)
//                    tBuildIrr(0, f1.f, f1.i.irrep, f1.f.user+g, cmpl, totalVector, 0);
//            }
            else
            {
                for( g = 0; g < EV ; g++)
                    tAddTw(f1.f, totalVector, 0, f1.f.user+g, cmpl);
            }
                    
            tCycleDecompostionGridOneMP(-2, f1.f, totalVector, 0, NULL,eigenVectors , cmpl, c1->rt.vCANON, part(f1.f,eigenVectors), c1->rt.powDecompose);
        }
        double norm = magnitude(f1.f, eigenVectors );
        if ( norm > c1->rt.TARGET ){
            printf("Normed from %f\n", norm );
            fflush(stdout);
            
            tScaleOne(f1.f, eigenVectors, 0, 1/norm);
            //testSA(f1.f,eigenVectors);
        }
        else
        {
            print(c1,f1,1,0,0,eigenVectors);
            return 1;
        }
        EV = 1;
        RdsSize = 1;
        fflush(stdout);
        tFilter(f1.f, EV, 0, eigenVectors);//classify
        printExpectationValues(f1.f, Ha, eigenVectors);
        fflush(stdout);
        print(c1,f1,1,0,1,eigenVectors);
        if ( c1->rt.runFlag == 0 )

        tEdges(f1.f, eigenVectors);

    }
    
    INT_TYPE flag;
    for ( iterator = 1 ; iterator < f1.i.Iterations ; iterator++){
        
        
       
            flag = 1;
        
        
        
        if ( ! tGreatDivideIteration(c1->i.shiftFlag, c1->i.shiftVector[iterator-1][0],c1->i.shiftVector[iterator-1][1],  f1.f,Iterator, 1,0,eigenVectors+RdsSize-EV,EV,2*EV,0)){
            RdsSize += EV;
            
            if(1){
                tFilter(f1.f, EV, !(!f1.i.filter )* f1.i.irrep, eigenVectors+RdsSize-EV);//filter
                printf ("Step \t%d\n", iterator);
                fflush(stdout);
                INT_TYPE iii ;
                next = 0;
                for ( iii = 0; iii < 1 ; iii++){
                    printf ( "\n Vector \t%d \t %d\n", iii+1, +RdsSize-EV+iii);
                //    tEigenCycle(1, f1.f, Ha, 1, RdsSize, eigenVectors, RdsSize, RdsSize, 1, 0, 0, nullName, twoBodyRitz);
                    printExpectationValues(f1.f, Ha, eigenVectors+RdsSize-EV+iii);
                    //assume EV = 1;
                    print(c1,f1,0,RdsSize-EV+iii,RdsSize-EV+iii+1,eigenVectors);
                   
                    if ( c1->rt.runFlag == 0 )
                        next = tEdges(f1.f , eigenVectors+RdsSize-EV+iii);
                    else
                        next = 0;
                    fflush(stdout);

                }
                if ( next ){
                    printf("advice:\t %d\n",next);
                    break;
                }
            }
        }else {
            break;
        }
        fflush(stdout);
    }
    fModel(&f1.f);

    return 0;
}

INT_TYPE decompose ( struct calculation *c1, struct field f1){
    //count canonical-rank...
    
    iModel(c1,&f1);
    
    if (0){
        inputFormat(f1.f, c1->name, hamiltonian, 1);
        tClear(f1.f,trainHamiltonian);
    } else {
        tEqua(f1.f, hamiltonian, 0, linear, 0);
    }
    tCycleDecompostionGridOneMP(-2, f1.f, hamiltonian, 0, NULL,trainHamiltonian , 0, c1->rt.CANON, part(f1.f,trainHamiltonian), c1->rt.powDecompose);
    printf("%f\n", traceOne(f1.f, linear, 0));
    ioStoreMatrix(f1.f, trainHamiltonian, 0, "trainLinear.matrix", 0);
    fModel(&f1.f);
    return 0;
}

INT_TYPE spitGauss ( struct calculation *c1, struct field f1){
    //load vectors...
    INT_TYPE fi,EV=0;
       f1.i.qFloor = countLinesFromFile(c1,f1,0,&f1.i.iRank, &f1.i.xRank);
       //count canonical-rank...

    switch(f1.i.body ){
        case one :
            f1.i.nStates =1+c1->i.gaussCount;
        case two:
            f1.i.nStates =1+c1->i.gaussCount*c1->i.gaussCount;
    }
       f1.i.Iterations = 1  ;
      
       iModel(c1,&f1);

       for ( fi =0 ; fi < f1.i.files ; fi++){
           tLoadEigenWeights (c1,f1, f1.i.fileList[fi], &EV,eigenVectors , 0);
       }
           if (EV == 0 ){
           return 1;
           }

    //categorize by body!

    
    INT_TYPE nn,G,g,g2,space,r,n,nnn,i,i2,i3,rank,l ;
    switch ( f1.i.body ){
        case one:
            for ( g = 0 ; g < c1->i.gaussCount ; g++)
                for ( space = 0; space < SPACE ; space++)
                    if ( f1.f.rose[space].body != nada){
                        
                        n = vectorLen(f1.f, space);
                        for ( r = 0; r < CanonicalRank(f1.f, eigenVectors, 0); r++)
                            streams(f1.f,eigenVectors,0,space)[g+n*r] = 0;


            }
//            for ( g = 0 ; g <= c1->i.gaussCount ; g++){
//                          if ( g == 0 )
//                          {
//                              print(c1,f1,1,0,1,eigenVectors+g);
//                          }else {
//                              tClear(f1.f   , eigenVectors+g);
//                              tId(f1.f, eigenVectors+g, 0);
//                              zero(f1.f,eigenVectors+g,0);
//
//                              for ( space = 0; space < SPACE ; space++)
//                                  if ( f1.f.rose[space].body != nada)
//                                      streams(f1.f,eigenVectors+g,0,space)[g-1] = 1;
//                              print(c1,f1,0,g,1+g,eigenVectors);
//                          }
//
//            }
            break;
            case two:
            for ( space = 0; space < SPACE ; space++)
                if ( f1.f.rose[space].body != nada){
                    n = vector1Len(f1.f, space);
                    nn = vectorLen(f1.f, space);

                    for ( i = 0 ; i < n ; i++)
                        for ( i2 = 0 ; i2 < n ; i2++)
                        {
                            for ( r = 0; r < CanonicalRank(f1.f, eigenVectors, 0); r++)
                                if ( i < c1->i.gaussCount || i2 < c1->i.gaussCount )
                                    streams(f1.f,eigenVectors,0,space)[i+n*i2+nn*r] = 0;
                        }
                
                }
//            G = c1->i.gaussCount;
//            for ( g = 0 ; g <= G ; g++)
//                for ( g2 = 0 ; g2 <= G ; g2++){
//                              if ( g == 0 && g2 == 0 )
//                              {
//                                  print(c1,f1,1,g+g2*G,1+g+g2*G,eigenVectors+g+g2*G);
//                              }else if ( g != 0 && g2 != 0 ) {
//                                  tClear(f1.f,eigenVectors+g+g2*G);
//                                  tId(f1.f,eigenVectors+g+g2*G,0);
//                                  zero(f1.f,eigenVectors+g+g2*G,0);
//
//                                  for (space = 0; space < SPACE ; space++)
//                                      if (f1.f.rose[space].body != nada){
//                                          n = vector1Len(f1.f, space);
//                                          nn = vectorLen(f1.f, space);
//                                          for ( r = 0; r < CanonicalRank(f1.f, eigenVectors, 0); r++)
//                                          streams(f1.f,eigenVectors+g+g2*G,0,space)[(g-1)+(g2-1)*n+nn*r] = 1;
//                                      }
//                                  print(c1,f1,0,g+g2*G,1+g+g2*G,eigenVectors);
//                              }else if ( g == 0 ) {
//                                  tClear(f1.f,eigenVectors+g+g2*G);
//                                  tId(f1.f,eigenVectors+g+g2*G,0);
//                                  zero(f1.f,eigenVectors+g+g2*G,0);
//                              }else if ( g2 == 0 ) {
//                                  tClear(f1.f,eigenVectors+g+g2*G);
//                                  tId(f1.f,eigenVectors+g+g2*G,0);
//                                  zero(f1.f,eigenVectors+g+g2*G,0);
//                              }
//
//
//                }
                break;
            case three:
            for ( space = 0; space < SPACE ; space++)
                if ( f1.f.rose[space].body != nada){
                    n = vector1Len(f1.f, space);
                    nnn = vectorLen(f1.f, space);

                    for ( i = 0 ; i < n ; i++)
                        for ( i2 = 0 ; i2 < n ; i2++)
                            for ( i3 = 0 ; i3 < n ; i3++)

                        {
                            for ( r = 0; r < CanonicalRank(f1.f, eigenVectors, 0); r++)
                                if ( i < c1->i.gaussCount || i2 < c1->i.gaussCount|| i3 < c1->i.gaussCount)
                                    streams(f1.f,eigenVectors,0,space)[i+n*i2+n*n*i3+nnn*r] = 0;
                        }
                
                }
            break;
    }
    print(c1,f1,0,0,1,eigenVectors);

    fModel(&f1.f);

    return 0;
}


INT_TYPE ritz( struct calculation * c1, struct field f1){
    size_t ms = MAXSTRING;
    char line0[MAXSTRING];
    char filename[MAXSTRING];    char str[SUPERMAXSTRING];


    char * line = line0;
    //c1->i.canonRank = 0;
    INT_TYPE fi,EV = 0,i,j,sp;
    double va;
    INT_TYPE lines = 0,flines = 0;

    f1.i.qFloor = countLinesFromFile(c1,f1,0,&f1.i.iRank, &f1.i.xRank);
  //  printf("-->%d\n", f1.i.iRank);
   // c1->i.iRank = c1->i.bRank;
    iModel(c1,&f1);
    
    countHam(c1,f1);
    
    for ( fi =0 ; fi < f1.i.files ; fi++)
        tLoadEigenWeights (c1,f1, f1.i.fileList[fi],&EV, f1.f.user,f1.i.collect);//UNUSUAL!!!
    if (EV == 0 ){
        printf ("no ritz vectors!\n");
        exit(0);
    }
    
    {
        INT_TYPE typer = 1;
        
    tEigenCycle(typer,c1->i.flipSignFlag , f1.f,Ha,CDT, f1.i.nStates, f1.f.user,EV,0, EV,0,1,eigenVectors,twoBodyRitz);
    }
    
    DCOMPLEX *V = (DCOMPLEX*)myStreams(f1.f,matrixHbuild,0);
    INT_TYPE iii,ii,stride = f1.f.maxEV;
    for ( iii = 0; iii < f1.i.nStates ; iii++){
        {
            FILE * outf ;
            sprintf(str, "%s-%d.vector",c1->name,iii+1);
            outf = fopen (str,"w");
            fclose(outf);
            sprintf(str, "%s-%d",c1->name,iii+1);
        }
        for ( ii = 0; ii < EV ; ii++)
            printVector(c1,f1.f,f1.f.tulip[f1.f.user+ii].value.title,str,f1.f.tulip[f1.f.user+ii].value.stage-1,f1.i.irrep, V+stride*iii+ii);
    }
    fModel(&f1.f);

    return 0;
}

//INT_TYPE svd ( struct calculation c, struct field f1){//E
//    if ( ! f1.i.filesVectorOperator){
//        printf("vectorOperator flag\n");
//        exit(1);
//    }
//
//    INT_TYPE EV = foundation(&c,f1);
//    printf("finished foundation\n");
//    fflush(stdout);
//    tEigenCycle(0,f1.f,Ha,'T', f1.i.nStates, f1.f.user,EV,0, EV,0,0,eigenVectors,twoBodyRitz);
//
//    INT_TYPE i,j=0;
//    for ( i = 0 ; i < f1.i.nStates ; i++ ){
//        if ( -myStreams(f1.f, twoBodyRitz, 0)[i] > c.rt.TARGET ){
//            printf("load %f\n", -myStreams(f1.f, twoBodyRitz, 0)[i]);
//            j++;
//        }
//    }
////    tEigenCycle(0,f1.f,Ha,'T', j, f1.f.user,EV,0, EV,0,4,eigenVectors,twoBodyRitz);
////
////    {
////           cblas_dscal(j, -1., myStreams(f1.f, twoBodyRitz, 0), 1);
////    }
//    return 0;
//}

INT_TYPE sumTo2(struct sinc_label f1,double scalar, INT_TYPE space,enum block bl, enum division mat,INT_TYPE ms, enum division sum,INT_TYPE spin){
    enum block bol;
    INT_TYPE I1,I2, I3, I4,r,N1;
    double value;
    for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
        for ( bol = tv1 ; bol <= tv2 ; bol++){
            {
                N1 = vector1Len(f1, space);
                Stream_Type * stream = streams(f1,sum,spin,space)+N1*N1*N1*N1*(CanonicalRank(f1, sum, spin)+r);
                for ( I1 = 0 ; I1 < N1 ; I1++)//body 0
                    for ( I2 = 0 ; I2 < N1 ; I2++)//body 0
                        for ( I3 = 0 ; I3 < N1 ; I3++)//body 1
                            for ( I4 = 0 ; I4 < N1 ; I4++)//body 1
                            if ( bol == bl)
                            {
                                if ( bol == tv1 ){
                                    value  = streams(f1,mat,ms,space)[ I1*N1+I2 + r*N1*N1 ] * delta(I3-I4);
                                }else if (bol == tv2){
                                    value  = streams(f1,mat,ms,space)[ I3*N1+I4 + r*N1*N1 ] * delta(I1-I2);
                                }
                                if ( space == 0 )
                                    value *= scalar;
                                stream[ (I1+I3*N1)+ ( I2+I4*N1)*N1*N1 ] = value;
                            }
            }
        }
    return 0;
}

INT_TYPE sumTo3(struct sinc_label f1,double scalar,INT_TYPE space,  enum block bl, enum division mat,INT_TYPE ms, enum division sum,INT_TYPE spin){

    INT_TYPE n2[SPACE];
    length(f1, sum,n2);

    if ( bodies ( f1, sum ) == three){

        if ( bodies ( f1, mat ) == one ){
            enum block bol;
            INT_TYPE I1,I2, I3, I4,I5,I6,r;
            INT_TYPE n1[SPACE];
            length1(f1,n1);
            double value;

            for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
                for ( bol = tv1 ; bol <= tv3 ; bol++){
                    {
                        Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*(CanonicalRank(f1, sum, spin)+r);
                        if ( bol == bl )
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                            {
                                                value = 0;
                                                if ( bol == tv1 ){
                                                    value  = streams(f1,mat,ms,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6);
                                                }else if ( bol == tv2 ) {
                                                    value  = streams(f1,mat,ms,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6);
                                                }else if ( bol == tv3 ) {
                                                    value  = streams(f1,mat,ms,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4);
                                                }
                                                if ( space == 0 )
                                                value *= scalar;
                                                stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]] = value;
                                            }
                    }
                }
        }else if ( bodies ( f1, mat ) == two ){
            enum block pair ;
            INT_TYPE I1,I2, I3, I4,I5,I6,r,ve;
            INT_TYPE n1[SPACE];

            length1(f1, n1);

            double value;

            for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
                for ( pair = e12 ; pair <=e23  ; pair++){
                    {
                        Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*(CanonicalRank(f1, sum, spin)+r);
                        if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
                            printf("part sum\n");
                            exit(0);
                        }
                        if ( pair == bl )
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)

                                            {
                                                value = 0;
                                                if ( pair == e12 ){
                                                    value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6);
                                                }else if ( pair == e13 ) {
                                                    value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4);
                                                }else if ( pair == e23 ) {
                                                    value  = streams(f1,mat,ms,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2);
                                                }
                                                ve = (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space];
                                                //                                                printf("%f %lld %lld\n", value,ve,n2[space]);
                                                //                                                fflush(stdout);
                                                if ( space == 0 )
                                                value *= scalar;
                                                stream[ ve ] = value;
                                                //                                                printf("x");
                                                //                                                fflush(stdout);

                                            }
                    }
                }
        }
        else {
            printf("Yo!");
            exit(0);
        }

    }
    return 0;
}








INT_TYPE sumTo4(struct sinc_label f1,double scalar,INT_TYPE space, enum block bl, enum division mat,INT_TYPE ms, enum division sum,INT_TYPE spin){

    INT_TYPE n2[SPACE];
    length(f1, sum,n2);
    if (bodies(f1,sum) == four && f1.rt->calcType == electronicStuctureCalculation){

    if ( bodies ( f1, mat ) == one ){
        enum block bol ;
        INT_TYPE I1,I2, I3, I4,I5,I6,I7,I8,r;
        INT_TYPE n1[SPACE];
        length1(f1,n1);
        double value;

        for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
            for ( bol = tv1 ; bol <= tv4 ; bol++){
                {
                    Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*(CanonicalRank(f1, sum, spin)+r);
                    if ( bol == bl )
                    for ( I1 = 0 ; I1 < n1[space] ; I1++)
                        for ( I2 = 0 ; I2 < n1[space] ; I2++)
                            for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                    for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                        for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                            for ( I7 = 0 ; I7 < n1[space] ; I7++)
                                                for ( I8 = 0 ; I8 < n1[space] ; I8++)
                                                        
                                                {
                                                    value = 0;
                                                    if ( bol == tv1 ){
                                                        value  = streams(f1,mat,ms,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6)*delta(I7-I8);
                                                    }else if ( bol == tv2 ) {
                                                        value  = streams(f1,mat,ms,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6)*delta(I7-I8);
                                                    }else if ( bol == tv3 ) {
                                                        value  = streams(f1,mat,ms,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I7-I8);
                                                    }else if ( bol == tv4 ) {
                                                        value  = streams(f1,mat,ms,space)[ I7*n1[space]+I8 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I5-I6);
                                                    }
                                                    if ( space == 0 )
                                                    value *= scalar;
                                                    stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
                                                }
                }
            }
    }else if ( bodies ( f1, mat ) == two ){
        enum block pair;
        INT_TYPE I1,I2, I3, I4,I5,I6,I7,I8,r;
        INT_TYPE n1[SPACE];
        length1(f1, n1);

        double value;

        for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
            for ( pair = e12 ; pair <= e34 ; pair++){
                {
                    Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*(CanonicalRank(f1, sum, spin)+r);
                    if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
                        printf("part sum\n");
                        exit(0);
                    }
                    if ( pair == bl )
                    for ( I1 = 0 ; I1 < n1[space] ; I1++)
                        for ( I2 = 0 ; I2 < n1[space] ; I2++)
                            for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                    for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                        for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                            for ( I7 = 0 ; I7 < n1[space] ; I7++)
                                                for ( I8 = 0 ; I8 < n1[space] ; I8++)

                                                {
                                                    //                                            0    e12,     1,3     2,4
                                                    //                                            1    e13,     1,5     2,6
                                                    //                                            2    e23,     3,5     4,6
                                                    //                                            3    e14,     1,7     2,8
                                                    //                                            4    e24,     3,7     4,8
                                                    //                                            5    e34      5,7     6,8


                                                    if ( pair == e12 ){
                                                        value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6)*delta(I7-I8);
                                                    }else if ( pair == e13 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I7-I8);
                                                    }else if ( pair == e23 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I7-I8);
                                                    }else if ( pair == e14 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I7) + (I2*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I5-I6);
                                                    }else if ( pair == e24 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I3*n1[space]+I7) + (I4*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I5-I6);
                                                    }else if ( pair == e34 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I5*n1[space]+I7) + (I6*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I3-I4);
                                                    }else {
                                                        printf ("rails!\n");
                                                        exit(0);
                                                    }
                                                    if ( space == 0 )
                                                    value *= scalar;
                                                    stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
                                                }
                }
            }
    }
    else {
        printf("Yo!");
        exit(0);
    }
    }
    return 0;
}

INT_TYPE report ( struct calculation c, struct field f1){

    if ( CanonicalRank(f1.f,interactionExchange,0) ){
        ioStoreMatrixScale(&f1,interactionExchange ,0,"interactionExchange.matrix",0);
    }
    if ( CanonicalRank(f1.f,interactionExchange,1) ){
        ioStoreMatrixScale(&f1,interactionExchange ,1,"interactionExchange.1.matrix",0);
    }

    if ( CanonicalRank(f1.f,interactionEwald,0) ){
        ioStoreMatrixScale(&f1,interactionEwald ,0,"interactionEwald.matrix",0);
    }
    if ( CanonicalRank(f1.f,interactionEwald,1) ){
        ioStoreMatrixScale(&f1,interactionEwald ,1,"interactionEwald.1.matrix",0);
    }

//    if ( CanonicalRank(f1.f,shortenPlus,0) )
//        ioStoreMatrix(f1.f,shortenPlus ,0,"shortenExchangePlus.matrix",0);
//
//    if ( CanonicalRank(f1.f,shortenMinus,0) )
//        ioStoreMatrix(f1.f,shortenMinus ,0,"shortenExchangeMinus.matrix",0);


    if ( CanonicalRank(f1.f,linear,0) ){
        ioStoreMatrix(f1.f,linear ,0,"linear.matrix",0);

    }

    if ( CanonicalRank(f1.f,linear,1) ){
        ioStoreMatrix(f1.f,linear ,1,"linear.1.matrix",0);

    }
    
    if ( CanonicalRank(f1.f, vectorMomentum , 0 ) ) {
        ioStoreMatrix(f1.f,vectorMomentum, 0,"vector.matrix",0);
    }
    if ( CanonicalRank(f1.f, vectorMomentum , 1 ) ) {
        ioStoreMatrix(f1.f,vectorMomentum, 1,"vector.1.matrix",0);
    }

    return 0;
}



INT_TYPE distill ( struct calculation c, struct field f1){
    
    iModel(&c, &f1);
    tClear(f1.f, hamiltonian);
    
        
        if ( c.rt.runFlag > 0 )
        {
            enum bodyType bootBodies = f1.f.rose[0].body;
            INT_TYPE N1;
            if (GAS == 1){

                    enum division in = interactionEwald;
                    enum division out = jelliumElectron;
                    INT_TYPE r, space,m,n,nl,Nl,n2,m2,cmpl;
                    tClear(f1.f, out);
                    tClear(f1.f, copy);
                    tId(f1.f, copy,0);
                    for ( cmpl = 0 ; cmpl < f1.f.cmpl; cmpl++){
                        for ( r = 0 ; r < CanonicalRank(f1.f, in, cmpl); r++){

                            for ( space = 0; space < SPACE ; space++)
                            {
                                nl = vector1Len(f1.f, space);
                                Nl = nl/2;
                                for ( n = 0; n < nl ; n++ )
                                    for ( m = 0 ; m < nl ; m++)
                                    {
                                        (streams(f1.f, copy, 0, space))[nl*m+n] = 0.;
                                        for ( n2 = 0; n2 < Nl ; n2++ )
                                            for ( m2 = 0 ; m2 < Nl ; m2++)
                                                (streams(f1.f, copy, 0, space))[nl*m+n] +=  (streams(f1.f, in, cmpl, space)+r*nl*nl*nl*nl)[nl*nl*(nl*m2+m)+(nl*n2+n)]/Nl;
                                    }
                            }
                            tAddTw(f1.f, out, cmpl, copy, 0);
                        }
                        tScaleOne(f1.f, out, cmpl, -(INT_TYPE)(bootBodies));
                        printf("jellium-%d %f\n", cmpl,traceOne(f1.f, out, cmpl));
                    }
                }
            

            {
                double offset=0.,sum=0.,sumt=0.,prod;
                enum division in = interactionEwald;
                INT_TYPE r, space,m,n,nl,Nl,n2,m2,cmpl;
                for ( cmpl = 0 ; cmpl < 1; cmpl++){
                    for ( r = 0 ; r < CanonicalRank(f1.f, in, cmpl); r++){
                        prod = 1.;
                        for ( space = 0; space < SPACE ; space++)
                        {
                            sum = 0.;
                            nl = vector1Len(f1.f, space);
                            Nl = nl/2;
                            for ( n = 0; n < Nl ; n++ )
                                for ( m = 0 ; m < Nl ; m++)
                                    for ( n2 = 0; n2 < Nl ; n2++ )
                                        for ( m2 = 0 ; m2 < Nl ; m2++)
                                        {
                                            sum +=  (streams(f1.f, in, cmpl, space)+r*nl*nl*nl*nl)[nl*nl*(nl*n+n2)+(nl*m+m2)]/(Nl*Nl);
                                        }
                            prod *= sum;
                        }
                        sumt += prod;
                    }
                }
                printf("jellium background\t%15.15f\n", sumt);
                switch ( bootBodies ) {
                    case one:
                        offset = (0.5)*sumt;
                        break;
                        //0
                        //0.5*2
                        // -1
                    case two:
                        offset = (2*0.5+1)*sumt/2.;
                        // 0
                        //=
                        //0.5 * 4 (FORM PLANES)--> 2 ewald + 2 constants
                        //1 (pair of +)
                        //1 (pair of - )
                        //-4 (together)
                        break;
                    case three:
                        offset = (3*0.5+3.)*sumt/3.;
                        break;
                        // 0
                        //=
                        //0.5 * 6 (FORM PLANES)--> 3 ewald + 3 constants
                        //3 (trio of +)
                        //3 (trio of - )
                        //-9 (together)

                        
                        
                        //1 ewald + 1 constant
                        //1 of +
                        //1 of -
                        //-3 JELLIUM
                        
                    case four:
                        offset = (4*0.5+6.)*sumt/4.;
                        break;
                    case five:
                        offset = (5*0.5+10.)*sumt/5.;
                        break;
                    case six:
                        offset = (6*0.5+16.)*sumt/6.;
                        break;

                }
                tClear(f1.f, copy);
                tId(f1.f, copy, 0);
                tScaleOne(f1.f, copy, 0, offset);
                tAddTw(f1.f, jelliumElectron, 0, copy, 0);
                printf("jelliumT-%d %f\n", 0,traceOne(f1.f, jelliumElectron, 0));

            }
            
            if (GAS==1 ){
                enum division in = interactionEwald;
                enum division out = intercellularSelfEwald;
                INT_TYPE r, space,m,n,nl,Nl,cmpl;
                tClear(f1.f, out);
                tClear(f1.f, copy);
                tId(f1.f, copy,0);
                for ( cmpl = 0 ; cmpl < f1.f.cmpl; cmpl++){
                    for ( r = 0 ; r < CanonicalRank(f1.f, in, cmpl); r++){
                        for ( space = 0; space < SPACE ; space++)
                        {
                            nl = vector1Len(f1.f, space);
                            Nl = nl/2;
                            for ( n = 0; n < nl ; n++ )
                                for ( m = 0 ; m < nl ; m++)
                                {
                                    (streams(f1.f, copy, 0, space))[nl*m+n] =  (streams(f1.f, in, cmpl, space)+r*nl*nl*nl*nl)[nl*nl*(nl*m+m)+(nl*n+n)];
                                }
                        }
                        tAddTw(f1.f, out, cmpl, copy, 0);
                    }
                    printf("interaction-%d %f\n", cmpl,traceOne(f1.f, out, cmpl));
                }
                tScaleOne(f1.f, out, 0, 0.5);
                tScaleOne(f1.f, out, 1, 0.5);

            }

            if (GAS == 1){
                enum division in = interactionExchange;
                enum division out = intracellularSelfEwald;
                INT_TYPE r, space,m,n,Nl,nl,cmpl;
                tClear(f1.f, out);
                tClear(f1.f, copy);
                tId(f1.f, copy,0);
                for ( cmpl = 0 ; cmpl < f1.f.cmpl; cmpl++){
                    for ( r = 0 ; r < CanonicalRank(f1.f, in, cmpl); r++){

                        for ( space = 0; space < SPACE ; space++)
                        {
                            nl = vector1Len(f1.f, space);
                            Nl = nl/2;
                            for ( n = 0; n < nl ; n++ )
                                for ( m = 0 ; m < nl ; m++)
                                {
                                    (streams(f1.f, copy, 0, space))[nl*m+n] =  (streams(f1.f, in, cmpl, space)+r*nl*nl*nl*nl)[nl*nl*(nl*m+m)+(nl*n+n)];
                                }
                        }
                        tAddTw(f1.f, out, cmpl, copy, 0);
                    }
                    printf("intraction-%d %f\n", cmpl,traceOne(f1.f, out, cmpl));
                }
                tScaleOne(f1.f, out, 0, -0.500);
                tScaleOne(f1.f, out, 1, -0.500);
            }
#ifndef APPLE
            if ( CanonicalRank(f1.f,jelliumElectron,0) ){
                ioStoreMatrix(f1.f,jelliumElectron ,0,"jelliumElectron.matrix",0);
            }

            if ( CanonicalRank(f1.f,jelliumElectron,1) ){
                ioStoreMatrix(f1.f,jelliumElectron ,1,"jelliumElectron.1.matrix",0);
            }

            if ( CanonicalRank(f1.f,intracellularSelfEwald,0) ){
                ioStoreMatrix(f1.f,intracellularSelfEwald ,0,"intracellularSelfEwald.matrix",0);
            }

            if ( CanonicalRank(f1.f,intracellularSelfEwald,1) ){
                ioStoreMatrix(f1.f,intracellularSelfEwald ,1,"intracellularSelfEwald.1.matrix",0);
            }

            if ( CanonicalRank(f1.f,intercellularSelfEwald,0) ){
                ioStoreMatrix(f1.f,intercellularSelfEwald ,0,"intercellularSelfEwald.matrix",0);
            }
            if ( CanonicalRank(f1.f,intercellularSelfEwald,1) ){
                ioStoreMatrix(f1.f,intercellularSelfEwald ,1,"intercellularSelfEwald.1.matrix",0);
            }
#endif
//            tClear(f1.f, interactionExchange);
//            tClear(f1.f, interactionEwald);

        }
    
    
    countHam(&c,f1);

    if ( allowQ(f1.f.rt, blockTrainHamiltonianBlock) && allowQ(f1.f.rt, blockHamiltonianBlock)&& allowQ(f1.f.rt, blockAllocateHamiltonianBlock)){
    enum division di;
    INT_TYPE spin,space ;
    for ( spin = 0 ; spin < f1.f.cmpl ; spin++){
        tClear(f1.f,hamiltonian);
        for ( di = Ha ; di!= nullName; di = f1.f.tulip[di].linkNext){
            if ( CanonicalRank(f1.f, f1.f.tulip[di].name, spin) ){
                if ( CanonicalRank(f1.f, hamiltonian, 0)+CanonicalRank(f1.f, f1.f.tulip[di].name, spin)  > part(f1.f, hamiltonian )){
                    printf("part sum\n");
                    exit(0);
                }
                        {//symmetric hamiltonian.
                            switch( bodies(f1.f, f1.f.tulip[di].name) ){
                                case one :
                                    switch ( f1.i.body ){
                                        case 1:
                                            tAddTw(f1.f, hamiltonian, 0, f1.f.tulip[di].name, spin);
                                            printf("add1 %d %d\n", f1.f.tulip[di].name,spin);
                                            break;
                                        case 2:
                                            for (space = 0 ; space < SPACE ; space++)
                                            if ( f1.f.rose[space].body != nada )
                                            sumTo2(f1.f,1.,space,f1.f.tulip[di].space[space].block, f1.f.tulip[di].name, spin, hamiltonian, 0);
                                            f1.f.tulip[hamiltonian].Current[0] += CanonicalRank(f1.f, f1.f.tulip[di].name, spin);

                                            printf("add1 %d %d\n", f1.f.tulip[di].name,spin);
                                            break;
                                        case 3:
                                            for (space = 0 ; space < SPACE ; space++)
                                            if ( f1.f.rose[space].body != nada )
                                            sumTo2(f1.f,1./2.,space,f1.f.tulip[di].space[space].block, f1.f.tulip[di].name, spin, hamiltonian, 0);
                                            f1.f.tulip[hamiltonian].Current[0] += CanonicalRank(f1.f, f1.f.tulip[di].name, spin);

                                            printf("add1 %d %d\n", f1.f.tulip[di].name,spin);
                                            break;
                                        case 4:
                                            for (space = 0 ; space < SPACE ; space++)
                                            if ( f1.f.rose[space].body != nada )
                                            sumTo2(f1.f,1./3.,space,f1.f.tulip[di].space[space].block, f1.f.tulip[di].name, spin, hamiltonian, 0);
                                            f1.f.tulip[hamiltonian].Current[0] += CanonicalRank(f1.f, f1.f.tulip[di].name, spin);

                                            printf("add1 %d %d\n", f1.f.tulip[di].name,spin);
                                            break;
                                    }
                                    break;
                                case two:
                                    switch ( f1.i.body ){//symmetric hamiltonian.
                                        case 2:
                                        case 3:
                                        case 4:
                                        case 5:
                                        case 6:
                                        {
                                            tAddTw(f1.f, hamiltonian, 0, f1.f.tulip[di].name, spin);
                                            printf("add2 %d %d\n", f1.f.tulip[di].name,spin);
                                        }
                                        case 1:
                                            break;
                                    }
                                break;
                            }
                }
                printf("<%f>%d\n", traceOne(f1.f, hamiltonian, 0),CanonicalRank(f1.f, hamiltonian, 0));
            }
        }
        tCycleDecompostionGridOneMP(-1, f1.f, hamiltonian, 0, NULL,trainHamiltonian  , spin, c.rt.CANON, part(f1.f,trainHamiltonian), c.rt.powDecompose);
    }

        
        
    }

    
    if ( allowQ(f1.f.rt, blockHamiltonianBlock)&&allowQ(f1.f.rt, blockTrainHamiltonianBlock))
    {

        if ( CanonicalRank(f1.f,trainHamiltonian,0) ){
            ioStoreMatrix(f1.f,trainHamiltonian ,0,"trainHamiltonian.matrix",0);
        }
        if ( CanonicalRank(f1.f,trainHamiltonian,1) ){
            ioStoreMatrix(f1.f,trainHamiltonian ,1,"trainHamiltonian.1.matrix",0);
        }
        
    }
        fModel(&f1.f);
    return 0;
}



int run (INT_TYPE argc , char * argv[]){
    argc--;//erase ./andromeda...
    argv++;
    struct calculation c;
    struct field f;
    
    
    
    if ( argc > 0 ){
        switch ( atoi( argv[0])){
            case -1 :
                argc--;
                argv++;
                //andromeda -1 inputFile
                //runs normally from file
                bootShell(argc, argv,&c,&f);
                break;

            case 1 :
                argc--;
                argv++;
                //andromeda 1 inputFile
                //spits out memory requirements only
                bootShell(argc, argv,&c,&f);
                c.i.RAMmax = 0;
                break;

            case 0 :
                //andromeda 0
                printf("----\nv7.7\n\n%s\n\n",getenv("LAUNCH"));
                exit(0);
        }

    }else {
        //andromeda inputFile
        bootShell(argc, argv,&c,&f);

   // test2();

    }
    
    INT_TYPE space,i,a,plusSize,nStatesTrans=0,nStatesFound=0 ,RdsSize = 0,totalIter = 0;
    FILE * out = stdout;
    struct runTime * rt = & c.rt;
    f.f.rt = rt;
#ifdef OMP
    if ( c.i.omp > MaxCore ){
        printf("lanes > MaxCore\n");
        c.i.omp = MaxCore;
    }
    if ( c.i.omp == -1 ){
#ifdef MKL
        if ( c.i.mkl < 1 )
        {
            printf("set parallel");
            exit(0);
        }
        
        c.i.omp = MaxCore/c.i.mkl;
#else
        c.i.omp = MaxCore;
#endif
    }

    rt->NLanes = c.i.omp;
#pragma omp parallel for private (i)
    for ( i = 0; i < MaxCore ; i++){
        rt->NSlot = omp_get_num_threads();
    }
    if ( rt->NLanes > rt->NSlot ){
        printf("decrease lanes or increase your number of OMP cores\n");
        rt->NLanes = rt->NSlot;
    }
    
#ifdef MKL
    if ( rt->NSlot < c.i.omp*c.i.mkl )
    {
        printf("not enough slots for omp\n" );
        c.i.omp = rt->NSlot/c.i.mkl;
    }
    rt->NParallel = c.i.mkl;
    printf("parallel \t %d\n", rt->NParallel);

#endif
    printf("lanes \t %d\n",  rt->NLanes);

#endif
    //0//...   //A//B//C//D//E
    if ( c.rt.phaseType == buildFoundation ){//0
#ifdef NBODY
        foundationM(&c,f);
#else
        foundation1(&c,f);
#endif
    }
    else if ( c.rt.phaseType == productKrylov ){//C
        krylov(&c,f);
    }
    else if ( c.rt.phaseType == solveRitz ){//A
        ritz(&c,f);
    }
    else if ( c.rt.phaseType == svdOperation ){//E
       // svd(c,f);
    }else if ( c.rt.phaseType == distillMatrix ){
        distill(c,f);
    } else if ( c.rt.phaseType == reportMatrix ){
        iModel(&c, &f);
        report(c,f);
        fModel(&f.f);
    }else if ( c.rt.phaseType == decomposeMatrix ){
        decompose ( &c, f);
    }else if ( c.rt.phaseType == gauss){
        spitGauss ( &c, f);
    }
    printf("\n\nFINIS.\n\n");
    return 0;
}
#if 1
    int main (INT_TYPE argc , char * argv[]){
        
#if 0
#ifdef APPLE
        for ( int n = 0 ; n < 10 ;n++){
            printf("%d\n",n);
        printf("%f\n",momentumIntegralInTrain(2, n, 1,eikonDiagonal, two));
        printf("%f\n",momentumIntegralInTrain(2, n, 1,eikonSemiDiagonal, two));
        printf("%f\n",momentumIntegralInTrain(2, n, 1,eikonOffDiagonal, two));
        }
        return 0;
#endif
#else
        return run(argc,argv);
#endif
    }
#else

int main (INT_TYPE argc , char * argv[]){
    struct general_index g ;
    g.bra.type = 1;
    g.bra.basis = GaussianBasisElement;
    g.bra.length = 2;
    g.bra.origin = 0;
    g.bra.index = 0;
    g.bra.index2 = 0;
    g.bra.grid = 9;

    
    g.ket.type = 1;
    g.ket.basis = SincBasisElement;
    g.ket.length = 1;
    g.ket.origin = 0.;
    g.ket.index = 0;
    g.ket.index2 = 0;
    g.ket.grid = 9;
    printf("%f", creal(BoB(g.bra,g.ket)));
    
    g.bra.type = 1;
    g.bra.basis = GaussianBasisElement;
    g.bra.length = 1;
    g.bra.origin = 0;
    g.bra.index = 0;
    g.bra.index2 = 0;
    g.bra.grid = 9;

    
    g.ket.type = 1;
    g.ket.basis = SincBasisElement;
    g.ket.length = 1;
    g.ket.origin = 0.;
    g.ket.index = 0;
    g.ket.index2 = 0;
    g.ket.grid = 9;
    printf("%f", cimag(Bd2B(g.bra,g.ket)));
    g.bra.type = 1;
    g.bra.basis = GaussianBasisElement;
    g.bra.length = 1;
    g.bra.origin = 0;
    g.bra.index = 0;
    g.bra.index2 = 0;
    g.bra.grid = 9;

    
    g.ket.type = 1;
    g.ket.basis = SincBasisElement;
    g.ket.length = 1;
    g.ket.origin = 0.;
    g.ket.index = 1;
    g.ket.index2 = 0;
    g.ket.grid = 9;
    printf("%f", cimag(Bd2B(g.ket,g.bra)));

    
    g.bra.type = 1;
    g.bra.basis = GaussianBasisElement;
    g.bra.length = 1;
    g.bra.origin = 0;
    g.bra.index = 0;
    g.bra.index2 = 0;
    g.bra.grid = 9;

    
    g.ket.type = 1;
    g.ket.basis = GaussianBasisElement;
    g.ket.length = 1;
    g.ket.origin = 0.;
    g.ket.index = 0;
    g.ket.index2 = 0;
    g.ket.grid = 9;
    printf("%f", cimag(Bd2B(g.bra,g.ket)));

    //printf("%f %f %f",creal( FGS(-2,2,&g)),cimag( FGS(-2,2,&g)),cabs( FGS(-2,2,&g)));

}
#endif
