/*
 *  jobs.c
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

#include "jobs.h"

INT_TYPE foundation1(struct calculation *c1, struct field f1){
    INT_TYPE EV;
    f1.i.Iterations = 1;
    if ( 1 ){
        iModel(c1,&f1);
        separateKinetic(f1.f, 0,kinetic, 1,electron);
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
    INT_TYPE EV;
    f1.i.Iterations = 1;
    if ( 1 ){
        iModel(c1,&f1);
        ioStoreMatrix(f1.f, quadCube, 0, "single.matrix", 1);
        switch(bodies(f1.f,eigen)){
            case one:
                tEqua(f1.f, eigen, 0, quadCube, 0);
                break;
            case two:
                tEqua(f1.f, eigen, 0, quadCube, 0);
                break;
            case three:
                sumTo3(f1.f, quadCube, 0,eigen, 0);
                tId(f1.f , eigen,0);
                tCycleDecompostionGridOneMP(-2, f1.f, build, 0, NULL,eigen , 0, c1->rt.vCANON, 1, 0);
                break;
            case four:
                sumTo4(f1.f, quadCube, 0,eigen, 0);
                tId(f1.f , eigen,0);
                tCycleDecompostionGridOneMP(-2, f1.f, build, 0, NULL,eigen , 0, c1->rt.vCANON, 1, 0);
                break;
            }
        tBootManyConstruction(c1,f1.f ,eigen);
        if ( ! tSortBoot(c1,f1.f,build) )
            EV =   tSlam(f1.f,f1.i.qFloor,f1.f.user,c1->i.level);
        else
            EV = 0;
        
    }
    fModel(&f1.f);

    return EV;
}



INT_TYPE krylov ( struct calculation *c1, struct field f1){
    INT_TYPE EV = 0,i,fi,next;


    f1.i.qFloor = countLinesFromFile(c1,f1,0,&f1.i.iRank, &f1.i.xRank);
    //count canonical-rank...

    
    f1.i.nStates =1  ;
    f1.i.nStates =f1.i.Iterations  ;
   
 iModel(c1,&f1);
    for ( fi =0 ; fi < f1.i.files ; fi++){
        tLoadEigenWeights (c1,f1, f1.i.fileList[fi], &EV,f1.f.user , f1.i.collect);
    }
        if (EV == 0 ){
        print(c1,f1,1,0,0,eigenVectors);
        return 1;
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
        
        
        
        if ( ! tGreatDivideIteration(flag, c1->i.shiftVector[iterator-1][0],c1->i.shiftVector[iterator-1][1],  f1.f,Iterator, 1,0,eigenVectors+RdsSize-EV,EV,2*EV,0)){
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
    inputFormat(f1.f, c1->name, hamiltonian, 1);
    
    tClear(f1.f,trainHamiltonian);
    
    tCycleDecompostionGridOneMP(-2, f1.f, hamiltonian, 0, NULL,trainHamiltonian , 0, c1->rt.CANON, part(f1.f,trainHamiltonian), c1->rt.powDecompose);

    FILE * fo = fopen("output.matrix", "w");
    outputFormat(f1.f, fo, trainHamiltonian, 0);
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
    

    for ( fi =0 ; fi < f1.i.files ; fi++)
        tLoadEigenWeights (c1,f1, f1.i.fileList[fi],&EV, f1.f.user,f1.i.collect);//UNUSUAL!!!
    if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
    }
    
    {
        INT_TYPE typer;
        if ( c1->i.shiftFlag )
            typer = -1;
        else
            typer = 1;
    
    tEigenCycle(typer,f1.f,Ha,CDT, f1.i.nStates, f1.f.user,EV,0, EV,0,1,eigenVectors,twoBodyRitz);
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
        for ( ii = 0; ii < EV ; ii++)  printVector(c1,f1.f,f1.f.tulip[f1.f.user+ii].value.title,str,f1.f.tulip[f1.f.user+ii].value.stage-1,f1.i.irrep, V+stride*iii+ii);
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

INT_TYPE sumTo2(struct sinc_label f1, enum division mat,INT_TYPE ms, enum division sum,INT_TYPE spin){
    
    INT_TYPE I1,I2, I3, I4,body,r,space,N1;
    double value;
    for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
        for ( body = 0 ; body < 2 ; body++){
            for ( space = 0; space < SPACE ; space++){
                N1 = vector1Len(f1, space);
                Stream_Type * stream = streams(f1,sum,spin,space)+N1*N1*N1*N1*CanonicalRank(f1, sum, spin);
                if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
                    printf("part sum\n");
                    exit(0);
                }
                for ( I1 = 0 ; I1 < N1 ; I1++)//body 0
                    for ( I2 = 0 ; I2 < N1 ; I2++)//body 0
                        for ( I3 = 0 ; I3 < N1 ; I3++)//body 1
                            for ( I4 = 0 ; I4 < N1 ; I4++)//body 1
                            {
                                if ( body == 0 ){
                                    value  = streams(f1,mat,ms,space)[ I1*N1+I2 + r*N1*N1 ] * delta(I3-I4);
                                }else {
                                    value  = streams(f1,mat,ms,space)[ I3*N1+I4 + r*N1*N1 ] * delta(I1-I2);
                                }
                                stream[ (I1+I3*N1)+ ( I2+I4*N1)*N1*N1 ] = value;
                            }
            }
            f1.tulip[sum].Current[spin]++;
            
            //  tAddTwo(f1, sum , quadCube);
        }
    return 0;
}

INT_TYPE sumTo3(struct sinc_label f1, enum division mat,INT_TYPE ms, enum division sum,INT_TYPE spin){

    INT_TYPE n2[SPACE];
    length(f1, sum,n2);

    if ( bodies ( f1, sum ) == three && f1.rt->calcType == electronicStuctureCalculation){

        if (bodies(f1, mat ) == three ){
            tAddTw(f1, sum ,spin, mat,ms );
        }else if ( bodies ( f1, mat ) == one ){
            INT_TYPE I1,I2, I3, I4,I5,I6,body,r,space;
            INT_TYPE n1[SPACE];
            length1(f1,n1);
            double value;

            for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
                for ( body = 0 ; body < 3 ; body++){
                    for ( space = 0; space < SPACE ; space++){
                        Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*CanonicalRank(f1, sum, spin);
                        if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
                            printf("part sum\n");
                            exit(0);
                        }

                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)

                                            {
                                                value = 0;
                                                if ( body == 0 ){
                                                    value  = streams(f1,mat,ms,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6);
                                                }else if ( body == 1 ) {
                                                    value  = streams(f1,mat,ms,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6);
                                                }else if ( body == 2 ) {
                                                    value  = streams(f1,mat,ms,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4);
                                                }
                                                stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]] = value;
                                            }
                    }
                    f1.tulip[sum].Current[spin]++;

                }
        }else if ( bodies ( f1, mat ) == two ){
            INT_TYPE I1,I2, I3, I4,I5,I6,pair,r,space,ve;
            INT_TYPE n1[SPACE];

            length1(f1, n1);

            double value;

            for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
                for ( pair = 0 ; pair < 3 ; pair++){
                    for ( space = 0; space < SPACE ; space++){
                        Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*CanonicalRank(f1, sum, spin);
                        if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
                            printf("part sum\n");
                            exit(0);
                        }
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)

                                            {
                                                value = 0;
                                                if ( pair == 0 ){
                                                    value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6);
                                                }else if ( pair == 1 ) {
                                                    value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4);
                                                }else if ( pair == 2 ) {
                                                    value  = streams(f1,mat,ms,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2);
                                                }
                                                ve = (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space];
                                                //                                                printf("%f %lld %lld\n", value,ve,n2[space]);
                                                //                                                fflush(stdout);
                                                stream[ ve ] = value;
                                                //                                                printf("x");
                                                //                                                fflush(stdout);

                                            }
                    }
                    f1.tulip[sum].Current[spin]++;
                }
        }
        else {
            printf("Yo!");
            exit(0);
        }

    }
    return 0;
}








INT_TYPE sumTo4(struct sinc_label f1, enum division mat,INT_TYPE ms, enum division sum,INT_TYPE spin){

    INT_TYPE n2[SPACE];
    length(f1, sum,n2);
    if (bodies(f1,sum) == four && f1.rt->calcType == electronicStuctureCalculation){

    if ( bodies ( f1, mat ) == one ){
        INT_TYPE I1,I2, I3, I4,I5,I6,I7,I8,body,r,space;
        INT_TYPE n1[SPACE];
        length1(f1,n1);
        double value;

        for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
            for ( body = 0 ; body < 4 ; body++){
                for ( space = 0; space < SPACE ; space++){
                    Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*CanonicalRank(f1, sum, spin);
                    if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
                        printf("part sum\n");
                        exit(0);
                    }
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
                                                    if ( body == 0 ){
                                                        value  = streams(f1,mat,ms,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6)*delta(I7-I8);
                                                    }else if ( body == 1 ) {
                                                        value  = streams(f1,mat,ms,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6)*delta(I7-I8);
                                                    }else if ( body == 2 ) {
                                                        value  = streams(f1,mat,ms,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I7-I8);
                                                    }else if ( body == 3 ) {
                                                        value  = streams(f1,mat,ms,space)[ I7*n1[space]+I8 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I5-I6);
                                                    }
                                                    stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
                                                }
                }
                f1.tulip[sum].Current[spin]++;

            }
    }else if ( bodies ( f1, mat ) == two ){
        INT_TYPE I1,I2, I3, I4,I5,I6,I7,I8,pair,r,space;
        INT_TYPE n1[SPACE];
        length1(f1, n1);

        double value;

        for ( r = 0; r < CanonicalRank(f1, mat , ms ); r++)
            for ( pair = 0 ; pair < 6 ; pair++){
                for ( space = 0; space < SPACE ; space++){
                    Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*CanonicalRank(f1, sum, spin);
                    if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
                        printf("part sum\n");
                        exit(0);
                    }
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


                                                    if ( pair == 0 ){
                                                        value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6)*delta(I7-I8);
                                                    }else if ( pair == 1 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I7-I8);
                                                    }else if ( pair == 2 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I7-I8);
                                                    }else if ( pair == 3 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I1*n1[space]+I7) + (I2*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I5-I6);
                                                    }else if ( pair == 4 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I3*n1[space]+I7) + (I4*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I5-I6);
                                                    }else if ( pair == 5 ) {
                                                        value  = streams(f1,mat,ms,space)[ (I5*n1[space]+I7) + (I6*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I3-I4);
                                                    }else {
                                                        printf ("rails!\n");
                                                        exit(0);
                                                    }
                                                    stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
                                                }
                }
                f1.tulip[sum].Current[spin]++;
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
        ioStoreMatrix(f1.f,interactionExchange ,0,"interactionExchange.matrix",0);
    }
    if ( CanonicalRank(f1.f,interactionExchange,1) ){
        ioStoreMatrix(f1.f,interactionExchange ,0,"interactionExchange.1.matrix",0);
    }


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


    if ( CanonicalRank(f1.f,interactionEwald,0) ){
        ioStoreMatrix(f1.f,interactionEwald ,0,"interactionEwald.matrix",0);
    }
    if ( CanonicalRank(f1.f,interactionEwald,1) ){
        ioStoreMatrix(f1.f,interactionEwald ,1,"interactionEwald.1.matrix",0);
    }

    if ( CanonicalRank(f1.f,shortenPlus,0) )
        ioStoreMatrix(f1.f,shortenPlus ,0,"shortenExchangePlus.matrix",0);

    if ( CanonicalRank(f1.f,shortenMinus,0) )
        ioStoreMatrix(f1.f,shortenMinus ,0,"shortenExchangeMinus.matrix",0);


    if ( CanonicalRank(f1.f,linear,0) ){
        ioStoreMatrix(f1.f,linear ,0,"linear.matrix",0);

    }
    if ( CanonicalRank(f1.f,linear,1) ){
        ioStoreMatrix(f1.f,linear ,1,"linear.1.matrix",0);

    }
    

    return 0;
}



INT_TYPE distill ( struct calculation c, struct field f1){
    double oneBodyFraction = 1.;
    switch( f1.i.body  ){
        case two:
            oneBodyFraction = 1/1.;
            break;
        case three:
            oneBodyFraction = 1/2.;
            break;
        case four:
            oneBodyFraction = 1/3.;
            break;
    }
    iModel(&c, &f1);
    tClear(f1.f, hamiltonian);
    
    
    if ( allowQ(f1.f.rt, blockTrainHamiltonianBlock) && allowQ(f1.f.rt, blockHamiltonianBlock)&& allowQ(f1.f.rt, blockTrainingHamiltonianBlock)){
        
        
        if ( f1.i.body >= two ){
            
            
            if ( c.rt.runFlag == 0 ){
                ioStoreMatrix(f1.f, hamiltonian, 0, "interactionExchange.matrix", 1);
                
                tScaleOne(f1.f, kinetic, 0, oneBodyFraction);
                sumTo2(f1.f, kinetic, 0, hamiltonian, 0);
                tScaleOne(f1.f, kinetic, 0, 1./oneBodyFraction);
                                
                tCycleDecompostionGridOneMP(-1, f1.f, hamiltonian, 0, NULL,trainHamiltonian  , 0, c.rt.CANON, part(f1.f,trainHamiltonian), c.rt.powDecompose);
                tClear(f1.f,hamiltonian);
                sortTerms(f1.f,trainHamiltonian,0,hamiltonian,0);
                tEqua(f1.f, trainHamiltonian,0, hamiltonian, 0);
                
            } else {
                ioStoreMatrix(f1.f, hamiltonian, 0, "interactionEwald.matrix", 1);

                tScaleOne(f1.f, linear, 0, oneBodyFraction);
                sumTo2(f1.f, linear, 0, hamiltonian, 0);
                tScaleOne(f1.f, linear, 0, 1./oneBodyFraction);
                
                tScaleOne(f1.f,intercellularSelfEwald, 0, oneBodyFraction);
                sumTo2(f1.f, intercellularSelfEwald, 0, hamiltonian, 0);
                tScaleOne(f1.f,intercellularSelfEwald, 0, 1./oneBodyFraction);
                
                tScaleOne(f1.f,intracellularSelfEwald, 0, oneBodyFraction);
                sumTo2(f1.f, intracellularSelfEwald, 0, hamiltonian, 0);
                tScaleOne(f1.f,intracellularSelfEwald, 0, 1./oneBodyFraction);
                
                tScaleOne(f1.f,jelliumElectron, 0, oneBodyFraction);
                sumTo2(f1.f, jelliumElectron, 0, hamiltonian, 0);
                tScaleOne(f1.f,jelliumElectron, 0, 1./oneBodyFraction);
                
                tScaleOne(f1.f, kinetic, 0, oneBodyFraction);
                sumTo2(f1.f, kinetic, 0, hamiltonian, 0);
                tScaleOne(f1.f, kinetic, 0, 1./oneBodyFraction);
                
                tCycleDecompostionGridOneMP(-2, f1.f, hamiltonian, 0, NULL,trainHamiltonian  , 0, c.rt.CANON, part(f1.f,trainHamiltonian), c.rt.powDecompose);
                tClear(f1.f,hamiltonian);
                sortTerms(f1.f,trainHamiltonian,0,hamiltonian,0);
                tEqua(f1.f, trainHamiltonian,0, hamiltonian, 0);
                
                tClear(f1.f,hamiltonian);
                
                //
                //**
                //
                ioStoreMatrix(f1.f, hamiltonian, 0, "interactionEwald.1.matrix", 1);
                tScaleOne(f1.f, kinetic, 1, oneBodyFraction);
                sumTo2(f1.f, kinetic, 1, hamiltonian, 0);
                tScaleOne(f1.f, kinetic, 1, 1./oneBodyFraction);
                
                tScaleOne(f1.f, linear, 1, oneBodyFraction);
                sumTo2(f1.f, linear, 1, hamiltonian, 0);
                tScaleOne(f1.f, linear, 1, 1./oneBodyFraction);
                
                tScaleOne(f1.f,intercellularSelfEwald, 1, oneBodyFraction);
                sumTo2(f1.f, intercellularSelfEwald, 1, hamiltonian, 0);
                tScaleOne(f1.f,intercellularSelfEwald, 1, 1./oneBodyFraction);
                
                tScaleOne(f1.f,intracellularSelfEwald, 1, oneBodyFraction);
                sumTo2(f1.f, intracellularSelfEwald, 1, hamiltonian, 0);
                tScaleOne(f1.f,intracellularSelfEwald, 1, 1./oneBodyFraction);
                
                tScaleOne(f1.f,jelliumElectron, 1, oneBodyFraction);
                sumTo2(f1.f, jelliumElectron, 1, hamiltonian, 0);
                tScaleOne(f1.f,jelliumElectron, 1, 1./oneBodyFraction);
                
                tCycleDecompostionGridOneMP(-1, f1.f, hamiltonian, 0, NULL,trainHamiltonian  , 1, c.rt.CANON, part(f1.f,trainHamiltonian), c.rt.powDecompose);
                tClear(f1.f,hamiltonian);
                sortTerms(f1.f,trainHamiltonian,1,hamiltonian,0);
                tEqua(f1.f, trainHamiltonian,1, hamiltonian, 0);
            }
        }else {
            if ( c.rt.runFlag == 0 ){
                
                if ( c.rt.calcType == electronicStuctureCalculation){
                    tAddTw(f1.f,hamiltonian,0,kinetic ,0);
                } else {
                    tAddTw(f1.f,hamiltonian,0,shortenPlus ,0);
                    tAddTw(f1.f,hamiltonian,0,shortenMinus ,0);
                    tAddTw(f1.f,hamiltonian,0,kinetic ,0);
                    tAddTw(f1.f,hamiltonian,0,kineticMass ,0);
                    tAddTw(f1.f,hamiltonian,0,protonRepulsion ,0);

                    
                    
                }
                tCycleDecompostionGridOneMP(-1, f1.f, hamiltonian, 0, NULL,trainHamiltonian  , 0, c.rt.CANON, part(f1.f,trainHamiltonian), c.rt.powDecompose);
                
                
                tClear(f1.f,hamiltonian);
                sortTerms(f1.f,trainHamiltonian,0,hamiltonian,0);
                tEqua(f1.f, trainHamiltonian,0, hamiltonian, 0);
                
            } else {

                tAddTw(f1.f,hamiltonian,0,kinetic ,0);
                tAddTw(f1.f,hamiltonian,0,intercellularSelfEwald ,0);
                tAddTw(f1.f,hamiltonian,0,intracellularSelfEwald ,0);
                tAddTw(f1.f,hamiltonian,0,jelliumElectron ,0);
                
                tCycleDecompostionGridOneMP(-1, f1.f, hamiltonian, 0, NULL,trainHamiltonian  , 0, c.rt.CANON, part(f1.f,trainHamiltonian), c.rt.powDecompose);
              //  printf("train : %d %d\n",CanonicalRank(f1.f,trainHamiltonian,0),CanonicalRank(f1.f,trainHamiltonian,1) );

                tClear(f1.f,hamiltonian);
                tAddTw(f1.f,hamiltonian,0,kinetic ,1);
                tAddTw(f1.f,hamiltonian,0,intercellularSelfEwald ,1);
                tAddTw(f1.f,hamiltonian,0,intracellularSelfEwald ,1);
                tAddTw(f1.f,hamiltonian,0,jelliumElectron ,1);
                
                tCycleDecompostionGridOneMP(-1, f1.f, hamiltonian, 0, NULL,trainHamiltonian  , 1, c.rt.CANON, part(f1.f,trainHamiltonian), c.rt.powDecompose);
               // printf("train : %d %d\n",CanonicalRank(f1.f,trainHamiltonian,0),CanonicalRank(f1.f,trainHamiltonian,1) );


                
            }
        }
    }
    if ( allowQ(f1.f.rt, blockHamiltonianBlock))
    {
        if ( CanonicalRank(f1.f,trainHamiltonian,0) ){
            ioStoreMatrix(f1.f,trainHamiltonian ,0,"trainHamiltonian.matrix",0);
            tId(f1.f , quadCube,0);
            tCycleDecompostionGridOneMP(-1, f1.f, trainHamiltonian, 0, NULL,quadCube  , 0, c.rt.CANON, 1, 0);
            ioStoreMatrix(f1.f,quadCube ,0,"single.matrix",0);

        }
        if ( CanonicalRank(f1.f,trainHamiltonian,1) ){
            ioStoreMatrix(f1.f,trainHamiltonian ,1,"trainHamiltonian.1.matrix",0);
            tId(f1.f , quadCube,0);
            tCycleDecompostionGridOneMP(-1, f1.f, trainHamiltonian, 1, NULL,quadCube  , 0, c.rt.CANON, 1,0);
            ioStoreMatrix(f1.f,quadCube ,0,"single.1.matrix",0);
        }
        
    }
        fModel(&f1.f);
    return 0;
}


#if 1

int main (INT_TYPE argc , char * argv[]){
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
                printf("----\nv7.5.4\n\n%s\n\n",getenv("LAUNCH"));
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
    
    INT_TYPE flag = 0;
#ifdef OMP
    if ( c.i.omp > MaxCore ){
        printf("lanes > MaxCore\n");
        c.i.omp = MaxCore;
        flag = 1;
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
        flag = 1;
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
        flag =1 ;
    }
    
#ifdef MKL
    if ( flag )
        c.i.mkl = 1;
    
    if ( rt->NSlot < c.i.mkl * c.i.omp )
    {
        printf("not enough slots for mkl*omp\n" );
        c.i.mkl = rt->NSlot/c.i.omp;
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
    }
    printf("\n\nFINIS.\n\n");
}

#endif
