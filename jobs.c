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

INT_TYPE foundation(struct calculation *c1, struct field f1){
    INT_TYPE EV;

    if ( 1 ){
        c1->i.irrep = 0;
        c1->i.twoBody.func.fn = nullFunction;
        iModel(c1,&f1);
        tBoot1Construction(c1,f1.f ,build);

//        if ( f1.i.irrep && f1.i.filter )
//            EV =   tCollect(f1.f,f1.i.irrep,f1.f.user,f1.i.qFloor ,1);
//        else
            EV =   tSlam(f1.f,f1.i.qFloor,f1.f.user,c1->i.level);
        tFilter(f1.f, EV,0, f1.f.user);//classify
    }else {
        exit(1);
    }
    
    INT_TYPE ii,flag = 1;;
    if ( c1->i.irrep )
        for ( ii = 0 ;ii < EV ;ii++){
            if ( f1.f.tulip[f1.f.user+ii].value.symmetry == f1.i.irrep ){
                print(c1,f1,flag,ii,ii+1 , f1.f.user);
                flag = 0;
            }
        }
    else
        print(c1,f1,1,0,EV , f1.f.user);
    
    fModel(&f1.f);
//    else{
//        iModel(c1);
//
////tBootManyConstruction(c1);
//        EV =   tCollect(f1,c1->i.irrep,f1.user,c1->i.qFloor ,1);
//        tFilter(f1, EV, 0, f1.user);//classify
//    }
    return EV;
}


INT_TYPE krylov ( struct calculation *c1, struct field f1){
    INT_TYPE EV = 0,i,fi;


    f1.i.qFloor = countLinesFromFile(f1,0,&f1.i.iRank);
    //count canonical-rank...

    
    f1.i.nStates =f1.i.qFloor* f1.i.Iterations  ;
    iModel(c1,&f1);
    for ( fi =0 ; fi < f1.i.files ; fi++){
        EV +=  tLoadEigenWeights (c1,f1, f1.i.fileList[fi], f1.f.user );//UNUSUAL!!!
    }
        if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
    }
    INT_TYPE RdsSize = EV,iterator=0;

    if(1){
        printf ("Step \t%d\n", iterator);
        INT_TYPE iii,sp ;
        for ( iii = 0; iii < EV ; iii++){
            printf ( "\n Vector \t%d \n", iii+1);
            for ( sp = 0 ; sp < spins(f1.f, f1.f.user); sp++)
                tCycleDecompostionGridOneMP(-1, f1.f, f1.f.user+iii, sp, NULL, eigenVectors+iii, sp, c1->rt.vCANON, part(f1.f,eigenVectors+iii), -1);
            
            tFilter(f1.f, EV, 0, eigenVectors+RdsSize-EV);//classify
            printExpectationValues(f1.f, Ha, eigenVectors+RdsSize-EV+iii);
            fflush(stdout);
            print(c1,f1,1,RdsSize-EV+iii,RdsSize-EV+iii+1,eigenVectors);
        }
    }

    
    for ( iterator = 1 ; iterator < f1.i.Iterations ; iterator++){
        RdsSize += tGreatDivideIteration(c1->i.shiftFlag, c1->i.realPart,  f1.f,Iterator, 1,0,eigenVectors+RdsSize-EV,EV,2*EV,0)-EV;
        if(1){
            tFilter(f1.f, EV, !(!c1->i.filter )* c1->i.irrep, eigenVectors+RdsSize-EV);//filter
            printf ("Step \t%d\n", iterator);
            INT_TYPE iii ;
            for ( iii = 0; iii < EV ; iii++){
                printf ( "\n Vector \t%d \n", iii+1);
                printExpectationValues(f1.f, Ha, eigenVectors+RdsSize-EV+iii);
                fflush(stdout);
                print(c1,f1,0,RdsSize-EV+iii,RdsSize-EV+iii+1,eigenVectors);
            }
        }
        fflush(stdout);
    }
    fModel(&f1.f);

    return 0;
}

INT_TYPE ritz( struct calculation * c1, struct field f1){
    size_t ms = MAXSTRING;
    char line0[MAXSTRING];
    char filename[MAXSTRING];    char str[MAXSTRING];


    char * line = line0;
    c1->i.canonRank = 0;
    INT_TYPE fi,EV = 0,i,j,sp;
    double va;
    INT_TYPE lines = 0,flines = 0;

    f1.i.qFloor = countLinesFromFile(f1,0,&f1.i.iRank);

   // c1->i.iRank = c1->i.bRank;
    iModel(c1,&f1);
    
    if ( CanonicalRank(f1.f,interactionExchange,0) )
        ioStoreMatrix(f1.f,interactionExchange ,0,"interactionExchange.matrix",0);
    if ( CanonicalRank(f1.f,interactionEwald,0) )
        ioStoreMatrix(f1.f,interactionEwald ,0,"interactionEwald.matrix",0);

    
        if ( CanonicalRank(f1.f,shortenPlus,0) )
            ioStoreMatrix(f1.f,shortenPlus ,0,"shortenExchangePlus.matrix",0);
        
        if ( CanonicalRank(f1.f,shortenMinus,0) )
            ioStoreMatrix(f1.f,shortenMinus ,0,"shortenExchangeMinus.matrix",0);

    
    if ( CanonicalRank(f1.f,linear,0) )
        ioStoreMatrix(f1.f,linear ,0,"linear.matrix",0);

    
    
    
#ifndef APPLE
    for ( fi =0 ; fi < f1.i.files ; fi++)
        EV +=  tLoadEigenWeights (c1,f1, f1.i.fileList[fi], f1.f.user+EV);//UNUSUAL!!!
    if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
    }
#endif
    tEigenCycle(f1.f,Ha,CDT, f1.i.nStates, f1.f.user,EV,0, EV,0,1,eigenVectors,twoBodyRitz);
    
    
    DCOMPLEX *V = (DCOMPLEX*)myStreams(f1.f,matrixHbuild,0);
    INT_TYPE iii,stride = f1.f.maxEV;
    char * token = filename;
    for ( iii = 0; iii < f1.i.nStates ; iii++){
        
        {
            FILE * outf ;
            sprintf(str, "%s-%d.vector",c1->name,iii+1);
            outf = fopen (str,"w");
            fclose(outf);
        }
        lines = 0;

        for ( fi =0 ; fi <f1.i.files ; fi++){
            flines = 0;
            FILE * fp = fopen(f1.i.fileList[fi],"r");
            if ( fp == NULL ) {
                printf("file?\n");
                exit(0);
            }
            getline(&line, &ms, fp);

            while(!feof(fp))
            {
                if (! comment(line))
                {
                    
                    token = strtok(line, "\"");
                    
                    sprintf(str, "%s-%d", c1->name, iii+1);
                    
                    {
                        INT_TYPE number,si;
                        char pa[MAXSTRING],*pa0;
                        pa0 = pa;
                        pa0 = strtok(NULL, "\"");
                        si = sscanf ( pa0, ",%d,",&number );
                        if ( si == 1  ) {
                            printVector(c1,f1.f,token,str,number-1,f1.i.irrep, V+stride*iii+lines);
                        }
                    }
                    lines++;
                }
                flines++;

                getline(&line, &ms, fp);

            }
            
            fclose(fp);
        }
    }
    fModel(&f1.f);

    return 0;
}

INT_TYPE decompose( struct calculation * c1 , struct field f1){

    INT_TYPE i,g,r,cmpl,cmpl2,rr,EV=0,fi;

    f1.i.qFloor = countLinesFromFile(f1,0,&f1.i.iRank);

    c1->i.twoBody.func.fn = nullFunction;
    c1->i.oneBody.func.fn = nullFunction;
    f1.i.nStates =1 ;
    iModel(c1,&f1);
    
#ifndef APPLE
    for ( fi =0 ; fi < f1.i.files ; fi++)
        EV +=  tLoadEigenWeights (c1, f1,f1.i.fileList[fi], f1.f.user+EV);//UNUSUAL!!!
#endif
    if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
    }
    
    
        for ( cmpl = 0; cmpl < spins(f1.f, f1.f.user) ; cmpl++){
            tClear(f1.f,totalVector);
            if ( c1->i.filter  && c1->i.irrep )
            {
                for( g = 0; g < f1.i.qFloor ; g++)
                    tBuildIrr(0, f1.f, c1->i.irrep, f1.f.user+g, cmpl, totalVector, 0);
            }else {
            for( g = 0; g < f1.i.qFloor ; g++)
                tAddTw(f1.f, totalVector, 0, f1.f.user+g, cmpl);
            }
            tCycleDecompostionGridOneMP(-2, f1.f, totalVector, 0, NULL,eigenVectors , cmpl, c1->rt.vCANON, part(f1.f,eigenVectors), c1->rt.powDecompose);
       //     tCycleDecompostionListOneMP(-1, f1, totalVector, 0, NULL,eigenVectors , cmpl, c1->rt.vCANON, part(f1,eigenVectors), 1.);
        }
        tFilter(f1.f, f1.i.nStates, 0, eigenVectors);//classify
    for ( i = 0; i < f1.i.nStates ; i++){
        tScale(f1.f, eigenVectors+i, 1./magnitude(f1.f, eigenVectors+i));
       // printf("printf %f\n", magnitude(f1, eigenVectors+i));
    }
    print(c1,f1,1,0,f1.i.nStates, eigenVectors);
    fModel(&f1.f);

    return 0;
}

#if 1

int main (INT_TYPE argc , char * argv[]){
    struct calculation c;
    struct field f;
    bootShell(argc, argv,&c,&f);
   // test2();
    
    INT_TYPE space,i,a,plusSize,nStatesTrans=0,nStatesFound=0 ,RdsSize = 0,totalIter = 0;
    FILE * out = stdout;
    struct runTime * rt = & c.rt;
    f.f.rt = rt;
    
    
#ifdef OMP
    if ( c.i.omp > MaxCore ){
        printf("lanes > MaxCore\n");
        exit(0);
    }
    rt->NLanes = c.i.omp;
#pragma omp parallel for private (i)
    for ( i = 0; i < MaxCore ; i++){
        rt->NSlot = omp_get_num_threads();
    }
    if ( rt->NLanes > rt->NSlot ){
        printf("lanes > available\n");
        exit(0);
    }
    
#ifdef MKL
    if ( rt->NSlot < c.i.mkl * c.i.omp )
    {
        printf("not enought slots for mkl*omp\n" );
        exit(0);
    }
    else
    {
        printf("mkl %d\n", c.i.mkl);
    }
#endif
    printf("lanes %lld\n", rt->NLanes);
    printf("slots %lld\n", rt->NSlot);
    fflush(stdout);
#endif
    
    //0//...   //A//B//C//D//E
    if ( c.rt.phaseType == buildFoundation ){//0
        foundation(&c,f);
    }
    else if ( c.rt.phaseType == productKrylov ){//C
        krylov(&c,f);
    }
    else if ( c.rt.phaseType == solveRitz ){//A
        ritz(&c,f);
    }
    else if ( c.rt.phaseType == decomposeTensor ){//B
        decompose(&c,f);
    }else if ( c.rt.phaseType == frameDensity ){//D
        if ( ! f.f.vectorOperator){
            printf("vectorOperator flag\n");
            exit(1);
        }
    //    c.i.OCSBflag = 0;

        INT_TYPE EV = foundation(&c,f);
        printf("finished foundation\n");
        fflush(stdout);
        tEigenCycle(f.f,Ha,CDT, f.i.nStates, f.f.user,EV,0, EV,0,0,eigenVectors,twoBodyRitz);
        
        INT_TYPE i,j=0;
        for ( i = 0 ; i < f.i.nStates ; i++ ){
            if ( -myStreams(f.f, twoBodyRitz, 0)[i] > c.rt.TARGET ){
                printf("load %f\n", -myStreams(f.f, twoBodyRitz, 0)[i]);
                j++;
            }
        }
        tEigenCycle(f.f,Ha,CDT, j, f.f.user,EV,0, EV,0,4,eigenVectors,twoBodyRitz);
        
        {
//            INT_TYPE i ;
//            double sum = 0;
//            for ( i = 0; i < j ; i++)
//                sum += -myStreams(f, twoBodyRitz, 0)[i];
            
            cblas_dscal(j, -1., myStreams(f.f, twoBodyRitz, 0), 1);
        }
        
        if ( 0 ){
            INT_TYPE ii;
            j = 1;
            myStreams(f.f, twoBodyRitz, 0)[0] = 1.;
            for ( space = 0; space < SPACE ; space++)
                if ( f.f.rose[space].body != nada   ){
                    for ( i = 0 ; i < vectorLen(f.f, space); i++)
                                if (  i < vectorLen(f.f, space)/2)
                                    streams(f.f, eigenVectors, 0, space)[i] = 1./sqrt(vectorLen(f.f, space)/2);
                                    else
                                        streams(f.f, eigenVectors, 0, space)[i] = 0.;
                    
                }
            f.f.tulip[eigenVectors].Current[0] = 1;
            
            f.f.tulip[diagonal2VectorA].Current[0] = 1;
            for ( space = 0; space < SPACE ; space++)
                if ( f.f.rose[space].body != nada   ){
                    for ( i = 0 ; i < vectorLen(f.f, space); i++)
                        for ( ii = 0 ; ii < vectorLen(f.f, space); ii++)
                            if ( ii < vectorLen(f.f, space)/2 && i < vectorLen(f.f, space)/2)
                                streams(f.f, diagonal2VectorA, 0, space)[ii*vectorLen(f.f,space)+i] = 1./(vectorLen(f.f, space)/2);
                            else
                                streams(f.f, diagonal2VectorA, 0, space)[ii*vectorLen(f.f,space)+i] = 0.;
                    
                }

            
            
            printf("mag %f\n", magnitude(f.f, diagonal2VectorA));
            FILE *out = fopen("uniform.1.0_mac","w");
            outputFormat(f.f, out, diagonal2VectorA, 0);
        }
        
        
        mySeparateEwaldCoulomb1(f.f,j,myStreams(f.f, twoBodyRitz, 0),eigenVectors, c.i.decomposeRankMatrix, interactionExchange,interactionEwald, shortenEwald, 1., 0, 0, electron);
        ioStoreMatrix(f.f, shortenEwald, 0, "shortenEwald.matrix", 0);
        
    }
//    else if ( c.rt.phaseType == collectKrylov ){//F
//        INT_TYPE on,EV = 0,i,fi,ev,Ve,Ven;
//        DCOMPLEX one = 1.;
//        char name[MAXSTRING];
//        c.i.qFloor = countLinesFromFile(&c,f,0);
//        c.i.iRank = c.i.bRank;
//        
//        c.i.nStates = c.i.qFloor ;
//        iModel(&c,f1);
//        for ( fi =0 ; fi < f1.files ; fi++){
//            EV =  tLoadEigenWeights (&c, f1, f1.fileList[fi], f1.user+EV);
//            
//        }
//        if (EV == 0 ){
//            printf ("ack!\n");
//            exit(0);
//        }
//        Ven = 0;
//        Ve = 0;
//        for ( ev = 0 ;ev < EV ; ev++){
//            Ve += tSelect(f1, Ve, 0, eigenVectors, f1.user + ev, 1);
//            // Ve =  tSASplit(f1, c.i.irrep, Ve,EV, eigenVectors,f1.user + ev);
//        }
//        printVector(&c,c.i.c.sinc,c.name,c.name,-1,0, &one);
//        for ( ev =0 ; ev < Ve ; ev++){
//            tFilename(c.name, ev+1, c.rt.body, 0, 0, name);
//            FILE * outVector = fopen(name, "w");
//            outputFormat(f1, outVector, eigenVectors+ev, 0);
//            fclose(outVector);
//            printVector(&c,c.i.c.sinc,c.name,c.name,ev,0, &one);
//        }
//    }else
//        if ( c.rt.phaseType == svdOperation ){//G
//        if ( ! f1.f.vectorOperator){
//            printf("vectorOperator flag\n");
//            exit(1);
//        }
//        INT_TYPE iterator,EV,iii;
//      //  c.i.OCSBflag = 0;
//        EV = foundation(&c);
//        RdsSize = EV;
//
//
//        tEigenCycle(f1.f,Ha,CDT, f1.i.nStates, f1.f.user,RdsSize,0, EV,0,4,eigenVectors,twoBodyRitz);
//        for ( iii = 0; iii < EV ; iii++)
//            print(&c,f1,0,iii,iii+1,eigenVectors);
//
//    }
    
}

#endif
