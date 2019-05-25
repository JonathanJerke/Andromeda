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

INT_TYPE foundation(struct calculation *c1){
    struct field * f1 = &c1->i.c;

    INT_TYPE EV;

    if ( !c1->i.OCSBflag ){
        c1->i.c.twoBody.func.fn = nullFunction;
        iModel(c1);
        tBoot1Construction(c1, build);
        EV =   tSlam(f1,c1->i.qFloor,f1->sinc.user,c1->i.level);
        tFilter(&c1->i.c, EV,0, f1->sinc.user);//classify
    }
    else{
        iModel(c1);

        tBootManyConstruction(c1);
        EV =   tCollect(f1,c1->i.irrep,f1->sinc.user,c1->i.qFloor ,1);
        tFilter(&c1->i.c, EV, 0, f1->sinc.user);//classify
    }
    return EV;
}


INT_TYPE krylov ( struct calculation *c1){
    struct field * f1 = &c1->i.c;
    INT_TYPE EV = 0,i,fi;


    c1->i.qFloor = countLinesFromFile(c1,0);
    i =c1->i.iRank;
    c1->i.iRank = c1->i.bRank;
    c1->i.bRank = i;
    
    c1->i.nTargets = c1->i.qFloor;
    c1->i.heliumFlag = c1->i.qFloor;
    c1->i.nStates =c1->i.qFloor  ;
    iModel(c1);
    for ( fi =0 ; fi < c1->mem.files ; fi++){
        EV +=  tLoadEigenWeights (c1, c1->mem.fileList[fi], eigenVectors);//UNUSUAL!!!
    }
        if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
    }
    INT_TYPE RdsSize = EV,iterator=0;

    if(1){
        printf ("Step \t%d\n", iterator);
        INT_TYPE iii ;
        for ( iii = 0; iii < EV ; iii++){
            printf ( "\n Vector \t%d \n", iii+1);
            tFilter(f1, EV, !(!c1->i.filter )* c1->i.irrep, eigenVectors+RdsSize-EV);//classify
            printExpectationValues(f1, Ha, eigenVectors+RdsSize-EV+iii);
            fflush(stdout);
            print(c1,1,RdsSize-EV+iii,RdsSize-EV+iii+1,eigenVectors);
        }
    }

    
    for ( iterator = 1 ; iterator < c1->i.Iterations ; iterator++){
        RdsSize += tGreatDivideIteration(c1->i.shiftFlag, c1->i.realPart,  f1,Iterator, 1,0,eigenVectors+RdsSize-EV,EV,2*EV,0)-EV;
        if(1){
            tFilter(f1, EV, !(!c1->i.filter )* c1->i.irrep, eigenVectors+RdsSize-EV);//classify
            printf ("Step \t%d\n", iterator);
            INT_TYPE iii ;
            for ( iii = 0; iii < EV ; iii++){
                printf ( "\n Vector \t%d \n", iii+1);
                printExpectationValues(f1, Ha, eigenVectors+RdsSize-EV+iii);
                fflush(stdout);
                print(c1,0,RdsSize-EV+iii,RdsSize-EV+iii+1,eigenVectors);
            }
        }
        fflush(stdout);
    }
    return 0;
}

INT_TYPE ritz( struct calculation * c1 ){
    size_t ms = MAXSTRING;
    char line0[MAXSTRING];
    char filename[MAXSTRING];    char str[MAXSTRING];


    char * line = line0;
    c1->i.canonRank = 0;
    struct field * f1 = &c1->i.c;
    INT_TYPE fi,EV = 0,i,j,sp;
    double va;
    INT_TYPE lines = 0,flines = 0;

    c1->i.qFloor = countLinesFromFile(c1,0);

   // c1->i.iRank = c1->i.bRank;
    iModel(c1);
    
    if ( CanonicalRank(f1,interactionExchange,0) )
        ioStoreMatrix(f1,interactionExchange ,0,"interactionExchange.matrix",0);
    if ( CanonicalRank(f1,interactionEwald,0) )
        ioStoreMatrix(f1,interactionEwald ,0,"interactionEwald.matrix",0);

    if ( c1->i.decomposeRankMatrix ){
        
        if ( CanonicalRank(f1,shortenPlus,0) )
            ioStoreMatrix(f1,shortenPlus ,0,"shortenExchangePlus.matrix",0);
        
        if ( CanonicalRank(f1,shortenMinus,0) )
            ioStoreMatrix(f1,shortenMinus ,0,"shortenExchangeMinus.matrix",0);

        
    }else {
        if ( CanonicalRank(f1,interactionExchangePlus,0) )
            ioStoreMatrix(f1,interactionExchangePlus ,0,"interactionExchangePlus.matrix",0);
        
        if ( CanonicalRank(f1,interactionExchangeMinus,0) )
            ioStoreMatrix(f1,interactionExchangeMinus ,0,"interactionExchangeMinus.matrix",0);
    }
    
    if ( CanonicalRank(f1,linear,0) )
        ioStoreMatrix(f1,linear ,0,"linear.matrix",0);

    
    
    
#ifndef APPLE
    for ( fi =0 ; fi < c1->mem.files ; fi++)
        EV +=  tLoadEigenWeights (c1, c1->mem.fileList[fi], f1->sinc.user+EV);//UNUSUAL!!!
    if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
    }
#endif
    tEigenCycle(&c1->i.c,Ha,'T', c1->i.nStates, f1->sinc.user,EV,0, EV,0,1,eigenVectors,twoBodyRitz);
    
    
    DCOMPLEX *V = (DCOMPLEX*)myStreams(&c1->i.c,matrixHbuild,0);
    INT_TYPE iii,stride = c1->i.c.sinc.maxEV;
    char * token = filename;
    for ( iii = 0; iii < c1->i.nStates ; iii++){
        
        {
            FILE * outf ;
            sprintf(str, "%s-%d.vector",c1->name,iii+1);
            outf = fopen (str,"w");
            fclose(outf);
        }
        lines = 0;

        for ( fi =0 ; fi < c1->mem.files ; fi++){
            flines = 0;
            FILE * fp = fopen(c1->mem.fileList[fi],"r");
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
                    
                    printVector(c1,token,str,flines,c1->i.irrep, V+stride*iii+lines);
                    //   printf("v%d s%d s%d\n", iii, stride,lines);
                    lines++;
                    flines++;
                }
                getline(&line, &ms, fp);

            }
            
            fclose(fp);
        }
    }
    
    return 0;
}

INT_TYPE decompose( struct calculation * c1 ){
    struct field * f1 = &c1->i.c;

    INT_TYPE i,g,r,cmpl,cmpl2,rr,EV=0,fi;

    c1->i.qFloor = countLinesFromFile(c1,0);

    c1->i.c.twoBody.func.fn = nullFunction;
    c1->i.c.oneBody.func.fn = nullFunction;
    c1->i.nTargets = 1;
    c1->i.heliumFlag = 1;
    c1->i.nStates =1 ;
    iModel(c1);
    
#ifndef APPLE
    for ( fi =0 ; fi < c1->mem.files ; fi++)
        EV +=  tLoadEigenWeights (c1, c1->mem.fileList[fi], f1->sinc.user+EV);//UNUSUAL!!!
#endif
    if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
    }
    
    
        for ( cmpl = 0; cmpl < spins(f1, f1->sinc.user) ; cmpl++){
            tClear(f1,totalVector);
            if ( c1->i.filter  && c1->i.irrep )
            {
                for( g = 0; g < c1->i.qFloor ; g++)
                    tBuildIrr(0, f1, c1->i.irrep, f1->sinc.user+g, cmpl, totalVector, 0);
            }else {
            for( g = 0; g < c1->i.qFloor ; g++)
                tAddTw(f1, totalVector, 0, f1->sinc.user+g, cmpl);
            }
            tCycleDecompostionGridOneMP(-2, f1, totalVector, 0, NULL,eigenVectors , cmpl, c1->rt.vCANON, part(f1,eigenVectors), c1->rt.powDecompose);
       //     tCycleDecompostionListOneMP(-1, f1, totalVector, 0, NULL,eigenVectors , cmpl, c1->rt.vCANON, part(f1,eigenVectors), 1.);
        }
        tFilter(f1, c1->i.nStates, !(!c1->i.filter )* c1->i.irrep, eigenVectors);//classify
    for ( i = 0; i < c1->i.nStates ; i++){
        tScale(f1, eigenVectors+i, 1./magnitude(f1, eigenVectors+i));
       // printf("printf %f\n", magnitude(f1, eigenVectors+i));
    }
    return 0;
}



int main (INT_TYPE argc , char * argv[]){
    struct calculation c = bootShell(argc, argv);
    
    INT_TYPE i,a,plusSize,nStatesTrans=0,nStatesFound=0 ,RdsSize = 0,totalIter = 0;
    FILE * out = stdout;
    struct runTime * rt = & c.rt;
    
    struct field *f1 = &(c.i.c);
    if ( c.i.iCharge > 0 )
        f1->Ne = c.i.iCharge;
    
    
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
    
    
    
    printf("\n\n\t  \tBox %d \t: Lattice  %1.3f \t\n\t\n",2* c.i.epi + 1,c.i.d);
    
    if ( SPACE > 3 )
        printf("\n\n\t  \tBox %d \t: Lattice  %1.3f \t\n\t\n",2* c.i.around + 1,c.i.D );
    
    //0//...   //A//B//C//D//E
    if ( c.rt.phaseType == buildFoundation ){//0
        c.i.nTargets = 0;
        c.i.heliumFlag = c.i.nTargets;
        c.i.nStates = c.i.heliumFlag;
        INT_TYPE EV= foundation(&c);
        print(&c,1,0,EV , f1->sinc.user);
    }
    else if ( c.rt.phaseType == productKrylov ){//C
        krylov(&c);
    }
    else if ( c.rt.phaseType == solveRitz ){//A
        ritz(&c);
    }
    else if ( c.rt.phaseType == decomposeTensor ){//B
        decompose(&c);
        print(&c,1,0,c.i.nStates, eigenVectors);
    }else if ( c.rt.phaseType == frameDensity ){//D
        if ( ! c.i.vectorOperatorFlag){
            printf("vectorOperator flag\n");
            exit(1);
        }
        c.i.OCSBflag = 0;

        foundation(&c);
        tEigenCycle(&c.i.c,Ha,'T', c.i.nStates, c.i.c.sinc.user,c.i.qFloor,0, c.i.qFloor,0,0,eigenVectors,twoBodyRitz);
        
        INT_TYPE i,j=0;
        for ( i = 0 ; i < c.i.nStates ; i++ ){
            if ( -myStreams(f1, twoBodyRitz, 0)[i] > c.rt.TARGET ){
                printf("load %f\n", -myStreams(f1, twoBodyRitz, 0)[i]);
                j++;
            }
        }
        tEigenCycle(&c.i.c,Ha,'T', j, c.i.c.sinc.user,c.i.qFloor,0, c.i.qFloor,0,4,eigenVectors,twoBodyRitz);
        
        {
//            INT_TYPE i ;
//            double sum = 0;
//            for ( i = 0; i < j ; i++)
//                sum += -myStreams(f1, twoBodyRitz, 0)[i];
            
            cblas_dscal(j, -1., myStreams(f1, twoBodyRitz, 0), 1);
        }
        
        
        mySeparateEwaldCoulomb1(&c.i.c,j,myStreams(f1, twoBodyRitz, 0),eigenVectors, c.i.decomposeRankMatrix, interactionExchange,interactionEwald, shortenEwald, 0, 0, 0, electron);
        ioStoreMatrix(&c.i.c, shortenEwald, 0, "shortenEwald.matrix", 0);
        
    }
    else if ( c.rt.phaseType == scaleBody ){//E
        if ( ! c.i.bodyFlag){
            printf("body flag\n");
            exit(1);
        }
        INT_TYPE fi,EV=0;
        c.i.qFloor = countLinesFromFile(&c,0);
        c.i.c.twoBody.func.fn = nullFunction;
        c.i.c.oneBody.func.fn = nullFunction;
        c.i.iRank = c.i.bRank;
        c.i.nTargets = 1;
        c.i.heliumFlag = 1;
        c.i.nStates =1 ;
        iModel(&c);
        
#ifndef APPLE
        for ( fi =0 ; fi < c.mem.files ; fi++)
            EV +=  tLoadEigenWeights (&c, c.mem.fileList[fi], c.i.c.sinc.user+EV);//UNUSUAL!!!
#endif
        if (EV == 0 ){
            printf ("ack!\n");
            exit(0);
        }
        
        tMap(&c);
    }
    else if ( c.rt.phaseType == collectKrylov ){//F
        struct field * f1 = &c.i.c;
        INT_TYPE on,EV = 0,i,fi,ev,Ve,Ven;
        DCOMPLEX one = 1.;
        char name[MAXSTRING];
        c.i.qFloor = countLinesFromFile(&c,0);
        c.i.iRank = c.i.bRank;
        
        c.i.nTargets = c.i.qFloor ;
        c.i.heliumFlag = c.i.nTargets;
        c.i.nStates = c.i.nTargets ;
        iModel(&c);
        for ( fi =0 ; fi < c.mem.files ; fi++){
            EV =  tLoadEigenWeights (&c, c.mem.fileList[fi], c.i.c.sinc.user+EV);
            
        }
        if (EV == 0 ){
            printf ("ack!\n");
            exit(0);
        }
        Ven = 0;
        Ve = 0;
        for ( ev = 0 ;ev < EV ; ev++){
            Ve += tSelect(f1, Ve, 0, eigenVectors, f1->sinc.user + ev, 1);
            // Ve =  tSASplit(f1, c.i.irrep, Ve,EV, eigenVectors,f1->sinc.user + ev);
        }
        printVector(&c,c.name,c.name,-1,0, &one);
        for ( ev =0 ; ev < Ve ; ev++){
            tFilename(c.name, ev+1, c.i.c.body, 0, 0, name);
            FILE * outVector = fopen(name, "w");
            outputFormat(f1, outVector, eigenVectors+ev, 0);
            fclose(outVector);
            printVector(&c,c.name,c.name,ev,0, &one);
        }
    }else if ( c.rt.phaseType == svdOperation ){//G
        if ( ! c.i.vectorOperatorFlag){
            printf("vectorOperator flag\n");
            exit(1);
        }
        INT_TYPE iterator,EV;
        struct calculation *c1 = &c;
        c.i.OCSBflag = 0;
        EV = foundation(&c);
        RdsSize = EV;
        tEigenCycle(&c.i.c,Ha,'T', c.i.nStates, c.i.c.sinc.user,RdsSize,0, EV,0,0,eigenVectors,twoBodyRitz);
        
    }
    
    fModel(&c);
}
