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
    iModel(c1);

    if ( c1->i.OCSBflag ){
       tBoot1Construction(c1, build);
        EV =   tCollect(f1,c1->i.irrep,f1->sinc.user,c1->i.qFloor ,1);
        tFilter(&c1->i.c, c1->i.qFloor, !(!c1->i.filter )* c1->i.irrep, eigenVectors);
    }
    else{
        tBootManyConstruction(c1);
        EV =   tCollect(f1,c1->i.irrep,f1->sinc.user,c1->i.qFloor ,1);
    }
    return 0;
}


INT_TYPE krylov ( struct calculation *c1){
    struct field * f1 = &c1->i.c;
    INT_TYPE EV = 0,i,fi;


    c1->i.qFloor = countLinesFromFile(c1,0);
    c1->i.iRank = c1->i.bRank;

    c1->i.nTargets = 1;
    c1->i.heliumFlag = 1;
    c1->i.nStates =1 ;
    iModel(c1);
    for ( fi =0 ; fi < c1->mem.files ; fi++){
        EV +=  tLoadEigenWeights (c1, c1->mem.fileList[fi], f1->sinc.user+EV);//UNUSUAL!!!
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
            printExpectationValues(f1, Ha, f1->sinc.user+RdsSize-EV+iii);
            tEdges(c1,f1->sinc.user+RdsSize-EV+iii);
            fflush(stdout);
        }
    }

    
    for ( iterator = 1 ; iterator < c1->i.Iterations ; iterator++){
        RdsSize += tGreatDivideIteration(c1->i.shiftFlag, c1->i.realPart,  f1,Iterator, 1,0,f1->sinc.user+RdsSize-EV,EV,2*EV,0)-EV;
        if(1){
            tFilter(f1, EV, !(!c1->i.filter )* c1->i.irrep, f1->sinc.user+RdsSize-EV);//HERE
            printf ("Step \t%d\n", iterator);
            INT_TYPE iii ;
            for ( iii = 0; iii < EV ; iii++){
                printf ( "\n Vector \t%d \n", iii+1);
                printExpectationValues(f1, Ha, f1->sinc.user+RdsSize-EV+iii);
                tEdges(c1,f1->sinc.user+RdsSize-EV+iii);
                fflush(stdout);
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

    struct field * f1 = &c1->i.c;
    INT_TYPE fi,EV = 0,i,j,sp;
    double va;
    INT_TYPE lines = 0,flines = 0;

    c1->i.qFloor = countLinesFromFile(c1,0);

    c1->i.iRank = c1->i.bRank;
    iModel(c1);
    
    if ( CanonicalRank(f1,interactionExchange,0) )
        ioStoreMatrix(f1,interactionExchange ,0,"interactionExchange.matrix",0);
    
    if ( c1->i.decomposeRankMatrix ){
        
        if ( CanonicalRank(f1,shortenPlus,0) )
            ioStoreMatrix(f1,shortenPlus ,0,"interactionExchangePlus.matrix",0);
        
        if ( CanonicalRank(f1,shortenMinus,0) )
            ioStoreMatrix(f1,shortenMinus ,0,"shortenExchangeMinus.matrix",0);

        
    }else {
        if ( CanonicalRank(f1,interactionExchangePlus,0) )
            ioStoreMatrix(f1,interactionExchangePlus ,0,"interactionExchangePlus.matrix",0);
        
        if ( CanonicalRank(f1,interactionExchangeMinus,0) )
            ioStoreMatrix(f1,interactionExchangeMinus ,0,"shortenExchangeMinus.matrix",0);
    }
    if ( PARTICLE == 2) {
        if ( CanonicalRank(f1,protonRepulsion,0) )
            ioStoreMatrix(f1,protonRepulsion ,0,"protonRepulsion.matrix",0);
    }
    
    if ( CanonicalRank(f1,linear,0) )
        ioStoreMatrix(f1,linear ,0,"linear.matrix",0);
    
    if ( CanonicalRank(f1,kinetic,0) )
        ioStoreMatrix(f1,kinetic ,0,"kinetic.matrix",0);

    
    
    
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

    INT_TYPE g,r,cmpl,cmpl2,rr,EV=0,fi;

    c1->i.qFloor = countLinesFromFile(c1,0);

    c1->i.c.twoBody.func.fn = nullFunction;
    c1->i.c.oneBody.func.fn = nullFunction;
    c1->i.iRank = c1->i.bRank;
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
            for( g = 0; g < c1->i.qFloor ; g++)
                tAddTw(f1, totalVector, 0, f1->sinc.user+g, cmpl);
            tCycleDecompostionListOneMP(-1, f1, totalVector, 0, NULL,eigenVectors , cmpl, c1->rt.vCANON, part(f1,eigenVectors), 1.);
        }
    tFilter(f1, c1->i.nStates, !(!c1->i.filter )* c1->i.irrep, eigenVectors);

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
    

    
    printf("\n\n\t  \tBox %d \t: Lattice  %1.3f \t\n\t| Attack :  %1.3f\t\t\t\t\n",2* c.i.epi + 1,c.i.d,c.i.attack);
    
    if ( SPACE > 3 )
        printf("\n\n\t  \tBox %d \t: Lattice  %1.3f \t\n\t| Attack :  %1.3f\t\t\t\t\n",2* c.i.around + 1,c.i.D,   1.0 );    
    
    //0//...   //A//B//C//D//E
    if ( c.rt.phaseType == buildFoundation ){//0
        c.i.nTargets = 0;
        c.i.heliumFlag = c.i.nTargets;
        c.i.nStates = c.i.heliumFlag;

        foundation(&c);
        print(&c,c.i.qFloor, f1->sinc.user);
    }
    else if ( c.rt.phaseType == productKrylov ){//C
        krylov(&c);
        print(&c,c.i.Iterations,f1->sinc.user);
    }
    else if ( c.rt.phaseType == solveRitz ){//A
        ritz(&c);
    }
    else if ( c.rt.phaseType == decomposeTensor ){//B
        decompose(&c);
        print(&c,c.i.nStates, eigenVectors);
    }else if ( c.rt.phaseType == frameDensity ){//D
        if ( ! c.i.vectorOperatorFlag){
            printf("vectorOperator flag\n");
            exit(1);
        }
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

        mySeparateEwaldCoulomb1(&c.i.c,j,eigenVectors, c.i.decomposeRankMatrix, interactionEwald, shortenEwald, 1./c.rt.body, 0, 1, electron);
        if ( c.i.decomposeRankMatrix )
            ioStoreMatrix(&c.i.c, shortenEwald, 0, "EwaldCoulomb.matrix", 0);
        else
            ioStoreMatrix(&c.i.c, interactionEwald, 0, "EwaldCoulomb.matrix", 0);

    }



    fModel(&c);

}
