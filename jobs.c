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
        c1->i.twoBody.func.fn = nullFunction;
        iModel(c1,&f1);
        tBoot1Construction(c1,f1.f ,build);
        tSortBoot(c1,f1.f,build);
        EV =   tSlam(f1.f,f1.i.qFloor,f1.f.user,c1->i.level);

//        tGreatDivideIteration(0, 0, f1.f, Ha, 1, 0, f1.f.user, 1, 2, 0);
//        tFilter(f1.f, EV,1, f1.f.user);//classify
        if ( OVERFLAG  || SPACE == 1 )
            tEigenCycle(0,f1.f,Ha,CDT, f1.i.nStates, f1.f.user,EV,0, EV,0,1,eigenVectors,twoBodyRitz);

    }else {
        exit(1);
    }
    
    INT_TYPE ii,flag = 1;;
        print(c1,f1,1,0,EV , f1.f.user);
    
    fModel(&f1.f);

    return EV;
}


INT_TYPE krylov ( struct calculation *c1, struct field f1){
    INT_TYPE EV = 0,i,fi;


    f1.i.qFloor = countLinesFromFile(c1,f1,0,&f1.i.iRank);
    //count canonical-rank...

    
    f1.i.nStates =1  ;
    f1.i.nStates =f1.i.Iterations  ;
   
 iModel(c1,&f1);
    for ( fi =0 ; fi < f1.i.files ; fi++){
        tLoadEigenWeights (c1,f1, f1.i.fileList[fi], &EV,f1.f.user , f1.i.collect);//UNUSUAL!!!
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
                for( g = 0; g < EV ; g++)
                    tBuild3Irr(0, f1.f, f1.i.irrep, f1.f.user+g, cmpl, totalVector, 0);
            }else {
                for( g = 0; g < EV ; g++)
                    tAddTw(f1.f, totalVector, 0, f1.f.user+g, cmpl);
            }
            double norm = magnitude(f1.f, totalVector );
            if ( norm > c1->rt.TARGET )
                printf("Normed from %f\n", norm );
            else
            {
                print(c1,f1,1,0,0,eigenVectors);
                return 1;
            }
            tScaleOne(f1.f, totalVector, 0, 1/norm);
            tCycleDecompostionGridOneMP(-2, f1.f, totalVector, 0, NULL,eigenVectors , cmpl, c1->rt.vCANON, part(f1.f,eigenVectors), c1->rt.powDecompose);
            //     tCycleDecompostionListOneMP(-1, f1, totalVector, 0, NULL,eigenVectors , cmpl, c1->rt.vCANON, part(f1,eigenVectors), 1.);
        }
        
        EV = 1;
        RdsSize = 1;
            tFilter(f1.f, EV, 0, eigenVectors);//classify
            printExpectationValues(f1.f, Ha, eigenVectors);
            fflush(stdout);
            print(c1,f1,1,0,1,eigenVectors);
        }
    

    
    for ( iterator = 1 ; iterator < f1.i.Iterations ; iterator++){
        RdsSize += tGreatDivideIteration(c1->i.shiftFlag, c1->i.realPart,  f1.f,Iterator, 1,0,eigenVectors+RdsSize-EV,EV,2*EV,0)-EV;
        if(1){
            tFilter(f1.f, EV, !(!f1.i.filter )* f1.i.irrep, eigenVectors+RdsSize-EV);//filter
            printf ("Step \t%d\n", iterator);
            fflush(stdout);
            INT_TYPE iii ;
            for ( iii = 0; iii < EV ; iii++){
                printf ( "\n Vector \t%d \t %d\n", iii+1, +RdsSize-EV+iii);
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
    //c1->i.canonRank = 0;
    INT_TYPE fi,EV = 0,i,j,sp;
    double va;
    INT_TYPE lines = 0,flines = 0;

    f1.i.qFloor = countLinesFromFile(c1,f1,0,&f1.i.iRank);
  //  printf("-->%d\n", f1.i.iRank);
   // c1->i.iRank = c1->i.bRank;
    iModel(c1,&f1);
    
    
    if ( CanonicalRank(f1.f,interactionExchange,0) ){
        ioStoreMatrix(f1.f,interactionExchange ,0,"interactionExchange.matrix",0);
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

INT_TYPE svd ( struct calculation c, struct field f1){//E
    if ( ! f1.i.filesVectorOperator){
        printf("vectorOperator flag\n");
        exit(1);
    }
    
    INT_TYPE EV = foundation(&c,f1);
    printf("finished foundation\n");
    fflush(stdout);
    tEigenCycle(0,f1.f,Ha,'T', f1.i.nStates, f1.f.user,EV,0, EV,0,0,eigenVectors,twoBodyRitz);
    
    INT_TYPE i,j=0;
    for ( i = 0 ; i < f1.i.nStates ; i++ ){
        if ( -myStreams(f1.f, twoBodyRitz, 0)[i] > c.rt.TARGET ){
            printf("load %f\n", -myStreams(f1.f, twoBodyRitz, 0)[i]);
            j++;
        }
    }
//    tEigenCycle(0,f1.f,Ha,'T', j, f1.f.user,EV,0, EV,0,4,eigenVectors,twoBodyRitz);
//    
//    {
//           cblas_dscal(j, -1., myStreams(f1.f, twoBodyRitz, 0), 1);
//    }
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
        printf("decrease lanes or increase your number of OMP cores\n");
        exit(0);
    }
    
#ifdef MKL
    if ( rt->NSlot < c.i.mkl * c.i.omp )
    {
        printf("not enough slots for mkl*omp\n" );
        exit(0);
    }
    else
    {
        rt->NParallel = c.i.mkl;
    }
#endif
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
    else if ( c.rt.phaseType == svdOperation ){//E
        ritz(&c,f);
    }

    printf("\n\nFINIS.\n\n");
}

#endif
