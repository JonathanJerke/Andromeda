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
        f1.i.irrep = 0;
        c1->i.twoBody.func.fn = nullFunction;
        iModel(c1,&f1);
        tBoot1Construction(c1,f1.f ,build);

        EV =   tSlam(f1.f,f1.i.qFloor,f1.f.user,c1->i.level);

//        tGreatDivideIteration(0, 0, f1.f, Ha, 1, 0, f1.f.user, 1, 2, 0);
//        tFilter(f1.f, EV,1, f1.f.user);//classify
        if ( OVERFLAG  || SPACE == 1 )
            tEigenCycle(f1.f,Ha,CDT, f1.i.nStates, f1.f.user,EV,0, EV,0,1,eigenVectors,twoBodyRitz);

    }else {
        exit(1);
    }
    
    INT_TYPE ii,flag = 1;;
    if ( f1.i.irrep )
        for ( ii = 0 ;ii < EV ;ii++){
            if ( f1.f.tulip[f1.f.user+ii].value.symmetry == f1.i.irrep ){
                print(c1,f1,flag,ii,ii+1 , f1.f.user);
                flag = 0;
            }
        }
    else
        print(c1,f1,1,0,EV , f1.f.user);
    
    fModel(&f1.f);

    return EV;
}


INT_TYPE krylov ( struct calculation *c1, struct field f1){
    INT_TYPE EV = 0,i,fi;


    f1.i.qFloor = countLinesFromFile(f1,0,&f1.i.iRank);
    //count canonical-rank...

    
    f1.i.nStates =1  ;
    f1.i.nStates =f1.i.Iterations  ;
   
 iModel(c1,&f1);
    for ( fi =0 ; fi < f1.i.files ; fi++){
        tLoadEigenWeights (c1,f1, f1.i.fileList[fi], &EV,f1.f.user , f1.i.collect);//UNUSUAL!!!
    }
        if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
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
            if ( f1.i.filter  && f1.i.irrep )
            {
                for( g = 0; g < EV ; g++)
                    tBuildIrr(0, f1.f, f1.i.irrep, f1.f.user+g, cmpl, totalVector, 0);
            }else {
                for( g = 0; g < EV ; g++)
                    tAddTw(f1.f, totalVector, 0, f1.f.user+g, cmpl);
            }
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

    
    
    
    for ( fi =0 ; fi < f1.i.files ; fi++)
        tLoadEigenWeights (c1,f1, f1.i.fileList[fi],&EV, f1.f.user,f1.i.collect);//UNUSUAL!!!
    if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
    }
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
    
    for ( fi =0 ; fi < f1.i.files ; fi++)
        tLoadEigenWeights (c1, f1,f1.i.fileList[fi], &EV,f1.f.user,f1.i.collect);//UNUSUAL!!!
    if (EV == 0 ){
        printf ("ack!\n");
        exit(0);
    }
    
    
        for ( cmpl = 0; cmpl < spins(f1.f, f1.f.user) ; cmpl++){
            tClear(f1.f,totalVector);
            if ( f1.i.filter  && f1.i.irrep )
            {
                for( g = 0; g < f1.i.qFloor ; g++)
                    tBuildIrr(0, f1.f, f1.i.irrep, f1.f.user+g, cmpl, totalVector, 0);
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

INT_TYPE frameEwald( struct calculation * c , struct field f)
{//D
    INT_TYPE space,EV ;

    INT_TYPE oi =1;
    f.i.nOperator = countLinesFromFile(f,1,&oi);

    if ( ! f.i.nOperator){
        printf("vectorOperator flag\n");
        exit(1);
    }
    
    {
        f.i.irrep = 0;
        iModel(c,&f);
        tBoot1Construction(c,f.f ,build);
        
        EV =   tSlam(f.f,f.i.qFloor,f.f.user,c->i.level);
    }
    
    
    printf("finished foundation\n");
    fflush(stdout);
    tEigenCycle(f.f,Ha,CDT, f.i.nStates, f.f.user,EV,0, EV,0,0,eigenVectors,twoBodyRitz);
    
    

    INT_TYPE i,j=0;
    double totalElectron = 0., consideredElectron = 0.;
    for ( i = 0 ; i < EV ; i++ ){
        if ( -myStreams(f.f, twoBodyRitz, 0)[i] > c->rt.EWALD ){
            printf("load-%d %f\n",i+1, myStreams(f.f, twoBodyRitz, 0)[i]);
            consideredElectron -= myStreams(f.f, twoBodyRitz, 0)[i];
            j++;//in order
        }
        totalElectron -= myStreams(f.f, twoBodyRitz, 0)[i];
    }
    printf ( " encapsulating %f probability\n considering %f probability using target %f\n",totalElectron, consideredElectron, c->rt.TARGET );
    
    
    tEigenCycle(f.f,Ha,CDT, j, f.f.user,EV,0, EV,0,4,eigenVectors,twoBodyRitz);
    print(c,f,1,0,j,f.f.user);

    if (1 ) {
        printf("normalizing considered probability to 1\n");
        cblas_dscal(EV, 1./consideredElectron, myStreams(f.f, twoBodyRitz, 0), 1);
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
    
    
    mySeparateEwaldCoulomb1(f.f,j,myStreams(f.f, twoBodyRitz, 0),eigenVectors, c->i.decomposeRankMatrix, interactionExchange, intracellularSelfEwald,-1., 0, 0, electron);
    ioStoreMatrix(f.f, intracellularSelfEwald, 0, "intracellularEwald.matrix", 0);
    
    mySeparateEwaldCoulomb1(f.f,j,myStreams(f.f, twoBodyRitz, 0),eigenVectors, c->i.decomposeRankMatrix,interactionEwald, intercellularSelfEwald, 1., 0, 0, electron);
    ioStoreMatrix(f.f, intercellularSelfEwald, 0, "intercellularEwald.matrix", 0);

    fModel(&f.f);

    return 0;
}


INT_TYPE collectSet( struct calculation * c , struct field f)
{
            INT_TYPE EV = 0,fi,ev,Ve,Ven;
            DCOMPLEX one = 1.;
            char name[MAXSTRING];
            f.i.qFloor = countLinesFromFile(f,0,&f.i.iRank);
            f.i.iRank = f.i.bRank;
    
            f.i.nStates = f.i.qFloor ;
            iModel(c,&f);
            for ( fi =0 ; fi < f.i.files ; fi++){
                tLoadEigenWeights (c, f, f.i.fileList[fi], &EV,f.f.user, f.i.collect);
    
            }
            if (EV == 0 ){
                printf ("ack!\n");
                exit(0);
            }
            Ven = 0;
            Ve = 0;
            for ( ev = 0 ;ev < EV ; ev++){
                Ve += tSelect(f.f, Ve, 0, eigenVectors, f.f.user + ev, 1);
            }
            printVector(c,f.f,c->name,c->name,-1,0, &one);
            for ( ev =0 ; ev < Ve ; ev++){
                tFilename(c->name, ev+1, f.i.body, 0, 0, name);
                FILE * outVector = fopen(name, "w");
                outputFormat(f.f, outVector, eigenVectors+ev, 0);
                fclose(outVector);
                printVector(c,f.f,c->name,c->name,ev,0, &one);
            }
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
    else if ( c.rt.phaseType == decomposeTensor ){//B
        decompose(&c,f);
    }else if ( c.rt.phaseType == frameDensity ){//D
        frameEwald (&c,f);
    }else if ( c.rt.phaseType == collectKrylov ){//F
        collectSet(&c,f);
    }
}

#endif
