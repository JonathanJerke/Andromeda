/**
 *  jobs.c
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

#include "jobs.h"

/**
 * Built a spherical foundation of one vector
 */
inta foundationS(  calculation *c1,   field f1){
    inta EV;
    f1.i.Iterations = 2;
    f1.i.qFloor = 2;
    f1.i.iRank = 3;
    f1.i.canonRank = 4;
    EV = f1.i.qFloor ;
    if ( 1 ){
        iModel(c1,&f1);
        
        tBoot(f1.f, f1.f.user, 0);
        tId(f1.f, f1.f.user, 0);

        printf("irrep%d\n",tClassify(f1.f, f1.f.user));
        tHXpY( f1.f, f1.f.user+1, f1.f.name[defSpiralMatrix(&f1.f,Iterator)].name, 0, f1.f.user, 1e-4, 1e-4, 1e-8, 1e-15, 5,1000, 5, 1);
        printExpectationValues(c1, f1.f, Iterator, totalVector);
        printf("irrep%d\n",tClassify(f1.f, totalVector));

        print(c1,f1,1,0,f1.i.qFloor , f1.f.user);
        fModel(&f1.f);
    }
    return EV;
}

/**
 *Build a complete digital foundation
 *
 *One vector per vector element
 */
inta foundationB(  calculation *c1,   field f1){
    inta EV;
    f1.i.Iterations = 1;
    inta space,i;
    f1.i.qFloor = 1;
    inta ssp;
    for ( space = 0 ;space< SPACE ; space++)
        if ( f1.f.canon[space].body != nada)
            f1.i.qFloor *= vectorLen(f1.f, space);
    
    EV = f1.i.qFloor ;
    if ( 1 ){
        iModel(c1,&f1);
        for ( i = 0; i < f1.i.qFloor ; i++)
        {
            ssp = 1;
            zero(f1.f, f1.f.user+i, 0);
            for  ( space =0; space < SPACE ; space++){
                if ( f1.f.canon[space].body != nada){
                    streams(f1.f,f1.f.user+i,0,space)[(i/ssp)%vectorLen(f1.f, space)] = 1;
                    ssp *=vectorLen(f1.f, space);
                }
            }
            f1.f.name[f1.f.user+i].Current[0] = 1;
        }
        print(c1,f1,1,0,f1.i.qFloor , f1.f.user);
        fModel(&f1.f);
    }
    return EV;
}

#if 0
/**
 *Test Permutation actions
 *
 */
double testPermutations (){
    calculation c2 = initCal();
    field f2 = initField();
    c2.i.numNames = 1000;
    c2.i.numVectors = 0;

    c2.i.termNumber = 0;
    c2.rt.NLanes = 1;
    f2.f.rt = &c2.rt;
    f2.f.rt->phaseType = productKrylov;
    
    f2.i.Iterations = 1;
    f2.i.files = 0;
    f2.i.filesVectorOperator = 0;
    f2.i.qFloor = 0;
    c2.i.lambda = 6;
    resetA(f2.f.rt);
    blockA(f2.f.rt, blockTrainVectorsblock);
    blockA(f2.f.rt, blockCopyBlock);
    blockA(f2.f.rt, blockMatrixElementsblock);
    blockA(f2.f.rt, blockPermutationsblock);
    blockA(f2.f.rt, blockParallelMultiplyblock);
    blockA(f2.f.rt, blockParallelMatrixElementblock);
    blockA(f2.f.rt, blockParallelPermuteblock);
    f2.i.body = two;
    f2.f.boot = noMatrices;
    
    
    return 0;
}
#endif

/**
 *Quad Traces
 *Tr [ mu.Transpose(mu) . ( a . Transpose( b ) ) ]
 *
 *outputs files labeled with mu,
 *each file is a matrix of a,b
 */
double quads ( calculation *c1, field f1){
    if ( ! allowQ(f1.f.rt,blockMatrixElementsblock)){
        printf("blockMatrixElementsblock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    inta EV = 0;
    inta fi,a,b,mu;
    char filename[MAXSTRING];
    {
        f1.i.nStates = countLinesFromFile(c1,f1,0,&f1.i.iRank, &f1.i.xRank);
        f1.i.qFloor = f1.i.nStates*f1.i.nStates;
        printf("load %d states\n",f1.i.nStates);
        iModel(c1,&f1);
        for ( fi = 0 ; fi < f1.i.files ; fi++){
            tLoadEigenWeights (c1,f1, f1.i.fileList[fi], &EV,eigenVectors , 0);
        }
        printf("loaded %d states\n",EV);

        for ( a = 0 ; a < f1.i.nStates ; a++){
            for ( b = 0 ; b < f1.i.nStates ; b++){
                tHXpY(f1.f, f1.f.user+f1.i.nStates*a+b, eigenVectors+a, 0, eigenVectors+b, f1.f.rt->TOLERANCE,f1.f.rt->relativeTOLERANCE,f1.f.rt->ALPHA,f1.f.rt->THRESHOLD,f1.f.rt->MAX_CYCLE,f1.f.rt->XCONDITION, part(f1.f,f1.f.user+f1.i.nStates*a+b), part(f1.f,f1.f.user+f1.i.nStates*a+b));
            }
        }
        
        for ( mu = 0 ; mu < f1.i.nStates ; mu++){
            for ( a = 0 ; a < f1.i.nStates ; a++){
                for ( b = 0 ; b < f1.i.nStates ; b++){
                    myStreams(f1.f,matrixSbuild,0)[f1.i.nStates*a+b] = pMatrixElement( f1.f, f1.f.user+f1.f.user*mu+mu, 0, nullOverlap, 0, f1.f.user+f1.f.user*a+b, 0);
                }
            }
            
            tFilename(c1->name, mu, 0, 0, 0, filename);
            ioArray(c1, f1, filename, f1.i.qFloor, (mea*)myStreams(f1.f,matrixSbuild,0), 0);
        }
        fModel(&f1.f);
    }
    return 0.;
}



/**
 *Multiply, decompose, and filter
 *
 *when totalVector is allocated, it will sum up all 'vector' inputs
 *when totalVector is not allocated, it will take 'solo' input
 *totalVector is blocked with 'blockMemory 1'
 */
double singlekrylov (   calculation *c1,   field f1){
    inta EV = 0;
#ifndef APPLE
    inta cmpl,g,fi;
    {
        f1.i.qFloor = countLinesFromFile(c1,f1,0,&f1.i.iRank, &f1.i.xRank);
        f1.i.nStates = 1;
        iModel(c1,&f1);
        for ( fi = 0 ; fi < f1.i.files ; fi++){
            tLoadEigenWeights (c1,f1, f1.i.fileList[fi], &EV,f1.f.user , f1.i.collect);
        }
        for ( cmpl = 0 ; cmpl < 1 ; cmpl++){
            tClear(f1.f, totalVector);
            for( g = 0; g < EV ; g++){
                tAddTw(f1.f, totalVector, 0, f1.f.user+g, cmpl);
            }
            CanonicalRankDecomposition( f1.f, NULL, totalVector, 0, eigenVectors, cmpl, c1->rt.TOLERANCE, c1->rt.relativeTOLERANCE, c1->rt.ALPHA, c1->rt.THRESHOLD, c1->rt.MAX_CYCLE,c1->rt.XCONDITION, part(f1.f,eigenVectors));
        }
        tClear(f1.f, totalVector);
    }
#else
    f1.i.iRank = 1;
    f1.i.qFloor = 1;
    f1.i.Iterations = 1;
    f1.i.nStates =  1;
    iModel(c1,&f1);
    for (int  i = 0 ; i < f1.i.canonRank ; i++)
        tId(f1.f,eigenVectors,0);
    EV = 1;
#endif
    division target,OpSpiral = defSpiralMatrix(&f1.f, Iterator);
    target = f1.f.name[OpSpiral].name;
    
    if ( f1.i.Iterations == 2 ){
        {
            double norm = sqrt(pMatrixElement(f1.f, eigenVectors ,0,nullOverlap,0,eigenVectors ,0));
            if ( norm > c1->rt.THRESHOLD ){
                printf("for multiply, Normed from %f\n", norm );
                fflush(stdout);
                tScaleOne(f1.f, eigenVectors, 0, 1./norm);
            }
        }
        tHXpY( f1.f,  eigenVectors, target, c1->i.shiftFlag , eigenVectors, f1.f.rt->TOLERANCE,f1.f.rt->relativeTOLERANCE,f1.f.rt->ALPHA,f1.f.rt->THRESHOLD,f1.f.rt->MAX_CYCLE,f1.f.rt->XCONDITION,  f1.f.name[eigenVectors].Partition,f1.f.name[eigenVectors].Partition*0.8);
    }
    if ( ((((f1.i.filter/4)%2)==1) * f1.i.irrep) ){
        tFilter(f1.f, 1, (((f1.i.filter/4)%2)==1) * f1.i.irrep, eigenVectors);
    }
    if (f1.i.Iterations == 1 ){
        double norm = sqrt(pMatrixElement(f1.f, eigenVectors ,0,nullOverlap,0,eigenVectors ,0));
        if ( norm > c1->rt.THRESHOLD ){
            printf("Normed from %f\n", norm );
            fflush(stdout);
            tScaleOne(f1.f, eigenVectors, 0, 1./norm);
        }
    }
    printExpectationValues(c1, f1.f, Iterator, eigenVectors);
    print(c1,f1,1,0,1,eigenVectors);
    fModel(&f1.f);
    return 0;
}

/**
 *Build ritz matrix, eigenSolve, print out
 *
 *Can parallelize build via OpIndex being non -1
*/
inta ritz(   calculation * c1,   field f1){
    char filename[MAXSTRING];    char str[SUPERMAXSTRING];

    inta fi,EV = 0;
    if ( ! allowQ(f1.f.rt,blockMatrixElementsblock)){
        printf("blockMatrixElementsblock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    f1.i.qFloor = countLinesFromFile(c1,f1,0,&f1.i.iRank, &f1.i.xRank);
    iModel(c1,&f1);
        
    for ( fi =0 ; fi < f1.i.files ; fi++){
        tLoadEigenWeights (c1,f1, f1.i.fileList[fi],&EV, f1.f.user,f1.i.collect);
    }
    if (EV == 0 ){
        printf ("no ritz vectors!\n");
        exit(0);
    }
    inta stride = f1.f.maxEV;
    if(c1->i.build){
        tBuildMatrix(0 , f1.f,Ha, f1.f.user,EV);
            
        if ( f1.i.OpIndex != -1 ){
        //WRITE OUT V
            tFilename(c1->name, 2, 0, 0, 0, filename);
            fflush(stdout);
            ioArray(c1, f1, filename, stride*stride, (mea*)myStreams(f1.f,matrixHbuild,0), 0);
            
            tFilename(c1->name, 1, 0, 0, 0, filename);
            fflush(stdout);
            ioArray(c1, f1, filename, stride*stride, (mea*)myStreams(f1.f,matrixSbuild,0), 0);
            
            
        }
    }else {
        myZero(f1.f, matrixHbuild, 0);
        myZero(f1.f, matrixSbuild, 0);

        //LOAD V
        for ( fi =0 ; fi < f1.i.matrices ; fi++){
            tFilename(f1.i.matrixList[fi], 2, 0, 0, 0, filename);
            ioArray(c1, f1, filename, stride*stride, (mea*)myStreams(f1.f,matrixHbuild,0)+stride*stride, 1);
            cblas_daxpy(stride*stride,1., myStreams(f1.f,matrixHbuild,0)+stride*stride, 1, myStreams(f1.f,matrixHbuild,0), 1);
        }
        tFilename(f1.i.matrixList[0], 1, 0, 0, 0, filename);
        ioArray(c1, f1, filename, stride*stride, (mea*)myStreams(f1.f,matrixSbuild,0), 1);
    }
    
    
    if ( f1.i.OpIndex <= 0 ){
        mea *V = (mea*)myStreams(f1.f,matrixHbuild,0);

        tSolveMatrix(1, f1.f, f1.i.nStates, f1.f.user, EV, twoBodyRitz);
        inta iii,ii,stride = f1.f.maxEV;
        for ( iii = 0; iii < f1.i.nStates ; iii++){
            {
                FILE * outf ;
                sprintf(str, "%s-%d.vector",c1->name,iii+1);
                outf = fopen (str,"w");
                fclose(outf);
                sprintf(str, "%s-%d",c1->name,iii+1);
            }
            for ( ii = 0; ii < EV ; ii++)
                printVector(c1,f1.f,f1.f.name[f1.f.user+ii].value.title,str,f1.f.name[f1.f.user+ii].value.stage-1,f1.i.irrep, V+stride*iii+ii);
        }
    }
    fModel(&f1.f);

    return 0;
}

/**
 *from prompt
 */
int run (inta argc , char * argv[]){
    argc--;///erase ./andromeda...
    argv++;
      calculation c;
      field f;

    
    if ( argc > 0 ){
        switch ( atoi( argv[0])){
            case -1 :
                argc--;
                argv++;
                ///andromeda -1 inputFile
                ///runs normally from file
                bootShell(argc, argv,&c,&f);
                break;

            case 1 :
                argc--;
                argv++;
                ///andromeda 1 inputFile
                ///spits out memory requirements only
                bootShell(argc, argv,&c,&f);
                c.i.RAMmax = 0;
                break;

            case 0 :
                ///andromeda 0
                printf("----\nv9.4\n\n%s\n\n",getenv("LAUNCH"));
                printf("cat file |  andromeda 1 \n\t\t--> MEMORY AND TERM output without committing\n");
                printf("cat file |  andromeda  \n");
                printf("andromeda -1 file \n");
                exit(0);
        }

    }else {
        ///andromeda inputFile
        bootShell(argc, argv,&c,&f);
    }
    defineCores(&c,&f);
    assignCores(f.f,1 );
    if ( c.rt.phaseType == buildFoundation ){//0
#ifdef SPHERE
        foundationS(&c,f);
#else
        foundationB(&c,f);
#endif
    }
    else if ( c.rt.phaseType == productKrylov ){
        singlekrylov(&c,f);
    }
    else if ( c.rt.phaseType == solveRitz ){
        ritz(&c,f);
    }
    else if ( c.rt.phaseType == buildQuad ){
        quads(&c,f);
    }

    ///Scripts look for this ...
    printf("\n\nFINIS.\n\n");
    return 0;
}


/**
 * basic main calling run
 */
int main (inta argc , char * argv[]){
        return run(argc,argv);
}
