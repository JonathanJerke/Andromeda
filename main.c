/*
 *  main.c
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

#include "main.h"

struct calculation bootShell (INT_TYPE argc , char * argv[]){
#ifndef APPLE

#ifdef GSL_LIB
    gsl_set_error_handler_off ();
#endif
    
    INT_TYPE broke;
    
    time_t t;
    /* Intializes random number generator */
    srand((unsigned) time(&t));
    
    struct calculation c1 = initCal();
    INT_TYPE i,c,EV,EV2,ER;
    
    FILE * in = stdin;
    char str[MAXSTRING];
    initCalculation(&c1);
    broke = readInput(&c1,in );
    if ( broke )
        exit(1);
    finalizeInit(&c1);
#else
    struct calculation c1 =initCal();
#endif

    return c1;
}

INT_TYPE print(struct calculation *c ){
    INT_TYPE irrep;
    struct field * f1 = & c->i.c;
    char str[MAXSTRING];
    INT_TYPE iii,jjj=1,cmpl;
    for ( irrep = 0 ; irrep <= 5 ; irrep++){
        jjj = 1;
        for ( iii = 0; iii < c->i.heliumFlag  ; iii++)
            if ( f1->sinc.tulip[eigenVectors+iii].symmetry  == irrep && ((! c->i.irrep)|| c->i.irrep == irrep)){
                
                printf("%dState%d:%d:,%d ,%1.15f, %d, %d , %1.1f,%1.15f\n", f1->body,jjj, f1->sinc.N1,jjj,f1->sinc.tulip[eigenVectors+iii].value,bodies(f1,eigenVectors+iii),irrep, deg(f1, irrep),f1->sinc.tulip[eigenVectors+iii].value2);
                if ( (c->i.outputFlag) % 2 == 1){
                    
                    for ( cmpl = 0 ; cmpl < spins(f1, eigenVectors+iii) ; cmpl++)
                    {
#ifndef APPLE
                        tFilename(c->name,jjj,bodies(f1, eigenVectors+iii) ,irrep, cmpl,str);

                        FILE * out = fopen ( str,"w" );
                        if ( out != NULL ){
                            outputFormat(&c->i.c, out, eigenVectors+iii,cmpl  );
                            fclose(out);
                        }
#endif
                    }
                }
                jjj++;
                
            }
    }
    if ( (c->i.outputFlag) % 2 == 1)
        tEdges(c);
    fflush(stdout);
    return 0;
}


INT_TYPE exec (struct calculation *c ){
    INT_TYPE splitFlag = 0,cycleFlag = 0;;
    if ( MaxCore > hardMaxCore ){
        printf("max core at this time is %d\n", hardMaxCore);
        exit(0);
    }
#ifdef splitTag
    if ( c->rt.body >= two)
    splitFlag = 1;
#endif
    struct calculation *c1 = malloc(sizeof(struct calculation));
    *c1 = *c;
    
    INT_TYPE a,plusSize,nStatesTrans=0,nStatesFound=0 ,RdsSize = 0,totalIter = 0;
    FILE * out = stdout;
    struct runTime * rt = & c1->rt;

    struct field *f1 = &(c1->i.c);
    if ( c1->i.iCharge > 0 )
        f1->Ne = c1->i.iCharge;
    
    
    
    double max1,max2=0;
    //fprintf(out , " lattice\t%1.3f\t%lld\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->rt.body);
    fprintf( out,   "  Number of electrons: %d\n  Number of Atoms: %d\n\n",c1->i.c.Ne,f1->Na);
    if ( c1->i.Angstroms )
    {
        printf("Angstroms\n");

    }else{
        printf("Bohrs\n");
    }
    for ( a = 1 ; a <= c1->i.c.Na ;a++){
        if ( c1->i.Angstroms )
            fprintf(out,"  %d %1.6f %1.6f %1.6f ",c1->i.c.atoms[a].label.Z,a0*c1->i.c.atoms[a].position[1],a0*c1->i.c.atoms[a].position[2],a0*c1->i.c.atoms[a].position[3]);
        else
            fprintf(out,"  %d %1.6f %1.6f %1.6f ",c1->i.c.atoms[a].label.Z,c1->i.c.atoms[a].position[1],c1->i.c.atoms[a].position[2],c1->i.c.atoms[a].position[3]);
        fprintf(out, "\n");
        max1 = max(max(c1->i.c.atoms[a].position[1],c1->i.c.atoms[a].position[2]),c1->i.c.atoms[a].position[3]);
        if ( max1 > max2 )
            max2 = max1;
        
        
    }
    
    INT_TYPE upCore,i;
    INT_TYPE nSlot ;
#ifdef OMP
    upCore = c1->i.omp;
#pragma omp parallel for private (i)

    for ( i = 0; i < MaxCore ; i++){
        nSlot = omp_get_num_threads();
    }
#else
    nSlot = 1;
#ifndef APPLE
#ifdef OMP
    rt->NCore = 1;
#endif
#endif
#endif

    
#ifdef OMP
    rt->NCore = upCore;
    rt->NSlot = nSlot;
    
    if ( rt->NSlot <= MaxCore ){
    }
    else{
        printf("MAX %lld\n", MaxCore);
        rt->NSlot = MaxCore;
    }
    if ( rt->NCore > rt->NSlot ){
        printf("lanes > available\n");
        rt->NCore = rt->NSlot;
    }
    printf("lanes %lld\n", rt->NCore);
    printf("slots %lld\n", rt->NSlot);
    fflush(stdout);

#endif
    
    if ( !(!c1->rt.printFlag)){
        c1->i.nStates = c->i.heliumFlag;
        iModel(c1);
        fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->rt.body);
#ifndef APPLE
        INT_TYPE lV ;
        if ( (c1->rt.printFlag/8)%2 ){
            lV = tLoadEigenWeightsWithConstraints ( c1, c1->mem.fileList,c1->mem.constraintFile);
        } else {
            lV = tLoadEigenWeights ( c1, c1->mem.fileList, eigenVectors);
        }
        if ( ( c1->rt.printFlag/4) % 2 ){
            printf("OPERATOR = KINETIC\n");
            ///here a text file input could be translated into a collection of operators.
            f1->sinc.tulip[printOperator].name = Ha;
            if ( bodies(f1, eigenVectors ) == two)
                f1->sinc.tulip[kinetic2].linkNext = nullName;
            else if ( bodies(f1, eigenVectors ) == three)
                f1->sinc.tulip[kinetic3].linkNext = nullName;
            else if ( bodies(f1, eigenVectors ) == four)
                f1->sinc.tulip[kinetic4].linkNext = nullName;
        }
        
        printOutput ( f1,lV);
        printVectorOutput ( f1,lV);
#endif
    } else {
        INT_TYPE xyz,cycle =0,iteration,iter=0,EV,type2;
        INT_TYPE space;
        enum division usz;
        if ( c1->i.sectors && ! c1->i.densityFlag){
            c1->i.heliumFlag = c->i.nTargets;
            c1->i.nStates = c1->i.heliumFlag;
            if ( splitFlag )
                c1->i.nStates *= tPaths(c1->rt.body , c1->i.irrep);
            c1->i.Iterations = 1;
            c1->rt.maxEntropy = 1;
            iModel(c1);
            if ( c1->rt.body == one )
                t1BodyConstruction(c1, build);
            else
                tSAboot(c1);
            usz = eigenVectors+c1->i.nStates;
            if ( tCollect(f1,c1->i.irrep,usz,c1->i.qFloor ,c->i.seekPower) != c1->i.qFloor){
                printf("could not muster \n");
                exit(0);
            }
            EV =c1->i.qFloor;
            cycleFlag = 1;
        }else if (! c1->i.sectors && ! c1->i.densityFlag  ){
            {
                INT_TYPE lines = 0;
                char ch ;
                FILE * fp = fopen(c1->mem.fileList,"r");
                if ( fp == NULL ) {
                    printf("file?\n");
                    exit(0);
                }
                while(!feof(fp))
                {
                    ch = fgetc(fp);
                    if(ch == '\n')
                    {
                        lines++;
                    }
                }
                fclose(fp);
                c->i.nTargets = lines;
            }
            c1->i.iRank = c->i.bRank;
            c1->i.qFloor = c->i.nTargets;
            c1->i.heliumFlag = c->i.nTargets;
            c1->i.nStates = c->i.nTargets;

            iModel(c1);

#ifndef APPLE
            if ( tLoadEigenWeights ( c1, c1->mem.fileList, eigenVectors) != c1->i.qFloor ){
                printf("set helium %d \n", c1->i.qFloor);
                exit(0);
            }
            tFilter(&c1->i.c, c1->i.nStates, 0, eigenVectors);
#endif
            usz = eigenVectors+c1->i.nStates;
            EV = xConstructFoundation (c1 , usz, c1->i.qFloor, c1,   eigenVectors,   c1->i.nStates ,1);

        }else if ( c1->i.sectors && c1->i.densityFlag ){
            {
                c1->i.heliumFlag = c->i.nTargets;
                c1->i.nStates = c1->i.heliumFlag;

                size_t ms = MAXSTRING;
                char input_line[MAXSTRING];
                char * mask = input_line;
                
                INT_TYPE lines = 0;
                char ch ;
                FILE * fp = fopen(c1->mem.densityName,"r");
                if ( fp == NULL ) {
                    printf("file?\n");
                    exit(0);
                }

                while(!feof(fp))
                {
                    ch = fgetc(fp);
                    if(ch == '\n')
                    {
                        lines++;
                    }
                }
                fclose(fp);
                c1->i.densityFlag = lines;
                printf("density Count %d\n",c1->i.densityFlag );

            }            
            iModel(c1);
            if ( c1->i.bodyDensity != bodies (&c1->i.c,c1->i.c.sinc.density) ){
                printf("body count!\n");
                exit(0);
                
            }
#ifndef APPLE
            if ( c1->i.densityFlag )
                printf("density Count %d\n",c1->i.densityFlag );
            
            if ( tLoadEigenWeights ( c1, c1->mem.densityName ,f1->sinc.density) != c1->i.densityFlag){
                printf("set helium %d \n", c1->i.densityFlag);
                exit(0);
            }

#endif
            if ( c1->rt.body == one )
                t1BodyConstruction(c1, build);
            else
                tSAboot(c1);

            usz = eigenVectors+c1->i.nStates;
            if ( tCollect(f1,c1->i.irrep,usz,c1->i.qFloor ,c->i.seekPower) != c1->i.qFloor){
                printf("could not muster \n");
                exit(0);
            }
            EV =c1->i.qFloor;
            cycleFlag = 1;

        }else {
            exit(0);
        }
        RdsSize = EV;
        c1->i.sectors = 0;
        for ( totalIter = 0,cycle = 0, iter = 1; cycle <= c->i.cycles ;totalIter++)
        {
            
            
            
            printf("Iter %d \n Cycle %d\n Iter %d\n cycleFlag %d\n\n", totalIter, cycle ,iter, cycleFlag);
            printf("EIGENVALUE %d @@ %d , %d \nUSZ %d @@ %d , %d \n", eigenVectors,c1->i.heliumFlag,c1->i.nStates, usz ,EV,RdsSize );

            fflush(stdout);
            
            

             {
                if ( tEigenCycle(&c1->i.c,Ha,'T', c1->i.heliumFlag, usz,RdsSize,0, EV,c1->i.irrep,  (3+!splitFlag),eigenVectors,twoBodyRitz))
                {
                    printf("usz:%d--ld %d %d\n", usz, RdsSize, EV);
                    
                    if (1 ){
                        printf("crap!");
                        fflush(stdout);
                        exit(0);
                    }
                }
                 tFilter(&c1->i.c, c1->i.nStates, !(!c1->i.filter )* c->i.irrep, eigenVectors);
                 
            }

            if (!  c->i.Iterations  )
                break;
            

            if ( (cycleFlag || iter  == c1->i.Iterations )){
                if ( cycle == c->i.cycles )
                    break;
                
                
                cycleFlag = 0;
                iter = 1;
                struct calculation *c2 = malloc(sizeof(struct calculation ));
                *c2 = *c1;
                c2->i.Iterations = c->i.Iterations;
                c2->rt.maxEntropy = c->rt.maxEntropy;
                c2->i.nStates = c->i.heliumFlag;
                c2->i.iRank = c1->i.bRank;
                
                
                c2->i.c.mem1->bootedMemory =0;
                c2->i.c.sinc.tulip = NULL;
                for ( space = 0; space <= SPACE ; space++)
                    c2->i.c.sinc.rose[space].stream = NULL;
                
                
                
                
                
                
                cycle++;
                c2->i.bRank++;
                INT_TYPE step = c1->i.cycleStep;
                c2->i.epi  += step;
                c2->i.d *= pow( (2.* c1->i.epi + 1.) /(2.*c2->i.epi + 1. ),c2->i.attack);
                
                c2->i.vectorMomentum *= (c1->i.d* (2.* c1->i.epi + 1.) )/(c1->i.d*(2.*c2->i.epi + 1. ));
                printf("attack %f %f -> %f\n",c2->i.attack,c1->i.d, c2->i.d);
                
                {
                    {
                        
                        {
                            nStatesTrans = 0;
                            nStatesFound = 0;
                            
                            INT_TYPE ii = 0,iii,irrep;
                            for ( irrep = 0 ; irrep <= 5 ; irrep++){
                                for ( iii = 0; iii < c1->i.nStates ; iii++)
                                {
                                    if (f1->sinc.tulip[eigenVectors+iii].symmetry == irrep && ((! c1->i.irrep) || irrep == c1->i.irrep) ){
                                        ii++;
                                        
                                        if ( ii <= c->i.nTargets){
                                            nStatesTrans = imax(nStatesTrans,iii+1);
                                            nStatesFound = ii;
                                        }
                                    }
                                }}
                            
                            c2->i.heliumFlag = nStatesFound;
                            c2->i.nStates = nStatesFound;
                        }
                        
                        
                        
                        INT_TYPE ii;
                        {
                            INT_TYPE ii,iii,iv=0;
                            for ( iii = 0; iii < nStatesTrans ; iii++){
                                ii = xConstructFoundation (c2 , 0,0, c1,eigenVectors+iii,1 ,1);
                                iv += ii;
                            }
                            printf("floor %d \n Found %d \n Trans %d\n", iv, nStatesFound , nStatesTrans);
                            fflush(stdout);
                            c2->i.qFloor = iv;
                            EV = iv;
                            RdsSize = EV;
                        }
                        
                        if ( c1->i.densityFlag )
                            c2->i.densityFlag = c1->i.densityFlag;
                        
                        iModel(c2);
                        xConstructFoundation (c2 , eigenVectors, c2->i.nStates, c1,eigenVectors,nStatesTrans ,1);
                        fModel(c1);

                        free(c1);
                        c1 = c2;
                        f1 = &c1->i.c;
                        usz = eigenVectors+c1->i.nStates;
                        
                        if ( c1->i.densityFlag ){

                        if ( c1->i.bodyDensity != bodies (&c1->i.c,c1->i.c.sinc.density) ){
                            printf("body count!\n");
                            exit(0);
                            
                        }
                            printf("density Count %d\n",c1->i.densityFlag );
                        
                        if ( tLoadEigenWeights ( c1, c1->mem.densityName ,f1->sinc.density) != c1->i.densityFlag){
                            printf("set helium %d \n", c1->i.densityFlag);
                            exit(0);
                        }

                        }
#ifdef OMP
#pragma omp parallel for private (ii) schedule(dynamic,1)
#endif
                        for ( ii = 0; ii < c1->i.nStates ; ii++){
                            tEqua(f1, usz+ii, 0, eigenVectors+ii , 0);
                            tEqua(f1, usz+ii, 1, eigenVectors+ii , 1);
//tClear(f1, eigenVectors+ii);
                            
                        }

                        EV =c1->i.nStates;
                        RdsSize = EV;
                    }
                }
                usz = eigenVectors + c1->i.nStates;

            }
            
                
                
            if ( iter < c1->i.Iterations)//extra catch-all
            
            {
                RdsSize += tGreatDivideIteration(f1,Ha, 1,0,usz+RdsSize-EV,EV,2*EV,0)-EV;

                while ( RdsSize > EV*(c1->i.lookBack+1)){
                    RdsSize -= EV;
                    usz += EV;
                    tClear(f1,matrixHbuild);
                    printf("forgetting...\n");
                    fflush(stdout);
                }
                
                
                if   ( tEigenCycle(&c1->i.c,Ha,'T', c1->i.heliumFlag, usz,RdsSize,0, EV, 0,0,eigenVectors,twoBodyRitz))
                {
                    RdsSize -= EV;
                    usz += EV;
                    if ( RdsSize == 0 )
                    {
                        printf("not good!\n");
                        fflush(stdout);
                        exit(0);
                    }
                    cycleFlag = 1;
                    
                    
                    INT_TYPE ii = 0,iii;
                    for ( iii = 0; iii < c1->i.heliumFlag ; iii++)
                    {
                        if (1){
                            fprintf(out,"Failure%d:%d:,%lld ,%1.15f\n", ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value);
                            ii++;
                        }
                    }
                    fflush(stdout);

                    
                    continue;//break from cycling/iterating...if failed LD test.
                }
                {
                    INT_TYPE iii ;
                    for ( iii = 0; iii < c1->i.heliumFlag ; iii++){
                        printf ( "\n\n Vector \t%d\n", iii+1);
                        printExpectationValues(f1, Ha, eigenVectors+iii);
                    }
                }
            }else {
                cycleFlag = 1;//not likely to happend unless errro.
                printf("woops! ... looks like you accidently attempted an extra iteration\n");
            }
            fprintf(out,"EV %d ; irrep %d ; offset %d\n", EV, c1->i.irrep,RdsSize-EV);
            fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->rt.body);
            fflush(stdout);
            
//            INT_TYPE ii = 0,iii;
//            for ( iii = 0; iii < c1->i.heliumFlag ; iii++)
//            {
//                if (1){
//                    fprintf(out,"Condition%d:%d:,%lld ,%1.15f\n", ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value);
//                    ii++;
//                }
//            }
//            fflush(stdout);
            iter++;

        }
        

    }
//    if ( splitFlag )
//    {//exit
//        INT_TYPE nP = tPerms(f1->body);
//        enum division usz;
//        INT_TYPE EV,space,sp,g,i,ii,iii = 0,rank=0,iv,sup;
//
//        struct calculation *c2 = malloc(sizeof(struct calculation ));
//        *c2 = *c1;
//        c2->i.Iterations = c->i.Iterations;
//        c2->rt.maxEntropy = c->rt.maxEntropy;
//        c2->i.nStates = c->i.heliumFlag;
//        c2->i.iRank = c1->i.bRank;
//
//
//        c2->i.c.mem1->bootedMemory =0;
//        c2->i.c.sinc.tulip = NULL;
//        for ( space = 0; space <= SPACE ; space++)
//            c2->i.c.sinc.rose[space].stream = NULL;
//
//
//
//
//
//        struct calculation *c3 = malloc(sizeof(struct calculation ));
//        *c3 = *c1;
//        c3->i.Iterations = 1;
//        c3->i.cycles = 0;
//        c3->i.nStates = c1->i.nStates;
//        if ( c->i.irrep == 3 )
//            c3->i.nStates *= 4;
//        if ( c->i.irrep > 3 )
//            c3->i.nStates *= 9;
//        c3->i.Iterations = 0;
//        c3->i.qFloor = 0;
//        c3->i.c.mem1->bootedMemory =0;
//        c3->i.c.sinc.tulip = NULL;
//        for ( space = 0; space <= SPACE ; space++)
//            c3->i.c.sinc.rose[space].stream = NULL;
//
//        iModel(c3);
//
//
//
//        for( g = 0; g < nSAG*nSAG*nSAG; g++)
//            c2->i.cSA[g] = 0;
//
//
//        for ( ii = 0; ii < c1->i.nStates ; ii++)
//        {
//            sup = tSizeUp(rank,f1,c1->i.irrep, eigenVectors+ii);
//            if ( sup ){
//                for ( sp = 0; sp  < spins (f1, eigenVectors+ii);sp++)
//                    xEqua(&c3->i.c, copyVector, sp, f1, eigenVectors+ii, sp);
//                if (tSelect(&c3->i.c,  iii,c3->i.irrep+nP, eigenVectors, copyVector,1)){
//                    c3->i.c.sinc.tulip[eigenVectors+iii].path  = tPath(f1, eigenVectors+ii );
//                    c2->i.cSA[tPath(f1, eigenVectors+ii ) ] += 1;
//                    iii++;
//                }
//            }
//        }
//        fModel(c1);
//        free(c1);
//
//        {
//            INT_TYPE nG = tSize(f1->body);
//            for( g = 0; g < nSAG*nSAG*nSAG; g++)
//                if ( c2->i.cSA[g] )
//                    printf("%d : final : %d:: (%d %d %d) --> %d\n",g, c2->i.cSA[g],c2->i.cSA[g]%nG+1 , (c2->i.cSA[g]/nG)%nG +1, (c2->i.cSA[g]/(nG*nG) % nG) +1, c->i.irrep);
//
//        }
//        c2->i.qFloor = iii;
//
//        c2->i.heliumFlag = c->i.nTargets;//number of states
//        c2->i.nStates = c2->i.heliumFlag;//number of paths ...
//        EV = c2->i.qFloor;
//
//        iModel(c2);
//        usz = eigenVectors+c2->i.nStates;
//        {
//            for ( i = 0; i < c2->i.qFloor ; i++)
//                for ( sp = 0; sp < spins(f1, usz); sp++)
//                    xEqua(&c2->i.c, usz+i, sp, &c3->i.c, eigenVectors+i, sp);
//
//        }
//        fModel(c3);
//        free(c3);
//
//        c1 = c2;
//        f1 = &c1->i.c;
//
//        RdsSize = EV;
//        {
//            if ( tEigenCycle(&c1->i.c,Ha,'T', c1->i.heliumFlag, eigenVectors + c1->i.nStates,RdsSize,0, EV,c1->i.irrep, 4,eigenVectors,twoBodyRitz))
//            {
//                printf("ld %d %d\n", RdsSize, EV);
//
//                if (1 ){
//                    printf("crap!");
//                    exit(0);
//                }
//            }
//            tFilter(&c1->i.c, c1->i.nStates, 0, eigenVectors);
//
//        }
//    }//exit
    
    *c = *c1;
    {
        INT_TYPE space ;
        c->i.c.sinc.tulip = c1->i.c.sinc.tulip;
        for ( space = 0; space <= SPACE ; space++){
            c->i.c.sinc.rose[space].stream =c1->i.c.sinc.rose[space].stream;
            c1->i.c.sinc.rose[space].stream = NULL;
        }
    }
    c1->mem.bootedMemory = 0;
    c1->i.c.sinc.tulip = NULL;
    fModel(c1);
    free(c1);
    return 0;
}
#if 1
int main (INT_TYPE argc , char * argv[]){
    
#if sincFlag

//    tTestSA(two, 2);
//    tTestSA(three, 3);
//    tTestSA(four, 5);
//    exit(0);
    
    struct calculation c = bootShell(argc, argv);
    exec(&c);
#ifndef splitTag
    print(&c);
#endif
    fModel(&c);
#else
    

    struct function_label func1 ;
    func1.interval = 1;
    func1.fn = Coulomb;//Yukawa, Erf
    func1.param[0] = 1.;//scalar out front
    func1.param[1] = 1.;//not used by elemCal...a quadrature parameter
    func1.param[2] = 1.;//parameter in Erf and Yukawa (scale/mass)
    func1.param[3] = 1.;//second a parameter, unused.
//    getDescription ( &func1 ,0.,stdout);

    INT_TYPE periodic = 0;//GTO's not yet coded for peridicity/.
    INT_TYPE i ;
    
#if 1
    if ( argc != 4 ){
        printf("L,Distance, Gaussian\n");
        exit(0);
        
    }
    
    
    INT_TYPE l = atoi(argv[1]);
    double x = atof(argv[2]);
    double b = atof(argv[3]);
#else
    INT_TYPE l =0;
    double x = 0;
    double b = 1;
#endif
    
    INT_TYPE acc = 1;//somehwo?
    
    INT_TYPE body = 2;
    struct general_2index g3[3];

    INT_TYPE l1=l,l2=l,l3=l,l4=l;
//    for ( l1 = 0 ; l1 <=5 ; l1++ )
//        for ( l2 = 0 ; l2 <= 5 ; l2++)
//            for (l3 = 0; l3<=5 ; l3++)
//                for (l4 = 0 ; l4<=5 ;l4++){
//
    {
                    {
                        
                        i = 0;
                        g3[i].powSpace = 0;
                        g3[i].gaussianAccelerationFlag = acc;
                        g3[i].realFlag = 1;
                        g3[i].point = 2;
                        g3[i].i[0].action = 0;//derivatives
                        g3[i].i[1].action = 0;//derivatives
                        g3[i].momentumShift = 0;
                        
                        g3[i].i[0].bra.basis = GaussianBasisElement;
                        g3[i].i[0].bra.length = b;//body 1 Matrix
                        g3[i].i[0].bra.origin = 0;
                        g3[i].i[0].bra.index = l1;
                        g3[i].i[0].bra.periodic = 0;
                        
                        g3[i].i[0].ket.basis = GaussianBasisElement;
                        g3[i].i[0].ket.length = b;
                        g3[i].i[0].ket.origin = 0 ;
                        g3[i].i[0].ket.index = l2 ;
                        g3[i].i[0].ket.periodic = 0;
                        
                        
                            g3[i].i[1].bra.basis = GaussianBasisElement;
                            g3[i].i[1].bra.length = 1;//body 1 Matrix
                            g3[i].i[1].bra.origin = 0;
                            g3[i].i[1].bra.index = l3;
                            g3[i].i[1].bra.periodic = 0;
                            
                            g3[i].i[1].ket.basis = GaussianBasisElement;
                            g3[i].i[1].ket.length = 1;
                            g3[i].i[1].ket.origin = 0;
                            g3[i].i[1].ket.index = l4;
                            g3[i].i[1].ket.periodic = 0;
                        g3[i].fl = & func1;
                        
                        
                        
                        
                        
                        if ( i == 0 )
                            g3[i].periodic = ( periodic ) % 2;
                        else if ( i == 1 )
                            g3[i].periodic = ( periodic / 2 ) % 2;
                        else if (i == 2 )
                            g3[i].periodic = ( periodic / 4) % 2;
                        else
                        {
                            printf("here\n");
                            exit(0);
                        }
                        g3[i].body = body;
                        
                        //                printf("%f %f\n", creal(aaGdnGdm(atoi(argv[3]),atoi(argv[4]), &g3[0].i[0])));
                        //                exit(0);
                        
                        
                        //
//                        DCOMPLEX va;
//                        for ( i =0; i< 1; i++)
//                            va = aaGGCGG(atof(argv[1]), &g3[0]);
//
//                        printf("%d\t%d\t%d\t%d\t %f\t%f\n",l1,l2,l3,l4, creal(va), cimag(va));
//                        //        exit(0);
                        //
                    }
                }
    {
        i = 1;
        g3[i].powSpace = 0;
        g3[i].gaussianAccelerationFlag = acc;
        g3[i].realFlag = 1;
        g3[i].point = 2;
        g3[i].i[0].action = 0;//derivatives
        g3[i].i[1].action = 0;//derivatives
        g3[i].momentumShift = 0;

        g3[i].i[0].bra.basis = GaussianBasisElement;
        g3[i].i[0].bra.length = b;//body 1 Matrix
        g3[i].i[0].bra.origin = x;
        g3[i].i[0].bra.index = l;
        g3[i].i[0].bra.periodic = 0;
        
        
        g3[i].i[0].ket.basis = GaussianBasisElement;
        g3[i].i[0].ket.length = b;
        g3[i].i[0].ket.origin = x;
        g3[i].i[0].ket.index = l;
        g3[i].i[0].ket.periodic = 0;
        
        g3[i].i[1].bra.basis = GaussianBasisElement;
        g3[i].i[1].bra.length = b;//body 1 Matrix
        g3[i].i[1].bra.origin = 0;
        g3[i].i[1].bra.index = l;
        g3[i].i[1].bra.periodic = 0;
        
        g3[i].i[1].ket.basis = GaussianBasisElement;
        g3[i].i[1].ket.length = b;
        g3[i].i[1].ket.origin = 0;
        g3[i].i[1].ket.index = l;
        g3[i].i[1].ket.periodic = 0;

        g3[i].fl = & func1;
        
        if ( i == 0 )
            g3[i].periodic = ( periodic ) % 2;
        else if ( i == 1 )
            g3[i].periodic = ( periodic / 2 ) % 2;
        else if (i == 2 )
            g3[i].periodic = ( periodic / 4) % 2;
        else
        {
            printf("here\n");
            exit(0);
        }
        g3[i].body = body;
    }
    {
        i = 2;
        g3[i].powSpace = 0;
        g3[i].gaussianAccelerationFlag = acc;
        g3[i].realFlag = 1;
        g3[i].point = 2;
        g3[i].i[0].action = 0;//derivatives
        g3[i].i[1].action = 0;//derivatives
        g3[i].momentumShift = 0;

        g3[i].i[0].bra.basis = GaussianBasisElement;
        g3[i].i[0].bra.length = b;//body 1 Matrix
        g3[i].i[0].bra.origin = x;
        g3[i].i[0].bra.index = l;
        g3[i].i[0].bra.periodic = 0;
        
        
        g3[i].i[0].ket.basis = GaussianBasisElement;
        g3[i].i[0].ket.length = b;
        g3[i].i[0].ket.origin = x;
        g3[i].i[0].ket.index = l;
        g3[i].i[0].ket.periodic = 0;
        
        g3[i].i[1].bra.basis = GaussianBasisElement;
        g3[i].i[1].bra.length = b;//body 1 Matrix
        g3[i].i[1].bra.origin = 0;
        g3[i].i[1].bra.index = l;
        g3[i].i[1].bra.periodic = 0;
        
        g3[i].i[1].ket.basis = GaussianBasisElement;
        g3[i].i[1].ket.length = b;
        g3[i].i[1].ket.origin = 0;
        g3[i].i[1].ket.index = l;
        g3[i].i[1].ket.periodic = 0;


        g3[i].fl = & func1;
        
        if ( i == 0 )
            g3[i].periodic = ( periodic ) % 2;
        else if ( i == 1 )
            g3[i].periodic = ( periodic / 2 ) % 2;
        else if (i == 2 )
            g3[i].periodic = ( periodic / 4) % 2;
        else
        {
            printf("here\n");
            fflush(stdout);
            exit(0);
        }
        g3[i].body = body;
    }

    {
        INT_TYPE Nt = 10000;
        time_t start_t, lapse_t;
        time(&start_t);

        
        printf("computing 10,000 matrix elements...\n\n");
        for ( i = 0; i < Nt ; i++){
           elementCal(1e-3,-1, g3);
        }
        
        printf("\nnorms of gaussians %f %f %f,  %f %f %f,  %f %f %f , %f %f %f\n", GTOnorm(g3[0].i[0].bra),GTOnorm(g3[0].i[0].ket) ,GTOnorm(g3[0].i[1].bra),GTOnorm(g3[0].i[1].ket)
               ,GTOnorm(g3[1].i[0].bra),GTOnorm(g3[1].i[0].ket) ,GTOnorm(g3[1].i[1].bra),GTOnorm(g3[1].i[1].ket)
               ,GTOnorm(g3[2].i[0].bra),GTOnorm(g3[2].i[0].ket) ,GTOnorm(g3[2].i[1].bra),GTOnorm(g3[2].i[1].ket));
        
        printf("%f\n", elementCal(1e-3,-1, g3)
               *GTOnorm(g3[0].i[0].bra)*GTOnorm(g3[0].i[0].ket) *GTOnorm(g3[0].i[1].bra)*GTOnorm(g3[0].i[1].ket)
               *GTOnorm(g3[1].i[0].bra)*GTOnorm(g3[1].i[0].ket) *GTOnorm(g3[1].i[1].bra)*GTOnorm(g3[1].i[1].ket)
               *GTOnorm(g3[2].i[0].bra)*GTOnorm(g3[2].i[0].ket) *GTOnorm(g3[2].i[1].bra)*GTOnorm(g3[2].i[1].ket)
               );
        time(&lapse_t);
        printf("\nFinal per \t %15.15f\n", difftime(lapse_t, start_t)/Nt);

    
    
}
#endif

}

#endif

INT_TYPE buildElectronFreeInteraction ( struct calculation * c1, enum division mat){
    double value;
    struct field * f1 = &c1->i.c;
    INT_TYPE g,space,r,i,j,n,m,N1 = f1->sinc.N1, N2 = f1->sinc.N1*f1->sinc.N1;

    if ( ! CanonicalRank(f1, interactionExchange, 0)){
        printf("canceled ee\n");
        return 0;
    }
    struct calculation c0 = *c1;
    for( g = 0; g < nSAG*nSAG*nSAG; g++)
        c0.i.cSA[g] = 0;

    
    c0.mem.bootedMemory = 0;
    c0.i.c.sinc.tulip = NULL;
    c0.i.c.sinc.rose[0].stream = NULL;
    c0.i.c.sinc.rose[1].stream = NULL;
    c0.i.c.sinc.rose[2].stream = NULL;
    c0.i.c.sinc.rose[3].stream = NULL;
    c0.i.nTargets = 1;
    c0.rt.printFlag = 0;
    c0.rt.body = one;
    c0.rt.bodyType = electron;
    c0.i.side = 4;
    c0.i.iCharge = 1;
    c0.i.irrep = 0;
    c0.i.filter = 0;
    c0.i.l2 = 10000000;
    c0.i.cycles = 0;
    c0.i.Iterations = 0;
    c0.i.nStates = 1;
    c0.i.heliumFlag = 1;
    INT_TYPE t = 0;
    if ( c0.rt.runFlag %2 )
        t++;
    if ( (c0.rt.runFlag/2) %2 )
        t++;
    if ( (c0.rt.runFlag/4) %2 )
        t++;
    if ( t == 3 )
        c0.i.qFloor = 18;
    else if ( t == 2)
        c0.i.qFloor = 6;
    else if ( t == 1)
        c0.i.qFloor = 3;

    c0.i.sectors = 1;
    c0.i.bRank = 1;//assume its a separable solution
    exec(&c0);
    tClear(f1, mat);
    //build jellium -electron from electron-electron
    for ( r = 0; r < CanonicalRank(f1, interactionExchange, 0); r++){
        for ( i = 0; i < N1 ; i++)
            for ( j = 0; j < N1; j++)
                for (space = 0; space < SPACE ;space++)
                    
                {
                    value = 0.;
                    for ( n = 0 ; n < N1 ; n++)
                        for ( m = 0 ; m < N1 ; m++)
                            value += streams(f1, interactionExchange, 0,space)[(N2*(N1*i+n)+(N1*j+m))+r*N2*N2]* streams(&c0.i.c,eigenVectors,0,space)[n]*streams(&c0.i.c,eigenVectors,0,space)[m];
                    streams(f1, mat,0,space)[N1*i+j+r*N2] =  value/(N1*N1*N1);
                }
    }
    double zero = 0.,sum,prod;
    
    for ( r = 0; r < CanonicalRank(f1, interactionExchange, 0); r++){
        prod = 1.;
        for (space = 0; space < SPACE ;space++){
            sum = 0.;
            for ( i = 0; i < N1 ; i++)
                for ( j = 0; j < N1; j++){
                    sum += streams(f1, mat,0,space)[N1*i+j+r*N2] *streams(&c0.i.c,eigenVectors,0,space)[i] *(streams(&c0.i.c,eigenVectors,0,space)[j]);
                }
            prod *= sum;
        }
        zero += prod;
    }
    double norm = 0.;
    prod = 1.;
    for (space = 0; space < SPACE ;space++){
        sum = 0.;
            for ( j = 0; j < N1; j++){
                sum += streams(&c0.i.c,eigenVectors,0,space)[j]*streams(&c0.i.c,eigenVectors,0,space)[j];
            }
        prod *= sum;
    }
    norm += prod;

    fModel(&c0);
        
    f1->sinc.tulip[mat].Current[0] = CanonicalRank(f1, interactionExchange, 0);
    
    tClear(f1, copy);
    tId(f1, copy,0);
    tScale(f1, copy,-traceOne(f1, mat, 0)/traceOne(f1, copy, 0));
    tAddTw(f1, mat, 0,copy,0);

    tClear(f1, quadCube);
    tId(f1, quadCube,0);
    tScale(f1, quadCube,-traceOne(f1, interactionExchange, 0)/traceOne(f1, quadCube, 0));
    tAddTw(f1, interactionExchange, 0,quadCube,0);

    printf ("background %1.15f\t :norm %f\n", zero,norm);
    
    
//    if ( bodies(f1, eigenVectors)==two){
    //        tClear(f1, copy);
    //        tId(f1, copy,0);
    //        tScale(f1, copy,-2*3.*traceOne(f1, mat, 0)/(N1*N1*N1)/sqr(f1->Ne));
    //        //1.5 = -2 * 2 + 1
    //        tAddTw(f1, mat,0, copy,0);

//    }
//    else     if ( bodies(f1, eigenVectors)==three){
//        tClear(f1, copy);
//        tId(f1, copy,0);
//        tScale(f1, copy,-3*6.*traceOne(f1, mat, 0)/(N1*N1*N1)/sqr(f1->Ne));
//        //2 = -3 * 3 + 3
//        tAddTw(f1, mat,0, copy,0);
//    }
//    else     if ( bodies(f1, eigenVectors)==four){
//        tClear(f1, copy);
//        tId(f1, copy,0);
//        tScale(f1, copy,-4*10.*traceOne(f1, mat, 0)/(N1*N1*N1)/sqr(f1->Ne));
//        //2.5 = -4 * 4 + 6
//        tAddTw(f1, mat, 0,copy,0);
//    }
    
    f1->sinc.tulip[mat].stop[0][0] = CanonicalRank(f1, mat,0);
    
    return 0;
}
