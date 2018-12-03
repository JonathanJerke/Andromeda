/*
 *  main.c
 *
 *
 *  Copyright 2018 Jonathan Jerke and Bill Poirier.
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
    INT_TYPE type;
    struct field * f1 = & c->i.c;
        char str[MAXSTRING];
        INT_TYPE iii,jjj=1,cmpl;
    for ( type = 1 ; type <= 24 ; type++){
        jjj = 1;
        for ( iii = 0; iii < c->i.heliumFlag  ; iii++)
            if ( f1->sinc.tulip[eigenVectors+iii].symmetry  == type && (! c->i.type|| c->i.type == type)){
                
                printf("%dState%d:%d:,%d ,%1.15f, %d, %d , %1.1f,%1.15f\n", f1->body,jjj, f1->sinc.N1,jjj,f1->sinc.tulip[eigenVectors+iii].value,bodies(f1,eigenVectors+iii),type, deg(f1, type),f1->sinc.tulip[eigenVectors+iii].value2);
                if ( (c->i.outputFlag) % 2 == 1){
                    
                    for ( cmpl = 0 ; cmpl < spins(f1, eigenVectors+iii) ; cmpl++)
                    {
                        sprintf(str,"%s.%lld.eigen-%lld.%lld_mac", c->name,jjj,type,cmpl);
#ifndef APPLE
                        
                        FILE * out = fopen ( str,"w" );
                        outputFormat(&c->i.c, out, eigenVectors+iii,cmpl  );
                        fclose(out);
#endif
                    }
                }
                jjj++;
                
            }
    }
    if ( (c->i.outputFlag) % 2 == 1)
        tEdges(c);

    return 0;
}


INT_TYPE exec (struct calculation *c ){
    INT_TYPE splitFlag = 0,cycleFlag = 0;;
    
#ifdef splitTag
    if ( c->i.body >= two)
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
    //fprintf(out , " lattice\t%1.3f\t%lld\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->i.body);
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
        fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->i.body);
#ifndef APPLE
        INT_TYPE lV ;
        if ( (c1->rt.printFlag/8)%2 ){
            lV = tLoadEigenWeightsWithConstraints ( c1, c1->mem.fileList,c1->mem.constraintFile);
        } else {
            lV = tLoadEigenWeights ( c1, c1->mem.fileList);
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
    }else {
        INT_TYPE cycle =0,iteration,iter=0,EV,type2;
        INT_TYPE space;
        enum division usz;
        if ( c1->i.sectors ){
            c1->i.heliumFlag = c->i.nTargets;
            c1->i.nStates = c1->i.heliumFlag;
            if ( splitFlag )
                c1->i.nStates *= tPaths(c1->i.body , c1->i.type);
            c1->i.Iterations = 1;
            c1->rt.maxEntropy = 1;
            iModel(c1);
            if (f1->body == one ){
                t1BodyConstruction ( c1, eigen,c1->i.l2);
            }
            else
            {
                
                assignCores(f1, 2);
                tNBodyConstruction ( c1, build,  eigen);

                INT_TYPE hits[nSAG],nG = tSize(f1->body),nP = tPerms(f1->body), v,rr,r2,info,r,rank=0,space,type,LN2[3],*n1 =vectorLen(f1, eigen);
                length(f1, eigen, LN2);
                double vo,va,vm;
                for ( type = 0 ; type <  nG ; type++){
                    
                    for ( r = 0 ; r < n1[0]; r++){

                        streams(f1,foundationStructure,0,0)[r+(type)*n1[0]] = 2e9;
                        streams(f1,foundationStructure,0,1)[r+(type)*n1[0]] = 2e9;
                        streams(f1,foundationStructure,0,2)[r+(type)*n1[0]] = 2e9;
                    }
                }
                for ( i = 0 ; i < nSAG ; i++)
                    hits[i] = 0;
#ifdef OMP
#pragma omp parallel for private (r,type,r2,rr,rank,space,va,vo) schedule(static,1)
#endif
                    
                    
                    
                for ( r = 0 ; r < (INT_TYPE)(pow(n1[0],c1->i.seekPower)); r++){
#ifdef OMP
                    rank = omp_get_thread_num();
#else
                    rank = 0;
#endif
                    
                        f1->sinc.tulip[diagonalVectorA].Current[rank]= 1;
                        for ( space = 0;  space < 3; space++)
                            cblas_dcopy(n1[space], streams(f1,eigen,0,space)+r*n1[space], 1, streams(f1, diagonalVectorA,rank,space), 1);
                        f1->sinc.tulip[permutationVector].Current[rank]= 0;
                        
                        tBuildIrr(rank, f1, -1, diagonalVectorA, rank, permutationVector, rank);
                        
                        for ( type = 0 ; type < nG ; type++){
                            if ( hits[type] < c1->i.side ){

                            zero(f1, diagonalVectorA, rank);
                            
                            for ( r2 = 0; r2 < CanonicalRank(f1, permutationVector, rank); r2++){
                                cblas_daxpy(n1[0], getter(f1->body, type+1, r2), streams(f1, permutationVector,rank,0)+r2*n1[0], 1, streams(f1,diagonalVectorA, rank, 0),1);
                                cblas_daxpy(n1[1], getter(f1->body, type+1, r2), streams(f1, permutationVector,rank,1)+r2*n1[1], 1, streams(f1,diagonalVectorA, rank, 1),1);
                                cblas_daxpy(n1[2], getter(f1->body, type+1, r2), streams(f1, permutationVector,rank,2)+r2*n1[2], 1, streams(f1,diagonalVectorA, rank, 2),1);
                                //zero?
                            }
                            
                            vo = inner(rank,f1, diagonalVectorA,rank);
                            if ( vo > 0.01){
                                tScaleOne(f1, diagonalVectorA,rank, 1./sqrt(vo));
                                
                                cblas_dcopy(n1[0], streams(f1, diagonalVectorA,rank,0), 1, streams(f1,build,0,0)+r*n1[0]+(type)*LN2[0], 1);
                                cblas_dcopy(n1[1], streams(f1, diagonalVectorA,rank,1), 1, streams(f1,build,0,1)+r*n1[1]+(type)*LN2[1], 1);
                                cblas_dcopy(n1[2], streams(f1, diagonalVectorA,rank,2), 1, streams(f1,build,0,2)+r*n1[2]+(type)*LN2[2], 1);
                                
                                
                                
                                for ( space = 0 ; space < SPACE ; space++){
                                
                                    cblas_dgemv(CblasColMajor, CblasNoTrans, n1[space], n1[space], 1.0, streams(f1,spam,0,space), n1[space], streams(f1, diagonalVectorA,rank,space), 1, 1., streams(f1, diagonalVector,rank,space), 1);
                                    streams(f1,foundationStructure,0,space)[r+(type)*n1[space]]  = cblas_ddot(n1[space], streams(f1, diagonalVectorA,rank,space), 1, streams(f1, diagonalVector,rank,space), 1);

                                }
                                hits[type]++;
                            }
                        }
                    }
                }
                for ( space = 0; space < SPACE ; space++)
                    for ( type = 0 ; type < nG ; type++){
                        rr = 0;
                        for ( r = 0 ; r < n1[space]; r++)
                            if (streams(f1,foundationStructure,0,space)[r+(type)*n1[space]] < 1e9 ){
                                if( r != rr ){
                                    cblas_dcopy(n1[space], streams(f1,build,0,space)+r*n1[space]+(type)*LN2[space], 1, streams(f1,build,0,space)+rr*n1[space]+(type)*LN2[space],1);
                                    streams(f1,foundationStructure,0,space)[rr+(type)*n1[space]] = streams(f1,foundationStructure,0,space)[r+(type)*n1[space]];
                                }
                                rr++;
                            }
                        for ( r = rr ; r < n1[space]; r++){
                            streams(f1,foundationStructure,0,space)[r+(type)*n1[space]] = 2e9;
                        }
                        
                    }

                
                for ( space = 0; space < SPACE ; space++){
                    va = 2e9;
                    for ( type = 0 ; type < nG ; type++)
                        for ( v = 0 ; v < n1[space] ; v++){
                            if ( va > streams(f1,foundationStructure,0,space)[v+(type)*n1[space]]&& 1e9 > streams(f1,foundationStructure,0,space)[v+(type)*n1[space]])
                                va = streams(f1,foundationStructure,0,space)[v+(type)*n1[space]];
                        }
                    for ( type = 0; type < nG ; type++)
                        for ( v = 0 ; v < n1[space] ; v++)
                            streams(f1,foundationStructure,0,space)[v+(type)*n1[space]] -= va;
                }
            }
            
            usz = eigenVectors+c1->i.nStates;
            if ( tCollect(f1,c1->i.type,usz,c1->i.qFloor ,c->i.seekPower) != c1->i.qFloor){
                printf("could not muster \n");
                exit(0);
            }
            EV =c1->i.qFloor;
            cycleFlag = 1;
        }else {

            {
                INT_TYPE lines = 0;
                char ch ;
                FILE * fp = fopen(c1->mem.fileList,"r");
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
            assignCores(f1, 2);

#ifndef APPLE
            if ( tLoadEigenWeights ( c1, c1->mem.fileList) != c1->i.qFloor ){
                printf("set helium %d \n", c1->i.qFloor);
                exit(0);
            }
#endif
            usz = eigenVectors+c1->i.nStates;
            EV = xConstructFoundation (c1 , usz, c1->i.qFloor, c1,   eigenVectors,   c1->i.nStates ,1);

        }
        RdsSize = EV;
        c1->i.sectors = 0;
        for ( totalIter = 0,cycle = 0, iter = 1; cycle <= c->i.cycles ;totalIter++)
        {
            
            
            
            printf("Iter %d \n Cycle %d\n Iter %d\n cycleFlag %d\n\n", totalIter, cycle ,iter, cycleFlag);
            printf("EIGENVALUE %d @@ %d , %d \nUSZ %d @@ %d , %d \n", eigenVectors,c1->i.heliumFlag,c1->i.nStates, usz ,EV,RdsSize );

            fflush(stdout);
            
            

             {
                if ( tEigenCycle(&c1->i.c,Ha,'T', c1->i.heliumFlag, usz,RdsSize,0, EV,c1->i.type,  (3+!splitFlag),eigenVectors,twoBodyRitz))
                {
                    printf("usz:%d--ld %d %d\n", usz, RdsSize, EV);
                    
                    if (1 ){
                        printf("crap!");
                        exit(0);
                    }
                }
                 tFilter(&c1->i.c, c1->i.nStates, !(!c1->i.filter )* c->i.type, eigenVectors);
                 
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
#ifndef splitTag
              //  print(c1);
#endif
                
                {
                    if ( splitFlag )
                    {
                        INT_TYPE nP = tPerms(f1->body);
                        struct calculation *c3 = malloc(sizeof(struct calculation ));
                        *c3 = *c1;
                        c3->i.Iterations = 1;
                        c3->i.cycles = 0;
                        c3->i.nStates = c1->i.nStates;
                        c3->i.Iterations = 0;
                        c3->i.qFloor = 7;
                        c3->i.c.mem1->bootedMemory =0;
                        c3->i.c.sinc.tulip = NULL;
                        for ( space = 0; space <= SPACE ; space++)
                            c3->i.c.sinc.rose[space].stream = NULL;
                        
                        iModel(c3);
                        
                        INT_TYPE sp,g,i,ii,iii = 0,rank=0,sup;
                        
                        
                        for( g = 0; g < nSAG*nSAG*nSAG; g++)
                            c2->i.cSA[g] = 0;
                        
                        
                        for ( ii = 0; ii < c1->i.nStates ; ii++)
                        {
                            sup = tSizeUp(rank,f1,c1->i.type, eigenVectors+ii);
                            if ( sup ){
                                for ( sp = 0; sp  < spins (f1, eigenVectors+ii);sp++)
                                    xEqua(&c3->i.c, copyVector, sp, f1, eigenVectors+ii, sp);
                                if (tSelect(&c3->i.c,  iii,c3->i.type+nP, eigenVectors, copyVector,1)){
                                    c3->i.c.sinc.tulip[eigenVectors+iii].path  = tPath(f1, eigenVectors+ii );
                                    c2->i.cSA[tPath(f1, eigenVectors+ii ) ] += 1;
                                    iii++;
                                }
                            }
                        }
                        
                        
                        {
                            printf("begin count\n");
                            INT_TYPE nG = tSize(f1->body);
                            for( g = 0; g < nSAG*nSAG*nSAG; g++)
                                if ( c2->i.cSA[g] )
                                    printf("%d : 00 : %d:: (%d %d %d) --> %d\n",g, c2->i.cSA[g],c2->i.cSA[g]%nG+1 , (c2->i.cSA[g]/nG)%nG +1, (c2->i.cSA[g]/(nG*nG) % nG) +1, c->i.type);
                            printf("end count\n");

                        }
                        c2->i.qFloor = iii;
                        
                        c2->i.heliumFlag = c->i.nTargets;//number of states
                        c2->i.nStates = tPaths(c->i.body , c->i.type)*c2->i.heliumFlag;//number of paths ...
                        EV = c2->i.qFloor;
                        
                        iModel(c2);
                        usz = eigenVectors+c2->i.nStates;
                        {
                            for ( i = 0; i < c2->i.qFloor ; i++){
                                for ( sp = 0; sp < spins(f1, usz); sp++){
                                    xEqua(&c2->i.c, usz+i, sp, &c3->i.c, eigenVectors+i, sp);
                                }
                                printf("c3 %i --> c2 %d\n", eigenVectors+i,usz+i);
                            }
                        }
                        fModel(c3);
                        free(c3);
                        fModel(c1);
                        free(c1);
                        c1 = c2;
                        f1 = &c1->i.c;
                        
                        RdsSize = EV;
                    }else {
                        
                        {
                            nStatesTrans = 0;
                            nStatesFound = 0;
                            
                            INT_TYPE ii = 0,iii,type;
                            for ( type = 1 ; type <= 24 ; type++){
                                ii= 0;
                                for ( iii = 0; iii < c1->i.nStates ; iii++)
                                {
                                    if (f1->sinc.tulip[eigenVectors+iii].symmetry == type && (! c1->i.type || type == c1->i.type) ){
                                        ii++;
                                        if ( ii <= c->i.nTargets){
                                            nStatesTrans = iii+1;
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
                            c2->i.qFloor = iv;
                            EV = iv;
                            RdsSize = EV;
                        }
                        
                        iModel(c2);
                        xConstructFoundation (c2 , eigenVectors, c2->i.nStates, c1,eigenVectors,nStatesTrans ,1);
                        fModel(c1);
                        free(c1);
                        c1 = c2;
                        f1 = &c1->i.c;
                        usz = eigenVectors+c1->i.nStates;
                        
#ifdef OMP
#pragma omp parallel for private (ii) schedule(dynamic,1)
#endif
                        for ( ii = 0; ii < c1->i.nStates ; ii++){
                            tEqua(f1, usz+ii, 0, eigenVectors+ii , 0);
                            tEqua(f1, usz+ii, 1, eigenVectors+ii , 1);
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
                }
                
                
                if   ( tEigenCycle(&c1->i.c,Ha,'T', c1->i.heliumFlag, usz,RdsSize,0, EV, 0,0,eigenVectors,twoBodyRitz))
                {
                    RdsSize -= EV;
                    usz += EV;
                    if ( RdsSize == 0 )
                    {
                        printf("not good!\n");
                        exit(0);
                    }
                    cycleFlag = 1;
                    continue;//break from cycling/iterating...if failed LD test.
                }

                
            }else {
                cycleFlag = 1;//not likely to happend unless errro.
                printf("woops! ... looks like you accidently attempted an extra iteration\n");
            }
            fprintf(out,"EV %d ; type %d ; offset %d\n", EV, c1->i.type,RdsSize-EV);
            fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->i.body);
            fflush(stdout);
            
            INT_TYPE ii = 0,iii;
            for ( iii = 0; iii < c1->i.heliumFlag ; iii++)
            {
                if (1){
                    fprintf(out,"Condition%d:%d:,%lld ,%1.15f\n", ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value);
                    ii++;
                }
            }
            fflush(stdout);
            iter++;

        }
        

    }
    if ( splitFlag )
    {//exit
        INT_TYPE nP = tPerms(f1->body);
        enum division usz;
        INT_TYPE EV,space,sp,g,i,ii,iii = 0,rank=0,iv,sup;
        
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

        
        
        
        
        struct calculation *c3 = malloc(sizeof(struct calculation ));
        *c3 = *c1;
        c3->i.Iterations = 1;
        c3->i.cycles = 0;
        c3->i.nStates = c1->i.nStates;
        if ( c->i.type == 3 )
            c3->i.nStates *= 4;
        if ( c->i.type > 3 )
            c3->i.nStates *= 9;
        c3->i.Iterations = 0;
        c3->i.qFloor = 0;
        c3->i.c.mem1->bootedMemory =0;
        c3->i.c.sinc.tulip = NULL;
        for ( space = 0; space <= SPACE ; space++)
            c3->i.c.sinc.rose[space].stream = NULL;
        
        iModel(c3);
        
        
        
        for( g = 0; g < nSAG*nSAG*nSAG; g++)
            c2->i.cSA[g] = 0;
        
        
        for ( ii = 0; ii < c1->i.nStates ; ii++)
        {
            sup = tSizeUp(rank,f1,c1->i.type, eigenVectors+ii);
            if ( sup ){
                for ( sp = 0; sp  < spins (f1, eigenVectors+ii);sp++)
                    xEqua(&c3->i.c, copyVector, sp, f1, eigenVectors+ii, sp);
                if (tSelect(&c3->i.c,  iii,c3->i.type+nP, eigenVectors, copyVector,1)){
                    c3->i.c.sinc.tulip[eigenVectors+iii].path  = tPath(f1, eigenVectors+ii );
                    c2->i.cSA[tPath(f1, eigenVectors+ii ) ] += 1;
                    iii++;
                }
            }
        }
        fModel(c1);
        free(c1);
        
        {
            INT_TYPE nG = tSize(f1->body);
            for( g = 0; g < nSAG*nSAG*nSAG; g++)
                if ( c2->i.cSA[g] )
                    printf("%d : final : %d:: (%d %d %d) --> %d\n",g, c2->i.cSA[g],c2->i.cSA[g]%nG+1 , (c2->i.cSA[g]/nG)%nG +1, (c2->i.cSA[g]/(nG*nG) % nG) +1, c->i.type);
            
        }
        c2->i.qFloor = iii;
        
        c2->i.heliumFlag = c->i.nTargets;//number of states
        c2->i.nStates = c2->i.heliumFlag;//number of paths ...
        EV = c2->i.qFloor;
        
        iModel(c2);
        usz = eigenVectors+c2->i.nStates;
        {
            for ( i = 0; i < c2->i.qFloor ; i++)
                for ( sp = 0; sp < spins(f1, usz); sp++)
                    xEqua(&c2->i.c, usz+i, sp, &c3->i.c, eigenVectors+i, sp);
            
        }
        fModel(c3);
        free(c3);
        
        c1 = c2;
        f1 = &c1->i.c;

        RdsSize = EV;
        {
            if ( tEigenCycle(&c1->i.c,Ha,'T', c1->i.heliumFlag, eigenVectors + c1->i.nStates,RdsSize,0, EV,c1->i.type, 4,eigenVectors,twoBodyRitz))
            {
                printf("ld %d %d\n", RdsSize, EV);
                
                if (1 ){
                    printf("crap!");
                    exit(0);
                }
            }
            tFilter(&c1->i.c, c1->i.nStates, 0, eigenVectors);
            
        }
    }//exit
    
    *c = *c1;
    c->i.c.sinc.tulip = c1->i.c.sinc.tulip;
    c->i.c.sinc.rose[0].stream =c1->i.c.sinc.rose[0].stream;
    c->i.c.sinc.rose[1].stream =c1->i.c.sinc.rose[1].stream;
    c->i.c.sinc.rose[2].stream =c1->i.c.sinc.rose[2].stream;
    c->i.c.sinc.rose[3].stream =c1->i.c.sinc.rose[3].stream;
    c1->mem.bootedMemory = 0;
    c1->i.c.sinc.tulip = NULL;
    c1->i.c.sinc.rose[0].stream = NULL;
    c1->i.c.sinc.rose[1].stream = NULL;
    c1->i.c.sinc.rose[2].stream = NULL;
    c1->i.c.sinc.rose[3].stream = NULL;
    fModel(c1);
    free(c1);
    return 0;
}
#ifndef APPLE
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
    
    func1.fn = Coulomb;//Yukawa, Erf
    func1.param[0] = 1.;//scalar out front
    func1.param[1] = 1.;//not used by elemCal...a quadrature parameter
    func1.param[2] = 1.;//parameter in Erf and Yukawa (scale/mass)
    func1.param[3] = 1.;//second a parameter, unused.
    getDescription ( &func1 ,0.,stdout);

    INT_TYPE periodic = 0;//not coded for peridicity/.
    INT_TYPE i ;
    
    
    if ( argc != 6 ){
        printf("L,BODY, Distance, Gaussian, accelerate\n");
        exit(0);
        
    }
    INT_TYPE l = atoi(argv[1]);
    INT_TYPE body = atoi(argv[2]);
    double x = atof(argv[3]);
    double b = atof(argv[4]);
    INT_TYPE acc = atoi(argv[5]);

    if ( acc ){
        if ( l > 1 ){
            printf("warning,  have not completed acceleration of L > 1 yet...");
            
        }
        if ( body == 1 ){
            printf("warning,  acceleration is focused on 2-body interactions");
        }
    }
    
    
    printf("ANGULAR %d \n %d-Body interaction", l , body);
    struct general_2index g3[3];
    {
        i = 0;
        g3[i].gaussianAccelerationFlag = acc;
        
        g3[i].i[0].b0 = b;//body 1 Matrix
        g3[i].i[0].x0 = x;
        g3[i].i[0].l0 = l;
        
        g3[i].i[0].b1 = b;
        g3[i].i[0].x1 = x;
        g3[i].i[0].l1 = l;
        
        g3[i].i[1].action = 0;
        g3[i].i[1].b0 = b;//body 2 matrix
        g3[i].i[1].x0 = 0;
        g3[i].i[1].l0 = l;
        
        g3[i].i[1].b1 = b;
        g3[i].i[1].x1 = 0;
        g3[i].i[1].l1 = l;

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
        i = 1;

        g3[i].gaussianAccelerationFlag = acc;

        g3[i].i[0].b0 = b;//body 1 Matrix
        g3[i].i[0].x0 = 0;
        g3[i].i[0].l0 = l;
        
        g3[i].i[0].b1 = b;
        g3[i].i[0].x1 = 0;
        g3[i].i[0].l1 = l;
        
        g3[i].i[1].action = 0;
        g3[i].i[1].b0 = b;//body 2 matrix
        g3[i].i[1].x0 = 0;
        g3[i].i[1].l0 = l;
        
        g3[i].i[1].b1 = b;
        g3[i].i[1].x1 = 0;
        g3[i].i[1].l1 = l;

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

        g3[i].gaussianAccelerationFlag = acc;

        g3[i].i[0].b0 = b;//body 1 Matrix
        g3[i].i[0].x0 = 0;
        g3[i].i[0].l0 = l;
        
        g3[i].i[0].b1 = b;
        g3[i].i[0].x1 = 0;
        g3[i].i[0].l1 = l;
        
        g3[i].i[1].action = 0;
        g3[i].i[1].b0 = b;//body 2 matrix
        g3[i].i[1].x0 = 0;
        g3[i].i[1].l0 = l;
        
        g3[i].i[1].b1 = b;
        g3[i].i[1].x1 = 0;
        g3[i].i[1].l1 = l;

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
//    g3[0].i[0].l0 = 2;
        printf(" %1.15f\n", elementCal(1e-3,-1, g3));

#endif

}

#endif

INT_TYPE buildElectronFreeInteraction ( struct calculation * c1, enum division mat){
    double value;
    struct field * f1 = &c1->i.c;
    INT_TYPE g,space,r,i,j,n,m,N1 = f1->sinc.N1, N2 = f1->sinc.N1*f1->sinc.N1;

    double coef = pow(f1->Ne,0.33333333333333);
    if ( ! CanonicalRank(f1, interactionExchange, 0))
        return 0;
    
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
    c0.i.body = 1;
    c0.i.iCharge = 1;
    c0.i.type = 1;
    c0.i.cycles = 0;
    c0.i.Iterations = 0;
    c0.i.nStates = 1;
    c0.i.heliumFlag = 1;
    c0.i.qFloor = f1->sinc.N1;
    c0.i.sectors = 1;
    c0.i.bRank = 1;//assume its a separable solution
    exec(&c0);
    for ( r = 0; r < CanonicalRank(f1, interactionExchange, 0); r++){
        
        for ( i = 0; i < N1 ; i++)
            for ( j = 0; j < N1; j++)
                for (space = 0; space < SPACE ;space++)
                    
                {
                    value = 0.;
                    for ( n = 0 ; n < N1 ; n++)
                        for ( m = 0 ; m < N1 ; m++)
                            value += streams(f1, interactionExchange, 0,space)[(N2*(N1*i+n)+(N1*j+m))+r*N2*N2]* streams(&c0.i.c,eigenVectors,0,space)[n]*streams(&c0.i.c,eigenVectors,0,space)[m];
                    streams(f1, mat,0,space)[N1*i+j+r*N2] = -coef*value;
                }
    }
    fModel(&c0);
    
    f1->sinc.tulip[mat].Current[0] = CanonicalRank(f1, interactionExchange, 0);
    
    if ( bodies(f1, eigenVectors)==two){
        tClear(f1, copy);
        tId(f1, copy,0);
        tScale(f1, copy,-3.*traceOne(f1, mat, 0)/(N1*N1*N1)/sqr(f1->Ne));
        //
        tAddTw(f1, mat, 0,copy,0);
    }
    else     if ( bodies(f1, eigenVectors)==three){
        tClear(f1, copy);
        tId(f1, copy,0);
        tScale(f1, copy,-6.*traceOne(f1, mat, 0)/(N1*N1*N1)/sqr(f1->Ne));
        //
        tAddTw(f1, mat,0, copy,0);
    }
    else     if ( bodies(f1, eigenVectors)==four){
        tClear(f1, copy);
        tId(f1, copy,0);
        tScale(f1, copy,-10.*traceOne(f1, mat, 0)/(N1*N1*N1)/sqr(f1->Ne));
        //
        tAddTw(f1, mat, 0,copy,0);
    }
    
    
    f1->sinc.tulip[mat].stop[0][0] = CanonicalRank(f1, mat,0);
    
    return 0;
}
