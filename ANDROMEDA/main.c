/*
 *  main.c
 *
 *
 *  Copyright 2018 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University
 *  and the Robert A. Welch Foundation.
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
    if ( (c->i.outputFlag) % 2 == 1){
        char str[MAXSTRING];
        INT_TYPE iii,jjj=1,cmpl;
        for ( type = 1 ; type <= 24 ; type++){
            jjj = 1;
            for ( iii = 0; iii < c->i.nStates  ; iii++)
                if ( f1->sinc.tulip[eigenVectors+iii].symmetry  == type && (! c->i.type|| c->i.type == type)){
                    
                    printf("%dState%d:%d:,%d ,%1.15f, %d, %d , %1.1f,%1.15f\n", f1->body,jjj, f1->sinc.N1,jjj,f1->sinc.tulip[eigenVectors+iii].value,bodies(f1,eigenVectors+iii),type, deg(f1, type),f1->sinc.tulip[eigenVectors+iii].value2);
                    
                    for ( cmpl = 0 ; cmpl < 2 ; cmpl++)
                    {
                        sprintf(str,"%s.%lld.eigen-%lld.%lld_mac", c->name,jjj,type,cmpl);
#ifndef APPLE

                        FILE * out = fopen ( str,"w" );
                        outputFormat(&c->i.c, out, eigenVectors+iii,cmpl  );
#endif
                    }
                    jjj++;
                }
        }
    }
    return 0;
}

> data <- transform ( data , WEIGHT = ifelse ( data$BODY == 4 , (ifelse ( data$CLASS == 4 , 1. ,0) + ifelse( data$CLASS == 3 , 0.5, 0 ) + ifelse( data$CLASS == 2 , 5., 0 ) ),0)+ ifelse ( data$BODY == 3 , (ifelse ( data$CLASS == 2 , 4. ,0) + ifelse( data$CLASS == 3 , 1., 0 ) ),0)+ ifelse ( data$BODY == 2 , (ifelse ( data$CLASS == 1 , 1. ,0) + ifelse( data$CLASS == 2 , 3., 0 ) ),0))
INT_TYPE exec (struct calculation *c ){
    
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
        INT_TYPE cycle =0,iteration,EV;
        INT_TYPE space;
        enum division usz;
        if ( c1->i.sectors ){
            c1->i.heliumFlag = c1->i.qFloor;
            c1->i.Iterations = 1;
            c1->i.iRank = 1;
            c1->rt.maxEntropy = 1;
            iModel(c1);
   
tNBodyConstruction ( c1, build,  eigen);

            usz = eigenVectors+c1->i.nStates;
            if ( tCollect(f1,0,usz,c1->i.qFloor) != c1->i.qFloor ){
                printf("could not muster \n");
                exit(0);
            }
            EV =c1->i.qFloor;
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
            iModel(c1);

#ifndef APPLE
            if ( tLoadEigenWeights ( c1, c1->mem.fileList) != c1->i.qFloor ){
                printf("set helium\n");
                exit(0);
            }
#endif
            usz = eigenVectors+c1->i.nStates;
            EV = xConstructFoundation (c1 , usz, c1->i.qFloor, c1,   eigenVectors,   c1->i.nStates ,1);
        }
        RdsSize = EV;

        c1->i.sectors = 0;
        c1->i.bRank = c->i.bRank;
        for ( totalIter = 0,cycle = 0; cycle < c->i.cycles ; )
        {
            if ( tEigenLoad(&c1->i.c,Ha,'T', c1->i.nStates, eigenVectors + c1->i.nStates,RdsSize,EV, 2,eigenVectors,twoBodyRitz))
            {
                printf("ld\n");
                exit(0);
            }
            {
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
                
            }
            
            
            
            
            
            struct calculation *c2 = malloc(sizeof(struct calculation ));
            *c2 = *c1;
            c2->i.Iterations = c->i.Iterations;
            c2->rt.maxEntropy = c->rt.maxEntropy;
            c2->i.heliumFlag = nStatesFound;
            c2->i.iRank = c1->i.bRank;
            if (! nStatesTrans ){
                printf("bail!\n");
                exit(0);
            }

            
            c2->i.c.mem1->bootedMemory =0;
            c2->i.c.sinc.tulip = NULL;
            for ( space = 0; space <= SPACE ; space++)
                c2->i.c.sinc.rose[space].stream = NULL;
            
            
            if ( totalIter  >= cycle *(c->i.Iterations-1) ){
                cycle++;
                c2->i.bRank++;
                INT_TYPE step = 1;
                c2->i.epi  += step;
                c2->i.d *= pow( (2.* c1->i.epi + 1.) /(2.*c2->i.epi + 1. ),c2->i.attack);
            
                c2->i.vectorMomentum *= (c1->i.d* (2.* c1->i.epi + 1.) )/(c1->i.d*(2.*c2->i.epi + 1. ));
                printf("attack %f %f -> %f\n",c2->i.attack,c1->i.d, c2->i.d);
                
                {
                    fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->i.body);
                    
                    INT_TYPE ii = 0,iii,type;
                    for ( type = 1 ; type <= 24 ; type++){
                        ii= 0;
                        for ( iii = 0; iii < c1->i.nStates ; iii++)
                        {
                            if (f1->sinc.tulip[eigenVectors+iii].symmetry == type && (! c1->i.type || type == c1->i.type) ){
                                printf("%dState%d:%d:,%d ,%1.15f, %d, %d , %1.1f,%1.15f\n", f1->body,ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value,bodies(f1,eigenVectors+iii),type, deg(f1, type),f1->sinc.tulip[eigenVectors+iii].value2);
                                ii++;
                            }
                        }}
                    
                    tEdges(c1);
                    fflush(stdout);
                }
            }
                {
                    INT_TYPE ii = 0,iii,type;
                    ii= 0;
                    for ( iii = 0; iii < c1->i.nStates ; iii++)
                    {
                        {
                            printf("%dREF%d:%d:,%d ,%1.15f, %d, %d , %1.1f,%1.15f\n", f1->body,ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value,bodies(f1,eigenVectors+iii),f1->sinc.tulip[eigenVectors+iii].symmetry, deg(f1, f1->sinc.tulip[eigenVectors+iii].symmetry),f1->sinc.tulip[eigenVectors+iii].value2);
                            ii++;
                        }
                    }
                }
                fflush(stdout);
                
            
            
            {
                INT_TYPE ii;
                {
                    INT_TYPE ii,iii,iv=0;
                    for ( iii = 0; iii < nStatesTrans ; iii++){
                        ii = xConstructFoundation (c2 , 0,0, c1,eigenVectors+iii,1 ,1);
                        iv += ii;
                    }
                    
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
                
            }

            
            while (totalIter  < cycle *(c->i.Iterations-1)){
                totalIter++;

                if ( tEigenLoad(&c1->i.c,Ha,'T', c1->i.nStates, usz,RdsSize,EV, 0,eigenVectors,twoBodyRitz))
                {
                    printf("ld\n");
                    
                    RdsSize -= EV;
                    
                    if (RdsSize == 0 ){
                        printf("ran out!\n");
                        exit(0);
                    }
                    break;
                }
 
                
                if (  0 ) {
                    plusSize = tMinorDivideIteration(f1,Ha, 1,0,usz+RdsSize-EV,EV,2*EV,c1->i.group);
                } else {
                    plusSize = tGreatDivideIteration(f1,Ha, 1,0,usz+RdsSize-EV,EV,2*EV,c1->i.group);
                }

                if ( ! plusSize )
                    break;
                
                RdsSize += plusSize - EV;

                tFilter(f1,EV, c1->i.type,usz+RdsSize-EV);

                fprintf(out,"EV %d ; type %d ; offset %d\n", EV, c1->i.type,RdsSize-EV);
                fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->i.body);

                INT_TYPE ii = 0,iii;
                for ( iii = 0; iii < imin(RdsSize,c1->i.nStates) ; iii++)
                {
                    if (1){
                        fprintf(out,"Condition%d:%d:,%lld ,%1.15f\n", ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value);
                        ii++;
                    }
                }
                fflush(stdout);
            }

        }
        
        if ( tEigenLoad(&c1->i.c,Ha,'T', c1->i.nStates, eigenVectors + c1->i.nStates,RdsSize,EV, 2,eigenVectors,twoBodyRitz))
        {
            printf("ld\n");
            exit(0);
        }

    }
    
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
    struct calculation c = bootShell(argc, argv);
    exec(&c);
    print(&c);
    fModel(&c);
}
#endif



INT_TYPE buildElectronFreeInteraction ( struct calculation * c1, enum division mat){
    double value;
    struct field * f1 = &c1->i.c;
    INT_TYPE space,r,i,j,n,m,N1 = f1->sinc.N1, N2 = f1->sinc.N1*f1->sinc.N1;

    double coef = pow(f1->Ne,0.33333333333333);
    if ( ! CanonicalRank(f1, interactionExchange, 0))
        return 0;
    
    struct calculation c0 = *c1;
    
    
    c0.mem.bootedMemory = 0;
    c0.i.c.sinc.tulip = NULL;
    c0.i.c.sinc.rose[0].stream = NULL;
    c0.i.c.sinc.rose[1].stream = NULL;
    c0.i.c.sinc.rose[2].stream = NULL;
    c0.i.c.sinc.rose[3].stream = NULL;

    c0.rt.printFlag = 0;
    c0.i.body = 1;
    c0.i.iCharge = 1;
    c0.i.type = 0;
    c0.i.cycles = 0;
    c0.i.Iterations = 1;
    c0.i.nStates = 1;
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
