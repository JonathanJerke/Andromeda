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
    finalizeInit(&c1);
#else
    struct calculation c1 =initCal();
#endif

    return c1;
}

INT_TYPE exec (struct calculation *c ){
    
    struct calculation *c1 = malloc(sizeof(struct calculation));
    *c1 = *c;
    
    INT_TYPE a;
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
    f1->twoBody.num = c1->i.canonRank;
    f1->oneBody.num = c1->i.canonRank;
    
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
        INT_TYPE cycle ,iteration,EV;
        INT_TYPE RdsSize = 0,space;
        enum division usz;
        if ( c1->i.sectors ){
            iModel(c1);
            tNBodyConstruction ( c1, build,  eigen);
            usz = eigenVectors+c1->i.nStates;
            if ( tCollect(f1,c1->i.type,usz,c1->i.qFloor) != c1->i.qFloor ){
                printf("could not muster \n");
                exit(0);
            }
            EV =c1->i.qFloor;
        }else {
            c1->i.qFloor = c1->i.nStates;
            iModel(c1);

#ifndef APPLE
            if ( tLoadEigenWeights ( c1, c1->mem.fileList) != c1->i.qFloor ){
                printf("set helium\n");
                exit(0);
            }
#endif
            usz = eigenVectors+c1->i.nStates;
            EV = xConstructFoundation (c1 , usz, c1->i.qFloor, c1,   eigenVectors,   c1->i.nStates ,c1->i.iRank);
        }
        RdsSize = EV;

        for ( iteration = 1 ; iteration < c1->i.Iterations; iteration++){
            if ( ! RdsSize ){
                printf("lost my marbles!\n");
                exit(0);
            }
            RdsSize += tGreatDivideIteration(f1,Ha, 1,0,usz+RdsSize-EV,EV,2*EV,0.0)-EV;
            tFilter(f1,EV, c1->i.type,usz+RdsSize-EV);
            while (tEigenLoad(&c1->i.c,Ha,c1->i.type, c1->i.nStates, usz,RdsSize,EV, 0,eigenVectors,twoBodyRitz)){
                INT_TYPE i,cmpl,mini = 999;
                for ( i = 0 ; i < EV ; i++){
                    for ( cmpl = 0; cmpl < 2 ;cmpl++)
                        if(f1->sinc.tulip[usz+RdsSize-EV+i].Current[cmpl]){
                            f1->sinc.tulip[usz+RdsSize-EV+i].Current[cmpl]--;
                        }
                    mini = imin(mini,CanonicalRank(f1, usz+RdsSize-EV+i, 0)+CanonicalRank(f1, usz+RdsSize-EV+i, 1));
                }
                if ( ! mini ){
                    iteration = c1->i.Iterations;
                    RdsSize -= EV;
                    break;
                }
            }
            {
                fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->i.body);
                
                INT_TYPE ii = 0,iii;
                for ( iii = 0; iii < imin(RdsSize,c1->i.nStates) ; iii++)
                {
                    if (1 ){
                        printf("Condition%d:%d:,%d ,%1.15f\n", ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value);
                        ii++;
                    }
                }
                fflush(stdout);
            }
        }
        tEigenLoad(&c1->i.c,Ha,c1->i.type, c1->i.nStates, usz,RdsSize,EV, 2,eigenVectors,twoBodyRitz);
        {
            fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->i.body);

            INT_TYPE ii = 0,iii,type;
            for ( type = 1 ; type <= 24 ; type++){
                ii= 0;
            for ( iii = 0; iii < imin(RdsSize,c1->i.nStates) ; iii++)
            {
                if (f1->sinc.tulip[eigenVectors+iii].symmetry == type && (! c1->i.type || type == c1->i.type) ){
                    printf("State%d:%d:,%d ,%1.15f, %d, %d , %1.1f\n", ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value,bodies(f1,eigenVectors+iii),type, deg(f1, type));
                    ii++;
                }
            }}
            
            tEdges(c1);
            fflush(stdout);
        }
        c1->i.sectors = 0;
        for ( cycle = 0; cycle < c1->i.cycles ; cycle++)
        {

            struct calculation *c2 = malloc(sizeof(struct calculation ));
            *c2 = *c1;
            c2->i.c.mem1->bootedMemory =0;
            c2->i.c.sinc.tulip = NULL;
            for ( space = 0; space <= SPACE ; space++)
                c2->i.c.sinc.rose[space].stream = NULL;
            {
                INT_TYPE step = 1;
                c2->i.epi  += step;
                c2->i.d *= pow( (2.* c1->i.epi + 1.) /(2.*c2->i.epi + 1. ),c2->i.attack);
                printf("attack %f %f -> %f\n",c2->i.attack,c1->i.d, c2->i.d);
                
                c2->i.vectorMomentum *= (c1->i.d* (2.* c1->i.epi + 1.) )/(c1->i.d*(2.*c2->i.epi + 1. ));

            }
            c2->i.iRank = c1->i.bRank;//ADDED
            c2->i.bRank += 1;
            {
                c2->i.qFloor = 0;
                   for ( i = 0 ; i < imin(RdsSize,c1->i.nStates) ; i++)
                       if ( f1->sinc.tulip[eigenVectors+i].symmetry == c1->i.type|| ! c1->i.type){
                           c2->i.qFloor++;
                       }
            }
            fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c2->i.d,2*c2->i.epi+1, c2->i.body);
            fflush(out);

            iModel(c2);
            enum division usz = eigenVectors+c2->i.nStates;
            EV = xConstructFoundation (c2 , usz, c2->i.qFloor, c1,   eigenVectors,   c1->i.nStates ,c2->i.iRank);
            fModel(c1);
            free(c1);
            c1 = c2;
            f1 = &c1->i.c;
            
            
            RdsSize= EV;
            for ( iteration = 1 ; iteration < c1->i.Iterations; iteration++){
                RdsSize += tGreatDivideIteration(f1,Ha, 1,0,usz+RdsSize-EV,EV,2*EV,0.0)-EV;
                tFilter(f1,EV, c1->i.type,usz+RdsSize-EV);
                
                while (tEigenLoad(&c1->i.c,Ha,c1->i.type, c1->i.nStates, usz,RdsSize,EV, 0,eigenVectors,twoBodyRitz)){
                    INT_TYPE i,cmpl,mini = 999;
                    for ( i = 0 ; i < EV ; i++){
                        for ( cmpl = 0; cmpl < 2 ;cmpl++)
                            if(f1->sinc.tulip[usz+RdsSize-EV+i].Current[cmpl]){
                                f1->sinc.tulip[usz+RdsSize-EV+i].Current[cmpl]--;
                            }
                        mini = imin(mini,CanonicalRank(f1, usz+RdsSize-EV+i, 0)+CanonicalRank(f1, usz+RdsSize-EV+i, 1));
                    }
                    if ( ! mini ){
                        iteration = c1->i.Iterations;
                        RdsSize -= EV;
                        break;
                    }
                }
                {
                    fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->i.body);
                    
                    INT_TYPE ii = 0,iii;
                    for ( iii = 0; iii < imin(RdsSize,c1->i.nStates) ; iii++)
                    {
                        if (1 ){
                            printf("Condition%d:%d:,%lld ,%1.15f\n", ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value);
                            ii++;
                        }
                    }
                    fflush(stdout);
                }
            }
            tEigenLoad(&c1->i.c,Ha,c1->i.type, c1->i.nStates, usz,RdsSize,EV, 2,eigenVectors,twoBodyRitz);
            {
                fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->i.body);
                
                INT_TYPE ii = 0,iii,type;
                for ( type = 1 ; type <= 24 ; type++){
                    ii= 0;
                    for ( iii = 0; iii < imin(RdsSize,c1->i.nStates) ; iii++)
                    {
                        if (f1->sinc.tulip[eigenVectors+iii].symmetry == type && (! c1->i.type || type == c1->i.type) ){
                            printf("State%d:%d:,%d ,%1.15f, %d, %d , %1.1f\n", ii+1, f1->sinc.N1,ii+1,f1->sinc.tulip[eigenVectors+iii].value,bodies(f1,eigenVectors+iii),type, deg(f1, type));
                            ii++;
                        }
                    }}
                
                tEdges(c1);
                fflush(stdout);
            }

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
    INT_TYPE type;
    struct field * f1 = & c.i.c;
    if ( (c.i.outputFlag) % 2 == 1){
        char str[MAXSTRING];
        INT_TYPE iii,jjj=1,cmpl;
        for ( type = 1 ; type <= 24 ; type++){
            jjj = 1;
            for ( iii = 0; iii < c.i.nStates  ; iii++)
                if ( f1->sinc.tulip[eigenVectors+iii].symmetry  == type && (! c.i.type|| c.i.type == type)){
                    for ( cmpl = 0 ; cmpl < 2 ; cmpl++)
                    {
                        sprintf(str,"%s.%lld.eigen-%lld.%lld_mac", c.name,jjj,type,cmpl);
                        FILE * out = fopen ( str,"w" );
                        outputFormat(&c.i.c, out, eigenVectors+iii,cmpl  );
                    }
                    jjj++;
                }
        }
    }
    fModel(&c);
}
#endif
