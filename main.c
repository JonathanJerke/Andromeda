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




INT_TYPE exec (struct calculation *c ){
    INT_TYPE splitFlag = 0,cycleFlag = 0;;
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
    fprintf( out,   "\t| %d Body \n\t| %d Components\t\t\n\t| %d %dD Electrons \t\n\t| %d Atoms \t\t\t\n",c1->rt.body,SPACE,c1->i.c.Ne,COMPONENT ,f1->Na);
    if ( c1->i.Angstroms )
    {
        printf("\t| Input Angstroms\n");

    }else{
        printf("\t| Input Units Bohrs\n\n\n");
    }
    printf ("*Geometry\n");
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
    printf(".Geometry\n\n");

    INT_TYPE i;
    INT_TYPE nSlot ;
#ifdef OMP
    if ( c1->i.omp > MaxCore ){
        printf("lanes > MaxCore\n");
        exit(0);
    }
    rt->NLanes = c1->i.omp;
#pragma omp parallel for private (i)
    for ( i = 0; i < MaxCore ; i++){
        nSlot = omp_get_num_threads();
    }
    rt->NSlot = nSlot;
    if ( rt->NLanes > rt->NSlot ){
        printf("lanes > available\n");
        exit(0);
    }

#ifdef MKL
   if ( rt->NSlot < c->i.mkl * c->i.omp )
   {
       printf("not enough slots for mkl*omp\n" );
       exit(0);
   }
    else
    {
        printf("mkl %d\n", c1->i.mkl);
    }
#endif
    printf("lanes %lld\n", rt->NLanes);
    printf("slots %lld\n", rt->NSlot);
    fflush(stdout);
#else
    nSlot = 1;
#ifndef APPLE
#endif
#endif

    
    if ( !(!c1->rt.printFlag)){
        c1->i.nStates = c->i.heliumFlag;
        iModel(c1);
        fprintf(out , " lattice\t%1.3f\t%d\t-Box\t body \t%d\n", c1->i.d,2*c1->i.epi+1, c1->rt.body);
#ifndef APPLE
        INT_TYPE lV ;
        if ( (c1->rt.printFlag/8)%2 ){
//lV = tLoadEigenWeightsWithConstraints ( c1, c1->mem.fileList,c1->mem.constraintFile);
        } else {
            lV = tLoadEigenWeights ( c1, c1->mem.fileList[0], eigenVectors);
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
        
       // printOutput ( f1,lV);
       // printVectorOutput ( f1,lV);
#endif
    } else {
        INT_TYPE xyz,cycle =0,iteration,iter=0,EV,type2;
        INT_TYPE space;
        enum division usz;
//        if ( c1->i.sectors && ! c1->i.vectorOperatorFlag){
//            c1->i.heliumFlag = c->i.nTargets;
//            c1->i.nStates = c1->i.heliumFlag;
//            c1->i.Iterations = 1;
//            c1->rt.maxEntropy = 1;
//            iModel(c1);
//
//
//            if ( c1->rt.body == one  &&0  ){
//                tBoot1Construction(c1, build);
//                EV =   tCollect(f1,0,f1->sinc.user,c1->i.qFloor ,1);
//            }
//            else{
//                tBootManyConstruction(c1);
//                EV =   tCollect(f1,0,f1->sinc.user,c1->i.qFloor ,1);
//            }
//            cycleFlag = 1;
//        }else if (! c1->i.sectors && ! c1->i.vectorOperatorFlag  ){
//            {
//                INT_TYPE lines = 0;
//                char ch ;
//                FILE * fp = fopen(c1->mem.fileList[0],"r");
//                if ( fp == NULL ) {
//                    printf("file?\n");
//                    exit(0);
//                }
//                while(!feof(fp))
//                {
//                    ch = fgetc(fp);
//                    if(ch == '\n')
//                    {
//                        lines++;
//                    }
//                }
//                fclose(fp);
//                c->i.nTargets = lines;
//            }
//            c1->i.iRank = c->i.bRank;
//            c1->i.qFloor = c->i.nTargets;
//            c1->i.heliumFlag = c->i.nTargets;
//            c1->i.nStates = c->i.nTargets;
//
//            iModel(c1);
//
//#ifndef APPLE
//            if ( tLoadEigenWeights ( c1, c1->mem.fileList[0], eigenVectors) != c1->i.qFloor ){
//                printf("set helium %d \n", c1->i.qFloor);
//                exit(0);
//            }
//            tFilter(&c1->i.c, c1->i.nStates, 0, eigenVectors);
//#endif
//            EV = xConstructFoundation (c1 , f1->sinc.user, c1->i.qFloor, c1,   eigenVectors,   c1->i.nStates ,1);
//
//        }else if ( c1->i.sectors && c1->i.vectorOperatorFlag ){
//            {
//                c1->i.heliumFlag = c->i.nTargets;
//                c1->i.nStates = c1->i.heliumFlag;
//                c1->i.vectorOperatorFlag = countLinesFromFile(c1, 1);
//            }
//            iModel(c1);
//
//            if ( c1->i.bodyVectorOperator != f1->sinc.tulip[f1->sinc.vectorOperator].space[0].body  ){
//                printf("body count!\n");
//                exit(0);
//
//            }
//#ifndef APPLE
//            if ( c1->i.vectorOperatorFlag )
//                printf("vectorOperator Count %d\n",c1->i.vectorOperatorFlag );
//
//            if ( tLoadEigenWeights ( c1, c1->mem.operatorName ,f1->sinc.vectorOperator) != c1->i.vectorOperatorFlag){
//                printf("set helium %d \n", c1->i.vectorOperatorFlag);
//                exit(0);
//            }
//
//#endif
//            if ( c1->rt.body == one &&0 )
//                tBoot1Construction(c1, build);
//            else
//                tBootManyConstruction(c1);
//
//            EV =   tCollect(f1,0,f1->sinc.user,c1->i.qFloor ,1);
//            cycleFlag = 1;
//
//        }else {
//            printf("sectors?");
//            exit(0);
//        }
        RdsSize = EV;
        c1->i.sectors = 0;
        for ( totalIter = 0,cycle = 0, iter = 1; cycle <= c->i.cycles ;totalIter++)
        {
            
            fprintf(out,"\n\n\t\t| Foundation \t: %d\t\t\n\t\t| Krylov Space \t: %d\t\t \n\t\t| Irrep \t\t: %d\t\t\t\n", EV,RdsSize, c1->i.irrep);
            fprintf(out , "\t\t| Body \t\t\t: %d\t\t\t\n\t\t| Box \t\t\t: %d\t\t\n\t\t| Lattice \t\t: %1.3f\t\t\n\n",c1->rt.body,2*c1->i.epi+1,c1->i.d );
            fflush(stdout);

            
            printf("\n\ttotal Count\t|  %d\n\tCycle\t\t|  %d\n\tIteration \t|  %d\n\n", totalIter, cycle ,iter);
            printf("\tenum numbers\n\t eigenVectors = %d | allocated %d \n\t Users begin = %d | allocated %d\n\t\t\n", eigenVectors, c1->i.heliumFlag, f1->sinc.user , EV   );
            
            fflush(stdout);
            
            
             {
                 struct name_label u  = f1->sinc.tulip[Ha];
                 struct name_label k  = f1->sinc.tulip[kinetic];
                 struct name_label v  = f1->sinc.tulip[linear];

                if ( tEigenCycle(&c1->i.c,Ha,'T', c1->i.heliumFlag, f1->sinc.user,RdsSize,0, EV,c1->i.irrep,  4,eigenVectors,twoBodyRitz))
                {
                    printf("usz:%d--ld %d %d\n", f1->sinc.user, RdsSize, EV);
                    
                    if (1 ){
                        printf("crap!");
                        fflush(stdout);
                        exit(0);
                    }
                }

                 tFilter(&c1->i.c, c1->i.nStates, !(!c1->i.filter )* c->i.irrep, eigenVectors);
                 {
                     INT_TYPE iii ;
                     for ( iii = 0; iii < c1->i.heliumFlag ; iii++){
                         printf ( "\n Vector \t%d\n", iii+1);
                         printExpectationValues(f1, Ha, eigenVectors+iii);
                 //        tEdges(c1,eigenVectors+iii);

                     }
                 }

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
                
                c2->i.around  += step;
                c2->i.D *= pow( ( 2*c1->i.around +1.) /(2*c2->i.around+1. ),1);

               // c2->i.vectorMomentum *= (c1->i.d* (2.* c1->i.epi + 1.) )/(c1->i.d*(2.*c2->i.epi + 1. ));
                printf("\n\n\t| Band Pass Staging \n\t| From \tBox %d \t: Lattice  %1.3f \t\n\t| To \tBox %d \t: Lattice %1.3f\t\t\n\t\n",2* c1->i.epi + 1,c1->i.d,2* c2->i.epi + 1,c2->i.d);
                
                if ( SPACE > COMPONENT )
                    printf("\n\n\t| Band Pass Staging \n\t| From \tBox %d \t: Lattice  %1.3f \t\n\t| To \tBox %d \t: Lattice %1.3f\t\t\n\t\n",2* c1->i.around + 1,c1->i.D,2* c2->i.around + 1,c2->i.D);


                {
                    {
                        
                        {
                            nStatesTrans = 0;
                            nStatesFound = 0;
                            
                            INT_TYPE ii = 0,iii,irrep;
                            for ( irrep = 0 ; irrep <= 5 ; irrep++){
                                for ( iii = 0; iii < c1->i.nStates ; iii++)
                                {
                                 //   if (f1->sinc.tulip[eigenVectors+iii].value.symmetry == irrep && ((! c1->i.irrep) || irrep == c1->i.irrep) )
                                    {
                                        ii++;
                                        
                                        if ( ii <= c->i.nTargets){
                                            nStatesTrans = imax(nStatesTrans,iii+1);
                                            nStatesFound = ii;
                                        }
                                    }
                                }
                                
                            }
                            
                            c2->i.heliumFlag = nStatesFound;
                            c2->i.nStates = nStatesFound;
                        }
                        
                        if (  nStatesFound == 0 )
                        {
                            printf("no states\n");
                            exit(0);
                        }
                        
                        INT_TYPE ii;
                        {
                            INT_TYPE ii,iii,iv=0;
                            for ( iii = 0; iii < nStatesTrans ; iii++){
                                ii = xConstructFoundation (c2 , 0,0, c1,eigenVectors+iii,1 ,1);
                                iv += ii;
                            }



                            fflush(stdout);
                            c2->i.qFloor = iv;
                            EV = iv;
                            RdsSize = EV;
                        }
                        
                        if ( c1->i.vectorOperatorFlag )
                            c2->i.vectorOperatorFlag = c1->i.vectorOperatorFlag;
                        
                        iModel(c2);

                        xConstructFoundation (c2 , eigenVectors, c2->i.nStates, c1,eigenVectors,nStatesTrans ,1);

                        fModel(c1);

                        free(c1);
                        c1 = c2;
                        f1 = &c1->i.c;
                        {
                            INT_TYPE iii ;
                            for ( iii = 0; iii < c1->i.heliumFlag ; iii++){
                                printf ( "\n Pass-Vector \t%d\n", iii+1);
                                printExpectationValues(f1, Ha, eigenVectors+iii);
                            }
                        }
//                        if ( c1->i.vectorOperatorFlag ){
//                            
//                            if ( c1->i.bodyVectorOperator != f1->sinc.tulip[f1->sinc.vectorOperator].space[0].body  ){
//                                printf("body count!\n");
//                                exit(0);
//                                
//                            }
//                            
//                            if ( tLoadEigenWeights ( c1, c1->mem.operatorName ,f1->sinc.vectorOperator) != c1->i.vectorOperatorFlag){
//                                printf("set helium %d \n", c1->i.vectorOperatorFlag);
//                                exit(0);
//                            }
//                            
//                        }
                        
                        INT_TYPE sp;
#ifdef OMP
#pragma omp parallel for private (ii,sp) schedule(dynamic,1)
#endif
                        for ( ii = 0; ii < c1->i.nStates ; ii++)
                            for ( sp = 0 ;sp < spins(f1, eigenVectors+ii);sp++){
                                tEqua(f1, f1->sinc.user+ii, sp, eigenVectors+ii , sp);
                            }

                        EV =c1->i.nStates;
                        RdsSize = EV;

                    }
                }

            }
            
            fprintf(out,"\n\n\t\t| Foundation \t: %d\t\t\n\t\t| Krylov Space \t: %d\t\t \n\t\t| Irrep \t\t: %d\t\t\t\n", EV,RdsSize, c1->i.irrep);
            fprintf(out , "\t\t| Body \t\t\t: %d\t\t\t\n\t\t| Box \t\t\t: %d\t\t\n\t\t| Lattice \t\t: %1.3f\t\t\n\n",c1->rt.body,2*c1->i.epi+1,c1->i.d );
            fflush(stdout);

                
            if ( iter < c1->i.Iterations)//extra catch-all
            
            {
                
                
                if   ( tEigenCycle(&c1->i.c,Ha,'T', c1->i.heliumFlag, f1->sinc.user,RdsSize,0, EV, 0,0,eigenVectors,twoBodyRitz))
                {
                    RdsSize -= EV;
                    f1->sinc.user += EV;
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
                            fprintf(out,"Failure%d:%d:,%lld ,%1.15f\n", ii+1, c1->i.epi*2+1,ii+1,f1->sinc.tulip[eigenVectors+iii].value);
                            ii++;
                        }
                    }
                    fflush(stdout);

                    
                    continue;//break from cycling/iterating...if failed LD test.
                }
                tFilter(f1, c1->i.nStates, !(!c1->i.filter )* c->i.irrep, f1->sinc.user);
                {
                    INT_TYPE iii ;
                    for ( iii = 0; iii < c1->i.heliumFlag ; iii++){
                        printf ( "\n Vector \t%d\n", iii+1);
                        printExpectationValues(f1, Ha, eigenVectors+iii);
                      //  tEdges(c1, eigenVectors+iii);

                    }
                }

                RdsSize += tGreatDivideIteration(f1,Iterator, 1,0,f1->sinc.user+RdsSize-EV,EV,2*EV,0)-EV;
                
                while ( RdsSize > EV*(c1->i.lookBack+1)){
                    RdsSize -= EV;
                    f1->sinc.user += EV;
                    tClear(f1,matrixHbuild);
                    printf("forgetting...\n");
                    fflush(stdout);
                }
                
            }else {
                cycleFlag = 1;//not likely to happend unless errro.
                printf("woops! ... looks like you accidently attempted an extra iteration\n");
            }
            
            INT_TYPE ii = 0,iii;
            for ( iii = 0; iii < c1->i.heliumFlag ; iii++)
            {
                if (1){
                    INT_TYPE irrep = f1->sinc.tulip[eigenVectors+iii].value.symmetry ;
                    fprintf(out,"%dCondition%d:,%1.15f, %d , %1.1f,%1.6f\n", bodies(f1,eigenVectors+iii),ii+1, f1->sinc.tulip[eigenVectors+iii].value.value,irrep, deg(f1, irrep),f1->sinc.tulip[eigenVectors+iii].value.value2);
                    ii++;
                }
            }
            fflush(stdout);
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
    
    
    
    
#define sincFlag 1

    
    
#if sincFlag

//    tTestSA(two, 2);
//    tTestSA(three, 3);
//    tTestSA(four, 5);
//    exit(0);
    
    struct calculation c = bootShell(argc, argv);
    exec(&c);
    print(&c,0,c.i.nStates,eigenVectors);
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
    if ( argc < 4 ){
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
    DCOMPLEX va;

    INT_TYPE acc = 1;//somehwo?
    
    INT_TYPE body = 2;
    struct general_2index g3[3];
    enum basisElementType bee = GaussianBasisElement;
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
                        
                        g3[i].i[0].bra.basis = bee;
                        g3[i].i[0].bra.length = b;//body 1 Matrix
                        g3[i].i[0].bra.origin = 0;
                        g3[i].i[0].bra.index = l1;
                        
                        g3[i].i[0].ket.basis = bee;
                        g3[i].i[0].ket.length = b;
                        g3[i].i[0].ket.origin = 0 ;
                        g3[i].i[0].ket.index = l2 ;
                        
                        if(1){
                            g3[i].i[1].bra.basis = bee;
                            g3[i].i[1].bra.length = 1;//body 1 Matrix
                            g3[i].i[1].bra.origin = 0;
                            g3[i].i[1].bra.index = l3;
                            
                            g3[i].i[1].ket.basis = bee;
                            g3[i].i[1].ket.length = 1;
                            g3[i].i[1].ket.origin = 0;
                            g3[i].i[1].ket.index = l4;
                        }else {
                            g3[i].i[1].bra.basis = DiracDeltaElement;
                            g3[i].i[1].bra.length = 1;//body 1 Matrix
                            g3[i].i[1].bra.origin = 0;
                            g3[i].i[1].bra.index = l3;
                            
                            g3[i].i[1].ket.basis = nullBasisElement ;
                            g3[i].i[1].ket.length = 1;
                            g3[i].i[1].ket.origin = 0;
                            g3[i].i[1].ket.index = l4;

                        }
                        
                        g3[i].fl = &func1;
                        
                        
                        
                        
                        //                printf("%f %f\n", creal(aaGdnGdm(atoi(argv[3]),atoi(argv[4]), &g3[0].i[0])));
                        //                exit(0);
                        
                        //
                         va  = collective(1., &g3[0]);
                        
                        
                        
//                        for ( i =0; i< 1; i++)
                  //          va = FGS(0,&g3[0].i[0]);
//
                 //       printf("%d\t%d\t%d\t%d\t %f\t%f\n",l1,l2,l3,l4, creal(va), cimag(va));
//                        //        exit(0);
                        //
                   //     exit(0);
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

        g3[i].i[0].bra.basis = bee;
        g3[i].i[0].bra.length = b;//body 1 Matrix
        g3[i].i[0].bra.origin = x;
        g3[i].i[0].bra.index = l;
        
        
        g3[i].i[0].ket.basis = bee;
        g3[i].i[0].ket.length = b;
        g3[i].i[0].ket.origin = x;
        g3[i].i[0].ket.index = l;
        
        g3[i].i[1].bra.basis = bee;
        g3[i].i[1].bra.length = b;//body 1 Matrix
        g3[i].i[1].bra.origin = 0;
        g3[i].i[1].bra.index = l;
        
        g3[i].i[1].ket.basis = bee;
        g3[i].i[1].ket.length = b;
        g3[i].i[1].ket.origin = 0;
        g3[i].i[1].ket.index = l;

        g3[i].fl = & func1;
        
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

        g3[i].i[0].bra.basis = bee;
        g3[i].i[0].bra.length = b;//body 1 Matrix
        g3[i].i[0].bra.origin = x;
        g3[i].i[0].bra.index = l;
        
        
        g3[i].i[0].ket.basis = bee;
        g3[i].i[0].ket.length = b;
        g3[i].i[0].ket.origin = x;
        g3[i].i[0].ket.index = l;
        
        g3[i].i[1].bra.basis = bee;
        g3[i].i[1].bra.length = b;//body 1 Matrix
        g3[i].i[1].bra.origin = 0;
        g3[i].i[1].bra.index = l;
        
        g3[i].i[1].ket.basis = bee;
        g3[i].i[1].ket.length = b;
        g3[i].i[1].ket.origin = 0;
        g3[i].i[1].ket.index = l;


        g3[i].fl = & func1;
        
    }

    {
        INT_TYPE Nt = 10000;
        time_t start_t, lapse_t;
        time(&start_t);

        
#if 1
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
    
#else
    printf("%f\n", elementCal(1e-3,-1, g3));
        exit(0);
#endif
    }
#endif

}

#endif

INT_TYPE buildElectronFreeInteraction ( struct calculation * c1, enum division mat){
    double value;
    struct field * f1 = &c1->i.c;
    INT_TYPE g,space,r,i,j,n,m,N1 = c1->i.epi*2+1;
    INT_TYPE N2 = N1*N1;

    if ( ! CanonicalRank(f1, interactionExchange, 0)){
        printf("canceled ee\n");
        return 0;
    }
    struct calculation c0 = *c1;
//    for( g = 0; g < nSAG*nSAG*nSAG; g++)
//        c0.i.cSA[g] = 0;

    
    c0.mem.bootedMemory = 0;
    c0.i.c.sinc.tulip = NULL;
    c0.i.c.sinc.rose[0].stream = NULL;
    c0.i.c.sinc.rose[1].stream = NULL;
    c0.i.c.sinc.rose[2].stream = NULL;
    c0.i.c.sinc.rose[3].stream = NULL;
    c0.i.nTargets = 1;
    c0.rt.printFlag = 0;
    c0.rt.calcType = electronicStuctureCalculation;
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
    
   // f1->sinc.tulip[mat].stop[0][0] = CanonicalRank(f1, mat,0);
    
    return 0;
}
