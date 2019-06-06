/*
 *  Model.c
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

#include "Model.h"

void resetExternal(struct calculation * i, INT_TYPE number, double scale ){
    if ( number == 1 ){
        i->i.c.Na = 1;
        i->i.c.atoms[1].label.Z = i->rt.body;
        i->i.c.atoms[1].position[1] = 0;
        i->i.c.atoms[1].position[2] = 0;
        i->i.c.atoms[1].position[3] = 0;
        
    } else
        if ( number== 2 ){
            i->i.c.Na = 2;

            i->i.c.atoms[1].label.Z = 1;
            i->i.c.atoms[1].position[1] = scale;
            i->i.c.atoms[1].position[2] = 0;
            i->i.c.atoms[1].position[3] = 0;
            i->i.c.atoms[2].label.Z = 1;
            i->i.c.atoms[2].position[1] = -scale;
            i->i.c.atoms[2].position[2] = 0;
            i->i.c.atoms[2].position[3] = 0;
            
        }
        else
            if ( number == 3 ){
                i->i.c.Na = 3;

                i->i.c.atoms[1].label.Z = 1;
                i->i.c.atoms[1].position[1] = 2*scale;
                i->i.c.atoms[1].position[2] = -scale;
                i->i.c.atoms[1].position[3] = -scale;
                i->i.c.atoms[2].label.Z = 1;
                i->i.c.atoms[2].position[1] = -scale;
                i->i.c.atoms[2].position[2] = -scale;
                i->i.c.atoms[2].position[3] = 2*scale;
                i->i.c.atoms[3].label.Z = 1;
                i->i.c.atoms[3].position[1] = -scale;
                i->i.c.atoms[3].position[2] = 2*scale;
                i->i.c.atoms[3].position[3] = -scale;
                
            }
            else
                if ( number == 4 ){
                    
                    i->i.c.Na = 4;

                    i->i.c.atoms[1].label.Z = 1;
                    i->i.c.atoms[1].position[1] = scale;
                    i->i.c.atoms[1].position[2] = 0;
                    i->i.c.atoms[1].position[3] = -scale/sqrt(2.);
                    i->i.c.atoms[2].label.Z = 1;
                    i->i.c.atoms[2].position[1] = -scale;
                    i->i.c.atoms[2].position[2] = 0;
                    i->i.c.atoms[2].position[3] = -scale/sqrt(2.);
                    i->i.c.atoms[3].label.Z = 1;
                    i->i.c.atoms[3].position[1] = 0;
                    i->i.c.atoms[3].position[2] = scale;
                    i->i.c.atoms[3].position[3] = scale/sqrt(2.);
                    i->i.c.atoms[4].label.Z = 1;
                    i->i.c.atoms[4].position[1] = 0;
                    i->i.c.atoms[4].position[2] = -scale;
                    i->i.c.atoms[4].position[3] = scale/sqrt(2.);
                    
                }
                else
                    if ( number == 5 ){
                        
                        i->i.c.Na = 1;
                        
                        i->i.c.atoms[1].label.Z = 2;
                        i->i.c.atoms[1].position[1] = 0;
                        i->i.c.atoms[1].position[2] = 0;
                        i->i.c.atoms[1].position[3] = 0;
                        
                    }

        i->i.charge = i->i.c.Na - i->rt.body;
    }


struct calculation initCal (void ) {
    INT_TYPE space;
    struct calculation i;
    i.rt.printFlag = 0;
    i.rt.boot = fullMatrices;
//    i.i.hartreeFockFlag = 0;
    i.i.springConstant = 0.;
    i.i.springFlag = 0;
    i.i.potentialFlag = 0;
    i.i.RAMmax = 4;
    i.rt.runFlag = 0;
    i.i.vectorMomentum = 0.;
    i.i.decomposeRankMatrix = 10;
    i.i.orgClamp = 2.;
    i.i.Angstroms = 0;
    i.i.cycles = 1;
    i.i.Iterations = 1;
    i.rt.printFlag = 0;
    i.i.heliumFlag = 1;

    i.rt.TARGET = 1e-2;
    i.rt.targetCondition = 1e-7;
    i.rt.ALPHA = 1e-8;
    i.rt.CANON = 1e-6;
    i.rt.vCANON = 1e-2;
    i.rt.TOL = 1e5;
    i.rt.maxEntropy = 1;
    i.i.level = 10000;
    i.i.complexType =2;
    i.i.sectors = 1;
    i.i.d = 1.;
    i.i.D = 2.;
    i.i.epi = 2;
        i.i.bRank = 4 ;
        i.i.cycleStep = 1;
        i.i.cycles = 25;
        i.i.Iterations = 12;
        i.i.lookBack = 12;
    i.i.level = 1;
        i.i.nStates =  3 ;
        i.i.nTargets = 3;
        i.i.qFloor = 50;
        
        i.i.iRank = 1;
        i.i.side = 5;
        i.i.l2 = 1000;
        
        i.i.attack = 0.66;
        
        
        
        i.rt.samples = 0;
        i.rt.monteCarlo = 0;
        i.i.turn = 1.;
        i.i.param1 = 1.;
        i.i.param2 = 1.;
        i.i.interval = 1;
        i.i.scalar = 1.;
        i.i.outputFlag = 0;
        
        i.i.vectorOperatorFlag = 0;
        i.i.dRank = 0;
    
        i.i.magFlag = 0;
        i.i.mag = 0.1;
        //THESE
        if ( SPACE == 3 ){
            i.rt.calcType = electronicStuctureCalculation;
            i.rt.body = two;
            i.i.irrep = 2;
            i.i.filter = 1;

        } else {
            i.rt.body = one;
            i.rt.calcType = clampProtonElectronCalculation;
            i.i.around = 7 ;
            i.i.D = 0.1;
            i.i.irrep = 0;

        }
    i.i.massElectron = 1.;
    i.i.massProton = 1836.15267245;
    i.i.massClampPair = 1836.15267245;
    resetExternal(&i, 1, 1);
    
    for( space = 0 ; space <= SPACE ; space++){
        i.i.c.sinc.rose[space].stream = NULL;
    }
    i.i.c.sinc.tulip = NULL;
    i.rt.phaseType = buildFoundation ;
    
    i.i.c.twoBody.func.fn = Coulomb;
    i.i.c.oneBody.func.fn = nullFunction;

    i.i.canonRank = 1000;
    i.i.c.twoBody.num = 25;
    i.i.c.twoBody.func.interval  = 0;
    i.i.c.twoBody.func.param[0]  = 1;
    i.i.c.twoBody.func.param[1]  = 1;
    i.i.c.twoBody.func.param[2]  = 1;

    i.i.c.oneBody.num = 25;
    i.i.c.oneBody.func.interval  = 0;
    i.i.c.oneBody.func.param[0]  = 1;
    i.i.c.oneBody.func.param[1]  = 1;
    i.i.c.oneBody.func.param[2]  = 1;
    i.i.iCharge = i.rt.body;
    i.i.c.body = i.rt.body;

    i.mem.bootedMemory = 0;
    return i;
}


struct calculation gas (void ) {
    struct calculation c = initCal();
    c.rt.runFlag = 3;
    return c;
}



INT_TYPE fModel ( struct calculation * c1 ){
    INT_TYPE i;
    struct field * f1  =&(c1->i.c);
    
    if ( c1->mem.bootedMemory ){
        for ( i = 0;i <= SPACE ; i++){
          //  fprintf(stdout, "SPACE %d\n", i);
           // fflush(stdout);
            free(f1->sinc.rose[i].stream);
            f1->sinc.rose[i].stream = NULL;
        }
        {
          //  fprintf(stdout, "TULIP\n" );
          //  fflush(stdout);

            free(f1->sinc.tulip);
            f1->sinc.tulip = NULL;
        }
    }
    return 0;
}


INT_TYPE iModel( struct calculation * c1){
    //c1->rt.printFlag = 1;
    struct name_label l2;;

    struct field *f1 = &(c1->i.c);
    {
        INT_TYPE space ;
        if ( f1->sinc.tulip != NULL ){
            printf("iModel\n");
            exit(0);
        }
        for ( space = 0; space <= SPACE  ; space++)
            if ( f1->sinc.rose[space].stream != NULL ){
                printf("iModel rose\n");

                exit(0);
            }
    }
    if(0){
        printf("\n\n*Parameters\n");
        printf("\t target\t\t\t%1.1f\n", -log(c1->rt.TARGET)/log(10));//quality of decomposition
        printf("\t tolerance\t\t%1.1f\n", log(c1->rt.TOL)/log(10));//max condition of foundation    s
        printf("\t threshold\t\t%1.1f\n", -log(c1->rt.CANON)/log(10));//matrix training standard
        printf("\t vectorThreshold\t%1.1f\n", -log(c1->rt.vCANON)/log(10));//vector training standard
        printf("\t condition\t\t%1.1f\n", -log(c1->rt.ALPHA )/log(10));//Beylkin parameter
        printf(".Parameters\n\n");
    }
    struct runTime * rt = &c1->rt;
    if ( c1->i.iCharge > 0 )
        f1->Ne = c1->i.iCharge;
    
    
        INT_TYPE space;
        {
            c1->mem.rt = rt;
            f1->mem1 = &c1->mem;
            
        }
        
        enum bodyType bootBodies = c1->rt.body;
        c1->i.c.body = bootBodies;
        INT_TYPE ra = tPerms(bootBodies),nG = c1->i.decomposeRankMatrix;//tSize(bootBodies);
        INT_TYPE N12 = c1->i.epi;
        INT_TYPE N1 =  2*N12+1;
        enum shape bootShape;
        INT_TYPE maxVector = imax(c1->i.decomposeRankMatrix, imax(c1->i.bRank,imax(1+c1->i.iRank,c1->i.Iterations+c1->i.iRank)));
        //rds defined in input.c
        
        bootShape = Cube;
        f1->sinc.d = c1->i.d;
        
        
        INT_TYPE FloorV = imax(0, c1->i.qFloor), CeilV = imax(0,0);
        INT_TYPE maxArray,EV,maxEV,NV = 0,FV = FloorV+CeilV ;
        INT_TYPE maxDensity = c1->i.vectorOperatorFlag;
        
        EV = FV;
        maxEV =EV*(imax(c1->i.Iterations,1));
        maxArray = imax(c1->i.nStates,maxEV);//slip Nb into spectra...
        
        f1->sinc.maxEV = maxArray;
      //  printf("states %d , maxEV %d\n", c1->i.nStates, maxEV);
        enum division vectorOperator  = eigenVectors +  c1->i.nStates+maxEV;
        f1->sinc.vectorOperator = vectorOperator;
        enum division end  = vectorOperator+maxDensity;
        //printf("end %d\n, vectorOp %d\n", end, vectorOperator);
        f1->sinc.end = end;
        f1->sinc.tulip = malloc ( (end+1) * sizeof(struct name_label));


        INT_TYPE matrixNumber =   c1->i.decomposeRankMatrix;
        INT_TYPE outVector  = imax(f1->oneBody.num , f1->twoBody.num )*maxVector;
        
        {//defaults
            INT_TYPE periodic = 1;
            //define vectors
            f1->sinc.rose[SPACE].component = nullComponent;
            f1->sinc.rose[SPACE].body = nada;
            f1->sinc.rose[SPACE].particle = nullParticle;
            if ( c1->rt.calcType == electronicStuctureCalculation){
                periodic = 1;
                for ( space = 0; space < SPACE ; space++){
                    f1->sinc.rose[space].particle = electron;
                    f1->sinc.rose[space].basis = SincBasisElement;
                    f1->sinc.rose[space].body = bootBodies;
                    f1->sinc.rose[space].component = spatialComponent1+space%COMPONENT + ((c1->rt.runFlag/periodic)%2)*3;
                    periodic *= 2;
                    f1->sinc.rose[space].lattice = c1->i.d;
                    f1->sinc.rose[space].count1Basis = N1;
                    if (f1->sinc.rose[space].component > 3  )
                        f1->sinc.rose[space].count1Basis *= 2;
                    f1->sinc.rose[space].origin = 0.;
                }
            } else  if ( c1->rt.calcType == clampProtonElectronCalculation){
                periodic = 1;
                for ( space = 0; space < COMPONENT ; space++){
                    f1->sinc.rose[space].particle = electron;
                    f1->sinc.rose[space].basis = SincBasisElement;
                    f1->sinc.rose[space].body = bootBodies;
                    f1->sinc.rose[space].component = spatialComponent1+space+ ((c1->rt.runFlag/periodic)%2)*3;
                    periodic *= 2;
                    f1->sinc.rose[space].lattice = c1->i.d;
                    f1->sinc.rose[space].count1Basis = N1;
                    f1->sinc.rose[space].origin = 0.;

                }
                periodic = 1;
                for ( space = COMPONENT; space <= COMPONENT ; space++){
                    f1->sinc.rose[space].particle = proton;
                    f1->sinc.rose[space].basis = SincBasisElement;
                    f1->sinc.rose[space].body = one;
                    f1->sinc.rose[space].component = spatialComponent1+space%COMPONENT+ ((c1->rt.runFlag/periodic)%2)*3;
                    periodic *= 2;
                    f1->sinc.rose[space].lattice = c1->i.D;
                    f1->sinc.rose[space].count1Basis = 2*c1->i.around+1;
                    f1->sinc.rose[space].origin = c1->i.orgClamp;
                }
                for ( space = COMPONENT+1; space < SPACE ; space++){
                    f1->sinc.rose[space].particle = proton;
                    f1->sinc.rose[space].basis = SincBasisElement;
                    f1->sinc.rose[space].body = nada;
                    f1->sinc.rose[space].component = spatialComponent1+space%COMPONENT+ ((c1->rt.runFlag/periodic)%2)*3;
                    periodic *= 2;
                    f1->sinc.rose[space].lattice = 0;
                    f1->sinc.rose[space].count1Basis = 0 ;
                    f1->sinc.rose[space].origin = 0.;
                }
                
            }
                //define vectors end
            {
                enum division label1;
                for ( label1 = 0 ;label1 <= end; label1++){
                    f1->sinc.tulip[label1].name = label1;
                    f1->sinc.tulip[label1].Partition = 0;
                    f1->sinc.tulip[label1].header = Cube;
                    f1->sinc.tulip[label1].spinor = real;

                    f1->sinc.tulip[label1].species = scalar;
                    f1->sinc.tulip[label1].linkNext = nullName;
                    f1->sinc.tulip[label1].memory = objectAllocation;
                    for ( space = 0; space <= SPACE ; space++)
                        f1->sinc.tulip[label1].space[space].Address = -1;
                    f1->sinc.tulip[label1].space[SPACE].block = id0;
                    for ( space = 0; space <= SPACE ; space++){
                        {
                            f1->sinc.tulip[label1].space[space].block = id0;//matrix prototype
                            f1->sinc.tulip[label1].space[space].body = nada;//matrix prototype
                        }
                    }
                    tClear(f1,label1);
                }
                
            }
        }//defaults
        fromBeginning(f1, kinetic, 0);
        f1->sinc.tulip[kinetic].Partition = COMPONENT;//
//        if ( c1->rt.runFlag > 0 )
//            f1->sinc.tulip[kinetic].spinor = cmpl;
        assignOneWithPointers(f1, kinetic,all);

        fromBeginning(f1, kineticMass, kinetic);
        f1->sinc.tulip[kineticMass].Partition = COMPONENT;//
        assignOneWithPointers(f1, kineticMass,all);
        struct name_label u = f1->sinc.tulip[kineticMass];

        fromBeginning(f1, protonRepulsion, kineticMass);
        f1->sinc.tulip[protonRepulsion].Partition = f1->oneBody.num;//
        assignOneWithPointers(f1, protonRepulsion,all);
        
        fromBeginning(f1, vectorMomentum, protonRepulsion);
        assignOneWithPointers(f1, vectorMomentum,electron);
        f1->sinc.tulip[vectorMomentum].spinor = cmpl;
        f1->sinc.tulip[vectorMomentum].Partition = COMPONENT * c1->i.springFlag+4*c1->i.magFlag;//

        fromBeginning(f1, harmonium, vectorMomentum);
        assignOneWithPointers(f1, harmonium,electron);
        f1->sinc.tulip[harmonium].Partition = 0;//
        
        fromBeginning(f1, X, harmonium);
        assignOneWithPointers (f1, X,electron);
        f1->sinc.tulip[X].Partition =  0 ;//make it a semi-local Gaussian * x
        
        fromBeginning(f1, linear, X);
        assignOneWithPointers (f1, linear,electron);
        f1->sinc.tulip[linear].Partition = f1->Na*f1->oneBody.num+f1->twoBody.num ;//

        fromBeginning(f1, overlap, linear);
        assignOneWithPointers (f1, overlap,electron);
        f1->sinc.tulip[overlap].Partition = 0;//

        {
            fromBeginning(f1, build, overlap);
            f1->sinc.tulip[build].Partition = c1->i.OCSBflag*!(!c1->i.sectors)*bootBodies*(6+(/*HERE 3*/part(f1, kinetic)+part(f1, vectorMomentum)) + matrixNumber);//easily reduce in cheaper ways!
            if ( bootBodies == two )
                f1->sinc.tulip[build].Partition += c1->i.OCSBflag*c1->i.sectors;
            else if ( bootBodies == three )
                f1->sinc.tulip[build].Partition += c1->i.OCSBflag*3*c1->i.sectors;
            else if ( bootBodies == four )
                f1->sinc.tulip[build].Partition += c1->i.OCSBflag*6*c1->i.sectors;
            
            f1->sinc.tulip[build].Partition = c1->i.OCSBflag*c1->i.sectors*imax(nG*c1->i.sectors,part(f1, build));//BUILD
            
            //f1->sinc.tulip[build].Partition = 0;
            // B2 :  2*(1+S)*3 + 1*sector +  Na*matrixNumber*2
            // B3 :  3*(1+S)*3 + 3*sector  + Na*matrixNumber*3
            f1->sinc.tulip[build].species = matrix;
            assignParticle(f1, build, all, bootBodies);
            INT_TYPE n1[SPACE];
            length1(f1,n1);
       //     printf("sectors %d\n", c1->i.sectors);
            for ( space = 0; space < SPACE ; space++)
            {
                if ( space == 0 )
                    fromBeginning(f1, bill1+space, build);
                else
                    fromBeginning(f1, bill1+space, bill1+space-1);
                
                f1->sinc.tulip[bill1+space].Partition = !(!c1->i.sectors)*pow(imin( n1[space], (bootBodies != one ) * c1->i.decomposeRankMatrix + (bootBodies == one )* n1[space]),c1->rt.body) * vectorLen(f1, space) ;
                f1->sinc.tulip[bill1+space].memory = bufferAllocation;
               // if ( c1->rt.runFlag )
               //     f1->sinc.tulip[bill1+space].spinor = cmpl;
                
            }
            fromBeginning(f1, eigen, bill1+SPACE-1);
            f1->sinc.tulip[eigen].Partition = c1->i.OCSBflag * !(!c1->i.sectors)*c1->i.decomposeRankMatrix;
            f1->sinc.tulip[eigen].species = matrix;
            assignParticle(f1, eigen, all, bootBodies);
            for ( space = 0; space < SPACE ; space++)
                if ( f1->sinc.rose[space].body != nada)
                    f1->sinc.tulip[eigen].space[space].block = tv1;
            
            
            {
                INT_TYPE di,cmpl;
                enum division last = eigen;
                INT_TYPE booting[7];
                for ( di = 0; di < 7 ; di++)
                    booting[di] = 0;
                
                
                if ( c1->i.vectorOperatorFlag )
                    
                {
                    enum bodyType bd;
                    INT_TYPE fi,lines = 0;
                    size_t ms = MAXSTRING;
                    char line0[MAXSTRING];
                    char name[MAXSTRING];
                    char *line = line0;
                    INT_TYPE FIT ;
                    FIT = c1->mem.filesVectorOperator ;
                    for ( fi =0 ; fi < FIT; fi++){
                        strcpy(name ,c1->mem.fileVectorOperator[fi]);
                        FILE * fp = fopen(name,"r");
                        if ( fp == NULL ) {
                            printf("file?\n");
                            exit(0);
                        }
                        getline(&line, &ms, fp);
                        while(!feof(fp))
                        {
                            if (! comment(line))
                            {
                                cmpl = 0;
                                tFromReadToFilename(NULL, line,  name, spins(f1,eigenVectors)-1,cmpl);
                              //  printf("%s\n", name);
                                fromBeginning(f1, vectorOperator+lines, last);
                                f1->sinc.tulip[vectorOperator+lines].Partition = inputFormat(f1, name, nullName, 2);
                                f1->sinc.tulip[vectorOperator+lines].species = outerVector;
                                for (space = 0; space < SPACE ; space++){
                                    bd = inputFormat(f1, name, nullName, 100+space/COMPONENT);
                                    f1->sinc.tulip[vectorOperator+lines].space[space].body = bd;
                                  //  printf("%d %d %d\n", lines+1, space, bd);//space/COMPONENT = particle
                                    booting[ bd - bootBodies ] = 1;

                                }
                                
                                last = vectorOperator+lines;
                                lines++;
                            }
                            getline(&line, &ms, fp);
                        }
                        if ( fi > MAXSTRING)
                        {
                            printf("too many files, increase MAXSTRING\n");
                            exit(0);
                        }
                        fclose(fp);
                    }
                }
                
                
//                for ( di = 1; di < 4 ; di++)
//                {
//                    fromBeginning(f1, complement+di, last);
//                    last = complement + di;
//                    f1->sinc.tulip[complement+ di].Partition = booting[di] ;//
//                    f1->sinc.tulip[complement+ di].spinor = parallel;
//                    f1->sinc.tulip[complement+ di].species = outerVector;
//                    assignParticle(f1, complement+di, electron, di);
//                    
//                    fromBeginning(f1,complementTwo+di , last);
//                    last = complementTwo + di;
//                    f1->sinc.tulip[complementTwo+ di].Partition = booting[di] ;//
//                    f1->sinc.tulip[complementTwo+ di].spinor = parallel;
//                    f1->sinc.tulip[complementTwo+ di].species = outerVector;
//                    assignParticle(f1, complementTwo+di, electron, di);
//
//                }
                
                fromBeginning(f1, eigenVectors, last);
            }
            
            struct name_label u ;
            struct name_label u2;


            if(1 /* !si*/){
                INT_TYPE di,d0=1;
#if VERBOSE
                printf("std USERs\n\n");
#endif
                f1->sinc.tulip[eigenVectors].Partition = c1->i.bRank;
                f1->sinc.tulip[eigenVectors].species = vector;
                f1->sinc.tulip[eigenVectors].spinor = c1->i.complexType;
                for ( di = 1 ; eigenVectors+di < vectorOperator; di++){
                    fromBeginning(f1, eigenVectors+di, eigenVectors+di-1);
                    f1->sinc.tulip[eigenVectors+di].spinor = c1->i.complexType;
                    if ( di < c1->i.nStates ){
                        f1->sinc.tulip[eigenVectors+di].Partition = c1->i.bRank;
                        d0++;
                    }
                    else if ( di < c1->i.nStates+maxEV){
                        {
                            f1->sinc.tulip[eigenVectors+di].Partition = ((di-d0)/EV)+c1->i.iRank;
                            NV += spins(f1,eigenVectors+di )*f1->sinc.tulip[eigenVectors+di].Partition;
                        }
                        
                    }else{
                        exit(1);
                    }
                        f1->sinc.tulip[eigenVectors+di].species = vector;
                    }
                fromBeginning(f1, diagonalVectorA, eigenVectors+di-1);
                ;
                 u = f1->sinc.tulip[eigenVectors];
                 u2 = f1->sinc.tulip[eigenVectors+di-1];
                f1->sinc.user = eigenVectors + d0;
            }
                
            

            f1->sinc.tulip[diagonalVectorA].spinor = parallel;
            f1->sinc.tulip[diagonalVectorA].Partition = !(!(c1->i.sectors));
            f1->sinc.tulip[diagonalVectorA].species = vector;
            
            fromBeginning(f1, diagonalVectorB, diagonalVectorA);
            f1->sinc.tulip[diagonalVectorB].spinor = parallel;
            f1->sinc.tulip[diagonalVectorB].Partition =  !(!(c1->i.sectors));
            f1->sinc.tulip[diagonalVectorB].species = vector;

        }
        
        INT_TYPE len[SPACE],mxlen=0;
        length(f1, eigenVectors, len);
        for ( space = 0 ; space < SPACE ; space++)
            if(len[space] > mxlen)
                mxlen = len[space];
        
        INT_TYPE mx1len=0;
        length1(f1, len);
        for ( space = 0 ; space < SPACE ; space++)
            if(len[space] > mx1len)
                mx1len = len[space];

        
        fromBeginning(f1,edgeElectronMatrix , diagonalVectorB);
        f1->sinc.tulip[edgeElectronMatrix].Partition = 1 ;//
        f1->sinc.tulip[edgeElectronMatrix].species = matrix;
        assignOneWithPointers(f1, edgeElectronMatrix, electron);
        
        fromBeginning(f1,edgeProtonMatrix,edgeElectronMatrix);
        f1->sinc.tulip[edgeProtonMatrix].Partition = 1 ;//
        f1->sinc.tulip[edgeProtonMatrix].species = matrix;
        assignOneWithPointers(f1, edgeProtonMatrix, proton);

        fromBeginning(f1,productVector,edgeProtonMatrix);
        f1->sinc.tulip[productVector].Partition = maxVector;
        f1->sinc.tulip[productVector].species = vector;
      //  f1->sinc.tulip[productVector].spinor = parallel;

        fromBeginning(f1,permutationVector,productVector);
        f1->sinc.tulip[permutationVector].Partition =  c1->i.filter *ra*maxVector;
        f1->sinc.tulip[permutationVector].species = vector;
        //f1->sinc.tulip[permutationVector].spinor = parallel;
        
        fromBeginning(f1,permutation2Vector,permutationVector);
        f1->sinc.tulip[permutation2Vector].Partition = c1->i.filter *maxVector;
        f1->sinc.tulip[permutation2Vector].species = vector;
//        f1->sinc.tulip[permutation2Vector].spinor = parallel;
    
        fromBeginning(f1,canonicalmvVector,permutation2Vector);
        f1->sinc.tulip[canonicalmvVector].Partition = 1;
        f1->sinc.tulip[canonicalmvVector].species = vector;
        f1->sinc.tulip[canonicalmvVector].spinor = parallel;

        fromBeginning(f1,canonicalmv2Vector,canonicalmvVector);
        f1->sinc.tulip[canonicalmv2Vector].Partition = 1;
        f1->sinc.tulip[canonicalmv2Vector].species = vector;
        f1->sinc.tulip[canonicalmv2Vector].spinor = parallel;

        fromBeginning(f1,canonicalmv3Vector,canonicalmv2Vector);
        f1->sinc.tulip[canonicalmv3Vector].Partition = 1;
        f1->sinc.tulip[canonicalmv3Vector].species = vector;
        f1->sinc.tulip[canonicalmv3Vector].spinor = parallel;

        fromBeginning(f1,canonicaldotVector,canonicalmv3Vector);
        f1->sinc.tulip[canonicaldotVector].Partition = 1;
        f1->sinc.tulip[canonicaldotVector].species = vector;;
        f1->sinc.tulip[canonicaldotVector].spinor = parallel;
        
        fromBeginning(f1,canonicaldot2Vector,canonicaldotVector);
        f1->sinc.tulip[canonicaldot2Vector].Partition = 1;
        f1->sinc.tulip[canonicaldot2Vector].species = vector;;
        f1->sinc.tulip[canonicaldot2Vector].spinor = parallel;
        
        fromBeginning(f1,canonicaldot3Vector,canonicaldot2Vector);
        f1->sinc.tulip[canonicaldot3Vector].Partition = 1;
        f1->sinc.tulip[canonicaldot3Vector].species = vector;;
        f1->sinc.tulip[canonicaldot3Vector].spinor = parallel;
        
        fromBeginning(f1,canonicalvvVector,canonicaldot3Vector);
        f1->sinc.tulip[canonicalvvVector].Partition = 0;
        f1->sinc.tulip[canonicalvvVector].species = vector;;
        f1->sinc.tulip[canonicalvvVector].spinor = parallel;
        
        fromBeginning(f1,canonicalvv2Vector,canonicalvvVector);
        f1->sinc.tulip[canonicalvv2Vector].Partition = 0;
        f1->sinc.tulip[canonicalvv2Vector].species = vector;;
        f1->sinc.tulip[canonicalvv2Vector].spinor = parallel;
        
        fromBeginning(f1,canonicalvv3Vector,canonicalvv2Vector);
        f1->sinc.tulip[canonicalvv3Vector].Partition = 0;
        f1->sinc.tulip[canonicalvv3Vector].species = vector;;
        f1->sinc.tulip[canonicalvv3Vector].spinor = parallel;

        fromBeginning(f1,canonicalmeVector,canonicalvv3Vector);
        f1->sinc.tulip[canonicalmeVector].Partition = 1;
        f1->sinc.tulip[canonicalmeVector].species = vector;;
        f1->sinc.tulip[canonicalmeVector].spinor = parallel;
        
        fromBeginning(f1,canonicalme2Vector,canonicalmeVector);
        f1->sinc.tulip[canonicalme2Vector].Partition = 0;
        f1->sinc.tulip[canonicalme2Vector].species = vector;;
        f1->sinc.tulip[canonicalme2Vector].spinor = parallel;
        
        fromBeginning(f1,canonicalme3Vector,canonicalme2Vector);
        f1->sinc.tulip[canonicalme3Vector].Partition = 0;
        f1->sinc.tulip[canonicalme3Vector].species = vector;;
        f1->sinc.tulip[canonicalme3Vector].spinor = parallel;

        fromBeginning(f1,copyVector,canonicalme3Vector);
        f1->sinc.tulip[copyVector].Partition = 0* maxVector;
        f1->sinc.tulip[copyVector].species = vector;
        f1->sinc.tulip[copyVector].spinor = parallel;

        fromBeginning(f1,copyTwoVector,copyVector);
        f1->sinc.tulip[copyTwoVector].Partition =   0*maxVector;//c1->i.bRank*c1->i.bRank;
        f1->sinc.tulip[copyTwoVector].species = vector;
        f1->sinc.tulip[copyTwoVector].spinor = parallel;

        fromBeginning(f1,copyThreeVector,copyTwoVector);
        f1->sinc.tulip[copyThreeVector].Partition =   0*maxVector;//c1->i.bRank*c1->i.bRank;
        f1->sinc.tulip[copyThreeVector].species = vector;
        f1->sinc.tulip[copyThreeVector].spinor = parallel;

        fromBeginning(f1,copyFourVector,copyThreeVector);
        f1->sinc.tulip[copyFourVector].Partition =    0*maxVector;//c1->i.bRank*c1->i.bRank;
        f1->sinc.tulip[copyFourVector].species = vector;
        f1->sinc.tulip[copyFourVector].spinor = parallel;
        
//        f1->sinc.tulip[oneVector].Address = fromBeginning(f1,copyFourVector);
//        f1->sinc.tulip[oneVector].Partition = !(!c1->rt.printFlag) * 24*imax(f1->oneBody.num, f1->twoBody.num);;
//        f1->sinc.tulip[oneVector].body = one;
//        f1->sinc.tulip[oneVector].species = vector;
//        f1->sinc.tulip[oneVector].header = Cube;
//        f1->sinc.tulip[oneVector].parallel = 2;
//
//        f1->sinc.tulip[twoVector].Address = fromBeginning(f1,oneVector);
//        f1->sinc.tulip[twoVector].Partition = !(!c1->rt.printFlag) * c1->i.bRank;;
//        f1->sinc.tulip[twoVector].body = two;
//        f1->sinc.tulip[twoVector].species = vector;
//        f1->sinc.tulip[twoVector].header = Cube;
//        f1->sinc.tulip[twoVector].parallel = 2;
//        f1->sinc.tulip[twoVector].purpose = Object;
////        f1->sinc.tulip[twoVector].symmetryType = nullSymmetry;
//        f1->sinc.tulip[twoVector].name = twoVector;
//
        fromBeginning(f1,totalVector,copyFourVector);
        f1->sinc.tulip[totalVector].Partition = (!( c1->rt.phaseType == buildFoundation ))*  maxVector*(c1->i.canonRank)*(1 + (ra*c1->i.filter));
        f1->sinc.tulip[totalVector].species = vector;
        //f1->sinc.tulip[totalVector].spinor = parallel;

        fromBeginning(f1,totalFuzzyVector,totalVector);
        f1->sinc.tulip[totalFuzzyVector].Partition = 0*(c1->i.canonRank);
        f1->sinc.tulip[totalFuzzyVector].species = vector;
        f1->sinc.tulip[totalFuzzyVector].spinor = parallel;

        fromBeginning(f1,diagonalCube,totalFuzzyVector);
        f1->sinc.tulip[diagonalCube].Partition = 1;
        f1->sinc.tulip[diagonalCube].species = matrix;
        f1->sinc.tulip[diagonalCube].spinor = parallel;
        assignParticle(f1, diagonalCube, all, one);
        
        fromBeginning(f1,diagonal1VectorA,diagonalCube);
        f1->sinc.tulip[diagonal1VectorA].Partition =1+!(!c1->rt.printFlag)+c1->i.bodyFlag;
        f1->sinc.tulip[diagonal1VectorA].spinor = parallel;
        f1->sinc.tulip[diagonal1VectorA].species = outerVector;
        assignParticle(f1, diagonal1VectorA, all , one);
        
        fromBeginning(f1,diagonal2VectorA,diagonal1VectorA);
        f1->sinc.tulip[diagonal2VectorA].Partition =c1->i.sectors+!(!c1->rt.printFlag)+c1->i.bodyFlag;
        f1->sinc.tulip[diagonal2VectorA].spinor = parallel;
        f1->sinc.tulip[diagonal2VectorA].species = outerVector;
        assignParticle(f1, diagonal2VectorA, all , two);

        fromBeginning(f1,diagonal2VectorB,diagonal2VectorA);
        f1->sinc.tulip[diagonal2VectorB].Partition =c1->i.sectors+!(!c1->rt.printFlag);
        f1->sinc.tulip[diagonal2VectorB].spinor = parallel;
        f1->sinc.tulip[diagonal2VectorB].species = outerVector;
        assignParticle(f1, diagonal2VectorB, all , two);

        fromBeginning(f1,diagonal1VectorB,diagonal2VectorB);
        f1->sinc.tulip[diagonal1VectorB].Partition =c1->i.sectors+!(!c1->rt.printFlag)+c1->i.bodyFlag;
        f1->sinc.tulip[diagonal1VectorB].spinor = parallel;
        f1->sinc.tulip[diagonal1VectorB].species = outerVector;
        assignParticle(f1, diagonal1VectorB, all , one);

        fromBeginning(f1,diagonal1VectorC,diagonal1VectorB);
        f1->sinc.tulip[diagonal1VectorC].Partition =c1->i.sectors+!(!c1->rt.printFlag)+c1->i.bodyFlag;
        f1->sinc.tulip[diagonal1VectorC].spinor = parallel;
        f1->sinc.tulip[diagonal1VectorC].species = outerVector;
        assignParticle(f1, diagonal1VectorC, all , one);

        fromBeginning(f1,diagonal1VectorD,diagonal1VectorC);
        f1->sinc.tulip[diagonal1VectorD].Partition = c1->i.sectors+!(!c1->rt.printFlag)+c1->i.bodyFlag;
        f1->sinc.tulip[diagonal1VectorD].spinor = parallel;
        f1->sinc.tulip[diagonal1VectorD].species = outerVector;
        assignParticle(f1, diagonal1VectorD, all , one);

        fromBeginning(f1,diagonal3VectorA,diagonal1VectorD);
        f1->sinc.tulip[diagonal3VectorA].Partition =0 *(c1->i.sectors+!(!c1->rt.printFlag)+c1->i.bodyFlag);
       // f1->sinc.tulip[diagonal3VectorA].spinor = parallel;
        f1->sinc.tulip[diagonal3VectorA].species = outerVector;
        assignParticle(f1, diagonal3VectorA, all , three);
        
        fromBeginning(f1,bandBasis,diagonal3VectorA);
        f1->sinc.tulip[bandBasis].Partition = mxlen;
        f1->sinc.tulip[bandBasis].memory = bufferAllocation;

        
        fromBeginning(f1,copy,bandBasis);
        f1->sinc.tulip[copy].Partition = imax(c1->i.bRank*c1->i.bRank,imax(mx1len*mx1len, matrixNumber+ imax(f1->oneBody.num,   imax(matrixNumber*f1->Na,outVector) )));
    if ( c1->rt.phaseType == frameDensity)
            f1->sinc.tulip[copy].spinor = parallel;
        f1->sinc.tulip[copy].species = matrix;
        assignParticle(f1, copy, all, one);
        for ( space = 0; space < SPACE ; space++)
            if ( f1->sinc.rose[space].body != nada)
                f1->sinc.tulip[copy].space[space].block = tv1;
        
        fromBeginning(f1,copyTwo,copy);
        f1->sinc.tulip[copyTwo].Partition = matrixNumber;
        f1->sinc.tulip[copyTwo].spinor = parallel;
        f1->sinc.tulip[copyTwo].species = matrix;
        assignParticle(f1, copyTwo, all, one);
        
        fromBeginning(f1,copyThree,copyTwo);
        f1->sinc.tulip[copyThree].Partition = 0*matrixNumber+ imax(f1->oneBody.num,   imax(matrixNumber*f1->Na,outVector) );
        f1->sinc.tulip[copyThree].spinor = parallel;
        f1->sinc.tulip[copyThree].species = matrix;
        assignParticle(f1, copyThree, all, one);
        
        fromBeginning(f1,squareTwo,copyThree);
        f1->sinc.tulip[squareTwo].Partition= c1->i.sectors;//BUILD
        f1->sinc.tulip[squareTwo].species = matrix;
        f1->sinc.tulip[squareTwo].header = Cube;
        assignParticle(f1, squareTwo, all, two);
        
        fromBeginning(f1,inversion,squareTwo);
        f1->sinc.tulip[inversion].Partition= 1;
        f1->sinc.tulip[inversion].species = matrix;
        f1->sinc.tulip[inversion].header = Cube;
        assignParticle(f1, inversion, all, one);
    
        
        fromBeginning(f1,canonicalBuffers,inversion);
        f1->sinc.tulip[canonicalBuffers].Partition = maxVector*maxVector+ imax(NV,maxVector*c1->i.canonRank)*maxVector;
        f1->sinc.tulip[canonicalBuffers].spinor = parallel;

        fromBeginning(f1,trackBuffer,canonicalBuffers);
        f1->sinc.tulip[trackBuffer].Partition = (2*maxVector+1)*maxVector;
        f1->sinc.tulip[trackBuffer].spinor = parallel;
        f1->sinc.tulip[trackBuffer].memory = bufferAllocation;

        fromBeginning(f1,guideBuffer,trackBuffer);
        f1->sinc.tulip[guideBuffer].Partition = c1->i.canonRank*maxVector*maxVector;
        f1->sinc.tulip[guideBuffer].spinor = parallel;
        f1->sinc.tulip[guideBuffer].memory = bufferAllocation;

        fromBeginning(f1,foundationStructure,guideBuffer);
        f1->sinc.tulip[foundationStructure].spinor = parallel;
        f1->sinc.tulip[foundationStructure].Partition = 1;
        f1->sinc.tulip[foundationStructure].species = vector;
        f1->sinc.tulip[foundationStructure].spinor = cmpl;//need two channels
        
        fromBeginning(f1,interactionExchange,foundationStructure);
            f1->sinc.tulip[interactionExchange].Partition = f1->twoBody.num*(( bootBodies > one )|| c1->rt.runFlag > 0);
            f1->sinc.tulip[interactionExchange].species = matrix;
            assignParticle(f1, interactionExchange, electron, two);
    {
        {
            enum block ee ;
            for (ee = e12 ; ee <= e34 ; ee++){
                f1->sinc.tulip[interaction12+ee-e12].name = interactionExchange;
                f1->sinc.tulip[interaction12+ee-e12].species = matrix;
                if ( ee < e34 )
                    f1->sinc.tulip[interaction12+ee-e12].linkNext = interaction12+ee-e12+1;
                for ( space = 0; space < SPACE ; space++)
                    if ( f1->sinc.rose[space].body >= two && f1->sinc.rose[space].particle == electron )
                        f1->sinc.tulip[interaction12+ee-e12].space[space].block = ee;
                f1->sinc.tulip[interaction12+ee-e12].Partition = f1->twoBody.num;
            }
        }
    }
        fromBeginning(f1,interactionExchangeB,interactionExchange);
        f1->sinc.tulip[interactionExchangeB].Partition = f1->twoBody.num*( bootBodies > one )* ( c1->rt.calcType == protonsElectronsCalculation);
        f1->sinc.tulip[interactionExchangeB].species = matrix;
        assignParticle(f1, interactionExchangeB, proton, two);
        {
        enum block ee ;
            for (ee = e12 ; ee <= e34 ; ee++){
                f1->sinc.tulip[interaction12B+ee-e12].name = interactionExchangeB;
                f1->sinc.tulip[interaction12B+ee-e12].species = matrix;
                if ( ee < e34 )
                    f1->sinc.tulip[interaction12B+ee-e12].linkNext = interaction12B+ee-e12+1;
                for ( space = 0; space < SPACE ; space++)
                    if ( f1->sinc.rose[space].body >= two && f1->sinc.rose[space].particle == proton ){
                        f1->sinc.tulip[interaction12B+ee-e12].space[space].block = ee;
                    }
                f1->sinc.tulip[interaction12B+ee-e12].Partition = f1->twoBody.num;
            }
        }
        
        fromBeginning(f1,interactionEwald,interactionExchangeB);
        f1->sinc.tulip[interactionEwald].Partition = f1->twoBody.num * (c1->rt.runFlag > 0 );
        f1->sinc.tulip[interactionEwald].species = matrix;
        //f1->sinc.tulip[interactionEwald].spinor = cmpl;
        assignParticle(f1, interactionEwald, electron, two);

    {
        enum block ee ;
        for (ee = e12 ; ee <= e34 ; ee++){
            f1->sinc.tulip[interaction12Ewald+ee-e12].name = interactionEwald;
            f1->sinc.tulip[interaction12Ewald+ee-e12].species = matrix;
            if ( ee < e34 )
                f1->sinc.tulip[interaction12Ewald+ee-e12].linkNext = interaction12Ewald+ee-e12+1;
            for ( space = 0; space < SPACE ; space++)
                if ( f1->sinc.rose[space].body >= two && f1->sinc.rose[space].particle == electron )
                    f1->sinc.tulip[interaction12Ewald+ee-e12].space[space].block = ee;
            f1->sinc.tulip[interaction12Ewald+ee-e12].Partition = f1->twoBody.num;
        }
    }
        fromBeginning(f1,shortenEwald,interactionEwald);
        f1->sinc.tulip[shortenEwald].Partition = f1->twoBody.num*2*c1->i.qFloor * c1->i.decomposeRankMatrix;
        f1->sinc.tulip[shortenEwald].species = matrix;
        assignOneWithPointers(f1, shortenEwald, electron);

        fromBeginning(f1,shortenPlus,shortenEwald);
        f1->sinc.tulip[shortenPlus].Partition = f1->twoBody.num*c1->i.decomposeRankMatrix*( c1->rt.calcType == clampProtonElectronCalculation );
        f1->sinc.tulip[shortenPlus].species = matrix;
        assignOneWithPointers(f1, shortenPlus, all);
        
        fromBeginning(f1,shortenMinus,shortenPlus);
        f1->sinc.tulip[shortenMinus].Partition = f1->twoBody.num*c1->i.decomposeRankMatrix*( c1->rt.calcType == clampProtonElectronCalculation );
        f1->sinc.tulip[shortenMinus].species = matrix;
        assignOneWithPointers(f1, shortenMinus, all);
        
        fromBeginning(f1,interactionExchangePlus,shortenMinus);
        f1->sinc.tulip[interactionExchangePlus].Partition =0* mx1len*mx1len*f1->twoBody.num*( c1->rt.calcType == clampProtonElectronCalculation );
        f1->sinc.tulip[interactionExchangePlus].species = matrix;
        assignOneWithPointers(f1, interactionExchangePlus,all);
        
        fromBeginning(f1,interactionExchangeMinus,interactionExchangePlus);
        f1->sinc.tulip[interactionExchangeMinus].Partition = 0*mx1len*mx1len*f1->twoBody.num*( c1->rt.calcType == clampProtonElectronCalculation );
        f1->sinc.tulip[interactionExchangeMinus].species = matrix;
        assignOneWithPointers(f1, interactionExchangeMinus,all);
        
        fromBeginning(f1,interactionTwoAcrossDimensions,interactionExchangeMinus);
        f1->sinc.tulip[interactionTwoAcrossDimensions].Partition = (!c1->i.decomposeRankMatrix)*mx1len*mx1len*f1->twoBody.num*( c1->rt.calcType == protonsElectronsCalculation );
        f1->sinc.tulip[interactionTwoAcrossDimensions].species = matrix;
        assignParticle(f1, interactionTwoAcrossDimensions, all, one);

        {
            enum block ee,eeB ;
            for (ee = tv1 ; ee <= tv4 ; ee++)
                for (eeB = tv1 ; eeB <= tv4 ; eeB++){
                    
                f1->sinc.tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].name = interactionTwoAcrossDimensions;
                f1->sinc.tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].species = matrix;
                if ( ee < e34 && eeB < e34 )
                    f1->sinc.tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].linkNext = interactionTwoAcrossDimensions+ee-tv1+eeB-tv1+1;
                for ( space = 0; space < SPACE ; space++)
                    if ( space < 3)
                        f1->sinc.tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].space[space].block = ee;
                    else if ( space < 6 )
                        f1->sinc.tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].space[space].block = eeB;

                f1->sinc.tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].Partition = f1->twoBody.num;
            }
        }
        
        fromBeginning(f1,shortTwoAcrossDimensions,interactionTwoAcrossDimensions);
        f1->sinc.tulip[shortTwoAcrossDimensions].Partition = f1->twoBody.num*c1->i.decomposeRankMatrix*f1->twoBody.num*( c1->rt.calcType == protonsElectronsCalculation );
        f1->sinc.tulip[shortTwoAcrossDimensions].species = matrix;
        assignParticle(f1, shortTwoAcrossDimensions, all, one);

        {
            enum block ee,eeB ;
            for (ee = tv1 ; ee <= tv4 ; ee++)
                for (eeB = tv1 ; eeB <= tv4 ; eeB++){
                    
                    f1->sinc.tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].name = shortTwoAcrossDimensions;
                    f1->sinc.tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].species = matrix;
                    if ( ee < e34 && eeB < e34 )
                        f1->sinc.tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].linkNext = shortTwoAcrossDimensions+ee-tv1+eeB-tv1+1;
                    for ( space = 0; space < SPACE ; space++)
                        if ( space < 3)
                            f1->sinc.tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].space[space].block = ee;
                        else if ( space < 6 )
                            f1->sinc.tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].space[space].block = eeB;

                    f1->sinc.tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].Partition = 0* c1->i.decomposeRankMatrix;
                }
        }

        
        fromBeginning(f1,quadCube,shortTwoAcrossDimensions);
        f1->sinc.tulip[quadCube].Partition =1;
        f1->sinc.tulip[quadCube].species = matrix;
        assignParticle(f1, quadCube, all, two);
    
        fromBeginning(f1,oneArray,quadCube);
        f1->sinc.tulip[oneArray].Partition = c1->i.M1*3+c1->i.M1*c1->i.M1;
        f1->sinc.tulip[oneArray].memory = bufferAllocation;

        fromBeginning(f1,threeArray,oneArray);
        f1->sinc.tulip[threeArray].Partition = 2*c1->i.M1*c1->i.M1*c1->i.M1;
        f1->sinc.tulip[threeArray].memory = bufferAllocation;

        fromBeginning(f1,oneBasis,threeArray);
        f1->sinc.tulip[oneBasis].Partition = c1->i.M1*N1;
        f1->sinc.tulip[oneBasis].memory = bufferAllocation;

        INT_TYPE vecLen = 1;
        for ( space = 0  ; space < SPACE ; space++)
            if ( vecLen < alloc(f1, eigenVectors, space))
                vecLen = alloc(f1, eigenVectors, space);
        
        
        fromBeginning(f1,tensorBuffers,oneBasis);
        f1->sinc.tulip[tensorBuffers].Partition = vecLen ;
        f1->sinc.tulip[tensorBuffers].spinor = parallel;
        f1->sinc.tulip[tensorBuffers].memory = bufferAllocation;

        fromBeginning(f1,tensorBuffers2,tensorBuffers);
        f1->sinc.tulip[tensorBuffers2].Partition = vecLen;
        f1->sinc.tulip[tensorBuffers2].spinor = parallel;
        f1->sinc.tulip[tensorBuffers2].memory = bufferAllocation;

        fromBeginning(f1,tensorBuffers3,tensorBuffers2);
        f1->sinc.tulip[tensorBuffers3].Partition = vecLen;
        f1->sinc.tulip[tensorBuffers3].spinor = parallel;
        f1->sinc.tulip[tensorBuffers3].memory = bufferAllocation;

        fromBeginning(f1,tensorBuffers4,tensorBuffers3);
        f1->sinc.tulip[tensorBuffers4].Partition = vecLen;
        f1->sinc.tulip[tensorBuffers4].spinor = parallel;
        f1->sinc.tulip[tensorBuffers4].memory = bufferAllocation;

        fromBeginning(f1,tensorBuffers5,tensorBuffers4);
        f1->sinc.tulip[tensorBuffers5].Partition = vecLen;
        f1->sinc.tulip[tensorBuffers5].spinor = parallel;
        f1->sinc.tulip[tensorBuffers5].memory = bufferAllocation;

        fromBeginning(f1,tensorBuffers6,tensorBuffers5);
        f1->sinc.tulip[tensorBuffers6].Partition = vecLen;
        f1->sinc.tulip[tensorBuffers6].spinor = parallel;
        f1->sinc.tulip[tensorBuffers6].memory = bufferAllocation;
        
        fromBeginning(f1,oneByOneBuffer,tensorBuffers6);
        f1->sinc.tulip[oneByOneBuffer].Partition = mx1len*mx1len*mx1len*mx1len* (c1->rt.calcType >= clampProtonElectronCalculation);
        f1->sinc.tulip[oneByOneBuffer].memory = bufferAllocation;
        
        fromBeginning(f1,canonicalBuffersB,oneByOneBuffer);
        f1->sinc.tulip[canonicalBuffersB].Partition = maxVector;
        f1->sinc.tulip[canonicalBuffersB].spinor = parallel;
        f1->sinc.tulip[canonicalBuffersB].memory = bufferAllocation;

        fromBeginning(f1,canonicalBuffersBM,canonicalBuffersB);//twobody
        f1->sinc.tulip[canonicalBuffersBM].Partition = mx1len*mx1len*mx1len*mx1len ;
        f1->sinc.tulip[canonicalBuffersBM].memory = bufferAllocation;

        fromBeginning(f1,canonicalBuffersC,canonicalBuffersBM);
        f1->sinc.tulip[canonicalBuffersC].Partition = (c1->rt.phaseType == frameDensity ) *   NV;
        f1->sinc.tulip[canonicalBuffersC].spinor = parallel;
        f1->sinc.tulip[canonicalBuffersC].memory = bufferAllocation;

        fromBeginning(f1,twoBodyRitz,canonicalBuffersC);
        f1->sinc.tulip[twoBodyRitz].Partition = maxArray;
        f1->sinc.tulip[twoBodyRitz].memory = bufferAllocation;

        fromBeginning(f1,conditionOverlapNumbers,twoBodyRitz);
        f1->sinc.tulip[conditionOverlapNumbers].Partition = maxArray;
        f1->sinc.tulip[conditionOverlapNumbers].memory = bufferAllocation;
        
        fromBeginning(f1,matrixHbuild,conditionOverlapNumbers);
    f1->sinc.tulip[matrixHbuild].Partition = ( c1->i.OCSBflag || c1->rt.phaseType == solveRitz|| c1->rt.phaseType == frameDensity ) *  imax(2*maxArray*maxArray, c1->i.sectors * 2*mxlen*mxlen);
        f1->sinc.tulip[matrixHbuild].memory = bufferAllocation;

        fromBeginning(f1,vectorHbuild,matrixHbuild);
        f1->sinc.tulip[vectorHbuild].Partition = 0;
        f1->sinc.tulip[vectorHbuild].memory = bufferAllocation;

        fromBeginning(f1,matrixSbuild,vectorHbuild);
        f1->sinc.tulip[matrixSbuild].Partition = 2*maxArray*maxArray;
        f1->sinc.tulip[matrixSbuild].memory = bufferAllocation;

        
        fromBeginning(f1,dsyBuffers,matrixSbuild);
#if 1
        f1->sinc.tulip[dsyBuffers].Partition = 2*8*(8*(imax(mxlen,maxEV))+72*c1->i.nStates*c1->i.nStates+ 8 * mxlen)+3*maxEV;
#else
        f1->sinc.tulip[dsyBuffers].Partition = maxVector*maxVector;
#endif
        f1->sinc.tulip[dsyBuffers].spinor = parallel;
        f1->sinc.tulip[dsyBuffers].memory = bufferAllocation;

        fromBeginning(f1,end,dsyBuffers);
        fromBeginning(f1,end,end);
        struct name_label e = f1->sinc.tulip[end];
        struct name_label k1 = f1->sinc.tulip[kinetic1];
        struct name_label k2 = f1->sinc.tulip[kinetic2];
        struct name_label it = f1->sinc.tulip[interaction12];
        struct name_label i = f1->sinc.tulip[interactionExchangePlus];
        struct name_label ip1 = f1->sinc.tulip[interaction1Plus];


        {
            printf("\t| SPACE \t:   Gb\t \n");
            double maxMem = 0.,currMem;
            for ( space = 0 ; space <= SPACE ; space++){
                currMem = (f1->sinc.tulip[end].space[space].Address)/(1000000000./(sizeof(Stream_Type)));
                
                
                if ( space < SPACE )
                    printf("\t| %d \t\t: \t%1.9f\n",space,currMem);
                else
                    printf("\t| my \t\t: \t%1.9f\n",currMem);

                fflush(stdout);

                maxMem += currMem;
            }
            if ( maxMem > c1->i.RAMmax ){
                printf("oops too much RAM required\n");
                fflush(stdout);

                exit(0);
            }
        }
        printf("\n\n\n");
        
        for ( space = 0; space <= SPACE ; space++){
            f1->sinc.rose[space].stream = malloc( (f1->sinc.tulip[end].space[space].Address)*sizeof(Stream_Type));
        }

        c1->mem.bootedMemory = 1;
        assignCores(f1, 2);

        INT_TYPE RdsSize;
        RdsSize = 0;

        
        if (c1->i.M1){
            INT_TYPE i,ii;
            INT_TYPE M1 = c1->i.M1,M12 = (M1-1)/2;
            double r = (double)(M1-1)/(double)(N1-1);
            double * u =myStreams(f1,oneBasis,0);
            for( i = 0; i < N1 ; i++)
                for ( ii = 0 ; ii < M1 ; ii++){
                    u[i*M1+ii] = Sinc(r, r*(i-N12)- (ii-M12))/sqrt(r);
                }
        }
        
    
        
        if(c1->rt.boot == fullMatrices){
            f1->sinc.tulip[Ha].species    = scalar;
            f1->sinc.tulip[Iterator].species    = scalar;
            if ( c1->rt.calcType == electronicStuctureCalculation  ){

                if ( f1->sinc.rose[0].component == periodicComponent1 ){
                    printf("func %d\n",f1->twoBody.func.fn);
                    if ( f1->twoBody.func.fn != nullFunction)
                    if ( ! ioStoreMatrix(&c1->i.c, shortenEwald, 0, "shortenEwald.matrix", 1)){
                        tClear(f1,shortenEwald);
                    }
                    if ( f1->twoBody.func.fn != nullFunction || c1->rt.phaseType == frameDensity)
                    if ( !  ioStoreMatrix(&c1->i.c, interactionEwald, 0, "interactionEwald.matrix",1) ){
                        if ( c1->rt.phaseType != solveRitz ){
                            printf("you can remove this barrier, but I would recommend you run ritz first\n");
                            exit(0);
                        }
                        mySeparateExactTwo(f1, interactionEwald, 1., 0, 1, electron);
                    }
                    tZeroSum(f1, interactionEwald, 0);
                    zero(f1, linear, 0);
                    if ( TEST4 )
                        buildElectronProtonInteraction(f1, linear, 0);
                    if ( f1->twoBody.func.fn != nullFunction|| c1->rt.phaseType == frameDensity)
                    if ( ! ioStoreMatrix(&c1->i.c, interactionExchange, 0, "interactionExchange.matrix",1) && c1->rt.phaseType == solveRitz ){
                        printf("skipping zero-cell intrainteraction\n");
                        if (0)
                            mySeparateExactTwo(f1, interactionExchange, 1., 0, 0, electron);
                    }
                    separateExternal(c1,linear,0,0,-1.0,-1,0,electron);
                }else {
//                    if ( f1->twoBody.func.fn != nullFunction)
                    if ( bootBodies > one ){
                        if ( ! ioStoreMatrix(&c1->i.c, interactionExchange, 0, "interactionExchange.matrix",1)){
                            mySeparateExactTwo(f1, interactionExchange, 1., 0, 0, electron);
                        }
                    }
                    
                    separateExternal(c1,linear,0,0,-1.0,-1,0,electron);
                    
                }
                separateKinetic(f1, 0,kinetic, c1->i.massElectron,electron);

                if ( c1->i.magFlag ){
                    INT_TYPE deriv[SPACE];
                    INT_TYPE power[SPACE];
                    
                    
                    deriv[0] = 1;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 0 ;
                    power[1] = 1;
                    power[2] = 0;
                    
                    separateDerivatives(f1, 0, vectorMomentum, power, deriv, 0.5*c1->i.mag, electron);
                    
                    deriv[0] = 0;
                    deriv[1] = 1;
                    deriv[2] = 0;
                    
                    power[0] = 1;
                    power[1] = 0;
                    power[2] = 0;
                    
                    separateDerivatives(f1,0, vectorMomentum, power, deriv, -0.5*c1->i.mag, electron);

                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 2;
                    power[1] = 0;
                    power[2] = 0;
                    
                    separateDerivatives(f1, 0, vectorMomentum, power, deriv, 0.25*c1->i.mag, electron);

                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 0;
                    power[1] = 2;
                    power[2] = 0;
                    
                    separateDerivatives(f1, 0, vectorMomentum, power, deriv, 0.25*c1->i.mag, electron);
                    

                }
                
                if ( c1->i.springFlag ){
                    INT_TYPE deriv[SPACE];
                    INT_TYPE power[SPACE];
                    
                    
                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 2 ;
                    power[1] = 0;
                    power[2] = 0;
                    
                    separateDerivatives(f1, 0, vectorMomentum, power, deriv, 0.5*c1->i.springConstant, electron);
                    
                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 0;
                    power[1] = 2;
                    power[2] = 0;
                    
                    separateDerivatives(f1,0, vectorMomentum, power, deriv,  0.5*c1->i.springConstant, electron);
                    
                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 0;
                    power[1] = 0;
                    power[2] = 2;
                    
                    separateDerivatives(f1, 0, vectorMomentum, power, deriv,  0.5*c1->i.springConstant, electron);
                    
                
                }
                

        } else if ( c1->rt.calcType == clampProtonElectronCalculation  ){
                if ( bootBodies > one ){
                    mySeparateExactTwo(f1,interactionExchange, 1. , 0,0,electron);
                }
            
                if ( !ioStoreMatrix(f1,shortenPlus ,0,"shortenExchangePlus.matrix",1) )
                    mySeparateExactOneByOne(f1, c1->i.decomposeRankMatrix,interactionExchangePlus, shortenPlus,-1., 1,c1->i.massClampPair/(c1->i.massProton + c1->i.massClampPair),electron, proton);
            
                if ( ! ioStoreMatrix(f1,shortenMinus ,0,"shortenExchangeMinus.matrix",1) )
                    mySeparateExactOneByOne(f1, c1->i.decomposeRankMatrix,interactionExchangeMinus, shortenMinus,-1., -1,c1->i.massProton/(c1->i.massProton + c1->i.massClampPair),electron, proton);

            struct name_label u = f1->sinc.tulip[shortenPlus];
            
            separateExternal(c1,protonRepulsion, 0,0,1.0,-1,0,proton);
            separateKinetic(f1, 0,kineticMass, c1->i.massProton*c1->i.massClampPair/(c1->i.massProton+c1->i.massClampPair),proton);

            separateKinetic(f1, 0,kinetic, c1->i.massElectron*(c1->i.massProton + c1->i.massClampPair )/(c1->i.massProton + c1->i.massClampPair + bootBodies*c1->i.massElectron ),electron);

            }
            else if ( c1->rt.calcType == protonsElectronsCalculation  ){

                if ( bootBodies > one ){
                    mySeparateExactTwo(f1,interactionExchange, 1. , 0,0,electron);
                    mySeparateExactTwo(f1,interactionExchangeB, 1. , 0,0,proton);
                //    tZeroSum(f1, interactionExchange, 0 );
                //    tZeroSum(f1, interactionExchangeB, 0 );
                }
                mySeparateExactOneByOne(f1,c1->i.decomposeRankMatrix,interactionTwoAcrossDimensions, shortTwoAcrossDimensions,-1., 1,1,electron, proton);
                separateKinetic(f1, 0,kineticMass, c1->i.massProton,proton);
                separateKinetic(f1, 0,kinetic, c1->i.massElectron,electron);

            }


        
    
        {
            if ( c1->i.vectorOperatorFlag ){

                f1->sinc.tulip[Ha].linkNext = vectorOperator;
                f1->sinc.tulip[Iterator].linkNext = vectorOperator;
                INT_TYPE fi,fe=0,ff;
                for ( fi =0 ; fi < c1->mem.filesVectorOperator ; fi++){
                    fe +=  tLoadEigenWeights (c1, c1->mem.fileVectorOperator[fi], f1->sinc.vectorOperator+fe);
                }
                
                for ( ff = 1; ff < fe ; ff++)
                    f1->sinc.tulip[vectorOperator+ff-1].linkNext = vectorOperator+ff;

            
                
                if ( fe > c1->i.vectorOperatorFlag )
                {
                    printf("failure to load vector Operators");
                    exit(0);
                }else {
                    c1->i.vectorOperatorFlag = fe;
                }
                
            }
            else
            if ( bootBodies == one && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->sinc.tulip[shorten1Ewald].linkNext = vectorMomentum1;
                f1->sinc.tulip[vectorMomentum1].linkNext = kinetic1;
                f1->sinc.tulip[kinetic1].linkNext = external1;
                f1->sinc.tulip[external1].linkNext = nullName;

                //active assignment
                f1->sinc.tulip[Ha].linkNext = shorten1Ewald;
                f1->sinc.tulip[Iterator].linkNext = shorten1Ewald;
                
            } else if ( bootBodies == two && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->sinc.tulip[shorten2Ewald].linkNext = vectorMomentum1;
                f1->sinc.tulip[vectorMomentum2].linkNext = kinetic1;
                if ( c1->rt.runFlag == 0 ){
                    f1->sinc.tulip[kinetic2].linkNext = external1;
                    f1->sinc.tulip[external2].linkNext = interaction12;
                    f1->sinc.tulip[interaction12].linkNext = nullName;

                }else {
                    f1->sinc.tulip[kinetic2].linkNext = external1;
                    f1->sinc.tulip[external2].linkNext = interaction12Ewald;
                    f1->sinc.tulip[interaction12Ewald].linkNext = nullName;
                }
                //active assignment
                f1->sinc.tulip[Ha].linkNext = shorten1Ewald;
                f1->sinc.tulip[Iterator].linkNext = shorten1Ewald;
                
            }else if ( bootBodies == three && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->sinc.tulip[shorten3Ewald].linkNext = vectorMomentum1;
                f1->sinc.tulip[vectorMomentum3].linkNext = kinetic1;

                if ( c1->rt.runFlag == 0 ){
                    f1->sinc.tulip[kinetic3].linkNext = external1;
                    f1->sinc.tulip[external3].linkNext = interaction12;
                    f1->sinc.tulip[interaction23].linkNext = nullName;
                }else {
                    f1->sinc.tulip[kinetic3].linkNext = external1;
                    f1->sinc.tulip[external3].linkNext = interaction12Ewald;
                    f1->sinc.tulip[interaction23Ewald].linkNext = nullName;
                }
                //active assignment
                    f1->sinc.tulip[Ha].linkNext = shorten1Ewald;
                    f1->sinc.tulip[Iterator].linkNext = shorten1Ewald;
                
            }
            
            else if ( bootBodies == four && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->sinc.tulip[shorten4Ewald].linkNext = vectorMomentum1;
                f1->sinc.tulip[vectorMomentum4].linkNext = kinetic1;

                if ( c1->rt.runFlag == 0 ){
                    f1->sinc.tulip[kinetic4].linkNext = external1;
                    f1->sinc.tulip[external4].linkNext = interaction12;
                    f1->sinc.tulip[interaction34].linkNext = nullName;
                }else {
                    f1->sinc.tulip[kinetic4].linkNext = external1;
                    f1->sinc.tulip[external4].linkNext = interaction12Ewald;
                    f1->sinc.tulip[interaction34Ewald].linkNext = nullName;
                }

                //active assignment
                f1->sinc.tulip[Ha].linkNext = shorten1Ewald;
                f1->sinc.tulip[Iterator].linkNext = shorten1Ewald;
                
            }
            else if (bootBodies == one && ( c1->rt.calcType == clampProtonElectronCalculation  )){
                struct name_label ip1 = f1->sinc.tulip[interaction1Plus];
                
                
                //paths do not matter unless assigned below...
                f1->sinc.tulip[interaction1Plus].linkNext =interaction1Minus;
                f1->sinc.tulip[interaction1Minus].linkNext =kinetic1;
                
                f1->sinc.tulip[shorten1Plus].linkNext =shorten1Minus;
                f1->sinc.tulip[shorten1Minus].linkNext =kinetic1;
                
                
                f1->sinc.tulip[kinetic1].linkNext = kineticMass1;
                f1->sinc.tulip[kineticMass1].linkNext = proton1;
                f1->sinc.tulip[proton1].linkNext = nullName;
                //active assignment
                if ( c1->i.decomposeRankMatrix ){
                    f1->sinc.tulip[Ha].linkNext = shorten1Plus;
                    f1->sinc.tulip[Iterator].linkNext = shorten1Plus;
                } else {
                    f1->sinc.tulip[Ha].linkNext = interaction1Plus;
                    f1->sinc.tulip[Iterator].linkNext = interaction1Plus;
                    
                }
            }
            else if (bootBodies == two && ( c1->rt.calcType == clampProtonElectronCalculation  )){
                struct name_label ip1 = f1->sinc.tulip[interaction1Plus];
                if ( c1->i.decomposeRankMatrix )
                    f1->sinc.tulip[interaction12].linkNext =shorten1Plus;
                else
                    f1->sinc.tulip[interaction12].linkNext =interaction1Plus;

                f1->sinc.tulip[interaction2Plus].linkNext =interaction1Minus;
                f1->sinc.tulip[interaction2Minus].linkNext =kinetic1;

                
                f1->sinc.tulip[shorten2Plus].linkNext =shorten1Minus;
                f1->sinc.tulip[shorten2Minus].linkNext =kinetic1;

                f1->sinc.tulip[kinetic1].linkNext = kinetic2;
                f1->sinc.tulip[kinetic2].linkNext = kineticMass1;
                f1->sinc.tulip[kineticMass1].linkNext = proton1;
                
                
                
                //active assignment
                f1->sinc.tulip[Ha].linkNext = interaction12;
                f1->sinc.tulip[Iterator].linkNext = interaction12;
                
            }
            else if (bootBodies == two && ( c1->rt.calcType == protonsElectronsCalculation  )){
                struct name_label ip1 = f1->sinc.tulip[interaction1Plus];
                f1->sinc.tulip[interaction12].linkNext =interaction12B;

                if ( c1->i.decomposeRankMatrix )
                    f1->sinc.tulip[interaction12B].linkNext =shortTAD11;
                else
                    f1->sinc.tulip[interaction12B].linkNext =interactionTAD11;
                
                f1->sinc.tulip[interactionTAD22].linkNext =kinetic1;
                f1->sinc.tulip[shortTAD22].linkNext =kinetic1;
                
                f1->sinc.tulip[kinetic2].linkNext = kineticMass1;
                f1->sinc.tulip[kineticMass2].linkNext = nullName;
                
                //active assignment
                f1->sinc.tulip[Ha].linkNext = interaction12;
                f1->sinc.tulip[Iterator].linkNext = interaction12;
                
            }

            
            
        }
        {
            INT_TYPE nn[SPACE],i,j;
            length1(f1,nn);
            for ( space = 0; space < SPACE ; space++)
                if ( f1->sinc.rose[space].body != nada )
                    for ( i = 0 ; i < nn[space];i++)
                        for ( j = 0; j < nn[space];j++)
                        {
                            if ( i+j+1 == nn[space] )
                                streams(f1, inversion,0,space)[i*nn[space]+j] = 1.;
                            else
                                streams(f1, inversion,0,space)[i*nn[space]+j] = 0.;
                            
                        }
            
        }
        }
    
   // tInnerTest(f1, kinetic, copy);
#if VERBOSE
    printf("boot complete\n");
#endif
    fflush(stdout);
    return 0;
}
