/**
 *  Model.c
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

#include "Model.h"


/**
 *initialize field parameters
 *
*/
  field initField (void ) {
    inta space;
      field i;
    i.i.body = nada;
    i.i.nStates = 0;
    i.i.qFloor = 0;
    i.i.nOperator = 0;
    i.i.OpIndex = 0;
    i.i.canonRank = 0;
    i.i.iRank = 0;
    i.i.xRank = 0;
    i.i.filter = 0;
    i.i.irrep = 0;
    i.i.collect = 0;
    i.i.body = one;
    i.i.Iterations = 0;
    i.f.boot = noMatrices;
    for( space = 0 ; space <= SPACE ; space++){
        i.f.canon[space].stream = NULL;
        i.f.canon[space].body = nada;//place holder...
        i.f.canon[space].label = 0;
    }
#ifdef APPLE
    space =0;
    i.f.canon[space].basis = SincBasisElement;
    i.f.canon[space].body = two;
    i.f.canon[space].component = 1;
    i.f.canon[space].count1Basis = 9;
    i.f.canon[space].space = 0;
    i.f.canon[space].label = 1;
    i.f.canon[space].particle[one].attack = 0.5;
    i.f.canon[space].particle[one].lattice = 2;
    i.f.canon[space].particle[one].origin = -4*2 ;
    i.f.canon[space].particle[two].attack = 0.5;
    i.f.canon[space].particle[two].lattice = 2;
    i.f.canon[space].particle[two].origin = -4*2 ;

    space =1;
    i.f.canon[space].basis = SincBasisElement;
    i.f.canon[space].body = two;
    i.f.canon[space].component = 1;
    i.f.canon[space].count1Basis = 9;
    i.f.canon[space].space = 0;
    i.f.canon[space].label = 1;
    i.f.canon[space].particle[one].attack = 0.5;
    i.f.canon[space].particle[one].lattice = 2;
    i.f.canon[space].particle[one].origin = -4*2 ;
    i.f.canon[space].particle[two].attack = 0.5;
    i.f.canon[space].particle[two].lattice = 2;
    i.f.canon[space].particle[two].origin = -4*2 ;
    
      space =2;
      i.f.canon[space].basis = SincBasisElement;
      i.f.canon[space].body = two;
      i.f.canon[space].component = 1;
      i.f.canon[space].count1Basis = 9;
      i.f.canon[space].space = 0;
      i.f.canon[space].label = 1;
      i.f.canon[space].particle[one].attack = 0.5;
      i.f.canon[space].particle[one].lattice = 2;
      i.f.canon[space].particle[one].origin = -4*2 ;
      i.f.canon[space].particle[two].attack = 0.5;
      i.f.canon[space].particle[two].lattice = 2;
      i.f.canon[space].particle[two].origin = -4*2 ;

#endif
//
//
//    for ( space = 1 ;space < 3 ;space++){
//        i.f.canon[space].attack = 0.5;
//        i.f.canon[space].basis = SincBasisElement;
//        i.f.canon[space].body = one;
//        i.f.canon[space].component = 2;
//        i.f.canon[space].count1Basis = 5;
//        i.f.canon[space].label = 2;
//        i.f.canon[space].lattice = 1.;
//        i.f.canon[space].space = space-1;
//
//    }
    i.f.bootedMemory = 0;
    i.f.name = NULL;
    i.i.files = 0;
    i.i.filesVectorOperator = 0;
    i.f.nullLabels.currLabel = 0;
    i.f.eikonLabels.currLabel = 0;
    i.f.nullLabels.maxLabel = 10000;
    i.f.eikonLabels.maxLabel = 100;

#ifdef APPLE


#endif
    
    
#ifdef OVERFLAG
        i.i.cmpl = cmpl;
        i.i.canonRank = 4;
        i.i.iRank = 1;
        i.i.nStates = 1;
        i.i.qFloor = 9*9*9;
        i.f.boot = fullMatrices;
        i.i.body = one;
        i.i.irrep = 1;
#else
    i.i.OpIndex =-1;
    i.i.cmpl = real;
    i.i.canonRank = 35;
    i.i.iRank = 1;
    i.i.nStates = 1;
    i.i.qFloor = 1;
    i.f.boot  = fullMatrices;
    i.i.irrep = 2;
#endif
    return i;
}

/**
 *initialize calculation parameters
 *
 */
calculation initCal (void ) {
      calculation i;
    i.i.build = 1;
    i.i.minIterationPrint = 0;
    i.i.termNumber = 0;
    i.i.shiftVector[0] = 0;
    i.i.shiftVector[1] = 1;
    i.rt.MAX_CYCLE = 10;
    i.rt.relativeTOLERANCE = 0.000000000000001;
    i.rt.THRESHOLD = 1e-12; 
    i.i.numNames = 1000;
    i.i.numVectors = 100;
#ifdef APPLE
    resetA(&i.rt);
    blockA(&i.rt, 3);
    blockA(&i.rt, 4);
    blockA(&i.rt, 5);
    i.i.RAMmax = 6;
    i.i.Angstroms = 0;
    i.rt.ALPHA = 1e-9;
    i.rt.TOLERANCE = 1e-7;
    i.rt.XCONDITION = 1e5;
    i.i.lambda = 2;
    
    
    i.rt.calcType = electronicStuctureCalculation;
    i.rt.phaseType = buildFoundation;
    i.i.Na = 0;
    i.rt.calcType = electronicStuctureCalculation;
    i.i.lambda = 200 ;
#else
    i.i.RAMmax = 1;
    i.i.Angstroms = 0;
    i.rt.TOLERANCE = 1e-7;
    i.rt.relativeTOLERANCE = 1e-7;
    i.rt.ALPHA = 1e-8;
    i.rt.THRESHOLD = 1e-13;
    i.rt.XCONDITION = 1e5;
    i.rt.calcType = nullCalculation;
    i.rt.phaseType = buildFoundation ;
#endif
    return i;
}

/**
 * Deallocation of resources
 *
 *Seen at the end of jobs
 *@param f1     The container only
 */
inta fModel (   sinc_label * f1){
    inta i;
    
    if ( f1->bootedMemory ){
        for ( i = 0;i <= SPACE ; i++)
            if ( f1->canon[i].body != nada || i == SPACE){
            free(f1->canon[i].stream);
            f1->canon[i].stream = NULL;

            if ( i < SPACE ){

            }
        }
        {
            free(f1->name);
            f1->name = NULL;
        }
        f1->bootedMemory = 0;
    }
    return 0;
}



/**
 *The allocation of resources
 *
 *The 'blockMemory' command follows blockMemoryType here.
 *This command will block aspects of the allocation.
 *
 *@param c1     parameters associated with runtime
 *@param f      the container of all allocated memory
 */
inta iModel(   calculation * c1,   field *f){
    f->f.nullLabels.maxLabel = c1->i.numNames;
    f->f.eikonLabels.maxLabel = c1->i.numVectors;

    inta splitOperator = 1 ;
    sinc_label *f1 = &f->f;

    {//SA++
        f->f.irrep = f->i.irrep;
        f->f.cmpl = f->i.cmpl;
    }//SA++
    
    if(f1->boot == fullMatrices){
        printf("\t tolerance\t\t10^-\t%1.1f\n", -log(c1->rt.TOLERANCE)/log(10));
        printf("\t relativeTolerance\t10^-\t%1.1f\n", -log(c1->rt.relativeTOLERANCE )/log(10));
        printf("\t condition\t\t10^-\t%1.1f\n", -log(c1->rt.ALPHA )/log(10));
        printf("\t threshold\t\t10^-\t%1.1f\n", -log(c1->rt.THRESHOLD)/log(10));
        printf("\t maxCycle\t\t\t%d\n", c1->rt.MAX_CYCLE);
        printf("\t maxCondition\t\t10^\t%1.1f\n", log(c1->rt.XCONDITION)/log(10));
    }
    
    
        runTime * rt = &c1->rt;
    
        inta space;
        {
            f->f.rt = rt;
        }
        
        bodyType bootBodies = f->i.body;
    
        ///simplify the allocations
        f->i.iRank = imax(f->i.iRank ,f->i.canonRank);
        f->i.xRank = imax( f->i.xRank,f->i.iRank*f->i.qFloor  );
        inta maxVector = imax(f->i.canonRank,1+f->i.iRank);
    
    

    inta FloorV = f->i.qFloor;
        inta maxArray,EV = FloorV,NV=0 ;
    
        maxArray = imax(f->i.nStates,EV);//slip Nb into spectra...
        
        f1->maxEV = EV;
          division vectorOperator  = eigenVectors + 1 +  f->i.nStates+EV;
        f1->vectorOperator = vectorOperator;
    {
        inta ir,ix;
        f->i.nOperator = countLinesFromFile(c1, *f, 1, &ir, &ix);
        if (f->i.nOperator ){
            printf("Operators\t %d\n", f->i.nOperator);
        }
        
        inta rx = 1;
        switch (bootBodies){
            case three:
                rx = 3;
                break;
            case four:
                rx = 6;
                break;
            case five:
                rx = 10;
                break;
            case six:
                rx = 15;
                break;
            default:
                break;
        }
        if ( !splitOperator )
            rx = 0;
    
          division end = vectorOperator+f->i.nOperator*(2*rx+1)+f1->nullLabels.maxLabel + f1->eikonLabels.maxLabel;
        //EIKONS
        //EIKONS
        f1->nullLabels.head = vectorOperator+f->i.nOperator*(2*rx+1);
        //EIKONS
        f1->eikonLabels.head  = f1->nullLabels.head +f1->nullLabels.maxLabel;

    
        f1->end = end;
        f1->name = malloc ( (end+1) * sizeof(  name_label));
    }
        
        {//defaults
            //define vectors end
            {
                inta c;
                  division label1;
                for ( label1 = 0 ;label1 <= f1->end; label1++){
                    f1->name[label1].name = label1;
                    f1->name[label1].Partition = 0;
                    f1->name[label1].spinor = f1->cmpl;
                    f1->name[label1].species = scalar;
                    f1->name[label1].linkNext = nullName;
                    f1->name[label1].chainNext = nullName;
                    f1->name[label1].multId = -label1;//unique entries do not multiply up.
                    f1->name[label1].loopNext = nullName;
                    f1->name[label1].memory = objectAllocation;
                   // f1->name[label1].operatorSignFlag = 0;
                    for ( space = 0; space <= SPACE ; space++)
                        f1->name[label1].space[space].Address = -1;
                    f1->name[label1].space[SPACE].block = id0;
                    for ( space = 0; space <= SPACE ; space++){
                        {
                            f1->name[label1].space[space].block = id0;//matrix prototype
                            f1->name[label1].space[space].body = nada;//matrix prototype
                            f1->name[label1].space[space].act = 1;
                            f1->name[label1].space[space].invert = 0;
                        }
                    }
                    tClear(*f1,label1);
                    for ( c = 0 ; c < MAX_CORE ; c++)
                        f1->name[label1].Begin[c] = 0;
                }
                
            }
        }//defaults
    
    
    
        {
            division last = 0;

            if ( f1->boot==fullMatrices ){
                {
                    inta i ;
                    for ( i = 1; i < BLOCK_COUNT ; i++){
                        if ( ! allowQ(f1->rt,i))
                            printf("\t blockMemory %d\n",i);
                    }
                }
                
                if ( f->i.OpIndex ){
                inta terms = defineTerms(c1, f1,Iterator,0);
                if ( terms)
                    printf("\n\n#spam %d\n\n",terms);
                }
            }
            
            
            
            
                inta cmpl;
                if ( f->i.filesVectorOperator )
                    
                {
                      bodyType bd;
                    inta fi,lines = 0,num;
                    int ms = MAXSTRING;
                    char line0[MAXSTRING];
                    char name[MAXSTRING];
                    char title[MAXSTRING];
                    char *line = line0;
                    inta FIT ;
                    FIT = f->i.filesVectorOperator ;
                    for ( fi =0 ; fi < FIT; fi++){
                        strcpy(name ,f->i.fileVectorOperator[fi]);
                        FILE * fp = fopen(name,"r");
                        if ( fp == NULL ) {
                            printf("file?\n");
                            exit(0);
                        }
                        fgets(line, ms, fp);
                        while(!feof(fp))
                        {
                            if (! comment(line))
                                {
                                    fromBeginning(*f1, vectorOperator+lines, last);
                                    inta part1 = 0;
                                      genusType genus ;
                                    for ( cmpl = f1->cmpl-1 ; cmpl >= 0 ; cmpl--){
                                        tFromReadToFilename(NULL, line,  name, f1->cmpl-1,cmpl,title,&num);
                                        part1 = imax(part1,inputFormat(*f1, name, nullName, 2));
                                    }//name = real component here.
                                    f1->name[vectorOperator+lines].Partition = part1;
                                    genus = inputFormat(*f1, name, nullName, 0);
                                    bd = inputFormat(*f1, name, nullName, 100);
                                    printf("%d %d %s\n", bd, genus,name);
                                    if ( (bd==two && genus == vector )|| (bd == one && genus == matrix )){
                                        printf("matrix-%dbody-Op\n",1);
                                        f1->name[vectorOperator+lines].species = matrix;
                                        for ( space = 0 ; space < 3/*particular*/ ; space++)
                                            f1->name[vectorOperator+lines].space[space].body = one;
                                    }else {
                                        printf("vector*vector-Op\n");

                                        f1->name[vectorOperator+lines].species = outerVector;
                                        for (space = 0; space < SPACE ; space++)
                                        if ( f1->canon[space].body != nada)
                                        {
                                            f1->name[vectorOperator+lines].space[space].body = bd;
                                         //   booting[ bd - bootBodies ] = 1;
                                        }
                                    }
                                    last = vectorOperator+lines;
                                    lines++;
                                }
                            fgets(line, ms, fp);
                        }
                        if ( fi > MAX_FILE)
                        {
                            printf("too many files, increase MAX_FILE\n");
                            exit(0);
                        }
                        fclose(fp);
                    }
                }
                
            
            
              name_label u ;
              name_label u2;


            if(1 /* !si*/){
                inta di,d0=0;
#if VERBOSE
                printf("std USERs\n\n");
#endif
                for ( di = 0 ; di < EV+f->i.nStates; di++){
                    fromBeginning(*f1, eigenVectors+di, last);
                    last = eigenVectors+di;
                    f1->name[eigenVectors+di].spinor =f1->cmpl;
                    if ( d0 < f->i.nStates ){
                        f1->name[eigenVectors+di].Partition = f->i.canonRank;
                        d0++;
                    }
                    else if ( eigenVectors+di <= eigenVectors+f->i.nStates+EV){
                        {
                            f1->name[eigenVectors+di].Partition = ((di-d0)/EV)+f->i.iRank;
                            NV += spins(*f1,eigenVectors+di )*f1->name[eigenVectors+di].Partition;
                        }
                        
                    }else{
                        exit(1);
                    }
                    f1->name[eigenVectors+di].species = vector;

                }
                fromBeginning(*f1, diagonalVectorA, last);
                u = f1->name[eigenVectors];
                u2 = f1->name[eigenVectors+di-1];
                f1->user = eigenVectors + d0;
            }
                
            

            f1->name[diagonalVectorA].spinor = parallel;
            f1->name[diagonalVectorA].Partition =(c1->rt.phaseType == buildFoundation);
            f1->name[diagonalVectorA].species = vector;
            
            fromBeginning(*f1, diagonalVectorB, diagonalVectorA);
            f1->name[diagonalVectorB].spinor = parallel;
            f1->name[diagonalVectorB].Partition = (c1->rt.phaseType == buildFoundation);
            f1->name[diagonalVectorB].species = vector;

        }
        inta len[SPACE];
        inta mxlen=0;
        inta mx1len=0;
        length1(*f1, len);
        for ( space = 0 ; space < SPACE ; space++)
            if ( f1->canon[space].body != nada ){
                len[space]+=f1->canon[space].count1Inc;
                if(len[space] > mx1len)
                    mx1len = len[space];
            }
        mxlen = Power(mx1len, bootBodies);
    
        fromBeginning(*f1,productVector,diagonalVectorB);
        f1->name[productVector].Partition =  0;
        f1->name[productVector].species = vector;
        f1->name[productVector].spinor = parallel;

    
        fromBeginning(*f1,scalarTemp,productVector);
        f1->name[scalarTemp].Partition =  ( f1->cmpl == 2 )*  maxVector*2;
        f1->name[scalarTemp].species = vector;
        f1->name[scalarTemp].spinor = real;

        fromBeginning(*f1,permutationVector,scalarTemp);
        f1->name[permutationVector].Partition =  allowQ(f1->rt, blockPermutationsblock);
        f1->name[permutationVector].species = vector;
        if ( allowQ(f1->rt, blockParallelPermuteblock) )
            f1->name[permutationVector].spinor = parallel;
        
        fromBeginning(*f1,permutation2Vector,permutationVector);
        f1->name[permutation2Vector].Partition = allowQ(f1->rt, blockPermutationsblock);
        f1->name[permutation2Vector].species = vector;
        if ( allowQ(f1->rt, blockParallelPermuteblock) )
            f1->name[permutation2Vector].spinor = parallel;
    
    
        fromBeginning(*f1,multiplyVector,permutation2Vector);
        f1->name[multiplyVector].Partition = 1;
        f1->name[multiplyVector].species = vector;
        if ( allowQ(f1->rt, blockParallelMultiplyblock) )
            f1->name[multiplyVector].spinor = parallel;

    
        fromBeginning(*f1,canonicalmvVector,multiplyVector);
        f1->name[canonicalmvVector].Partition = 1*0;
        f1->name[canonicalmvVector].species = vector;
        f1->name[canonicalmvVector].spinor = parallel;

        fromBeginning(*f1,canonicalmv2Vector,canonicalmvVector);
        f1->name[canonicalmv2Vector].Partition = 1*0;
        f1->name[canonicalmv2Vector].species = vector;

        fromBeginning(*f1,canonicalmv3Vector,canonicalmv2Vector);
        f1->name[canonicalmv3Vector].Partition = 2;
        f1->name[canonicalmv3Vector].species = vector;
        if ( allowQ(f1->rt, blockParallelMultiplyblock) )
            f1->name[canonicalmv3Vector].spinor = parallel;

        fromBeginning(*f1,canonicaldotVector,canonicalmv3Vector);
        f1->name[canonicaldotVector].Partition = allowQ(f1->rt, blockPermutationsblock);
        f1->name[canonicaldotVector].species = vector;;
        if ( allowQ(f1->rt, blockParallelPermuteblock) )
            f1->name[canonicaldotVector].spinor = parallel;
        
        fromBeginning(*f1,canonicaldot2Vector,canonicaldotVector);
        f1->name[canonicaldot2Vector].Partition = allowQ(f1->rt, blockPermutationsblock);
        f1->name[canonicaldot2Vector].species = vector;;
        if ( allowQ(f1->rt, blockParallelPermuteblock) )
            f1->name[canonicaldot2Vector].spinor = parallel;
        
        fromBeginning(*f1,canonicaldot3Vector,canonicaldot2Vector);
        f1->name[canonicaldot3Vector].Partition = 1*0;
        f1->name[canonicaldot3Vector].species = vector;;
        f1->name[canonicaldot3Vector].spinor = parallel;
        
        fromBeginning(*f1,canonicalvvVector,canonicaldot3Vector);
        f1->name[canonicalvvVector].Partition = 1*0;
        f1->name[canonicalvvVector].species = vector;;
        f1->name[canonicalvvVector].spinor = parallel;
        
        fromBeginning(*f1,canonicalvv2Vector,canonicalvvVector);
        f1->name[canonicalvv2Vector].Partition = 1*0;
        f1->name[canonicalvv2Vector].species = vector;;
        f1->name[canonicalvv2Vector].spinor = parallel;
        
        fromBeginning(*f1,canonicalvv3Vector,canonicalvv2Vector);
        f1->name[canonicalvv3Vector].Partition = 1*0;
        f1->name[canonicalvv3Vector].species = vector;;
        f1->name[canonicalvv3Vector].spinor = parallel;

        fromBeginning(*f1,canonicalmeVector,canonicalvv3Vector);
        f1->name[canonicalmeVector].Partition = 1;
        f1->name[canonicalmeVector].species = vector;;
        if (allowQ(f1->rt, blockParallelMatrixElementblock))
            f1->name[canonicalmeVector].spinor = parallel;
        
        fromBeginning(*f1,copyVector,canonicalmeVector);
        if ( allowQ(&c1->rt, blockCopyBlock))
            f1->name[copyVector].Partition = maxVector;
        f1->name[copyVector].species = vector;

        fromBeginning(*f1,copyTwoVector,copyVector);
        if ( allowQ(&c1->rt, blockCopyBlock))
            f1->name[copyTwoVector].Partition =   maxVector;
        f1->name[copyTwoVector].species = vector;

        fromBeginning(*f1,copyThreeVector,copyTwoVector);
        if ( allowQ(&c1->rt, blockCopyBlock))
            f1->name[copyThreeVector].Partition =   maxVector;
        f1->name[copyThreeVector].species = vector;

        fromBeginning(*f1,copyFourVector,copyThreeVector);
        if ( allowQ(&c1->rt, blockCopyBlock))
            f1->name[copyFourVector].Partition =    maxVector;
        f1->name[copyFourVector].species = vector;
        
        fromBeginning(*f1,totalVector,copyFourVector);
  
    
    {
        inta ra = 1;
        if ( f->i.filter && f->i.irrep )
        switch (bootBodies){
            case two:
                ra = 2;
                break;
            case three:
                ra = 6;
                break;
            case four:
                ra = 24;
                break;
            default:
                break;
    }
    
    
    if ( allowQ(f1->rt,blockTotalVectorBlock) ){
            f1->name[totalVector].Partition =  imax( imax(f->i.xRank, c1->i.lambda * maxVector), ra * maxVector) ;
    }
    }
        f1->name[totalVector].species = vector;
                    
        fromBeginning(*f1,bandBasis,totalVector);
        f1->name[bandBasis].Partition = allowQ(f1->rt,blockTransferBasisblock)*(4*mx1len*mx1len+2*mxlen);
        f1->name[bandBasis].memory = bufferAllocation;

        inta maxOriginRank =  part(*f1,totalVector);
        inta maxTrainRank  = maxVector;
        inta flag = allowQ(f1->rt, blockTrainVectorsblock) ;
        fromBeginning(*f1,canonicalBuffers,bandBasis);
        f1->name[canonicalBuffers].Partition = flag*(maxTrainRank*maxTrainRank+ maxOriginRank*maxTrainRank+maxTrainRank) ;
        f1->name[canonicalBuffers].species = scalar;
    
        fromBeginning(*f1,canonicalBuffersAster,canonicalBuffers);
        f1->name[canonicalBuffersAster].Partition = flag*(2*maxOriginRank+ 2*maxOriginRank*maxOriginRank) ;
        f1->name[canonicalBuffersAster].memory = bufferAllocation;
        f1->name[canonicalBuffersAster].species = scalar;
    
        fromBeginning(*f1,trackBuffer,canonicalBuffersAster);
        f1->name[trackBuffer].Partition = 2*flag*(maxTrainRank*maxTrainRank);
        f1->name[trackBuffer].memory = bufferAllocation;
        f1->name[trackBuffer].species = scalar;

        fromBeginning(*f1,guideBuffer,trackBuffer);
        f1->name[guideBuffer].Partition = flag*(maxOriginRank*maxTrainRank);
        f1->name[guideBuffer].memory = bufferAllocation;
        f1->name[guideBuffer].species = scalar;

        fromBeginning(*f1,canonicalBuffersB,guideBuffer);
        f1->name[canonicalBuffersB].Partition = (allowQ(f1->rt, blockTrainVectorsblock))* (maxTrainRank+maxOriginRank);
        f1->name[canonicalBuffersB].spinor = parallel;
        f1->name[canonicalBuffersB].memory = bufferAllocation;
        f1->name[canonicalBuffersB].species = scalar;

        fromBeginning(*f1,canonicalVector,canonicalBuffersB);
        f1->name[canonicalVector].Partition = 1;
        f1->name[canonicalVector].spinor = parallel;
        if ( allowQ(f1->rt, blockTrainVectorsblock))
            f1->name[canonicalVector].species = vector;
        assignParticle(*f1, canonicalVector, all, two);
    
       ///EIKON
    fromBeginning(*f1,eikonBuffer,canonicalVector);
    f1->name[eikonBuffer].Partition = !(!f1->eikonLabels.maxLabel);
    f1->name[eikonBuffer].species = eikon;
    assignParticle(*f1, eikonBuffer, all, two);
        
    
    {
        inta ii;
          division prev = eikonBuffer;
        for ( ii = 0 ; ii < f1->eikonLabels.maxLabel ; ii++)
            {
                  division label1 = f1->eikonLabels.head+ii;
                f1->name[label1].Partition = 1;
                assignParticle(*f1, label1, 0, two);
                f1->name[label1].species = eikon;
                fromBeginning(*f1,label1,prev);
                prev = label1;
            }        
        
        fromBeginning(*f1,twoBodyRitz,prev);
    }
        f1->name[twoBodyRitz].Partition =  allowQ(f1->rt, blockMatrixElementsblock) * maxArray;
        f1->name[twoBodyRitz].spinor = real;
        f1->name[twoBodyRitz].memory = bufferAllocation;
        
        fromBeginning(*f1,matrixHbuild,twoBodyRitz);
    f1->name[matrixHbuild].Partition = allowQ(f1->rt, blockMatrixElementsblock) *  (maxArray*maxArray);
   if ( f->i.cmpl == cmpl )
       f1->name[matrixHbuild].Partition *= 4;
        f1->name[matrixHbuild].spinor = real;
        f1->name[matrixHbuild].memory = bufferAllocation;

        fromBeginning(*f1,vectorHbuild,matrixHbuild);
        f1->name[vectorHbuild].Partition = 0;
        f1->name[vectorHbuild].memory = bufferAllocation;

        fromBeginning(*f1,matrixSbuild,vectorHbuild);
        f1->name[matrixSbuild].Partition =  allowQ(f1->rt, blockMatrixElementsblock) * (maxArray*maxArray);
    if ( f->i.cmpl == cmpl )
        f1->name[matrixSbuild].Partition *= 4;

        f1->name[matrixSbuild].spinor = real;
        f1->name[matrixSbuild].memory = bufferAllocation;

        fromBeginning(*f1,dsyBuffers,matrixSbuild);
#ifdef APPLE
        ///this count needs work
        f1->name[dsyBuffers].Partition = 1000;
#endif
        f1->name[dsyBuffers].spinor = parallel;
        f1->name[dsyBuffers].memory = bufferAllocation;

        fromBeginning(*f1,f1->end,dsyBuffers);
        fromBeginning(*f1,f1->end,f1->end);

    
        {
              if ( f1->boot == fullMatrices ){
            printf(" b\tp\ty\t1\t<\t>\n");
              }
            double maxMem = 0.,currMem;
              bodyType body;
            for ( space = 0 ; space <= SPACE ; space++){
                currMem = (f1->name[f1->end].space[space].Address)/(1000000000./(sizeof(floata)));
             
                if ( f1->boot == fullMatrices ){
                    if ( f1->canon[space].body != nada ){
                        for ( body = one ; body <= f1->canon[space].body ; body++){
                            printf(" %d\t%d\t%d\t%d",body,f1->canon[space].label, f1->canon[space].basis, f1->canon[space].count1Basis);
                        if ( f1->canon[space].basis == 1 )
                            printf("\t%1.1f\t%1.1f\t\n",f1->canon[space].particle[body].origin, f1->canon[space].particle[body].origin+(f1->canon[space].particle[body].lattice *(f1->canon[space].count1Basis-1)));
                        else  if ( f1->canon[space].basis == 4 )
                            printf("\t0\t%d\n",f1->canon[space].count1Basis-1);
                        else
                            printf("\n");
                        }
                    }
                }
                if ( currMem >= 0 )
                    maxMem += currMem;
                else {
                    exit(0);
                }
            }
              if ( f1->boot == fullMatrices )
            printf("\n\nmaxRAM %1.3f\n\n",maxMem);
            if ( maxMem > c1->i.RAMmax ){
                printf("oops too much RAM required\n %f > %d\n",maxMem, c1->i.RAMmax);
                fflush(stdout);

                exit(0);
            }
        }
    
        for ( space = 0; space <= SPACE ; space++){
            f1->canon[space].stream = malloc( (f1->name[f1->end].space[space].Address)*sizeof(floata));
        }

        f1->bootedMemory = 1;
        assignCores(*f1, 1);

        inta RdsSize;
        RdsSize = 0;
    
    
    f1->name[Ha].name = nullName;
    f1->name[Ha].linkNext = Iterator;
    f1->name[Iterator].name = nullName;

    
    if ( f1->boot==fullMatrices )
        defineTerms(c1, f1,Iterator,f->i.OpIndex);
#if VERBOSE
    linkDetails(*f1, Iterator);
#endif
    
    return 0;
}
