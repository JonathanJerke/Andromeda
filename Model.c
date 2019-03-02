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
        i->i.c.atoms[1].label.Z = 1;
        i->i.c.atoms[1].position[1] = 0;
        i->i.c.atoms[1].position[2] = 0;
        i->i.c.atoms[1].position[3] = 0;
        
    } else
        if ( number== 2 ){
            i->i.c.Na = 2;

            i->i.c.atoms[1].label.Z = 1;
            i->i.c.atoms[1].position[1] = 2*scale;
            i->i.c.atoms[1].position[2] = 0;
            i->i.c.atoms[1].position[3] = 0;
            i->i.c.atoms[2].label.Z = 1;
            i->i.c.atoms[2].position[1] = -2*scale;
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
        i->i.charge = i->i.c.Na - i->rt.body;
    }


struct calculation initCal (void ) {
    INT_TYPE space;
    struct calculation i;
    i.i.irrep = 0;
    i.i.hartreeFockFlag = 0;
    i.i.springConstant = 0.;
    i.i.springFlag = 0;
    i.i.potentialFlag = 0;
    i.i.RAMmax = 1;
    i.rt.runFlag = 0;
    i.i.bRank = 5;
    i.i.vectorMomentum = 0.;
    i.i.decomposeRankMatrix = 3;
    i.i.iRank = 1;
    i.i.qFloor = 0;
    i.i.Angstroms = 0;
    i.i.cycles = 1;
    i.i.Iterations = 1;
    i.rt.printFlag = 0;
    i.i.heliumFlag = 1;
    i.rt.body = 1;
    i.i.iCharge = i.rt.body;
    
    i.rt.TARGET = 1e-3;
    i.rt.ALPHA = 1e-6;
    i.rt.CANON = 1e-7;
    i.rt.vCANON = 1e-6;
    i.rt.TOL = 1e-6;
//    i.rt.CONVERGENCE = 1e-6;
//    i.rt.vCONVERGENCE = 1e-6;
    i.rt.maxEntropy = 1e-6;
    
    
//    i.p.iCondition = 5;
//    i.p.iTolerance = 1;
//    i.p.iThreshold = 2;
//    i.p.iTarget = 5;
//    i.p.iConvergence = 1;
//    i.p.iEntropy = 1;
//    i.p.vectorThreshold = 1;
//    i.p.vectorConvergence = 1;
    
    i.i.sectors = 1;
    i.i.d = 1;
    i.i.epi = 4;
    i.i.M1 = 4;
    resetExternal(&i, 1, 2.);
    
    for( space = 0 ; space <= SPACE ; space++){
        i.i.c.sinc.rose[space].stream = NULL;
    }
    i.i.c.sinc.tulip = NULL;
    
    
    i.i.c.twoBody.func.fn = nullFunction;
    i.i.c.oneBody.func.fn = Pseudo;

    i.i.canonRank = 15;
    i.i.c.twoBody.num = i.i.canonRank;
    i.i.c.twoBody.func.interval  = 0;
    i.i.c.twoBody.func.param[0]  = 1;
    i.i.c.twoBody.func.param[1]  = 1;
    i.i.c.twoBody.func.param[2]  = 1;

    i.i.c.oneBody.num = i.i.canonRank;
    i.i.c.oneBody.func.interval  = 0;
    i.i.c.oneBody.func.param[0]  = 1;
    i.i.c.oneBody.func.param[1]  = 1;
    i.i.c.oneBody.func.param[2]  = 1;

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
            fprintf(stdout, "SPACE %d\n", i);
            fflush(stdout);
            free(f1->sinc.rose[i].stream);
            f1->sinc.rose[i].stream = NULL;
        }
        {
            fprintf(stdout, "TULIP\n" );
            fflush(stdout);

            free(f1->sinc.tulip);
            f1->sinc.tulip = NULL;
        }
    }
    return 0;
}


INT_TYPE iModel( struct calculation * c1){
    //c1->rt.printFlag = 1;
    
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
    

    
    struct runTime * rt = &c1->rt;
    //rt->eigenTime = 0.;
    //rt->lanczosTime = 0.;
    //rt->buildTime = 0.;
//    for ( c = 0 ; c < MaxCore ; c++)
//        for ( i = 0; i < MaxCycle ; i++){
//            rt->timeBeylkin[c][i] = 0.;
//            rt->timeGEMV[c][i] = 0.;
//            rt->distBeylkin[c][i] = 0.;
//            rt->count[c][i] = 0;
//        }
    
    if ( c1->i.iCharge > 0 )
        f1->Ne = c1->i.iCharge;
    
    
    {
        INT_TYPE space;
        struct field * f1  =&(c1->i.c);
        {
            c1->mem.rt = rt;
            f1->mem1 = &c1->mem;
            
        }
        
        enum body bootBodies = c1->rt.body;
        enum bodyType bootType = c1->rt.bodyType;
        c1->i.c.body = bootBodies;
        INT_TYPE ra = tPerms(bootBodies),nG = tSize(bootBodies);
        c1->i.paths = tPaths(bootBodies, c1->i.irrep);
        printf("'states' %d\n helium %d \n floor %d\n", c1->i.nStates,c1->i.heliumFlag,c1->i.qFloor);
        INT_TYPE N12 = c1->i.epi;
        f1->sinc.N1 = 2*N12+1;
        INT_TYPE N1 = f1->sinc.N1;

        
        if (bootType == electron ){
            INT_TYPE p,s;
            for ( p = 1 ; p < 7 ; p++)
                for ( s = 0 ; s < SPACE ; s++)
                    f1->sinc.dims[p][s] = N1;
            
        }else if ( bootType == h2plus){
            INT_TYPE p,s;
            for ( p = 1 ; p < 7 ; p++)
                for ( s = 0 ; s < SPACE ; s++)
                    f1->sinc.dims[p][s] = N1;
        }
        
        
        f1->sinc.Basis[0] = N1;
        f1->sinc.Basis[1] = N1;
        f1->sinc.Basis[2] = N1;
        
        f1->sinc.Basis2[0] = f1->sinc.Basis[0]*f1->sinc.Basis[0];
        f1->sinc.Basis2[1] = f1->sinc.Basis[1]*f1->sinc.Basis[1];
        f1->sinc.Basis2[2] = f1->sinc.Basis[2]*f1->sinc.Basis[2];
        
        f1->sinc.Basis3[0] = f1->sinc.Basis[0]*f1->sinc.Basis[0]*f1->sinc.Basis[0];
        f1->sinc.Basis3[1] = f1->sinc.Basis[1]*f1->sinc.Basis[1]*f1->sinc.Basis[1];
        f1->sinc.Basis3[2] = f1->sinc.Basis[2]*f1->sinc.Basis[2]*f1->sinc.Basis[2];
        
        f1->sinc.Basis4[0] = f1->sinc.Basis[0]*f1->sinc.Basis[0]*f1->sinc.Basis[0]*f1->sinc.Basis[0];
        f1->sinc.Basis4[1] = f1->sinc.Basis[1]*f1->sinc.Basis[1]*f1->sinc.Basis[1]*f1->sinc.Basis[1];
        f1->sinc.Basis4[2] = f1->sinc.Basis[2]*f1->sinc.Basis[2]*f1->sinc.Basis[2]*f1->sinc.Basis[2];
        
        enum shape bootShape;
        INT_TYPE maxVector = imax(c1->i.decomposeRankMatrix, imax(c1->i.bRank,imax(1+c1->i.iRank,c1->i.Iterations+c1->i.iRank)));
        
        //rds defined in input.c
        
        bootShape = Cube;
        if ( 0)
            bootShape = Rds;
        f1->sinc.d = c1->i.d;
        
        
        INT_TYPE FloorV = imax(0, c1->i.qFloor), CeilV = imax(0,0);
        INT_TYPE maxArray,EV,maxEV,NV = 0,FV = FloorV+CeilV ;
        INT_TYPE maxDensity = c1->i.densityFlag;
        
        EV = FV;
        maxEV =EV*(imax(c1->i.Iterations,1));
        maxArray = imax(c1->i.nStates,maxEV);//slip Nb into spectra...
        
        f1->sinc.maxEV = maxArray;
        
        enum division density  = eigenVectors +  c1->i.nStates+maxEV;
        f1->sinc.density = density;
        enum division end  = density+maxDensity;
        f1->sinc.end = end;
        f1->sinc.tulip = malloc ( (end+1) * sizeof(struct name_label));


        INT_TYPE matrixNumber =   c1->i.decomposeRankMatrix;
        INT_TYPE outVector  = imax(f1->oneBody.num , f1->twoBody.num )*maxVector;
        
        {
         //   initPointerTensors(f1);
            enum division label1;
            for ( label1 = 0 ;label1 <= end; label1++){
                f1->sinc.tulip[label1].name = label1;
                f1->sinc.tulip[label1].header = Cube;
                f1->sinc.tulip[label1].spinor = cmpl;
                f1->sinc.tulip[label1].purpose = Object;
                f1->sinc.tulip[label1].memory = threeObject;
                f1->sinc.tulip[label1].species = scalar;
                f1->sinc.tulip[label1].linkNext = nullName;
                f1->sinc.tulip[label1].blockType = tv1;
                f1->sinc.tulip[label1].NBody = bootBodies;
                f1->sinc.tulip[label1].TBody = electron;
                f1->sinc.tulip[label1].Partition = 0;
                f1->sinc.tulip[label1].parallel = 0;
                tClear(f1,label1);
                
            }
        }
        
        if(1){
            printf("TARGET \t\t%f\n", log(c1->rt.TARGET)/log(10));//quality of decomposition
            printf("TOL \t\t%f\n", log(c1->rt.TOL)/log(10));//max condition of foundation vectors
            printf("CANON \t\t%f\n", log(c1->rt.CANON)/log(10));//matrix training standard
            printf("vCANON \t\t%f\n", log(c1->rt.vCANON)/log(10));//vector training standard
        //    printf("CONVERGENCE \t%f\n", log(c1->rt.CONVERGENCE )/log(10));
        //    printf("vCONVERGENCE\t %f\n", log(c1->rt.vCONVERGENCE )/log(10));
            printf("ALPHA\t\t %f\n", log(c1->rt.ALPHA )/log(10));//Beylkin parameter

        }
        
//  f1->sinc.tulip[kinetic].spinor = cmpl;
        f1->sinc.tulip[kinetic].Address =0;
        f1->sinc.tulip[kinetic].Partition = SPACE;//
        f1->sinc.tulip[kinetic].NBody = one;
        f1->sinc.tulip[kinetic].species = matrix;
        f1->sinc.tulip[kinetic].header = Cube;
//        f1->sinc.tulip[kinetic].symmetryType = nullSymmetry;
        f1->sinc.tulip[kinetic].name = kinetic;
        
 //       f1->sinc.tulip[kineticMass].spinor = cmpl;
        f1->sinc.tulip[kineticMass].Address =fromBegining(f1,kinetic);
        f1->sinc.tulip[kineticMass].Partition = SPACE;//
        f1->sinc.tulip[kineticMass].NBody = one;
        f1->sinc.tulip[kineticMass].species = matrix;
        f1->sinc.tulip[kineticMass].header = Cube;
//        f1->sinc.tulip[kineticMass].symmetryType = nullSymmetry;
        f1->sinc.tulip[kineticMass].name = kineticMass;

        f1->sinc.tulip[harmonium].spinor = none;
        f1->sinc.tulip[harmonium].Address =fromBegining(f1,kineticMass);
        f1->sinc.tulip[harmonium].Partition = c1->i.springFlag*SPACE;//
        f1->sinc.tulip[harmonium].NBody = one;
        f1->sinc.tulip[harmonium].species = matrix;
        f1->sinc.tulip[harmonium].header = Cube;
//        f1->sinc.tulip[harmonium].symmetryType = nullSymmetry;
        f1->sinc.tulip[harmonium].name = harmonium;
        
        f1->sinc.tulip[X].spinor = none;
        f1->sinc.tulip[X].Address = fromBegining(f1,harmonium);
        f1->sinc.tulip[X].Partition =  1 ;
        f1->sinc.tulip[X].NBody = one;
        f1->sinc.tulip[X].species = matrix;
        f1->sinc.tulip[X].header = Cube;
        f1->sinc.tulip[X].purpose = Object;
//        f1->sinc.tulip[X].symmetryType = nullSymmetry;
        f1->sinc.tulip[X].name = X;
        
        
        f1->sinc.tulip[linear].spinor = none;
        f1->sinc.tulip[linear].Address = fromBegining(f1,X);
        f1->sinc.tulip[linear].Partition = (f1->Na*f1->oneBody.num + f1->twoBody.num + 1 );//
        f1->sinc.tulip[linear].NBody = one;
        f1->sinc.tulip[linear].species = matrix;
        f1->sinc.tulip[linear].header = Cube;
//        f1->sinc.tulip[linear].symmetryType = nullSymmetry;
        f1->sinc.tulip[linear].name = linear;
        
        f1->sinc.tulip[overlap].spinor = none;
        f1->sinc.tulip[overlap].Address = fromBegining(f1,linear);
        f1->sinc.tulip[overlap].Partition = 0;//
        f1->sinc.tulip[overlap].species = matrix;
        f1->sinc.tulip[overlap].header = bootShape;
//        f1->sinc.tulip[overlap].symmetryType = nullSymmetry;
        f1->sinc.tulip[overlap].name = overlap;
        
        
        {
            f1->sinc.tulip[build].spinor = none;
            f1->sinc.tulip[build].Address = fromBegining(f1,overlap);
            
            enum body runBodies = bootBodies;
//            if ( bootBodies == one ){
//                runBodies = one;
//            }else if ( bootBodies == three){
//                runBodies = three;
//            }
            
            f1->sinc.tulip[build].Partition = c1->i.sectors*runBodies*((part(f1, kinetic)+part(f1, harmonium)) + matrixNumber * f1->Na);//easily reduce in cheaper ways!
            if ( runBodies == two )
                f1->sinc.tulip[build].Partition += c1->i.sectors+ ( bootType > electron );
            else if ( runBodies == three )
                f1->sinc.tulip[build].Partition += 3*c1->i.sectors;
            else if ( runBodies == four )
                f1->sinc.tulip[build].Partition += 6*c1->i.sectors;
            
            f1->sinc.tulip[build].Partition = imax(nG*c1->i.sectors,part(f1, build));
            
            // B2 :  2*(1+S)*3 + 1*sector +  Na*matrixNumber*2
            // B3 :  3*(1+S)*3 + 3*sector  + Na*matrixNumber*3
            f1->sinc.tulip[build].NBody = runBodies;
            f1->sinc.tulip[build].species = matrix;
            f1->sinc.tulip[build].header = Cube;
            f1->sinc.tulip[build].purpose = Object;
//            f1->sinc.tulip[build].symmetryType = nullSymmetry;
            f1->sinc.tulip[build].name = build;
            
            f1->sinc.tulip[eigen].Address = fromBegining(f1,build);
            f1->sinc.tulip[eigen].Partition = c1->i.sectors;
            f1->sinc.tulip[eigen].NBody = runBodies;
            f1->sinc.tulip[eigen].species = matrix;
            f1->sinc.tulip[eigen].header = Cube;
            f1->sinc.tulip[eigen].purpose = tObject;
//            f1->sinc.tulip[eigen].symmetryType = nullSymmetry;
            f1->sinc.tulip[eigen].name = eigen;
        
            
            INT_TYPE si = 0,g;
            for ( g = 0; g < nSAG*nSAG*nSAG ; g++)
                si += c1->i.cSA[g];
            
            if ( c1->i.densityFlag ){
                INT_TYPE di;
                printf("std DENSITY vectors with %d ranks\n\n", c1->i.bRank);
                f1->sinc.tulip[density].Address = fromBegining(f1,eigen);
                f1->sinc.tulip[density].Partition = c1->i.dRank;
                printf("density rakn %d\n", c1->i.dRank);
                f1->sinc.tulip[density].species = vector;
                f1->sinc.tulip[density].NBody = c1->i.bodyDensity;
                f1->sinc.tulip[density].header = Cube;
                f1->sinc.tulip[density].purpose = Object;
                f1->sinc.tulip[density].name = density;
                
                for ( di = 1 ; density+di < end; di++){
                    f1->sinc.tulip[density+di-1].linkNext = density+di;// linked up!
                    f1->sinc.tulip[density+di].Address = fromBegining(f1,density+di-1);
                    f1->sinc.tulip[density+di].Partition = c1->i.dRank ;
                    f1->sinc.tulip[density+di].header = Cube;
                    f1->sinc.tulip[density+di].species = vector;
                    f1->sinc.tulip[density+di].NBody = c1->i.bodyDensity;
                    f1->sinc.tulip[density+di].purpose = Object;
                    f1->sinc.tulip[density+di].name = density+di;
                }
                f1->sinc.tulip[eigenVectors].Address = fromBegining(f1,density+di-1);
                
            }else{
                f1->sinc.tulip[eigenVectors].Address = fromBegining(f1,eigen);
            }
            
            
            if(1 /* !si*/){
                INT_TYPE di,d0=1;
                printf("std USERs\n\n");
//                f1->sinc.tulip[eigenVectors].Address = fromBegining(f1,density);
                f1->sinc.tulip[eigenVectors].Partition = c1->i.bRank;
                f1->sinc.tulip[eigenVectors].species = vector;
                f1->sinc.tulip[eigenVectors].header = Cube;//READ AND WRITE | => CUBE
                f1->sinc.tulip[eigenVectors].purpose = tObject;
                f1->sinc.tulip[eigenVectors].name = eigenVectors;
                
                for ( di = 1 ; eigenVectors+di < density; di++){
                    f1->sinc.tulip[eigenVectors+di].Address = fromBegining(f1,eigenVectors+di-1);
                    if ( di < c1->i.nStates ){
                        f1->sinc.tulip[eigenVectors+di].Partition = c1->i.bRank;
                        d0++;
                    }
                    else if ( di < c1->i.nStates+maxEV){
                        {

                            INT_TYPE nextRank = ((di-d0)/EV)+c1->i.iRank;
                            if ( (c1->i.iRank > 12) ){
                                if ( !(((di-d0)/EV) % 2) ){
                                   }
                                else
                                    nextRank--;
                            }
                            
                            f1->sinc.tulip[eigenVectors+di].Partition = nextRank;
                          //  printf("new %d %d\n",eigenVectors+di, part(f1, eigenVectors+di) );
                          //  f1->sinc.tulip[eigenVectors+di].spinor = none;

                            NV += spins(f1,eigenVectors+di )*f1->sinc.tulip[eigenVectors+di].Partition;
                        }
                        
                    }
                        f1->sinc.tulip[eigenVectors+di].header = Cube;
                        f1->sinc.tulip[eigenVectors+di].species = vector;
                        f1->sinc.tulip[eigenVectors+di].purpose = tObject;
                        f1->sinc.tulip[eigenVectors+di].name = eigenVectors+di;
                        
                    }
                    f1->sinc.tulip[diagonalVectorA].Address = fromBegining(f1,eigenVectors+di-1);
                    
            }else
            {
                INT_TYPE di,d0=1,ii,jj,kk,h,i,last,fpath,v;
                printf("path USERs\n\n");
                fflush(stdout);
                exit(0);
                f1->sinc.tulip[eigenVectors].Address = fromBegining(f1,density);
                f1->sinc.tulip[eigenVectors].Partition = c1->i.bRank;
                f1->sinc.tulip[eigenVectors].species = vector;
                f1->sinc.tulip[eigenVectors].header = Cube;//READ AND WRITE | => CUBE
                f1->sinc.tulip[eigenVectors].purpose = tObject;
                f1->sinc.tulip[eigenVectors].name = eigenVectors;
                
                for ( di = 1 ; di < c1->i.nStates; di++){
                    f1->sinc.tulip[eigenVectors+di].Address = fromBegining(f1,eigenVectors+di-1);
                    f1->sinc.tulip[eigenVectors+di].Partition = c1->i.bRank;
                    f1->sinc.tulip[eigenVectors+di].species = vector;
                    f1->sinc.tulip[eigenVectors+di].header = Cube;//READ AND WRITE | => CUBE
                    f1->sinc.tulip[eigenVectors+di].purpose = tObject;
                    f1->sinc.tulip[eigenVectors+di].name = eigenVectors+di;
                    d0++;
                }
                
                last = d0 -1;
                fpath = 0;
                v=0;
                for ( v = 0 ; v < nG*nG*nG ; v++ )
                {
                    ii = v % nG;
                    jj = (v/nG)%nG;
                    kk = (v/(nG*nG))%nG;
                    if ( tIR(bootBodies, ii, jj, kk, c1->i.irrep)){
                        for ( h = 0; h < c1->i.cSA[v] ;h++)
                            for ( i = 0; i < c1->i.Iterations ; i++)//ORDER IN MEMORY
                                
                            {
                                enum division nm = eigenVectors+d0+ h+fpath + i*si;//ORDER IN USR
                                
  //                             printf("new %d rank %d last %d path %d \n",nm, i+c1->i.iRank,eigenVectors+last,v);//
//                               fflush(stdout);
                                f1->sinc.tulip[nm].Address = fromBegining(f1,eigenVectors+last);
                              //   f1->sinc.tulip[nm].spinor = none;

                                last = nm-eigenVectors;
                                f1->sinc.tulip[nm].Partition = i+c1->i.iRank;
                                NV += spins(f1,nm)*(i+c1->i.iRank);
                               // printf("spines %d\n",spins(f1,nm));
                                f1->sinc.tulip[nm].header = Cube;
                                f1->sinc.tulip[nm].species = vector;
                                f1->sinc.tulip[nm].purpose = tObject;
                                f1->sinc.tulip[nm].name = nm;
                                f1->sinc.tulip[nm].path = v;
                            }
                        fpath+=c1->i.cSA[v];
                    }
                    f1->sinc.tulip[diagonalVectorA].Address = fromBegining(f1,eigenVectors+last );
                    
                }
            }
                
            
            
            
            f1->sinc.tulip[diagonalVectorA].Partition = !(!(c1->i.sectors));
            f1->sinc.tulip[diagonalVectorA].species = vector;
            f1->sinc.tulip[diagonalVectorA].parallel = 2;
            f1->sinc.tulip[diagonalVectorA].header = Cube;
            f1->sinc.tulip[diagonalVectorA].NBody = runBodies;
            f1->sinc.tulip[diagonalVectorA].purpose = tObject;
//            f1->sinc.tulip[diagonalVectorA].symmetryType = nullSymmetry;
            f1->sinc.tulip[diagonalVectorA].name = diagonalVectorA;
            
            f1->sinc.tulip[diagonalVectorB].Address = fromBegining(f1,diagonalVectorA);
            f1->sinc.tulip[diagonalVectorB].Partition =  !(!(c1->i.sectors));
            f1->sinc.tulip[diagonalVectorB].species = vector;
            f1->sinc.tulip[diagonalVectorB].parallel = 2;
            f1->sinc.tulip[diagonalVectorB].header = Cube;
            f1->sinc.tulip[diagonalVectorB].NBody = runBodies;
            f1->sinc.tulip[diagonalVectorB].purpose = tObject;
//            f1->sinc.tulip[diagonalVectorB].symmetryType = nullSymmetry;
            f1->sinc.tulip[diagonalVectorB].name = diagonalVectorB;

        }
        
        f1->sinc.tulip[diagonalVector].Address = fromBegining(f1,diagonalVectorB);
        f1->sinc.tulip[diagonalVector].Partition =  maxVector;
        f1->sinc.tulip[diagonalVector].species = vector;
        f1->sinc.tulip[diagonalVector].header = Cube;
        f1->sinc.tulip[diagonalVector].parallel = 2;
        f1->sinc.tulip[diagonalVector].purpose = tObject;
//        f1->sinc.tulip[diagonalVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[diagonalVector].name = diagonalVector;

        
        f1->sinc.tulip[entropyVector].Address = fromBegining(f1,diagonalVector);
        f1->sinc.tulip[entropyVector].spinor = none;
        f1->sinc.tulip[entropyVector].Partition = 0;//
        if ( 0 )
            f1->sinc.tulip[entropyVector].NBody = 0;
        f1->sinc.tulip[entropyVector].species = vector;
        f1->sinc.tulip[entropyVector].header = Cube;
        f1->sinc.tulip[entropyVector].purpose = Object;
//        f1->sinc.tulip[entropyVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[entropyVector].name = entropyVector;
        
        f1->sinc.tulip[entropyUnit].spinor = none;
        f1->sinc.tulip[entropyUnit].Address = fromBegining(f1,entropyVector);
        f1->sinc.tulip[entropyUnit].Partition = !(!(part(f1, entropyVector)))*c1->i.sectors;//
        f1->sinc.tulip[entropyUnit].NBody = bodies(f1, entropyVector);
        f1->sinc.tulip[entropyUnit].species = vector;
        f1->sinc.tulip[entropyUnit].header = Cube;
        f1->sinc.tulip[entropyUnit].purpose = tObject;
//        f1->sinc.tulip[entropyUnit].symmetryType = nullSymmetry;
        f1->sinc.tulip[entropyUnit].name = entropyUnit;
        
        f1->sinc.tulip[complement].Address = fromBegining(f1,entropyUnit);
        f1->sinc.tulip[complement].Partition = !(!c1->i.densityFlag ) * c1->i.dRank * c1->i.canonRank;//
        f1->sinc.tulip[complement].parallel = 1;
        f1->sinc.tulip[complement].species = vector;
        if ( c1->i.densityFlag && (c1->i.bodyDensity - bootBodies) > 0)
            f1->sinc.tulip[complement].NBody = (c1->i.bodyDensity - bootBodies);
        f1->sinc.tulip[complement].header = Cube;
        f1->sinc.tulip[complement].purpose = Object;
        f1->sinc.tulip[complement].name = complement;
        
        f1->sinc.tulip[complementTwo].Address = fromBegining(f1,entropyUnit);
        f1->sinc.tulip[complementTwo].Partition = !(!c1->i.densityFlag ) * c1->i.dRank * c1->i.canonRank;//
        f1->sinc.tulip[complementTwo].parallel = 1;
        f1->sinc.tulip[complementTwo].species = vector;
        if ( c1->i.densityFlag && (c1->i.bodyDensity - bootBodies) > 0)
            f1->sinc.tulip[complementTwo].NBody = (c1->i.bodyDensity - bootBodies);
        f1->sinc.tulip[complementTwo].header = Cube;
        f1->sinc.tulip[complementTwo].purpose = Object;
        f1->sinc.tulip[complementTwo].name = complementTwo;

        
        f1->sinc.tulip[MomentumDot].spinor = none;
        f1->sinc.tulip[MomentumDot].Address = fromBegining(f1,complementTwo);
        f1->sinc.tulip[MomentumDot].Partition =0 ;//
        f1->sinc.tulip[MomentumDot].NBody = two;
        f1->sinc.tulip[MomentumDot].species = matrix;
        f1->sinc.tulip[MomentumDot].header = Cube;
        f1->sinc.tulip[MomentumDot].purpose = Object;
//        f1->sinc.tulip[MomentumDot].symmetryType = nullSymmetry;
        f1->sinc.tulip[MomentumDot].name = MomentumDot;
        
        f1->sinc.tulip[edgeMatrix].spinor = none;
        f1->sinc.tulip[edgeMatrix].Address = fromBegining(f1,MomentumDot);
        f1->sinc.tulip[edgeMatrix].Partition = 6 ;//
        f1->sinc.tulip[edgeMatrix].NBody = one;
        f1->sinc.tulip[edgeMatrix].species = matrix;
        f1->sinc.tulip[edgeMatrix].header = Cube;
        f1->sinc.tulip[edgeMatrix].purpose = Object;
//        f1->sinc.tulip[edgeMatrix].symmetryType = nullSymmetry;
        f1->sinc.tulip[edgeMatrix].name = edgeMatrix;
        
        f1->sinc.tulip[nullVector].Address = fromBegining(f1,edgeMatrix);
        f1->sinc.tulip[nullVector].Partition = 0;
        f1->sinc.tulip[nullVector].species = scalar;
        f1->sinc.tulip[nullVector].header = Cube;
        f1->sinc.tulip[nullVector].purpose = ptObject;
        f1->sinc.tulip[nullVector].name = nullVector;
        
        f1->sinc.tulip[nullMatrix].Address = fromBegining(f1,nullVector);
        f1->sinc.tulip[nullMatrix].Partition = 0;
        f1->sinc.tulip[nullMatrix].species = scalar;
        f1->sinc.tulip[nullMatrix].header = Cube;
        f1->sinc.tulip[nullMatrix].purpose = ptObject;
        f1->sinc.tulip[nullMatrix].name = nullMatrix;
        
        f1->sinc.tulip[vectorCubeBuffers].Address = fromBegining(f1,nullMatrix);
        f1->sinc.tulip[vectorCubeBuffers].Partition = 0;
        f1->sinc.tulip[vectorCubeBuffers].parallel = 2;
        f1->sinc.tulip[vectorCubeBuffers].species = vector;
        f1->sinc.tulip[vectorCubeBuffers].header = Cube;//WARNING
        f1->sinc.tulip[vectorCubeBuffers].purpose = tObject;
        f1->sinc.tulip[vectorCubeBuffers].name = vectorCubeBuffers;
        
        f1->sinc.tulip[secondVector].Address = fromBegining(f1,vectorCubeBuffers);
        f1->sinc.tulip[secondVector].Partition = 0*(!(!c1->rt.printFlag)) * maxVector + ( ( !(!c1->i.sectors)))* 24;//(!(!Nb))*outVector ;
        f1->sinc.tulip[secondVector].species = vector;
        f1->sinc.tulip[secondVector].header = Cube;
        f1->sinc.tulip[secondVector].parallel = 1;
        f1->sinc.tulip[secondVector].purpose = tObject;
//        f1->sinc.tulip[secondVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[secondVector].name = secondVector;
                
        f1->sinc.tulip[edges].Partition = 0;
        f1->sinc.tulip[edges].Address = fromBegining(f1,secondVector);
        f1->sinc.tulip[edges].species = vector;
        f1->sinc.tulip[edges].header = Cube;
        f1->sinc.tulip[edges].purpose = Object;
        f1->sinc.tulip[edges].name = edges;
        
        f1->sinc.tulip[basis].spinor = none;
        f1->sinc.tulip[basis].Address = fromBegining(f1,edges);
        f1->sinc.tulip[basis].Partition = 0;
        f1->sinc.tulip[basis].species = matrix;
        f1->sinc.tulip[basis].header = RdsBasis;
        f1->sinc.tulip[basis].purpose = Object;
//        f1->sinc.tulip[basis].symmetryType = nullSymmetry;
        f1->sinc.tulip[basis].name = basis;
        
        f1->sinc.tulip[bandBasis].spinor = none;
        f1->sinc.tulip[bandBasis].Address = fromBegining(f1,basis);
        f1->sinc.tulip[bandBasis].Partition = 0;
        f1->sinc.tulip[bandBasis].species = matrix;
        f1->sinc.tulip[bandBasis].header = BandBasis;
        f1->sinc.tulip[bandBasis].purpose = Object;
//        f1->sinc.tulip[bandBasis].symmetryType = nullSymmetry;
        f1->sinc.tulip[bandBasis].name = bandBasis;
        
        
        INT_TYPE largestRank = 0;
        INT_TYPE largestRankMatrix = 0;
        INT_TYPE maxTrainingRank = 0;
        INT_TYPE maxVectorTrain = 0;
        INT_TYPE maxMatrixTrain = 0;
        INT_TYPE maxQuarticTrain = 0;
        
        
        
        enum division label;
        for ( label = 0; label < end ; label++){
            if ( part(f1, label ) > largestRankMatrix &&species(f1,label ) == matrix && purpose(f1, label) == Object)
                largestRankMatrix = part(f1, label);
            if ( part(f1, label ) > maxTrainingRank && purpose(f1, label) == tObject)
                maxTrainingRank = part(f1, label);
            if (part(f1, label ) > maxVectorTrain && species(f1,label ) == vector&& purpose(f1, label) == tObject)
                maxVectorTrain = part(f1, label);
            if ( part(f1, label ) > maxMatrixTrain &&species(f1,label ) == matrix&& purpose(f1, label) == tObject)
                maxMatrixTrain = 0*part(f1, label);
            if ( part(f1, label ) > maxQuarticTrain &&species(f1,label ) == quartic&& purpose(f1, label) == tObject)
                maxQuarticTrain = 0*part(f1, label);
            
        }
        
        
        f1->sinc.tulip[productVector].Address = fromBegining(f1,bandBasis);
        f1->sinc.tulip[productVector].Partition = (c1->i.canonRank)*maxVector;
        f1->sinc.tulip[productVector].species = vector;
        f1->sinc.tulip[productVector].header = Cube;
        f1->sinc.tulip[productVector].parallel = 2;
        f1->sinc.tulip[productVector].purpose = tObject;
//        f1->sinc.tulip[productVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[productVector].name = productVector;
        
        f1->sinc.tulip[permutationVector].Address = fromBegining(f1,productVector);
        f1->sinc.tulip[permutationVector].Partition = ra*maxVector;
        f1->sinc.tulip[permutationVector].species = vector;
        f1->sinc.tulip[permutationVector].header = Cube;
        f1->sinc.tulip[permutationVector].parallel = 2;
        f1->sinc.tulip[permutationVector].purpose = Object;
//        f1->sinc.tulip[permutationVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[permutationVector].name = permutationVector;
        
        f1->sinc.tulip[copyVector].Address = fromBegining(f1,permutationVector);
        f1->sinc.tulip[copyVector].Partition = maxVector;//c1->i.bRank*c1->i.bRank;
        f1->sinc.tulip[copyVector].species = vector;
        f1->sinc.tulip[copyVector].header = Cube;
        if ( !(!c1->rt.printFlag) && bootBodies > two){
            f1->sinc.tulip[copyVector].NBody = two;
        }
        f1->sinc.tulip[copyVector].parallel = 2;
        f1->sinc.tulip[copyVector].purpose = tObject;
//        f1->sinc.tulip[copyVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[copyVector].name = copyVector;
        
        
        f1->sinc.tulip[copyTwoVector].Address = fromBegining(f1,copyVector);
        f1->sinc.tulip[copyTwoVector].Partition =   maxVector;//c1->i.bRank*c1->i.bRank;
        f1->sinc.tulip[copyTwoVector].species = vector;
        f1->sinc.tulip[copyTwoVector].header = Cube;
        if ( !(!c1->rt.printFlag)  && bootBodies > two){
            f1->sinc.tulip[copyTwoVector].NBody = two;
        }
        f1->sinc.tulip[copyTwoVector].parallel = 1;
        f1->sinc.tulip[copyTwoVector].purpose = tObject;
//        f1->sinc.tulip[copyTwoVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[copyTwoVector].name = copyTwoVector;
        
        f1->sinc.tulip[copyThreeVector].Address = fromBegining(f1,copyTwoVector);
        f1->sinc.tulip[copyThreeVector].Partition = maxVector;//c1->i.bRank*c1->i.bRank;
        f1->sinc.tulip[copyThreeVector].species = vector;
        f1->sinc.tulip[copyThreeVector].header = Cube;
        f1->sinc.tulip[copyThreeVector].parallel = 1;
        f1->sinc.tulip[copyThreeVector].purpose = tObject;
        f1->sinc.tulip[copyThreeVector].name = copyThreeVector;
        
        f1->sinc.tulip[copyFourVector].Address = fromBegining(f1,copyThreeVector);
        f1->sinc.tulip[copyFourVector].Partition = maxVector;//c1->i.bRank*c1->i.bRank;
        f1->sinc.tulip[copyFourVector].species = vector;
        f1->sinc.tulip[copyFourVector].header = Cube;
        f1->sinc.tulip[copyFourVector].parallel = 1;
        f1->sinc.tulip[copyFourVector].purpose = tObject;
        f1->sinc.tulip[copyFourVector].name = copyFourVector;
        
        f1->sinc.tulip[oneVector].Address = fromBegining(f1,copyFourVector);
        f1->sinc.tulip[oneVector].Partition = !(!c1->rt.printFlag) * 24*imax(f1->oneBody.num, f1->twoBody.num);;
        f1->sinc.tulip[oneVector].NBody = one;
        f1->sinc.tulip[oneVector].species = vector;
        f1->sinc.tulip[oneVector].header = Cube;
        f1->sinc.tulip[oneVector].parallel = 2;
        f1->sinc.tulip[oneVector].purpose = Object;
//        f1->sinc.tulip[oneVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[oneVector].name = oneVector;
        
        f1->sinc.tulip[twoVector].Address = fromBegining(f1,oneVector);
        f1->sinc.tulip[twoVector].Partition = !(!c1->rt.printFlag) * c1->i.bRank;;
        f1->sinc.tulip[twoVector].NBody = two;
        f1->sinc.tulip[twoVector].species = vector;
        f1->sinc.tulip[twoVector].header = Cube;
        f1->sinc.tulip[twoVector].parallel = 2;
        f1->sinc.tulip[twoVector].purpose = Object;
//        f1->sinc.tulip[twoVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[twoVector].name = twoVector;

        
        f1->sinc.tulip[squareVector].Address = fromBegining(f1,twoVector);
        f1->sinc.tulip[squareVector].Partition = 0*maxVector*imax(ra,c1->i.canonRank);
        f1->sinc.tulip[squareVector].species = vector;
        f1->sinc.tulip[squareVector].header = Cube;
        f1->sinc.tulip[squareVector].parallel = 2;
        f1->sinc.tulip[squareVector].purpose = tObject;
        f1->sinc.tulip[squareVector].name = squareVector;
        
        f1->sinc.tulip[totalVector].Address = fromBegining(f1,squareVector);
        f1->sinc.tulip[totalVector].Partition = (c1->i.canonRank)*part(f1,copyVector);
        f1->sinc.tulip[totalVector].species = vector;
        f1->sinc.tulip[totalVector].header = Cube;
        f1->sinc.tulip[totalVector].parallel = 1;//easily reduce this size...
        f1->sinc.tulip[totalVector].purpose = Object;
        f1->sinc.tulip[totalVector].name = totalVector;
        
        f1->sinc.tulip[diagonal].Address = fromBegining(f1,totalVector);
        f1->sinc.tulip[diagonal].Partition = 1;
        f1->sinc.tulip[diagonal].species = matrix;
        f1->sinc.tulip[diagonal].NBody = one;
        f1->sinc.tulip[diagonal].header = Cube;
        f1->sinc.tulip[diagonal].purpose = tObject;
//        f1->sinc.tulip[diagonal].symmetryType = nullSymmetry;
        f1->sinc.tulip[diagonal].name = diagonal;
        
        f1->sinc.tulip[diagonalCube].spinor = none;
        f1->sinc.tulip[diagonalCube].Address = fromBegining(f1,diagonal);
        f1->sinc.tulip[diagonalCube].Partition = 1;
        f1->sinc.tulip[diagonalCube].NBody = one;
        f1->sinc.tulip[diagonalCube].species = matrix;
        f1->sinc.tulip[diagonalCube].header = Cube;
        f1->sinc.tulip[diagonalCube].purpose = tObject;
//        f1->sinc.tulip[diagonalCube].symmetryType = nullSymmetry;
        f1->sinc.tulip[diagonalCube].name = diagonalCube;
        
        f1->sinc.tulip[diagonal1VectorA].spinor = none;
        f1->sinc.tulip[diagonal1VectorA].Address = fromBegining(f1,diagonalCube);
        f1->sinc.tulip[diagonal1VectorA].Partition =1+!(!c1->rt.printFlag);
        f1->sinc.tulip[diagonal1VectorA].parallel = 2;
        f1->sinc.tulip[diagonal1VectorA].species = vector;
        f1->sinc.tulip[diagonal1VectorA].header = Cube;
        f1->sinc.tulip[diagonal1VectorA].NBody = one;
        f1->sinc.tulip[diagonal1VectorA].purpose = tObject;
//        f1->sinc.tulip[diagonal1VectorA].symmetryType = nullSymmetry;
        f1->sinc.tulip[diagonal1VectorA].name = diagonal1VectorA;
                
        f1->sinc.tulip[diagonal2VectorA].spinor = none;
        f1->sinc.tulip[diagonal2VectorA].Address = fromBegining(f1,diagonal1VectorA);
        f1->sinc.tulip[diagonal2VectorA].Partition =1+!(!c1->rt.printFlag);
        f1->sinc.tulip[diagonal2VectorA].parallel = 2;
        f1->sinc.tulip[diagonal2VectorA].species = vector;
        f1->sinc.tulip[diagonal2VectorA].header = Cube;
        f1->sinc.tulip[diagonal2VectorA].NBody = two;
        f1->sinc.tulip[diagonal2VectorA].purpose = tObject;
//        f1->sinc.tulip[diagonal2VectorA].symmetryType = nullSymmetry;
        f1->sinc.tulip[diagonal2VectorA].name = diagonal2VectorA;
        
        f1->sinc.tulip[diagonal2VectorB].spinor = none;
        f1->sinc.tulip[diagonal2VectorB].Address = fromBegining(f1,diagonal2VectorA);
        f1->sinc.tulip[diagonal2VectorB].Partition =1+!(!c1->rt.printFlag);
        f1->sinc.tulip[diagonal2VectorB].parallel = 2;
        f1->sinc.tulip[diagonal2VectorB].species = vector;
        f1->sinc.tulip[diagonal2VectorB].header = Cube;
        f1->sinc.tulip[diagonal2VectorB].NBody = two;
        f1->sinc.tulip[diagonal2VectorB].purpose = tObject;
//        f1->sinc.tulip[diagonal2VectorB].symmetryType = nullSymmetry;
        f1->sinc.tulip[diagonal2VectorB].name = diagonal2VectorB;
        
        f1->sinc.tulip[diagonal1VectorB].spinor = none;
        f1->sinc.tulip[diagonal1VectorB].Address = fromBegining(f1,diagonal2VectorB);
        f1->sinc.tulip[diagonal1VectorB].Partition =1+!(!c1->rt.printFlag);
        f1->sinc.tulip[diagonal1VectorB].parallel = 2;
        f1->sinc.tulip[diagonal1VectorB].species = vector;
        f1->sinc.tulip[diagonal1VectorB].header = Cube;
        f1->sinc.tulip[diagonal1VectorB].NBody = one;
        f1->sinc.tulip[diagonal1VectorB].purpose = tObject;
//        f1->sinc.tulip[diagonal1VectorB].symmetryType = nullSymmetry;
        f1->sinc.tulip[diagonal1VectorB].name = diagonal1VectorB;
        //setBoxLimits(f1, diagonal1VectorB,N1);
        
        f1->sinc.tulip[diagonal1VectorC].spinor = none;
        f1->sinc.tulip[diagonal1VectorC].Address = fromBegining(f1,diagonal1VectorB);
        f1->sinc.tulip[diagonal1VectorC].Partition =1+!(!c1->rt.printFlag);
        f1->sinc.tulip[diagonal1VectorC].parallel = 2;
        f1->sinc.tulip[diagonal1VectorC].species = vector;
        f1->sinc.tulip[diagonal1VectorC].header = Cube;
        f1->sinc.tulip[diagonal1VectorC].NBody = one;
        f1->sinc.tulip[diagonal1VectorC].purpose = tObject;
//        f1->sinc.tulip[diagonal1VectorC].symmetryType = nullSymmetry;
        f1->sinc.tulip[diagonal1VectorC].name = diagonal1VectorC;
        //setBoxLimits(f1, diagonal1VectorC,N1);
        
        f1->sinc.tulip[diagonal1VectorD].spinor = none;
        f1->sinc.tulip[diagonal1VectorD].Address = fromBegining(f1,diagonal1VectorC);
        f1->sinc.tulip[diagonal1VectorD].Partition = 1+!(!c1->rt.printFlag);
        f1->sinc.tulip[diagonal1VectorD].parallel = 2;
        f1->sinc.tulip[diagonal1VectorD].species = vector;
        f1->sinc.tulip[diagonal1VectorD].header = Cube;
        f1->sinc.tulip[diagonal1VectorD].NBody = one;
        f1->sinc.tulip[diagonal1VectorD].purpose = tObject;
//        f1->sinc.tulip[diagonal1VectorD].symmetryType = nullSymmetry;
        f1->sinc.tulip[diagonal1VectorD].name = diagonal1VectorD;
        //setBoxLimits(f1, diagonal1VectorD,N1);
        
        f1->sinc.tulip[seconds].spinor = none;
        f1->sinc.tulip[seconds].Address = fromBegining(f1,diagonal1VectorD);
        f1->sinc.tulip[seconds].Partition = 0 ;
        f1->sinc.tulip[seconds].parallel = 0;
        f1->sinc.tulip[seconds].species = matrix;
        f1->sinc.tulip[seconds].header = Cube;
        f1->sinc.tulip[seconds].NBody = one;
        f1->sinc.tulip[seconds].purpose = tObject;
//        f1->sinc.tulip[seconds].symmetryType = nullSymmetry;
        f1->sinc.tulip[seconds].name = seconds;
        
        f1->sinc.tulip[copy].Address = fromBegining(f1,seconds);
        f1->sinc.tulip[copy].Partition = matrixNumber+ imax(f1->oneBody.num,   imax(matrixNumber*f1->Na,outVector) );
        f1->sinc.tulip[copy].parallel = 1;
        f1->sinc.tulip[copy].species = matrix;
        f1->sinc.tulip[copy].header = Cube;
        f1->sinc.tulip[copy].NBody = one;
        f1->sinc.tulip[copy].purpose = tObject;
//        f1->sinc.tulip[copy].symmetryType = nullSymmetry;
        f1->sinc.tulip[copy].name = copy;
        
        f1->sinc.tulip[copyTwo].Address = fromBegining(f1,copy);
        f1->sinc.tulip[copyTwo].Partition = matrixNumber ;
        f1->sinc.tulip[copyTwo].parallel = 1;
        f1->sinc.tulip[copyTwo].species = matrix;
        f1->sinc.tulip[copyTwo].header = Cube;
        f1->sinc.tulip[copyTwo].NBody = one;
        if ( c1->rt.printFlag ){
            f1->sinc.tulip[copyTwo].NBody = two;
            f1->sinc.tulip[copyTwo].Partition = 2*c1->i.bRank * c1->i.bRank;
        }
        f1->sinc.tulip[copyTwo].purpose = tObject;
//        f1->sinc.tulip[copyTwo].symmetryType = nullSymmetry;
        f1->sinc.tulip[copyTwo].name = copyTwo;
        
        
        f1->sinc.tulip[copyThree].Address = fromBegining(f1,copyTwo);
        f1->sinc.tulip[copyThree].Partition = 0 ;
        f1->sinc.tulip[copyThree].parallel = 1;
        f1->sinc.tulip[copyThree].species = matrix;
        f1->sinc.tulip[copyThree].header = Cube;
        f1->sinc.tulip[copyThree].NBody = one;
        if ( c1->rt.printFlag ){
            f1->sinc.tulip[copyThree].NBody = two;
            f1->sinc.tulip[copyThree].Partition = c1->i.decomposeRankMatrix;
        }
        f1->sinc.tulip[copyThree].purpose = tObject;
//        f1->sinc.tulip[copyThree].symmetryType = nullSymmetry;
        f1->sinc.tulip[copyThree].name = copyThree;
        
        f1->sinc.tulip[square].Address = fromBegining(f1,copyThree);
        f1->sinc.tulip[square].Partition=   4*maxVector*maxVector;
        f1->sinc.tulip[square].species = matrix;
        f1->sinc.tulip[square].header = Cube;
        f1->sinc.tulip[square].NBody = one;
        f1->sinc.tulip[square].parallel = 0;
        f1->sinc.tulip[square].purpose = Object;
//        f1->sinc.tulip[square].symmetryType = nullSymmetry;
        f1->sinc.tulip[square].name = square;
        
        
        f1->sinc.tulip[squareTwo].spinor = none;
        f1->sinc.tulip[squareTwo].Address = fromBegining(f1,square);
        f1->sinc.tulip[squareTwo].Partition= 1+!(!c1->rt.printFlag) * c1->i.nStates* c1->i.decomposeRankMatrix;
        f1->sinc.tulip[squareTwo].species = matrix;
        f1->sinc.tulip[squareTwo].header = Cube;
        f1->sinc.tulip[squareTwo].NBody = two;
        f1->sinc.tulip[squareTwo].parallel = 0;
        f1->sinc.tulip[squareTwo].purpose = tObject;
//        f1->sinc.tulip[squareTwo].symmetryType = nullSymmetry;
        f1->sinc.tulip[squareTwo].name = squareTwo;
        
        
        for ( label = 0; label < end ; label++)
            if ( part(f1, label ) > largestRank && purpose(f1, label) == Object )
                largestRank = part(f1, label);
        
        f1->sinc.tulip[trainVector].Address = fromBegining(f1,squareTwo);
        f1->sinc.tulip[trainVector].Partition = (bootBodies == one ) *maxVector;
        f1->sinc.tulip[trainVector].parallel = 1;
        f1->sinc.tulip[trainVector].NBody = one;
        f1->sinc.tulip[trainVector].species = vector;
        f1->sinc.tulip[trainVector].header = Cube;
        f1->sinc.tulip[trainVector].purpose = tObject;
//        f1->sinc.tulip[trainVector].symmetryType = nullSymmetry;
        f1->sinc.tulip[trainVector].name = trainVector;
        
        f1->sinc.tulip[trainVector2].Address = fromBegining(f1,trainVector);
        f1->sinc.tulip[trainVector2].Partition = maxVector;
        f1->sinc.tulip[trainVector2].parallel = 1;
        f1->sinc.tulip[trainVector2].NBody = two;
        f1->sinc.tulip[trainVector2].species = vector;
        f1->sinc.tulip[trainVector2].header = Cube;
        f1->sinc.tulip[trainVector2].purpose = tObject;
//        f1->sinc.tulip[trainVector2].symmetryType = nullSymmetry;
        f1->sinc.tulip[trainVector2].name = trainVector2;
        
        f1->sinc.tulip[trainVector3].Address = fromBegining(f1,trainVector2);
        f1->sinc.tulip[trainVector3].Partition = ((bootBodies == three ) *maxVector ) + part(f1, entropyUnit) ;
        f1->sinc.tulip[trainVector3].parallel = 1;
        f1->sinc.tulip[trainVector3].NBody = three;
        f1->sinc.tulip[trainVector3].species = vector;
        f1->sinc.tulip[trainVector3].header = Cube;
        f1->sinc.tulip[trainVector3].purpose = tObject;
//        f1->sinc.tulip[trainVector3].symmetryType = nullSymmetry;
        f1->sinc.tulip[trainVector3].name = trainVector3;
        
        f1->sinc.tulip[trainVector4].Address = fromBegining(f1,trainVector3);
        f1->sinc.tulip[trainVector4].Partition = (bootBodies == four ) *maxVector;
        f1->sinc.tulip[trainVector4].parallel = 1;
        f1->sinc.tulip[trainVector4].NBody = four;
        f1->sinc.tulip[trainVector4].species = vector;
        f1->sinc.tulip[trainVector4].header = Cube;
        f1->sinc.tulip[trainVector4].purpose = tObject;
//        f1->sinc.tulip[trainVector4].symmetryType = nullSymmetry;
        f1->sinc.tulip[trainVector4].name = trainVector4;
        
        f1->sinc.tulip[trainMatrix].Address = fromBegining(f1,trainVector4);
        f1->sinc.tulip[trainMatrix].Partition = imax(matrixNumber*f1->Na,c1->i.sectors);
        f1->sinc.tulip[trainMatrix].parallel = 0;
        f1->sinc.tulip[trainMatrix].species = matrix;
        f1->sinc.tulip[trainMatrix].NBody = one;
        f1->sinc.tulip[trainMatrix].header = Cube;
        f1->sinc.tulip[trainMatrix].purpose = Object;
//        f1->sinc.tulip[trainMatrix].symmetryType = nullSymmetry;
        f1->sinc.tulip[trainMatrix].name = trainMatrix;
        
        f1->sinc.tulip[trainMatrix2].Address = fromBegining(f1,trainMatrix);
        f1->sinc.tulip[trainMatrix2].Partition =  c1->i.sectors+!(!c1->rt.printFlag)*imax(c1->i.decomposeRankMatrix,c1->i.bRank);
        f1->sinc.tulip[trainMatrix2].parallel = 0;
        f1->sinc.tulip[trainMatrix2].species = matrix;
        f1->sinc.tulip[trainMatrix2].NBody = two;
        f1->sinc.tulip[trainMatrix2].header = Cube;
        f1->sinc.tulip[trainMatrix2].purpose = Object;
//        f1->sinc.tulip[trainMatrix2].symmetryType = nullSymmetry;
        f1->sinc.tulip[trainMatrix2].name = trainMatrix2;
        
        f1->sinc.tulip[trainMatrix3].Address = fromBegining(f1,trainMatrix2);
        f1->sinc.tulip[trainMatrix3].Partition = (bootBodies == three)*c1->i.sectors;
        f1->sinc.tulip[trainMatrix3].parallel = 0;
        f1->sinc.tulip[trainMatrix3].species = matrix;
        f1->sinc.tulip[trainMatrix3].NBody = three;
        f1->sinc.tulip[trainMatrix3].header = Cube;
        f1->sinc.tulip[trainMatrix3].purpose = Object;
//        f1->sinc.tulip[trainMatrix3].symmetryType = nullSymmetry;
        f1->sinc.tulip[trainMatrix3].name = trainMatrix3;
        
        f1->sinc.tulip[trainMatrix4].Address = fromBegining(f1,trainMatrix3);
        f1->sinc.tulip[trainMatrix4].Partition = (bootBodies == four )*c1->i.sectors;
        f1->sinc.tulip[trainMatrix4].parallel = 0;
        f1->sinc.tulip[trainMatrix4].species = matrix;
        f1->sinc.tulip[trainMatrix4].NBody = four;
        f1->sinc.tulip[trainMatrix4].header = Cube;
        f1->sinc.tulip[trainMatrix4].purpose = Object;
//        f1->sinc.tulip[trainMatrix4].symmetryType = nullSymmetry;
        f1->sinc.tulip[trainMatrix4].name = trainMatrix4;
        
        f1->sinc.tulip[trainQuartic].Address = fromBegining(f1,trainMatrix4);
        f1->sinc.tulip[trainQuartic].Partition =0;
        f1->sinc.tulip[trainQuartic].parallel = 1;
        f1->sinc.tulip[trainQuartic].species = quartic;
        f1->sinc.tulip[trainQuartic].header = Cube;
        f1->sinc.tulip[trainQuartic].purpose = Object;
//        f1->sinc.tulip[trainQuartic].symmetryType = nullSymmetry;
        f1->sinc.tulip[trainQuartic].name = trainQuartic;
        
      
        f1->sinc.tulip[highBallVector].Address =  fromBegining(f1,trainQuartic);
        f1->sinc.tulip[highBallVector].Partition =maxVector*(( c1->rt.printFlag / 8 )%2);
        f1->sinc.tulip[highBallVector].species = vector;
        f1->sinc.tulip[highBallVector].NBody = four;
        f1->sinc.tulip[highBallVector].header = Cube;
        f1->sinc.tulip[highBallVector].purpose = Object;
        f1->sinc.tulip[highBallVector].name = highBallVector;
        
        f1->sinc.tulip[highBallVector2].Address =  fromBegining(f1,highBallVector);
        f1->sinc.tulip[highBallVector2].Partition =maxVector*(( c1->rt.printFlag / 8 )%2);
        f1->sinc.tulip[highBallVector2].species = vector;
        f1->sinc.tulip[highBallVector2].NBody = three;
        f1->sinc.tulip[highBallVector2].header = Cube;
        f1->sinc.tulip[highBallVector2].purpose = Object;
        f1->sinc.tulip[highBallVector2].name = highBallVector2;

        f1->sinc.tulip[highBallVector3].Address =  fromBegining(f1,highBallVector2);
        f1->sinc.tulip[highBallVector3].Partition =maxVector*(( c1->rt.printFlag / 8 )%2);
        f1->sinc.tulip[highBallVector3].species = vector;
        f1->sinc.tulip[highBallVector3].NBody = two;
        f1->sinc.tulip[highBallVector3].header = Cube;
        f1->sinc.tulip[highBallVector3].purpose = Object;
        f1->sinc.tulip[highBallVector3].name = highBallVector3;

        f1->sinc.tulip[highBallVector4].Address =  fromBegining(f1,highBallVector3);
        f1->sinc.tulip[highBallVector4].Partition =maxVector*(( c1->rt.printFlag / 8 )%2);
        f1->sinc.tulip[highBallVector4].species = vector;
        f1->sinc.tulip[highBallVector4].NBody = one;
        f1->sinc.tulip[highBallVector4].header = Cube;
        f1->sinc.tulip[highBallVector4].purpose = Object;
        f1->sinc.tulip[highBallVector4].name = highBallVector4;

        
        {
            INT_TYPE len[3];
            length(f1, eigenVectors, len);
            
            f1->sinc.tulip[canonicalBuffers].spinor = none;
            f1->sinc.tulip[canonicalBuffers].parallel = 1;
            f1->sinc.tulip[canonicalBuffers].Address = fromBegining(f1,highBallVector4);
            f1->sinc.tulip[canonicalBuffers].Partition = maxVector*maxVector+ imax(len[0],imax(NV,part(f1,totalVector)))*maxVector;
            f1->sinc.tulip[canonicalBuffers].species = scalar;
            f1->sinc.tulip[canonicalBuffers].purpose = Object;
            //        f1->sinc.tulip[canonicalBuffers].symmetryType = nullSymmetry;
            f1->sinc.tulip[canonicalBuffers].name = canonicalBuffers;
        }
        
        f1->sinc.tulip[foundationStructure].Address =  fromBegining(f1,canonicalBuffers);
        f1->sinc.tulip[foundationStructure].Partition = nG;
        f1->sinc.tulip[foundationStructure].NBody = bootBodies;
        f1->sinc.tulip[foundationStructure].parallel = 1;
        f1->sinc.tulip[foundationStructure].species = vector;
        f1->sinc.tulip[foundationStructure].header = Cube;
        f1->sinc.tulip[foundationStructure].purpose = Object;
        f1->sinc.tulip[foundationStructure].name = foundationStructure;
        
        f1->sinc.tulip[foundationEquals].Address =  fromBegining(f1,foundationStructure);
        f1->sinc.tulip[foundationEquals].Partition = 1;
        f1->sinc.tulip[foundationEquals].parallel = 1;
        f1->sinc.tulip[foundationEquals].species = vector;
        f1->sinc.tulip[foundationEquals].header = Cube;
        f1->sinc.tulip[foundationEquals].purpose = Object;
        f1->sinc.tulip[foundationEquals].name = foundationEquals;

        {
            f1->sinc.tulip[interactionDirect].Address = fromBegining(f1,foundationEquals);
            f1->sinc.tulip[interactionDirect].spinor = none;
            f1->sinc.tulip[interactionDirect].Partition = (!(!c1->i.hartreeFockFlag))*f1->twoBody.num*( bootBodies > one );
            f1->sinc.tulip[interactionDirect].species = matrix;
            f1->sinc.tulip[interactionDirect].NBody = two;
            f1->sinc.tulip[interactionDirect].header = Cube;
            f1->sinc.tulip[interactionDirect].purpose = Object;
            f1->sinc.tulip[interactionDirect].name = interactionDirect;
            
            f1->sinc.tulip[interactionExchange].Address = fromBegining(f1,interactionDirect);
            f1->sinc.tulip[interactionExchange].spinor = none;
            f1->sinc.tulip[interactionExchange].Partition = f1->twoBody.num*( bootBodies > one );
            f1->sinc.tulip[interactionExchange].species = matrix;
            f1->sinc.tulip[interactionExchange].NBody = two;
            f1->sinc.tulip[interactionExchange].header = Cube;
            f1->sinc.tulip[interactionExchange].purpose = Object;
            f1->sinc.tulip[interactionExchange].name = interactionExchange;
            
            f1->sinc.tulip[interactionExchangePlus].Address = fromBegining(f1,interactionExchange);
            f1->sinc.tulip[interactionExchangePlus].spinor = none;
            f1->sinc.tulip[interactionExchangePlus].Partition = f1->twoBody.num*( bootBodies >= two && bootType > electron );
            f1->sinc.tulip[interactionExchangePlus].species = matrix;
            f1->sinc.tulip[interactionExchangePlus].NBody = two;
            f1->sinc.tulip[interactionExchangePlus].header = Cube;
            f1->sinc.tulip[interactionExchangePlus].purpose = Object;
            f1->sinc.tulip[interactionExchangePlus].name = interactionExchangePlus;

            
            f1->sinc.tulip[quadCube].spinor = none;
            f1->sinc.tulip[quadCube].Address = fromBegining(f1,interactionExchangePlus);
            f1->sinc.tulip[quadCube].Partition =1;
            f1->sinc.tulip[quadCube].species = matrix;
            f1->sinc.tulip[quadCube].header = Cube;
            f1->sinc.tulip[quadCube].NBody = two;
            f1->sinc.tulip[quadCube].purpose = Object;
            f1->sinc.tulip[quadCube].name = quadCube;
            
            f1->sinc.tulip[quad2Cube].spinor = none;
            f1->sinc.tulip[quad2Cube].Address = fromBegining(f1,quadCube);
            f1->sinc.tulip[quad2Cube].Partition =0;
            f1->sinc.tulip[quad2Cube].species = matrix;
            f1->sinc.tulip[quad2Cube].header = Cube;
            f1->sinc.tulip[quad2Cube].NBody = two;
            f1->sinc.tulip[quad2Cube].purpose = Object;
            f1->sinc.tulip[quad2Cube].name = quad2Cube;
            
            f1->sinc.tulip[quad].spinor = none;
            f1->sinc.tulip[quad].Address = fromBegining(f1,quad2Cube);
            f1->sinc.tulip[quad].Partition = 0;
            f1->sinc.tulip[quad].species = matrix;
            f1->sinc.tulip[quad].header = bootShape;
            f1->sinc.tulip[quad].NBody = two;
            f1->sinc.tulip[quad].purpose = Object;
            f1->sinc.tulip[quad].name = quad;
            
            f1->sinc.tulip[end].Partition = 0;
            f1->sinc.tulip[end].Address = fromBegining(f1,quad);
            f1->sinc.tulip[end].species = scalar;
            f1->sinc.tulip[end].purpose = Object;
            f1->sinc.tulip[end].name = end;
            
            f1->sinc.tulip[basisBuffers].parallel = 0;        //for parallel basis transforming
            f1->sinc.tulip[basisBuffers].species = matrix;
            if ( bootBodies == one )
                f1->sinc.tulip[basisBuffers].NBody = two;//for one body matrices...
            
            f1->sinc.tulip[basisBuffers].myAddress = 0;
            f1->sinc.tulip[basisBuffers].Partition = 0;
            f1->sinc.tulip[basisBuffers].header = Cube;
            f1->sinc.tulip[basisBuffers].memory = oneObject;
            f1->sinc.tulip[basisBuffers].name = basisBuffers;
        }
    
        f1->sinc.tulip[oneArray].parallel = 0;
        f1->sinc.tulip[oneArray].myAddress = fromMyBegining(f1,basisBuffers);
        f1->sinc.tulip[oneArray].Partition = c1->i.M1*3+c1->i.M1*c1->i.M1;
        f1->sinc.tulip[oneArray].species = scalar;
        f1->sinc.tulip[oneArray].memory = oneObject;
        f1->sinc.tulip[oneArray].name = oneArray;

        f1->sinc.tulip[threeArray].parallel = 0;
        f1->sinc.tulip[threeArray].myAddress = fromMyBegining(f1,oneArray);
        f1->sinc.tulip[threeArray].Partition = 2*c1->i.M1*c1->i.M1*c1->i.M1;
        f1->sinc.tulip[threeArray].species = scalar;
        f1->sinc.tulip[threeArray].memory = oneObject;
        f1->sinc.tulip[threeArray].name = threeArray;

        f1->sinc.tulip[oneBasis].parallel = 0;
        f1->sinc.tulip[oneBasis].myAddress = fromMyBegining(f1,threeArray);
        f1->sinc.tulip[oneBasis].Partition = c1->i.M1*N1;
        f1->sinc.tulip[oneBasis].species = scalar;
        f1->sinc.tulip[oneBasis].memory = oneObject;
        f1->sinc.tulip[oneBasis].name = oneBasis;
        
        f1->sinc.tulip[tensorBuffers].spinor = none;
        f1->sinc.tulip[tensorBuffers].myAddress =  fromMyBegining(f1,oneBasis);
        f1->sinc.tulip[tensorBuffers].Partition = 1;
        f1->sinc.tulip[tensorBuffers].parallel = 1;
        f1->sinc.tulip[tensorBuffers].species = vector;
        f1->sinc.tulip[tensorBuffers].header = Cube;
        f1->sinc.tulip[tensorBuffers].purpose = Object;
        f1->sinc.tulip[tensorBuffers].memory = oneObject;
        f1->sinc.tulip[tensorBuffers].name = tensorBuffers;
        
        f1->sinc.tulip[tensorBuffers2].spinor = none;
        f1->sinc.tulip[tensorBuffers2].myAddress =  fromMyBegining(f1,tensorBuffers);
        f1->sinc.tulip[tensorBuffers2].Partition = 1;
        f1->sinc.tulip[tensorBuffers2].parallel = 1;
        f1->sinc.tulip[tensorBuffers2].species = vector;
        f1->sinc.tulip[tensorBuffers2].header = Cube;
        f1->sinc.tulip[tensorBuffers2].purpose = Object;
        f1->sinc.tulip[tensorBuffers2].memory = oneObject;
        f1->sinc.tulip[tensorBuffers2].name = tensorBuffers2;

        f1->sinc.tulip[tensorBuffers3].spinor = none;
        f1->sinc.tulip[tensorBuffers3].myAddress =  fromMyBegining(f1,tensorBuffers2);
        f1->sinc.tulip[tensorBuffers3].Partition = 1;
        f1->sinc.tulip[tensorBuffers3].parallel = 1;
        f1->sinc.tulip[tensorBuffers3].species = vector;
        f1->sinc.tulip[tensorBuffers3].header = Cube;
        f1->sinc.tulip[tensorBuffers3].purpose = Object;
        f1->sinc.tulip[tensorBuffers3].memory = oneObject;
        f1->sinc.tulip[tensorBuffers3].name = tensorBuffers3;

        f1->sinc.tulip[tensorBuffers4].spinor = none;
        f1->sinc.tulip[tensorBuffers4].myAddress =  fromMyBegining(f1,tensorBuffers3);
        f1->sinc.tulip[tensorBuffers4].Partition = 1;
        f1->sinc.tulip[tensorBuffers4].parallel = 1;
        f1->sinc.tulip[tensorBuffers4].species = vector;
        f1->sinc.tulip[tensorBuffers4].header = Cube;
        f1->sinc.tulip[tensorBuffers4].purpose = Object;
        f1->sinc.tulip[tensorBuffers4].memory = oneObject;
        f1->sinc.tulip[tensorBuffers4].name = tensorBuffers4;

        f1->sinc.tulip[tensorBuffers5].spinor = none;
        f1->sinc.tulip[tensorBuffers5].myAddress =  fromMyBegining(f1,tensorBuffers4);
        f1->sinc.tulip[tensorBuffers5].Partition = 1;
        f1->sinc.tulip[tensorBuffers5].parallel = 1;
        f1->sinc.tulip[tensorBuffers5].species = vector;
        f1->sinc.tulip[tensorBuffers5].header = Cube;
        f1->sinc.tulip[tensorBuffers5].purpose = Object;
        f1->sinc.tulip[tensorBuffers5].memory = oneObject;
        f1->sinc.tulip[tensorBuffers5].name = tensorBuffers5;
        
        f1->sinc.tulip[tensorBuffers6].spinor = none;
        f1->sinc.tulip[tensorBuffers6].myAddress =  fromMyBegining(f1,tensorBuffers5);
        f1->sinc.tulip[tensorBuffers6].Partition = 1;
        f1->sinc.tulip[tensorBuffers6].parallel = 1;
        f1->sinc.tulip[tensorBuffers6].species = vector;
        f1->sinc.tulip[tensorBuffers6].header = Cube;
        f1->sinc.tulip[tensorBuffers6].purpose = Object;
        f1->sinc.tulip[tensorBuffers6].memory = oneObject;
        f1->sinc.tulip[tensorBuffers6].name = tensorBuffers6;
        
        f1->sinc.tulip[canonicalBuffersB].spinor = none;
        f1->sinc.tulip[canonicalBuffersB].parallel = 1;
        f1->sinc.tulip[canonicalBuffersB].myAddress = fromMyBegining(f1,tensorBuffers6);
        f1->sinc.tulip[canonicalBuffersB].Partition = c1->i.canonRank*imax(1,maxVector);
        f1->sinc.tulip[canonicalBuffersB].species = vector;
        f1->sinc.tulip[canonicalBuffersB].memory = oneObject;
        f1->sinc.tulip[canonicalBuffersB].name = canonicalBuffersB;
        
        f1->sinc.tulip[canonicalBuffersBX].spinor = none;
        f1->sinc.tulip[canonicalBuffersBX].parallel = 0;
        f1->sinc.tulip[canonicalBuffersBX].myAddress = fromMyBegining(f1,canonicalBuffersB);
        f1->sinc.tulip[canonicalBuffersBX].Partition = 0*part(f1, entropyUnit);
        f1->sinc.tulip[canonicalBuffersBX].NBody = bodies (f1, entropyVector);
        f1->sinc.tulip[canonicalBuffersBX].species = vector;
        f1->sinc.tulip[canonicalBuffersBX].memory = oneObject;
        f1->sinc.tulip[canonicalBuffersBX].name = canonicalBuffersBX;
        
        f1->sinc.tulip[canonicalBuffersBM].spinor = none;
        f1->sinc.tulip[canonicalBuffersBM].parallel = 0;
        f1->sinc.tulip[canonicalBuffersBM].myAddress = fromMyBegining(f1,canonicalBuffersBX);
        f1->sinc.tulip[canonicalBuffersBM].Partition = 1;
        f1->sinc.tulip[canonicalBuffersBM].species = matrix;
        f1->sinc.tulip[canonicalBuffersBM].memory = oneObject;
        f1->sinc.tulip[canonicalBuffersBM].name = canonicalBuffersBM;
        
        f1->sinc.tulip[canonicalBuffersC].spinor = none;
        f1->sinc.tulip[canonicalBuffersC].parallel = 1;
        f1->sinc.tulip[canonicalBuffersC].myAddress = fromMyBegining(f1,canonicalBuffersBM);
        f1->sinc.tulip[canonicalBuffersC].Partition = NV;
        f1->sinc.tulip[canonicalBuffersC].species = scalar;
        f1->sinc.tulip[canonicalBuffersC].memory = oneObject;
        f1->sinc.tulip[canonicalBuffersC].name = canonicalBuffersC;
        
        f1->sinc.tulip[canonicalBuffersD].spinor = none;
        f1->sinc.tulip[canonicalBuffersD].parallel = 1;
        f1->sinc.tulip[canonicalBuffersD].myAddress = fromMyBegining(f1,canonicalBuffersC);
        f1->sinc.tulip[canonicalBuffersD].Partition = 0*NV*c1->i.bRank;
        f1->sinc.tulip[canonicalBuffersD].species = scalar;
        f1->sinc.tulip[canonicalBuffersD].memory = oneObject;
        f1->sinc.tulip[canonicalBuffersD].name = canonicalBuffersD;
        
        f1->sinc.tulip[twoBodyRitz].spinor = none;
        f1->sinc.tulip[twoBodyRitz].myAddress = fromMyBegining(f1,canonicalBuffersD);
        f1->sinc.tulip[twoBodyRitz].Partition = maxArray;
        f1->sinc.tulip[twoBodyRitz].parallel = 0;
        f1->sinc.tulip[twoBodyRitz].species = scalar;
        f1->sinc.tulip[twoBodyRitz].header = Cube;
        f1->sinc.tulip[twoBodyRitz].memory = oneObject;
        f1->sinc.tulip[twoBodyRitz].name = twoBodyRitz;
        
        f1->sinc.tulip[conditionOverlapNumbers].spinor = none;
        f1->sinc.tulip[conditionOverlapNumbers].myAddress = fromMyBegining(f1,twoBodyRitz);
        f1->sinc.tulip[conditionOverlapNumbers].Partition = maxArray;
        f1->sinc.tulip[conditionOverlapNumbers].parallel = 0;
        f1->sinc.tulip[conditionOverlapNumbers].species = scalar;
        f1->sinc.tulip[conditionOverlapNumbers].header = Cube;
        f1->sinc.tulip[conditionOverlapNumbers].memory = oneObject;
        f1->sinc.tulip[conditionOverlapNumbers].name = conditionOverlapNumbers;

        f1->sinc.tulip[twoBodyProjector].spinor = none;
        f1->sinc.tulip[twoBodyProjector].myAddress = fromMyBegining(f1,conditionOverlapNumbers);
        f1->sinc.tulip[twoBodyProjector].Partition = 0*maxEV;
        f1->sinc.tulip[twoBodyProjector].parallel = 0;
        f1->sinc.tulip[twoBodyProjector].species = scalar;
        f1->sinc.tulip[twoBodyProjector].header = Cube;
        f1->sinc.tulip[twoBodyProjector].memory = oneObject;
        f1->sinc.tulip[twoBodyProjector].name = twoBodyProjector;
        
        f1->sinc.tulip[matrixHbuild].spinor = none;
        f1->sinc.tulip[matrixHbuild].myAddress = fromMyBegining(f1,twoBodyProjector);
        f1->sinc.tulip[matrixHbuild].Partition = imax(6*maxArray*maxArray, 2*vectorLen(f1, eigenVectors)[0]*vectorLen(f1, eigenVectors)[0]  );
        f1->sinc.tulip[matrixHbuild].parallel = 0;
        f1->sinc.tulip[matrixHbuild].species = scalar;
        f1->sinc.tulip[matrixHbuild].header = Cube;
        f1->sinc.tulip[matrixHbuild].memory = oneObject;
        f1->sinc.tulip[matrixHbuild].name = matrixHbuild;
        
        f1->sinc.tulip[vectorHbuild].spinor = none;
        f1->sinc.tulip[vectorHbuild].myAddress = fromMyBegining(f1,matrixHbuild);
        f1->sinc.tulip[vectorHbuild].Partition = maxArray;
        f1->sinc.tulip[vectorHbuild].parallel = 0;
        f1->sinc.tulip[vectorHbuild].species = scalar;
        f1->sinc.tulip[vectorHbuild].header = Cube;
        f1->sinc.tulip[vectorHbuild].memory = oneObject;
        f1->sinc.tulip[vectorHbuild].name = vectorHbuild;
        
        
        f1->sinc.tulip[matrixSbuild].spinor = none;
        f1->sinc.tulip[matrixSbuild].myAddress = fromMyBegining(f1,vectorHbuild);
        f1->sinc.tulip[matrixSbuild].Partition = 4*maxArray*maxArray;
        f1->sinc.tulip[matrixSbuild].parallel = 0;
        f1->sinc.tulip[matrixSbuild].species = scalar;
        f1->sinc.tulip[matrixSbuild].header = Cube;
        f1->sinc.tulip[matrixSbuild].memory = oneObject;
        f1->sinc.tulip[matrixSbuild].name = matrixSbuild;
        
        f1->sinc.tulip[dsyBuffers].spinor = none;
        f1->sinc.tulip[dsyBuffers].myAddress = fromMyBegining(f1,matrixSbuild);
#ifdef APPLE
        f1->sinc.tulip[dsyBuffers].Partition = 8*(8*(imax(vectorLen(f1,eigenVectors)[0],maxEV))+72*c1->i.nStates*c1->i.nStates+ 8 * vectorLen(f1, squareVector)[0])+3*maxEV;
#else
        f1->sinc.tulip[dsyBuffers].Partition = maxVector*maxVector;
#endif
        f1->sinc.tulip[dsyBuffers].parallel = 1;
        f1->sinc.tulip[dsyBuffers].species = scalar;
        f1->sinc.tulip[dsyBuffers].header = Cube;
        f1->sinc.tulip[dsyBuffers].memory = oneObject;
        f1->sinc.tulip[dsyBuffers].name = dsyBuffers;
        
        f1->sinc.tulip[end].myAddress = fromMyBegining(f1,dsyBuffers);

        if ( (3*fromBegining(f1,end )+fromMyBegining(f1,end ))/(1000000000./(sizeof(Stream_Type))) > c1->i.RAMmax ){
            printf("oops too much RAM required\n");
            printf("\n\tG%f",(3*fromBegining(f1,end ))/(1000000000./(sizeof(Stream_Type))));
            printf("\tG%f\n",(fromMyBegining(f1,end ))/(1000000000./(sizeof(Stream_Type))));
            exit(0);
        }
        printf("\n\tG%f",(3*fromBegining(f1,end ))/(1000000000./(sizeof(Stream_Type))));
        printf("\tG%f\n",(fromMyBegining(f1,end ))/(1000000000./(sizeof(Stream_Type))));

        
        
        fflush(stdout);
        for ( space = 0; space < SPACE ; space++){
            f1->sinc.rose[space].stream = malloc( fromBegining(f1,end )*sizeof(Stream_Type));
        }

        f1->sinc.rose[SPACE].stream = malloc( fromMyBegining(f1,end )*sizeof(Stream_Type));
        c1->mem.bootedMemory = 1;
        assignCores(f1, 2);

        tClear(f1, linear);

        INT_TYPE RdsSize;
        RdsSize = 0;

        
        {
            INT_TYPE i,ii;
            INT_TYPE M1 = c1->i.M1,M12 = (M1-1)/2;
            double r = (double)(M1-1)/(double)(N1-1);
            double * u =myStreams(f1,oneBasis,0);
            for( i = 0; i < N1 ; i++)
                for ( ii = 0 ; ii < M1 ; ii++){
                    u[i*M1+ii] = Sinc(r, r*(i-N12)- (ii-M12))/sqrt(r);
                }
        }
        
        
        tEnd(f1, edgeMatrix, 0, 0);
        tEnd(f1, edgeMatrix, 0, 1);
        tEnd(f1, edgeMatrix, 0, 2);
        tAlt(f1, edgeMatrix, 0, 0);
        tAlt(f1, edgeMatrix, 0, 1);
        tAlt(f1, edgeMatrix, 0, 2);

        
        {
            
            if ( bootType == electron  ){
                if ( bootBodies > one ){
                    mySeparateExactTwo(f1,c1->rt.runFlag, 1. , 0,0,1,2);
                    if ( (c1->rt.runFlag == 7 && SPACE ==3 ) || (c1->rt.runFlag == 3 && SPACE == 2 )){
                        
                        buildElectronProtonInteraction(f1, linear,0);
                        tZeroSum(f1, interactionExchange, 0 );
                        tZeroSum(f1, linear, 0 );

                    }
                }
                if (f1->Na != 0 ){
                    separateExternal(c1,c1->rt.runFlag,0,1.0,4,0,1);
                }
                if ( c1->i.springFlag ){
                    double vspr[3];
                    vspr[0] = c1->i.springFlag;
                    vspr[1] = c1->i.springFlag;
                    vspr[2] = c1->i.springFlag;

                    separateHarmonicExternal(c1,c1->rt.runFlag,1.,vspr,0,1);
                }
                separateKinetic(f1, c1->rt.runFlag,kinetic, 1.0,1);
            
            } else if ( bootType == h2plus ){
                separateExternal(c1,c1->rt.runFlag,0,1.0,4,0,1);
                mySeparateExactTwo(f1,c1->rt.runFlag, 1. , 0,0,1,2);
                mySeparateExactTwo(f1,c1->rt.runFlag, 1. , 0,1,1,2);
                separateKinetic(f1, c1->rt.runFlag,kineticMass, 1836.15267245/8.  ,1);
                separateKinetic(f1, c1->rt.runFlag,kinetic, 0.9997277656208493,2);
            }
        }
    
        {
            
            if ( c1->i.densityFlag ){
                f1->sinc.tulip[Ha].linkNext = density;
                f1->sinc.tulip[Iterator].linkNext = density;
            }
            else
            if ( bootBodies == one && bootType == electron){
                
                //FOR LINEAR
                f1->sinc.tulip[external1].spinor = none;
                f1->sinc.tulip[external1].NBody = one;
                f1->sinc.tulip[external1].species = matrix;
                f1->sinc.tulip[external1].blockType = tv1;
                f1->sinc.tulip[external1].name = linear;
                f1->sinc.tulip[external1].purpose = ptObject;
                f1->sinc.tulip[external1].Current[0] = part(f1, linear);
                f1->sinc.tulip[external1].ptRank[0] = 0;
                
                if ( c1->i.springFlag )
                    f1->sinc.tulip[kinetic1].linkNext = harmonium1;
                else
                    f1->sinc.tulip[kinetic1].linkNext = external1;
                f1->sinc.tulip[kinetic1].NBody = one;
                f1->sinc.tulip[kinetic1].species = matrix;
                f1->sinc.tulip[kinetic1].blockType = tv1;
                f1->sinc.tulip[kinetic1].name = kinetic;
                f1->sinc.tulip[kinetic1].purpose = ptObject;
                f1->sinc.tulip[kinetic1].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic1].ptRank[0] = 0;
                f1->sinc.tulip[kinetic1].Current[1] = 0;
                f1->sinc.tulip[kinetic1].ptRank[1] = 0;

                
                
                f1->sinc.tulip[harmonium1].spinor = none;
                f1->sinc.tulip[harmonium1].linkNext = external1;
                f1->sinc.tulip[harmonium1].NBody = one;
                f1->sinc.tulip[harmonium1].species = matrix;
                f1->sinc.tulip[harmonium1].blockType = tv1;
                f1->sinc.tulip[harmonium1].name = harmonium;
                f1->sinc.tulip[harmonium1].purpose = ptObject;
                f1->sinc.tulip[harmonium1].Current[0] = part(f1, harmonium1);
                f1->sinc.tulip[harmonium1].ptRank[0] = 0;
                
                f1->sinc.tulip[X1].spinor = none;
                f1->sinc.tulip[X1].linkNext = kinetic1;
                f1->sinc.tulip[X1].NBody = one;
                f1->sinc.tulip[X1].species = matrix;
                f1->sinc.tulip[X1].blockType = tv1;
                f1->sinc.tulip[X1].name = X;
                f1->sinc.tulip[X1].purpose = ptObject;
                f1->sinc.tulip[X1].Current[0] = 1;
                f1->sinc.tulip[X1].ptRank[0] = 0;

                //active assignment
                    f1->sinc.tulip[Ha].linkNext = kinetic1;
                
                    f1->sinc.tulip[Iterator].linkNext = kinetic1;
                
                
            } else if ( bootBodies == two && bootType == electron){
                
                
                
                f1->sinc.tulip[interaction12].species = matrix;
                f1->sinc.tulip[interaction12].name = interactionExchange;
                if ( c1->i.hartreeFockFlag )
                    f1->sinc.tulip[interaction12].linkNext = hartree;
                f1->sinc.tulip[interaction12].blockType = e12;
                f1->sinc.tulip[interaction12].NBody = two;
               // f1->sinc.tulip[interaction12].memory = threeObject;
                f1->sinc.tulip[interaction12].Current[0] = CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction12].ptRank[0] = 0;
                f1->sinc.tulip[interaction12].header = Cube;
                
                f1->sinc.tulip[hartree].species = matrix;
                f1->sinc.tulip[hartree].name = interactionDirect;
                f1->sinc.tulip[hartree].blockType = e12;
                f1->sinc.tulip[hartree].NBody = two;
               // f1->sinc.tulip[hartree].memory = threeObject;
                f1->sinc.tulip[hartree].Current[0] = (!(!c1->i.hartreeFockFlag))*CanonicalRank(f1, interactionDirect,0);//only run if SET BY ExactOne...
                f1->sinc.tulip[hartree].ptRank[0] = 0;
                f1->sinc.tulip[hartree].header = Cube;
                
                
                f1->sinc.tulip[Ha].Partition = 0;//
                f1->sinc.tulip[Ha].species = matrix;
                f1->sinc.tulip[Ha].header = bootShape;
                f1->sinc.tulip[Ha].name = Ha;
                
                
                f1->sinc.tulip[harmonium1].spinor = none;
                f1->sinc.tulip[harmonium1].linkNext = harmonium2;
                f1->sinc.tulip[harmonium1].NBody = one;
                f1->sinc.tulip[harmonium1].species = matrix;
                f1->sinc.tulip[harmonium1].blockType = tv1;
                f1->sinc.tulip[harmonium1].name = harmonium;
                f1->sinc.tulip[harmonium1].purpose = ptObject;
                f1->sinc.tulip[harmonium1].Current[0] = part(f1, harmonium1);
                f1->sinc.tulip[harmonium1].ptRank[0] = 0;
                
                
                f1->sinc.tulip[harmonium2].spinor = none;
                f1->sinc.tulip[harmonium2].linkNext = external1;
                f1->sinc.tulip[harmonium2].NBody = one;
                f1->sinc.tulip[harmonium2].species = matrix;
                f1->sinc.tulip[harmonium2].blockType = tv2;
                f1->sinc.tulip[harmonium2].name = harmonium;
                f1->sinc.tulip[harmonium2].purpose = ptObject;
                f1->sinc.tulip[harmonium2].Current[0] = part(f1, harmonium2);
                f1->sinc.tulip[harmonium2].ptRank[0] = 0;
                
                
                
                
                //FOR LINEAR
                f1->sinc.tulip[external1].spinor = none;
                f1->sinc.tulip[external1].linkNext = external2;
                f1->sinc.tulip[external1].NBody = one;
                f1->sinc.tulip[external1].species = matrix;
                f1->sinc.tulip[external1].blockType = tv1;
                f1->sinc.tulip[external1].name = linear;
                f1->sinc.tulip[external1].purpose = ptObject;
                f1->sinc.tulip[external1].Current[0] = part(f1, linear);
                f1->sinc.tulip[external1].ptRank[0] = 0;
                
                f1->sinc.tulip[kinetic1].linkNext = kinetic2;
                f1->sinc.tulip[kinetic1].NBody = one;
                f1->sinc.tulip[kinetic1].species = matrix;
                f1->sinc.tulip[kinetic1].blockType = tv1;
                f1->sinc.tulip[kinetic1].name = kinetic;
                f1->sinc.tulip[kinetic1].purpose = ptObject;
                f1->sinc.tulip[kinetic1].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic1].ptRank[0] = 0;
                f1->sinc.tulip[kinetic1].Current[1] = 0;
                f1->sinc.tulip[kinetic1].ptRank[1] = 0;

                //FOR LINEAR
                f1->sinc.tulip[external2].spinor = none;
                f1->sinc.tulip[external2].linkNext = interaction12;
                f1->sinc.tulip[external2].NBody = one;
                f1->sinc.tulip[external2].species = matrix;
                f1->sinc.tulip[external2].name = linear;
                f1->sinc.tulip[external2].blockType = tv2;
                f1->sinc.tulip[external2].purpose = ptObject;
                f1->sinc.tulip[external2].Current[0] = part(f1, linear);
                f1->sinc.tulip[external2].ptRank[0] = 0;
                
                if ( c1->i.springFlag )
                    f1->sinc.tulip[kinetic2].linkNext = harmonium1;
                else
                    f1->sinc.tulip[kinetic2].linkNext = external1;
                f1->sinc.tulip[kinetic2].NBody = one;
                f1->sinc.tulip[kinetic2].species = matrix;
                f1->sinc.tulip[kinetic2].blockType = tv2;
                f1->sinc.tulip[kinetic2].name = kinetic;
                f1->sinc.tulip[kinetic2].purpose = ptObject;
                f1->sinc.tulip[kinetic2].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic2].ptRank[0] = 0;
                f1->sinc.tulip[kinetic2].Current[1] = 0;
                f1->sinc.tulip[kinetic2].ptRank[1] = 0;

                f1->sinc.tulip[X1].spinor = none;
                f1->sinc.tulip[X1].linkNext = X2;
                f1->sinc.tulip[X1].NBody = one;
                f1->sinc.tulip[X1].species = matrix;
                f1->sinc.tulip[X1].blockType = tv1;
                f1->sinc.tulip[X1].name = X;
                f1->sinc.tulip[X1].purpose = ptObject;
                f1->sinc.tulip[X1].Current[0] = 1;
                f1->sinc.tulip[X1].ptRank[0] = 0;

                f1->sinc.tulip[X2].spinor = none;
                f1->sinc.tulip[X2].linkNext = kinetic1;
                f1->sinc.tulip[X2].NBody = one;
                f1->sinc.tulip[X2].species = matrix;
                f1->sinc.tulip[X2].blockType = tv2;
                f1->sinc.tulip[X2].name = X;
                f1->sinc.tulip[X2].purpose = ptObject;
                f1->sinc.tulip[X2].Current[0] = 1;
                f1->sinc.tulip[X2].ptRank[0] = 0;

                
                //active assignment
                    f1->sinc.tulip[Ha].linkNext = kinetic1;
                
                    f1->sinc.tulip[Iterator].linkNext = kinetic1;
                
            }else if ( bootBodies == three && bootType == electron){
                
                f1->sinc.tulip[interaction12].species = matrix;
                f1->sinc.tulip[interaction12].name = interactionExchange;
                f1->sinc.tulip[interaction12].blockType = e12;
                f1->sinc.tulip[interaction12].NBody = two;
                f1->sinc.tulip[interaction12].Current[0] =CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction12].ptRank[0] = 0;
                f1->sinc.tulip[interaction12].linkNext = interaction13;
                f1->sinc.tulip[interaction12].header = Cube;
                
                f1->sinc.tulip[interaction13].species = matrix;
                f1->sinc.tulip[interaction13].name = interactionExchange;
                f1->sinc.tulip[interaction13].blockType = e13;
                f1->sinc.tulip[interaction13].NBody = two;
                f1->sinc.tulip[interaction13].Current[0] =CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction13].ptRank[0] = 0;
                f1->sinc.tulip[interaction13].linkNext = interaction23;
                f1->sinc.tulip[interaction13].header = Cube;
                
                f1->sinc.tulip[interaction23].species = matrix;
                f1->sinc.tulip[interaction23].name = interactionExchange;
                f1->sinc.tulip[interaction23].blockType = e23;
                f1->sinc.tulip[interaction23].NBody = two;
                f1->sinc.tulip[interaction23].Current[0] = CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction23].ptRank[0] = 0;
                f1->sinc.tulip[interaction23].header = Cube;
                
                f1->sinc.tulip[harmonium1].spinor = none;
                f1->sinc.tulip[harmonium1].linkNext = harmonium2;
                f1->sinc.tulip[harmonium1].NBody = one;
                f1->sinc.tulip[harmonium1].species = matrix;
                f1->sinc.tulip[harmonium1].blockType = tv1;
                f1->sinc.tulip[harmonium1].name = harmonium;
                f1->sinc.tulip[harmonium1].purpose = ptObject;
                f1->sinc.tulip[harmonium1].Current[0] = part(f1, harmonium1);
                f1->sinc.tulip[harmonium1].ptRank[0] = 0;
                
                f1->sinc.tulip[harmonium2].spinor = none;
                f1->sinc.tulip[harmonium2].linkNext = harmonium3;
                f1->sinc.tulip[harmonium2].NBody = one;
                f1->sinc.tulip[harmonium2].species = matrix;
                f1->sinc.tulip[harmonium2].blockType = tv2;
                f1->sinc.tulip[harmonium2].name = harmonium;
                f1->sinc.tulip[harmonium2].purpose = ptObject;
                f1->sinc.tulip[harmonium2].Current[0] = part(f1, harmonium2);
                f1->sinc.tulip[harmonium2].ptRank[0] = 0;
                
                f1->sinc.tulip[harmonium3].spinor = none;
                f1->sinc.tulip[harmonium3].linkNext = external1;
                f1->sinc.tulip[harmonium3].NBody = one;
                f1->sinc.tulip[harmonium3].species = matrix;
                f1->sinc.tulip[harmonium3].blockType = tv3;
                f1->sinc.tulip[harmonium3].name = harmonium;
                f1->sinc.tulip[harmonium3].purpose = ptObject;
                f1->sinc.tulip[harmonium3].Current[0] = part(f1, harmonium2);
                f1->sinc.tulip[harmonium3].ptRank[0] = 0;
                
                f1->sinc.tulip[Ha].Partition = 0;//
                f1->sinc.tulip[Ha].linkNext = external1   ;
                f1->sinc.tulip[Ha].species = matrix;
                f1->sinc.tulip[Ha].header = bootShape;
                f1->sinc.tulip[Ha].name = Ha;
                
                
                //FOR LINEAR
                f1->sinc.tulip[external1].spinor = none;
                f1->sinc.tulip[external1].linkNext = external2;
                f1->sinc.tulip[external1].NBody = one;
                f1->sinc.tulip[external1].species = matrix;
                f1->sinc.tulip[external1].blockType = tv1;
                f1->sinc.tulip[external1].name = linear;
                f1->sinc.tulip[external1].purpose = ptObject;
                f1->sinc.tulip[external1].Current[0] = part(f1, linear);
                f1->sinc.tulip[external1].ptRank[0] = 0;
                
                f1->sinc.tulip[kinetic1].linkNext = kinetic2;
                f1->sinc.tulip[kinetic1].NBody = one;
                f1->sinc.tulip[kinetic1].species = matrix;
                f1->sinc.tulip[kinetic1].blockType = tv1;
                f1->sinc.tulip[kinetic1].name = kinetic;
                f1->sinc.tulip[kinetic1].purpose = ptObject;
                f1->sinc.tulip[kinetic1].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic1].ptRank[0] = 0;
                f1->sinc.tulip[kinetic1].Current[1] = 0;
                f1->sinc.tulip[kinetic1].ptRank[1] = 0;

                //FOR LINEAR
                f1->sinc.tulip[external2].spinor = none;
                f1->sinc.tulip[external2].linkNext = external3;
                f1->sinc.tulip[external2].NBody = one;
                f1->sinc.tulip[external2].species = matrix;
                f1->sinc.tulip[external2].name = linear;
                f1->sinc.tulip[external2].blockType = tv2;
                f1->sinc.tulip[external2].purpose = ptObject;
                f1->sinc.tulip[external2].Current[0] = part(f1, linear);
                f1->sinc.tulip[external2].ptRank[0] = 0;
                
                f1->sinc.tulip[kinetic2].linkNext = kinetic3;
                f1->sinc.tulip[kinetic2].NBody = one;
                f1->sinc.tulip[kinetic2].species = matrix;
                f1->sinc.tulip[kinetic2].blockType = tv2;
                f1->sinc.tulip[kinetic2].name = kinetic;
                f1->sinc.tulip[kinetic2].purpose = ptObject;
                f1->sinc.tulip[kinetic2].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic2].ptRank[0] = 0;
                f1->sinc.tulip[kinetic2].Current[1] = 0;
                f1->sinc.tulip[kinetic2].ptRank[1] = 0;

                f1->sinc.tulip[external3].spinor = none;
                f1->sinc.tulip[external3].linkNext = interaction12;
                f1->sinc.tulip[external3].NBody = one;
                f1->sinc.tulip[external3].species = matrix;
                f1->sinc.tulip[external3].name = linear;
                f1->sinc.tulip[external3].blockType = tv3;
                f1->sinc.tulip[external3].purpose = ptObject;
                f1->sinc.tulip[external3].Current[0] = part(f1, linear);
                f1->sinc.tulip[external3].ptRank[0] = 0;
                
                f1->sinc.tulip[kinetic3].NBody = one;
                if ( c1->i.springFlag )
                    f1->sinc.tulip[kinetic3].linkNext = harmonium1;
                else
                    f1->sinc.tulip[kinetic3].linkNext = external1;
                f1->sinc.tulip[kinetic3].species = matrix;
                f1->sinc.tulip[kinetic3].blockType = tv3;
                f1->sinc.tulip[kinetic3].name = kinetic;
                f1->sinc.tulip[kinetic3].purpose = ptObject;
                f1->sinc.tulip[kinetic3].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic3].ptRank[0] = 0;
                f1->sinc.tulip[kinetic3].Current[1] = 0;
                f1->sinc.tulip[kinetic3].ptRank[1] = 0;

                
                //active assignment
                    f1->sinc.tulip[Ha].linkNext = kinetic1;
                
                    f1->sinc.tulip[Iterator].linkNext = kinetic1;
                
            }
            
            else if ( bootBodies == four && bootType == electron){
                
                f1->sinc.tulip[interaction12].species = matrix;
                f1->sinc.tulip[interaction12].name = interactionExchange;
                f1->sinc.tulip[interaction12].blockType = e12;
                f1->sinc.tulip[interaction12].NBody = two;
                f1->sinc.tulip[interaction12].Current[0] = CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction12].ptRank[0] = 0;
                f1->sinc.tulip[interaction12].linkNext = interaction13;
                f1->sinc.tulip[interaction12].header = Cube;
                
                f1->sinc.tulip[interaction13].species = matrix;
                f1->sinc.tulip[interaction13].name = interactionExchange;
                f1->sinc.tulip[interaction13].blockType = e13;
                f1->sinc.tulip[interaction13].NBody = two;
                f1->sinc.tulip[interaction13].Current[0] = CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction13].ptRank[0] = 0;
                f1->sinc.tulip[interaction13].linkNext = interaction23;
                f1->sinc.tulip[interaction13].header = Cube;
                
                f1->sinc.tulip[interaction23].species = matrix;
                f1->sinc.tulip[interaction23].name = interactionExchange;
                f1->sinc.tulip[interaction23].blockType = e23;
                f1->sinc.tulip[interaction23].NBody = two;
                f1->sinc.tulip[interaction23].Current[0] = CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction23].ptRank[0] = 0;
                f1->sinc.tulip[interaction23].linkNext = interaction14;
                f1->sinc.tulip[interaction23].header = Cube;
                
                f1->sinc.tulip[interaction14].species = matrix;
                f1->sinc.tulip[interaction14].name = interactionExchange;
                f1->sinc.tulip[interaction14].blockType = e14;
                f1->sinc.tulip[interaction14].NBody = two;
                f1->sinc.tulip[interaction14].Current[0] = CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction14].ptRank[0] = 0;
                f1->sinc.tulip[interaction14].linkNext = interaction24;
                f1->sinc.tulip[interaction14].header = Cube;
                
                f1->sinc.tulip[interaction24].species = matrix;
                f1->sinc.tulip[interaction24].name = interactionExchange;
                f1->sinc.tulip[interaction24].blockType = e24;
                f1->sinc.tulip[interaction24].NBody = two;
                f1->sinc.tulip[interaction24].Current[0] = CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction24].ptRank[0] = 0;
                f1->sinc.tulip[interaction24].linkNext = interaction34;
                f1->sinc.tulip[interaction24].header = Cube;
                
                f1->sinc.tulip[interaction34].species = matrix;
                f1->sinc.tulip[interaction34].name = interactionExchange;
                f1->sinc.tulip[interaction34].blockType = e34;
                f1->sinc.tulip[interaction34].NBody = two;
                f1->sinc.tulip[interaction34].Current[0] = CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction34].ptRank[0] = 0;
                f1->sinc.tulip[interaction34].header = Cube;
                
                f1->sinc.tulip[harmonium1].spinor = none;
                f1->sinc.tulip[harmonium1].linkNext = harmonium2;
                f1->sinc.tulip[harmonium1].NBody = one;
                f1->sinc.tulip[harmonium1].species = matrix;
                f1->sinc.tulip[harmonium1].blockType = tv1;
                f1->sinc.tulip[harmonium1].name = harmonium;
                f1->sinc.tulip[harmonium1].purpose = ptObject;
                f1->sinc.tulip[harmonium1].Current[0] = part(f1, harmonium1);
                f1->sinc.tulip[harmonium1].ptRank[0] = 0;
                
                f1->sinc.tulip[harmonium2].spinor = none;
                f1->sinc.tulip[harmonium2].linkNext = harmonium3;
                f1->sinc.tulip[harmonium2].NBody = one;
                f1->sinc.tulip[harmonium2].species = matrix;
                f1->sinc.tulip[harmonium2].blockType = tv2;
                f1->sinc.tulip[harmonium2].name = harmonium;
                f1->sinc.tulip[harmonium2].purpose = ptObject;
                f1->sinc.tulip[harmonium2].Current[0] = part(f1, harmonium2);
                f1->sinc.tulip[harmonium2].ptRank[0] = 0;
                
                f1->sinc.tulip[harmonium3].spinor = none;
                f1->sinc.tulip[harmonium3].linkNext = harmonium4;
                f1->sinc.tulip[harmonium3].NBody = one;
                f1->sinc.tulip[harmonium3].species = matrix;
                f1->sinc.tulip[harmonium3].blockType = tv3;
                f1->sinc.tulip[harmonium3].name = harmonium;
                f1->sinc.tulip[harmonium3].purpose = ptObject;
                f1->sinc.tulip[harmonium3].Current[0] = part(f1, harmonium2);
                f1->sinc.tulip[harmonium3].ptRank[0] = 0;
                
                f1->sinc.tulip[harmonium4].spinor = none;
                f1->sinc.tulip[harmonium4].linkNext = external1;
                f1->sinc.tulip[harmonium4].NBody = one;
                f1->sinc.tulip[harmonium4].species = matrix;
                f1->sinc.tulip[harmonium4].blockType = tv4;
                f1->sinc.tulip[harmonium4].name = harmonium;
                f1->sinc.tulip[harmonium4].purpose = ptObject;
                f1->sinc.tulip[harmonium4].Current[0] = part(f1, harmonium4);
                f1->sinc.tulip[harmonium4].ptRank[0] = 0;
                
                
                f1->sinc.tulip[Ha].Partition = 0;//
                f1->sinc.tulip[Ha].linkNext = external1   ;
                f1->sinc.tulip[Ha].species = matrix;
                f1->sinc.tulip[Ha].header = bootShape;
                f1->sinc.tulip[Ha].name = Ha;
                
                
                //FOR LINEAR
                f1->sinc.tulip[external1].spinor = none;
                f1->sinc.tulip[external1].linkNext = external2;
                f1->sinc.tulip[external1].NBody = one;
                f1->sinc.tulip[external1].species = matrix;
                f1->sinc.tulip[external1].blockType = tv1;
                f1->sinc.tulip[external1].name = linear;
                f1->sinc.tulip[external1].purpose = ptObject;
                f1->sinc.tulip[external1].Current[0] = part(f1, linear);
                f1->sinc.tulip[external1].ptRank[0] = 0;
                
                f1->sinc.tulip[kinetic1].linkNext = kinetic2;
                f1->sinc.tulip[kinetic1].NBody = one;
                f1->sinc.tulip[kinetic1].species = matrix;
                f1->sinc.tulip[kinetic1].blockType = tv1;
                f1->sinc.tulip[kinetic1].name = kinetic;
                f1->sinc.tulip[kinetic1].purpose = ptObject;
                f1->sinc.tulip[kinetic1].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic1].ptRank[0] = 0;
                f1->sinc.tulip[kinetic1].Current[1] = 0;
                f1->sinc.tulip[kinetic1].ptRank[1] = 0;

                //FOR LINEAR
                f1->sinc.tulip[external2].spinor = none;
                f1->sinc.tulip[external2].linkNext = external3;
                f1->sinc.tulip[external2].NBody = one;
                f1->sinc.tulip[external2].species = matrix;
                f1->sinc.tulip[external2].name = linear;
                f1->sinc.tulip[external2].blockType = tv2;
                f1->sinc.tulip[external2].purpose = ptObject;
                f1->sinc.tulip[external2].Current[0] = part(f1, linear);
                f1->sinc.tulip[external2].ptRank[0] = 0;
                
                f1->sinc.tulip[kinetic2].linkNext = kinetic3;
                f1->sinc.tulip[kinetic2].NBody = one;
                f1->sinc.tulip[kinetic2].species = matrix;
                f1->sinc.tulip[kinetic2].blockType = tv2;
                f1->sinc.tulip[kinetic2].name = kinetic;
                f1->sinc.tulip[kinetic2].purpose = ptObject;
                f1->sinc.tulip[kinetic2].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic2].ptRank[0] = 0;
                f1->sinc.tulip[kinetic2].Current[1] = 0;
                f1->sinc.tulip[kinetic2].ptRank[1] = 0;
                
                f1->sinc.tulip[external3].spinor = none;
                f1->sinc.tulip[external3].linkNext = external4;
                f1->sinc.tulip[external3].NBody = one;
                f1->sinc.tulip[external3].species = matrix;
                f1->sinc.tulip[external3].name = linear;
                f1->sinc.tulip[external3].blockType = tv3;
                f1->sinc.tulip[external3].purpose = ptObject;
                f1->sinc.tulip[external3].Current[0] = part(f1, linear);
                f1->sinc.tulip[external3].ptRank[0] = 0;
                
                f1->sinc.tulip[kinetic3].linkNext = kinetic4;
                f1->sinc.tulip[kinetic3].NBody = one;
                f1->sinc.tulip[kinetic3].species = matrix;
                f1->sinc.tulip[kinetic3].blockType = tv3;
                f1->sinc.tulip[kinetic3].name = kinetic;
                f1->sinc.tulip[kinetic3].purpose = ptObject;
                f1->sinc.tulip[kinetic3].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic3].ptRank[0] = 0;
                f1->sinc.tulip[kinetic3].Current[1] = 0;
                f1->sinc.tulip[kinetic3].ptRank[1] = 0;

                f1->sinc.tulip[external4].spinor = none;
                f1->sinc.tulip[external4].linkNext = interaction12;
                f1->sinc.tulip[external4].NBody = one;
                f1->sinc.tulip[external4].species = matrix;
                f1->sinc.tulip[external4].name = linear;
                f1->sinc.tulip[external4].blockType = tv4;
                f1->sinc.tulip[external4].purpose = ptObject;
                f1->sinc.tulip[external4].Current[0] = part(f1, linear);
                f1->sinc.tulip[external4].ptRank[0] = 0;
                
                if ( c1->i.springFlag )
                    f1->sinc.tulip[kinetic4].linkNext = harmonium1;
                else
                    f1->sinc.tulip[kinetic4].linkNext = external1;
                f1->sinc.tulip[kinetic4].NBody = one;
                f1->sinc.tulip[kinetic4].species = matrix;
                f1->sinc.tulip[kinetic4].blockType = tv4;
                f1->sinc.tulip[kinetic4].name = kinetic;
                f1->sinc.tulip[kinetic4].purpose = ptObject;
                f1->sinc.tulip[kinetic4].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic4].ptRank[0] = 0;
                f1->sinc.tulip[kinetic4].Current[1] = 0;
                f1->sinc.tulip[kinetic4].ptRank[1] = 0;

                //active assignment
                    f1->sinc.tulip[Ha].linkNext = kinetic1;
                
                    f1->sinc.tulip[Iterator].linkNext = external1;
                
            }
            else if (bootType == h2plus && bootBodies == two){
                
                
                
                f1->sinc.tulip[interaction12].species = matrix;
                f1->sinc.tulip[interaction12].name = interactionExchange;
                f1->sinc.tulip[interaction12].linkNext = interaction12Plus;
                f1->sinc.tulip[interaction12].blockType = e12;
                f1->sinc.tulip[interaction12].NBody = two;
                f1->sinc.tulip[interaction12].Current[0] = CanonicalRank(f1, interactionExchange,0);
                f1->sinc.tulip[interaction12].ptRank[0] = 0;
                f1->sinc.tulip[interaction12].header = Cube;
                
                f1->sinc.tulip[interaction12Plus].species = matrix;
                f1->sinc.tulip[interaction12Plus].name = interactionExchangePlus;
                f1->sinc.tulip[interaction12Plus].linkNext = harmonium1;
                f1->sinc.tulip[interaction12Plus].blockType = e12;
                f1->sinc.tulip[interaction12Plus].NBody = two;
                f1->sinc.tulip[interaction12Plus].Current[0] = CanonicalRank(f1, interactionExchangePlus,0);
                f1->sinc.tulip[interaction12Plus].ptRank[0] = 0;
                f1->sinc.tulip[interaction12Plus].header = Cube;
                
                f1->sinc.tulip[harmonium1].spinor = none;
                f1->sinc.tulip[harmonium1].NBody = one;
                f1->sinc.tulip[harmonium1].species = matrix;
                f1->sinc.tulip[harmonium1].blockType = tv1;
                f1->sinc.tulip[harmonium1].name = harmonium;
                f1->sinc.tulip[harmonium1].purpose = ptObject;
                f1->sinc.tulip[harmonium1].Current[0] = part(f1, harmonium1)-1;
                f1->sinc.tulip[harmonium1].ptRank[0] = 1;//nock out one dimension!!!...make it a cylindrical R axis


                f1->sinc.tulip[Ha].Partition = 0;//
                f1->sinc.tulip[Ha].species = matrix;
                f1->sinc.tulip[Ha].TBody = h2plus;
                f1->sinc.tulip[Ha].header = bootShape;
                f1->sinc.tulip[Ha].name = Ha;
                
                
                //FOR LINEAR
                f1->sinc.tulip[external1].spinor = none;
                f1->sinc.tulip[external1].linkNext = interaction12;
                f1->sinc.tulip[external1].NBody = one;
                f1->sinc.tulip[external1].TBody = proton;
                f1->sinc.tulip[external1].species = matrix;
                f1->sinc.tulip[external1].blockType = tv1;
                f1->sinc.tulip[external1].name = linear;
                f1->sinc.tulip[external1].purpose = ptObject;
                f1->sinc.tulip[external1].Current[0] = part(f1, linear);
                f1->sinc.tulip[external1].ptRank[0] = 0;
                
                f1->sinc.tulip[kinetic1].linkNext = kinetic2;
                f1->sinc.tulip[kinetic1].NBody = one;
                f1->sinc.tulip[kinetic1].TBody = proton;
                f1->sinc.tulip[kinetic1].species = matrix;
                f1->sinc.tulip[kinetic1].blockType = tv1;
                f1->sinc.tulip[kinetic1].name = kineticMass;
                f1->sinc.tulip[kinetic1].purpose = ptObject;
                f1->sinc.tulip[kinetic1].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic1].ptRank[0] = 0;
                f1->sinc.tulip[kinetic1].Current[1] = 0;
                f1->sinc.tulip[kinetic1].ptRank[1] = 0;
                
                //FOR LINEAR
                
//                if ( c1->i.springFlag )
//                    f1->sinc.tulip[kinetic2].linkNext = harmonium1;
//                else
                
                f1->sinc.tulip[kinetic2].linkNext = external1;
                f1->sinc.tulip[kinetic2].NBody = one;
                f1->sinc.tulip[kinetic2].species = matrix;
                f1->sinc.tulip[kinetic2].blockType = tv2;
                f1->sinc.tulip[kinetic2].name = kinetic;
                f1->sinc.tulip[kinetic2].purpose = ptObject;
                f1->sinc.tulip[kinetic2].Current[0] = part(f1, kinetic);
                f1->sinc.tulip[kinetic2].ptRank[0] = 0;
                f1->sinc.tulip[kinetic2].Current[1] = 0;
                f1->sinc.tulip[kinetic2].ptRank[1] = 0;
                
                
                //active assignment
                f1->sinc.tulip[Ha].linkNext = kinetic1;
                
                f1->sinc.tulip[Iterator].linkNext = kinetic1;
                
            }
        }
        
    }

    printf("boot complete\n");
    fflush(stdout);
//    if ( 0)
//        tInnerTest(f1, kinetic, copy);
    return 0;
}
