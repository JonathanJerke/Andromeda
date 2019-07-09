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
        i->i.Na = 1;
        i->i.atoms[1].label.Z =1;
        i->i.atoms[1].position[1] = 0;
        i->i.atoms[1].position[2] = 0;
        i->i.atoms[1].position[3] = 0;
        
    } else
        if ( number== 2 ){
            i->i.Na = 2;

            i->i.atoms[1].label.Z = 1;
            i->i.atoms[1].position[1] = scale;
            i->i.atoms[1].position[2] = 0;
            i->i.atoms[1].position[3] = 0;
            i->i.atoms[2].label.Z = 1;
            i->i.atoms[2].position[1] = -scale;
            i->i.atoms[2].position[2] = 0;
            i->i.atoms[2].position[3] = 0;
            
        }
        else
            if ( number == 3 ){
                i->i.Na = 3;

                i->i.atoms[1].label.Z = 1;
                i->i.atoms[1].position[1] = 2*scale;
                i->i.atoms[1].position[2] = -scale;
                i->i.atoms[1].position[3] = -scale;
                i->i.atoms[2].label.Z = 1;
                i->i.atoms[2].position[1] = -scale;
                i->i.atoms[2].position[2] = -scale;
                i->i.atoms[2].position[3] = 2*scale;
                i->i.atoms[3].label.Z = 1;
                i->i.atoms[3].position[1] = -scale;
                i->i.atoms[3].position[2] = 2*scale;
                i->i.atoms[3].position[3] = -scale;
                
            }
            else
                if ( number == 4 ){
                    
                    i->i.Na = 4;

                    i->i.atoms[1].label.Z = 1;
                    i->i.atoms[1].position[1] = scale;
                    i->i.atoms[1].position[2] = 0;
                    i->i.atoms[1].position[3] = -scale/sqrt(2.);
                    i->i.atoms[2].label.Z = 1;
                    i->i.atoms[2].position[1] = -scale;
                    i->i.atoms[2].position[2] = 0;
                    i->i.atoms[2].position[3] = -scale/sqrt(2.);
                    i->i.atoms[3].label.Z = 1;
                    i->i.atoms[3].position[1] = 0;
                    i->i.atoms[3].position[2] = scale;
                    i->i.atoms[3].position[3] = scale/sqrt(2.);
                    i->i.atoms[4].label.Z = 1;
                    i->i.atoms[4].position[1] = 0;
                    i->i.atoms[4].position[2] = -scale;
                    i->i.atoms[4].position[3] = scale/sqrt(2.);
                    
                }
                else
                    if ( number == 5 ){
                        
                        i->i.Na = 1;
                        
                        i->i.atoms[1].label.Z = 2;
                        i->i.atoms[1].position[1] = 0;
                        i->i.atoms[1].position[2] = 0;
                        i->i.atoms[1].position[3] = 0;
                        
                    }

    }
struct field initField (void ) {
    INT_TYPE space;
    struct field i;
    i.i.attack = 0.;
    i.i.body = nada;
    
    i.i.D = 1.0;
    i.i.around = 0;

    i.i.d = 1.0;
    i.i.epi = 0;

    i.i.nStates = 0;
    i.i.qFloor = 0;
    i.i.nOperator = 0;

    
    i.i.bRank = 0;
    i.i.iRank = 0;
    i.i.filter = 0;
    i.i.irrep = 0;
    i.i.collect = 0;
    
    i.i.Iterations = 0;
    i.f.boot = noMatrices;
    for( space = 0 ; space <= SPACE ; space++){
        i.f.rose[space].stream = NULL;
        i.f.rose[space].basisList = NULL;
    }
    i.f.bootedMemory = 0;
    i.f.tulip = NULL;
    i.i.files = 0;
    i.i.filesVectorOperator = 0;
    
#ifdef APPLE
    i.i.d = 2.;

    i.i.bRank = 4;
    i.i.iRank = 2;
    i.i.nStates = 1;
    i.i.qFloor = 27;
    i.i.filter = 1;
    i.i.irrep = 1;
    i.f.boot = fullMatrices;
    i.i.body = one;
    i.i.epi = 4;
    i.i.qFloor = 46;
    i.i.body = three;
    
#endif
    return i;
}
struct calculation initCal (void ) {
    struct calculation i;
    

#ifdef APPLE
   // i.i.OCSBflag = 0;
    i.i.springConstant = 0.;
    i.i.springFlag = 0;
    i.i.RAMmax = 1;
    i.rt.runFlag = 0;
    i.i.vectorMomentum = 0.;
    i.i.decomposeRankMatrix = 1;
    i.i.orgClamp = 2.;
    i.i.Angstroms = 0;
    i.rt.TARGET = 1e-6;
    i.rt.targetCondition = 1e-7;
    i.rt.ALPHA = 1e-8;
    i.rt.CANON = 1e-7;
    i.rt.vCANON = 1e-6;
    i.rt.TOL = 1e5;
    i.rt.maxEntropy = 1;
    i.i.complexType =2;
    i.i.level = 0.1;
    
    
    i.rt.powDecompose = 2;
    
        i.i.turn = 1.;
        i.i.param1 = 1.;
        i.i.param2 = 1.;
        i.i.interval = 1;
        i.i.scalar = 1.;
    
        i.i.magFlag = 0;
        i.i.mag = 0.1;
        //THESE
        if ( SPACE == 3 ){
            i.rt.calcType = electronicStuctureCalculation;
            i.rt.runFlag = 0;
        } else {

            i.rt.calcType = clampProtonElectronCalculation;
            i.i.minClamp = 0.5;
            i.i.maxClamp = 4;
            
            i.i.orgClamp = 0.5*(i.i.minClamp+i.i.maxClamp);


        }
    i.i.massElectron = 1.;
    i.i.massProton = 1836.15267245;
    i.i.massClampPair = 1836.15267245;
    resetExternal(&i, 1, 1);

    i.rt.phaseType = buildFoundation ;
    
    i.i.twoBody.func.fn = nullFunction;
    i.i.oneBody.func.fn = Coulomb;
    //i.i.springFlag = 1;
    i.i.springConstant = 0.25;
    i.i.canonRank = 50;
    i.i.twoBody.num = 0;
    i.i.twoBody.func.interval  = 0;
    i.i.twoBody.func.param[0]  = 1;
    i.i.twoBody.func.param[1]  = 1;
    i.i.twoBody.func.param[2]  = 1;

    i.i.oneBody.num = 50;
    i.i.oneBody.func.interval  = 0;
    i.i.oneBody.func.param[0]  = 1.;
    i.i.oneBody.func.param[1]  = 1;
    i.i.oneBody.func.param[2]  = 1;
#else
    i.i.springConstant = 0.;
    i.i.springFlag = 0;
    i.i.RAMmax = 1;
    i.rt.runFlag = 0;
    i.i.vectorMomentum = 0.;
    i.i.decomposeRankMatrix = 1;
    i.i.orgClamp = 2.;
    i.i.Angstroms = 0;
    i.rt.TARGET = 1e-6;
    i.rt.targetCondition = 1e-7;
    i.rt.ALPHA = 1e-8;
    i.rt.CANON = 1e-7;
    i.rt.vCANON = 1e-6;
    i.rt.TOL = 1e5;
    i.rt.maxEntropy = 1;
    i.i.level = 10000;
    i.i.complexType =2;
    i.i.level = 1;
    i.rt.powDecompose = 2;
    i.i.turn = 1.;
    i.i.param1 = 1.;
    i.i.param2 = 1.;
    i.i.interval = 1;
    i.i.scalar = 1.;
    i.i.magFlag = 0;
    i.i.mag = 0.1;
    //THESE
    i.rt.calcType = nullCalculation;
    i.rt.runFlag = 0;
    i.i.massElectron = 1.;
    i.i.massProton = 1836.15267245;
    i.i.massClampPair = 1836.15267245;
    
    i.rt.phaseType = buildFoundation ;
    i.i.twoBody.func.fn = nullFunction;
    i.i.oneBody.func.fn = nullFunction;
    i.i.canonRank = 50;
    i.i.twoBody.num = 0;
    i.i.twoBody.func.interval  = 0;
    i.i.twoBody.func.param[0]  = 1;
    i.i.twoBody.func.param[1]  = 1;
    i.i.twoBody.func.param[2]  = 1;
    
    i.i.oneBody.num = 50;
    i.i.oneBody.func.interval  = 0;
    i.i.oneBody.func.param[0]  = 1.;
    i.i.oneBody.func.param[1]  = 1;
    i.i.oneBody.func.param[2]  = 1;

#endif
    return i;
}

INT_TYPE fModel ( struct sinc_label * f1){
    INT_TYPE i;
    
    if ( f1->bootedMemory ){
        for ( i = 0;i <= SPACE ; i++){
            free(f1->rose[i].stream);
            f1->rose[i].stream = NULL;

            if ( i < SPACE ){
                free(f1->rose[i].basisList);
                f1->rose[i].basisList = NULL;
            }
        }
        {
            free(f1->tulip);
            f1->tulip = NULL;
        }
    }
    return 0;
}

//INT_TYPE multModel( struct calculation * c1, INT_TYPE nv,struct field *v, struct field *f){
//    struct sinc_label *f1 = &f->f;
//    INT_TYPE len, bl,i,periodic = 1,space, N1 =0,N2 =0;
//    
//    for ( i = 0; i < nv ; i++){
//        N1 += v[i].i.epi*2+1;
//        N2 += v[i].i.around*2+1;
//    }
//    
//    //define vectors
//    f1->rose[SPACE].component = nullComponent;
//    f1->rose[SPACE].body = nada;
//    f1->rose[SPACE].particle = nullParticle;
//    if ( c1->rt.calcType == electronicStuctureCalculation){
//        periodic = 1;
//        for ( space = 0; space < SPACE ; space++){
//            f1->rose[space].particle = electron;
//            f1->rose[space].component = spatialComponent1+space%COMPONENT + ((c1->rt.runFlag/periodic)%2)*3;
//            periodic *= 2;
//            
//            len = N1;
//            
//    
//            f1->rose[space].count1Basis = len;
//            if (f1->rose[space].component > 3  )
//                f1->rose[space].count1Basis *= 2;
//            f1->rose[space].origin = 0.;
//            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*N1);
//            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
//                f1->rose[space].basisList[bl] = defineBasis(nullNote, spatialComponent1+space%COMPONENT + ((c1->rt.runFlag/periodic)%2)*3, SincBasisElement, f->i.d, 0.,             f1->rose[space].count1Basis);
//            
//        }
//    } else  if ( c1->rt.calcType == clampProtonElectronCalculation){
//        periodic = 1;
//        for ( space = 0; space < COMPONENT ; space++){
//            f1->rose[space].particle = electron;
//            f1->rose[space].component = spatialComponent1+space+ ((c1->rt.runFlag/periodic)%2)*3;
//            periodic *= 2;
//            
//            f1->rose[space].basis = nullBasisElement;
//            f1->rose[space].body = nada;
//            f1->rose[space].lattice = 0.;
//            
//            f1->rose[space].count1Basis = N1;
//            f1->rose[space].origin = 0.;
//            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
//            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
//                f1->rose[space].basisList[bl] = grabBasis(f1, space % COMPONENT, space / COMPONENT, bl);
//
//        }
//        periodic = 1;
//        for ( space = COMPONENT; space <= COMPONENT ; space++){
//            f1->rose[space].particle = proton;
//            f1->rose[space].component = spatialComponent1+space%COMPONENT+ ((c1->rt.runFlag/periodic)%2)*3;
//            periodic *= 2;
//            
//            f1->rose[space].basis = nullBasisElement;
//            f1->rose[space].body = nada;
//            f1->rose[space].lattice = 0.;
//            
//            f1->rose[space].count1Basis = N2;
//            f1->rose[space].origin = c1->i.orgClamp;
//            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
//            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
//                f1->rose[space].basisList[bl] = grabBasis(f1, space % COMPONENT, space / COMPONENT, bl);
//
//        }
//        for ( space = COMPONENT+1; space < SPACE ; space++){
//            f1->rose[space].particle = proton;
//            f1->rose[space].component = spatialComponent1+space%COMPONENT+ ((c1->rt.runFlag/periodic)%2)*3;
//            periodic *= 2;
//            
//            
//            f1->rose[space].lattice = 0;
//            f1->rose[space].basis = nullBasisElement;
//            f1->rose[space].body = nada;
//
//            f1->rose[space].count1Basis = 0 ;
//            f1->rose[space].origin = 0.;
//            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
//            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
//                f1->rose[space].basisList[bl] = grabBasis(f1, space % COMPONENT, space / COMPONENT, bl);
//
//        }
//        
//    }
//    return N1;
//}

INT_TYPE singleSincModel( struct calculation * c1, struct field * f){
    struct sinc_label *f1 = &f->f;
    
    INT_TYPE bl,periodic = 1,space, N1 = f->i.epi*2+1,N2=2*f->i.around+1;
    //define vectors
    f1->rose[SPACE].component = nullComponent;
    f1->rose[SPACE].body = nada;
    f1->rose[SPACE].particle = nullParticle;
    if ( c1->rt.calcType == electronicStuctureCalculation){
        periodic = 1;
        for ( space = 0; space < SPACE ; space++){
            f1->rose[space].particle = electron;
            f1->rose[space].basis = SincBasisElement;
            f1->rose[space].body = f->i.body;
            f1->rose[space].component = spatialComponent1+space%COMPONENT + ((c1->rt.runFlag/periodic)%2)*3;
            if (  f1->rose[space].component > 3){
                printf("\nperiodic boundaries %d\n",f1->rose[space].component);
            }
            periodic *= 2;
            f1->rose[space].count1Basis = N1;
            if (f1->rose[space].component > 3  )
                f1->rose[space].count1Basis *= 2;
            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
                f1->rose[space].basisList[bl] = defineSincBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, f->i.d, 0.,             f1->rose[space].count1Basis,bl);

        }
    } else  if ( c1->rt.calcType == clampProtonElectronCalculation){
        periodic = 1;
        for ( space = 0; space < COMPONENT ; space++){
            f1->rose[space].particle = electron;
            f1->rose[space].basis = SincBasisElement;
            f1->rose[space].body = f->i.body;
            f1->rose[space].component = spatialComponent1+space+ ((c1->rt.runFlag/periodic)%2)*3;
            periodic *= 2;
            if (  f1->rose[space].component > 3){
                printf("\nperiodic boundaries %d\n",f1->rose[space].component);
            }

            f1->rose[space].count1Basis = N1;
            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
                f1->rose[space].basisList[bl] = defineSincBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, f->i.d, 0.,             f1->rose[space].count1Basis,bl);


        }
        periodic = 1;
        for ( space = COMPONENT; space <= COMPONENT ; space++){
            f1->rose[space].particle = proton;
            f1->rose[space].basis = SincBasisElement;
            f1->rose[space].body = one;
            f1->rose[space].component = spatialComponent1+space%COMPONENT+ ((c1->rt.runFlag/periodic)%2)*3;
            periodic *= 2;
            if (  f1->rose[space].component > 3){
                printf("\nperiodic boundaries %d\n",f1->rose[space].component);
            }

            f1->rose[space].count1Basis = N2;
            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
                f1->rose[space].basisList[bl] = defineSincBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, f->i.D, 0.,             f1->rose[space].count1Basis,bl);


        }
        for ( space = COMPONENT+1; space < SPACE ; space++){
            f1->rose[space].particle = proton;
            f1->rose[space].basis = SincBasisElement;
            f1->rose[space].body = nada;
            f1->rose[space].component = spatialComponent1+space%COMPONENT+ ((c1->rt.runFlag/periodic)%2)*3;
            periodic *= 2;
            f1->rose[space].count1Basis = 0;
            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*1);
            
        }
        
    }
    return N1;
}

INT_TYPE singleGaussModel( struct calculation * c1, struct field * f){
    struct sinc_label *f1 = &f->f;
    double ble[] = {0.2,1.,5.};
    INT_TYPE bl,space;
    //define vectors
    f1->rose[SPACE].component = nullComponent;
    f1->rose[SPACE].body = nada;
    f1->rose[SPACE].particle = nullParticle;
    if ( c1->rt.calcType == electronicStuctureCalculation){
        for ( space = 0; space < SPACE ; space++){
            f1->rose[space].particle = electron;
            f1->rose[space].basis = GaussianBasisElement;
            f1->rose[space].body = 1;
            f1->rose[space].component = spatialComponent1+space%COMPONENT ;
            f1->rose[space].count1Basis = 3;//change
            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
                f1->rose[space].basisList[bl] = defineGaussBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, ble[bl], 0.,             f1->rose[space].count1Basis,0);
            
        }
    }
    return 1;
}

INT_TYPE iModel( struct calculation * c1, struct field *f){
    struct name_label l2;
    if ( OVERFLAG )
        singleGaussModel(c1, f);
    else
        singleSincModel(c1, f);

    struct sinc_label *f1 = &f->f;

    if(1){
        printf("\n\n*Parameters\n");
        printf("\t target\t\t\t%1.1f\n", -log(c1->rt.TARGET)/log(10));//quality of decomposition
        printf("\t tolerance\t\t%1.1f\n", log(c1->rt.TOL)/log(10));//max condition of foundation    s
        printf("\t threshold\t\t%1.1f\n", -log(c1->rt.CANON)/log(10));//matrix training standard
        printf("\t vectorThreshold\t%1.1f\n", -log(c1->rt.vCANON)/log(10));//vector training standard
        printf("\t condition\t\t%1.1f\n", -log(c1->rt.ALPHA )/log(10));//Beylkin parameter
        printf(".Parameters\n\n");
        
        printf("\n\n\t  \tBox %d \t: Lattice  %1.3f \t\n\t\n",2* f->i.epi + 1,f->i.d);
        if ( SPACE > 3 )
            printf("\n\n\t  \tBox %d \t: Latte  %1.3f \t\n\t\n",2* f->i.around + 1,f->i.D );

    }
    struct runTime * rt = &c1->rt;
    
        INT_TYPE space;
        {
            f->f.rt = rt;
        }
        
        enum bodyType bootBodies = f->i.body;
        INT_TYPE ra = tPerms(bootBodies);//tSize(bootBodies);
        INT_TYPE N12 = f->i.epi;
        INT_TYPE N1 =  2*N12+1;
        enum shape bootShape;
        INT_TYPE maxVector = imax(0*c1->i.decomposeRankMatrix, imax(f->i.bRank,imax(1+f->i.iRank,f->i.Iterations+f->i.iRank)));
        //rds defined in input.c
        
        bootShape = Cube;
        

    INT_TYPE FloorV = imax(0, f->i.qFloor), CeilV = imax(0,0);
        INT_TYPE maxArray,EV,maxEV,NV = 0,FV = FloorV+CeilV ;
    
        EV = FV;
    maxEV =EV*(imax(f->i.Iterations,1));
    maxArray = imax(f->i.nStates,maxEV);//slip Nb into spectra...
        
        f1->maxEV = maxArray;
        enum division vectorOperator  = eigenVectors + 1 +  f->i.nStates+maxEV;
        f1->vectorOperator = vectorOperator;
        enum division end  = vectorOperator+f->i.nOperator;
        f1->end = end;
        f1->tulip = malloc ( (end+1) * sizeof(struct name_label));


        INT_TYPE outVector  = imax(c1->i.oneBody.num , c1->i.twoBody.num )*maxVector;
        
        {//defaults
            //define vectors end
            {
                enum division label1;
                for ( label1 = 0 ;label1 <= end; label1++){
                    f1->tulip[label1].name = label1;
                    f1->tulip[label1].Partition = 0;
                    f1->tulip[label1].header = Cube;
                    f1->tulip[label1].spinor = real;

                    f1->tulip[label1].species = scalar;
                    f1->tulip[label1].linkNext = nullName;
                    f1->tulip[label1].memory = objectAllocation;
                    for ( space = 0; space <= SPACE ; space++)
                        f1->tulip[label1].space[space].Address = -1;
                    f1->tulip[label1].space[SPACE].block = id0;
                    for ( space = 0; space <= SPACE ; space++){
                        {
                            f1->tulip[label1].space[space].block = id0;//matrix prototype
                            f1->tulip[label1].space[space].body = nada;//matrix prototype
                        }
                    }
                    tClear(*f1,label1);
                }
                
            }
        }//defaults
        fromBeginning(*f1, kinetic, 0);
        f1->tulip[kinetic].Partition = COMPONENT;//
//        if ( c1->rt.runFlag > 0 )
//            f1->tulip[kinetic].spinor = cmpl;
        assignOneWithPointers(*f1, kinetic,all);

        fromBeginning(*f1, kineticMass, kinetic);
        f1->tulip[kineticMass].Partition = COMPONENT;//
        assignOneWithPointers(*f1, kineticMass,all);
        struct name_label u = f1->tulip[kineticMass];

        fromBeginning(*f1, protonRepulsion, kineticMass);
        f1->tulip[protonRepulsion].Partition = c1->i.oneBody.num;//
        assignOneWithPointers(*f1, protonRepulsion,all);
        
        fromBeginning(*f1, vectorMomentum, protonRepulsion);
        assignOneWithPointers(*f1, vectorMomentum,electron);
        f1->tulip[vectorMomentum].spinor = cmpl;
        f1->tulip[vectorMomentum].Partition = COMPONENT * c1->i.springFlag+4*c1->i.magFlag;//

        fromBeginning(*f1, harmonium, vectorMomentum);
        assignOneWithPointers(*f1, harmonium,electron);
        f1->tulip[harmonium].Partition = 0;//
        
        fromBeginning(*f1, X, harmonium);
        assignOneWithPointers (*f1, X,electron);
        f1->tulip[X].Partition =  0 ;//make it a semi-local Gaussian * x
        
        fromBeginning(*f1, linear, X);
        assignOneWithPointers (*f1, linear,electron);
        f1->tulip[linear].Partition = c1->i.Na*c1->i.oneBody.num+c1->i.twoBody.num ;//

        fromBeginning(*f1, overlap, linear);
        assignOneWithPointers (*f1, overlap,all);
    for ( space = 0 ; space < SPACE ; space++)
        f1->tulip[overlap].Partition = 1;//

        {
            fromBeginning(*f1, build, overlap);
//            f1->tulip[build].Partition = f->i.OCSBflag*!(!f->i.sectors)*bootBodies*(6+(/*HERE 3*/part(*f1, kinetic)+part(*f1, vectorMomentum)) + matrixNumber);//easily reduce in cheaper ways!
//            if ( bootBodies == two )
//                f1->tulip[build].Partition += f->i.OCSBflag*f->i.sectors;
//            else if ( bootBodies == three )
//                f1->tulip[build].Partition += f->i.OCSBflag*3*f->i.sectors;
//            else if ( bootBodies == four )
//                f1->tulip[build].Partition += f->i.OCSBflag*6*f->i.sectors;
            
//            f1->tulip[build].Partition = f->i.OCSBflag*f->i.sectors*imax(nG*f->i.sectors,part(*f1, build));//BUILD
            
            //f1->tulip[build].Partition = 0;
            // B2 :  2*(1+S)*3 + 1*sector +  Na*matrixNumber*2
            // B3 :  3*(1+S)*3 + 3*sector  + Na*matrixNumber*3
            f1->tulip[build].species = matrix;
            assignParticle(*f1, build, all, bootBodies);
       //     printf("sectors %d\n", f->i.sectors);
            for ( space = 0; space < SPACE ; space++)
            {
                if ( space == 0 )
                    fromBeginning(*f1, bill1+space, build);
                else
                    fromBeginning(*f1, bill1+space, bill1+space-1);
                
                f1->tulip[bill1+space].Partition = (c1->rt.phaseType == buildFoundation) *vectorLen(*f1, space)*vectorLen(*f1, space) ;
                f1->tulip[bill1+space].memory = bufferAllocation;
               // if ( c1->rt.runFlag )
               //     f1->tulip[bill1+space].spinor = cmpl;
                
            }
            fromBeginning(*f1, eigen, bill1+SPACE-1);
           // f1->tulip[eigen].Partition = f->i.OCSBflag * !(!f->i.sectors)*c1->i.decomposeRankMatrix;
            f1->tulip[eigen].species = matrix;
            assignParticle(*f1, eigen, all, bootBodies);
            for ( space = 0; space < SPACE ; space++)
                if ( f1->rose[space].body != nada)
                    f1->tulip[eigen].space[space].block = tv1;
            
            
            {
                INT_TYPE di,cmpl;
                enum division last = eigen;
                INT_TYPE booting[7];
                for ( di = 0; di < 7 ; di++)
                    booting[di] = 0;
                
                
                if ( f->i.filesVectorOperator )
                    
                {
                    enum bodyType bd;
                    INT_TYPE fi,lines = 0;
                    size_t ms = MAXSTRING;
                    char line0[MAXSTRING];
                    char name[MAXSTRING];
                    char *line = line0;
                    INT_TYPE FIT ;
                    FIT = f->i.filesVectorOperator ;
                    for ( fi =0 ; fi < FIT; fi++){
                        strcpy(name ,f->i.fileVectorOperator[fi]);
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
                                tFromReadToFilename(NULL, line,  name, spins(*f1,eigenVectors)-1,cmpl);
                              //  printf("%s\n", name);
                                fromBeginning(*f1, vectorOperator+lines, last);
                                f1->tulip[vectorOperator+lines].Partition = inputFormat(*f1, name, nullName, 2);
                                f1->tulip[vectorOperator+lines].species = outerVector;
                                for (space = 0; space < SPACE ; space++){
                                    bd = inputFormat(*f1, name, nullName, 100+space/COMPONENT);
                                    f1->tulip[vectorOperator+lines].space[space].body = bd;
                                    printf("%d %d %d\n", lines+1, space, bd);//space/COMPONENT = particle
                                    booting[ bd - bootBodies ] = 1;

                                }
                                
                                last = vectorOperator+lines;
                                lines++;
                            }
                            getline(&line, &ms, fp);
                        }
                        if ( fi > MAXFILE)
                        {
                            printf("too many files, increase MAXSTRING\n");
                            exit(0);
                        }
                        fclose(fp);
                    }
                }
                
                
//                for ( di = 1; di < 4 ; di++)
//                {
//                    fromBeginning(*f1, complement+di, last);
//                    last = complement + di;
//                    f1->tulip[complement+ di].Partition = booting[di] ;//
//                    f1->tulip[complement+ di].spinor = parallel;
//                    f1->tulip[complement+ di].species = outerVector;
//                    assignParticle(*f1, complement+di, electron, di);
//                    
//                    fromBeginning(*f1,complementTwo+di , last);
//                    last = complementTwo + di;
//                    f1->tulip[complementTwo+ di].Partition = booting[di] ;//
//                    f1->tulip[complementTwo+ di].spinor = parallel;
//                    f1->tulip[complementTwo+ di].species = outerVector;
//                    assignParticle(*f1, complementTwo+di, electron, di);
//
//                }
                
                fromBeginning(*f1, eigenVectors, last);
            }
            
            struct name_label u ;
            struct name_label u2;


            if(1 /* !si*/){
                enum division last = eigenVectors;
                INT_TYPE di,d0=0;
#if VERBOSE
                printf("std USERs\n\n");
#endif
                for ( di = 0 ; di < maxEV+f->i.nStates; di++){
                    if ( di ){
                        fromBeginning(*f1, eigenVectors+di, last);
                        last = eigenVectors+di;
                    }
                    f1->tulip[eigenVectors+di].spinor =c1->i.complexType;
                    if ( d0 < f->i.nStates ){
                        f1->tulip[eigenVectors+di].Partition = f->i.bRank;
                        d0++;
                    }
                    else if ( eigenVectors+di <= eigenVectors+f->i.nStates+maxEV){
                        {
                            f1->tulip[eigenVectors+di].Partition = ((di-d0)/EV)+f->i.iRank;
                            NV += spins(*f1,eigenVectors+di )*f1->tulip[eigenVectors+di].Partition;
                        }
                        
                    }else{
                        exit(1);
                    }
                    f1->tulip[eigenVectors+di].species = vector;
                }
                fromBeginning(*f1, diagonalVectorA, last);
                ;
                 u = f1->tulip[eigenVectors];
                 u2 = f1->tulip[eigenVectors+di-1];
                f1->user = eigenVectors + d0;
            }
                
            

            f1->tulip[diagonalVectorA].spinor = parallel;
            f1->tulip[diagonalVectorA].Partition =(c1->rt.phaseType == buildFoundation);
            f1->tulip[diagonalVectorA].species = vector;
            
            fromBeginning(*f1, diagonalVectorB, diagonalVectorA);
            f1->tulip[diagonalVectorB].spinor = parallel;
            f1->tulip[diagonalVectorB].Partition = (c1->rt.phaseType == buildFoundation);
            f1->tulip[diagonalVectorB].species = vector;

        }
        
        INT_TYPE len[SPACE],mxlen=0;
        length(*f1, eigenVectors, len);
        for ( space = 0 ; space < SPACE ; space++)
            if(len[space] > mxlen)
                mxlen = len[space];
        
        INT_TYPE mx1len=0;
        length1(*f1, len);
        for ( space = 0 ; space < SPACE ; space++)
            if(len[space] > mx1len)
                mx1len = len[space];

        
        fromBeginning(*f1,edgeElectronMatrix , diagonalVectorB);
        f1->tulip[edgeElectronMatrix].Partition = 1 ;//
        f1->tulip[edgeElectronMatrix].species = matrix;
        assignOneWithPointers(*f1, edgeElectronMatrix, electron);
        
        fromBeginning(*f1,edgeProtonMatrix,edgeElectronMatrix);
        f1->tulip[edgeProtonMatrix].Partition = 1 ;//
        f1->tulip[edgeProtonMatrix].species = matrix;
        assignOneWithPointers(*f1, edgeProtonMatrix, proton);

        fromBeginning(*f1,productVector,edgeProtonMatrix);
        f1->tulip[productVector].Partition =  1;
        f1->tulip[productVector].species = vector;
        f1->tulip[productVector].spinor = parallel;

        fromBeginning(*f1,permutationVector,productVector);
        f1->tulip[permutationVector].Partition =  f->i.filter *ra*maxVector;
        f1->tulip[permutationVector].species = vector;
        //f1->tulip[permutationVector].spinor = parallel;
        
        fromBeginning(*f1,permutation2Vector,permutationVector);
        f1->tulip[permutation2Vector].Partition = f->i.filter *ra*maxVector;
        f1->tulip[permutation2Vector].species = vector;
        f1->tulip[permutation2Vector].spinor = parallel;
    
        fromBeginning(*f1,canonicalmvVector,permutation2Vector);
        f1->tulip[canonicalmvVector].Partition = 1;
        f1->tulip[canonicalmvVector].species = vector;
        f1->tulip[canonicalmvVector].spinor = parallel;

        fromBeginning(*f1,canonicalmv2Vector,canonicalmvVector);
        f1->tulip[canonicalmv2Vector].Partition = 1;
        f1->tulip[canonicalmv2Vector].species = vector;
        f1->tulip[canonicalmv2Vector].spinor = parallel;

        fromBeginning(*f1,canonicalmv3Vector,canonicalmv2Vector);
        f1->tulip[canonicalmv3Vector].Partition = 1;
        f1->tulip[canonicalmv3Vector].species = vector;
        f1->tulip[canonicalmv3Vector].spinor = parallel;

        fromBeginning(*f1,canonicaldotVector,canonicalmv3Vector);
        f1->tulip[canonicaldotVector].Partition = 1;
        f1->tulip[canonicaldotVector].species = vector;;
        f1->tulip[canonicaldotVector].spinor = parallel;
        
        fromBeginning(*f1,canonicaldot2Vector,canonicaldotVector);
        f1->tulip[canonicaldot2Vector].Partition = 1;
        f1->tulip[canonicaldot2Vector].species = vector;;
        f1->tulip[canonicaldot2Vector].spinor = parallel;
        
        fromBeginning(*f1,canonicaldot3Vector,canonicaldot2Vector);
        f1->tulip[canonicaldot3Vector].Partition = 1;
        f1->tulip[canonicaldot3Vector].species = vector;;
        f1->tulip[canonicaldot3Vector].spinor = parallel;
        
        fromBeginning(*f1,canonicalvvVector,canonicaldot3Vector);
        f1->tulip[canonicalvvVector].Partition = 0;
        f1->tulip[canonicalvvVector].species = vector;;
        f1->tulip[canonicalvvVector].spinor = parallel;
        
        fromBeginning(*f1,canonicalvv2Vector,canonicalvvVector);
        f1->tulip[canonicalvv2Vector].Partition = 0;
        f1->tulip[canonicalvv2Vector].species = vector;;
        f1->tulip[canonicalvv2Vector].spinor = parallel;
        
        fromBeginning(*f1,canonicalvv3Vector,canonicalvv2Vector);
        f1->tulip[canonicalvv3Vector].Partition = 0;
        f1->tulip[canonicalvv3Vector].species = vector;;
        f1->tulip[canonicalvv3Vector].spinor = parallel;

        fromBeginning(*f1,canonicalmeVector,canonicalvv3Vector);
        f1->tulip[canonicalmeVector].Partition = 1;
        f1->tulip[canonicalmeVector].species = vector;;
        f1->tulip[canonicalmeVector].spinor = parallel;
        
        fromBeginning(*f1,canonicalme2Vector,canonicalmeVector);
        f1->tulip[canonicalme2Vector].Partition = 0;
        f1->tulip[canonicalme2Vector].species = vector;;
        f1->tulip[canonicalme2Vector].spinor = parallel;
        
        fromBeginning(*f1,canonicalme3Vector,canonicalme2Vector);
        f1->tulip[canonicalme3Vector].Partition = 0;
        f1->tulip[canonicalme3Vector].species = vector;;
        f1->tulip[canonicalme3Vector].spinor = parallel;

        fromBeginning(*f1,copyVector,canonicalme3Vector);
        f1->tulip[copyVector].Partition = maxVector;
        f1->tulip[copyVector].species = vector;
        f1->tulip[copyVector].spinor = parallel;

        fromBeginning(*f1,copyTwoVector,copyVector);
        f1->tulip[copyTwoVector].Partition =   0*maxVector;//f->i.bRank*f->i.bRank;
        f1->tulip[copyTwoVector].species = vector;
        f1->tulip[copyTwoVector].spinor = parallel;

        fromBeginning(*f1,copyThreeVector,copyTwoVector);
        f1->tulip[copyThreeVector].Partition =   0*maxVector;//f->i.bRank*f->i.bRank;
        f1->tulip[copyThreeVector].species = vector;
        f1->tulip[copyThreeVector].spinor = parallel;

        fromBeginning(*f1,copyFourVector,copyThreeVector);
        f1->tulip[copyFourVector].Partition =    0*maxVector;//f->i.bRank*f->i.bRank;
        f1->tulip[copyFourVector].species = vector;
        f1->tulip[copyFourVector].spinor = parallel;
        
//        f1->tulip[oneVector].Address = fromBeginning(*f1,copyFourVector);
//        f1->tulip[oneVector].Partition = !(!c1->rt.printFlag) * 24*imax(*f1->oneBody.num, f1->twoBody.num);;
//        f1->tulip[oneVector].body = one;
//        f1->tulip[oneVector].species = vector;
//        f1->tulip[oneVector].header = Cube;
//        f1->tulip[oneVector].parallel = 2;
//
//        f1->tulip[twoVector].Address = fromBeginning(*f1,oneVector);
//        f1->tulip[twoVector].Partition = !(!c1->rt.printFlag) * f->i.bRank;;
//        f1->tulip[twoVector].body = two;
//        f1->tulip[twoVector].species = vector;
//        f1->tulip[twoVector].header = Cube;
//        f1->tulip[twoVector].parallel = 2;
//        f1->tulip[twoVector].purpose = Object;
////        f1->tulip[twoVector].symmetryType = nullSymmetry;
//        f1->tulip[twoVector].name = twoVector;
//
        fromBeginning(*f1,totalVector,copyFourVector);
        f1->tulip[totalVector].Partition = (!( c1->rt.phaseType == buildFoundation )|| 1)*  f->i.bRank*(c1->i.canonRank);
        f1->tulip[totalVector].species = vector;
        //f1->tulip[totalVector].spinor = parallel;

        fromBeginning(*f1,totalFuzzyVector,totalVector);
        f1->tulip[totalFuzzyVector].Partition = 0*(c1->i.canonRank);
        f1->tulip[totalFuzzyVector].species = vector;
        f1->tulip[totalFuzzyVector].spinor = parallel;

        fromBeginning(*f1,diagonalCube,totalFuzzyVector);
        f1->tulip[diagonalCube].Partition = 1;
        f1->tulip[diagonalCube].species = matrix;
        f1->tulip[diagonalCube].spinor = parallel;
        assignParticle(*f1, diagonalCube, all, one);
        
        fromBeginning(*f1,diagonal1VectorA,diagonalCube);
        f1->tulip[diagonal1VectorA].Partition =1;
        f1->tulip[diagonal1VectorA].spinor = parallel;
        f1->tulip[diagonal1VectorA].species = outerVector;
        assignParticle(*f1, diagonal1VectorA, all , one);
        
        fromBeginning(*f1,diagonal2VectorA,diagonal1VectorA);
        f1->tulip[diagonal2VectorA].Partition =(c1->rt.phaseType == buildFoundation);
        f1->tulip[diagonal2VectorA].spinor = parallel;
        f1->tulip[diagonal2VectorA].species = outerVector;
        assignParticle(*f1, diagonal2VectorA, all , two);

        fromBeginning(*f1,diagonal2VectorB,diagonal2VectorA);
        f1->tulip[diagonal2VectorB].Partition =(c1->rt.phaseType == buildFoundation);
        f1->tulip[diagonal2VectorB].spinor = parallel;
        f1->tulip[diagonal2VectorB].species = outerVector;
        assignParticle(*f1, diagonal2VectorB, all , two);

        fromBeginning(*f1,diagonal1VectorB,diagonal2VectorB);
        f1->tulip[diagonal1VectorB].Partition =(c1->rt.phaseType == buildFoundation);
        f1->tulip[diagonal1VectorB].spinor = parallel;
        f1->tulip[diagonal1VectorB].species = outerVector;
        assignParticle(*f1, diagonal1VectorB, all , one);

        fromBeginning(*f1,diagonal1VectorC,diagonal1VectorB);
        f1->tulip[diagonal1VectorC].Partition =(c1->rt.phaseType == buildFoundation);
        f1->tulip[diagonal1VectorC].spinor = parallel;
        f1->tulip[diagonal1VectorC].species = outerVector;
        assignParticle(*f1, diagonal1VectorC, all , one);

        fromBeginning(*f1,diagonal1VectorD,diagonal1VectorC);
        f1->tulip[diagonal1VectorD].Partition = (c1->rt.phaseType == buildFoundation);
        f1->tulip[diagonal1VectorD].spinor = parallel;
        f1->tulip[diagonal1VectorD].species = outerVector;
        assignParticle(*f1, diagonal1VectorD, all , one);

        fromBeginning(*f1,diagonal3VectorA,diagonal1VectorD);
        f1->tulip[diagonal3VectorA].Partition =(c1->rt.phaseType == buildFoundation);
       // f1->tulip[diagonal3VectorA].spinor = parallel;
        f1->tulip[diagonal3VectorA].species = outerVector;
        assignParticle(*f1, diagonal3VectorA, all , three);
        
        fromBeginning(*f1,bandBasis,diagonal3VectorA);
        f1->tulip[bandBasis].Partition = mxlen;
        f1->tulip[bandBasis].memory = bufferAllocation;

        
        fromBeginning(*f1,copy,bandBasis);
        f1->tulip[copy].Partition = imax(f->i.bRank*f->i.bRank,imax(mx1len*mx1len, c1->i.decomposeRankMatrix+ imax(c1->i.oneBody.num,   imax(c1->i.decomposeRankMatrix*c1->i.Na,outVector) )));
    if ( c1->rt.phaseType == frameDensity){
    f1->tulip[copy].spinor = parallel;
    }
        f1->tulip[copy].species = matrix;
        assignParticle(*f1, copy, all, one);
        for ( space = 0; space < SPACE ; space++)
            if ( f1->rose[space].body != nada)
                f1->tulip[copy].space[space].block = tv1;
        
        fromBeginning(*f1,copyTwo,copy);
        f1->tulip[copyTwo].Partition = c1->i.decomposeRankMatrix;
    if ( c1->rt.phaseType == frameDensity || c1->rt.phaseType == svdOperation ){
        f1->tulip[copyTwo].Partition = 10*2*c1->i.twoBody.num * f->i.qFloor * f->i.iRank *f->i.iRank;
    }
        f1->tulip[copyTwo].species = matrix;
        assignParticle(*f1, copyTwo, all, one);
        
        fromBeginning(*f1,copyThree,copyTwo);
        f1->tulip[copyThree].Partition = 0*c1->i.decomposeRankMatrix+ imax(c1->i.oneBody.num,   imax(c1->i.decomposeRankMatrix*c1->i.Na,outVector) );
        f1->tulip[copyThree].spinor = parallel;
        f1->tulip[copyThree].species = matrix;
        assignParticle(*f1, copyThree, all, one);
        
        fromBeginning(*f1,squareTwo,copyThree);
        f1->tulip[squareTwo].Partition= 0*(c1->rt.phaseType == buildFoundation);//BUILD
        f1->tulip[squareTwo].species = matrix;
        f1->tulip[squareTwo].header = Cube;
        assignParticle(*f1, squareTwo, all, two);
        
        fromBeginning(*f1,inversion,squareTwo);
        f1->tulip[inversion].Partition= 1;
        f1->tulip[inversion].species = matrix;
        f1->tulip[inversion].header = Cube;
        assignParticle(*f1, inversion, all, one);
    
        
        fromBeginning(*f1,canonicalBuffers,inversion);
        f1->tulip[canonicalBuffers].Partition = maxVector*maxVector+ imax(NV,maxVector*c1->i.canonRank)*maxVector;
        f1->tulip[canonicalBuffers].spinor = parallel;

        fromBeginning(*f1,trackBuffer,canonicalBuffers);
        f1->tulip[trackBuffer].Partition = (2*maxVector+1)*maxVector;
        f1->tulip[trackBuffer].spinor = parallel;
        f1->tulip[trackBuffer].memory = bufferAllocation;

        fromBeginning(*f1,guideBuffer,trackBuffer);
        f1->tulip[guideBuffer].Partition = c1->i.canonRank*maxVector*maxVector;
        f1->tulip[guideBuffer].spinor = parallel;
        f1->tulip[guideBuffer].memory = bufferAllocation;

        fromBeginning(*f1,foundationStructure,guideBuffer);
        f1->tulip[foundationStructure].spinor = parallel;
        f1->tulip[foundationStructure].Partition = 1;
        f1->tulip[foundationStructure].species = vector;
        f1->tulip[foundationStructure].spinor = cmpl;//need two channels
        
        fromBeginning(*f1,interactionExchange,foundationStructure);
            f1->tulip[interactionExchange].Partition = c1->i.twoBody.num*(( bootBodies > one )|| c1->rt.runFlag > 0);
            f1->tulip[interactionExchange].species = matrix;
            assignParticle(*f1, interactionExchange, electron, two);
    {
        {
            enum block ee ;
            for (ee = e12 ; ee <= e34 ; ee++){
                f1->tulip[interaction12+ee-e12].name = interactionExchange;
                f1->tulip[interaction12+ee-e12].species = matrix;
                if ( ee < e34 )
                    f1->tulip[interaction12+ee-e12].linkNext = interaction12+ee-e12+1;
                for ( space = 0; space < SPACE ; space++)
                    if ( f1->rose[space].body >= two && f1->rose[space].particle == electron )
                        f1->tulip[interaction12+ee-e12].space[space].block = ee;
            }
        }
    }
        fromBeginning(*f1,interactionExchangeB,interactionExchange);
        f1->tulip[interactionExchangeB].Partition = c1->i.twoBody.num*( bootBodies > one )* ( c1->rt.calcType == protonsElectronsCalculation);
        f1->tulip[interactionExchangeB].species = matrix;
        assignParticle(*f1, interactionExchangeB, proton, two);
        {
        enum block ee ;
            for (ee = e12 ; ee <= e34 ; ee++){
                f1->tulip[interaction12B+ee-e12].name = interactionExchangeB;
                f1->tulip[interaction12B+ee-e12].species = matrix;
                if ( ee < e34 )
                    f1->tulip[interaction12B+ee-e12].linkNext = interaction12B+ee-e12+1;
                for ( space = 0; space < SPACE ; space++)
                    if ( f1->rose[space].body >= two && f1->rose[space].particle == proton ){
                        f1->tulip[interaction12B+ee-e12].space[space].block = ee;
                    }
            }
        }
        
        fromBeginning(*f1,interactionEwald,interactionExchangeB);
        f1->tulip[interactionEwald].Partition = c1->i.twoBody.num * (c1->rt.runFlag > 0 );
        f1->tulip[interactionEwald].species = matrix;
        assignParticle(*f1, interactionEwald, electron, two);

    {
        enum block ee ;
        for (ee = e12 ; ee <= e34 ; ee++){
            f1->tulip[interaction12Ewald+ee-e12].name = interactionEwald;
            f1->tulip[interaction12Ewald+ee-e12].species = matrix;
            if ( ee < e34 )
                f1->tulip[interaction12Ewald+ee-e12].linkNext = interaction12Ewald+ee-e12+1;
            for ( space = 0; space < SPACE ; space++)
                if ( f1->rose[space].body >= two && f1->rose[space].particle == electron )
                    f1->tulip[interaction12Ewald+ee-e12].space[space].block = ee;
        }
    }
        fromBeginning(*f1,intracellularSelfEwald,interactionEwald);
        f1->tulip[intracellularSelfEwald].Partition = c1->i.decomposeRankMatrix;
        f1->tulip[intracellularSelfEwald].species = matrix;
        assignOneWithPointers(*f1, intracellularSelfEwald, electron);
    
        fromBeginning(*f1,intercellularSelfEwald,intracellularSelfEwald);
        f1->tulip[intercellularSelfEwald].Partition = c1->i.decomposeRankMatrix;
        f1->tulip[intercellularSelfEwald].species = matrix;
        assignOneWithPointers(*f1, intercellularSelfEwald, electron);
    
        fromBeginning(*f1,shortenPlus,intercellularSelfEwald);
        f1->tulip[shortenPlus].Partition = c1->i.twoBody.num*c1->i.decomposeRankMatrix*( c1->rt.calcType == clampProtonElectronCalculation );
        f1->tulip[shortenPlus].species = matrix;
        assignOneWithPointers(*f1, shortenPlus, all);
        
        fromBeginning(*f1,shortenMinus,shortenPlus);
        f1->tulip[shortenMinus].Partition = c1->i.twoBody.num*c1->i.decomposeRankMatrix*( c1->rt.calcType == clampProtonElectronCalculation );
        f1->tulip[shortenMinus].species = matrix;
        assignOneWithPointers(*f1, shortenMinus, all);
        
        fromBeginning(*f1,interactionExchangePlus,shortenMinus);
        f1->tulip[interactionExchangePlus].Partition =0* mx1len*mx1len*c1->i.twoBody.num*( c1->rt.calcType == clampProtonElectronCalculation );
        f1->tulip[interactionExchangePlus].species = matrix;
        assignOneWithPointers(*f1, interactionExchangePlus,all);
        
        fromBeginning(*f1,interactionExchangeMinus,interactionExchangePlus);
        f1->tulip[interactionExchangeMinus].Partition = 0*mx1len*mx1len*c1->i.twoBody.num*( c1->rt.calcType == clampProtonElectronCalculation );
        f1->tulip[interactionExchangeMinus].species = matrix;
        assignOneWithPointers(*f1, interactionExchangeMinus,all);
        
        fromBeginning(*f1,interactionTwoAcrossDimensions,interactionExchangeMinus);
        f1->tulip[interactionTwoAcrossDimensions].Partition = (!c1->i.decomposeRankMatrix)*mx1len*mx1len*c1->i.twoBody.num*( c1->rt.calcType == protonsElectronsCalculation );
        f1->tulip[interactionTwoAcrossDimensions].species = matrix;
        assignParticle(*f1, interactionTwoAcrossDimensions, all, one);

        {
            enum block ee,eeB ;
            for (ee = tv1 ; ee <= tv4 ; ee++)
                for (eeB = tv1 ; eeB <= tv4 ; eeB++){
                    
                f1->tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].name = interactionTwoAcrossDimensions;
                f1->tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].species = matrix;
                if ( ee < e34 && eeB < e34 )
                    f1->tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].linkNext = interactionTwoAcrossDimensions+ee-tv1+eeB-tv1+1;
                for ( space = 0; space < SPACE ; space++)
                    if ( space < 3)
                        f1->tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].space[space].block = ee;
                    else if ( space < 6 )
                        f1->tulip[interactionTwoAcrossDimensions+ee-tv1+eeB-tv1].space[space].block = eeB;
            }
        }
        
        fromBeginning(*f1,shortTwoAcrossDimensions,interactionTwoAcrossDimensions);
        f1->tulip[shortTwoAcrossDimensions].Partition = c1->i.twoBody.num*c1->i.decomposeRankMatrix*c1->i.twoBody.num*( c1->rt.calcType == protonsElectronsCalculation );
        f1->tulip[shortTwoAcrossDimensions].species = matrix;
        assignParticle(*f1, shortTwoAcrossDimensions, all, one);

        {
            enum block ee,eeB ;
            for (ee = tv1 ; ee <= tv4 ; ee++)
                for (eeB = tv1 ; eeB <= tv4 ; eeB++){
                    
                    f1->tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].name = shortTwoAcrossDimensions;
                    f1->tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].species = matrix;
                    if ( ee < e34 && eeB < e34 )
                        f1->tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].linkNext = shortTwoAcrossDimensions+ee-tv1+eeB-tv1+1;
                    for ( space = 0; space < SPACE ; space++)
                        if ( space < 3)
                            f1->tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].space[space].block = ee;
                        else if ( space < 6 )
                            f1->tulip[shortTwoAcrossDimensions+ee-tv1+eeB-tv1].space[space].block = eeB;
                }
        }

        
        fromBeginning(*f1,quadCube,shortTwoAcrossDimensions);
        f1->tulip[quadCube].Partition =1;
        f1->tulip[quadCube].species = matrix;
        assignParticle(*f1, quadCube, all, two);
    
        fromBeginning(*f1,oneArray,quadCube);
//        f1->tulip[oneArray].Partition = f->i.M1*3+f->i.M1*f->i.M1;
        f1->tulip[oneArray].memory = bufferAllocation;

        fromBeginning(*f1,threeArray,oneArray);
//        f1->tulip[threeArray].Partition = 2*f->i.M1*f->i.M1*f->i.M1;
        f1->tulip[threeArray].memory = bufferAllocation;

        fromBeginning(*f1,oneBasis,threeArray);
//        f1->tulip[oneBasis].Partition = f->i.M1*N1;
        f1->tulip[oneBasis].memory = bufferAllocation;

        INT_TYPE vecLen = 1;
        for ( space = 0  ; space < SPACE ; space++)
            if ( vecLen < alloc(*f1, eigenVectors, space))
                vecLen = alloc(*f1, eigenVectors, space);
        
        
        fromBeginning(*f1,tensorBuffers,oneBasis);
        f1->tulip[tensorBuffers].Partition = vecLen ;
        f1->tulip[tensorBuffers].spinor = parallel;
        f1->tulip[tensorBuffers].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers2,tensorBuffers);
        f1->tulip[tensorBuffers2].Partition = vecLen;
        f1->tulip[tensorBuffers2].spinor = parallel;
        f1->tulip[tensorBuffers2].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers3,tensorBuffers2);
        f1->tulip[tensorBuffers3].Partition = vecLen;
        f1->tulip[tensorBuffers3].spinor = parallel;
        f1->tulip[tensorBuffers3].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers4,tensorBuffers3);
        f1->tulip[tensorBuffers4].Partition = vecLen;
        f1->tulip[tensorBuffers4].spinor = parallel;
        f1->tulip[tensorBuffers4].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers5,tensorBuffers4);
        f1->tulip[tensorBuffers5].Partition = vecLen;
        f1->tulip[tensorBuffers5].spinor = parallel;
        f1->tulip[tensorBuffers5].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers6,tensorBuffers5);
        f1->tulip[tensorBuffers6].Partition = vecLen;
        f1->tulip[tensorBuffers6].spinor = parallel;
        f1->tulip[tensorBuffers6].memory = bufferAllocation;
        
        fromBeginning(*f1,oneByOneBuffer,tensorBuffers6);
        f1->tulip[oneByOneBuffer].Partition = mx1len*mx1len*mx1len*mx1len* (c1->rt.calcType >= clampProtonElectronCalculation);
        f1->tulip[oneByOneBuffer].memory = bufferAllocation;
        
        fromBeginning(*f1,canonicalBuffersB,oneByOneBuffer);
        f1->tulip[canonicalBuffersB].Partition = maxVector;
        f1->tulip[canonicalBuffersB].spinor = parallel;
        f1->tulip[canonicalBuffersB].memory = bufferAllocation;

        fromBeginning(*f1,canonicalBuffersBM,canonicalBuffersB);//twobody
        f1->tulip[canonicalBuffersBM].Partition = mx1len*mx1len*mx1len*mx1len ;
        f1->tulip[canonicalBuffersBM].memory = bufferAllocation;
    if ( c1->rt.phaseType == frameDensity    )
        f1->tulip[canonicalBuffersBM].spinor = parallel;

    
        fromBeginning(*f1,canonicalBuffersC,canonicalBuffersBM);
        f1->tulip[canonicalBuffersC].Partition = (c1->rt.phaseType == frameDensity ) *   NV;
        f1->tulip[canonicalBuffersC].spinor = parallel;
        f1->tulip[canonicalBuffersC].memory = bufferAllocation;

        fromBeginning(*f1,twoBodyRitz,canonicalBuffersC);
        f1->tulip[twoBodyRitz].Partition = maxArray;
        f1->tulip[twoBodyRitz].memory = bufferAllocation;

        fromBeginning(*f1,conditionOverlapNumbers,twoBodyRitz);
        f1->tulip[conditionOverlapNumbers].Partition = maxArray;
        f1->tulip[conditionOverlapNumbers].memory = bufferAllocation;
        
        fromBeginning(*f1,matrixHbuild,conditionOverlapNumbers);
    f1->tulip[matrixHbuild].Partition = (  c1->rt.phaseType == solveRitz|| c1->rt.phaseType == frameDensity ) *  (2*maxArray*maxArray);
        f1->tulip[matrixHbuild].memory = bufferAllocation;

        fromBeginning(*f1,vectorHbuild,matrixHbuild);
        f1->tulip[vectorHbuild].Partition = 0;
        f1->tulip[vectorHbuild].memory = bufferAllocation;

        fromBeginning(*f1,matrixSbuild,vectorHbuild);
        f1->tulip[matrixSbuild].Partition = 2*maxArray*maxArray;
        f1->tulip[matrixSbuild].memory = bufferAllocation;

    fromBeginning(*f1,square,matrixSbuild);
    f1->tulip[square].Partition=   4*maxVector*maxVector;
    f1->tulip[square].species = matrix;
    assignParticle(*f1, square, all, one);
    
    
    
        fromBeginning(*f1,dsyBuffers,square);
#if 1
        f1->tulip[dsyBuffers].Partition = 2*8*(8*(imax(mxlen,maxEV))+72*f->i.nStates*f->i.nStates+ 8 * mxlen)+3*maxEV;
#else
        f1->tulip[dsyBuffers].Partition = maxVector*maxVector;
#endif
        f1->tulip[dsyBuffers].spinor = parallel;
        f1->tulip[dsyBuffers].memory = bufferAllocation;

        fromBeginning(*f1,end,dsyBuffers);
        fromBeginning(*f1,end,end);
        struct name_label e = f1->tulip[end];
        struct name_label k1 = f1->tulip[kinetic1];
        struct name_label k2 = f1->tulip[kinetic2];
        struct name_label it = f1->tulip[interaction12];
        struct name_label i = f1->tulip[interactionExchangePlus];
        struct name_label ip1 = f1->tulip[interaction1Plus];

    
        {
            printf("\t| SPACE \t:   Gb\t \n");
            double maxMem = 0.,currMem;
            for ( space = 0 ; space <= SPACE ; space++){
                currMem = (f1->tulip[end].space[space].Address)/(1000000000./(sizeof(Stream_Type)));
                
                
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
            f1->rose[space].stream = malloc( (f1->tulip[end].space[space].Address)*sizeof(Stream_Type));
        }

        f1->bootedMemory = 1;
        assignCores(*f1, 1);

        INT_TYPE RdsSize;
        RdsSize = 0;

        
//        if (f->i.M1){
//            INT_TYPE i,ii;
//            INT_TYPE M1 = f->i.M1,M12 = (M1-1)/2;
//            double r = (double)(M1-1)/(double)(N1-1);
//            double * u =myStreams(*f1,oneBasis,0);
//            for( i = 0; i < N1 ; i++)
//                for ( ii = 0 ; ii < M1 ; ii++){
//                    u[i*M1+ii] = Sinc(r, r*(i-N12)- (ii-M12))/sqrt(r);
//                }
//        }
        
    
        
        if(f1->boot == fullMatrices){
            f1->tulip[Ha].species    = scalar;
            f1->tulip[Iterator].species    = scalar;
            if ( c1->rt.calcType == electronicStuctureCalculation  ){
                separateKinetic(*f1, 0,kinetic, c1->i.massElectron,electron);
                separateOverlap(*f1, 0,overlap, 0,all);

                if ( f1->rose[0].component == periodicComponent1 ){
                    printf("func %d\n",c1->i.twoBody.func.fn);
                    if ( c1->i.twoBody.func.fn != nullFunction)
                    if ( ! ioStoreMatrix(*f1, intracellularSelfEwald, 0, "intracellularEwald.matrix", 1)){
                        printf("running without ewald-self interaction\n");
                        tClear(*f1,intracellularSelfEwald);
                    }
                    if ( ! ioStoreMatrix(*f1, intercellularSelfEwald, 0, "intercellularEwald.matrix", 1)){
                        printf("running without ewald-self interaction\n");
                        tClear(*f1,intercellularSelfEwald);
                    }
                    if (c1->i.twoBody.func.fn != nullFunction )
                    if ( !  ioStoreMatrix(*f1, interactionEwald, 0, "interactionEwald.matrix",1) ){
                        if (! (c1->rt.phaseType == frameDensity || c1->rt.phaseType == solveRitz) ){
                            printf("you can remove this barrier, but I would recommend you run ritz or frame first\n");
                            exit(0);
                        }
                        mySeparateExactTwo(*f1, c1->i.twoBody,interactionEwald, 1., 0, 1, electron);
                    }
//                    if ( TEST4 )
//                        buildElectronProtonInteraction(*f1, linear, 0);
                    if ( c1->i.twoBody.func.fn != nullFunction&& c1->rt.phaseType == frameDensity)
                    if ( ! ioStoreMatrix(*f1, interactionExchange, 0, "interactionExchange.matrix",1) ){
                        
                        if (c1->rt.phaseType != frameDensity ){
                            printf("you can remove this barrier, but I would recommend you run frame first\n");
                            exit(0);
                        }

//                        printf("skipping zero-cell intrainteraction\n");
//                        if (0)
                        
                        mySeparateExactTwo(*f1,c1->i.twoBody, interactionExchange, 1., 0, 0, electron);
                        ioStoreMatrix(*f1, interactionExchange, 0, "interactionExchange.matrix", 0);
                    }
                    separateExternal(c1,*f1,linear,0,0,-1.0,-1,0,electron);
                }else {
//                    if ( f1->twoBody.func.fn != nullFunction)
                    if ( bootBodies > one ){
                        if ( c1->i.twoBody.func.fn != nullFunction )
                            if(  ! ioStoreMatrix(*f1, interactionExchange, 0, "interactionExchange.matrix",1)){
                            if ( c1->rt.phaseType != solveRitz ){
                                printf("you can remove this barrier, but I would recommend you run ritz first\n");
                                exit(0);
                            }

                            mySeparateExactTwo(*f1,c1->i.twoBody, interactionExchange, 1., 0, 0, electron);
                        }
                    }
                    
                    separateExternal(c1,*f1,linear,0,0,-1.0,-1,0,electron);
                    
                }

                if ( c1->i.magFlag ){
                    INT_TYPE deriv[SPACE];
                    INT_TYPE power[SPACE];
                    
                    
                    deriv[0] = 1;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 0 ;
                    power[1] = 1;
                    power[2] = 0;
                    
                    separateDerivatives(*f1, 0, vectorMomentum, power, deriv, 0.5*c1->i.mag, electron);
                    
                    deriv[0] = 0;
                    deriv[1] = 1;
                    deriv[2] = 0;
                    
                    power[0] = 1;
                    power[1] = 0;
                    power[2] = 0;
                    
                    separateDerivatives(*f1,0, vectorMomentum, power, deriv, -0.5*c1->i.mag, electron);

                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 2;
                    power[1] = 0;
                    power[2] = 0;
                    
                    separateDerivatives(*f1, 0, vectorMomentum, power, deriv, 0.25*c1->i.mag, electron);

                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 0;
                    power[1] = 2;
                    power[2] = 0;
                    
                    separateDerivatives(*f1, 0, vectorMomentum, power, deriv, 0.25*c1->i.mag, electron);
                    

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
                    
                    separateDerivatives(*f1, 0, vectorMomentum, power, deriv, 0.5*c1->i.springConstant, electron);
                    
                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 0;
                    power[1] = 2;
                    power[2] = 0;
                    
                    separateDerivatives(*f1,0, vectorMomentum, power, deriv,  0.5*c1->i.springConstant, electron);
                    
                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 0;
                    power[1] = 0;
                    power[2] = 2;
                    
                    separateDerivatives(*f1, 0, vectorMomentum, power, deriv,  0.5*c1->i.springConstant, electron);
                    
                
                }
                

        } else if ( c1->rt.calcType == clampProtonElectronCalculation  ){
                if ( bootBodies > one ){
                    mySeparateExactTwo(*f1,c1->i.twoBody,interactionExchange, 1. , 0,0,electron);
                }
            
            if ( c1->i.twoBody.func.fn != nullFunction && !ioStoreMatrix(*f1,shortenPlus ,0,"shortenExchangePlus.matrix",1) ){
                if ( c1->rt.phaseType != solveRitz ){
                    printf("you can remove this barrier, but I would recommend you run ritz first\n");
                    exit(0);
                }

                    mySeparateExactOneByOne(*f1,c1->i.twoBody, c1->i.decomposeRankMatrix,interactionExchangePlus, shortenPlus,-1., 1,c1->i.massClampPair/(c1->i.massProton + c1->i.massClampPair),electron, proton);
            }
            if ( c1->i.twoBody.func.fn != nullFunction&& ! ioStoreMatrix(*f1,shortenMinus ,0,"shortenExchangeMinus.matrix",1) ){
                if ( c1->rt.phaseType != solveRitz ){
                    printf("you can remove this barrier, but I would recommend you run ritz first\n");
                    exit(0);
                }

                    mySeparateExactOneByOne(*f1,c1->i.twoBody, c1->i.decomposeRankMatrix,interactionExchangeMinus, shortenMinus,-1., -1,c1->i.massProton/(c1->i.massProton + c1->i.massClampPair),electron, proton);

            }
            struct name_label u = f1->tulip[shortenPlus];
            
            if ( c1->i.oneBody.func.fn != nullFunction&&c1->rt.phaseType != buildFoundation )
                separateExternal(c1,*f1,protonRepulsion, 0,0,0.5,-1,0,proton);
            separateKinetic(*f1, 0,kineticMass,4* c1->i.massProton*c1->i.massClampPair/(c1->i.massProton+c1->i.massClampPair),proton);

            separateKinetic(*f1, 0,kinetic, c1->i.massElectron*(c1->i.massProton + c1->i.massClampPair )/(c1->i.massProton + c1->i.massClampPair + bootBodies*c1->i.massElectron ),electron);

            }
            else if ( c1->rt.calcType == protonsElectronsCalculation  ){

                if ( bootBodies > one ){
                    mySeparateExactTwo(*f1,c1->i.twoBody,interactionExchange, 1. , 0,0,electron);
                    mySeparateExactTwo(*f1,c1->i.twoBody,interactionExchangeB, 1. , 0,0,proton);
                //    tZeroSum(*f1, interactionExchange, 0 );
                //    tZeroSum(*f1, interactionExchangeB, 0 );
                }
                mySeparateExactOneByOne(*f1,c1->i.twoBody,c1->i.decomposeRankMatrix,interactionTwoAcrossDimensions, shortTwoAcrossDimensions,-1., 1,1,electron, proton);
                separateKinetic(*f1, 0,kineticMass, c1->i.massProton,proton);
                separateKinetic(*f1, 0,kinetic, c1->i.massElectron,electron);

            }


        
    
        {
            if ( f->i.nOperator ){

                f1->tulip[Ha].linkNext = vectorOperator;
                f1->tulip[Iterator].linkNext = vectorOperator;
                INT_TYPE fi,fe=0,ff;
                for ( fi =0 ; fi < f->i.filesVectorOperator ; fi++){
                    tLoadEigenWeights (c1,*f,f->i.fileVectorOperator[fi], &fe,f1->vectorOperator,0);
                }
                
                for ( ff = 1; ff < fe ; ff++)
                    f1->tulip[vectorOperator+ff-1].linkNext = vectorOperator+ff;

            
                
                if ( fe > f->i.nOperator )
                {
                    printf("failure to load vector Operators");
                    exit(0);
                }else {
                    f1->vectorOperator = fe;
                }
                
            }
            else
            if ( bootBodies == one && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->tulip[intracellularSelf1Ewald].linkNext = intercellularSelf1Ewald;
                f1->tulip[intercellularSelf1Ewald].linkNext = vectorMomentum1;
                f1->tulip[vectorMomentum1].linkNext = kinetic1;
                f1->tulip[kinetic1].linkNext = external1;
                f1->tulip[external1].linkNext = nullName;

                //active assignment
                f1->tulip[Ha].linkNext = intracellularSelf1Ewald;
                f1->tulip[Iterator].linkNext = intracellularSelf1Ewald;
                
            } else if ( bootBodies == two && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->tulip[intracellularSelf2Ewald].linkNext = intercellularSelf1Ewald;
                f1->tulip[intercellularSelf2Ewald].linkNext = vectorMomentum1;
                f1->tulip[vectorMomentum2].linkNext = kinetic1;
                f1->tulip[kinetic2].linkNext = external1;

                if ( c1->rt.runFlag == 0 ){
                    f1->tulip[external2].linkNext = interaction12;
                    f1->tulip[interaction12].linkNext = nullName;

                }else {
                    f1->tulip[external2].linkNext = interaction12Ewald;
                    f1->tulip[interaction12Ewald].linkNext = nullName;
                }
                //active assignment
                f1->tulip[Ha].linkNext = intracellularSelf1Ewald;
                f1->tulip[Iterator].linkNext = intracellularSelf1Ewald;
                
            }else if ( bootBodies == three && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->tulip[intracellularSelf3Ewald].linkNext = intercellularSelf1Ewald;
                f1->tulip[intercellularSelf3Ewald].linkNext = vectorMomentum1;

                f1->tulip[vectorMomentum3].linkNext = kinetic1;
                f1->tulip[kinetic3].linkNext = external1;

                if ( c1->rt.runFlag == 0 ){
                    f1->tulip[external3].linkNext = interaction12;
                    f1->tulip[interaction23].linkNext = nullName;
                }else {
                    f1->tulip[external3].linkNext = interaction12Ewald;
                    f1->tulip[interaction23Ewald].linkNext = nullName;
                }
                //active assignment
                    f1->tulip[Ha].linkNext = intracellularSelf1Ewald;
                    f1->tulip[Iterator].linkNext = intracellularSelf1Ewald;
                
            }
            
            else if ( bootBodies == four && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->tulip[intracellularSelf4Ewald].linkNext = intercellularSelf1Ewald;
                f1->tulip[intercellularSelf4Ewald].linkNext = vectorMomentum1;

                f1->tulip[vectorMomentum4].linkNext = kinetic1;
                f1->tulip[kinetic4].linkNext = external1;

                if ( c1->rt.runFlag == 0 ){
                    f1->tulip[external4].linkNext = interaction12;
                    f1->tulip[interaction34].linkNext = nullName;
                }else {
                    f1->tulip[external4].linkNext = interaction12Ewald;
                    f1->tulip[interaction34Ewald].linkNext = nullName;
                }

                //active assignment
                f1->tulip[Ha].linkNext = intracellularSelf1Ewald;
                f1->tulip[Iterator].linkNext = intracellularSelf1Ewald;
                
            }
            else if (bootBodies == one && ( c1->rt.calcType == clampProtonElectronCalculation  )){
                //paths do not matter unless assigned below...
                f1->tulip[shorten1Plus].linkNext =shorten1Minus;
                f1->tulip[shorten1Minus].linkNext =kinetic1;
                
                f1->tulip[kinetic1].linkNext = kineticMass1;
                f1->tulip[kineticMass1].linkNext = proton1;
                f1->tulip[proton1].linkNext = nullName;
                //active assignment
                f1->tulip[Ha].linkNext = shorten1Plus;
                f1->tulip[Iterator].linkNext = shorten1Plus;
            }
            else if (bootBodies == two && ( c1->rt.calcType == clampProtonElectronCalculation  )){
                f1->tulip[interaction12].linkNext =interaction1Plus;
                f1->tulip[interaction2Plus].linkNext =interaction1Minus;
                f1->tulip[interaction2Minus].linkNext =kinetic1;

                f1->tulip[kinetic2].linkNext = kineticMass1;
                f1->tulip[kineticMass1].linkNext = proton1;
                
                
                struct name_label u = f1->tulip[kinetic1];

                //active assignment
                f1->tulip[Ha].linkNext = interaction12;
                f1->tulip[Iterator].linkNext = interaction12;
                
            }
//            else if (bootBodies == two && ( c1->rt.calcType == protonsElectronsCalculation  )){
//                f1->tulip[interaction12].linkNext =interaction1Plus;
//
//                f1->tulip[interaction3Plus].linkNext =interaction1Minus;
//                f1->tulip[interaction3Minus].linkNext =kinetic1;
//                
//                f1->tulip[kinetic3].linkNext = kineticMass1;
//                f1->tulip[kineticMass1].linkNext = proton1;
//
//                //active assignment
//                f1->tulip[Ha].linkNext = interaction12;
//                f1->tulip[Iterator].linkNext = interaction12;
//                
//            }

            
            
        }
        {
            INT_TYPE nn[SPACE],i,j;
            length1(*f1,nn);
            for ( space = 0; space < SPACE ; space++)
                if ( f1->rose[space].body != nada )
                    for ( i = 0 ; i < nn[space];i++)
                        for ( j = 0; j < nn[space];j++)
                        {
                            if ( i+j+1 == nn[space] )
                                streams(*f1, inversion,0,space)[i*nn[space]+j] = 1.;
                            else
                                streams(*f1, inversion,0,space)[i*nn[space]+j] = 0.;
                            
                        }
            
        }
        }
    
   // tInnerTest(*f1, kinetic, copy);
#if VERBOSE
    printf("boot complete\n");
#endif
    fflush(stdout);
    return 0;
}
