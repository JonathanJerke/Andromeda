/*
 *  Model.c
 *
 *
 *  Copyright 2020 Jonathan Jerke and Bill Poirier.
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

    i.i.OpIndex = 0;
    i.i.bRank = 0;
    i.i.iRank = 0;
    i.i.xRank = 0;
    i.i.filter = 0;
    i.i.irrep = 0;
    i.i.cat = 0;
    i.i.collect = 0;
    
    i.i.Iterations = 0;
    i.f.boot = noMatrices;
    for( space = 0 ; space <= SPACE ; space++){
        i.f.rose[space].stream = NULL;
     //   i.f.rose[space].basisList = NULL;
    }
    i.f.bootedMemory = 0;
    i.f.tulip = NULL;
    i.i.files = 0;
    i.i.file2Fill= 0;
    i.i.filesVectorOperator = 0;
    i.f.nullLabels.currLabel = 0;
    i.f.eikonLabels.currLabel = 0;
    i.f.nullLabels.maxLabel = MAXNAMES;
    i.f.eikonLabels.maxLabel = MAXNAMES;

#ifdef OVERFLAG
        i.i.cmpl = cmpl;
        i.i.bRank = 4;
        i.i.iRank = 1;
        i.i.nStates = 1;
        i.i.qFloor = 9*9*9;
        i.i.filter = 0;
        i.f.boot = fullMatrices;
        i.i.body = one;
        i.i.irrep = 1;
        i.i.cat  = 1;
        i.i.epi = 2;
        i.i.d = 1;
#else
    i.i.OpIndex =1;
    i.i.qFloor = 1;

    i.i.d = 0.35 ;
    i.i.D = 0.1*2;
    i.i.cmpl = real;
    i.i.bRank = 35;
    i.i.iRank = 1;
    i.i.nStates = 1;
    i.i.qFloor = 0;
    i.i.filter = 0;
    i.f.boot  = fullMatrices;
    i.i.body  = one;
    i.i.irrep = 0;
    i.i.cat   = 0;
    i.i.epi   = 4;
    i.i.around=0;
    i.i.eikonFlag = 1;
    i.i.chainFlag = 1;

#endif
    return i;
}
struct calculation initCal (void ) {
    struct calculation i;
    i.i.minIterationPrint = 0;
    i.i.flipSignFlag = 0;
    i.i.shiftFlag = 0;
    i.i.termNumber = 0;
    INT_TYPE n;
    for ( n= 0 ; n < 100 ; n++){
        i.i.shiftVector[n][0] = -50;
        i.i.shiftVector[n][1] = 1;
    }
    i.rt.GAMMA = 100;
    i.rt.BETA = 0.0001;
#ifdef APPLE
    resetA(&i.rt);
    blockA(&i.rt, 1);
    blockA(&i.rt, 2);
    blockA(&i.rt, 3);
    blockA(&i.rt, 4);
    blockA(&i.rt, 5);
    blockA(&i.rt, 6);
   // blockA(&i.rt, 7);
    blockA(&i.rt, 8);
    //blockA(&i.rt, 9);
    blockA(&i.rt, 10);
    blockA(&i.rt, 11);
    blockA(&i.rt, 12);
    blockA(&i.rt, 13);

#ifdef CHROME
    i.i.chromaticRank = 1000;
    i.i.chromos = 0.9999;
    i.i.chroma = 1;
#endif
    //i.i.OCSBflag = 0;
    i.i.springConstant = 0.;
    i.i.springFlag = 0;
    i.i.RAMmax = 6;
    i.rt.runFlag = 0;
    i.i.vectorMomentum = 0.;
    i.i.decomposeRankMatrix = 0;
    i.i.orgClamp = 2.;
    i.i.Angstroms = 0;
    i.rt.TARGET = 1e-4;
    i.rt.ALPHA = 1e-9;
    i.rt.CANON = 1e-7;
    i.rt.vCANON = 1e-3;
    i.rt.TOL = 1e5;
    i.rt.maxEntropy = 1;
    i.i.level = 3100;
//    if ( SPACE == 1 ){
//        i.i.M1 = 0;
//        i.i.Na = 0;
//        i.i.level = 100;
//    }else
    {
    i.i.M1 = 0;
    }
    
    i.rt.powDecompose = 3;
    
        i.i.magFlag = 0;
        i.i.mag = 0.1;
        //THESE
    
    i.rt.calcType = electronicStuctureCalculation;
    i.rt.runFlag = 0;
    i.rt.phaseType = singleKrylov;
    i.i.gaussCount = 0;
    //i.i.Na = 1;

        if ( SPACE == 3 || SPACE == 1 ){
            i.rt.calcType = electronicStuctureCalculation;
            i.rt.runFlag = 0;
        } else if ( SPACE == 6 ){

            i.rt.calcType = clampProtonElectronCalculation;
            i.i.orgClamp = 2;
        }else if ( SPACE == 1 ){
            i.rt.calcType = electronicStuctureCalculation;
            i.rt.runFlag = 0;
            i.rt.phaseType = buildFoundation;

        }
    i.i.massElectron = 1.;
    i.i.massProton = 1836.15267245;
    i.i.massClampPair = 1836.15267245;
    resetExternal(&i, 1, 1);
    
    i.i.twoBody.func.fn = Coulomb;
    i.i.oneBody.func.fn = Coulomb;
    //i.i.springFlag = 1;
    i.i.springConstant = 0.25;
    i.i.canonRank = 100 ;
    i.i.twoBody.num = 30;
    i.i.twoBody.func.contr = 2;
    i.i.twoBody.func.interval  = 1;
    i.i.twoBody.func.param[0]  = 1;
    i.i.twoBody.func.param[1]  = 1;
    i.i.twoBody.func.param[2]  = 2;

    i.i.oneBody.func.contr = 1;
    i.i.oneBody.num = 7;
    i.i.oneBody.func.interval  = 1;
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
    i.rt.ALPHA = 1e-8;
    i.rt.CANON = 1e-7;
    i.rt.vCANON = 1e-6;
    i.rt.TOL = 1e5;
    i.rt.maxEntropy = 1;
    i.i.level = 10000;
    i.i.level = 1;
    i.rt.powDecompose = 3;
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
    i.i.twoBody.func.contr = 0;
    i.i.twoBody.func.interval  = 0;
    i.i.twoBody.func.param[0]  = 1;
    i.i.twoBody.func.param[1]  = 1;
    i.i.twoBody.func.param[2]  = 1;
    
    i.i.oneBody.num = 50;
    i.i.oneBody.func.contr = 0;
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
         //       free(f1->rose[i].basisList);
           //     f1->rose[i].basisList = NULL;
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
                //printf("\nperiodic boundaries %d\n",f1->rose[space].component);
            }
            periodic *= 2;
            f1->rose[space].count1Basis = N1;
            f1->rose[space].lattice = f->i.d;
            if (f1->rose[space].component > 3  )
                f1->rose[space].count1Basis *= 2;
          //  f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);

//            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
//                f1->rose[space].basisList[bl] = defineSincBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, f->i.d, 0.,             f1->rose[space].count1Basis,bl);

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
               // printf("\nperiodic boundaries %d\n",f1->rose[space].component);
            }

            f1->rose[space].count1Basis = N1;
//            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
//            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
//                f1->rose[space].basisList[bl] = defineSincBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, f->i.d, 0.  ,             f1->rose[space].count1Basis,bl);


        }
        periodic = 1;
        for ( space = COMPONENT; space <= COMPONENT ; space++){
            f1->rose[space].particle = proton;
            f1->rose[space].basis = SincBasisElement;
            f1->rose[space].body = one;
            f1->rose[space].component = spatialComponent1+space%COMPONENT+ ((c1->rt.runFlag/periodic)%2)*3;
            periodic *= 2;
            if (  f1->rose[space].component > 3){
               // printf("\nperiodic boundaries %d\n",f1->rose[space].component);
            }

            f1->rose[space].count1Basis = N2;
//            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
//            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
//                f1->rose[space].basisList[bl] = defineSincBasis(nullNote, f1->rose[space].component, f1->rose[space].basis, f->i.D, c1->i.orgClamp,             f1->rose[space].count1Basis,bl);


        }
        for ( space = COMPONENT+1; space < SPACE ; space++){
            f1->rose[space].particle = proton;
            f1->rose[space].basis = SincBasisElement;
            f1->rose[space].body = nada;
            f1->rose[space].component = spatialComponent1+space%COMPONENT+ ((c1->rt.runFlag/periodic)%2)*3;
            periodic *= 2;
            f1->rose[space].count1Basis = 0;
//            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*1);
            
        }
        
    }
    return N1;
}

#ifdef PURITY
INT_TYPE definePurity(struct sinc_label *f1, INT_TYPE R,INT_TYPE can, enum division head){
    if ( head == nullName ){
        return 3*R*R;
    }
    else {
        enum division last = head;
        INT_TYPE r;
        f1->temp = f1->purity +R*R;
        f1->purityOverlap = f1->purity+2*R*R;

        for ( r = 0 ; r < 3*R*R ; r++){
   fromBeginning(*f1,f1->purity+ r,last);
            last = f1->purity + r;
            if ( r < 2*R*R )
                f1->tulip[f1->purity+ r].Partition = can;
            else
                f1->tulip[f1->purity+ r].Partition = 2;
            f1->tulip[f1->purity+ r].species = matrix;
            assignParticle(*f1, f1->purity+ r, electron, one);
        }
        f1->purityCanon = R;
    }
    return 3*R*R;
}
#else

INT_TYPE definePurity(struct sinc_label *f1, INT_TYPE R,INT_TYPE can, enum division head){
    return 0;
}
#endif


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
            f1->rose[space].body = 1;
            f1->rose[space].component = spatialComponent1+space%COMPONENT ;
            f1->rose[space].count1Basis = 3;//change
//            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
//            for ( bl = 0 ; bl < f1->rose[space].count1Basis  ; bl++)
//                f1->rose[space].basisList[bl] = defineGaussBasis(nullNote, f1->rose[space].component, GaussianBasisElement, ble[bl], 0.,             f1->rose[space].count1Basis,0);
            
        }
    }
    return 1;
}

INT_TYPE singleTestModel( struct calculation * c1, struct field * f){
    struct sinc_label *f1 = &f->f;
    INT_TYPE bl,B1 = c1->i.gaussCount,periodic = 1,space, N1 = f->i.epi*2+1;
    if ( c1->i.gaussCount > 10){
        printf("gau");
        exit(0);
    }
    if ( B1 ){
        printf("GaussianSinc-Target 4\n");
        printf("ACKTUNG!!!!\n\n,  only one set of %d GTOs and on only the origin!\n\nACKTUNG!!!!\n\n",B1);
    }
    double ble[] = {2.8, 11.2, 44.8, 179.2, 716.8, 2867.2, 11468.8, 45875.2, 183501.,734003.};
    //
    //define vectors
    f1->rose[SPACE].component = nullComponent;
    f1->rose[SPACE].body = nada;
    f1->rose[SPACE].particle = nullParticle;
        for ( space = 0; space < SPACE ; space++){
            f1->rose[space].particle = electron;
            f1->rose[space].body = f->i.body;
            f1->rose[space].component = spatialComponent1+space%COMPONENT ;
            f1->rose[space].count1Basis = N1+B1;
//            f1->rose[space].basisList = malloc(sizeof(struct basisElement)*f1->rose[space].count1Basis);
//            for ( bl = 0; bl < B1  ; bl++){
//                f1->rose[space].basisList[bl] = defineGaussBasis(nullNote, f1->rose[space].component, GaussianBasisElement, ble[bl]/(f->i.d)/( f->i.d ), 0.,             f1->rose[space].count1Basis,0);
//            }
//            for ( bl = 0 ; bl < N1  ; bl++){
//                f1->rose[space].basisList[bl+B1] = defineSincBasis(nullNote, f1->rose[space].component, SincBasisElement, f->i.d, 0.,             f1->rose[space].count1Basis-B1,bl);
//            }
        }
    
return 1;
}


INT_TYPE iModel( struct calculation * c1, struct field *f){
    enum spinType c;
    f->f.eikonFlag = f->i.eikonFlag;
    f->f.chainFlag = f->i.chainFlag;
    INT_TYPE splitOperator = 1 ;
#ifdef GAUSSIANSINC
        singleTestModel(c1, f);

#else
        singleSincModel(c1, f);
#endif
    struct sinc_label *f1 = &f->f;

    {//SA++
    //    f->f.cat = f->i.cat;
        f->f.irrep = f->i.irrep;
        f->f.cmpl = f->i.cmpl;
    }//SA++
    
    if(f1->boot == fullMatrices){
        printf("\n\n*Parameters\n");
        printf("\t target\t\t\t%1.1f\n", -log(c1->rt.TARGET)/log(10));//quality of decomposition
        printf("\t tolerance\t\t%1.1f\n", log(c1->rt.TOL)/log(10));//max condition of foundation    s
        printf("\t threshold\t\t%1.1f\n", -log(c1->rt.CANON)/log(10));//matrix training standard
        printf("\t vectorThreshold\t%1.1f\n", -log(c1->rt.vCANON)/log(10));//vector training standard
        printf("\t condition\t\t%1.1f\n", -log(c1->rt.ALPHA )/log(10));//Beylkin parameter
        printf("\t relativeThreshold\t\t%1.1f\n", -log(c1->rt.BETA )/log(10));//
        printf("\t minALS\t\t%1.1f\n", c1->rt.GAMMA);//

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
        INT_TYPE maxVector = imax(0*c1->i.decomposeRankMatrix, imax(f->i.bRank,1+f->i.iRank));
        //rds defined in input.c
        
        bootShape = Cube+c1->i.gaussCount;
        

    INT_TYPE FloorV = imax(0, f->i.qFloor), CeilV = imax(0,0);
        INT_TYPE maxArray,EV,maxEV,NV = 0,FV = FloorV+CeilV ;
    
        EV = FV;
    maxEV =EV*(imax(f->i.Iterations,1));
    maxArray = imax(f->i.nStates,maxEV);//slip Nb into spectra...
        
        f1->maxEV = maxArray;
        enum division vectorOperator  = eigenVectors + 1 +  f->i.nStates+maxEV;
        f1->vectorOperator = vectorOperator;
    {
        INT_TYPE ir,ix;
        f->i.nOperator = countLinesFromFile(c1, *f, 1, &ir, &ix);
        if (f->i.nOperator ){
            printf("Operators\t %d\n", f->i.nOperator);
        }
        
        INT_TYPE rx = 1;
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
        }
        if ( !splitOperator )
            rx = 0;
    
        enum division end = vectorOperator+f->i.nOperator*(2*rx+1)+f1->nullLabels.maxLabel + f1->eikonLabels.maxLabel;
        //EIKONS
        //EIKONS
        f1->nullLabels.head = vectorOperator+f->i.nOperator*(2*rx+1);
        //EIKONS
        f1->eikonLabels.head  = f1->nullLabels.head +f1->nullLabels.maxLabel;

    
        f1->end = end;
        f1->tulip = malloc ( (end+1) * sizeof(struct name_label));
    }
        INT_TYPE outVector  = imax( c1->i.Na,1) * c1->i.decomposeRankMatrix *maxVector;
        
        {//defaults
            //define vectors end
            {
                INT_TYPE c;
                enum division label1;
                for ( label1 = 0 ;label1 <= f1->end; label1++){
                    f1->tulip[label1].name = label1;
                    f1->tulip[label1].Partition = 0;
                    f1->tulip[label1].header = bootShape;
                    f1->tulip[label1].spinor = f1->cmpl;
                    f1->tulip[label1].species = scalar;
                    f1->tulip[label1].linkNext = nullName;
                    f1->tulip[label1].chainNext = nullName;
                    f1->tulip[label1].loopNext = nullName;
                    f1->tulip[label1].memory = objectAllocation;
                   // f1->tulip[label1].operatorSignFlag = 0;
                    for ( space = 0; space <= SPACE ; space++)
                        f1->tulip[label1].space[space].Address = -1;
                    f1->tulip[label1].space[SPACE].block = id0;
                    for ( space = 0; space <= SPACE ; space++){
                        {
                            f1->tulip[label1].space[space].mapTo = space;
                            f1->tulip[label1].space[space].block = id0;//matrix prototype
                            f1->tulip[label1].space[space].body = nada;//matrix prototype
                            f1->tulip[label1].space[space].act = 1;
                            f1->tulip[label1].space[space].invert = 0;
                        }
                    }
                    tClear(*f1,label1);
                    for ( c = 0 ; c < MaxCore ; c++)
                        f1->tulip[label1].Begin[c] = 0;
                }
                
            }
        }//defaults
    if(0){
        INT_TYPE i ;
        for ( i = 0; i <= 6 ; i++){
            printf("%d %d\n",i, allowQ(f1->rt,i));
        }
    }
    struct name_label *u = &f->f.tulip[intracellularSelfEwald];
    {
        
        {
            enum division prev = 0;

            if ( f->i.OpIndex ){
            INT_TYPE terms = defineTerms(c1, f1,Iterator,0);
            if ( terms){
                printf("yo, there are %d terms,\n",terms);
            if ( maxVector % terms != 0 ){
                printf("increasing basisRank by %d\n" ,terms - maxVector % terms );
                maxVector += terms - maxVector % terms;
            }
            }
            INT_TYPE term;
            
            
            
           //krylov            for (term = 0; term < terms*terms; term++){
            if ( c1->rt.phaseType == productKrylov && f->i.Iterations > 1)
            for (term = 0; term < terms*terms ; term++){
                enum division brat = anotherLabel(f1, all, nada);
                f1->tulip[brat].Partition = maxVector;//HERE
                    
                f1->tulip[brat].spinor = f1->cmpl;
                if ( c1->rt.phaseType == productKrylov )
                    f1->tulip[brat].species = vector;
                else
                if ( c1->rt.phaseType == singleKrylov)
                {
                    f1->tulip[brat].species = outerVector;
                    assignParticle(*f1, brat, all, f->i.bodyTo);
                }
                fromBeginning(*f1, brat, prev);
                if ( ! prev )
                    f1->bra = brat;
                prev = brat;
            }
        
            }
            fromBeginning(*f1, kinetic1, prev);

        }
        
        
        
        assignOneWithPointers(*f1, kinetic, all);
        f1->tulip[kinetic1].Partition = allowQ(f1->rt,blockKineticEnergyBlock)*COMPONENT;//
        f1->tulip[kinetic1].name = kinetic1;
        assignParticle(*f1, kinetic1, all, one);
        tClear(*f1, kinetic1);
        
        fromBeginning(*f1, kinetic2, kinetic1);
        f1->tulip[kinetic2].Partition = allowQ(f1->rt,blockKineticEnergyBlock)*COMPONENT;//
        f1->tulip[kinetic2].name = kinetic2;
        assignParticle(*f1, kinetic2, all, one);
        tClear(*f1, kinetic2);

        fromBeginning(*f1, kinetic3, kinetic2);
        f1->tulip[kinetic3].Partition = allowQ(f1->rt,blockKineticEnergyBlock)*COMPONENT;//
        f1->tulip[kinetic3].name = kinetic3;
        assignParticle(*f1, kinetic3, all, one);
        tClear(*f1, kinetic3);

        fromBeginning(*f1, kinetic4, kinetic3);
        f1->tulip[kinetic4].Partition = allowQ(f1->rt,blockKineticEnergyBlock)*COMPONENT;//
        f1->tulip[kinetic4].name = kinetic4;
        assignParticle(*f1, kinetic4, all, one);
        tClear(*f1, kinetic4);
        
        assignOneWithPointers(*f1, kinetic_2, all);
        fromBeginning(*f1, kinetic1_2, kinetic4);
        f1->tulip[kinetic1_2].Partition = allowQ(f1->rt,blockKineticEnergyBlock)*COMPONENT;//
        f1->tulip[kinetic1_2].name = kinetic1;
        assignParticle(*f1, kinetic1_2, all, one);
        tClear(*f1, kinetic1_2);
        
        fromBeginning(*f1, kinetic2_2, kinetic1_2);
        f1->tulip[kinetic2_2].Partition = allowQ(f1->rt,blockKineticEnergyBlock)*COMPONENT;//
        f1->tulip[kinetic2_2].name = kinetic2_2;
        assignParticle(*f1, kinetic2_2, all, one);
        tClear(*f1, kinetic2_2);

        fromBeginning(*f1, kinetic3_2, kinetic2_2);
        f1->tulip[kinetic3_2].Partition = allowQ(f1->rt,blockKineticEnergyBlock)*COMPONENT;//
        f1->tulip[kinetic3_2].name = kinetic3_2;
        assignParticle(*f1, kinetic3_2, all, one);
        tClear(*f1, kinetic3_2);

        fromBeginning(*f1, kinetic4_2, kinetic3_2);
        f1->tulip[kinetic4_2].Partition = allowQ(f1->rt,blockKineticEnergyBlock)*COMPONENT;//
        f1->tulip[kinetic4_2].name = kinetic4_2;
        assignParticle(*f1, kinetic4_2, all, one);
        tClear(*f1, kinetic4_2);

        
        
        fromBeginning(*f1, kineticMass, kinetic4_2);
        f1->tulip[kineticMass].Partition = allowQ(f1->rt,blockKineticEnergyBlock);//
        assignOneWithPointers(*f1, kineticMass,all);
    }
    fromBeginning(*f1, hamiltonian, kineticMass);
        f1->tulip[hamiltonian].Partition = allowQ(f1->rt,blockAllocateHamiltonianBlock)*c1->i.canonRank;//
        f1->tulip[hamiltonian].species = matrix;
        f1->tulip[hamiltonian].spinor = real;

    if ( bootBodies == one ){
        assignParticle(*f1, hamiltonian, all, one);
        //no active matrice-blocks...so its just an assign memory.
        //NEED-TO ALLOW CONTROL ON hamiltonian body number for 7.7
    }else {
        assignParticle(*f1, hamiltonian, electron, two);
    }
    
    
        fromBeginning(*f1, trainHamiltonian, hamiltonian);
        f1->tulip[trainHamiltonian].Partition = allowQ(f1->rt,blockTrainHamiltonianBlock)*c1->i.decomposeRankMatrix;//
        f1->tulip[trainHamiltonian].species = matrix;
    
        if ( bootBodies == one ){
                assignOneWithPointers(*f1, trainHamiltonian, all);
            //NEED-TO ALLOW CONTROL ON trainHamiltonian body number for 7.7

        }
    else{
        assignParticle(*f1, trainHamiltonian, electron, two);
        
        {
            {
                enum block ee ;
                for (ee = e12 ; ee <= e34 ; ee++){
                    f1->tulip[trainHamiltonian+ee-e12+1].name = trainHamiltonian;
                    f1->tulip[trainHamiltonian+ee-e12+1].spinor = spins(*f1,trainHamiltonian);
                    f1->tulip[trainHamiltonian+ee-e12+1].species = matrix;
                    if ( ee < e34 )
                        f1->tulip[trainHamiltonian+ee-e12+1].linkNext = trainHamiltonian+ee-e12+2;
                    for ( space = 0; space < SPACE ; space++)
                        if ( f1->rose[space].body >= two && f1->rose[space].particle == electron )
                            f1->tulip[trainHamiltonian+ee-e12+1].space[space].block = ee;
                }
            }
        }
    }
        fromBeginning(*f1, protonRepulsion, trainHamiltonian);
        f1->tulip[protonRepulsion].Partition = allowQ(f1->rt,blockHamiltonianBlock)*c1->i.oneBody.num;//
        assignOneWithPointers(*f1, protonRepulsion,all);
        
        fromBeginning(*f1, vectorMomentum, protonRepulsion);
        f1->tulip[vectorMomentum].spinor = cmpl;
        assignOneWithPointers(*f1, vectorMomentum,electron);
        f1->tulip[vectorMomentum].Partition = allowQ(f1->rt,blockHamiltonianBlock)*(COMPONENT* c1->i.springFlag+4*c1->i.magFlag);//

        fromBeginning(*f1, harmonium, vectorMomentum);
        assignParticle(*f1, harmonium,electron,one);
        f1->tulip[harmonium].Partition = 0;//
        
        fromBeginning(*f1, X, harmonium);
        assignOneWithPointers (*f1, X,electron);
        f1->tulip[X].Partition =  0 ;//make it a semi-local Gaussian * x
        
        fromBeginning(*f1, linear, X);
      
    if (! f1->chainFlag)
    assignOneWithPointers (*f1, linear,electron);
if (f1->eikonFlag)
    f1->tulip[linear].Partition =  0;
//else
    //    f1->tulip[linear].Partition = buildExternalPotential(c1, f1,1.,nullName,0,1,id0,1,electron ,0,real)+(GAS==0)*COMPONENT  ;//

    fromBeginning(*f1, overlap, linear);
        f1->tulip[overlap].species = matrix;
#ifdef GAUSSIANSINC
        f1->tulip[overlap].Partition = 1;//
#endif
        f1->tulip[overlap].header = bootShape;
    assignOneWithPointers(*f1, overlap, all);

            fromBeginning(*f1, overlapTwo,overlap);
            f1->tulip[overlapTwo].species = matrix;
    #ifdef GAUSSIANSINC
            f1->tulip[overlapTwo].Partition = 1;//
    #endif
            f1->tulip[overlapTwo].header = bootShape;
            assignParticle(*f1, overlapTwo, all, two);

    
    
        {
            fromBeginning(*f1, build, overlapTwo);
            
            {
                        
                INT_TYPE ra = 0;
                       
                       switch ( bootBodies ){
                           case three:
                               ra = 3;
                               break;
                           case four :
                               ra = 6;
                               break;
                       }
                f1->tulip[build].Partition = allowQ(f1->rt, blockfoundationMblock)*ra;
                
            }
            f1->tulip[build].spinor = real;;
            f1->tulip[build].species = matrix;
            assignParticle(*f1, build, all, bootBodies);

            for ( space = 0; space < SPACE ; space++)
            {
                if ( space == 0 )
                    fromBeginning(*f1, bill1+space, build);
                else
                    fromBeginning(*f1, bill1+space, bill1+space-1);
                
                f1->tulip[bill1+space].Partition = allowQ(f1->rt,blockFoundationBlock)*vectorLen(*f1, space)*vectorLen(*f1, space) ;
                INT_TYPE p = f1->tulip[bill1+space].Partition;
              //  f1->tulip[bill1+space].spinor = real;
                f1->tulip[bill1+space].memory = bufferAllocation;
                
            }
            fromBeginning(*f1, eigen, bill1+SPACE-1);
            f1->tulip[eigen].Partition = allowQ(f1->rt, blockfoundationMblock);
            f1->tulip[eigen].species = matrix;
            f1->tulip[eigen].spinor = f1->cmpl;
            assignParticle(*f1, eigen, all, bootBodies);
            
            {
                INT_TYPE di,cmpl;
                enum division last = eigen;
            //    INT_TYPE booting[7];
            //    for ( di = 0; di < 7 ; di++)
            //        booting[di] = 0;
                
                
                if ( f->i.filesVectorOperator )
                    
                {
                    enum bodyType bd;
                    INT_TYPE fi,lines = 0,num;
                    size_t ms = MAXSTRING;
                    char line0[MAXSTRING];
                    char name[MAXSTRING];
                    char title[MAXSTRING];
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
                                    fromBeginning(*f1, vectorOperator+lines, last);
                                    INT_TYPE part1 = 0;
                                    enum genus genus ;
                                    for ( cmpl = f1->cmpl-1 ; cmpl >= 0 ; cmpl--){
                                        tFromReadToFilename(NULL, line,  name, f1->cmpl-1,cmpl,title,&num);
                                        part1 = imax(part1,inputFormat(*f1, name, nullName, 2));
                                    }//name = real component here.
                                    f1->tulip[vectorOperator+lines].Partition = part1;
                                    genus = inputFormat(*f1, name, nullName, 0);
                                    bd = inputFormat(*f1, name, nullName, 100+0/COMPONENT);
                                    printf("%d %d %s\n", bd, genus,name);
                                    if ( (bd==two && genus == vector )|| (bd == one && genus == matrix )){
                                        printf("matrix-%dbody-Op\n",1);
                                        f1->tulip[vectorOperator+lines].species = matrix;
                                        for ( space = 0 ; space < 3/*particular*/ ; space++)
                                            f1->tulip[vectorOperator+lines].space[space].body = one;
                                    }else {
                                        printf("vector*vector-Op\n");

                                        f1->tulip[vectorOperator+lines].species = outerVector;
                                        for (space = 0; space < SPACE ; space++){
                                            f1->tulip[vectorOperator+lines].space[space].body = bd;
                                         //   booting[ bd - bootBodies ] = 1;
                                        }
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
                    f1->tulip[eigenVectors+di].spinor =f1->cmpl;
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
        f1->tulip[edgeElectronMatrix].spinor = real;
        assignOneWithPointers(*f1, edgeElectronMatrix, electron);
        
        fromBeginning(*f1,edgeHamiltonian, edgeElectronMatrix );
        f1->tulip[edgeHamiltonian].Partition = COMPONENT ;//
        f1->tulip[edgeHamiltonian].species = matrix;
        f1->tulip[edgeHamiltonian].spinor = real;
        assignOneWithPointers(*f1, edgeHamiltonian, electron);

        fromBeginning(*f1,edgeProtonMatrix,edgeHamiltonian);
        f1->tulip[edgeProtonMatrix].Partition = 1 ;//
        f1->tulip[edgeProtonMatrix].species = matrix;
        f1->tulip[edgeProtonMatrix].spinor = real;
        assignOneWithPointers(*f1, edgeProtonMatrix, proton);

        fromBeginning(*f1,productVector,edgeProtonMatrix);
        f1->tulip[productVector].Partition =  1;
        f1->tulip[productVector].species = vector;
        f1->tulip[productVector].spinor = parallel;

    
        fromBeginning(*f1,scalarTemp,productVector);
        f1->tulip[scalarTemp].Partition =  ( f1->cmpl == 2 )*  maxVector*2;
        f1->tulip[scalarTemp].species = vector;
        f1->tulip[scalarTemp].spinor = real;

    
        fromBeginning(*f1,permutationVector,scalarTemp);
        f1->tulip[permutationVector].Partition =  maxVector;
        f1->tulip[permutationVector].species = vector;
        f1->tulip[permutationVector].spinor = parallel;
        
        fromBeginning(*f1,permutation2Vector,permutationVector);
        f1->tulip[permutation2Vector].Partition = maxVector;
        f1->tulip[permutation2Vector].species = vector;
        f1->tulip[permutation2Vector].spinor = parallel;
    
    fromBeginning(*f1,northoKet,permutation2Vector);
#ifdef GAUSSIANSINC
    f1->tulip[northoKet].Partition = 1;
#endif
    f1->tulip[northoKet].species = vector;
    f1->tulip[northoKet].spinor = parallel;

    
    
        fromBeginning(*f1,canonicalmvVector,northoKet);
        f1->tulip[canonicalmvVector].Partition = 1;
        f1->tulip[canonicalmvVector].species = vector;
        f1->tulip[canonicalmvVector].spinor = parallel;

        fromBeginning(*f1,canonicalmv2Vector,canonicalmvVector);
        f1->tulip[canonicalmv2Vector].Partition = 1;
        f1->tulip[canonicalmv2Vector].species = vector;
        f1->tulip[canonicalmv2Vector].spinor = parallel;

        fromBeginning(*f1,canonicalmv3Vector,canonicalmv2Vector);
        f1->tulip[canonicalmv3Vector].Partition = 2;
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
        f1->tulip[canonicalvvVector].Partition = 1*0;
        f1->tulip[canonicalvvVector].species = vector;;
        f1->tulip[canonicalvvVector].spinor = parallel;
        
        fromBeginning(*f1,canonicalvv2Vector,canonicalvvVector);
        f1->tulip[canonicalvv2Vector].Partition = 1*0;
        f1->tulip[canonicalvv2Vector].species = vector;;
        f1->tulip[canonicalvv2Vector].spinor = parallel;
        
        fromBeginning(*f1,canonicalvv3Vector,canonicalvv2Vector);
        f1->tulip[canonicalvv3Vector].Partition = 1*0;
        f1->tulip[canonicalvv3Vector].species = vector;;
        f1->tulip[canonicalvv3Vector].spinor = parallel;

        fromBeginning(*f1,canonicalmeVector,canonicalvv3Vector);
        f1->tulip[canonicalmeVector].Partition = 1;
        f1->tulip[canonicalmeVector].species = vector;;
        f1->tulip[canonicalmeVector].spinor = parallel;
        
        fromBeginning(*f1,canonicalme2Vector,canonicalmeVector);
#ifdef GAUSSIANSINC
        f1->tulip[canonicalme2Vector].Partition = 1;
#endif
    f1->tulip[canonicalme2Vector].species = vector;;
        f1->tulip[canonicalme2Vector].spinor = parallel;
        
        fromBeginning(*f1,canonicalme3Vector,canonicalme2Vector);
#ifdef GAUSSIANSINC
        f1->tulip[canonicalme3Vector].Partition = 1;
#endif
        f1->tulip[canonicalme3Vector].species = vector;;
        f1->tulip[canonicalme3Vector].spinor = parallel;

        fromBeginning(*f1,copyVector,canonicalme3Vector);
        f1->tulip[copyVector].Partition = maxVector;
     //   if ( ! f1->cat )
       //     f1->tulip[copyVector].Partition *= ra;
        f1->tulip[copyVector].species = vector;
        //f1->tulip[copyVector].spinor = parallel;

        fromBeginning(*f1,copyTwoVector,copyVector);
        f1->tulip[copyTwoVector].Partition =   0;//f->i.bRank*f->i.bRank;
       // if ( ! f1->cat )
        //    f1->tulip[copyTwoVector].Partition *= ra;
        f1->tulip[copyTwoVector].species = vector;
        //f1->tulip[copyTwoVector].spinor = parallel;

        fromBeginning(*f1,copyThreeVector,copyTwoVector);
        f1->tulip[copyThreeVector].Partition =   maxVector*0;//f->i.bRank*f->i.bRank;
        f1->tulip[copyThreeVector].species = vector;
        f1->tulip[copyThreeVector].spinor = parallel;

        fromBeginning(*f1,copyFourVector,copyThreeVector);
        f1->tulip[copyFourVector].Partition =    maxVector*0;//f->i.bRank*f->i.bRank;
        f1->tulip[copyFourVector].species = vector;
        f1->tulip[copyFourVector].spinor = parallel;
        
        fromBeginning(*f1,oneVector,copyFourVector);
if (f1->eikonFlag)
        f1->tulip[oneVector].Partition = 0;
else
    f1->tulip[oneVector].Partition = 0;//buildExternalPotential(c1, f1,1.,0,1,id0,1,nullName,electron ,0,real);
        f1->tulip[oneVector].species = outerVector;
       // f1->tulip[oneVector].spinor = parallel;
        assignParticle(*f1, oneVector, all, one);

        fromBeginning(*f1,twoVector,oneVector);
        f1->tulip[twoVector].Partition = 0 * f->i.bRank;;
        f1->tulip[twoVector].species = outVector;
        f1->tulip[twoVector].spinor = parallel;
        assignParticle(*f1, oneVector, all, two);

        fromBeginning(*f1,totalVector,twoVector);
  
    if ( allowQ(f1->rt,blockTotalVectorBlock) ){
            f1->tulip[totalVector].Partition =  imax( f->i.xRank, c1->i.canonRank * maxVector) ;
        
    }
    if ( c1->rt.phaseType == productKrylov   )
        f1->tulip[totalVector].species = vector;
    else
    if ( c1->rt.phaseType == singleKrylov ){//not implemented to effect yet...
        f1->tulip[totalVector].species = outerVector;
        assignParticle(*f1, totalVector, all,f->i.bodyTo);
    }
    if (  f1->rt->phaseType == productKrylov )
        f1->tulip[totalVector].spinor = parallel;

    
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
        f1->tulip[bandBasis].Partition = 4*mx1len*mx1len+2*mxlen;
        f1->tulip[bandBasis].memory = bufferAllocation;

    

    
        fromBeginning(*f1,copy,bandBasis);
        f1->tulip[copy].Partition = 1+(c1->rt.phaseType == buildFoundation)*(c1->rt.phaseType == svdOperation)* imax(maxVector*f->i.bRank,imax(mx1len*mx1len, c1->i.decomposeRankMatrix+ imax(c1->i.oneBody.num,   imax(c1->i.decomposeRankMatrix*c1->i.Na,outVector) )));
        if ( c1->rt.phaseType == svdOperation){
                f1->tulip[copy].spinor = parallel;
        }
        f1->tulip[copy].species = matrix;
        assignParticle(*f1, copy, all, one);
        for ( space = 0; space < SPACE ; space++)
            if ( f1->rose[space].body != nada)
                f1->tulip[copy].space[space].block = tv1;
        
        fromBeginning(*f1,copyTwo,copy);
        f1->tulip[copyTwo].Partition = (c1->rt.phaseType == buildFoundation)*c1->i.decomposeRankMatrix;
    if (  c1->rt.phaseType == svdOperation ){
        f1->tulip[copyTwo].Partition = 10*2*c1->i.twoBody.num * f->i.qFloor * f->i.iRank *f->i.iRank;
    }
        f1->tulip[copyTwo].species = matrix;
        assignParticle(*f1, copyTwo, all, one);
        
        fromBeginning(*f1,copyThree,copyTwo);
        f1->tulip[copyThree].Partition = 0*c1->i.decomposeRankMatrix+ 0* imax(c1->i.oneBody.num,   imax(c1->i.decomposeRankMatrix*c1->i.Na,outVector) );
        f1->tulip[copyThree].spinor = parallel;
        f1->tulip[copyThree].species = matrix;
        assignParticle(*f1, copyThree, all, one);
        
        fromBeginning(*f1,squareTwo,copyThree);
        f1->tulip[squareTwo].Partition= 0*(c1->rt.phaseType == buildFoundation);//BUILD
        f1->tulip[squareTwo].species = matrix;
        assignParticle(*f1, squareTwo, all, two);
    
    fromBeginning(*f1,tempOneMatrix,squareTwo);
    f1->tulip[tempOneMatrix].Partition= allowQ(f1->rt,blockBuildHamiltonianBlock)*(c1->rt.calcType == clampProtonElectronCalculation )*c1->i.twoBody.num*vector1Len(*f1, 0)*vector1Len(*f1, 0);//BUILD
    f1->tulip[tempOneMatrix].species = matrix;
    assignParticle(*f1, tempOneMatrix, all, one);
//
    
    //for new interaction builder
    //no need to double train

//    fromBeginning(*f1,tempTwoMatrix,tempOneMatrix);
//    f1->tulip[tempTwoMatrix].Partition= allowQ(f1->rt,blockBuildHamiltonianBlock)*buildPairWisePotential(c1, *f1,1.,nullName, electron, 0,real);//BUILD
//    f1->tulip[tempTwoMatrix].species = matrix;
//    f1->tulip[tempTwoMatrix].header = Cube;
//    assignParticle(*f1, tempTwoMatrix, all, two);

    
        fromBeginning(*f1,inversion,tempOneMatrix);
#ifdef GAUSSIANSINC
        f1->tulip[inversion].Partition= 1;
#endif
        f1->tulip[inversion].species = matrix;
    assignOneWithPointers(*f1, inversion, all);
    
    
            fromBeginning(*f1,inversionTwo,inversion);
    #ifdef GAUSSIANSINC
            f1->tulip[inversionTwo].Partition= 1;
    #endif
            f1->tulip[inversionTwo].species = matrix;
            assignParticle(*f1, inversionTwo, all, two);

    
    {
        fromBeginning(*f1,bufferChromatic,inversionTwo);
        
#ifdef CHROME
        //canonical rank of each training piece.
        f1->chroma = c1->i.chroma;
        //threshold initial
        f1->chromos = c1->i.chromos;
        f1->chromous = c1->i.chromous;
        f1->chromaticStep  = 1.03;
#endif

        f1->tulip[bufferChromatic].spinor = parallel;
#ifdef CHROME
        f1->tulip[bufferChromatic].Partition = c1->i.chromaticRank;
#endif

        if (f1->rt->phaseType == reportMatrix|| f1->rt->phaseType == distillMatrix  || f1->rt->phaseType == decomposeMatrix){
            f1->tulip[bufferChromatic].species = matrix;
            if ( bootBodies >= two ){
                assignParticle(*f1, bufferChromatic, all, two);
            }
            else{
                assignParticle(*f1, bufferChromatic, all, one);
            }
        }else {//train vectors
            f1->tulip[bufferChromatic].species = vector;
        }

    
    }
    
    {
        INT_TYPE maxOriginRank = imax( NV*allowQ(f1->rt, blockEigenDecomposeBlock) , part(*f1,totalVector));
        if ( c1->rt.phaseType == buildFoundation)
            maxOriginRank =c1->i.decomposeRankMatrix;
        INT_TYPE maxTrainRank;
        if (c1->rt.phaseType == distillMatrix || c1->rt.phaseType == reportMatrix || c1->rt.phaseType == decomposeMatrix)//could decrease!!!
            maxTrainRank = c1->i.decomposeRankMatrix;
        else
            maxTrainRank = maxVector;
        INT_TYPE flag = 1;
       // printf("maxOrigin %d\nmaxTrain %d\n", maxOriginRank,maxTrainRank);
        fromBeginning(*f1,canonicalBuffers,bufferChromatic);
        f1->tulip[canonicalBuffers].Partition = flag*(maxTrainRank*maxTrainRank+ maxOriginRank*maxTrainRank+maxTrainRank) ;
        f1->tulip[canonicalBuffers].spinor = parallel;
        f1->tulip[canonicalBuffers].species = scalar;

        fromBeginning(*f1,trackBuffer,canonicalBuffers);
        f1->tulip[trackBuffer].Partition = 2*flag*(maxTrainRank*maxTrainRank);
        f1->tulip[trackBuffer].spinor = parallel;
        f1->tulip[trackBuffer].memory = bufferAllocation;
        f1->tulip[trackBuffer].species = scalar;

        fromBeginning(*f1,guideBuffer,trackBuffer);
        f1->tulip[guideBuffer].Partition = flag*(maxOriginRank*maxTrainRank);
        f1->tulip[guideBuffer].spinor = parallel;
        f1->tulip[guideBuffer].memory = bufferAllocation;
        f1->tulip[guideBuffer].species = scalar;

        fromBeginning(*f1,canonicalBuffersB,guideBuffer);
        f1->tulip[canonicalBuffersB].Partition = (allowQ(f1->rt, blockTrainMatricesblock)||allowQ(f1->rt, blockTrainVectorsblock))* (maxTrainRank+maxOriginRank);
        f1->tulip[canonicalBuffersB].spinor = parallel;
        f1->tulip[canonicalBuffersB].memory = bufferAllocation;
        f1->tulip[canonicalBuffersB].species = scalar;

        
#ifdef BUFFERSOLVE
        //no need for TOTAL ALS SOLVES
        flag = 0;
        
#else
        //need total ALS SOLVES
        flag =0  ;
#endif
        if ( c1->rt.calcType == clampProtonElectronCalculation && c1->rt.phaseType == reportMatrix ){
            maxOriginRank = imax( maxOriginRank, c1->i.twoBody.num * N1*N1);
        }
        else{
            if ( c1->rt.phaseType == distillMatrix || c1->rt.phaseType == decomposeMatrix)
                maxOriginRank = part(*f1,hamiltonian);
            else
                maxOriginRank = part(*f1,totalVector);
        }
        
        fromBeginning(*f1,canonicalBuffers0,canonicalBuffersB);
        f1->tulip[canonicalBuffers0].Partition = flag * maxTrainRank*maxTrainRank+ maxOriginRank*maxTrainRank ;
        
        fromBeginning(*f1,trackBuffer0,canonicalBuffers0);
        f1->tulip[trackBuffer0].Partition =  flag * (2*maxTrainRank+1)*maxTrainRank;
        f1->tulip[trackBuffer0].spinor = real;
        f1->tulip[trackBuffer0].memory = bufferAllocation;
        
        fromBeginning(*f1,guideBuffer0,trackBuffer0);
        f1->tulip[guideBuffer0].Partition = flag * maxOriginRank*maxTrainRank;
        f1->tulip[guideBuffer0].spinor = real;
        f1->tulip[guideBuffer0].memory = bufferAllocation;


        
    }
        fromBeginning(*f1,foundationStructure,guideBuffer0);
        f1->tulip[foundationStructure].spinor = parallel;
        f1->tulip[foundationStructure].Partition = 1;
        f1->tulip[foundationStructure].species = vector;
        f1->tulip[foundationStructure].spinor = cmpl;//need two channels
        
    
    
        fromBeginning(*f1,interactionExchange,foundationStructure);
if (f1->eikonFlag)
        f1->tulip[interactionExchange].Partition = 0;
    else
        f1->tulip[interactionExchange].Partition = allowQ(f1->rt, blockSeparateTwoBodyBlock)*(  allowQ(f1->rt,blockHamiltonianBlock))*50+allowQ(f1->rt, blockResolveTwoBodyBlock)*part(*f1,trainHamiltonian);
            f1->tulip[interactionExchange].species = matrix;
            f1->tulip[interactionExchange].spinor = f1->cmpl;
            assignParticle(*f1, interactionExchange, electron, two);
    
    
    {
        {
            enum block ee ;
            for (ee = e12 ; ee <= e34 ; ee++){
                f1->tulip[interaction12+ee-e12].name = interactionExchange;
                f1->tulip[interaction12+ee-e12].spinor = spins(*f1,interactionExchange);
                f1->tulip[interaction12+ee-e12].species = matrix;
                if ( ee < e34 ){
                    if (! f1->chainFlag)
                        f1->tulip[interaction12+ee-e12].linkNext = interaction12+ee-e12+1;
                }
                for ( space = 0; space < SPACE ; space++)
                    if ( f1->rose[space].body >= two && f1->rose[space].particle == electron )
                        f1->tulip[interaction12+ee-e12].space[space].block = ee;
            }
        }
    }
    
    fromBeginning(*f1,interactionDirect,interactionExchange);
    //7.7
    f1->tulip[interactionDirect].Partition = allowQ(f1->rt, blockResolveTwoBodyBlock)*part(*f1, interactionExchange);
    f1->tulip[interactionDirect].species = matrix;
    f1->tulip[interactionDirect].spinor = f1->cmpl;
    assignParticle(*f1, interactionDirect, electron, two);
    
    {
        {
            enum block ee ;
            for (ee = e12 ; ee <= e34 ; ee++){
                f1->tulip[interactionD12+ee-e12].name = interactionDirect;
                f1->tulip[interactionD12+ee-e12].spinor = spins(*f1,interactionDirect);
                f1->tulip[interactionD12+ee-e12].species = matrix;
                if ( ee < e34 )
                    f1->tulip[interactionD12+ee-e12].linkNext = interactionD12+ee-e12+1;
                for ( space = 0; space < SPACE ; space++)
                    if ( f1->rose[space].body >= two && f1->rose[space].particle == electron )
                        f1->tulip[interactionD12+ee-e12].space[space].block = ee;
            }
        }
    }
    
    
        fromBeginning(*f1,interactionExchangeB,interactionDirect);
        f1->tulip[interactionExchangeB].Partition =  allowQ(f1->rt,blockHamiltonianBlock)*c1->i.twoBody.num*( bootBodies > one )* ( c1->rt.calcType == protonsElectronsCalculation);
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
if (f1->eikonFlag)
    f1->tulip[interactionEwald].Partition =  0;
else
    f1->tulip[interactionEwald].Partition =  allowQ(f1->rt, blockSeparateTwoBodyBlock)* allowQ(f1->rt,blockHamiltonianBlock)*50 * (c1->rt.runFlag > 0 );

        f1->tulip[interactionEwald].species = matrix;
        f1->tulip[interactionEwald].spinor = f1->cmpl;
        assignParticle(*f1, interactionEwald, electron, two);

    {
        enum block ee ;
        for (ee = e12 ; ee <= e34 ; ee++){
            f1->tulip[interaction12Ewald+ee-e12].name = interactionEwald;
            f1->tulip[interaction12Ewald+ee-e12].species = matrix;
            f1->tulip[interaction12Ewald+ee-e12].spinor = spins(*f1,interactionEwald );
            if ( ee < e34 )
                f1->tulip[interaction12Ewald+ee-e12].linkNext = interaction12Ewald+ee-e12+1;
            for ( space = 0; space < SPACE ; space++)
                if ( f1->rose[space].body >= two && f1->rose[space].particle == electron )
                    f1->tulip[interaction12Ewald+ee-e12].space[space].block = ee;
        }
    }
    
        fromBeginning(*f1,intracellularSelfEwald,interactionEwald);
#ifndef EIKON
        f1->tulip[intracellularSelfEwald].Partition =  allowQ(f1->rt,blockHamiltonianBlock)*buildPairWisePotential(c1, f1,1.,nullName, electron, 0,real)* (c1->rt.runFlag > 0 );
#endif
        f1->tulip[intracellularSelfEwald].species = matrix;
        assignOneWithPointers(*f1, intracellularSelfEwald, electron);
    
        fromBeginning(*f1,intercellularSelfEwald,intracellularSelfEwald);
    #ifndef EIKON

        f1->tulip[intercellularSelfEwald].Partition =  allowQ(f1->rt,blockHamiltonianBlock)* buildPairWisePotential(c1, f1,1.,nullName, electron, 0,real)* (c1->rt.runFlag > 0 );
#endif
        f1->tulip[intercellularSelfEwald].species = matrix;
        assignOneWithPointers(*f1, intercellularSelfEwald, electron);
    
        fromBeginning(*f1,jelliumElectron,intercellularSelfEwald);
    #ifndef EIKON

        f1->tulip[jelliumElectron].Partition =  allowQ(f1->rt,blockHamiltonianBlock)*( buildPairWisePotential(c1, f1,1.,nullName, electron, 0,real)+1)* (c1->rt.runFlag > 0 );
#endif
        f1->tulip[jelliumElectron].species = matrix;
        assignOneWithPointers(*f1, jelliumElectron, electron);

        fromBeginning(*f1,shortenPlus,jelliumElectron);
        f1->tulip[shortenPlus].Partition = allowQ(f1->rt,blockHamiltonianBlock)*c1->i.decomposeRankMatrix*( c1->rt.calcType == clampProtonElectronCalculation )/2;
        f1->tulip[shortenPlus].species = matrix;
        assignOneWithPointers(*f1, shortenPlus, all);
        
        fromBeginning(*f1,shortenMinus,shortenPlus);
        f1->tulip[shortenMinus].Partition = allowQ(f1->rt,blockHamiltonianBlock)*c1->i.decomposeRankMatrix*( c1->rt.calcType == clampProtonElectronCalculation )/2;
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
        f1->tulip[shortTwoAcrossDimensions].Partition = 0*c1->i.twoBody.num*c1->i.decomposeRankMatrix*c1->i.twoBody.num*( c1->rt.calcType == protonsElectronsCalculation );
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
        f1->tulip[quadCube].Partition = allowQ(f1->rt,blockSeparateTwoBodyBlock)||allowQ(f1->rt,blockBuildHamiltonianBlock) || allowQ(f1->rt,blockHamiltonianBlock)||allowQ(f1->rt,blockAllocateHamiltonianBlock) || allowQ(f1->rt,blockfoundationMblock);
        f1->tulip[quadCube].species = matrix;
        f1->tulip[quadCube].spinor = real;
        assignParticle(*f1, quadCube, all, two);
   
    
    
    
        //EIKON
        //EIKON
    fromBeginning(*f1,eikonBuffer,quadCube);
    f1->tulip[eikonBuffer].Partition = 1;//beta independent
    f1->tulip[eikonBuffer].species = eikon;
    if ( bootBodies > one )
        assignParticle(*f1, eikonBuffer, all, two);
    else
        assignParticle(*f1, eikonBuffer, all, one);
        
        
    
    
    {
        INT_TYPE ii;
        enum division prev = eikonBuffer;
        for ( ii = 0 ; ii < f1->eikonLabels.maxLabel ; ii++)
            {
                enum division label1 = f1->eikonLabels.head+ii;
                f1->tulip[label1].Partition = 1;
                if ( bootBodies > one )
                    assignParticle(*f1, label1, all, two);
                else
                    assignParticle(*f1, label1, all, one);
                f1->tulip[label1].species = eikon;
                fromBeginning(*f1,label1,prev);
                prev = label1;
            }        
        //EIKON
        //EIKON

        
        fromBeginning(*f1,oneArray,prev);
    }
        f1->tulip[oneArray].Partition = c1->i.M1*3+c1->i.M1*c1->i.M1;
        f1->tulip[oneArray].memory = bufferAllocation;

        fromBeginning(*f1,threeArray,oneArray);
        f1->tulip[threeArray].Partition = 2*c1->i.M1*c1->i.M1*c1->i.M1;
        f1->tulip[threeArray].memory = bufferAllocation;

        fromBeginning(*f1,oneBasis,threeArray);
        f1->tulip[oneBasis].Partition = c1->i.M1*N1;
        f1->tulip[oneBasis].memory = bufferAllocation;

        INT_TYPE vecLen = 1;
        for ( space = 0  ; space < SPACE ; space++)
            if ( vecLen < alloc(*f1, eigenVectors, space))
                vecLen = alloc(*f1, eigenVectors, space);
        
        
        fromBeginning(*f1,tensorBuffers,oneBasis);
        f1->tulip[tensorBuffers].Partition = vecLen*0 ;
        f1->tulip[tensorBuffers].spinor = parallel;
        f1->tulip[tensorBuffers].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers2,tensorBuffers);
        f1->tulip[tensorBuffers2].Partition = vecLen*0;
        f1->tulip[tensorBuffers2].spinor = parallel;
        f1->tulip[tensorBuffers2].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers3,tensorBuffers2);
        f1->tulip[tensorBuffers3].Partition = vecLen*0;
        f1->tulip[tensorBuffers3].spinor = parallel;
        f1->tulip[tensorBuffers3].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers4,tensorBuffers3);
        f1->tulip[tensorBuffers4].Partition = vecLen*0;
        f1->tulip[tensorBuffers4].spinor = parallel;
        f1->tulip[tensorBuffers4].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers5,tensorBuffers4);
        f1->tulip[tensorBuffers5].Partition = vecLen*0;
        f1->tulip[tensorBuffers5].spinor = parallel;
        f1->tulip[tensorBuffers5].memory = bufferAllocation;

        fromBeginning(*f1,tensorBuffers6,tensorBuffers5);
        f1->tulip[tensorBuffers6].Partition = vecLen*0;
        f1->tulip[tensorBuffers6].spinor = parallel;
        f1->tulip[tensorBuffers6].memory = bufferAllocation;
        
        fromBeginning(*f1,oneByOneBuffer,tensorBuffers6);
        f1->tulip[oneByOneBuffer].Partition = mx1len*mx1len*mx1len*mx1len* (c1->rt.calcType >= clampProtonElectronCalculation)*allowQ(f1->rt,blockBuildHamiltonianBlock);
        f1->tulip[oneByOneBuffer].spinor = real;
        f1->tulip[oneByOneBuffer].memory = bufferAllocation;
    
        fromBeginning(*f1,canonicalBuffersC,oneByOneBuffer);
        f1->tulip[canonicalBuffersC].Partition = allowQ(f1->rt,blockEigenDecomposeBlock) *   NV;
        f1->tulip[canonicalBuffersC].spinor = parallel;
        f1->tulip[canonicalBuffersC].memory = bufferAllocation;

        fromBeginning(*f1,twoBodyRitz,canonicalBuffersC);
        f1->tulip[twoBodyRitz].Partition = maxArray;
        f1->tulip[twoBodyRitz].spinor = real;
        f1->tulip[twoBodyRitz].memory = bufferAllocation;

        fromBeginning(*f1,conditionOverlapNumbers,twoBodyRitz);
        f1->tulip[conditionOverlapNumbers].Partition = maxArray*0;
        f1->tulip[conditionOverlapNumbers].memory = bufferAllocation;
        
        fromBeginning(*f1,matrixHbuild,conditionOverlapNumbers);
    f1->tulip[matrixHbuild].Partition = (  c1->rt.phaseType == solveRitz|| c1->rt.phaseType == svdOperation ) *  (maxArray*maxArray)+( c1->rt.phaseType == buildFoundation)* mxlen*mxlen;
    f1->tulip[matrixHbuild].Partition *= 4;
        f1->tulip[matrixHbuild].spinor = real;
        f1->tulip[matrixHbuild].memory = bufferAllocation;

        fromBeginning(*f1,vectorHbuild,matrixHbuild);
        f1->tulip[vectorHbuild].Partition = 0;
        f1->tulip[vectorHbuild].memory = bufferAllocation;

        fromBeginning(*f1,matrixSbuild,vectorHbuild);
        f1->tulip[matrixSbuild].Partition = ( c1->rt.phaseType == buildFoundation ||  c1->rt.phaseType == solveRitz|| c1->rt.phaseType == svdOperation )*(maxArray*maxArray)+( c1->rt.phaseType == buildFoundation)* mxlen*mxlen;
    f1->tulip[matrixSbuild].Partition *= 4;

        f1->tulip[matrixSbuild].spinor = real;
        f1->tulip[matrixSbuild].memory = bufferAllocation;

    fromBeginning(*f1,square,matrixSbuild);
    f1->tulip[square].Partition=   0*4*maxVector*maxVector;
    f1->tulip[square].species = matrix;
    assignParticle(*f1, square, all, one);
    
    
    
        fromBeginning(*f1,dsyBuffers,square);
        f1->tulip[dsyBuffers].Partition = 1000+2*8*(8*(imax(mxlen,maxEV))+72*f->i.nStates*f->i.nStates+ 8 * mxlen)+3*maxEV;
    //unsure.
//    if ( PARTICLE == 1 )
//        f1->tulip[dsyBuffers].Partition = maxVector*maxVector;
        f1->tulip[dsyBuffers].spinor = parallel;
        f1->tulip[dsyBuffers].memory = bufferAllocation;

#ifdef PURITY
 INT_TYPE RS=  definePurity(f1, 3, 12, dsyBuffers);
    fromBeginning(*f1,end,f1->purity+RS);//

#else
    fromBeginning(*f1,f1->end,dsyBuffers);//
#endif
   
        fromBeginning(*f1,f1->end,f1->end);
        struct name_label e = f1->tulip[f1->end];
        struct name_label k1 = f1->tulip[kinetic1];
        struct name_label k2 = f1->tulip[kinetic2];
        struct name_label it = f1->tulip[interaction12];
        struct name_label i = f1->tulip[interactionExchangePlus];
        struct name_label ip1 = f1->tulip[interaction1Plus];

    
        {
            double maxMem = 0.,currMem;
            for ( space = 0 ; space <= SPACE ; space++){
                currMem = (f1->tulip[f1->end].space[space].Address)/(1000000000./(sizeof(Stream_Type)));
             
                if ( f1->boot == fullMatrices ){
//                    printf("\t| SPACE \t:   Gb\t \n");
                    if ( space < SPACE )
                        printf("\t| %d \t\t: \t%1.9f\n",space,currMem);
                    else
                        printf("\t| my \t\t: \t%1.9f\n",currMem);
                    fflush(stdout);
                }
                if ( currMem >= 0 )
                    maxMem += currMem;
                else {
                    exit(0);
                }
            }
            if ( maxMem > c1->i.RAMmax ){
                printf("oops too much RAM required\n");
                fflush(stdout);

                exit(0);
            }
          //  printf("\n\n\n");
        }
    
        for ( space = 0; space <= SPACE ; space++){
            f1->rose[space].stream = malloc( (f1->tulip[f1->end].space[space].Address)*sizeof(Stream_Type));
        }

        f1->bootedMemory = 1;
        assignCores(*f1, 1);

        INT_TYPE RdsSize;
        RdsSize = 0;

        
        if (c1->i.M1){
            INT_TYPE i,ii;
            INT_TYPE M1 = c1->i.M1,M12 = (M1-1)/2;
            double r = (double)(M1-1)/(double)(N1-1);
            double * u =myStreams(*f1,oneBasis,0);
            for( i = 0; i < N1 ; i++)
                for ( ii = 0 ; ii < M1 ; ii++){
                    u[i*M1+ii] = Sinc(r, r*(i-N12)- (ii-M12))/sqrt(r);
                }
        }
    
    
#ifdef GAUSSIANSINC
    separateOverlap(*f1, overlap,overlapTwo, inversion,inversionTwo);
#endif
    
    if ( allowQ(f1->rt,blockKineticEnergyBlock) ) {
        separateKinetic(*f1, 0,kinetic1, c1->i.massElectron,all);
        separateKinetic(*f1, 0,kinetic2, c1->i.massElectron,all);
        separateKinetic(*f1, 0,kinetic3, c1->i.massElectron,all);
//        separateKinetic(*f1, 0,kinetic1_2, c1->i.massElectron,all);
   //     separateKinetic(*f1, 0,kinetic2_2, c1->i.massElectron,all);
       // separateKinetic(*f1, 0,kinetic3_2, c1->i.massElectron,all);
    }
    f1->tulip[Ha].linkNext = Iterator;
    defineTerms(c1, f1,Iterator,f->i.OpIndex);
    
        
#if 0
#if 1
    if (bootBodies == one ){
        f1->tulip[external1].species = matrix;
        f1->tulip[external1].linkNext = nullName;
        f1->tulip[external1].chainNext = nullName;
        f1->tulip[Ha].linkNext = kinetic1;
        f1->tulip[Iterator].linkNext = kinetic1;
        f1->tulip[kinetic1].linkNext = external1;

        buildExternalPotential(c1, f1,tv1,external1,electron,!(!c1->rt.runFlag),real);

    }else if (1){
        switch(2){
            case 0:
        f1->tulip[external1].species = matrix;
        f1->tulip[external2].species = matrix;
                f1->tulip[interaction12].species = matrix;

    f1->tulip[kinetic1].linkNext = kinetic2;
   f1->tulip[kinetic2].linkNext = interaction12;
    f1->tulip[interaction12].linkNext = nullName;

    
    f1->tulip[Ha].linkNext = kinetic1;
    f1->tulip[Iterator].linkNext = kinetic1;

    buildPairWisePotential(c1, f1,e12,interaction12,electron, 0,real);
    buildExternalPotential(c1, f1,tv1,kinetic1,electron,!(!c1->rt.runFlag),real);
    buildExternalPotential(c1, f1,tv2,kinetic2,electron,!(!c1->rt.runFlag),real);
                break;
            case 1:
                     f1->tulip[external1].species = matrix;
                     f1->tulip[external2].species = matrix;

                 f1->tulip[kinetic1].chainNext = kinetic2;
                f1->tulip[kinetic1].linkNext = external1;
                 f1->tulip[external1].linkNext = external2;
                     f1->tulip[external2].linkNext = interaction12;

                 
                 f1->tulip[Ha].linkNext = kinetic1;
                 f1->tulip[Iterator].linkNext = kinetic1;

                 buildPairWisePotential(c1, f1,e12,interaction12,electron, 0,real);
                 buildExternalPotential(c1, f1,tv1,external1,electron,!(!c1->rt.runFlag),real);
                 buildExternalPotential(c1, f1,tv2,external2,electron,!(!c1->rt.runFlag),real);
                break;
            case 2:
                    f1->tulip[external1].species = matrix;

                    
                                    f1->tulip[kinetic1].linkNext = kinetic2;
                                    f1->tulip[kinetic2].linkNext = kinetic3;
                                    f1->tulip[kinetic3].linkNext = interaction12;
                                        
                //HERE
                                    
                                    f1->tulip[Ha].linkNext = kinetic1;
                                    f1->tulip[Iterator].linkNext = kinetic1;
                buildPairWisePotential(c1, f1,e12,interaction12,electron, 0,real);
                buildPairWisePotential(c1, f1,e13,interaction12,electron, 0,real);
                buildPairWisePotential(c1, f1,e23,interaction12,electron, 0,real);

                                        buildExternalPotential(c1, f1,tv1,kinetic1,electron,!(!c1->rt.runFlag),real);
                                        buildExternalPotential(c1, f1,tv2,kinetic2,electron,!(!c1->rt.runFlag),real);
                                        buildExternalPotential(c1, f1,tv3,kinetic3,electron,!(!c1->rt.runFlag),real);

                break;
        }
    }
    
    for( int i = 0 ; i < 1 ;i++) {
       INT_TYPE i;
        tClear ( *f1,copyTwoVector);
        zero(*f1, copyTwoVector, 0);
       for ( i = 0; i < part(*f1,copyTwoVector);i++)
           tId(*f1,copyTwoVector,0);
    //   tId(*f1,copyTwoVector,0);

        enum division spiralOp  =  defSpiralMatrix(f1, Ha);
        tLesserDivide(0, 0., 1., *f1,defSpiralVector(f1, spiralOp, copyVector) ,spiralOp,defSpiralVector(f1, spiralOp, copyTwoVector));
        printExpectationValues(*f1,copyVector, Ha, copyVector);

    }
    printf("\n\nFINIS.\n\n");
    if(0){
        tClear ( *f1,copyVector);
        zero(*f1, copyVector, 0);
        tId(*f1,copyVector,0);

        tHXpX( 0,*f1, Ha, 0, 0., 1., 0., copyVector, 1e-8, part(*f1,copyVector),0);
        printExpectationValues(*f1,copyVector, Ha, copyVector);

    }
    //JUST NEED TO ZERO FROM CURRENT -> PARTION ON EACH subVector Rank structure.
#if 0
    f1->tulip[copyVector].Current[0] = 1;
    streams(*f1,copyVector,0,0)[(N1-1)/2] = 1;
    streams(*f1,copyVector,0,1)[(N1-1)/2] = 1;
    streams(*f1,copyVector,0,2)[(N1-1)/2] = 1;
    streams(*f1,copyVector,0,0)[(N1-1)/2+1] = -1;
    streams(*f1,copyVector,0,1)[(N1-1)/2+1] = -1;
    streams(*f1,copyVector,0,2)[(N1-1)/2+1] = -1;
    streams(*f1,copyVector,0,0)[(N1-1)/2-1] = -1;
    streams(*f1,copyVector,0,1)[(N1-1)/2-1] = -1;
    streams(*f1,copyVector,0,2)[(N1-1)/2-1] = -1;
#endif
    exit(0);
#endif
#endif
//if one has memory to build and wants Hamiltonian...
         if(allowQ(f1->rt,blockBuildHamiltonianBlock)* allowQ(f1->rt,blockHamiltonianBlock))
        {
            if ( c1->rt.calcType == electronicStuctureCalculation  ){
               // separateKinetic(*f1, 0,kinetic, c1->i.massElectron,electron);
                if ( c1->i.Na ){
#ifndef EIKON
                    INT_TYPE flag ;
                    flag =   ioStoreMatrix(*f1,linear ,0,"linear.matrix",1);
                    if ( f1->cmpl == cmpl)
                        flag = flag && ioStoreMatrix(*f1,linear ,1,"linear.1.matrix",1);
#endif
                    
                    for ( c = real ; c <= spins (*f1, linear) ; c++){
                         //  buildExternalPotential(c1, f1,1.,0,1,id0,1,linear,electron,!(!c1->rt.runFlag),c);
                    }
                }
                if ( f1->rose[0].component == periodicComponent1 ){
                    if (c1->i.twoBody.func.fn != nullFunction ){
#ifndef APPLE
                            INT_TYPE flag = ioStoreMatrixScale(f,interactionEwald ,0,"interactionEwald.matrix",1);
                            if ( f1->cmpl == cmpl)
                                flag = flag &&   ioStoreMatrixScale(f,interactionEwald ,1,"interactionEwald.1.matrix",1);
                                if ( ! flag )
#endif
                                    if ( allowQ(f1->rt, blockSeparateTwoBodyBlock ) || f1->eikonFlag)
                                    {
                                        //for ( c = real ; c <= spins (*f1, interactionEwald) ; c++)
                                        //    buildPairWisePotential(c1, f1,1.,0.,1,id0,1,interactionEwald,electron, 1,c);
                            }
                    }
                    if ( c1->i.twoBody.func.fn != nullFunction){
                        //only diagonal!
#ifndef APPLE

                                INT_TYPE flag = ioStoreMatrixScale(f,interactionExchange ,0,"interactionExchange.matrix",1);
                               if ( f1->cmpl == cmpl)
                                   flag = flag &&   ioStoreMatrixScale(f,interactionExchange ,1,"interactionExchange.1.matrix",1);
                     //           if ( ! flag )
#endif
                                //    for ( c = real ; c <= spins (*f1, interactionExchange) ; c++)
                               //         buildPairWisePotential(c1, f1,1.,0,1,id0,1,interactionExchange,electron, 2/*diagonal*/,c);
                        }
                }else {//non-periodic
                    if ( bootBodies > one ){
                        
#ifndef EIKON
                        
                        if ( allowQ(f1->rt, blockSeparateTwoBodyBlock))
                        {
                                    
                                            INT_TYPE flag ;
                                            flag =   ioStoreMatrixScale(f,interactionExchange ,0,"interactionExchange.matrix",1);
                                            if ( f1->cmpl == cmpl)
                                                flag = flag && ioStoreMatrixScale(f,interactionExchange ,1,"interactionExchange.1.matrix",1);
                                            if ( ! flag )
                                                for ( c = real ; c <= spins (*f1, interactionExchange) ; c++){
                                                    buildPairWisePotential(c1, f1,interactionExchange,electron,0,c);
                                                }
                                
                        }
                    
#else
                    for ( c = real ; c <= spins (*f1, interactionExchange) ; c++){
                 //       buildPairWisePotential(c1, f1,1.,0,1,id0,1,interactionExchange,electron,0,c);
                    }

                    
#endif
                    
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
                
                if ( c1->i.springFlag && allowQ(f1->rt, blockKineticEnergyBlock)){
                    INT_TYPE deriv[SPACE];
                    INT_TYPE power[SPACE];
                    
                    
                    deriv[0] = 0;
                    deriv[1] = 0;
                    deriv[2] = 0;
                    
                    power[0] = 2 ;
                    power[1] = 0;
                    power[2] = 0;
                    printf("k %f\n", c1->i.springConstant);
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
                }
            } else if ( c1->rt.calcType == clampProtonElectronCalculation  ){
                
                if ( bootBodies > one ){
                    if ( allowQ(f1->rt, blockSeparateTwoBodyBlock )){
                        if ( c1->i.twoBody.func.fn != nullFunction ){
                            INT_TYPE flag = ioStoreMatrixScale(f,interactionExchange ,0,"interactionExchange.matrix",1);
                           if ( f1->cmpl == cmpl)

                               flag = flag &&   ioStoreMatrixScale(f,interactionExchange ,1,"interactionExchange.1.matrix",1);

                        
                    
                        
                            if ( ! flag )
                                for ( c = real ; c <= spins (*f1, interactionExchange) ; c++){
                                 //   buildPairWisePotential(c1, f1,1.,0,1,id0,1,interactionExchange,electron,0,c);
                                    }
                    }
                }
                }
                
                mySeparateExactOneByOne(*f1,c1->i.twoBody, part(*f1,shortenPlus),interactionExchangePlus, shortenPlus,-1., 1,c1->i.massClampPair/(c1->i.massClampPair+c1->i.massProton),electron, proton);
                
                mySeparateExactOneByOne(*f1,c1->i.twoBody, part(*f1,shortenMinus),interactionExchangeMinus, shortenMinus,-1., -1,c1->i.massProton/(c1->i.massClampPair+c1->i.massProton),electron, pair);
                
                
                
                if ( c1->i.oneBody.func.fn != nullFunction&&c1->rt.phaseType != buildFoundation )
                    separateExternal(c1,*f1,protonRepulsion, 0,0,1.0,-1,0,proton);
               // separateKinetic(*f1, 0,kineticMass,c1->i.massClampPair*c1->i.massProton/(c1->i.massProton + c1->i.massClampPair),proton);
                separateKinetic(*f1, 0,kinetic, c1->i.massElectron*(c1->i.massProton + c1->i.massClampPair)/(c1->i.massElectron+c1->i.massProton + c1->i.massClampPair),electron);
                
            }
        
            //else if one wants the Hamiltonian
        }else     if(allowQ(f1->rt,blockHamiltonianBlock) && f->i.filesVectorOperator == 0 )
        {
            if ( c1->rt.calcType == electronicStuctureCalculation  ){
               // separateKinetic(*f1, 0,kinetic, c1->i.massElectron,electron);

                if ( c1->i.Na )
                    if ( c1->i.oneBody.func.fn != nullFunction )
                        if(   ! ioStoreMatrix(*f1, linear, 0, "linear.matrix",1)){
                            for ( c = real ; c <= real ; c++){
                                //    buildExternalPotential(c1, f1,1.,0,1,id0,1,linear,electron,!(!c1->rt.runFlag),c);
                            }
                        }
                if ( c1->i.springFlag ){
                    ioStoreMatrix( *f1, vectorMomentum, 0 , "vector.matrix", 1 ) ;
                    if ( f1->cmpl == cmpl)
                        ioStoreMatrix( *f1, vectorMomentum, 1 , "vector.1.matrix", 1 ) ;
                }

                if ( f1->rose[0].component == periodicComponent1 ){
                    ioStoreMatrix(*f1, intracellularSelfEwald, 0, "intracellularSelfEwald.matrix",1);
                    ioStoreMatrix(*f1, intercellularSelfEwald, 0,"intercellularSelfEwald.matrix",1);
                    ioStoreMatrix(*f1, jelliumElectron, 0, "jelliumElectron.matrix",1);
                    if ( allowQ(f1->rt, blockSeparateTwoBodyBlock)){
                        ioStoreMatrixScale(f, interactionEwald, 0, "interactionEwald.matrix",1);
                        ioStoreMatrixScale(f, interactionExchange, 0, "interactionExchange.matrix",1);
                    }
                    
                        if ( f1->cmpl == cmpl ){
                            ioStoreMatrix(*f1, intracellularSelfEwald, 1, "intracellularSelfEwald.1.matrix",1);
                            ioStoreMatrix(*f1, intercellularSelfEwald, 1, "intercellularSelfEwald.1.matrix",1);
                            ioStoreMatrix(*f1, jelliumElectron, 1, "jelliumElectron.1.matrix",1);
                    
                            if ( allowQ(f1->rt, blockSeparateTwoBodyBlock)){
                                ioStoreMatrixScale(f, interactionEwald, 1, "interactionEwald.1.matrix",1);
                                ioStoreMatrixScale(f, interactionExchange, 1, "interactionExchange.1.matrix",1);
                            }
                        }
//                    if (!flag ){
//                        printf("ewald terms absent");
//
//                        exit(0);
//                    }
                }else{
                    if ( bootBodies > one )
                        if ( allowQ(f1->rt, blockSeparateTwoBodyBlock))
                        {
                            if ( c1->i.twoBody.func.fn != nullFunction )
                                if(   ! ioStoreMatrixScale(f, interactionExchange, 0, "interactionExchange.matrix",1)){
                                    printf("exchange absent");
                                    
                                    exit(0);
                                }
                        }
                }
                
            } else if ( c1->rt.calcType == clampProtonElectronCalculation  ){
                if ( bootBodies > one )
                    if ( c1->i.twoBody.func.fn != nullFunction )
                        if ( allowQ(f1->rt, blockSeparateTwoBodyBlock))
                        if(  ! ioStoreMatrixScale(f, interactionExchange, 0, "interactionExchange.matrix",1)){
                            printf("exchange absent");

                            exit(0);
                        }
                if ( c1->i.twoBody.func.fn != nullFunction && !ioStoreMatrix(*f1,shortenPlus ,0,"shortenExchangePlus.matrix",1) ){
                    printf("failed to load PLUS\n");
                    exit(0);
                }
                
                if ( c1->i.twoBody.func.fn != nullFunction&& ! ioStoreMatrix(*f1,shortenMinus ,0,"shortenExchangeMinus.matrix",1) ){
                    printf("failed to load MINUS\n");
                    exit(0);
                }
                
                
                if ( c1->i.oneBody.func.fn != nullFunction&&c1->rt.phaseType != buildFoundation )
                    separateExternal(c1,*f1,protonRepulsion, 0,0,1.0,-1,0,proton);
              //  separateKinetic(*f1, 0,kineticMass,c1->i.massClampPair*c1->i.massProton/(c1->i.massProton + c1->i.massClampPair),proton);
                separateKinetic(*f1, 0,kinetic, c1->i.massElectron*(c1->i.massProton + c1->i.massClampPair)/(c1->i.massElectron+c1->i.massProton + c1->i.massClampPair),electron);
                
            }

        }
    
    
        //if One wants to link of Ha as Hamiltonian
    if (allowQ(f1->rt,blockHamiltonianBlock))
        {
          
                if ( bootBodies == one && ( c1->rt.calcType == electronicStuctureCalculation  )){


                    f1->tulip[jellium1Electron].linkNext = intracellularSelf1Ewald;
                    f1->tulip[intracellularSelf1Ewald].linkNext = intercellularSelf1Ewald;
                    f1->tulip[intercellularSelf1Ewald].linkNext = vectorMomentum1;
                    f1->tulip[vectorMomentum1].linkNext = kinetic1;
                    f1->tulip[kinetic1].linkNext = external1;
                    f1->tulip[external1].linkNext = nullName;
                    
                    //active assignment
                    f1->tulip[Ha].linkNext = jellium1Electron;
                    f1->tulip[Iterator].linkNext = jellium1Electron;
                } else if ( bootBodies == two && ( c1->rt.calcType == electronicStuctureCalculation  )){

                    f1->tulip[jellium2Electron].linkNext = intracellularSelf1Ewald;
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
                    f1->tulip[Ha].linkNext = jellium1Electron;
                    f1->tulip[Iterator].linkNext = jellium1Electron;
                    
        
                }else if ( bootBodies == three && ( c1->rt.calcType == electronicStuctureCalculation  )){

                    if ( !f1->chainFlag ){

                    f1->tulip[jellium3Electron].linkNext = intracellularSelf1Ewald;
                    
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
                    f1->tulip[Ha].linkNext = jellium1Electron;
                    f1->tulip[Iterator].linkNext = jellium1Electron;
                    
                    } else {
                    f1->tulip[kinetic3].chainNext = interaction12;
                    f1->tulip[interaction23].chainNext = nullName;
                    f1->tulip[interaction23].linkNext = external1;
                    f1->tulip[external3].linkNext = nullName;
                    
                    f1->tulip[Ha].linkNext = kinetic1;
                    f1->tulip[Iterator].linkNext = kinetic1;

                    }
                }
            
                else if ( bootBodies == four && ( c1->rt.calcType == electronicStuctureCalculation  )){
                    f1->tulip[jellium4Electron].linkNext = intracellularSelf1Ewald;
                    
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
                    f1->tulip[Ha].linkNext = jellium1Electron;
                    f1->tulip[Iterator].linkNext = jellium1Electron;
                    
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
    //if oen wants to link up as trained Hamiltonian
        else if(allowQ(f1->rt,blockTrainHamiltonianBlock)) {
        f1->tulip[Ha].linkNext = h12;
        f1->tulip[Iterator].linkNext = h12;
        if ( c1->i.Na > 0 )
            printf("Trained: Correlated Electron Gas + (External)\n");
        else
            printf("Trained: Correlated Electron Gas\n");

        fflush(stdout);
        ioStoreMatrix(*f1, trainHamiltonian, 0, "trainHamiltonian.matrix",1);
        if ( f->i.cmpl == cmpl )
            ioStoreMatrix(*f1, trainHamiltonian, 1, "trainHamiltonian.1.matrix",1);

        if(allowQ(f1->rt,blockResolveTwoBodyBlock)) {
            printf("exchanging trainHamiltonian...\n");
            f1->tulip[Ha].linkNext = interactionD12;
            f1->tulip[Iterator].linkNext = interactionD12;
            if ( bootBodies == two && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->tulip[interactionD12].linkNext = nullName;
            } else if ( bootBodies == three && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->tulip[interactionD23].linkNext = nullName;
            }
            else if ( bootBodies == four && ( c1->rt.calcType == electronicStuctureCalculation  )){
                f1->tulip[interactionD34].linkNext = nullName;
            }
            exchangeToDirect(*f1, trainHamiltonian, interactionDirect);

        }else {
        
        //really want to keep linear separate and uncarved...
        if ( c1->i.Na && c1->rt.calcType != clampProtonElectronCalculation )
            if ( c1->i.oneBody.func.fn != nullFunction ){
                tClear(*f1, linear);
                if(   ! ioStoreMatrix(*f1, linear, 0, "trainLinear.matrix",1))
                    if ( ! ioStoreMatrix(*f1,linear ,0,"linear.matrix",1)) {
                        printf("failed to load either external or distilled External\n");
                        exit(0);//if not already present, go back and build it
                    }
                if ( f1->cmpl == cmpl )
                    if(   ! ioStoreMatrix(*f1, linear, 1, "trainLinear.1.matrix",1))
                        if ( ! ioStoreMatrix(*f1,linear ,1,"linear.matrix",1)){
                            printf("failed to load either external or distilled External\n");
                            exit(0);//if not already present, go back and build it
                        }

            }
        
        
        
        
        
            switch ( bootBodies ){
            case one:
                f1->tulip[h12].linkNext = external1;
                f1->tulip[external1].linkNext = kinetic1;
                f1->tulip[kinetic1].linkNext  = nullName;
                break;
            case two:
                f1->tulip[h12].linkNext = external1;
                f1->tulip[external2].linkNext = kinetic1;
                f1->tulip[kinetic2].linkNext  = nullName;
                break;
            case three:
                f1->tulip[h23].linkNext = external1;
                f1->tulip[external3].linkNext = kinetic1;
                f1->tulip[kinetic3].linkNext  = nullName;

                break;
            case four:
                f1->tulip[h34].linkNext = external1;
                f1->tulip[external4].linkNext = kinetic1;
                f1->tulip[kinetic4].linkNext  = nullName;

                break;
            }
        }
        }else   if ( f->i.nOperator ){
                      INT_TYPE fi,fe=0,ff;
                      for ( fi =0 ; fi < f->i.filesVectorOperator ; fi++){
                          tLoadEigenWeights (c1,*f,f->i.fileVectorOperator[fi], &fe,f1->vectorOperator,0);
                      }
                      
        f1->tulip[Ha].linkNext = h12;
        f1->tulip[Iterator].linkNext = h12;
        switch ( bootBodies ){
            case one:
                f1->tulip[h12].linkNext = external1;
                f1->tulip[external1].linkNext = kinetic1;
                f1->tulip[kinetic1].linkNext  = nullName;
                break;
            case two:
                f1->tulip[h12].linkNext = external1;
                f1->tulip[external2].linkNext = kinetic1;
                f1->tulip[kinetic2].linkNext  = nullName;
                break;
            case three:
                f1->tulip[h23].linkNext = external1;
                f1->tulip[external3].linkNext = kinetic1;
                f1->tulip[kinetic3].linkNext  = nullName;

                break;
            case four:
                f1->tulip[h34].linkNext = external1;
                f1->tulip[external4].linkNext = kinetic1;
                f1->tulip[kinetic4].linkNext  = nullName;

                break;
            }
           // assignSplit( *f1, h12 ,f->i.nOperator, f1->vectorOperator, f1->vectorOperator+f->i.nOperator );

            
            //really want to keep linear separate and uncarved...
            if ( c1->i.Na && c1->rt.calcType != clampProtonElectronCalculation )
                if ( c1->i.oneBody.func.fn != nullFunction ){
                    tClear(*f1, linear);
                    if(   ! ioStoreMatrix(*f1, linear, 0, "trainLinear.matrix",1))
                        if ( ! ioStoreMatrix(*f1,linear ,0,"linear.matrix",1)) {
                            printf("failed to load either external or distilled External\n");
                            exit(0);//if not already present, go back and build it
                        }
                    if ( f1->cmpl == cmpl )
                        if(   ! ioStoreMatrix(*f1, linear, 1, "trainLinear.1.matrix",1))
                            if ( ! ioStoreMatrix(*f1,linear ,1,"linear.matrix",1)){
                                printf("failed to load either external or distilled External\n");
                                exit(0);//if not already present, go back and build it
                            }

                }

        }
//        else {
//            if ( bootBodies == one && ( c1->rt.calcType == electronicStuctureCalculation  )&& f1->chainFlag){
//                    f1->tulip[external1].species = matrix;
//
//                f1->tulip[Ha].linkNext = kinetic1;
//                f1->tulip[Iterator].linkNext = kinetic1;
//                f1->tulip[kinetic1].linkNext = nullName;
//                tClear(*f1,kinetic1);
//                buildKinetic(c1, f1, 1., 1,tv1, kinetic1, electron, 0, real);
//
//                buildExternalPotential(c1, f1,1.,1,tv1,kinetic1,electron,!(!c1->rt.runFlag),real);
//
//                if (c1->i.springFlag)
//                    buildSHO(c1, f1, c1->i.springConstant,1, tv1, kinetic1, electron, 0, real);
//
//            }
//            if ( bootBodies == two && ( c1->rt.calcType == electronicStuctureCalculation  )&& f1->chainFlag){
//
//                if( 1 ){
//                    f1->tulip[kinetic1].chainNext  = kinetic1_2;
//                    f1->tulip[kinetic1].linkNext  = kinetic2;
//                    f1->tulip[kinetic2].chainNext  = kinetic2_2;
//                    f1->tulip[kinetic2].linkNext  = external1;
//
//
//                    if ( 1 ){
//                    f1->tulip[external1].chainNext  = external1_2;
//                    f1->tulip[external1].linkNext  = external2;
//                    f1->tulip[external2].chainNext  =external2_2;
//                    f1->tulip[external2].linkNext  = interaction12;
//                    f1->tulip[interaction12].chainNext = interaction12_2;
//                    } else{
//                        f1->tulip[external1].linkNext  = external2;
//
//                        f1->tulip[external2].linkNext  = interaction12;
//
//                    }
//                }else {
//
//                        f1->tulip[kinetic1].linkNext  = kinetic1_2;
//                        f1->tulip[kinetic1_2].linkNext  = kinetic2;
//                        f1->tulip[kinetic2].linkNext  = kinetic2_2;
//                        f1->tulip[kinetic2_2].linkNext  = external1;
//                    //    f1->tulip[external1].chainNext  = external1_2;
//                       f1->tulip[external1].linkNext  = external2;
//                //    f1->tulip[external2].chainNext  = external2_2;
//                    f1->tulip[external2].linkNext = interaction12;
//                    f1->tulip[interaction12].chainNext = interaction12_2;
//
//
//                }
//
////
////
////                                    f1->tulip[kinetic1].linkNext = kinetic2;
////
////            //HERE
//
//
//
//
////
//                                    f1->tulip[Ha].linkNext = kinetic1;
//                                    f1->tulip[Iterator].linkNext = kinetic1;
//                                    buildPairWisePotential(c1, f1,1.0,1,e12,interaction12,electron, 0,real);
//                                 //   buildPairWisePotential(c1, f1,-0.5,2,e12,interaction12,electron, 0,real);
//
//                buildExternalPotential(c1, f1,1.,1,tv1,external1,electron,!(!c1->rt.runFlag),real);
//                               //  buildExternalPotential(c1, f1,0.5,2,tv1,external1,electron,!(!c1->rt.runFlag),real);
//                buildExternalPotential(c1, f1,1.,tv2,1,external1,electron,!(!c1->rt.runFlag),real);
//                //  buildExternalPotential(c1, f1,0.5,tv2,1,external2,electron,!(!c1->rt.runFlag),real);
//
//                if ( c1->i.springFlag ){
//                    buildSHO(c1, f1, c1->i.springConstant,1, tv1, external2, electron, 0, real);
//                    buildSHO(c1, f1, c1->i.springConstant,1, tv2, external2, electron, 0, real);
//                }
//
//
//
//                if ((1)){
//                    buildKinetic(c1, f1, 1.,1, tv1, kinetic1, electron, 0, real);
//                    buildKinetic(c1, f1, 1.,1, tv2, kinetic1, electron, 0, real);
////                                    buildKinetic(c1, f1, 0.5,1, tv2, kinetic2, electron, 0, real);
////                                    buildKinetic(c1, f1, 0.5,2, tv2, kinetic2_2, electron, 0, real);
//                }
//
//            }else           if ( bootBodies == three && ( c1->rt.calcType == electronicStuctureCalculation  )&& f1->chainFlag){
//
//
//
//                              //      f1->tulip[kinetic1].linkNext = kinetic2;
//                              //      f1->tulip[kinetic2].linkNext = kinetic3;
//                f1->tulip[interaction12].linkNext = interaction13;
//                f1->tulip[interaction13].linkNext = interaction23;
//                //HERE
//                                    f1->tulip[Ha].linkNext = interaction12;
//                                    f1->tulip[Iterator].linkNext = interaction12;
//                buildPairWisePotential(c1, f1,1.,1,e12,interaction12,electron, 0,real);
//                buildPairWisePotential(c1, f1,1.,1,e13,interaction13,electron, 0,real);
//                buildPairWisePotential(c1, f1,1.,1,e23,interaction23,electron, 0,real);
//
//
//
//
//
//                buildExternalPotential(c1, f1,0.5,1,tv1,interaction12,electron,!(!c1->rt.runFlag),real);
//                buildExternalPotential(c1, f1,0.5,1,tv1,interaction13,electron,!(!c1->rt.runFlag),real);
//                buildKinetic(c1, f1, 0.5,1, tv1, interaction12, electron, 0, real);
//                buildKinetic(c1, f1, 0.5,1, tv1, interaction13, electron, 0, real);
//
//
//                buildExternalPotential(c1, f1,0.5,1,tv2,interaction12,electron,!(!c1->rt.runFlag),real);
//                buildExternalPotential(c1, f1,0.5,1,tv2,interaction23,electron,!(!c1->rt.runFlag),real);
//                buildKinetic(c1, f1, 0.5,1, tv2, interaction12, electron, 0, real);
//                buildKinetic(c1, f1, 0.5,1, tv2, interaction23, electron, 0, real);
//
//
//
//
//                buildExternalPotential(c1, f1,0.5,1,tv3,interaction13,electron,!(!c1->rt.runFlag),real);
//                buildExternalPotential(c1, f1,0.5,1,tv3,interaction23,electron,!(!c1->rt.runFlag),real);
//                buildKinetic(c1, f1, 0.5,1, tv3, interaction13, electron, 0, real);
//                buildKinetic(c1, f1, 0.5,1, tv3, interaction23, electron, 0, real);
//
//
//                if (0){
//                    double scalar = -9/f->i.d/f->i.d;
//                buildConstant(c1, f1, 2.*scalar,1, tv1, Iterator, electron, 0, real);
//                buildConstant(c1, f1, 2.*scalar,1, tv2, Iterator, electron, 0, real);
//                buildConstant(c1, f1, 2.*scalar,1, tv3, Iterator, electron, 0, real);
//                buildConstant(c1, f1, -0.5*scalar,4, tv1,Iterator, electron, 0, real);
//                buildConstant(c1, f1, -0.5*scalar,4, tv2, Iterator, electron, 0, real);
//                buildConstant(c1, f1, -0.5*scalar,4, tv3, Iterator, electron, 0, real);
//                buildConstant(c1, f1, -0.5*scalar,5, tv1, Iterator, electron, 0, real);
//                buildConstant(c1, f1, -0.5*scalar,5, tv2, Iterator, electron, 0, real);
//                buildConstant(c1, f1, -0.5*scalar,5, tv3, Iterator, electron, 0, real);
//                }
//                {
//                    double scalar = 1.;
////                    buildKinetic(c1, f1, 1.,1, tv1, kinetic1, electron, 0, real);
////                    buildKinetic(c1, f1, 1.,1, tv2, kinetic2, electron, 0, real);
////                    buildKinetic(c1, f1, 1.,1, tv3, kinetic3, electron, 0, real);
//                    if ( c1->i.springFlag){
//                    buildSHO(c1, f1, c1->i.springConstant/2,1, tv1, interaction12, electron, 0, real);
//                        buildSHO(c1, f1, c1->i.springConstant/2,1, tv2, interaction12, electron, 0, real);
//
//
//
//                        buildSHO(c1, f1, c1->i.springConstant/2,1, tv1, interaction13, electron, 0, real);
//                    buildSHO(c1, f1, c1->i.springConstant/2,1, tv3, interaction13, electron, 0, real);
//
//
//                        buildSHO(c1, f1, c1->i.springConstant/2,1, tv2, interaction23, electron, 0, real);
//                    buildSHO(c1, f1, c1->i.springConstant/2,1, tv3, interaction23, electron, 0, real);
//                    }
////                    buildKinetic(c1, f1, 2.*scalar,1, tv1, kinetic1, electron, 0, real);
////                    buildKinetic(c1, f1, 2.*scalar,1, tv2, kinetic1, electron, 0, real);
////                    buildKinetic(c1, f1, 2.*scalar,1, tv3, kinetic1, electron, 0, real);
////                    buildKinetic(c1, f1, -0.5*scalar,4, tv1, kinetic1, electron, 0, real);
////                    buildKinetic(c1, f1, -0.5*scalar,4, tv2, kinetic1, electron, 0, real);
////                    buildKinetic(c1, f1, -0.5*scalar,4, tv3, kinetic1, electron, 0, real);
////                    buildKinetic(c1, f1, -0.5*scalar,5, tv1, kinetic1, electron, 0, real);
////                    buildKinetic(c1, f1, -0.5*scalar,5, tv2, kinetic1, electron, 0, real);
////                    buildKinetic(c1, f1, -0.5*scalar,5, tv3, kinetic1, electron, 0, real);
//
//                }
//        }
//        }
    //7.7
#if VERBOSE
    printf("boot complete\n");
#endif
    fflush(stdout);
    return 0;
}