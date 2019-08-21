/*
 *  input.c
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

#include "input.h"

INT_TYPE add_atom( struct input *f1, INT_TYPE l, double x1,double y1,double z1 ){
    INT_TYPE n = f1->Na;
    INT_TYPE Z1,label ;
    if ( n+1 == MAXATOM ){
        fprintf(stderr,"TOO MANY ATOMS");
        return -1;
    }
    f1->Na++;
    n++;
    f1->atoms[n].position[1] = x1;
    f1->atoms[n].position[2] = y1;
    f1->atoms[n].position[3] = z1;
    
    Z1 = ((l)%MAX_PROTOTYPE_ATOMS);
    label = labs(l)/MAX_PROTOTYPE_ATOMS;
    
    f1->atoms[n].label.Z = (Z1);
    f1->atoms[n].label.label = label;
    
    return 0;
};


INT_TYPE centerMass ( struct calculation * c1){
    INT_TYPE alpha;
    INT_TYPE a;
    INT_TYPE l;
    double cm[] = {0.,0.,0.,0.};
    double wcm[] = {0.,0.,0.,0.};
    
    
    
    for ( alpha = 1 ; alpha <= c1->i.Na ; alpha++){
        for( l = 1 ; l <= 3 ; l++){
            cm[l] += c1->i.atoms[alpha].label.Z*c1->i.atoms[alpha].position[l];
            wcm[l] += c1->i.atoms[alpha].label.Z;
        }
    }
    
    cm[1] /= wcm[1];
    cm[2] /= wcm[2];
    cm[3] /= wcm[3];
    
    for ( a = 1 ; a <= c1->i.Na ;a++){
        for ( l = 1 ; l<= 3 ; l++){
            c1->i.atoms[a].position[l] -= cm[l];
            //            printf("%f\n",cm[l]);
        }
    }
    
    return 0;
}

INT_TYPE comment (const char * line ){
   return (NULL != strstr(line, "#"));
}

INT_TYPE control ( const char * line ){
    INT_TYPE m,s;
    INT_TYPE Ncom = 6;
    char *list_control []= {"#",
        "Body",
        "Parameters",       "Geometry",
        "Translation" ,       "Basis",
        "InputOutput",
    };

    char *list_line [] = {".","*"};
    char line0 [MAXSTRING];
    char * line1 ;

        for (m = 1 ; m <= Ncom ; m++){
            
            for ( s = 0; s < 2 ; s++){
                
                sprintf(line0, "%s%s", list_line[s],list_control[m] );
                line1 = strstr( line, line0 );
                if ( line1 != NULL ){
                    return (2*s-1)*m;
                }
            }
        }
 
        return 0;
    
}


INT_TYPE getParam ( struct calculation * c,struct input_label *f1, const char * input_line ){
    INT_TYPE i,d,ivalue;
    char test_line [MAXSTRING];
    double value;
    INT_TYPE NINT_TYPE = 119;
    char *list_INT_TYPE []= {"#",
        "LOST1","maxCycle" , "spinor", "charge","fineStr",//5
        "process", "NB", "MB", "percentFull","general",//10
        "center","xTranslate","yTranslate","zTranslate","postCalc",//15
        "goK","goV","goC","goX","goS",//20
        "iGold","LOST", "LOST2" ,"core","canon",//25
        "pseudo","minDIIS","LOST","weylet", "nWeylet",//30
        "mWeylet","helium","correlation","initRank","LOST3",//35
        "spinBlocks","LOST5", "LOST6", "maxLevelShift","diis",//40
        "minSPC","maxEV","inverseQuad","maxSPC","maxDIIS",//45
        "intDIIS", "trace","basisRank","LUMO","foundation",//50
        "Angstroms","train", "ceilFlag","iSymm","type",//55
        "zone","eZone","cycles","weightRank","printConvergence",//60
        "runFlag","LOST10","annulus","exclusion","runType",//65
        "SpinSqr","setRange","range","golden","bandStage",//70
        "boot","jCycle","fermi","signature","filter",//75
        "minRank","skipBuild","printLevel","stack","lanes",//80
        "sectors","body","LOST100","rds1","rds2",//85
        "rds3","interactionOne","interactionTwo","oCycle","interactionZero",//90
        "breakBody","interval","RAM","monteCarlo","samples",//95
        "hartreeFock","basisStage","iterations","collect","states",//100
        "length","XHA","lookBack","step","theory",//105
        "configuration","densityRank","densityBody","parallel","phase",//1101
        "around","cmpl","clampStage","OCSB","decompose",
        "shiftNO","matrix","catalog","increment"
    };
    INT_TYPE NDOUBLE = 76;
    char *list_DOUBLE []= {"#",
        "lattice","mix", "aoDirectDensity","aoExchangeDensity", "LOST" ,//1-5
        "xB", "yB", "zB", "xyRange" , "zRange",//6-10
        "XX", "scfTolerance","boundTolerance","cycleTolerance","oExternal",//11-15
        "LOST","LOST", "XXz" ,"LOST4", "aWeylet" ,//16-20
        "bWeylet","levelLevel","tolerance","threshold","target",//21-25
        "convergence","external","vectorThreshold","buildThreshold","maxPi",//26-30
        "EMPTY","EMPTY", "vectorConvergence","powState", "powVacuum",//35
        "mass1", "mass2","Charge1", "Charge2","Charge3",//40
        "beta","EMPTY10","ceilValue","floorValue","electronGasDensity",//45
        "shift","EMPTIER","crystal","jelliumRadius","spring",//50
        "REMOVEREMOVE", "maxDomain", "parcel","minDomain","param",
        "entropy","attack","scalar","turn","augment",
        "linearDependence","condition","seek","width","latte",
        "magnetismZ","clampMin","clampMax","electronMass","protonMass",
        "pairMass","gamma0","ewald","levelScale","scaleVectorThreshold",
        "scaleTarget"
    };
    
    for ( i = 1 ; i <= NINT_TYPE ; i++){
        
        if ( strstr( input_line, list_INT_TYPE[i])!= NULL){
#ifdef APPLE
            sscanf(input_line,"%s %d", test_line,&ivalue);
#else
            sscanf(input_line,"%s %lld", test_line,&ivalue);
#endif
            
            switch ( i ){
                case 1 :
                    if ( ivalue > 0 )
//                    {
//                        c->p.iTolerance = ivalue ;
//                        return i;
//                    } else
                        return 0;
                case 2 :
                    if ( ivalue > 0) {
             //       c->i.Ncycle = ivalue;
                    return i;
                    } else
                        return 0;
                case 3 :
                 //   c->i.spinor = ivalue;
                    return i;
                    
//                case 4 :
//                    c->i.charge = ivalue;
//                    return i;
                    
//                case 5 :
//                    if ( ivalue == 0 || ivalue == 1) {
//                    c->i.fineStructure = ivalue;
//                    return i;
//                    } else{
//                        return 0;
//                    }
//                case 6:
//                    if ( ivalue > 0 ){
//                        c->i.nProcDivison = ivalue;
//                    return i;
//                    } else
//                        return 0;
//                case 7 :
//                    if ( ivalue > 0){
//                    c->i.NB = ivalue;
//                    return i;
//                    } else
//                        return 0;
//                case 8 :
//                    if ( ivalue > 0){
//                        c->i.MB = ivalue;
//                        return i;
//                    } else
//                        return 0;
                    
//                case 9 :
//                    if ( ivalue >0 && ivalue <= 100 ){
//                    c->i.percentFull = ivalue;
//                    return i;
//                    }else
//                        return 0;
////                case 10 :
//                    if ( ivalue == 0 || ivalue == 1){
//                    c->i.General = ivalue;
//                    return i;
//                    } else
//                        return 0;
                    
                case 11 :
                    if ( ivalue == 0 || ivalue == 1 || ivalue == 2){
                   // c->i.centerMass = ivalue;
                        
                        //if ( c->i.centerMass == 2 )
                        //    centerPosition( c);
                        if ( ivalue == 1 )
                            centerMass( c);

                        
                        
                     return i;
                    } else
                        return 0;
//                case 12 :
//                    c->i.translate[1] = ivalue;
//                    return i;
//                   
//                case 13 :
//                    c->i.translate[2] = ivalue;
//                    return i;
//
//                case 14 :
//                    c->i.translate[3] = ivalue;
//                    return i;

//                case 15 :
//                    if ( ivalue >= 0 ){
//                        c->i.postCalc = ivalue;
//                        return i;
//                    } else
//                        return 0;
//                case 16 :
//                    if ( ivalue == 0 || ivalue == 1){
//
//                    c->i.goKinetic = ivalue;
//                    return i;
//                    } else return 0;
//                    
//                case 17 :
//                    if ( ivalue == 0 || ivalue == 1){
//                        
//                        c->i.goExternal = ivalue;
//                        return i;
//                    } else return 0;
//
//                case 18 :
//                    if ( ivalue == 0 || ivalue == 1){
//                        
//                        c->i.goDirect = ivalue;
//                        return i;
//                    } else return 0;
//
//                                   case 19 :
//                    if ( ivalue == 0 || ivalue == 1){
//                        
//                        c->i.goExchange = ivalue;
//                        return i;
//                    } else return 0;
//
                   
//                case 20 :
//                    if ( ivalue == 0 || ivalue == 1){
//                        
//                        c->i.goSpin = ivalue;
//                        return i;
//                    } else return 0;
            
//                case 21 :
//                    
//                        c->i.iGold = ivalue;
//                        return i;
                    
                case 22 :
                    
                    //c->p.iCondition = ivalue;
                    return i;
                case 23 :
                   // c->p.iThreshold = ivalue;
                    return i;
                case 24 :
                //   c->i.qCore  = ivalue;
                    return i;
                case 25 :
                    c->i.canonRank  = ivalue;
                    return i;
                case 27 :
                  //  c->i.minDIIS = ivalue;
                    return i;
                case 28 :
                 //   c->i.iCharge = ivalue;
               //     return i;
                case 29:
                  //  c->i.weyletFlag = ivalue;
                    return i;
//                case 30:
//                    c->i.nWeylet = ivalue;
//                    return i;
//                case 31:
//                    c->i.mWeylet = ivalue;
//                    return i;
                case 32:
//                    c->i.heliumFlag = ivalue;
                    return i;
//                case 33:
//                    c->i.correlationFlag = ivalue;
//                    return i;
                case 34:
                    f1->iRank = ivalue;
                    return i;
                case 35:
                ///    c->p.iTarget = ivalue;
                    return i;
                case 36:
                  //  c->p.iSpin = ivalue;
                    return i;
                case 37 :
//                    c->p.iXTolerance = ivalue;
                    return i;
                case 38 :
//                    c->p.iXThreshold = ivalue;
                    return i;
                case 39 :
                //    c->i.maxLevelShift = ivalue;
                    return i;
                case 40 :
                 //   c->i.diis = ivalue;
                    return i;
//                case 41 :
//                    c->i.minSPC = ivalue;
//                    return i;
//                case 42 :
//                    c->i.maxEV = ivalue;
//                    return i;
//                case 43 :
//                    c->i.ar = ivalue;
//                    return i;
//                case 44 :
//                    c->i.maxSPC = ivalue;
//                    return i;
                case 45 :
                //    c->i.maxDIIS = ivalue;
                    return i;
                case 46 :
                  //  c->i.intDIIS = ivalue;
                    return i;
//                case 47 :
//                    c->p.iTrace = ivalue;
//                    return i;
                case 48 :
                    f1->bRank = ivalue;
                    return i;
                case 49 :
              //      c->i.lumos = ivalue;
                    return i;
                case 50 :
                    f1->qFloor = ivalue;
               //     f1->sectors = !(!ivalue);
                    return i;
                case 51 :
                    c->i.Angstroms = ivalue;
                    return i;
                case 52 :
               //     c->i.trainFlag = ivalue;
                    return i;
                case 53 :
                //    c->i.qCeil = ivalue;
                    return i;
                case 54 :
              //      c->i.inversionSymm = ivalue;
                    return i;
                case 55 :
                    f1->irrep = ivalue;
                    return i;
                case 56 :
                //    c->i.zone = ivalue;
                    return i;
                case 57 :
                 //   c->i.eZone = ivalue;
                    return i;
                case 58 :
//                    c->i.cycles = ivalue;
                    return i;
                case 59 :
                 //   c->i.wRank  = ivalue;
                    return i;
                case 60 :
                 //   c->i.checkFlag = ivalue;
                    return i;
                case 61 :
                    c->rt.runFlag = ivalue;
                    return i;
//                case 62 :
//                    c->i.resolveState = ivalue;
//                    return i;
                case 63 :
               //     c->i.annulus = ivalue;
                    return i;
                case 64 :
               //     c->p.iExclusion = ivalue;
                    return i;
                case 65 :
              //      c->i.runType = ivalue;
                    return i;
                case 66 :
               //     c->i.SpinSqr = ivalue;
                    return i;
                case 67 :

                    f1->d /= (2*ivalue+1)/(2*f1->epi+1);
                    f1->epi = ivalue;
                    
                    
                    
                    
                    
                    
                    
                    return i;
                    
                case 68 :

                    f1->epi = ivalue;
                    
                    
                    return i;

                    
                case 69 :
                    //c->i.runCharacter = ivalue;
                    return i;
                    
                    
                case 70:
                    //expand range by number given.
                    f1->d *= pow( (2.* f1->epi + 1.) /(2.*f1->epi + 2*ivalue + 1),f1->attack);
                    f1->epi  += ivalue;
                    
                    

                    return i;

                case 71:
         //           c->i.bootRestriction = ivalue;
                    return i;
             
                case 72:
                 //   c->i.Jcycle = ivalue;
                    return i;

                case 73:
                //    c->i.fermi = ivalue;
                    return i;

                case 74:
               //     c->i.singleSignature = ivalue;
                    return i;

                case 75:
                    f1->filter = ivalue;
                    return i;
                    
                case 76:
            //        c->i.minRank = ivalue;
                    return i;

                case 77:
               //     c->i.overBuildFlag = ivalue;
                    return i;

                case 78 :
//                    if ( c->rt.printFlag ){
//
//                        c->rt.printFlag = ivalue;
//
//                        return i;
//                    }
                    return 0;
                    
                case 79:
//                    if ( c->i.iCharge){
//                        if ( c->i.iCharge == 2 ){
//                            c->i.heliumFlag = ivalue * (ivalue-1)/2;
//                            c->i.iCharge = 2 * ivalue;
//                        }else if ( c->i.iCharge == 3 ){
//                            c->i.heliumFlag = ivalue * (ivalue-1)* (ivalue-2)/6;
//                            c->i.iCharge = 3 * ivalue;
//                        }
//
//
//
                    
                        
                        
                        
//                        return i;
//                    }
                    return 0;
                    
                case 80:
                    if ( ivalue < 0 )
                        return 0;
#ifdef OMP
                    c->i.omp = ivalue;
#endif
                    return i;

                case 81:
                    if ( ivalue < 0 )
                        return 0;
                  //  f1->sectors = ivalue;
                    return i;
                    
                case 82:
                    if ( ivalue == 0 )
                        f1->body = nada;
                  else   if ( ivalue == 1 )
                        f1->body = one;
                    else if ( ivalue == 2 )
                        f1->body = two;
                    else if ( ivalue == 3 )
                        f1->body = three;
                    else if ( ivalue == 4 )
                        f1->body = four;
                    else
                    {
                        return 0;
                    }
                    return i;

                case 83:
                    if ( ivalue < 0 )
                        return 0;
                   // c->i.traffic = ivalue;
                    return i;

                case 84:
                    if ( ivalue < 0 )
                        return 0;
                    
//                    if ( c->rt.body == one ){
//                        c->i.c.rds.Basis[0] = ivalue;
//                        c->i.c.rds.flag = 1;
//                    }else if ( c->rt.body == two ){
//                        c->i.c.rds.Basis2[0] = ivalue;
//                        c->i.c.rds.flag = 2;
//                    }else if ( c->rt.body == three ){
//                        c->i.c.rds.Basis3[0] = ivalue;
//                        c->i.c.rds.flag = 3;
//                    }else if ( c->rt.body == four ){
//                        c->i.c.rds.Basis4[0] = ivalue;
//                        c->i.c.rds.flag = 4;
//                    }
                    return i;
                case 85:
                    if ( ivalue < 0 )
                        return 0;
//                    if ( c->rt.body == one ){
//                        c->i.c.rds.Basis[1] = ivalue;
//                        c->i.c.rds.flag = 1;
//
//                    }else if ( c->rt.body == two ){
//                        c->i.c.rds.Basis2[1] = ivalue;
//                        c->i.c.rds.flag = 2;
//
//                    }else if ( c->rt.body == three ){
//                        c->i.c.rds.Basis3[1] = ivalue;
//                        c->i.c.rds.flag = 3;
//
//                    }else if ( c->rt.body == four ){
//                        c->i.c.rds.Basis4[1] = ivalue;
//                        c->i.c.rds.flag = 4;
//
//                    }
                    return i;
                case 86:
                    if ( ivalue < 0 )
                        return 0;
//                    if ( c->rt.body == one ){
//                        c->i.c.rds.Basis[2] = ivalue;
//                        c->i.c.rds.flag = 1;
//
//                    }else if ( c->rt.body == two ){
//                        c->i.c.rds.Basis2[2] = ivalue;
//                        c->i.c.rds.flag = 2;
//
//                    }else if ( c->rt.body == three ){
//                        c->i.c.rds.Basis3[2] = ivalue;
//                        c->i.c.rds.flag = 3;
//
//                    }else if ( c->rt.body == four ){
//                        c->i.c.rds.Basis4[2] = ivalue;
//                        c->i.c.rds.flag = 4;
//
//                    }
                    return i;
                
                case 87:
                        c->i.oneBody.func.interval = c->i.interval;
                        c->i.oneBody.func.fn = ivalue;
                        c->i.oneBody.func.param[0] = c->i.scalar;
                        c->i.oneBody.func.param[1] = c->i.turn;
                        c->i.oneBody.func.param[2] = c->i.param1;
                        c->i.oneBody.func.param[3] = c->i.param2;

                        return i;

                    
                case 88:
                        c->i.twoBody.func.interval = c->i.interval;
                        c->i.twoBody.func.fn = ivalue;
                        c->i.twoBody.func.param[0] = c->i.scalar;
                        c->i.twoBody.func.param[1] = c->i.turn;
                        c->i.twoBody.func.param[2] = c->i.param1;
                        c->i.twoBody.func.param[3] = c->i.param2;

                        return i;
                case 89:
                   // c->i.OutCycle = ivalue;
                    return i;

                case 90:
//                    c->i.c.zroBody.func.interval = c->i.interval;
//                    c->i.c.zroBody.func.fn = ivalue;
//                    c->i.c.zroBody.func.param[0] = c->i.param0;
//                    return i;
                case 91:
                 //   c->i.entropyFlag = ivalue;
                    return i;
                case 92:
                    c->i.interval = ivalue;
                    return i;
                case 93:
                    c->i.RAMmax = ivalue;
                    return i;
                case 94:
//                    c->rt.monteCarlo = ivalue;
//                    return i;
                case 95:
                 //   c->rt.samples = ivalue;
                    return i;
                case 96:
//                    c->i.hartreeFockFlag = ivalue;
                    return i;
                case 97:
                    f1->bRank += ivalue;
                    return i;
                case 98:
                   f1->Iterations = ivalue;
                    return i;
                case 99:
                    f1->collect = ivalue;
                    return i;
                case 100:
                    f1->nStates = ivalue;
                    return i;
                case 101:
               //     c->i.l2 = ivalue;
              //      return i;
//                case 102:
//                    c->i.side = ivalue;
//                    return i;
                case 103:
                 //   c->i.lookBack = ivalue;
                    return i;
                case 104:
                   // c->i.cycleStep = ivalue;
                    return i;
                case 105:
                   // c->i.theory = ivalue;
                    return i;
                case 106:
                    if ( ivalue == 0 )
                        c->rt.calcType = electronicStuctureCalculation;
                    else if ( ivalue == 1 )
                        c->rt.calcType = clampProtonElectronCalculation;
                    else
                    {
                        return 0;
                    }
                    return i;
                case 107:
                    return i;
                case 108:
                    return 0;
                case 109:
#ifdef MKL
                    c->i.mkl = ivalue;
#endif
                    return i;
                case 110:
                    c->rt.phaseType = ivalue;
                    return i;
                    
                case 111:
                    c->rt.calcType = 2;
                    c->i.orgClamp = 0.5*(c->i.minClamp+c->i.maxClamp);
                    f1->around = floor((c->i.orgClamp-c->i.minClamp)/f1->D);
                    return i;
                    
                case 112:
                    f1->cmpl = ivalue;
                    return i;
                    
                case 113:
                    f1->D *= pow( (2.* f1->around + 1.) /(2.*f1->around + 2*ivalue + 1),1.);
                    f1->around  += ivalue;
                    return i;
                    
                case 114:
            //        c->i.OCSBflag = ivalue;
                    return i;
                case 115:
                    c->rt.powDecompose = ivalue;
                    return i;
                case 116:
                    c->i.shiftFlag = 0;
                    return i;
                case 117:
                    c->i.decomposeRankMatrix = ivalue;
                    return i;
                case 118:
                    f1->cat = ivalue;
                    return i;
                case 119:
                    c->i.decomposeRankMatrix += ivalue;
                    return i;

            }
        
        }
    
    }


    for ( d = 1; d <= NDOUBLE ; d++){
        
        sprintf(test_line ,"%s", list_DOUBLE[d]);
        if ( strstr( input_line, test_line)!= NULL){
            sscanf(input_line,"%s %lf", test_line, &value);
            
            switch ( d ){
                case 1 :
                    if ( value > 0 ){
                        f1->d = value;
                        return d;
                    } else {
                        f1->d = 1;
                        f1->epi = 0;
                        return d;
                    }
                case 2 :
                 //   c->i.mix = value;
                    if ( value < 0 ){
                 //       c->i.mixState = 0;
                    }
                    else{
                     //   c->i.mixState = 1;
                    }
                    
                    return d;
                    
                    //                case 3 :
                    //                    c->i.ddOffSpace = value;
                    //                    return d;
                    //
                    //                case 4 :
                    //                    c->i.dxOffSpace = value;
                    //                    return d;
                    //
                    //                case 5 :
                    //                    c->i.hdOffSpace = value;
                    //                    return d;
                    
                    //                case 6 :
                    //                    c->i.c.bf[1] = value;
                    //                    return d;
                    //
                    //                case 7 :
                    //                    c->i.c.bf[2] = value;
                    //                    return d;
                    //
                    //                case 8 :
                    //                    c->i.c.bf[3] = value;
                    //                    return d;
                    
                case 9 :
                    // c->i.xy_epi = value;
                    // return d;
                    
                case 10:
                   f1->epi = value;
                    return d;
                    
                case 11 :
                    //                    c->i.xy_epi = value;
                    //                    c->i.z_epi = value;
                    // return 0;
                    
                    //                case 12:
                    //                    c->i.scfTolerance = value;
                    //                    return d;
                    //
                    //                case 13:
                    //                    c->i.evTolerance = value;
                    //                    return d;
                    //                case 14:
                    //                    c->i.sp2Tolerance = value;
                    //                    return d;
                    //                case 15 :
                    //                    c->i.xOffSpace = value;
                    //                    return d;
                    //                case 16 :
                    //                    c->i.hxOffSpace = value;
                    //                    return d;
                    //                case 17 :
                    //                    c->i.xxOffSpace = value;
                    //                    return d;
                case 18 :
                    f1->epi *= f1->d / value ;
                    f1->d = value;
                    return d;
                case 19 :
                    //   c->p.sTarget = value;
                    return d;
                    //                case 20 :
                    //                    c->i.aWeylet = value;
                    //                    return d;
                    //                case 21 :
                    //                    c->i.bWeylet = value;
                    //                    return d;
                case 22 :
                    c->i.level = value;
                    return d;
                case 23 :
                    c->rt.TOL = pow(10., value);
                    return d;
                case 24 :
                    c->rt.CANON = pow(0.1, value);
                    return d;
                case 25 :
                    c->rt.TARGET = pow(0.1,value);
                    return d;
                case 26 :
//                    c->rt.CONVERGENCE = c->rt.TARGET*pow(0.1, value);
//                    return d;
                case 27 :
                   // c->p.iExternal = value;
                    return d;
                case 28 :
                    c->rt.vCANON = pow(0.1, value);
                    return d;
                case 29 :
                  //  c->p.iBuild = value;
                    return d;
                case 30 :
              //      c->i.Pi = value;
                    return d;
                case 31 :
                //    c->i.cPi = value;
                    return d;
                    //                case 32 :
                    //                    c->i.fPi = value;
                    //                    return d;
                case 33 :
//                    c->rt.vCONVERGENCE = c->rt.TARGET*pow(0.1, value);
//                    return d;
                case 34 :
                //    c->i.powState = value;
                    return d;
                case 35 :
                ///    c->i.powVacuum = value;
                    return d;
                case 36 :
                //    c->i.mass1 = value;
                    return d;
                case 37 :
                //    c->i.mass2 = value;
                    return d;
                case 38 :
                //    c->i.charge1 = value;
                    return d;
                case 39 :
                 //   c->i.charge2 = value;
                    return d;
                case 40 :
                //    c->i.charge3 = value;
                    return d;
                case 41 :
                 //   c->i.beta = value;
                    return d;
//                case 42 :
//                    c->i.Aspect = value;
//                    return d;
                case 43 :
                //    c->i.valueCeil = value;
                    return d;
                case 44 :
                //    c->i.valueFloor = value;
                    return d;
                case 45 :
                   f1->d = pow( 4. * pi /3. *f1->body ,0.33333333333333333333)* value / (2.*f1->epi+1.);
                    //value = Rs

                    //4pi/3 Rs^3 = pow(d * N1 ,3.0)/f1->Ne
                    return d;
                case 46 :
                {
                    INT_TYPE i;
                    for ( i = 0; i < 100; i++){
                        c->i.shiftVector[i][0] = 1.;
                        c->i.shiftVector[i][1] = value;
                    }
                }
                    c->i.shiftFlag = 1;
                    c->i.realPart = value;
                    return d;
                case 47 :
                 //   c->i.shift = value * pi * pi / c->i.d / c->i.d ;
                    return d;
                case 48 :
                    printf("Crystal Momentum -universal-t%f\n", value);

                    value /= f1->d*( 2.*f1->epi+1.);
                    c->i.vectorMomentum = value;
                    return d;
                    
                case 49:
                    if ( value <= 0 )
                        return 0;
                    else{
                    //    c->i.jelliumSphereRadius = value;
//                        return d;
                    }
                    
                case 50:
                    c->i.springConstant = value;
                    c->i.springFlag = 1;
                    return d;


                case 52:
                  //  c->i.domain = value;
                    return d;

                case 53:
                 //   c->i.eps = value;
                    return d;
                case 54:
                 //   c->i.domainMin = value;
                    return d;
                case 55:
                    c->i.param1 = value;
                    return d;
                case 56:
                    c->rt.maxEntropy = pow(0.1,value);
                    return d;
                case 57:
                      f1->attack = value;
                    return d;
                case 58:
                    c->i.scalar = value;
                    return d;
                case 59:
                    c->i.turn = value;
                    return d;
                case 60:
                    c->i.param2 = value;
                    return d;
                case 61:
                   // c->rt.condition = pow(10., value);
                    return d;
                case 62:
                   // c->rt.targetCondition =  pow(0.1,value);
                    c->rt.ALPHA = pow(0.1,value);
                    return d;
                case 63:
                   // c->i.seekPower = value;
                    return d;
                case 64:
                    c->i.springFlag = 1;
                    c->i.springConstant = 1./(f1->d*value)/(f1->d*value);
//                    printf("spring %f (%f,%f)\n",c->i.springConstant,f1->d,value);
                    return d;
                case 65:
                    f1->D = value;
                    return d;
                case 66:
                    c->i.mag = value;
                    c->i.magFlag = 1;
                    return d;
                case 67:
                    c->i.minClamp = value;
                    return d;
                case 68:
                    c->i.maxClamp = value;
                    return d;
                case 69:
                    c->i.massElectron = value;
                    return d;
                case 70:
                    c->i.massProton = value;
                    return d;
                case 71:
                    c->i.massClampPair = value;
                    return d;
                case 72:
                    //Jacek Karwowski 06/17/2019
                    c->i.turn = pow(120.*(1.+sqrt(value/(20.*2*c->i.param1*c->i.param1*(c->i.massProton * c->i.massClampPair/(c->i.massProton + c->i.massClampPair) )*c->i.scalar))) ,1./6);
                    return d;
                case 73:
                    c->rt.EWALD = pow(0.1, value);
                    return d;
                case 74:
                    c->i.level = value/sqr( f1->d );
                    return d;
                case 75:
                    c->rt.vCANON *= value;
                    return d;
                case 76:
                    c->rt.TARGET *= value;
                    return d;
            }

        }
        
    }
    return 0;
}

INT_TYPE getGeometry(struct calculation * c, const char * input_line ){
    
    double x,y,z;
    INT_TYPE l,Z,label;
    
#ifdef APPLE
        if ( 4 == sscanf ( input_line, "%d %lf %lf %lf ", &l, &x,&y,&z ))
#else
        if ( 4 == sscanf ( input_line, "%lld %lf %lf %lf ", &l, &x,&y,&z ))
#endif
    {
        
        
        
        if ( l != 0 ){
            
            Z = ((l)%MAX_PROTOTYPE_ATOMS);
#ifdef APPLE
            label = abs(l)/MAX_PROTOTYPE_ATOMS;

#else
            label = llabs(l)/MAX_PROTOTYPE_ATOMS;

#endif
            
            if ( llabs(Z) >= MAX_PROTOTYPE_ATOMS || label >= MAX_LABEL )
            {
                return 1;
            }
            
            
            
            
            if ( c->i.Angstroms )
                add_atom ( &(c->i), l,x/a0,y/a0,z/a0 );
            else
                add_atom ( &(c->i), l,x,y,z );
//            
//            if ( Z > 12 ){
//                (c->i.c).atoms[(c->i.c).Na].label.iZ = Z -10;
//            }
//            else if ( Z > 4 ){
//                (c->i.c).atoms[(c->i.c).Na].label.iZ = Z -2;
//            } else {
//                (c->i.c).atoms[(c->i.c).Na].label.iZ = Z;
                
            }
//        if ( 1)  {
//                INT_TYPE Cc;
//                Cc = 0;
//                if ( Z > 2)
//                    Cc += 2;
//                if ( Z > 10)
//                    Cc += 8;
//                if ( Z > 18)
//                    Cc += 8;
//                if ( Z > 36){
//                    Cc += 18;
//                    exit(0);
//                }
//                if ( Z > 54){
//                    Cc += 18;
//                    exit(0);
//                }
//                (c->i.c).atoms[(c->i.c).Na].label.Cc = Cc;
//                
//            }

            
            
            return 0;
        }
        else {
            
            double u[] = {x,y,z};
            
            rotateGeometry(c,u);
            
            
            return 0;

        }
        
        
        
        
     
    
        return 1;
    
}


INT_TYPE rotateGeometry (struct calculation * c, double * u  ){
    double mag = cblas_dnrm2(3,u,1);
    INT_TYPE alpha;
    cblas_dscal(3,1./mag, u,1);
    mag *= 2*pi;
    
    printf("%f %f* %f\n",mag, cos(mag),sin(mag));
    double id [] = {1,0,0,  0,1,0,  0,0,1};
    double ux [] = {0, u[2],-u[1],  -u[2],0,u[0],   u[1],-u[0],0};
    double uu [] = {u[0]*u[0], u[0]*u[1],u[0]*u[2],   u[1]*u[0], u[1]*u[1], u[1]*u[2],   u[2]*u[0],u[2]*u[1],u[2]*u[2]};
    double R [9],x[3],y[3];
    cblas_dcopy(9,id,1,R,1);
    cblas_dscal(9,cos(mag),R,1);
    cblas_daxpy(9,sin(mag),ux,1,R,1);
    cblas_daxpy(9,1-cos(mag),uu,1,R,1);
    
    for ( alpha = 1; alpha <= c->i.Na ; alpha++)
    {
        x [0] = c->i.atoms[alpha].position[1];
        x [1] = c->i.atoms[alpha].position[2];
        x [2] = c->i.atoms[alpha].position[3];
        
        cblas_dgemv(CblasColMajor, CblasNoTrans,3,3,1.,R,3,x,1,0.,y,1);
        
        c->i.atoms[alpha].position[1] = y[0];
        c->i.atoms[alpha].position[2] = y[1];
        c->i.atoms[alpha].position[3] = y[2];

    }
    
    return 0;
    
}

INT_TYPE modGeometry(struct calculation * c, const char * input_line ){
    
    double x,y,z;
    INT_TYPE l;
    INT_TYPE input ;
#ifdef APPLE
    input =  sscanf ( input_line, "%d %lf %lf %lf ", &l, &x,&y,&z );

#else
   input =  sscanf ( input_line, "%lld %lf %lf %lf ", &l, &x,&y,&z );
#endif
    if ( input == 4 ){
        if ( l > 0 && l <= c->i.Na ){
            if ( c->i.Angstroms )
            {
                c->i.atoms[l].position[1] += x/a0;
                c->i.atoms[l].position[2] += y/a0;
                c->i.atoms[l].position[3] += z/a0;
            } else {
                {
                    c->i.atoms[l].position[1] += x;
                    c->i.atoms[l].position[2] += y;
                    c->i.atoms[l].position[3] += z;
                }
                
                
            }
            return 0;
        }
        else
            return 1;
    }
    return 0;
}

INT_TYPE intervalGeometry(struct calculation * c, const char * input_line ){
    
    double x,y,z;
    double x2,y2,z2;
    double x3,y3,z3;
    double lambda;
    INT_TYPE l,Z ,label;
    INT_TYPE input ;
#ifdef APPLE
    
    input =  sscanf ( input_line, "%lf %d %lf %lf %lf %lf %lf %lf ", &lambda, &l, &x,&y,&z, &x2,&y2,&z2 );

#else
    input =  sscanf ( input_line, "%lf %lld %lf %lf %lf %lf %lf %lf ", &lambda, &l, &x,&y,&z, &x2,&y2,&z2 );

#endif
    if ( input == 8 ){
        
        Z = ((l)%MAX_PROTOTYPE_ATOMS);
        
#ifdef APPLE
        label = abs(l)/MAX_PROTOTYPE_ATOMS;

#else
        label = llabs(l)/MAX_PROTOTYPE_ATOMS;
#endif
        
        if ( llabs(Z) >= MAX_PROTOTYPE_ATOMS || label >= MAX_LABEL )
        {
            return 1;
        }
        
        {
            
            x3 = (x+ lambda*(x2-x));
            y3 = (y+ lambda*(y2-y));
            z3 = (z+ lambda*(z2-z));
            if ( c->i.Angstroms )
                add_atom ( &(c->i), l,x3/a0,y3/a0,z3/a0 );
            else
                add_atom ( &(c->i), l,x3,y3,z3 );

        }
        
        
        return 0;
    }
    else
        return 1;

}


INT_TYPE getInputOutput(struct calculation * c,struct input_label * f1, const char * input_line ){
    INT_TYPE io;
    INT_TYPE Nio = 33;
    char test_line [MAXSTRING];
    char *list_IO[] = {"#",
        "densityIn","hartreeIn", "densityOut" ,//3
        "hartreeOut", "LOST" , "pAtomicTSW",//6
        "pWhen", "pAlt" , "read" ,//9
        "chdir","mkdir","externalIn",//12
        "externalOut","pMEM","pPotential",//15
        "xCoreMEM","xTotalMEM","pPercentFull",//18
        "pGold","pTime","nameDensityOut",//21
        "nameHartreeOut","Eigen" ,"set",//24
        "component","byHand","Spec",//27
        "vector","operator","print",//30
        "body","shift","twist"//33
    };
    char filename[MAXSTRING];
    FILE * mid;
    
    for( io = 1 ; io <= Nio ; io++){
        if ( strstr( input_line, list_IO [io])!=NULL){
            sscanf(input_line,"%s %s", test_line,  filename);
            
//            if ( strstr(filename , "*" ) != NULL )
//                sprintf(filename,"%s", c->name);
        switch ( io ){
            
//            case 1:
//                
//            
//                {
//                
//                    c->i.cpf = 1;
////                    sprintf(c->i.cpfName, "%s.%s",filename, DENSITY_SUFFIX);
//                }
//                return io;
//                
//            case 2 :
//                
//                {
//                    c->i.cpf = 2;
////                    sprintf(c->i.cpfName, "%s.%s",filename, HARTREE_SUFFIX);
//                }
//
//
//                return io;

            case 3 :
                

//                sprintf(c->i.nameOutDensity,"%s.%s", c->name,DENSITY_SUFFIX);
//                sprintf(c->i.oDensity,"%s",filename);

                
                
                return io;

                
            case 4 :
                
//                sprintf(c->i.nameOutHartree,"%s.%s", c->name,HARTREE_SUFFIX);
//                sprintf(c->i.oHartree,"%s",filename);
                return io;
             
            case 5:
                
               
//                sprintf(c->outFileName,"%s.%s", filename,OUTPUT_SUFFIX);
            
                return io;
                
            
            case 6:
                
//                sprintf(c->i.pTSW ,"%s", filename);
//                return io;
                
                //add outputs like prints... and Atomic core TSW.
                //linear plots...   define line segement and calculate
                
            case 7:
                
//                sprintf(c->i.pCalculation, "%s",filename);
//
//                return io;

            case 8:
//                
//                sprintf(c->i.pAlt,"%s", filename);
//                return io;
                
            case 9:
                
                mid = NULL;
                mid = fopen ( filename, "r");
                if( mid == NULL )
                    return 0;
                strcpy(test_line,c->name);
                if ( readInput( c, f1,mid ) )
                    return 0;
                strcpy(c->name, test_line);//does not inheret name of run
                fclose( mid);
                return io;
                
            case 10:
                
           //     if ( !chdir( filename ) )
            //        return 0;
            //    else
                    return io;
            case 11:
                
            //    if ( !mkdir( filename , 0777) )
            ///        return 0;
            ///    else
                    return io;
                
            case 12:
                
                
//                sprintf(c->i.cpfName, "%s.%s",filename, EXTERNAL_SUFFIX);
//                if ( c->i.cpf / 4 == 0 )
//                    c->i.cpf += 4;
                return io;
                
            case 13 :
                
//                sprintf(c->i.oufName, "%s.%s",filename,EXTERNAL_SUFFIX);
//                if ( c->i.ouf / 4 == 0 )
//                    c->i.ouf += 4;
//
                return io;

//            case 14:
//                
//                sprintf(c->i.pMEM,"%s", filename);
//                return io;
                
//            case 15:
//                
//                sprintf(c->i.pPotential ,"%s", filename);
//                return io;
//                
//            case 16:
//                
//                sprintf(c->i.xCoreMEM,"%s",filename);
//                return io;
//                
//            case 17:
//                
//                sprintf(c->i.xTotalMEM,"%s",filename);
//                return io;
//             
//                
//            case 18:
//                
//                sprintf(c->i.pPercentFull,"%s", filename);
//                
//                return io;
//            case 19:
//                
//                sprintf(c->i.pGold,"%s", filename);
//                
//                return io;
//            case 20:
//                
//                sprintf(c->i.pTime,"%s", filename);
//                
//                return io;
//                
//            case 21:
//             //   sprintf(c->i.nameOutDensity,"%s.%s", filename,DENSITY_SUFFIX);
//                return io;
//                
//            case 22:
//               // sprintf(c->i.nameOutHartree,"%s.%s", filename,HARTREE_SUFFIX);
//                return io;
            case 23:
            
                return io;
            case 24:{
                return io;
            }
            case 25:
            {
                return io;

            }
                
            case 26:
            {
               // c->i.entropyFlag = ! c->i.entropyFlag;
//                c->i.buildLength = 1;
//                sprintf(c->i.buildName,"%s", filename);
                return io;
            }
                
                
                
            case 27:
            {
                struct dirent *de;  // Pointer for directory entry
                
                // opendir() returns a pointer of DIR type.
                DIR *dr = opendir(".");
                
                if (dr == NULL)  // opendir returns NULL if couldn't open directory
                {
                    printf("Could not open current directory" );
                    return 0;
                }
                
                // Refer http://pubs.opengroup.org/onlinepubs/7990989775/xsh/readdir.html
                // for readdir()
                while ((de = readdir(dr)) != NULL){
                    if ( strstr(de->d_name,".spec.") != NULL&& strstr(de->d_name,filename)!= NULL)
                        printf("%s\n", de->d_name);
                }
                closedir(dr);
                return 0;
            }
                
            case 28:
            {
                sprintf(f1->fileList[f1->files++],"%s.vector", filename);
                if (  strstr(filename,c->name) != NULL){
                    printf(" cannot name inputs same as outputs\n");
                    printf("%st %s\n", filename, c->name);
                    exit(1);
                }
//                if ( (c->rt.printFlag/2 ) % 2 == 0 )
//                    c->rt.printFlag += 2;//vector
//
                return io;
            }
            case 29:
            {
                sprintf(f1->fileVectorOperator[f1->filesVectorOperator++],"%s.vector", filename);
                return io;
            }
            case 30:
            {
//                c->i.outputFlag = 1;
                return io;
            }
            case 32:
            {
                FILE * in = fopen(filename,"r");
                INT_TYPE i,ii=0;
                size_t ms = MAXSTRING;
                char input_line[MAXSTRING];
                char * mask = input_line;
                
                do{
                    getline(&mask, &ms, in);
                    if ( 2 == sscanf(input_line,"%lf,%lf", &c->i.shiftVector[ii][0],&c->i.shiftVector[ii][1]) )
                        ii++;
                }while ( ! feof(in ));
                fclose(in);
                
//                c->i.bodyFlag = 1;
//                return io;
                return io;

            }
            case 33:
            {
                double twistedD [] = {0., -0.5, -0.666667, -0.75, -0.8, -0.833333, -0.857143, -0.875,
                    -0.888889, -0.9, -0.909091, -0.916667, -0.923077, -0.928571,
                    -0.933333, -0.9375, -0.941176, -0.944444, -0.947368, -0.95,
                    -0.952381, -0.954545, -0.956522, -0.958333, -0.96, -0.961538,
                    -0.962963, -0.964286, -0.965517, -0.966667, -0.967742, -0.96875,
                    -0.969697, -0.970588, -0.971429, -0.972222, -0.972973, -0.973684,
                    -0.974359, -0.975, -0.97561, -0.97619, -0.976744, -0.977273,
                    -0.977778, -0.978261, -0.978723, -0.979167, -0.979592, -0.98,
                    -0.980392, -0.980769, -0.981132, -0.981481, -0.981818, -0.982143,
                    -0.982456, -0.982759, -0.983051, -0.983333, -0.983607, -0.983871,
                    -0.984127, -0.984375, -0.984615, -0.984848, -0.985075, -0.985294,
                    -0.985507, -0.985714, -0.985915, -0.986111, -0.986301, -0.986486,
                    -0.986667, -0.986842, -0.987013, -0.987179, -0.987342, -0.9875,
                    -0.987654, -0.987805, -0.987952, -0.988095, -0.988235, -0.988372,
                    -0.988506, -0.988636, -0.988764, -0.988889, -0.989011, -0.98913,
                    -0.989247, -0.989362, -0.989474, -0.989583, -0.989691, -0.989796,
                    -0.989899, -0.99, -0.990099};
                double twistedC [] = {1., 1.5, 1.66667, 1.75, 1.8, 1.83333, 1.85714, 1.875, 1.88889, 1.9,
                    1.90909, 1.91667, 1.92308, 1.92857, 1.93333, 1.9375, 1.94118,
                    1.94444, 1.94737, 1.95, 1.95238, 1.95455, 1.95652, 1.95833, 1.96,
                    1.96154, 1.96296, 1.96429, 1.96552, 1.96667, 1.96774, 1.96875,
                    1.9697, 1.97059, 1.97143, 1.97222, 1.97297, 1.97368, 1.97436, 1.975,
                    1.97561, 1.97619, 1.97674, 1.97727, 1.97778, 1.97826, 1.97872,
                    1.97917, 1.97959, 1.98, 1.98039, 1.98077, 1.98113, 1.98148, 1.98182,
                    1.98214, 1.98246, 1.98276, 1.98305, 1.98333, 1.98361, 1.98387,
                    1.98413, 1.98438, 1.98462, 1.98485, 1.98507, 1.98529, 1.98551,
                    1.98571, 1.98592, 1.98611, 1.9863, 1.98649, 1.98667, 1.98684,
                    1.98701, 1.98718, 1.98734, 1.9875, 1.98765, 1.9878, 1.98795, 1.9881,
                    1.98824, 1.98837, 1.98851, 1.98864, 1.98876, 1.98889, 1.98901,
                    1.98913, 1.98925, 1.98936, 1.98947, 1.98958, 1.98969, 1.9898, 1.9899,
                    1.99, 1.9901};

                INT_TYPE i;
                for ( i= 0; i < 100 ; i++){
                    c->i.shiftVector[i][0] = twistedD[i];
                    c->i.shiftVector[i][1] = twistedC[i];
                }
                return io;
            }
        }
        }
    }
    return 0;
}



//default includes write to outFileName = name.kappa

INT_TYPE readInput(struct calculation *c , struct input_label * f1, FILE * in){
    size_t ms = MAXSTRING;
    size_t read;
    INT_TYPE state,com;
    INT_TYPE run = 0;
    char input_line[MAXSTRING];
    char * mask = input_line;
    char line0[MAXSTRING];
    state = 0;
    INT_TYPE temp;
    do {
        read =  getline( &mask,&ms, in );
        if (! comment( input_line) && strlen(input_line) > 1  ){
            com = control( input_line);
            if ( strstr ( c->name , "print")){
                printf("%s", input_line);
            }
            
            if ( com == 1 && run == 0 ){
                run = 1;
                sscanf(input_line, "%s %s", line0 , c->name);
            }
            if ( com != 1 && run == 0){
                return 1;
            }else if ( com == 1 && run == 0 ){
                run = 1;
            }else if ( com == -1 && run ==1 ){
                run = 2;
            }
            
            if (run == 1 )
            {
                
                    
            if ( com != 1 && com != 0 ){
                
                if ( com > 0 ){
                    if ( state  )
                        return -1;
                    else
                        state = com;
                }
                else //com < 0
                {
                    if ( state ){
                        if ( state == - com  )
                            state = 0;
                        else
                            return -1;
                    } else{
                        return -1;
                    }
                }
                }
            else {
               if ( state == 2 ){
                   if (! getParam(c,f1,input_line)){
                       printf("%s", input_line);
                       return state;
                   }
                }
                else if ( state == 3 ){
                    if (getGeometry(c,input_line)) {
                        printf("%s", input_line);
                        return state;
                    }
                }
                else if ( state == 4 ){
                    if (modGeometry(c,input_line)) {
                        printf("%s", input_line);
                        return state;
                    }
                }
                else if ( state == 5 ){
                
                }
                else if ( state == 6 ){
                    temp =   getInputOutput(c,f1,input_line);
                    //printf("*%lld\n", temp);
                    if (!temp) {
                        printf("%s", input_line);
                        return state;
                    }
                }
       //         else if ( state == 7 ){
            //        if (getInitialGeneral(c,input_line)) return state;
            //    }
//                else if ( state == 8 ){
//                    if ( getBuildParameters(c,input_line)) return state;
//                }
                else if ( state == 9 ){
                    if ( intervalGeometry(c,input_line)) return state;
                }
                else if ( state == 12 ){
#if 0

                    if ( getCore(c,input_line)) return state;
#endif
                }
           //     else if ( state == 13 ){
            //        if ( generateRandom(c,input_line)) return state;
           //     }
//                else if ( state == 10 ){
//                    if ( defineNaturalBasis(c,input_line)) return state;
//                }
//                else if ( state == 11 ){
//                    if ( assignNaturalBasis(c,input_line)) return state;
//                }
            }
               

        }
    }
    }  while ( run != 2 );
    if ( run != 2 )
        return 1;
    else {
   // printf("exit\n");
    return 0;
    }
}

INT_TYPE initCalculation(struct calculation * c ){
    INT_TYPE g,space;
    c->i.level = 1e9;
    c->i.massElectron = 1.;
    c->i.massProton = 1836.15267245;
    c->i.massClampPair = 1836.15267245;
    c->rt.powDecompose = 2;
    c->i.shiftFlag = 0;
    c->i.complexType = 1;//real =1 , cmpl = 2
    c->i.RAMmax = 0;//Gb  needs updating
    c->i.springFlag = 0;
    c->i.magFlag = 0;
    c->i.M1 = 0;
    c->i.Na = 0;
 //   c->i.OCSBflag = 0;
//    for ( g = 0; g < nSAG*nSAG*nSAG ; g++)
//        c->i.cSA[g] = 0;
//#ifdef PARAMETER_PATH
//    FILE * same;
//    char filename[MAXSTRING];
//    sprintf(*&filename, "%s/%s.%s", PARAMETER_PATH, "default","alpha" );
//    same = fopen (filename, "r");
//    if (  readInput(c,same) )
//       exit(0);
//    fclose(same);
//#endif
    return 0;
}

INT_TYPE estSize ( INT_TYPE interval ){
    
    if ( interval == 0 ) {
        return 14;
    }
    if ( interval == 1 ) {
        return 22;
    }
    if ( interval == 2 ){
        return 50;
    }
    if ( interval == 3 ){
        return 134;
    }
    return 0;
}



INT_TYPE finalizeInit(struct calculation * c ){
    c->i.oneBody.num =estSize(3);
    c->i.twoBody.num =estSize(2);
    return 0;
}


INT_TYPE bootShell (INT_TYPE argc , char * argv[],struct calculation * c1, struct field *f){
#ifndef APPLE
    argc--;
    
    
    INT_TYPE broke;
    
    time_t t;
    /* Intializes random number generator */
    srand((unsigned) time(&t));
    
    INT_TYPE i,c,EV,EV2,ER;
    FILE * in = stdin;
    char str[MAXSTRING];
    *f = initField();
    *c1 = initCal();
    initCalculation(c1);

//
//    if ( argc > 0 && argc < MAXFIELD ){
//
//        for ( i = argc ; i > 0 ; i++){
//            broke = readInput(&c1,list[i],in );
//            if ( broke )
//                exit(1);
//        }
//        printf("erasing previous geometries");
//        c1.i.Na = 0;
//    }

    //    
    {
        broke = readInput(c1,&f->i,in);
        if ( broke )
            exit(1);
    }
    finalizeInit(c1);
    f->f.boot = fullMatrices;

#else
    *c1 =initCal();
    *f = initField();
#endif
#ifdef GSL_LIB
        gsl_set_error_handler_off ();
#endif
        
        return 1;
}
