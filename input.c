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

INT_TYPE add_atom( struct field *f1, INT_TYPE l, double x1,double y1,double z1 ){
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
    
    
    
    for ( alpha = 1 ; alpha <= c1->i.c.Na ; alpha++){
        for( l = 1 ; l <= 3 ; l++){
            cm[l] += c1->i.c.atoms[alpha].label.Z*c1->i.c.atoms[alpha].position[l];
            wcm[l] += c1->i.c.atoms[alpha].label.Z;
        }
    }
    
    cm[1] /= wcm[1];
    cm[2] /= wcm[2];
    cm[3] /= wcm[3];
    
    for ( a = 1 ; a <= c1->i.c.Na ;a++){
        for ( l = 1 ; l<= 3 ; l++){
            c1->i.c.atoms[a].position[l] -= cm[l];
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


INT_TYPE getParam ( struct calculation * c, const char * input_line ){
    INT_TYPE i,d,ivalue;
    char test_line [MAXSTRING];
    double value;
    INT_TYPE NINT_TYPE = 116;
    char *list_INT_TYPE []= {"#",
        "LOST1","maxCycle" , "spinor", "charge","fineStr",//5
        "process", "NB", "MB", "percentFull","general",//10
        "center","xTranslate","yTranslate","zTranslate","postCalc",//15
        "goK","goV","goC","goX","goS",//20
        "iGold","LOST", "LOST2" ,"core","canon",//25
        "pseudo","minDIIS","iCharge","weylet", "nWeylet",//30
        "mWeylet","helium","correlation","initRank","LOST3",//35
        "spinBlocks","LOST5", "LOST6", "maxLevelShift","diis",//40
        "minSPC","maxEV","inverseQuad","maxSPC","maxDIIS",//45
        "intDIIS", "trace","basisRank","LUMO","foundation",//50
        "Angstroms","train", "ceilFlag","iSymm","type",//55
        "zone","eZone","cycles","weightRank","printConvergence",//60
        "runFlag","LOST10","annulus","exclusion","runType",//65
        "SpinSqr","setRange","range","golden","bandStage",//70
        "matrix","jCycle","fermi","signature","filter",//75
        "minRank","skipBuild","printLevel","stack","lanes",//80
        "sectors","body","LOST100","rds1","rds2",//85
        "rds3","interactionOne","interactionTwo","oCycle","interactionZero",//90
        "breakBody","interval","RAM","monteCarlo","samples",//95
        "hartreeFock","basisStage","iterations","group","states",//100
        "length","side","lookBack","step","theory",
        "configuration","densityRank","densityBody","parallel","phase",
        "around","cmpl","clampStage","OCSB","decompose",
        "shiftNO"
    };
    INT_TYPE NDOUBLE = 71;
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
        "pairMass"
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
                    
                case 4 :
                    c->i.charge = ivalue;
                    return i;
                    
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
                    c->i.iCharge = ivalue;
                    return i;
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
                    c->i.heliumFlag = ivalue;
                    return i;
//                case 33:
//                    c->i.correlationFlag = ivalue;
//                    return i;
                case 34:
                    c->i.iRank = ivalue;
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
                    c->i.bRank = ivalue;
                    return i;
                case 49 :
              //      c->i.lumos = ivalue;
                    return i;
                case 50 :
                    c->i.qFloor = ivalue;
                    c->i.sectors = !(!ivalue);
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
                    c->i.irrep = ivalue;
                    return i;
                case 56 :
                //    c->i.zone = ivalue;
                    return i;
                case 57 :
                 //   c->i.eZone = ivalue;
                    return i;
                case 58 :
                    c->i.cycles = ivalue;
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

                    c->i.d /= (2*ivalue+1)/(2*c->i.epi+1);
                    c->i.epi = ivalue;
                    
                    
                    
                    
                    
                    
                    
                    return i;
                    
                case 68 :

                    c->i.epi = ivalue;
                    
                    
                    return i;

                    
                case 69 :
                    //c->i.runCharacter = ivalue;
                    return i;
                    
                    
                case 70:
                    //expand range by number given.
                    c->i.d *= pow( (2.* c->i.epi + 1.) /(2.*c->i.epi + 2*ivalue + 1),c->i.attack);
                    c->i.epi  += ivalue;
                    
                    

                    return i;

                case 71:
                    c->i.decomposeRankMatrix = ivalue;
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
                    c->i.filter = ivalue;
                    return i;
                    
                case 76:
            //        c->i.minRank = ivalue;
                    return i;

                case 77:
               //     c->i.overBuildFlag = ivalue;
                    return i;

                case 78 :
                    if ( c->rt.printFlag ){
                        
                        c->rt.printFlag = ivalue;

                        return i;
                    }
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
                    c->i.sectors = ivalue;
                    return i;
                    
                case 82:
                    if ( ivalue == 0 )
                        c->rt.body = nada;
                  else   if ( ivalue == 1 )
                        c->rt.body = one;
                    else if ( ivalue == 2 )
                        c->rt.body = two;
                    else if ( ivalue == 3 )
                        c->rt.body = three;
                    else if ( ivalue == 4 )
                        c->rt.body = four;
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
                        c->i.c.oneBody.func.interval = c->i.interval;
                        c->i.c.oneBody.func.fn = ivalue;
                        c->i.c.oneBody.func.param[0] = c->i.scalar;
                        c->i.c.oneBody.func.param[1] = c->i.turn;
                        c->i.c.oneBody.func.param[2] = c->i.param1;
                        c->i.c.oneBody.func.param[3] = c->i.param2;

                        return i;

                    
                case 88:
                        c->i.c.twoBody.func.interval = c->i.interval;
                        c->i.c.twoBody.func.fn = ivalue;
                        c->i.c.twoBody.func.param[0] = c->i.scalar;
                        c->i.c.twoBody.func.param[1] = c->i.turn;
                        c->i.c.twoBody.func.param[2] = c->i.param1;
                        c->i.c.twoBody.func.param[3] = c->i.param2;

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
                    c->rt.monteCarlo = ivalue;
                    return i;
                case 95:
                    c->rt.samples = ivalue;
                    return i;
                case 96:
//                    c->i.hartreeFockFlag = ivalue;
                    return i;
                case 97:
                    c->i.bRank += ivalue;
                    return i;
                case 98:
                    c->i.Iterations = ivalue;
                    return i;
                case 99:
            //        c->i.group = ivalue;
                    return i;
                case 100:
                    c->i.nTargets = ivalue;
                    c->i.heliumFlag = ivalue;
                    c->i.nStates = ivalue;
                    return i;
                case 101:
                    c->i.l2 = ivalue;
                    return i;
                case 102:
                    c->i.side = ivalue;
                    return i;
                case 103:
                    c->i.lookBack = ivalue;
                    return i;
                case 104:
                    c->i.cycleStep = ivalue;
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
                    c->i.dRank = ivalue;
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
                    printf("%f - %f in spacing %f\n", c->i.minClamp, c->i.maxClamp, c->i.D);
                    c->i.orgClamp = 0.5*(c->i.minClamp+c->i.maxClamp);
                    c->i.around = floor((c->i.orgClamp-c->i.minClamp)/c->i.D);
                    printf("grid %d\n", 2*c->i.around +1);
                    return i;
                    
                case 112:
                    c->i.complexType = ivalue;
                    return i;
                    
                case 113:
                    c->i.D *= pow( (2.* c->i.around + 1.) /(2.*c->i.around + 2*ivalue + 1),1.);
                    c->i.around  += ivalue;
                    return i;
                    
                case 114:
                    c->i.OCSBflag = ivalue;
                    return i;
                case 115:
                    c->rt.powDecompose = ivalue;
                    return i;
                case 116:
                    c->i.shiftFlag = 0;
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
                        c->i.d = value;
                        return d;
                    } else {
                        c->i.d = 1;
                        c->i.epi = 0;
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
                    c->i.epi = value;
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
                    c->i.epi *= c->i.d / value ;
                    c->i.d = value;
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
                    c->i.d = pow( 4. * pi /3. *c->rt.body ,0.33333333333333333333)* value / (2.*c->i.epi+1.);
                    //value = Rs

                    //4pi/3 Rs^3 = pow(d * N1 ,3.0)/f1->Ne
                    return d;
                case 46 :
                    c->i.shiftFlag = 1;
                    c->i.realPart = value;
                    return d;
                case 47 :
                 //   c->i.shift = value * pi * pi / c->i.d / c->i.d ;
                    return d;
                case 48 :
                    printf("Crystal Momentum -universal- \t%f\n", value);

                    value /= c->i.d*( 2.*c->i.epi+1.);
                    c->i.vectorMomentum = value;
                    return d;
                    
                case 49:
                    if ( value <= 0 )
                        return 0;
                    else{
                    //    c->i.jelliumSphereRadius = value;
                        c->i.potentialFlag = 1;
                        return d;
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
                      c->i.attack = value;
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
                    c->i.springConstant = 1./(c->i.d*value)/(c->i.d*value);
                    printf("spring %f (%f,%f)\n",c->i.springConstant,c->i.d,value);
                    return d;
                case 65:
                    c->i.D = value;
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
                add_atom ( &(c->i.c), l,x/a0,y/a0,z/a0 );
            else
                add_atom ( &(c->i.c), l,x,y,z );
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
    
    for ( alpha = 1; alpha <= c->i.c.Na ; alpha++)
    {
        x [0] = c->i.c.atoms[alpha].position[1];
        x [1] = c->i.c.atoms[alpha].position[2];
        x [2] = c->i.c.atoms[alpha].position[3];
        
        cblas_dgemv(CblasColMajor, CblasNoTrans,3,3,1.,R,3,x,1,0.,y,1);
        
        c->i.c.atoms[alpha].position[1] = y[0];
        c->i.c.atoms[alpha].position[2] = y[1];
        c->i.c.atoms[alpha].position[3] = y[2];

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
        if ( l > 0 && l <= c->i.c.Na ){
            if ( c->i.Angstroms )
            {
                c->i.c.atoms[l].position[1] += x/a0;
                c->i.c.atoms[l].position[2] += y/a0;
                c->i.c.atoms[l].position[3] += z/a0;
            } else {
                {
                    c->i.c.atoms[l].position[1] += x;
                    c->i.c.atoms[l].position[2] += y;
                    c->i.c.atoms[l].position[3] += z;
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
                add_atom ( &(c->i.c), l,x3/a0,y3/a0,z3/a0 );
            else
                add_atom ( &(c->i.c), l,x3,y3,z3 );

        }
        
        
        return 0;
    }
    else
        return 1;

}


INT_TYPE getInputOutput(struct calculation * c, const char * input_line ){
    INT_TYPE io;
    INT_TYPE Nio = 31;
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
        "body"
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
                if ( readInput( c, mid ) )
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
                sprintf(c->mem.fileList[c->mem.files++],"%s.vector", filename);
                if ( (c->rt.printFlag/2 ) % 2 == 0 )
                    c->rt.printFlag += 2;//vector

                return io;
            }
            case 29:
            {
                sprintf(c->mem.fileVectorOperator[c->mem.filesVectorOperator++],"%s.vector", filename);
                if ( c->i.vectorOperatorFlag % 2 == 0 )
                    c->i.vectorOperatorFlag += 1;//vectorOperator
                return io;
            }
            case 30:
            {
                c->i.outputFlag = 1;
                return io;
            }
            case 31:
            {
                c->i.bodyFlag = 1;
                return io;
            }

        }
        }
    }
    return 0;
}



//default includes write to outFileName = name.kappa

INT_TYPE readInput(struct calculation *c , FILE * in){
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
                   if (! getParam(c,input_line)){
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
                    temp =   getInputOutput(c,input_line);
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
    INT_TYPE g;
    c->i.level = 1e9;
    c->i.massElectron = 1.;
    c->i.massProton = 1836.15267245;
    c->i.massClampPair = 1836.15267245;
    c->rt.boot = fullMatrices;
    c->rt.powDecompose = 1;
    c->i.shiftFlag = 0;
    c->i.bodyFlag = 0;
    c->i.complexType = 1;//real =1 , cmpl = 2
    c->i.RAMmax = 0;//Gb  needs updating
    c->rt.printFlag = 0;
    c->i.potentialFlag = 0;
    c->i.vectorOperatorFlag = 0;
    c->i.springFlag = 0;
    c->i.outputFlag = 0;
    c->i.magFlag = 0;
    c->i.M1 = 0;
    c->i.c.Na = 0;
    c->mem.files = 0;
    c->mem.filesVectorOperator  = 0;
    c->i.OCSBflag = 0;
//    for ( g = 0; g < nSAG*nSAG*nSAG ; g++)
//        c->i.cSA[g] = 0;
#ifdef PARAMETER_PATH
    FILE * same;
    char filename[MAXSTRING];
    sprintf(*&filename, "%s/%s.%s", PARAMETER_PATH, "default","alpha" );
    same = fopen (filename, "r");
    if (  readInput(c,same) )
       exit(0);
    fclose(same);
#endif
    return 0;
}



INT_TYPE finalizeInit(struct calculation * c ){
    //count up atoms and electrons//
    INT_TYPE a,nc;
    nc = 0;
    
    for ( a = 1 ; a <= c->i.c.Na ; a++ ){
        {
            nc += c->i.c.atoms[a].label.Z;
        }
    }
    c->i.c.Ne = nc-c->i.charge;//set charges
    
    
    if ( c->i.c.oneBody.func.interval <= 1 )
        c->i.c.oneBody.num =23;
    else {
        printf("interval 1\n");
        exit(9);
    }
    if ( c->i.c.twoBody.func.interval <= 1 )
        c->i.c.twoBody.num =23;
    else {
        printf("interval 1\n");

        exit(9);
    }
    c->i.vectorOperatorFlag  =  countLinesFromFile(c,1);
    return 0;
}


struct calculation bootShell (INT_TYPE argc , char * argv[]){
#ifndef APPLE
    
    
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
#ifdef GSL_LIB
        gsl_set_error_handler_off ();
#endif
        
        return c1;
}
