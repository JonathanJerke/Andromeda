/**
 *  input.c
 *
 *
 *  Copyright 2021 Jonathan Jerke and Bill Poirier.
 *  Ongoing support for this program is coordinated through quantumgalaxies.org.
 *  We acknowledge the generous support of Texas Tech University,
 *  the Robert A. Welch Foundation, and the Army Research Office.
 *
 
 *   *   This file is part of Andromeda.
 
 *   *   Andromeda is free software: you can redistribute it and/or modify
 *   *   it under the terms of the GNU General Public License as published by
 *   *   the Free Software Foundation, either version 3 of the License.
 
 *   *   Andromeda is distributed in the hope that it will be useful,
 *   *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   *   GNU General Public License for more details.
 
 *   *   You should have received a copy of the GNU General Public License
 *   *   along with Andromeda.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "input.h"

inta resetA(   runTime *rt){
    inta i ;
    for ( i = 0; i < BLOCK_COUNT ; i++)
        rt->memBlock[i] = passBlock;
    return 1;
}

inta allowQ(   runTime *rt,   blockMemoryType a ){
    inta i,flag = 0;
    for ( i = 0; i < BLOCK_COUNT ; i++)
        flag += (rt->memBlock[i]==a);
    return (flag == 0) ;
}

inta blockA(   runTime *rt,   blockMemoryType a ){
    inta iii=0;
    while ( rt->memBlock[iii] != passBlock ){
        iii++;
        if ( iii == BLOCK_COUNT ){
            printf("add BLOCK_COUNT\n");
            exit(0);
        }
    }
    rt->memBlock[iii] = a;
    return 1;
}

inta add_atom(   input *f1, inta l, double x1,double y1,double z1 ){
    inta n = f1->Na;
    inta Z1 ;
    if ( n+1 == MAXATOM ){
        fprintf(stderr,"TOO MANY ATOMS");
        return -1;
    }
    f1->Na++;
    n++;
    f1->atoms[n].position[1] = x1;
    f1->atoms[n].position[2] = y1;
    f1->atoms[n].position[3] = z1;
    
    Z1 = ((l));
    
    f1->atoms[n].Z = (Z1);
    
    return 0;
};


inta centerMass (   calculation * c1){
    inta alpha;
    inta a;
    inta l;
    double cm[] = {0.,0.,0.,0.};
    double wcm[] = {0.,0.,0.,0.};
    
    
    
    for ( alpha = 1 ; alpha <= c1->i.Na ; alpha++){
        for( l = 1 ; l <= COMPONENT ; l++){
            cm[l] += c1->i.atoms[alpha].Z*c1->i.atoms[alpha].position[l];
            wcm[l] += c1->i.atoms[alpha].Z;
        }
    }
    
    cm[1] /= wcm[1];
    cm[2] /= wcm[2];
    cm[3] /= wcm[3];
    
    for ( a = 1 ; a <= c1->i.Na ;a++){
        for ( l = 1 ; l<= COMPONENT ; l++){
            c1->i.atoms[a].position[l] -= cm[l];
            //            printf("%f\n",cm[l]);
        }
    }
    
    return 0;
}

inta comment (const char * line ){
   return (NULL != strstr(line, "#"));
}

inta control ( const char * line ){
    inta m,s;
    inta Ncom = 8;
    char *list_control []= {"#",
        "Body",
        "Parameters",       "Geometry",
        "Translation" ,     "Basis",
        "InputOutput",      "Term" ,
        "Dim"
    };

    char *list_line [] = {".","*","?"};
    char line0 [MAXSTRING];
    char * line1 ;

        for (m = 1 ; m <= Ncom ; m++){
            
            for ( s = 0; s < 3 ; s++){
                
                sprintf(line0, "%s%s", list_line[s],list_control[m] );
                line1 = strstr( line, line0 );
                if ( line1 != NULL ){
                    if ( s == 2 ){
                        return 10000;
                    }
                    return (2*s-1)*m;
                }
            }
        }
 
        return 0;
    
}
inta getParam (   calculation * c,  input_label *f1, const char * input_line ){
    inta i,d,ivalue;
    char test_line [MAXSTRING];
    double value;

    
    inta NINT_TYPE = 32;
    char *list_INT_TYPE []= {"#",
        "lambda",
        "initRank",
        "canonRank",///changed name from basisRank
        "foundation",
        "Angstroms",
        "irrep",///changed name from type
        "filter",
        "lanes",
        "RAM",
        "canonStage",///changed name from basisStage
        "collect",
        "states",
        "configuration",
        "parallel",
        "phase",
        "cmpl",
        "shiftNO",
        "blockMemory",
        "blockReset",
        "switchGeometry",
        "Above",
        "maxCycle",
        "spam",
        "names",
        "eikons",
        "body",
        "iterations",
        "dynamic",
        "levelFoundation",
        "seed",
        "iocsb",
        "nocsb",
    };
    inta NDOUBLE = 8;
    char *list_DOUBLE []= {"#",
        "maxCondition",
        "tolerance",
        "relativeTolerance",
        "shift",
        "condition",
        "power",
        "threshold",
        "widthFoundation"
    };
    
    for ( i = 1 ; i <= NINT_TYPE ; i++){
        
        if ( strstr( input_line, list_INT_TYPE[i])!= NULL){
#ifdef BIT_INT
            sscanf(input_line,"%s %d", test_line,&ivalue);
#endif
            
#ifdef BIT_LONG
            sscanf(input_line,"%s %ld", test_line,&ivalue);
#endif

#ifdef MKL
            sscanf(input_line,"%s %lld", test_line,&ivalue);
#endif
            
            switch ( i ){
                case 1:
                    c->i.Lambda  = ivalue;
                    return i;
                case 2:
                    f1->iRank = ivalue;
                    return i;
                case 3 :
                    f1->canonRank = ivalue;
                    return i;
                case 4 :
                    f1->qFloor = ivalue;
                    return i;
                case 5 :
                    c->i.Angstroms = ivalue;
                    return i;
                case 6 :
                    f1->irrep = ivalue;
                    return i;
                case 7:
                    f1->filter = ivalue;
                    return i;
                case 8:
#ifdef OMP
                    c->i.omp = ivalue;
#endif
                    return i;
                case 9:
                    c->i.RAMmax = ivalue;
                    return i;
                case 10:
                    f1->canonRank += ivalue;
                    return i;
                case 11:
                    f1->collect = ivalue;
                    return i;
                case 12:
                    f1->nStates = ivalue;
                    return i;
                case 13:
                    if ( ivalue == 0 )
                        c->rt.calcType = electronicStuctureCalculation;
                    else
                    {
                        return 0;
                    }
                    return i;
                case 14:
#ifdef MKL
                    c->i.mkl = ivalue;
#endif
                    return i;
                case 15:
                    c->rt.phaseType = ivalue;
                    return i;
                case 16:
                    f1->cmpl = ivalue;
                    return i;
                case 17:
                    c->i.shiftFlag = ivalue;
                    return i;
                case 18:
                    blockA(&c->rt,ivalue);
                    return i;
                case 19:
                    resetA(&c->rt);
                    return i;
                case 20:
                    c->i.Na *= -1;
                    return i;
                case 21:
                    c->i.minIterationPrint = ivalue;
                    return i;
                case 22:
                    c->rt.MAX_CYCLE = ivalue ;
                    return i;
                case 23:
                    f1->OpIndex = ivalue ;
                    return i;
                case 24:
                    c->i.numNames = ivalue ;
                    return i;
                case 25:
                    c->i.numVectors = ivalue ;
                    return i;
                case 26:
                    ///for use on SA only
                    f1->body = ivalue;
                    return i;
                case 27:
                    f1->Iterations = ivalue;
                    return i;
                case 28:
                    c->rt.dynamic = ivalue;
                    return i;
                case 29:
                    c->i.SymmetrizedGaussianLevel = ivalue;
                    return i;
                case 30:
                    srand(ivalue);
                    return i;
                case 31:
                    c->i.iocsb = ivalue;
                    return i;
                case 32:
                    c->i.nocsb = ivalue;
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
                    c->rt.XCONDITION = pow(10., value);
                    return d;
                case 2 :
                    c->rt.TOLERANCE = pow(0.1, value);
                    return d;
                case 3 :
                    c->rt.relativeTOLERANCE = pow(10.,-value);
                    return d;
                case 4:
                {
              //      c->i.shiftVector[0] = 1.;
                //    c->i.shiftVector[1] = value;
                    c->i.shiftFlag = 1;
                }
                    return d;
                case 5:
                    c->rt.ALPHA = pow(0.1,value);
                    return d;
                case 6:
                {
            //        c->i.shiftVector[0] = 0;
              //      c->i.shiftVector[1] = value;
                    c->i.shiftFlag = 0;
                    return d;
                }
                case 7:
                {
                    c->rt.THRESHOLD = pow(0.1,value);
                    return d;
                }
                case 8:
                    c->i.SymmetrizedGaussianWidth = value;
                    return d;
            }

        }
        
    }
    return 0;
}

inta getGeometry(  calculation * c, const char * input_line ){
    
    double x,y,z;
    inta l,Z;
    
#ifdef BIT_INT
        if ( 4 == sscanf ( input_line, "%d %lf %lf %lf ", &l, &x,&y,&z ))
#endif
            
#ifdef BIT_LONG
        if ( 4 == sscanf ( input_line, "%ld %lf %lf %lf ", &l, &x,&y,&z ))
#endif

#ifdef MKL
        if ( 4 == sscanf ( input_line, "%lld %lf %lf %lf ", &l, &x,&y,&z ))
#endif
    {
        
        
        
        if ( l != 0 ){
            
            Z = ((l));
            
                
        
            
            if ( c->i.Angstroms )
                add_atom ( &(c->i), l,x/a0,y/a0,z/a0 );
            else
                add_atom ( &(c->i), l,x,y,z );
            }
            return 0;
        }
        else {
            
            double u[] = {x,y,z};
            
            rotateGeometry(c,u);
            
            
            return 0;

        }
        
        
        
        
     
    
        return 1;
    
}


inta rotateGeometry (  calculation * c, double * u  ){
    double mag = cblas_dnrm2(3,u,1);
    inta alpha;
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

inta modGeometry(  calculation * c, const char * input_line ){
    
    double x,y,z;
    inta l;
    inta input ;
   #ifdef BIT_INT
       input =  sscanf ( input_line, "%d %lf %lf %lf ", &l, &x,&y,&z );
   #endif
    
       #ifdef BIT_LONG
           input =  sscanf ( input_line, "%ld %lf %lf %lf ", &l, &x,&y,&z );
       #endif
       #ifdef MKL
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

inta getInputOutput(  calculation * c,  field * f, const char * input_line ){
      input_label * f1 = &f->i;
    inta io;
    inta Nio = 6;
    char test_line [MAXSTRING];
    char *list_IO[] = {"#",
        "read" ,
        "vector",
        "operator",
        "control",
        "path",
        "sum"
    };
    char filename[MAXSTRING];
    FILE * mid;
    
    for( io = 1 ; io <= Nio ; io++){
        if ( strstr( input_line, list_IO [io])!=NULL){
            sscanf(input_line,"%s %s", test_line,  filename);
            
//            if ( strstr(filename , "*" ) != NULL )
//                sprintf(filename,"%s", c->name);
        switch ( io ){
            
            case 1:
                
                mid = NULL;
                mid = fopen ( filename, "r");
                if( mid == NULL )
                    return 0;
                strcpy(test_line,c->name);
                if ( readInput( c, f,mid ) )
                    return 0;
                strcpy(c->name, test_line);//does not inheret name of run
                fclose( mid);
                return io;
                
           
            case 2:
            {
                sprintf(f1->fileList[f1->files++],"%s.vector", filename);
                if (  strstr(filename,c->name) != NULL){
                    printf(" cannot name inputs same as outputs\n");
                    printf("%st %s\n", filename, c->name);
                    exit(1);
                }
                return io;
            }
            case 3:
            {
                sprintf(f1->fileVectorOperator[f1->filesVectorOperator++],"%s.vector", filename);
                return io;
            }
            case 4:
                mid = NULL;
            {
                char file2[SUPERMAXSTRING];
                sprintf(file2,"%s/control/%s/%s",getenv("LAUNCH"),c->i.controlPath, filename);
                mid = fopen ( file2, "r");
            }
                if( mid == NULL )
                    return 0;
                strcpy(test_line,c->name);
                if ( readInput( c, f,mid ) )
                    return 0;
                strcpy(c->name, test_line);//does not inheret name of run
                fclose( mid);
                return io;
            case 5:
                strcpy(c->i.controlPath,filename);
                return io;
            case 6:
                {
                    sprintf(f1->matrixList[f1->matrices++],"%s", filename);
                    if (  strstr(filename,c->name) != NULL){
                        printf(" cannot name inputs same as outputs\n");
                        printf("%st %s\n", filename, c->name);
                        exit(1);
                    }
                    c->i.build = 0;
                    return io;
                }
        }
    }
    }
    return 0;
}




inta getTermDefinitions(  calculation * c, const char * input_line ){
    static char filename[SUPERMAXSTRING];
    static inta atom = 1;
    static inta embed = 0;
    static inta flagScalar = 0;
    static blockType bl = id0;
    static inta act = 1;
    static inta newTerm = 1;
    static inta invert = 0;
    static inta Interval = 1;
    static inta momentumInterval = 217;
    static inta contr = 2;
    static double adjustOne = 1.;
    static double scalar = 1;
    static double turn = 1;
    static double param1 = 1;
    static double param2 = 1;
    static double mBeta = 0;
    static double xBeta = 10;
    static inta bra = 0;
    static inta ket = 0;
    static inta particle1 = 1;
    static inta particle2 = 1;
    static inta funcType = 3;
    char test_line [MAXSTRING];
    inta i,d,ivalue,space;
    metric_label metric;
    double value;    
    inta io;
    inta Nio = 11;
    char *list_IO[] = {"#",
        "constant","linear","spring",
        "deriv","kinetic","clamp",
        "element","oneBody","twoBody",
        "diagonal",
        
        "filename"
    };
    char input[MAXSTRING];
    
    for( io = 1 ; io <= Nio ; io++){
        if ( strstr( input_line, list_IO [io])!=NULL){
                sscanf(input_line,"%s %s", test_line,  input);
                if ( io <= 10 ){
                    c->i.terms[c->i.termNumber].type = io;
                    c->i.terms[c->i.termNumber].act = act;
                    c->i.terms[c->i.termNumber].bl = bl;
                c->i.terms[c->i.termNumber].embed = embed;
                    c->i.terms[c->i.termNumber].scalar = scalar;
                    c->i.terms[c->i.termNumber].invert = invert;
                    strcpy(c->i.terms[c->i.termNumber].desc,input);
                    c->i.terms[c->i.termNumber].headFlag = newTerm;
                    c->i.terms[c->i.termNumber].adjustOne = adjustOne;
                c->i.terms[c->i.termNumber].bra = bra;
                c->i.terms[c->i.termNumber].ket = ket;
                sprintf( c->i.terms[c->i.termNumber].filename,"%s",filename);
                    strcpy(c->i.terms[c->i.termNumber].filename,filename);
                c->i.terms[c->i.termNumber].atom     = atom;
                c->i.terms[c->i.termNumber].label[0] = particle1;
                c->i.terms[c->i.termNumber].label[1] = particle2;

                    newTerm = 0;
                metric.fn.fn = funcType;
                metric.fn.param[0] = 1.;
                metric.fn.param[1] = turn;
                metric.fn.param[2] = param1;
                metric.fn.param[3] = param2;
                metric.fn.interval = Interval;
                metric.fn.momentumInterval = momentumInterval;
                metric.fn.contr    = contr;
                metric.beta[0]     = mBeta;
                //add gaussian!
                if ( !flagScalar ){
                    if ( xBeta > 0 ){
                        metric.beta[1] = xBeta;
                        metric.metric = interval;
                    } else {
                        metric.metric = semiIndefinite;
                    }
                }else {
                    if ( xBeta > 0 ){
                        metric.beta[1] = xBeta;
                        metric.metric = pureInterval;
                    } else {
                        metric.metric = pureSemiIndefinite;
                    }
                }
                if ( Interval == 1 ){
                    //only mBeta matters
                    metric.metric = dirac;
                }
                
                for ( space = 0; space < SPACE ; space++){
                    metric.pow[space] = 0;
                    metric.deriv[space] = 0;
                }
                if ( io >= 8 ){
                    c->i.terms[c->i.termNumber].mu = metric;
                }
                c->i.termNumber++;
                return io;
                    
                } else switch(io){
                    case 11:
                        sprintf(filename, "%s.vector",input);
                        return io;
                }
            }
                
        }
        
    inta NINT_TYPE = 19;
    char *list_INT_TYPE []= {"#",
        "invert","block","act","newTerm","buffer",
        "interval","offDiagonals","funcType","atom","axis",
        "bra","ket","flags","embed","revise",
        "reset","axis1","axis2","momentumInterval"
    };
    inta NDOUBLE = 10;
    char *list_DOUBLE []= {"#",
        "placeholder","scalar","funcTurn","funcParam1","funcParam2",
        "adjustOne","minBeta","maxBeta","infBeta","soloBeta"
    };
    
    for ( i = 1 ; i <= NINT_TYPE ; i++){
        if ( strstr( input_line, list_INT_TYPE[i])!= NULL){
#ifdef BIT_INT
                sscanf(input_line,"%s %d", test_line,&ivalue);
#endif
                
#ifdef BIT_LONG
                sscanf(input_line,"%s %ld", test_line,&ivalue);
#endif

#ifdef MKL
                sscanf(input_line,"%s %lld", test_line,&ivalue);
#endif
                
                switch ( i ){
                    case 1:
                        invert = ivalue;
                        return i;
                    case 2:
                        bl = ivalue;
                        return i;
                    case 3:
                        act = ivalue ;
                        return i;
                    case 4:
                        newTerm = ivalue ;
                        return i;
                    case 5:
                        return i;
                    case 6 :
                        Interval = ivalue;
                        return i;
                    case 7:
                        contr = ivalue;
                        return i;
                    case 8:
                        funcType = ivalue;
                        return i;
                    case 9:
                        atom = ivalue;
                        return i;
                    case 10:
                        particle1 = ivalue;
                        particle2 = ivalue;
                        return i;
                    case 11:
                        bra = ivalue;
                        return i;
                    case 12:
                        ket = ivalue;
                        return i;
                    case 13:
                        flagScalar = ivalue;
                        return i;
                    case 14:
                        embed = ivalue;
                        return i;
                    case 15:
                        c->i.terms[ivalue].scalar *= scalar;
                        return i;
                    case 16:
                        c->i.termNumber = 0;
                        return i;
                    case 17:
                        particle1 = ivalue;
                        return i;
                    case 18:
                        particle2 = ivalue;
                        return i;
                    case 19:
                        momentumInterval = ivalue;
                        return i;
                }
//
            }
        }
    
    
        for ( d = 1; d <= NDOUBLE ; d++){
            
            sprintf(test_line ,"%s", list_DOUBLE[d]);
            if ( strstr( input_line, test_line)!= NULL){
                sscanf(input_line,"%s %lf", test_line, &value);
                
                switch ( d ){
                    case 1:
                       // funcScalar = value;
                      //  return d;
                    case 2:
                        scalar = value;
                        return d;
                    case 3:
                        turn = value;
                        return d;
                    case 4:
                        param1 = value;
                        return d;
                    case 5:
                        param2 = value;
                        return d;
                    case 6:
                        adjustOne = value;
                        return d;
                    case 7:
                        mBeta = value;
                        return d;
                    case 8:
                        xBeta = value;
                        return d;
                    case 9:
                        xBeta = -1;
                        return d;
                    case 10:
                        mBeta = value;
                        Interval = 1;
                        funcType = 8;
                        return d;
                }
            }
        }
    return 0;
};


//inta getWaveFunction(  calculation * c, const char * input_line ){
//    static   Parti parti ;
//    static inta component=0;
//    static inta particle = 1;
//    static inta index = 0;//0..
//    char test_line [MAXSTRING];
//    inta i,d,ivalue;
//    double value;
//    inta io;
//    inta Nio = 2;
//    char *list_IO[] = {"#",
//        "out","user"
//    };
//
//    char input[MAXSTRING];
//    for( io = 1 ; io <= Nio ; io++){
//        if ( strstr( input_line, list_IO [io])!=NULL){
//            switch(io){
//                case 1://out-->  join particles
//                    c->basis[index].tion[particle] = parti;
//                    c->basis[index].mask = 1;
//                    return io;
//                case 2:
//                    sprintf(parti.cle[component].file ,"%s", input);
//                    if (  strstr(input,c->name) != NULL){
//                        printf(" cannot name inputs same as outputs\n");
//                        printf("%st %s\n", filename, c->name);
//                        exit(1);
//                    }
//                    return io;
//            }
//        }
//    }
//    inta NINT_TYPE = 6;
//    char *list_INT_TYPE []= {"#",
//        "particle","component","kinetic",
//        "sho","user","index"
//        };
//    inta NDOUBLE = 1;
//    char *list_DOUBLE []= {"#",
//        "x0"
//        };
//
//    for ( i = 1 ; i <= NINT_TYPE ; i++){
//
//        if ( strstr( input_line, list_INT_TYPE[i])!= NULL){
//    #ifdef BIT_INT
//            sscanf(input_line,"%s %d", test_line,&ivalue);
//    #endif
//
//    #ifdef BIT_LONG
//            sscanf(input_line,"%s %ld", test_line,&ivalue);
//    #endif
//
//    #ifdef MKL
//            sscanf(input_line,"%s %lld", test_line,&ivalue);
//    #endif
//
//            switch ( i ){
//                case 1:
//                    particle  = ivalue;
//                    return i;
//                case 2:
//                    component = ivalue;
//                    return i;
//                case 3:
//                    parti.cle[component].funcType = 1;//kinetic
//                    parti.cle[component].state = ivalue;
//                    return i;
//                case 4:
//                    parti.cle[component].funcType = 2;//sho
//                    parti.cle[component].state = ivalue;
//                    return i;
//                case 5:
//                    parti.cle[component].funcType = 0;//user
//                    parti.cle[component].state = ivalue;
//                    return i;
//                case 6:
//                    index = ivalue;
//            }
//        }
//    }
//
//
//    for ( d = 1; d <= NDOUBLE ; d++){
//
//        sprintf(test_line ,"%s", list_DOUBLE[d]);
//        if ( strstr( input_line, test_line)!= NULL){
//            sscanf(input_line,"%s %lf", test_line, &value);
//            switch ( d ){
//                case 1:
//                    parti.cle[component].x0 = value;
//            }
//        }
//    }
//    return 0;
//}
    
inta getDimensionalDefinitions(struct calculation * c,struct field * f, const char * input_line ){
    static inta count1Basis = 1;
    static inta count1Inc = 1;
    static   bodyType body = nada;
    static inta space = 0;
    static inta label = 1;
    static double attack1 = 0.5;
    static double attack2 = 0.5;
    static double attack3 = 0.5;
    static double attack4 = 0.5;

    static   basisElementType basisType = SincBasisElement;
    static double lattice1 = 1.0;
    static double lattice2 = 1.0;
    static double lattice3 = 1.0;
    static double lattice4 = 1.0;

    static double origin1 = 0.;
    static double origin2 = 0.;
    static double origin3 = 0.;
    static double origin4 = 0.;

    static double anchor1 = 0.5;
    static double anchor2 = 0.5;
    static double anchor3 = 0.5;
    static double anchor4 = 0.5;

        char test_line [MAXSTRING];
        inta i,d,ivalue,dim,components = 0;
        double value;
      bodyType body1;
            inta io;
            inta Nio = 3;
            char *list_IO[] = {"#",
                "oneComponent",
                "twoComponent",
                "threeComponent"
            };
            char filename[MAXSTRING];
            
            for( io = 1 ; io <= Nio ; io++){
                if ( strstr( input_line, list_IO [io])!=NULL){
                    sscanf(input_line,"%s %s", test_line,  filename);
                    
                    components = 0;
                    
                    switch ( io ) {
                        case 1:
                            components = 1;
                            ///  allocate 1 space.
                            break;
                        case 2:
                            components = 2;
                            break;
                            ///  allocate 2 space.
                        case 3:
                            components = 3;
                            break;
                            ///  allocate 3 space.
                    }
                    for ( dim = c->i.dimNumber ; dim < c->i.dimNumber + components ; dim++){
                        f->f.canon[dim].basis = basisType;
                        f->f.canon[dim].body = body;
                        f->f.canon[dim].count1Basis = count1Basis;
                        f->f.canon[dim].count1Inc = count1Inc;
                    //    f->f.canon[dim].component = components;
                        f->f.canon[dim].space = dim - space;
                        f->f.canon[dim].particle[one].anchor = anchor1;
                        f->f.canon[dim].particle[two].anchor = anchor2;
                        f->f.canon[dim].particle[three].anchor = anchor3;
                        f->f.canon[dim].particle[four].anchor = anchor4;

                        f->f.canon[dim].particle[one].lattice = lattice1;
                        f->f.canon[dim].particle[two].lattice = lattice2;
                        f->f.canon[dim].particle[three].lattice = lattice3;
                        f->f.canon[dim].particle[four].lattice = lattice4;

                        f->f.canon[dim].label = label;
                        f->f.canon[dim].particle[one].attack = attack1;
                        f->f.canon[dim].particle[two].attack = attack2;
                        f->f.canon[dim].particle[three].attack = attack3;
                        f->f.canon[dim].particle[four].attack = attack4;

                        f->f.canon[dim].particle[one].origin =   lattice1*floor(origin1/lattice1)-lattice1*floor(count1Basis*anchor1);
                        f->f.canon[dim].particle[two].origin =   lattice2*floor(origin2/lattice2)-lattice2*floor(count1Basis*anchor2);
                        f->f.canon[dim].particle[three].origin = lattice3*floor(origin3/lattice3)-lattice3*floor(count1Basis*anchor3);
                        f->f.canon[dim].particle[four].origin  = lattice4*floor(origin4/lattice4)-lattice4*floor(count1Basis*anchor4);

                        if ( dim >= SPACE ){
                            printf("add space\n");
                            return 0;
                        }
                    }
                    if ( components ){
                        label++;
                        c->i.dimNumber += components;
                    }
                    return components;
                }
                    
            }
            
        inta NINT_TYPE = 6;
            char *list_INT_TYPE []= {"#",
                "count1Basis","body","basis",
                "count1Stage","revise","count1Inc"
            };
            inta NDOUBLE = 20;
            char *list_DOUBLE []= {"#",
               "latticeAll","latticeOne","latticeTwo","latticeThree"
                ,"attackAll","attackOne","attackTwo","attackThree"
                ,"originAll","originOne","originTwo","originThree"
                ,"anchorAll","anchorOne","anchorTwo","anchorThree"
                ,"latticeFour","attackFour","originFour","anchorFour"
            };
            
            for ( i = 1 ; i <= NINT_TYPE ; i++){
                
                if ( strstr( input_line, list_INT_TYPE[i])!= NULL){
        #ifdef BIT_INT
                    sscanf(input_line,"%s %d", test_line,&ivalue);
        #endif
                    
        #ifdef BIT_LONG
                    sscanf(input_line,"%s %ld", test_line,&ivalue);
        #endif

        #ifdef MKL
                    sscanf(input_line,"%s %lld", test_line,&ivalue);
        #endif
                    
                    switch ( i ){
                        case 1:
                            count1Basis = ivalue;
                            return i;
                        case 2:
                             if ( ivalue == 0 )
                                  body = nada;
                             else if ( ivalue == 1 )
                                  body = one;
                             else if ( ivalue == 2 )
                                  body = two;
                             else if ( ivalue == 3 )
                                  body = three;
                             else if ( ivalue == 4 )
                                  body = four;
                            return i;
                        case 3:
                                 if ( ivalue == 0 )
                                      basisType  = nullBasisElement;
                                 else if ( ivalue == 1 )
                                      basisType  = SincBasisElement;
                                 else if ( ivalue == 2 )
                                      basisType  = PeriodicSincBasisElement;
                                 else if ( ivalue == 3 )
                                      basisType  = GaussianBasisElement;
                                 else if ( ivalue == 4 )
                                      basisType  = DiracDeltaElement;
                                 else if ( ivalue == 5 )
                                      basisType  = StateBasisElement;
                            return i;
                        case 4:
                            for ( dim = 0; dim < SPACE ; dim++)
                                if ( f->f.canon[dim].body != nada)
                                    if ( f->f.canon[dim].basis == SincBasisElement || f->f.canon[dim].basis == PeriodicSincBasisElement){
                                        ///Use count1Inc to update step increments
                                        ///this has been done to reduce memory in basis-transfers
                                        if (ivalue > 1 )
                                            ivalue = 1;
                                        
                                        
                                        for (body1 = one ; body1 <= f->f.canon[dim].body; body1++){
                                            f->f.canon[dim].particle[body1].origin += +f->f.canon[dim].particle[body1].lattice*floor(f->f.canon[dim].count1Basis*f->f.canon[dim].particle[body1].anchor);
                                            ///remove previous origin

                                            f->f.canon[dim].particle[body1].lattice *= pow( (f->f.canon[dim].count1Basis-1)*1./(f->f.canon[dim].count1Basis + ivalue*f->f.canon[dim].count1Inc -1),f->f.canon[dim].particle[body1].attack );
                                            ///attack lattice
                                            
                                            f->f.canon[dim].particle[body1].origin -= +f->f.canon[dim].particle[body1].lattice*floor((f->f.canon[dim].count1Basis+ivalue*f->f.canon[dim].count1Inc)*f->f.canon[dim].particle[body1].anchor);
                                            ///reintroduce previous origin

                                        }

                                        f->f.canon[dim].count1Basis += ivalue*f->f.canon[dim].count1Inc;

                                    }
                            return i;
                        case 5:
                            for ( dim = 0 ; dim < space ; dim++)
                                if ( ivalue == f->f.canon[dim].label ){
                                    f->f.canon[dim].particle[one].origin   += lattice1*floor(origin1/lattice1);
                                    f->f.canon[dim].particle[two].origin   += lattice2*floor(origin2/lattice2);
                                    f->f.canon[dim].particle[three].origin += lattice3*floor(origin3/lattice3);
                                    f->f.canon[dim].particle[four].origin  += lattice4*floor(origin4/lattice4);

                                    f->f.canon[dim].count1Inc = count1Inc;
                                    f->f.canon[dim].particle[one].attack = attack1;
                                    f->f.canon[dim].particle[two].attack = attack2;
                                    f->f.canon[dim].particle[three].attack = attack3;
                                    f->f.canon[dim].particle[four].attack = attack4;

                            }
                            return i;
                        case 6:
                            count1Inc = ivalue;
                            return i;
                    }
                }
            }
        
        
            for ( d = 1; d <= NDOUBLE ; d++){
                
                sprintf(test_line ,"%s", list_DOUBLE[d]);
                if ( strstr( input_line, test_line)!= NULL){
                    sscanf(input_line,"%s %lf", test_line, &value);
                    switch ( d ){
                        case 1:
                            lattice1 = value;
                            lattice2 = value;
                            lattice3 = value;
                            lattice4 = value;
                            return d;
                        case 2:
                            lattice1 = value;
                            return d;
                        case 3:
                            lattice2 = value;
                            return d;
                        case 4:
                            lattice3 = value;
                            return d;
                        case 5:
                            attack1 = value;
                            attack2 = value;
                            attack3 = value;
                            attack4 = value;
                            return d;
                        case 6:
                            attack1 = value;
                            return d;
                        case 7:
                            attack2 = value;
                            return d;
                        case 8:
                            attack3 = value;
                            return d;
                        case 9:
                            origin1 = value;
                            origin2 = value;
                            origin3 = value;
                            origin4 = value;
                            return d;
                        case 10:
                            origin1 = value;
                            return d;
                        case 11:
                            origin2 = value;
                            return d;
                        case 12:
                            origin3 = value;
                            return d;
                        case 13:
                            anchor1 = value;
                            anchor2 = value;
                            anchor3 = value;
                            anchor4 = value;
                            return d;
                        case 14:
                            anchor1 = value;
                            return d;
                        case 15:
                            anchor2 = value;
                            return d;
                        case 16:
                            anchor3 = value;
                            return d;
                        case 17:
                            lattice4 = value;
                            return d;
                        case 18:
                            attack4 = value;
                            return d;
                        case 19:
                            origin4 = value;
                            return d;
                        case 20:
                            anchor4 = value;
                            return d;
                    }
                }
            }
        return 0;
    }
    

inta readInput(struct calculation *c , struct field * f, FILE * in){
    struct input_label * f1 = &f->i;
    int ms = MAXSTRING;
    inta state,com;
    inta run = 0;
    char input_line[MAXSTRING];
    char line0[MAXSTRING];
    state = 0;
    inta temp;
    do {
        fgets( input_line , ms, in );
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
                
                    
            if ( com != 1 && com != 0 && com != 10000){
                if ( com == 10000 ){
                    if ( state == 10000 )
                        state = 0;
                    else if ( state == 0 )
                        state = 10000;
                    else
                        return 0;
                }else if ( com > 0 ){
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
                    temp =   getInputOutput(c,f,input_line);
                    //printf("*%lld\n", temp);
                    if (!temp) {
                        printf("%s", input_line);
                        return state;
                    }
                }
                else if ( state == 7 ){
                  if (!getTermDefinitions(c,input_line)) return state;
                }
                else if ( state == 8 ){
                    if (! getDimensionalDefinitions(c,f,input_line)) return state;
                }
            }
               

        }
    }
    }  while ( run != 2 );
    if ( run != 2 )
        return 1;
    else {
        return 0;
    }
}

inta initCalculation(struct calculation * c ){
    inta i;
    c->i.shiftFlag = 0;
    c->i.RAMmax = 0;//Gb  needs updating
    c->i.Na = 0;
    for ( i = 0 ; i < BLOCK_COUNT ; i++)
        c->rt.memBlock[i] = passBlock;

  //  c->i.shiftVector[0] = 0;
  //  c->i.shiftVector[1] = 1;

    return 0;
}

inta finalizeInit(struct calculation * c ){
    return 0;
}


inta bootShell (inta argc , char * argv[], calculation * c1,  field *f){
#ifndef APPLE
    FILE * in ;

    if ( argc >= 1 ){
        in = fopen(argv[0],"r");
    } else {
        in = stdin;
    }
    
    inta broke;
    
    
    inta i,c,EV,EV2,ER;
    char str[MAXSTRING];
    initField(f);
    initCal(c1);
    initCalculation(c1);
    {
        broke = readInput(c1,f,in);
        if ( broke )
            exit(1);
    }
    finalizeInit(c1);
    f->f.boot = fullMatrices;
    if ( argc >= 1 ){
        fclose(in);
    }
#else
    initCal(c1);
    initField(f);
#endif
#ifdef GSL_LIB
        gsl_set_error_handler_off ();
#endif
        
        return 1;
}


inta readShell (inta argc , char * argv[], calculation * c1,  field *f){
    FILE * in ;
    inta broke;
    if ( argc >= 1 ){
        in = fopen(argv[0],"r");
    } else {
        in = stdin;
    }
    
    {
        broke = readInput(c1,f,in);
        if ( broke )
            exit(1);
    }
    
    fclose(in);
    return 1;
}
