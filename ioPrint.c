/*
 *  ioPrint.c
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

#include "ioPrint.h"

INT_TYPE printVector (struct calculation *c, char * name,char * vectorName,  INT_TYPE iv, INT_TYPE irrep,DCOMPLEX * vector){
    if ( vector == NULL )
        return 1;
    struct field * f1 = &c->i.c;
    INT_TYPE iii;
    char str [MAXSTRING];
    FILE * outf ;
    if ( iv == -1 ){
        sprintf(str, "%s.vector",vectorName);
        outf = fopen (str,"w");
        fclose(outf);
        return -1;
    }
    
    
    sprintf(str, "%s.vector",vectorName);
    outf = fopen (str,"a");
    for ( iii = iv; iii <= iv  ; iii++)
    {
        if ( cabs(*vector)> f1->mem1->rt->TARGET){
            if( spins(f1, eigenVectors )== 2)
                fprintf(outf, "\"%s\",%d,%15.15f,%15.15f\n",name, iii+1,creal(*vector),cimag(*vector));
            else
                fprintf(outf, "\"%s\",%d,%15.15f\n", name, iii+1,creal(*vector));
        }
    }
    fclose(outf);
    return 0;
}

INT_TYPE print(struct calculation *c , INT_TYPE reset,INT_TYPE mv, INT_TYPE lv,enum division eigenVectors){
    INT_TYPE irrep;
    struct field * f1 = & c->i.c;
    INT_TYPE iii,jjj=1,cmpl;
    char str [MAXSTRING];
    DCOMPLEX one = 1.;
    if ( reset ) {
        FILE * outf ;
        sprintf(str, "%s.vector",c->name);
        outf = fopen (str,"w");
        fclose(outf);
    }
        for ( iii = mv; iii < lv  ; iii++)
          //  if( (! c->i.irrep || f1->sinc.tulip[eigenVectors+iii].value.symmetry  == irrep)&& irrep == c->i.irrep)
        {
                tEdges(c , eigenVectors+iii);
            irrep = f1->sinc.tulip[eigenVectors+iii].value.symmetry;
                printf("State%d:%d:,%d ,%1.15f, %d, %d , %1.1f,%1.15f\n", iii+1, c->i.epi*2+1,iii+1,f1->sinc.tulip[eigenVectors+iii].value.value,bodies(f1,eigenVectors+iii),irrep, deg(f1, irrep),f1->sinc.tulip[eigenVectors+iii].value.value2);
                
            //    if ( (c->i.outputFlag) % 2 == 1){
                    printVector(c, c->name,c->name, iii,irrep, &one);
                    for ( cmpl = 0 ; cmpl < spins(f1, eigenVectors+iii) ; cmpl++)
                    {
#ifndef APPLE
                        tFilename(c->name,iii+1,bodies(f1, eigenVectors+iii) ,irrep, cmpl,str);
                        
                        FILE * out = fopen ( str,"w" );
                        if ( out != NULL ){
                            outputFormat(&c->i.c, out, eigenVectors+iii,cmpl  );
                            fclose(out);
                        }
#endif
                    }
               // }
            
            }

    fflush(stdout);
    return 0;
}


INT_TYPE ioStoreMatrix(struct field * f1, enum division op, INT_TYPE spin, char * filename, INT_TYPE ioIn ){
    INT_TYPE matchFlag = 0,tempFlag=1,space;
    //check if previous is acceptable...
    //0 genus
    //2 ranks
    //4 length
    //7 second length
    // 5 d
    // 6 D
    if (   access( filename, F_OK ) != -1 ){

        if (  inputFormat(f1, filename, nullName, 0) == 2 ){
            for ( space = 0;space < SPACE ; space++)
                if ( f1->sinc.rose[space].body != nada ){
                if ( f1->sinc.tulip[op].space[space].body != inputFormat(f1, filename, nullName, 100+space/COMPONENT)){
                  //  printf("body");
                    tempFlag = 0;
                }
                if ( f1->sinc.rose[space].count1Basis != inputFormat(f1, filename, nullName, 200+space/COMPONENT)){
                 //   printf("count");
                 //   printf("%d--%d != %d \n",space,f1->sinc.rose[space].count1Basis ,inputFormat(f1, filename, nullName, 200+space/COMPONENT) );
                    tempFlag = 0;
                }

            }
            if ( tempFlag && part(f1, op ) >= inputFormat(f1, filename, nullName, 2)){
                        matchFlag = 1;
                
            }else {
             //   printf("ranks");
            }
        }
        else
            return 0;
    }
    
    
    {
        if ( ioIn == 0 ){
            if ( ! matchFlag){
                //print out.
                FILE* out = fopen(filename,"w");
                outputFormat(f1, out, op, spin);
                fclose(out);
            }
            return 1;
        }else if ( ioIn ==1 && matchFlag ){
                inputFormat(f1, filename, op, 1);
                return 1;
        }
    }

    
    return 0;
}

void outputFormat(struct field  * f1, FILE * out, enum division output ,INT_TYPE spin){
    INT_TYPE  parts,p1,flag2,flag3,flag4,r,l,space, M[SPACE];
    length(f1, output, M);
    if ( header (f1, output ) != Cube )
    {
        printf("outputFormat: WARNING,Non Cubic output\n");
    }

    fprintf(out,"component = %d\n", COMPONENT);
    enum genus g = species(f1, output);
    if (g == 3 )
        g = 1;
    
    fprintf(out,"genus = %d \n", g );
    fprintf(out,"header = %d\n", header(f1, output));
    fprintf(out,"canonRank = %d\n",CanonicalRank(f1, output,spin ) );
    fprintf(out,"spin = %d\n", spin);
    fprintf(out,"symmetry = %d\n", f1->sinc.tulip[output].value.symmetry);
    parts=0;
    while ( parts*COMPONENT < SPACE )
        parts++;
    fprintf(out,"particles = %d\n",parts);
    for ( p1 = 0 ; p1 < parts;p1++){
        if ( species(f1, output)  == vector )
            fprintf(out,"%cbody = %d\n", 'A'+p1,f1->sinc.rose[COMPONENT * p1].body);
        else if ( species(f1, output)  == matrix || species(f1,output) == outerVector )
            fprintf(out,"%cbody = %d\n", 'A'+p1,f1->sinc.tulip[output].space[COMPONENT * p1].body);

            
        fprintf(out,"%ccount1Basis = %d\n", 'A'+p1, f1->sinc.rose[COMPONENT * p1].count1Basis);
    }
    
    fprintf (out,"Tensor = {\n");
    flag2 = 0;
    
    for ( r = 0; r < CanonicalRank(f1, output, spin ) ; r++){
        if ( ! flag2 ){
            flag2 =1 ;
        }else
            fprintf (out,",\n");
        fprintf (out,"{\n");
        flag3 = 0;
        
        
        for ( space = 0; space < SPACE ; space++)
            if ( f1->sinc.rose[space].body != nada)
            {
            
            if ( ! flag3 ){
                flag3 = 1 ;
            }else fprintf (out,",\n");
            
            fprintf (out,"{\n");
            flag4 = 0;
            
            for ( l = 0 ; l < M[space] ;l++){
                if ( ! flag4 ){
                    flag4 = 1 ;
                }else
                    fprintf (out,",\n");
                fprintf (out,"%1.15lf\n", streams(f1,output,spin,space)[r*M[space]+l]);
                
            }
            fprintf (out,"\n}");
        }
        fprintf (out,"\n}");
    }
    fprintf (out,"\n}");
}

void tFilename (char * cycleName, INT_TYPE count, INT_TYPE body ,INT_TYPE IRREP, INT_TYPE cmpl, char * filename){
    sprintf(filename,"%s.%d.%d_mac",cycleName,count,cmpl);
   // printf("%s\n", filename);
}


DCOMPLEX tFromReadToFilename (char * cycleName, char * read , char * filename,INT_TYPE cmplFlag, INT_TYPE cmpl){
    double Occ,iOcc=0;
    INT_TYPE number,si,str0;
    char tokens[2][MAXSTRING];
    char * token = &*tokens[0];
    char * pa = &* tokens[1];
    token = strtok(read, "\"");
    /* walk through other tokens */
    pa = strtok(NULL, "\"");
    if ( ! cmplFlag )
        si = sscanf ( pa, ",%d,%lf",&number, &Occ );
    else{
        si = sscanf ( pa, ",%d,%lf,%lf",&number, &Occ,&iOcc );
        si--;
    }
    if ( si == 2 || si == 3  ) {
        tFilename(token, number, 0, 0, cmpl, filename);
        return Occ+iOcc;
    }else {
        return 0.;
    }
}

#if 1
INT_TYPE inputFormat(struct field * f1,char * name,  enum division buffer, INT_TYPE input){
    size_t maxRead = MAXSTRING;
    char input_line [maxRead];
    double value,lvalue;
    char * inputPt= input_line;;
    char c;
    //    struct calculation * c2 = malloc( sizeof(struct calculation));
    //    initCalculation(c2);
    INT_TYPE head, genus;
    INT_TYPE Nbody[SPACE],parts,p1,comp,i,M[SPACE],r1,r,space,flag2,flag3,flag4,l,sp,sy,N1[SPACE];

    
    FILE * in = NULL;
    in = fopen(name, "r");
    if ( in == NULL ){
        return 0;
        printf("failed to load %s\n", name);
        
    }
    //    broke = readInput(c2,in);
    //    finalizeInit(c2);
    
    //  getline(&inputPt,&maxRead, in   );
    //  getline(&inputPt,&maxRead, in   );
    getline(&inputPt,&maxRead, in   );
    INT_TYPE ct = 0;
    
#ifndef MKL
//    fprintf(out,"component = %d\n", COMPONENT);
//    fprintf(out,"genus = %d \n", species(f1, output) );
//    fprintf(out,"header = %d\n", header(f1, output));
//    fprintf(out,"body = %d\n", bodies(f1, output));
//    fprintf(out,"canonRank = %d\n",CanonicalRank(f1, output,spin ) );
//    fprintf(out,"spin = %d\n", spin);
//    fprintf(out,"symmetry = %d\n", f1->sinc.tulip[output].value.symmetry);
//    fprintf(out,"count1Basis = %d\n",  f1->sinc.rose[0].count1Basis);
//    fprintf(out,"lattice %f\n",f1->sinc.rose[0].lattice);
    sscanf(inputPt,"component = %d",&comp);
    getline(&inputPt,&maxRead, in   );
    if ( comp != COMPONENT )
    {
        printf("component \n");
        exit(0);
    }
    sscanf(inputPt, "genus = %d", &genus );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "header = %d", &head );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "canonRank = %d", &r1 );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "spin = %d", &sp );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "symmetry = %d", &sy );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "particles = %d", &parts );
    getline(&inputPt,&maxRead, in   );
    
    for ( p1 = 0 ; p1 < parts ; p1++){

        sscanf(inputPt, "%cbody = %d", &c,&Nbody[p1] );
        getline(&inputPt,&maxRead, in   );
        sscanf(inputPt, "%ccount1Basis = %d",&c, &N1[p1] );
        getline(&inputPt,&maxRead, in   );
        for ( i = 0; i < COMPONENT ; i++)
            M[p1*COMPONENT + i] = pow(N1[p1],Nbody[p1] * genus );

    }
        flag2 = 0;
     //   printf("header = %d\ngenus = %d %d %d\n%d %d %d %d \n",head, genus,parts,  r1 , sp, Nbody[0], N1[0],M[0]);
#else
    sscanf(inputPt,"component = %lld",&comp);
    getline(&inputPt,&maxRead, in   );
    if ( comp != COMPONENT )
    {
        printf("component \n");
        exit(0);
    }

    sscanf(inputPt, "genus = %lld", &genus );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "header = %lld", &head );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "canonRank = %lld", &r1 );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "spin = %lld", &sp );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "symmetry = %lld", &sy );
    getline(&inputPt,&maxRead, in   );
    
    sscanf(inputPt, "particles = %lld", &parts );
    getline(&inputPt,&maxRead, in   );
    
    for ( p1 = 0 ; p1 < parts ; p1++){
        sscanf(inputPt, "%cbody = %lld",&c, &Nbody[p1] );
        getline(&inputPt,&maxRead, in   );
        sscanf(inputPt, "%ccount1Basis = %lld",&c, &N1[p1] );
        getline(&inputPt,&maxRead, in   );
        
        for ( i = 0; i < COMPONENT ; i++)
            M[p1*COMPONENT + i] = Power(N1[p1],Nbody[p1] * genus );
    }
        //   printf("header = %lld\ngenus = %lld %lld %f\n%lld %lld %lld %lld %lld \n",head, genus,Nbody, D,  r1 , sp, M[0], M[1],M[2]);
    flag2 = 0;

#endif
  
    if ( input == 0 ){
        fclose(in);
        return genus;
    }
    if ( 100 <= input && input < 200 ){
        fclose(in);
        return Nbody[input - 100];
    }
    if ( 200 <= input && input < 300 ){
        fclose(in);
        return N1[input - 200];
    }

    if ( input == 2 ){
        fclose(in);
        return r1;
    }
    
    if ( r1 > part(f1, buffer)){
        printf("increase part %d %d\n", r1,part(f1,buffer));
        fflush(stdout);
        exit(0);
    }

    f1->sinc.tulip[buffer].value.symmetry = sy;
    
    if ( r1 > part(f1, buffer ) ){
        printf("io error\n");
        exit(0);
    }

    for ( r = 0; r < r1 ; r++){
//        if ( ! flag2 ){
//            flag2 =1 ;
//        }else
//            getline(&inputPt,&maxRead, in   );
        getline(&inputPt,&maxRead, in   );
        flag3 = 0;

        for ( space = 0; space < SPACE ; space++)
            if ( f1->sinc.rose[space].body != nada)
{
//            if ( ! flag3 ){
//                flag3 = 1 ;
//            }else     getline(&inputPt,&maxRead, in   );

            getline(&inputPt,&maxRead, in   );
            flag4 = 0;

            for ( l = 0 ; l < M[space] ;l++){
//                if ( ! flag4 ){
//                    flag4 = 1 ;
//                }else
//                    getline(&inputPt,&maxRead, in   );

                getline(&inputPt,&maxRead, in   );
                if ( ! sscanf (inputPt,"%lf\n",&( value)) ){
                    printf("--%lld :%s:\n",ct,inputPt);
                    exit(0);
                };

                ct++;
//                if ( value == lvalue )
//                    printf("%lld %1.15f %lld %lld\n", ct,value,r,l);
                lvalue = value;
                streams(f1, buffer, sp, space)[ M[space] * r+ l ] = value;
                getline(&inputPt,&maxRead, in   );
            //    printf(":%s:\n",inputPt);

            }
            getline(&inputPt,&maxRead, in   );
        }
        getline(&inputPt,&maxRead, in   );
        //   printf(">%lld %f\n",r, cblas_dnrm2( M[0],streams(f1, buffer, sp, 0)+ M[1] *r , 1  )*cblas_dnrm2( M[1],streams(f1, buffer, sp, 1)+ M[1] *r , 1  )*cblas_dnrm2( M[2],streams(f1, buffer, sp, 2)+ M[2] *r , 1  ));
    }
    getline(&inputPt,&maxRead, in   );
    
    
    f1->sinc.tulip[buffer].Current[sp] = r1;
    
    
    
    if ( input == 1 ){
        fclose(in);
        return ct;
    }
    
    
    return 0;
}
#endif



//INT_TYPE printOutput ( struct field * f1,INT_TYPE number){
//    INT_TYPE i,NC,numC;
//    char printCom[MAXSTRING];
//    FILE * printOut = stdout;
//    double radius,res;
//    assignCores(f1,1);
//    if ( header == NULL ){
//        printf("file\n");
//        exit(0);
//    }
//    size_t ms = MAXSTRING,read;
//    char input_line[MAXSTRING];
//    char * mask = input_line;
//    char name[MAXSTRING];
//    if ( f1->mem1->rt->printFlag % 2 ){
//    
//    }else {
//        return 0;
//    }
//    numC = f1->mem1->rt->monteCarlo;
//    NC = f1->mem1->rt->samples;
//    //Rs^3 4pi/3 = Volme / Ne
//    double Volume = 1;
//    INT_TYPE space;
//    for ( space = 0 ;space < SPACE ; space++)
//        Volume *= vectorLen(f1,space)*lattice(f1, space);
//    
//    
//    double Rs = pow(Volume * 3./4./pi/f1->body,0.333333333);
//    
//    {
//#ifdef OMP
//#pragma omp parallel for private (i,radius,res) schedule(dynamic,1)
//#endif
//        for ( i=0 ; i <= NC ; i++ )
//        {
//            radius = pow(Volume,1/3.)*(i*1./NC);
//            if ( f1->mem1->rt->runFlag == 7 )
//                radius *= 3;
//            else if (  ! i )
//                continue;
//            res =  tComputeRadialPlot( f1 ,number,f1->mem1->rt->runFlag,radius ,numC);
//            if ( f1->mem1->rt->runFlag == 7 ){
//                res *= Volume/4./pi;
//            }
//            else {
//            }
//            fprintf(printOut,"%d,%d,%d,%d,%12.6f,%12.6f,%12.16f\n",bodies(f1,eigenVectors),(f1->mem1->rt->runFlag)%2,(f1->mem1->rt->runFlag/2)%2, (f1->mem1->rt->runFlag/4)%2,radius,radius/Rs,res);
//            fflush(printOut);
//        }
//    }
//    return NC;
//}

//INT_TYPE printVectorOutput ( struct field * f1,INT_TYPE number){
//    INT_TYPE i,j,k,l,NC,numC;
//    char printCom[MAXSTRING];
//    FILE * printOut = stdout;
//    double displacement[3],res;
//    assignCores(f1,1);
//    //  FILE * header = fopen(f1->mem1->fileList,"r");
//    if ( header == NULL ){
//        printf("file\n");
//        exit(0);
//    }
//    size_t ms = MAXSTRING,read;
//    char input_line[MAXSTRING];
//    char * mask = input_line;
//    char name[MAXSTRING];
//    //  getline(&mask,&ms,header);
//    if ( (f1->mem1->rt->printFlag/2)% 2  ){
//    }
//    else {
//        return 0;
//    }
//    numC = f1->mem1->rt->monteCarlo;
//    NC = f1->mem1->rt->samples;
//    //   fclose(header);
//    {
//#ifdef OMP
//#pragma omp parallel for private (l,i,j,k,displacement,res) schedule(dynamic,1)
//#endif
//        
//        for ( l = 0; l < NC*NC*NC ; l++)
//        {
//            i = (l % NC) - (NC-1)/2;
//            j = ((l/NC) % NC) - (NC-1)/2;
//            k = ((l/(NC*NC)) % NC) - (NC-1)/2;
//
//            displacement[0] = f1->sinc.N1*f1->sinc.d*(i*(1./NC));
//            displacement[1] = f1->sinc.N1*f1->sinc.d*(j*(1./NC));
//            displacement[2] = f1->sinc.N1*f1->sinc.d*(k*(1./NC));
//
//#ifdef Bill
//            if ( i != 0 )
//                continue;
//#endif
//            
//            res =  tComputeVectorPlot( f1 ,number,f1->mem1->rt->runFlag,displacement ,numC);
//            fprintf(printOut,"%d,%d,%d,%d,%12.6f,%12.6f,%12.6f,%d,%d,%d,%12.16f\n",bodies(f1,eigenVectors),(f1->mem1->rt->runFlag)%2,(f1->mem1->rt->runFlag/2)%2, (f1->mem1->rt->runFlag/4)%2, displacement[0],displacement[1],displacement[2],i,j,k,res);
//            fflush(printOut);
//        }
//    }
//    return NC;
//}



//double tComputeRadialPlot(struct field * f1,INT_TYPE number,  INT_TYPE class,  double radius,INT_TYPE numC )
//{
//    double res= 0.,err=1.;
//
//    INT_TYPE N1 = f1->sinc.N1;
//    INT_TYPE N12 = (N1-1)/2;
//    double xl[5],xu[5];
//    INT_TYPE space,num=0;
//
//
//
//
//        struct fieldArray fA ;
//        fA.f1 = f1;
//        fA.number = number;
//
//        xl[0] = 0.;
//        xu[0] = pi;
//
//        xl[1] = 0.;
//        xu[1] = 2*pi;
//
//
//
//        fA.class = class;
//
//
//
//
//        space = 0;
//        if ( class % 2 ){
//
//            //        xl[space+2] = - pi/f1->sinc.d;
//            //        xu[space+2] =   pi/f1->sinc.d;
//            xl[space+2] = - N12*f1->sinc.d;
//            xu[space+2] =   N12*f1->sinc.d;
//
//        }
//        else {
//
//            xl[space+2] = - N12*f1->sinc.d;
//            xu[space+2] =   N12*f1->sinc.d;
//
//        }
//        space = 1;
//        if ( (class/2) % 2 ){
//
//            //        xl[space+2] = - pi/f1->sinc.d;
//            //        xu[space+2] =   pi/f1->sinc.d;
//            xl[space+2] = - N12*f1->sinc.d;
//            xu[space+2] =   N12*f1->sinc.d;
//
//        }
//        else {
//
//            xl[space+2] = - N12*f1->sinc.d;
//            xu[space+2] =   N12*f1->sinc.d;
//
//        }
//        space = 2;
//        if ( (class/4) % 2 ){
//
//            //        xl[space+2] = - pi/f1->sinc.d;
//            //        xu[space+2] =   pi/f1->sinc.d;
//            xl[space+2] = - N12*f1->sinc.d;
//            xu[space+2] =   N12*f1->sinc.d;
//
//        }
//        else {
//
//            xl[space+2] = - N12*f1->sinc.d;
//            xu[space+2] =   N12*f1->sinc.d;
//
//        }
//
//
//        INT_TYPE rank = 0,dim = 5;
//
//        if ( bodies (f1, eigenVectors ) == one )
//            dim = 2;
//        else
//            dim  = 5;
//#ifdef OMP
//        rank = omp_get_thread_num();
//        gsl_monte_function G;
//        G.f = &evaluateDensityBracket;
//        G.dim = dim;
//        G.params = &fA;
//
//        INT_TYPE iter = 0;
//        const gsl_rng_type * T;
//        gsl_rng  * r;
//        ((struct fieldArray*)(G.params))->radius[rank] = radius;
//        gsl_monte_vegas_state * s;
//        gsl_rng_env_setup();
//        T = gsl_rng_default;
//        r = gsl_rng_alloc(T);
//        s = gsl_monte_vegas_alloc(G.dim);
//
//        gsl_monte_vegas_integrate ( & G , xl, xu, dim , numC, r,s,&res, &err );
//        do
//        {
//            gsl_monte_vegas_integrate ( & G , xl, xu, dim , numC, r,s,&res, &err );
//            //    printf("%f %f %f\n", res,err, gsl_monte_vegas_chisq (s));
//        } while (fabs ((gsl_monte_vegas_chisq (s) - 1.0) > 0.1 || max(1.,fabs(res)) > f1->mem1->rt->TARGET)&& num++ < 100);
//
//        gsl_monte_vegas_free(s);
//        gsl_rng_free(r);
//#else
//
//        //    {
//        //        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
//        //
//        //        gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
//        //                                   &res, &err);
//        //        display_results ("vegas warm-up", res, err);
//        //
//        //        printf ("converging...\n");
//        //
//        //        do
//        //        {
//        //            gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
//        //                                       &res, &err);
//        //            printf ("result = % .6f sigma = % .6f "
//        //                    "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
//        //        }
//        //        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//        //
//        //        display_results ("vegas final", res, err);
//        //
//        //        gsl_monte_vegas_free (s);
//        //    }
//#endif
//
//        return res;
//}


//double tComputeVectorPlot(struct field * f1,INT_TYPE number,  INT_TYPE class,  double *displacement,INT_TYPE numC )
//{
//    double res= 0.,err=1.;
//
//    INT_TYPE N1 = f1->sinc.N1;
//    INT_TYPE N12 = (N1-1)/2;
//    double xl[5],xu[5];
//    INT_TYPE space,num=0;
//
//    struct fieldArray fA ;
//    fA.f1 = f1;
//    fA.number = number;
//    fA.class = class;
//    space = 0;
//
//    if ( bodies(f1, eigenVectors ) == one ){
//        res = evaluateVectorBracket(displacement, 3,(void*)(&fA));
//    }
//    else {
//
//
//
//        if ( class % 2 ){
//            xl[space] = - N12*f1->sinc.d;
//            xu[space] =   N12*f1->sinc.d;
//            //
//            //        xl[space] = - pi/f1->sinc.d;
//            //        xu[space] =   pi/f1->sinc.d;
//        }
//        else {
//
//            xl[space] = - N12*f1->sinc.d;
//            xu[space] =   N12*f1->sinc.d;
//
//        }
//        space = 1;
//        if ( (class/2) % 2 ){
//            xl[space] = - N12*f1->sinc.d;
//            xu[space] =   N12*f1->sinc.d;
//
//            //        xl[space] = - pi/f1->sinc.d;
//            //        xu[space] =   pi/f1->sinc.d;
//        }
//        else {
//
//            xl[space] = - N12*f1->sinc.d;
//            xu[space] =   N12*f1->sinc.d;
//
//        }
//        space = 2;
//        if ( (class/4) % 2 ){
//            xl[space] = - N12*f1->sinc.d;
//            xu[space] =   N12*f1->sinc.d;
//
//            //        xl[space] = - pi/f1->sinc.d;
//            //        xu[space] =   pi/f1->sinc.d;
//        }
//        else {
//
//            xl[space] = - N12*f1->sinc.d;
//            xu[space] =   N12*f1->sinc.d;
//
//        }
//
//
//
//        INT_TYPE rank = 0,dim = 3;
//#ifndef APPLE
//#ifdef OMP
//        rank = omp_get_thread_num();
//#else
//        rank = 0;
//#endif
//        gsl_monte_function G;
//        G.f = &evaluateVectorBracket;
//        G.dim = 3;
//        G.params = &fA;
//
//        INT_TYPE iter = 0;
//        const gsl_rng_type * T;
//        gsl_rng  * r;
//        ((struct fieldArray*)(G.params))->displacement[rank][0] = displacement[0];
//        ((struct fieldArray*)(G.params))->displacement[rank][1] = displacement[1];
//        ((struct fieldArray*)(G.params))->displacement[rank][2] = displacement[2];
//        gsl_monte_vegas_state * s;
//        gsl_rng_env_setup();
//        T = gsl_rng_default;
//        r = gsl_rng_alloc(T);
//        s = gsl_monte_vegas_alloc(G.dim);
//
//#ifdef BILL
//        double x0[3];
//        x0[0] = 0.;
//        x0[1] = 0.;
//        x0[2] = 0.;
//        res = evaluateVectorBracket(x0,3,&fA );
//#else
//
//        gsl_monte_vegas_integrate ( & G , xl, xu, dim , numC, r,s,&res, &err );
//        do
//        {
//            gsl_monte_vegas_integrate ( & G , xl, xu, dim , numC, r,s,&res, &err );
//            // printf("%f %f %f\n", res,err, gsl_monte_vegas_chisq (s));
//        } while ((fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.1 || err/max(1.,fabs(res))  > f1->mem1->rt->TARGET) && num++ < 100);
//
//        gsl_monte_vegas_free(s);
//        gsl_rng_free(r);
//#endif
//#else
//
//        //    {
//        //        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
//        //
//        //        gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
//        //                                   &res, &err);
//        //        display_results ("vegas warm-up", res, err);
//        //
//        //        printf ("converging...\n");
//        //
//        //        do
//        //        {
//        //            gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
//        //                                       &res, &err);
//        //            printf ("result = % .6f sigma = % .6f "
//        //                    "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
//        //        }
//        //        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//        //
//        //        display_results ("vegas final", res, err);
//        //
//        //        gsl_monte_vegas_free (s);
//        //    }
//#endif
//    }
//    return res;
//}

INT_TYPE tFillBasis(Stream_Type ** pt/*3 vectors*/, double * coordinates/*3 numbers*/, INT_TYPE class,INT_TYPE N1,double lattice){
    INT_TYPE space,h,N12 = (N1-1)/2,aPeriodic;
    for ( space = 0 ; space < SPACE ;space++)
    {
        
        if ( space == 0 ){
            aPeriodic = class%2;
        }
        else if ( space == 1 ){
            aPeriodic = (class/2)%2;
        }
        else
            aPeriodic = (class/4)%2;
        
        if ( ! aPeriodic ){
            double * position = coordinates;
            for ( h = 0 ; h < N1 ; h++)
                pt[space][h] = Sinc(lattice,position[space]-lattice*(h-N12))/sqrt(lattice);
        }
//        else
//        {
//                double periodicPosition;
//                periodicPosition = coordinates[space];
//                for ( h = 0 ; h < N1 ; h++)
//                    pt[space][h] = periodicSinc(lattice,periodicPosition-lattice*(h-N12),N1)/sqrt(lattice);
//
//        }
            
//        {
//            double phase ;
//            double * waveNumber = coordinates;
//            if ( class == 1 )
//                phase = 0.;
//            else
//                phase = 0.5;
//            for ( h = 0 ; h < N1 ; h++){
//                pt[space][h] = cos(waveNumber[space]*lattice*(h+phase))/sqrt(N1);
//            }
//            for ( h = 0 ; h < N1 ; h++)
//                if ( space == 0)
//                    pt[space][h] = -sin(waveNumber[space]*lattice*(h+phase))/sqrt(N1);
//                else
//                    pt[space][h] = cos(waveNumber[space]*lattice*(h+phase))/sqrt(N1);
//
//            for ( h = 0 ; h < N1 ; h++)
//                if ( space == 1)
//                    pt[space][h] = -sin(waveNumber[space]*lattice*(h+phase))/sqrt(N1);
//                else
//                    pt[space][h] = cos(waveNumber[space]*lattice*(h+phase))/sqrt(N1);
//
//            for ( h = 0 ; h < N1 ; h++)
//                if ( space == 2)
//                    pt[space][h] = -sin(waveNumber[space]*lattice*(h+phase))/sqrt(N1);
//                else
//                    pt[space][h] = cos(waveNumber[space]*lattice*(h+phase))/sqrt(N1);
//
//        }
        
    }
            return 1;
    
}

//double evaluateDensityBracket( double x [], size_t dim , void * params ){
//    double y[3],sum=0.;;
//    INT_TYPE space,info,cl,cmpl;
//    Stream_Type *pt[SPACE];
//    INT_TYPE rank = 0;
//#ifdef OMP
//    rank = omp_get_thread_num();
//#endif
//    enum division wavefunction,product = squareVector;
//    struct field *f1 = ((struct fieldArray *)(params))->f1;
//    double radius = ((struct fieldArray*)(params))->radius[rank];
//    INT_TYPE class = ((struct fieldArray*)(params))->class;
//    INT_TYPE number = ((struct fieldArray*)(params))->number;
//
//    if ( dim == 5 ){
//        
//        {
//            {
//                for ( space = 0; space < SPACE ; space++)
//                    pt[space] = streams(f1, diagonal1VectorA,rank,space);
//                y[0] = x[0+2];
//                y[1] = x[1+2];
//                y[2] = x[2+2];
//                
//                f1->sinc.tulip[diagonal1VectorA].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
//                
//                
//                for ( space = 0; space < SPACE ; space++)
//                    pt[space] = streams(f1, diagonal1VectorB,rank,space);
//                y[0] = cos(x[0])*           radius+x[0+2];
//                y[1] = sin(x[0])*cos(x[1])* radius+x[1+2];
//                y[2] = sin(x[0])*sin(x[1])* radius+x[2+2];
//                
//                f1->sinc.tulip[diagonal1VectorB].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
//                
//                
//                
//                for ( cmpl = 0; cmpl < 2 ; cmpl++)
//                for ( wavefunction = eigenVectors ; wavefunction < eigenVectors+number ; wavefunction++){
////                    if ( f1->sinc.tulip[printOperator].linkNext != nullName || CanonicalRank(f1, printOperator, 0) || CanonicalRank(f1, printOperator, 1)){
////                        tEqua(f1, product,rank, wavefunction,);
////                        tHXpX(rank, f1, printOperator, 0, 1., 0., product, f1->mem1->rt->CONVERGENCE, part(f1,product));
////                    }else {
////                        product = wavefunction;
////                    }
//
//                    
//                    if ( bodies ( f1, wavefunction ) == two ){
//                        //oneVector :  basisRank  and ONE
//                        tMultiplyMP(rank, &info,f1, 1.0, -1, oneVector , rank, 'N', wavefunction, cmpl, 'N', diagonal1VectorB, rank);
//                        sum+= sqr(tMultiplyMP(rank,&info, f1, 1.0, -1, nullVector , 0, 'T', oneVector, rank, 'N', diagonal1VectorA, rank));
//                    }    if ( bodies ( f1, wavefunction ) == three ){
//                        //oneVector : basisRank and BODY - TWO
//                        //diagonal2VectorA : 1 and TWO
//                        f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
//                        tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
//                        tMultiplyMP(rank, &info,f1, 1.0, -1, oneVector , rank, 'N', wavefunction, cmpl, 'N', diagonal2VectorA, rank);
//                        sum += tMultiplyMP(rank, &info,f1, 1.0, -1, nullVector , 0, 'T', oneVector, rank, 'N', oneVector, rank);
//                    }   else   if ( bodies ( f1, wavefunction ) == four ){
//                        //oneVector : basisRank and BODY - TWO
//                        //diagonal2VectorA : 1 and TWO
//                        f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
//                        tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
//                        tMultiplyMP(rank, &info,f1, 1.0, -1, twoVector , rank, 'N', wavefunction, cmpl, 'N', diagonal2VectorA, rank);
//                        sum += tMultiplyMP(rank,&info, f1, 1.0, -1, nullVector , 0, 'T', twoVector, rank, 'N', twoVector, rank);
//                    }
//                }
//            }
//        }
//    }
//    else if (dim == 2 ){
//        {
//            for ( space = 0; space < SPACE ; space++)
//                pt[space] = streams(f1, diagonal1VectorA,rank,space);
//            y[0] = cos(x[0])*           radius;
//            y[1] = sin(x[0])*cos(x[1])* radius;
//            y[2] = sin(x[0])*sin(x[1])* radius;
//            
//            f1->sinc.tulip[diagonal1VectorA].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
//            for ( cmpl = 0; cmpl < 2 ; cmpl++)
//
//            for ( wavefunction = eigenVectors ; wavefunction < eigenVectors+number ; wavefunction++){
////                if ( f1->sinc.tulip[printOperator].linkNext != nullName || CanonicalRank(f1, printOperator, 0) || CanonicalRank(f1, printOperator, 1)){
////                    tHXpX(rank, f1, printOperator, 0, 1., 0., product, f1->mem1->rt->CONVERGENCE, part(f1,product));
////                }else {
////                    product = wavefunction;
////                }
//                
//                sum += sqr(tMultiplyMP(rank, &info, f1, 1.0, -1, nullVector, rank, 'T', wavefunction, cmpl, 'N', diagonal1VectorA, rank));
//            }
//        }
//    }
//    if( ! class)
//        return sqr(radius)*sin(x[0])*sum;
//    else
//        return sin(x[0])*sum;
//
//}


//double evaluateVectorBracket( double x [], size_t dim , void * params ){
//    double y[3],sum=0.;;
//    INT_TYPE space,info,cl,cmpl;
//    Stream_Type *pt[SPACE];
//    INT_TYPE rank = 0;
//#ifdef OMP
//    rank = omp_get_thread_num();
//#endif
//    double *displacement = ((struct fieldArray*)(params))->displacement[rank];
//    enum division wavefunction,product = squareVector;
//    struct field *f1 = ((struct fieldArray *)(params))->f1;
//    INT_TYPE class = ((struct fieldArray*)(params))->class;
//    INT_TYPE number = ((struct fieldArray*)(params))->number;
//
//    if ( dim == 3 ){
//        
//        {
//            {
//                for ( space = 0; space < SPACE ; space++)
//                    pt[space] = streams(f1, diagonal1VectorA,rank,space);
//                y[0] = x[0];
//                y[1] = x[1];
//                y[2] = x[2];
//                
//                f1->sinc.tulip[diagonal1VectorA].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
//                
//                
//                for ( space = 0; space < SPACE ; space++)
//                    pt[space] = streams(f1, diagonal1VectorB,rank,space);
//                y[0] = displacement[0]+x[0];
//                y[1] = displacement[1]+x[1];
//                y[2] = displacement[2]+x[2];
//                
//                f1->sinc.tulip[diagonal1VectorB].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
//                
//                for ( cmpl = 0; cmpl < 2 ; cmpl++)
//
//                for ( wavefunction = eigenVectors ; wavefunction < eigenVectors+number ; wavefunction++){
////                    if ( f1->sinc.tulip[printOperator].linkNext != nullName || CanonicalRank(f1, printOperator, 0) || CanonicalRank(f1, printOperator, 1)){
////                        tHXpX(rank, f1, printOperator, 0, 1., 0., product, f1->mem1->rt->CONVERGENCE, part(f1,product));
////                    }else {
////                        product = wavefunction;
////                    }
//
//                    if ( bodies ( f1, wavefunction ) == one ){
//                        sum += sqr(tMultiplyMP(rank, &info,f1, 1.0, -1, nullVector , 0, 'T', wavefunction, cmpl, 'N', diagonal1VectorA, rank));
//                    } else
//                    if ( bodies ( f1, wavefunction ) == two ){
//                        //oneVector :  basisRank  and ONE
//                        tMultiplyMP(rank, &info,f1, 1.0, -1, oneVector , rank, 'N', wavefunction, cmpl, 'N', diagonal1VectorB, rank);
//                        sum+= sqr(tMultiplyMP(rank,&info, f1, 1.0, -1, nullVector , 0, 'N', oneVector, rank, 'N', diagonal1VectorA, rank));
//                    }  else   if ( bodies ( f1, wavefunction ) == three ){
//                        //oneVector : basisRank and BODY - TWO
//                        //diagonal2VectorA : 1 and TWO
//                        f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
//                        tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
//                        tMultiplyMP(rank, &info,f1, 1.0, -1, oneVector , rank, 'N', wavefunction, cmpl, 'N', diagonal2VectorA, rank);
//                        sum += tMultiplyMP(rank, &info,f1, 1.0, -1, nullVector , 0, 'T', oneVector, rank, 'N', oneVector, rank);
//                    }   else   if ( bodies ( f1, wavefunction ) == four ){
//                        //oneVector : basisRank and BODY - TWO
//                        //diagonal2VectorA : 1 and TWO
//                        f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
//                        tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
//                        tMultiplyMP(rank, &info,f1, 1.0, -1, twoVector , rank, 'N', wavefunction, cmpl, 'N', diagonal2VectorA, rank);
//                        sum += tMultiplyMP(rank,&info, f1, 1.0, -1, nullVector , 0, 'T', twoVector, rank, 'N', twoVector, rank);
//                    }
//                }
//            }
//        }
//    }
//
//    return sum/pow(f1->sinc.d,SPACE);
//    
//}

INT_TYPE tLoadEigenWeights (struct calculation * c1, char * filename, enum division inputVectors){
    INT_TYPE space,ct = 0,ct2,number,class,weight,cmpl;
    FILE * in = NULL;
    in = fopen(filename, "r");
    struct field *f1 = &c1->i.c;
    if ( in == NULL ){
        printf("file of occupations is missing\n");
        exit(0);
    }
    INT_TYPE flagLoad;
    DCOMPLEX ov;
    size_t ms = MAXSTRING;
    char input_line[MAXSTRING];
    char * mask = input_line;
    DCOMPLEX Occ;
    char name[MAXSTRING];
    while (1){
        if (  getline(&mask,&ms,in) > 0 ){
            if ( (!comment(input_line)) && (strlen(input_line) > 1) ){
                Occ = 0.;
                flagLoad = 0;
                for ( cmpl = 0; cmpl < spins(&c1->i.c, inputVectors); cmpl++)
                {
                    Occ = tFromReadToFilename(NULL, input_line,  name, spins(f1,eigenVectors)-1,cmpl);
                    if ( cabs(Occ) > c1->rt.TARGET){
                        {
                            f1->sinc.tulip[inputVectors+ct].Current[cmpl] = 0;

                            struct calculation c2 = *c1;
                            c2.rt.body = inputFormat(f1, name, nullName, 100);
                            c2.i.c.mem1->bootedMemory =0;
                            c2.i.c.sinc.tulip = NULL;
                            c2.i.sectors = 0;
                            c2.i.decomposeRankMatrix = 1;
                            for ( space = 0; space <= SPACE ; space++)
                                c2.i.c.sinc.rose[space].stream = NULL;
                            c2.rt.boot = noMatrices;
                            c2.i.nStates= 0;
                            c2.i.vectorOperatorFlag= 0;
                            c2.i.springFlag= 0;
                            c2.i.heliumFlag = 0;
                            c2.i.nTargets = 0;
                            c2.i.qFloor = 0;
#ifdef OMP
                            c2.rt.NLanes = spins(f1, eigenVectors);
#endif
#ifdef MKL
                            c2.rt.NParallel = 1;
#endif
                            c2.i.canonRank = imax(c1->i.bRank,2+inputFormat(f1, name, nullName, 2));
                            c2.i.c.oneBody.func.fn = nullFunction;
                            c2.i.c.twoBody.func.fn = nullFunction;
                            if ( SPACE > COMPONENT ){
                                if ( c2.rt.runFlag )
                                    c2.i.around = (inputFormat(f1, name, nullName, 201)/2-1)/2;
                                else
                                    c2.i.around = (inputFormat(f1, name, nullName, 201)-1)/2;
                                c2.i.D = c1->i.D * pow( (2.* c1->i.around + 1.) /(2.*c2.i.around + 1),1.);
                            }
                            if ( c2.rt.runFlag )
                                c2.i.epi =(inputFormat(f1, name, nullName, 200)/2-1)/2;
                            else
                                c2.i.epi =(inputFormat(f1, name, nullName, 200)-1)/2;

                            c2.i.d = c1->i.d* pow( (2.* c1->i.epi + 1.) /(2.*c2.i.epi + 1),c2.i.attack);

                            
                            printf("\n\n\t  \tBox %d \t: Lattice  %1.3f \t\n\t\n",2* c2.i.epi + 1,c2.i.d);
                            
                            if ( SPACE > 3 )
                                printf("\n\n\t  \tBox %d \t: Lattice  %1.3f \t\n\t\n",2* c2.i.around + 1,c2.i.D);
                            

                            
                            
                            iModel(&c2);
                            c2.i.c.sinc.tulip[totalVector].Current[0] = 0;
                            inputFormat(&c2.i.c, name, totalVector,1);
                            xEqua(f1,inputVectors+ct, cmpl, &c2.i.c, totalVector,0);
                            f1->sinc.tulip[inputVectors+ct].value.symmetry = c2.i.c.sinc.tulip[totalVector].value.symmetry;

                            printf("%s\tSA%d\n", name,f1->sinc.tulip[inputVectors+ct].value.symmetry);
                            flagLoad = 1;
                            fModel(&c2);
                        }
                        
                        
                    }
                }
                if ( ! flagLoad )
                    continue;
                if (( ct > c1->i.vectorOperatorFlag&& inputVectors == f1->sinc.vectorOperator )|| ( ct > c1->i.nStates&& inputVectors == eigenVectors  )){
                    printf("maxed out buffer of states\n");
                    exit(0);
                }
                if ( cimag(Occ) != 0. ){
                    printf("fixme\n");
                    exit(1);
                }

                tScale(&c1->i.c, inputVectors+ct, creal(Occ));
                
                ov = magnitude(f1, inputVectors+ct);

                if ( inputVectors >= c1->i.c.sinc.vectorOperator){
                    printf("Density %f\n",sqr(creal(ov)));
                    if ( sqr(cabs(ov)) > c1->rt.TARGET)
                        ct++;
                    else
                        printf("skipped\n");

                }
                else{
                    printf("Norm %f\n",(creal(ov)));
                    if ( (cabs(ov)) > c1->rt.TARGET)
                        ct++;
                    else
                        printf("skipped\n");
                }
            }
        }else {
            break;
        }
    }
    
    return ct;
};
