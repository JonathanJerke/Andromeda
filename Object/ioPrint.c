/**
 *  ioPrint.c
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

#include "ioPrint.h"

inta printVector (  calculation *c,  sinc_label f1, char * name,char * vectorName,  inta iv, inta irrep,mea * vector){
    if ( vector == NULL )
        return 1;
    inta iii;
    char str [SUPERMAXSTRING];
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
#ifdef COMPLEXME
        {
            if( f1.cmpl == 2)
                fprintf(outf, "\"%s\",%d,%15.15f,%15.15f\n",name, iii+1,creal(*vector),cimag(*vector));
            else
                fprintf(outf, "\"%s\",%d,%15.15f\n", name, iii+1,creal(*vector));
        }
#else
        fprintf(outf, "\"%s\",%d,%15.15f\n", name, iii+1,(*vector));
#endif
    }
    fclose(outf);
    return 0;
}

inta printOut(  calculation *c ,   field *f1,inta reset, inta lv,  division vector){
    inta irrep;
    inta jjj=1,cmpl;
    char str [SUPERMAXSTRING];
#ifndef APPLE
    mea one = 1.;
    if ( reset ) {
        FILE * outf ;
        sprintf(str, "%s.vector",c->name);
        outf = fopen (str,"w");
        fclose(outf);
    }
        {
                
                    printf("State%d: %1.15f, body%d, irrep%d, %1.15f\n", lv+1,f1->f.name[vector].value.value,bodies(f1->f,vector),f1->f.name[vector].value.symmetry, f1->f.name[vector].value.value2);
            
                    printVector(c,f1->f, c->name,c->name, lv,irrep, &one);
                    for ( cmpl = 0 ; cmpl < spins(f1->f, vector) ; cmpl++)
                    {
#ifndef APPLE
                        tFilename(c->name,lv+1,bodies(f1->f, vector) ,irrep, cmpl,str);
                        
                        
#ifdef writeHDF5
                        {
                            inta space;
                        for ( space = 0; space < SPACE ; space++)
                            if ( f1->f.canon[space].body != nada)
                                writeFast(f1->f, str, space, vector,cmpl);
                        }
#else
                        FILE * out = NULL;
                        out = fopen ( str,"w" );
                        if ( out != NULL ){
                            outputFormat(f1->f, out, vector,cmpl  );
                            fclose(out);
                        }
#endif
#endif
                    }
            
            }

    fflush(stdout);
#endif
    return 0;
}




void outputFormat(  sinc_label f1, FILE * out, division output ,inta spin){
    inta  dims,parts,p1,flag2,flag3,flag4,r,l,space, M[SPACE],prev;
    length(f1, output, M);

    dims = 0;
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.canon[space].body != nada)
            dims++;
    
    
    fprintf(out,"component = %d\n", dims);
      genusType g = species(f1, output);
    if (g == 3 )
        g = 1;
    
    fprintf(out,"genus = %d \n", g );
    fprintf(out,"header = %d\n", header(f1, output));
    fprintf(out,"canonRank = %d\n",CanonicalRank(f1, output,spin ) );
    fprintf(out,"spin = %d\n", spin);
    fprintf(out,"symmetry = %d\n", f1.name[output].value.symmetry);
    parts=1;
    prev = 0;
    for ( space= 1 ; space < SPACE; space++ )
        if ( f1.canon[space].body != nada )
        if (f1.canon[prev].label != f1.canon[space].label){
            parts++;
            prev = space;
        }
    fprintf(out,"particles = %d\n",parts);
    prev = 0;
    for ( p1 = 0 ; p1 < parts;p1++)
        for ( space = 0; space < SPACE ; space++)
            if ( f1.canon[space].body != nada)
                if (f1.canon[space].label == p1+1){
                    if ( species(f1, output)  == vector )
                        fprintf(out,"%cbody = %d\n", 'A'+p1,f1.canon[space].body);
                    else if ( species(f1, output)  == matrix || species(f1,output) == outerVector )
                        fprintf(out,"%cbody = %d\n", 'A'+p1,f1.name[output].space[space].body);
                    fprintf(out,"%ccount1Basis = %d\n", 'A'+p1, f1.canon[space].count1Basis);
                   // fprintf(out,"%caxes = %d\n", 'A'+p1, f1.canon[space].component);
                    break;
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
            if ( f1.canon[space].body != nada)
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

void tFilename (char * cycleName, inta count, inta body ,inta IRREP, inta cmpl, char * filename){
    sprintf(filename,"%s.%d.%d_mac",cycleName,count,cmpl);
}


mea tFromReadToFilename (char * cycleName, char * read , char * filename,inta cmplFlag, inta cmpl,char * title, inta *number){
    double Occ,iOcc=0;
    inta si;
    char tokens[2][MAXSTRING];
    char * token = &*tokens[0];
    char * pa = &* tokens[1];
    token = strtok(read, "\"");
    strcpy( title, token);
    /* walk through other tokens */
    pa = strtok(NULL, "\"");
    if ( ! cmplFlag )
        si = sscanf ( pa, ",%d,%lf",number, &Occ );
    else{
        si = sscanf ( pa, ",%d,%lf,%lf",number, &Occ,&iOcc );
        si--;
    }
    if ( si == 2  ) {
        tFilename(token, *number, 0, 0, cmpl, filename);
      //  printf("%s\n", filename);
        return Occ+iOcc;
    }else {
        printf("failed to load\n");
        exit(0);
        return 0.;
    }
}

inta ioArray(  calculation *c1,   field f,char * name,inta N1, floata * array, inta ioIn){
    inta space,i;
    field f2 ;
    initField(&f2);
    calculation c2;
    c2 = *c1;
    c2.i.numNames = 0;
    c2.i.numVectors = 0;
    c2.i.termNumber = 0;
    c2.rt.NLanes = 1;
    f2.f.rt = &c2.rt;
    f2.f.rt->phaseType = productKrylov;
    f2.i = f.i;
    f2.i.iRank = 1;
    f2.i.xRank = 1;
    c2.i.Lambda = 0;
    f2.i.Iterations = 1;
    f2.i.files = 0;
    f2.i.filesVectorOperator = 0;
    f2.i.qFloor = 0;
    resetA(f2.f.rt);
    blockA(f2.f.rt, blockTotalVectorBlock);
    blockA(f2.f.rt, blockTrainVectorsblock);
    blockA(f2.f.rt, blockCopyBlock);
    blockA(f2.f.rt, blockTransferBasisblock);
    blockA(f2.f.rt, blockMatrixElementsblock);
    blockA(f2.f.rt, blockPermutationsblock);
    blockA(f2.f.rt, blockParallelMultiplyblock);
    blockA(f2.f.rt, blockParallelMatrixElementblock);
    blockA(f2.f.rt, blockParallelPermuteblock);
    blockA(f2.f.rt, blockTotalVectorParallelBlock);
    blockA(f2.f.rt, blockComponentblock);
    blockA(f2.f.rt, blockDiagonalMatrixblock);
    blockA(f2.f.rt, blockCompressionBlock);

    f2.i.body = one;

    f2.f.boot = noMatrices;
    
    f2.i.canonRank  = 1;
    f2.i.nStates = 1;
    
    for ( space = 1 ; space < SPACE ; space++)
        f2.f.canon[space].body = nada;

    f2.f.canon[0].body = one;
    f2.f.canon[0].count1Basis = N1;
    f2.f.canon[0].basis = 4;
    //f2.f.canon[0].component = 1;
    f2.f.canon[0].label = 1;

    iModel(&c2,&f2);
    
    f2.f.name[eigenVectors].Current[0] = 1;
    if ( ioIn ){
        //IN
        inputFormat(f2.f, name, eigenVectors,1);
        floata *pt = streams(f2.f,eigenVectors,0,0);
        for ( i = 0 ; i < N1  ; i++)
            array[i] = pt[i];//column major matrix.
        //IN
    } else {
        //OUT.
        floata *pt = streams(f2.f,eigenVectors,0,0);
        for ( i = 0 ; i < N1  ; i++)
            pt[i] = array[i];//column major matrix.
        
#ifdef writeHDF5
        writeFast(f2.f, name, 0, eigenVectors,0);
#else
        FILE* file = fopen(name,"w");
        outputFormat(f2.f, file, eigenVectors, 0);
        fclose(file);
#endif
        //OUT
    }
    fModel(&f2.f);
    
    return 0;
}

inta inputFormat(  sinc_label f1,char * name,    division buffer, inta input){
#ifdef readHDF5
    inta space ,part = 1,prev = 0;
    if ( input == 0 )
        return readFast(f1,name, 0,0,buffer,0,0);
    if ( input == 2 )
        return readFast(f1,name, 2,0,buffer,0,0);
    if ( 100 <= input && input < 200 )
        return readFast(f1,name, 3,input-100,buffer,0,0);
    if ( 200 <= input  )
        return readFast(f1,name, 4,input-200,buffer,0,0);

    if ( input == 1 ){
        for ( space = 0; space < SPACE ; space++)
            if ( f1.canon[space].body != nada)
                readFast(f1,name,1,space,buffer,0,space);
        f1.name[buffer].Current[0] = readFast(f1,name, 2,0,buffer,0,0);
    }
    return 0;
#else
    int maxRead = MAXSTRING;
    char input_line [maxRead];
    double value,lvalue;
    char * inputPt= input_line;;
    char c;
    //      calculation * c2 = malloc( sizeof(  calculation));
    //    initCalculation(c2);
    inta head, genus;
    inta Nbody[SPACE],comps[SPACE],parts,p1,comp,i,M[SPACE],r1,r,space,flag2,flag3,flag4,l,sp,sy,N1[SPACE];

    FILE * in = NULL;
    in = fopen(name, "r");
    if ( in == NULL ){
        return 0;
    }
    fgets(inputPt,maxRead, in   );
    inta ct = 0;
    
#ifndef MKL
    sscanf(inputPt,"component = %d",&comp);
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "genus = %d", &genus );
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "header = %d", &head );
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "canonRank = %d", &r1 );
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "spin = %d", &sp );
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "symmetry = %d", &sy );
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "particles = %d", &parts );
    fgets(inputPt,maxRead, in   );
    space = 0;
    for ( p1 = 0 ; p1 < parts ; p1++){

        sscanf(inputPt, "%cbody = %d", &c,&Nbody[p1] );
        fgets(inputPt,maxRead, in   );
        sscanf(inputPt, "%ccount1Basis = %d",&c, &N1[p1] );
        fgets(inputPt,maxRead, in   );
        if ( 2 != sscanf(inputPt, "%caxes = %d", &c,&comps[p1] )){
                comps[p1] = COMPONENT;
            }else {
                fgets(inputPt,maxRead, in   );
            }
        for ( i = 0; i < comps[p1] ; i++)
            M[space++] = pow(N1[p1],Nbody[p1] * genus );
    
    }
        flag2 = 0;
     //   printf("header = %d\ngenus = %d %d %d\n%d %d %d %d \n",head, genus,parts,  r1 , sp, Nbody[0], N1[0],M[0]);
#else
    sscanf(inputPt,"component = %lld",&comp);
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "genus = %lld", &genus );
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "header = %lld", &head );
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "canonRank = %lld", &r1 );
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "spin = %lld", &sp );
    fgets(inputPt,maxRead, in   );
    sscanf(inputPt, "symmetry = %lld", &sy );
    fgets(inputPt,maxRead, in   );
    
    sscanf(inputPt, "particles = %lld", &parts );
    fgets(inputPt,maxRead, in   );
    space = 0;
    for ( p1 = 0 ; p1 < parts ; p1++){

        sscanf(inputPt, "%cbody = %lld", &c,&Nbody[p1] );
        fgets(inputPt,maxRead, in   );
        sscanf(inputPt, "%ccount1Basis = %lld",&c, &N1[p1] );
        fgets(inputPt,maxRead, in   );
        if ( 2 != sscanf(inputPt, "%caxes = %lld", &c,&comps[p1] )){
            comps[p1] = COMPONENT;
        }else {
            fgets(inputPt,maxRead, in   );
        }
        for ( i = 0; i < comps[p1] ; i++)
            M[space++] = pow(N1[p1],Nbody[p1] * genus );

    }
   // printf("header = %lld\ngenus = %lld %lld\n%lld %lld %lld %lld %lld \n",head, genus,Nbody[0],  r1 , sp, M[0], M[1],M[2]);
   // fflush(stdout);

    flag2 = 0;

#endif
    inta dims;
    dims = 0;
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.canon[space].body != nada)
            dims++;
    if ( comp != dims ){
        printf("dimensional twist! %d %d", comp,dims);
        exit(0);
    }
  
    if ( input == -1 ){
        fclose(in);
        return head;
    }

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
    if ( 300 <= input && input < 400 ){
        fclose(in);
        return comps[input - 300];
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
//    for ( i = 0 ; i < SPACE ; i++)
//        if ( f1.canon[i].body != nada)
//            printf("%d M %d\n",i, M[i]);
//    fflush(stdout);
    
    
    f1.name[buffer].value.symmetry = sy;
    
    if ( r1 > part(f1, buffer ) ){
        printf("io error\n");
        exit(0);
    }
    if ( sp >= spins(f1, buffer ) )
        sp = 0;

    for ( r = 0; r < r1 ; r++){
        fgets(inputPt,maxRead, in   );
        flag3 = 0;

        for ( space = 0; space < SPACE ; space++)
            if ( f1.canon[space].body != nada)
{
            fgets(inputPt,maxRead, in   );
            flag4 = 0;

            for ( l = 0 ; l < M[space] ;l++){
                fgets(inputPt,maxRead, in   );
                if ( ! sscanf (inputPt,"%lf\n",&( value)) ){
                    printf("--%d :%s: %d\n",ct,inputPt,space);
                    exit(0);
                };

                ct++;
                lvalue = value;
                streams(f1, buffer, sp, space)[ M[space] * r+ l ] = value;
                fgets(inputPt,maxRead, in   );

            }
            fgets(inputPt,maxRead, in   );
        }
        fgets(inputPt,maxRead, in   );
    }
    f1.name[buffer].Current[sp] = r1;
    
    
    
    if ( input == 1 ){
        fclose(in);
        return ct;
    }
    

    return 0;
#endif
}


inta tLoadEigenWeights (  calculation * c1,   field *f,char * filename, inta *ct,  division inputVectors, inta collect){
      sinc_label f1 = f->f;
    inta space,cmpl;
    FILE * in = NULL;
    in = fopen(filename, "r");
    if ( in == NULL ){
        printf("file of occupations is missing\n");
        exit(0);
    }
      bodyType body;
    inta flagLoad,stage;
    DCOMPLEX ov;
    int ms = MAXSTRING;
    char input_line[MAXSTRING];
    char input_line2[MAXSTRING];
    char * mask = input_line;
    DCOMPLEX Occ;
    char name[MAXSTRING];
    while (1){
        if (  fgets(mask,ms,in) > 0 ){
            if ( (!comment(input_line)) && (strlen(input_line) > 1) ){
                Occ = 0.;
                flagLoad = 0;
                for ( cmpl = 0; cmpl < spins(f1, inputVectors); cmpl++)
                {
                    strcpy(input_line2 , input_line);
                    Occ = tFromReadToFilename(NULL, input_line2,  name, spins(f1,eigenVectors)-1,cmpl,f1.name[inputVectors+*ct].value.title,&stage);
                    if ( cabs(Occ) > c1->rt.THRESHOLD){
                        
                            f1.name[inputVectors+*ct].Current[cmpl] = 0;
                        field f2;
                        initField(&f2);
                            calculation c2;
                            c2 = *c1;
                            c2.i.numNames = 100+ 3*inputFormat(f1, name, nullName, 2);
                            c2.i.numVectors = 0;

                            c2.i.termNumber = 0;
                            c2.rt.NLanes = 1;
                            f2.f.rt = &c2.rt;
                            f2.f.rt->phaseType = productKrylov;
                            f2.i = f->i;

                            f2.i.Iterations = 1;
                            f2.i.files = 0;
                            f2.i.filesVectorOperator = 0;
                            f2.i.qFloor = 0;
                            c2.i.Lambda = 6;
                            resetA(f2.f.rt);
                            if ( (f->i.filter/2)%2 == 0 ){
                                blockA(f2.f.rt, blockTotalVectorBlock);
                                blockA(f2.f.rt, blockTrainVectorsblock);
                                blockA(f2.f.rt, blockPermutationsblock);
                            }
                            blockA(f2.f.rt, blockCopyBlock);
                            blockA(f2.f.rt, blockParallelMultiplyblock);
                            blockA(f2.f.rt, blockParallelMatrixElementblock);
                            blockA(f2.f.rt, blockTotalVectorParallelBlock);

                            f2.i.body = inputFormat(f1,name, nullName, 100);
                        
                            f2.f.boot = noMatrices;
                            
                            ///below it *May* filter-decompose
                            f2.i.canonRank  = imax(part(f1,inputVectors+*ct), inputFormat(f1, name, nullName, 2));
                        
                            f2.i.nStates = 1;
                            for (space = 0 ;space < SPACE ; space++){
                                
                                ///previous stage
                                f2.f.canon[space] = f1.canon[space];
                                f2.f.canon[space].stream = NULL;
                                if ( f2.f.canon[space].body != nada )
                                if ( f1.canon[space].basis == SincBasisElement  )
                                {
#ifdef readHDF5
                                    f2.f.canon[space].count1Basis =((inputFormat(f1, name, nullName,200+space)));

#else
                                    f2.f.canon[space].count1Basis =((inputFormat(f1, name, nullName, 200+f2.f.canon[space].label-1)));
#endif
                                    
                                    
                                    
                                    
                                    for (body = one ; body <= f2.f.canon[space].body; body++)
                                        {
                                            f2.f.canon[space].particle[body].origin +=  +f1.canon[space].particle[body].lattice*floor(f1.canon[space].count1Basis*f1.canon[space].particle[body].anchor);
                                            ///remove previous origin
                                            
                                            
                                            f2.f.canon[space].particle[body].lattice  = f1.canon[space].particle[body].lattice * pow( f1.canon[space].count1Basis*1./f2.f.canon[space].count1Basis,f1.canon[space].particle[body].attack);
                                            ///attack lattice
                                            
                                            f2.f.canon[space].particle[body].origin -= +f2.f.canon[space].particle[body].lattice*floor(f2.f.canon[space].count1Basis*f2.f.canon[space].particle[body].anchor);
                                            ///reintroduce origin
                                        }
                                    ///new stage..
                                }
                            }
                            iModel(&c2,&f2);
                            inputFormat(f2.f, name, eigenVectors,1);
                            
                        {
                            inta sp;
                            ///filter with +2,  will filter input vector.
                            if ( (((f->i.filter/2)%2)==1)*f->f.irrep ) {
                                if ( ! allowQ(f2.f.rt,blockTotalVectorBlock)){
                                    printf("blockTotalVectorBlock Allow!\n");
                                    fflush(stdout);
                                    exit(0);
                                }

                                for ( sp = 0; sp < spins(f1, eigenVectors);sp++){
                                    f2.f.name[totalVector].Current[0] = 0;
                                    tBuildIrr(0, f2.f, f->f.irrep, eigenVectors, sp, totalVector, 0);
                                    CanonicalRankDecomposition( f2.f, NULL,totalVector, 0,eigenVectors,sp, f1.rt->TOLERANCE,f1.rt->relativeTOLERANCE, f1.rt->ALPHA,f1.rt->THRESHOLD, f1.rt->MAX_CYCLE,f1.rt->XCONDITION, part(f2.f,eigenVectors),f2.f.rt->dynamic);
                                }
                            }
                        }
                            if ( collect ){
                                xEqua(f1,inputVectors+*ct, cmpl , f2.f, eigenVectors,0);
                                flagLoad = 0;
                                if ( tSelect(f1, *ct, 0, inputVectors, 1) ) {
                                    
                                    if ( (((f->i.filter)%2)==1)*f->f.irrep ){
                                        f1.name[inputVectors+*ct].value.symmetry = tClassify( f1, inputVectors+*ct);
                                        printf("%s\tcollect-SA%d\n", name,  f1.name[inputVectors+*ct].value.symmetry);
                                        fflush(stdout);
                                    }
                                    else {
                                        f1.name[inputVectors+*ct].value.symmetry = 0;
                                    }
                                    if ( ((f->i.filter%2)==1) && (f1.name[inputVectors+*ct].value.symmetry != f->f.irrep) && f->f.irrep)
                                        flagLoad = 0;
                                    else
                                        flagLoad = 1;

                                }
                            } else {
                                xEqua(f1,inputVectors+*ct, cmpl, f2.f, eigenVectors,0);

                                if ( (((f->i.filter)%2)==1)*f->f.irrep ){
                                    f1.name[inputVectors+*ct].value.symmetry = tClassify( f1, inputVectors+*ct);
                                    printf("%s\tSA%d\n", name,  f1.name[inputVectors+*ct].value.symmetry);
                                    fflush(stdout);
                                }
                                else {
                                    f1.name[inputVectors+*ct].value.symmetry = 0;
                                }
                                if ( ((f->i.filter%2)==1) && (f1.name[inputVectors+*ct].value.symmetry != f->f.irrep) && f->f.irrep)
                                    flagLoad = 0;
                                else
                                    flagLoad = 1;

                            }
                            f1.name[inputVectors+*ct].value.stage = stage;
                            fModel(&f2.f);
                        }
                        
                        
                    }
                


                if (( *ct > f->i.nOperator && inputVectors == f1.vectorOperator )){
                    printf("maxed out buffer of states\n");
                    exit(0);
                }
                if ( cimag(Occ) != 0. ){
                    printf("fix me\n");
                    exit(1);
                }

                if ( flagLoad ) {
                    ov = 1;
                    if ( species(f1, inputVectors+*ct) == vector){
                        tScale(f1, inputVectors+*ct, creal(Occ));//error... need to scale real and complex separately!!!

                        printf("vector --");
                        fflush(stdout);
                        if ( allowQ (f1.rt, blockPrintStuffblock))
                            ov = sqrt(tMatrixElements(0,f1, inputVectors+*ct ,0,nullOverlap,0,inputVectors+*ct ,0));
                    }
                    else{
                        tScale(f1, inputVectors+*ct, fabs(creal(Occ)));
                        //setting Occ >0 so that split operators will not lose their signage via operatorSignFlag
                        
                        if ( allowQ (f1.rt, blockPrintStuffblock)){
                            
                            if ( creal(Occ) < 0. ){
                                printf("matrix -neg-");
                                ov = traceOne(f1, inputVectors+*ct, 0);

                            }else {
                                printf("matrix -pos-");
                                ov = traceOne(f1, inputVectors+*ct, 0);

                            }
                            if ( spins(f1, inputVectors+*ct)> 1 ){
                                ov = I*traceOne(f1, inputVectors+*ct, 1);
                            }
                        }
                    }
                    
                    if ( inputVectors >= f1.vectorOperator){
                        if ( allowQ (f1.rt, blockPrintStuffblock))
                            printf("Operator %f\n",cabs(ov));
                        (*ct)++;
                    }
                    else{
                        if ( allowQ (f1.rt, blockPrintStuffblock))
                            printf("Norm %f\n",(creal(ov)));
                        (*ct)++;
                    }
                }
            }
        }else {
            break;
        }
    }
    fflush(stdout);
    return 0;
};

#ifdef writeHDF5
#ifndef WRITE_FAST
/**
 *An IO solution for big systems
 *@param f1          container
 *@param filename  char*
 *@param space to write
 *@param label  the destination of the input
 *@param spin to write
 *
 *http://web.mit.edu/fwtools_v3.1.0/www/Intro/IntroExamples.html
*/
inta writeFast( sinc_label f1, char * filename, inta space, division label ,inta spin){
    
    hid_t       file;                        /* handles */
    hid_t       dataset;
    hid_t       filespace;
    hid_t       dataspace;
    hid_t       datatype;
    hid_t       memspace;
    hid_t       aid2;
    hid_t       attr1,attr2,attr3,attr4;
    hid_t       ret;
    hid_t       ret1;
    hsize_t     dims[1];                     /* dataset */

    herr_t      status, status_n;
       
    
    int canonRank = CanonicalRank(f1,label,spin),genus=1,particle,body = f1.canon[space].body ,count1 = vector1Len(f1,space);

    
    
#ifdef MODULARIZE_OUTPUT
    
    char tokens[3][MAXSTRING];
    char * stage = &*(tokens[0]);
    char * phase = &*(tokens[1]);
    char * remainder = &*(tokens[2]);
    char str[SUPERMAXSTRING];
    char destroy [SUPERMAXSTRING];
    strcpy(destroy, filename);
    const char * pstr;

    stage = strtok(destroy, "/");
    phase = strtok(NULL, ".");
    remainder = strtok(NULL, "_");

    sprintf(str , "%s.%s.%s.%d.%d", stage,phase,remainder,space,spin);
    pstr = &str[0];
    
    fflush(stdout);
    file = -1;
    
    while ( file < 0 ){
        if ( !strcmp(phase,"D") ){
            file = H5Fopen("D", H5F_ACC_RDWR, H5P_DEFAULT);
        }else {
            char fileout[MAXSTRING];
            sprintf(fileout,"%s/%s", stage, "T");
            file = H5Fopen(fileout, H5F_ACC_RDWR, H5P_DEFAULT);
        }
        sleep(rand() % SLEEP_DURATION +SLEEP_DURATION );
    }
#else
    char str[6];
    const char * pstr;
    sprintf(str,"%3d-%1d",space,spin);
    pstr = &str[0];

    /*
     * Open the file and the dataset.
     */
    file = -1;
    
    while ( file < 0 ){
        if ( !space )
            file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        else
            file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
        sleep(rand() % SLEEP_DURATION +SLEEP_DURATION );
    }
#endif
    
    
    dims[0] = canonRank*vectorLen(f1,space);
    dataspace = H5Screate_simple(1, dims, NULL);
        
    if ( ! H5Lexists(file,pstr,H5P_DEFAULT))
        dataset = H5Dcreate2(file, pstr, H5T_NATIVE_DOUBLE, dataspace,
                        H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    else
        dataset = H5Dopen2(file, pstr,H5P_DEFAULT);

    
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
              H5P_DEFAULT, streams(f1,label,spin,space));
    
    
    aid2  = H5Screate(H5S_SCALAR);
    
    if ( H5Aexists(dataset, "genus"))
            H5Adelete(dataset , "genus");

    attr1 = H5Acreate2(dataset, "genus", H5T_NATIVE_INT, aid2,H5P_DEFAULT,
                      H5P_DEFAULT);
    ret = H5Awrite(attr1, H5T_NATIVE_INT, &genus);
    H5Aclose(attr1);
    if ( H5Aexists(dataset, "canonRank"))
            H5Adelete(dataset , "canonRank");

    attr2 = H5Acreate2(dataset, "canonRank", H5T_NATIVE_INT, aid2,H5P_DEFAULT,
                      H5P_DEFAULT);
    ret = H5Awrite(attr2, H5T_NATIVE_INT, &canonRank);
    H5Aclose(attr2);

    if ( H5Aexists(dataset, "count1Basis"))
            H5Adelete(dataset , "count1Basis");

    attr3 = H5Acreate2(dataset, "count1Basis", H5T_NATIVE_INT, aid2,H5P_DEFAULT,
                      H5P_DEFAULT);
    ret = H5Awrite(attr3, H5T_NATIVE_INT, &count1);
    H5Aclose(attr3);

    if ( H5Aexists(dataset, "body"))
            H5Adelete(dataset , "body");

    attr4 = H5Acreate2(dataset, "body", H5T_NATIVE_INT, aid2,H5P_DEFAULT,
                      H5P_DEFAULT);
    ret = H5Awrite(attr4, H5T_NATIVE_INT, &body);
    H5Aclose(attr4);
    H5Sclose(aid2);


    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file);

}

#else

/**
 *An IO solution for big systems, chunked
 *@param f1          container
 *@param filename  char*
 *@param space to write
 *@param label  the destination of the input
 *@param spin to write
 *
 *http://web.mit.edu/fwtools_v3.1.0/www/Intro/IntroExamples.html
*/
inta writeFast( sinc_label f1, char * filename, inta space, division label ,inta spin){
    
    hid_t       file;                        /* handles */
    hid_t       dataset;
    hid_t       filespace;
    hid_t       dataspace;
    hid_t       datatype;
    hid_t       memspace;
    hid_t       aid2;
    hid_t       attr1,attr2,attr3,attr4;
    hid_t       ret;
    hid_t       ret1;
    hsize_t     dims[2];                     /* dataset */

    herr_t      status, status_n;
       
    
    int s,canonRank = CanonicalRank(f1,label,spin),genus=1,particle,body = f1.canon[space].body ,count1 = vector1Len(f1,space);

#ifdef MODULARIZE_OUTPUT
    
    char tokens[3][MAXSTRING];
    char * stage = &*(tokens[0]);
    char * phase = &*(tokens[1]);
    char * remainder = &*(tokens[2]);
    char str[SUPERMAXSTRING];
    char destroy [SUPERMAXSTRING];
    strcpy(destroy, filename);
    const char * pstr;

    stage = strtok(destroy, "/");
    phase = strtok(NULL, ".");
    remainder = strtok(NULL, "_");

    sprintf(str , "%s.%s.%s.%d.%d", stage,phase,remainder,space,spin);
    pstr = &str[0];
    
    fflush(stdout);
    file = -1;
    while ( file < 0 ){
        if ( !strcmp(phase,"D") ){
            file = H5Fopen("D", H5F_ACC_RDWR, H5P_DEFAULT);
        }else {
            char fileout[MAXSTRING];
            sprintf(fileout,"%s/%s", stage, "T");
            file = H5Fopen(fileout, H5F_ACC_RDWR, H5P_DEFAULT);
        }
        sleep(rand() % SLEEP_DURATION +SLEEP_DURATION );
    }
#else
    char str[6];
    const char * pstr;
    sprintf(str,"%3d-%1d",space,spin);
    pstr = &str[0];

    /*
     * Open the file and the dataset.
     */
    while ( file < 0 ){
        if ( !space )
            file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        else
            file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
        sleep(rand() % SLEEP_DURATION +SLEEP_DURATION );
    }
#endif


    dims[0] = canonRank;
    dims[1] = vectorLen(f1,space);

    /*
     * Open the file and the dataset.
     */    
    
    dataspace = H5Screate_simple(2, dims, NULL);
    //dapl!!!
    if ( ! H5Lexists(file,pstr,H5P_DEFAULT))
        dataset = H5Dcreate2(file, pstr, H5T_NATIVE_DOUBLE, dataspace,
                        H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    else
        dataset = H5Dopen2(file, pstr,H5P_DEFAULT);

    
    
    double *ptr[dims[0]];
    for ( s = 0 ; s < dims[0] ; s++)
        ptr[s] = streams(f1,label,spin,space)+s*dims[1];
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
              H5P_DEFAULT, ptr[0]);
    
    
    
    aid2  = H5Screate(H5S_SCALAR);
    if ( H5Aexists(dataset, "genus"))
            H5Adelete(dataset , "genus");
    
    attr1 = H5Acreate2(dataset, "genus", H5T_NATIVE_INT, aid2,H5P_DEFAULT,
                      H5P_DEFAULT);
    ret = H5Awrite(attr1, H5T_NATIVE_INT, &genus);
    H5Aclose(attr1);

    
    if ( H5Aexists(dataset, "canonRank"))
            H5Adelete(dataset , "canonRank");

    attr2 = H5Acreate2(dataset, "canonRank", H5T_NATIVE_INT, aid2,H5P_DEFAULT,
                      H5P_DEFAULT);
    ret = H5Awrite(attr2, H5T_NATIVE_INT, &canonRank);
    H5Aclose(attr2);

    
    if ( H5Aexists(dataset, "count1Basis"))
            H5Adelete(dataset , "count1Basis");

    
    attr3 = H5Acreate2(dataset, "count1Basis", H5T_NATIVE_INT, aid2,H5P_DEFAULT,
                      H5P_DEFAULT);
    ret = H5Awrite(attr3, H5T_NATIVE_INT, &count1);
    H5Aclose(attr3);
    
    if ( H5Aexists(dataset, "body"))
            H5Adelete(dataset , "body");



    attr4 = H5Acreate2(dataset, "body", H5T_NATIVE_INT, aid2,H5P_DEFAULT,
                      H5P_DEFAULT);
    ret = H5Awrite(attr4, H5T_NATIVE_INT, &body);
    H5Aclose(attr4);
    H5Sclose(aid2);


    H5Sclose(dataspace);
    H5Dclose(dataset);
    printf("pre-close\n");
    fflush(stdout);
    H5Fclose(file);
    printf("close\n");
    fflush(stdout);

}

#endif
#endif

#ifdef readHDF5
#ifndef READ_FAST

/**
 *An IO solution for big systems
 *@param f1          container
 *@param filename  char*
 *@param command to program, may read various attributes or file
 *0 output genus
 *1 read vector
 *2 output CanonRank
 *3 body
 *4 count1Basis
 *@param space to read from disk
 *@param label  the destination of the input
 *@param space2 to place into label
*/
inta readFast( sinc_label f1, char * filename, inta command, inta space, division label ,inta spin, inta space2){
    
    hid_t       file;                        /* handles */
    hid_t       dataset;
    hid_t       filespace;
    hid_t       attr;
    hid_t       ret;
    hid_t       memspace;
    hsize_t     dims[1];                     /* dataset */

    herr_t      status, status_n;
   
    int canonRank,genus,particle,body,count1;
    /*
     * Open the file and the dataset.
     */
#ifdef MODULARIZE_INPUT
    char tokens[3][MAXSTRING];
    char * stage = &*(tokens[0]);
    char * phase = &*(tokens[1]);
    char * remainder = &*(tokens[2]);
    char str[SUPERMAXSTRING];
    char destroy [SUPERMAXSTRING];
    strcpy(destroy, filename);
    const char * pstr;

    stage = strtok(destroy, "/");
    phase = strtok(NULL, ".");
    remainder = strtok(NULL, "_");

    sprintf(str , "%s.%s.%s.%d.%d", stage,phase,remainder,space,spin);
    pstr = &str[0];
    
    fflush(stdout);
    
    file = -1;
    while ( file < 0 ){
        if ( !strcmp(phase,"D") ){
            file = H5Fopen("D", H5F_ACC_RDWR, H5P_DEFAULT);
        }else {
            char fileout[MAXSTRING];
            sprintf(fileout,"%s/%s", stage, "T");
            file = H5Fopen(fileout, H5F_ACC_RDWR, H5P_DEFAULT);
        }
        sleep(rand() % SLEEP_DURATION +SLEEP_DURATION );
    }

#ifdef BACKWARDS
    if ( ! H5Lexists(file,pstr,H5P_DEFAULT)){
        H5Fclose(file);
        return readBackwards(f1,filename,command,space,label,spin,space2);
    }
#endif

#else
    char str[6];
    const char * pstr;
    file = -1;
    while ( file < 0 ){
        file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        sleep(rand() % SLEEP_DURATION +SLEEP_DURATION );
    }
    sprintf(str,"%3d-%1d",space,spin);
    pstr = &str[0];
#endif
    ///
        ///
        ///
    {
        dataset = H5Dopen(file, pstr, H5P_DEFAULT);
        
        if ( command == 0 ){
            attr = H5Aopen_name(dataset,"genus");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &genus);
            ///close.
            H5Aclose(attr);
            H5Dclose(dataset);
            H5Fclose(file);
            
            return genus;
        }
        
        if ( command == 2 ){
            attr = H5Aopen_name(dataset,"canonRank");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &canonRank);
            ///close.
            H5Aclose(attr);
            H5Dclose(dataset);
            H5Fclose(file);

            return canonRank;
        }
                
        if ( command == 3 ){
            attr = H5Aopen_name(dataset,"body");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &body);
            ///close.
            H5Aclose(attr);
            H5Dclose(dataset);
            H5Fclose(file);

            return body;
        }

        if ( command == 4 ){
            attr = H5Aopen_name(dataset,"count1Basis");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &count1);
            ///close.
            H5Aclose(attr);

            H5Dclose(dataset);
            H5Fclose(file);

            return count1;
        }
        attr = H5Aopen_name(dataset,"canonRank");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &canonRank);
        H5Aclose(attr);

        attr = H5Aopen_name(dataset,"genus");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &genus);
        H5Aclose(attr);

        attr = H5Aopen_name(dataset,"body");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &body);
        H5Aclose(attr);

        attr = H5Aopen_name(dataset,"count1Basis");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &count1);
        H5Aclose(attr);

        ///close.
        dims[0] = canonRank*pow(count1,body * genus );


        filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
        memspace = H5Screate_simple(1,dims,NULL);
         
        /*
        * Read dataset
        */
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,H5P_DEFAULT, streams(f1,label,spin,space2) );
        
        H5Sclose(filespace);
        
        H5Sclose(memspace);
        
        H5Dclose(dataset);

     }///spin,space

    f1.name[label].Current[spin] = canonRank;
    H5Fclose(file);

    return 0;
}

#else
/**
 *An IO solution for big systems, makes an 2D image.  Can accept older 1D images.
 *
 *@param f1          container
 *@param filename  char*
 *@param command to program, may read various attributes or file
 *0 output genus
 *1 read vector
 *2 output CanonRank
 *3 body
 *4 count1Basis
 *@param space to read from disk
 *@param label  the destination of the input
 *@param space2 to place into label
*/
inta readFast( sinc_label f1, char * filename, inta command, inta space, division label ,inta spin, inta space2){
    
    hid_t       file;                        /* handles */
    hid_t       dataset;
    hid_t       filespace;
    hid_t       attr;
    hid_t       ret;
    hid_t       memspace;
    hsize_t     dims[2];                     /* dataset */
    herr_t      status, status_n;
   
    
    
    int type2,s,canonRank,genus,particle,body,count1;
#ifdef MODULARIZE_INPUT
    char tokens[3][MAXSTRING];
    char * stage = &*(tokens[0]);
    char * phase = &*(tokens[1]);
    char * remainder = &*(tokens[2]);
    char str[SUPERMAXSTRING];
    char destroy [SUPERMAXSTRING];
    strcpy(destroy, filename);
    const char * pstr;

    stage = strtok(destroy, "/");
    phase = strtok(NULL, ".");
    remainder = strtok(NULL, "_");

    sprintf(str , "%s.%s.%s.%d.%d", stage,phase,remainder,space,spin);
    pstr = &str[0];
    
    fflush(stdout);
    file = -1;
    while ( file < 0 ){
        if ( !strcmp(phase,"D") ){
            file = H5Fopen("D", H5F_ACC_RDWR, H5P_DEFAULT);
        }else {
            char fileout[MAXSTRING];
            sprintf(fileout,"%s/%s", stage, "T");
            file = H5Fopen(fileout, H5F_ACC_RDWR, H5P_DEFAULT);
        }
        sleep(rand() % SLEEP_DURATION +SLEEP_DURATION );
    }

#ifdef BACKWARDS
    if ( ! H5Lexists(file,pstr,H5P_DEFAULT)){
        H5Fclose(file);
        return readBackwards(f1,filename,command,space,label,spin,space2);
    }
#endif
    
    
    
#else
    char str2[8];

    char str[6];
    const char * pstr;
    file = -1;
    while ( file < 0 ){
        file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        sleep(rand() % SLEEP_DURATION +SLEEP_DURATION );
    }
    sprintf(str,"%3d-2- %1d",space,spin);
    pstr = &str[0];
#endif

    {///
        ///first attempt to open -2-
        ///
        type2 = 1;
        dataset = H5Dopen(file, pstr, H5P_DEFAULT);
                
        
        if ( command == 0 ){
            attr = H5Aopen_name(dataset,"genus");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &genus);
            ///close.
            H5Aclose(attr);
            H5Dclose(dataset);
            H5Fclose(file);
            
            return genus;
        }
        
        if ( command == 2 ){
            attr = H5Aopen_name(dataset,"canonRank");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &canonRank);
            ///close.
            H5Aclose(attr);
            H5Dclose(dataset);
            H5Fclose(file);

            return canonRank;
        }
                
        if ( command == 3 ){
            attr = H5Aopen_name(dataset,"body");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &body);
            ///close.
            H5Aclose(attr);
            H5Dclose(dataset);
            H5Fclose(file);

            return body;
        }

        if ( command == 4 ){
            attr = H5Aopen_name(dataset,"count1Basis");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &count1);
            ///close.
            H5Aclose(attr);

            H5Dclose(dataset);
            H5Fclose(file);

            return count1;
        }
        attr = H5Aopen_name(dataset,"canonRank");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &canonRank);
        H5Aclose(attr);

        attr = H5Aopen_name(dataset,"genus");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &genus);
        H5Aclose(attr);

        attr = H5Aopen_name(dataset,"body");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &body);
        H5Aclose(attr);

        attr = H5Aopen_name(dataset,"count1Basis");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &count1);
        H5Aclose(attr);

        
        {
            ///close.
            dims[0] = canonRank;
            dims[1] = pow(count1, body);
             
            filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
            memspace = H5Screate_simple(2,dims,NULL);
            
         
            /*
            * Define 2D image 
            */

            double *ptr[dims[0]];
            for ( s = 0 ; s < dims[0] ; s++)
                ptr[s] = streams(f1,label,spin,space)+s*dims[1];
            /*
            * Read dataset
            */

            status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,H5P_DEFAULT, ptr[0] );
                   
         
         

        }
        
        H5Sclose(memspace);
        H5Sclose(filespace);

        H5Dclose(dataset);

     }///spin,space

    f1.name[label].Current[spin] = canonRank;
    H5Fclose(file);

    return 0;
}
#endif
#endif


#ifdef readHDF5

/**
 *An IO solution for big systems
 *@param f1          container
 *@param filename  char*
 *@param command to program, may read various attributes or file
 *0 output genus
 *1 read vector
 *2 output CanonRank
 *3 body
 *4 count1Basis
 *@param space to read from disk
 *@param label  the destination of the input
 *@param space2 to place into label
*/
inta readBackwards( sinc_label f1, char * filename, inta command, inta space, division label ,inta spin, inta space2){
    
    hid_t       file;                        /* handles */
    hid_t       dataset;
    hid_t       filespace;
    hid_t       attr;
    hid_t       ret;
    hid_t       memspace;
    hsize_t     dims[1];                     /* dataset */

    herr_t      status, status_n;
   
    int canonRank,genus,particle,body,count1;
    /*
     * Open the file and the dataset.
     */
    char str[6];
    const char * pstr;

    file = -1;
    while ( file < 0 ){
        file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        sleep(rand() % SLEEP_DURATION +SLEEP_DURATION );
    }
    sprintf(str,"%3d-%1d",space,spin);
    pstr = &str[0];

    {
        dataset = H5Dopen(file, pstr, H5P_DEFAULT);
        
        if ( command == 0 ){
            attr = H5Aopen_name(dataset,"genus");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &genus);
            ///close.
            H5Aclose(attr);
            H5Dclose(dataset);
            H5Fclose(file);
            
            return genus;
        }
        
        if ( command == 2 ){
            attr = H5Aopen_name(dataset,"canonRank");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &canonRank);
            ///close.
            H5Aclose(attr);
            H5Dclose(dataset);
            H5Fclose(file);

            return canonRank;
        }
                
        if ( command == 3 ){
            attr = H5Aopen_name(dataset,"body");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &body);
            ///close.
            H5Aclose(attr);
            H5Dclose(dataset);
            H5Fclose(file);

            return body;
        }

        if ( command == 4 ){
            attr = H5Aopen_name(dataset,"count1Basis");
            ret  = H5Aread(attr, H5T_NATIVE_INT, &count1);
            ///close.
            H5Aclose(attr);

            H5Dclose(dataset);
            H5Fclose(file);

            return count1;
        }
        attr = H5Aopen_name(dataset,"canonRank");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &canonRank);
        H5Aclose(attr);

        attr = H5Aopen_name(dataset,"genus");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &genus);
        H5Aclose(attr);

        attr = H5Aopen_name(dataset,"body");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &body);
        H5Aclose(attr);

        attr = H5Aopen_name(dataset,"count1Basis");
        ret  = H5Aread(attr, H5T_NATIVE_INT, &count1);
        H5Aclose(attr);

        ///close.
        dims[0] = canonRank*pow(count1,body * genus );


        filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
        memspace = H5Screate_simple(1,dims,NULL);
         
        /*
        * Read dataset
        */
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,H5P_DEFAULT, streams(f1,label,spin,space2) );
        
        H5Sclose(filespace);
        
        H5Sclose(memspace);
        
        H5Dclose(dataset);

     }///spin,space

    f1.name[label].Current[spin] = canonRank;
    H5Fclose(file);

    return 0;
}
#endif
