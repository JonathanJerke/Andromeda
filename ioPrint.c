/*
 *  ioPrint.c
 *
 *
 *  Copyright 2018 Jonathan Jerke and Bill Poirier.
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

void outputFormat(struct field  * f1, FILE * out, enum division output ,INT_TYPE spin){
    INT_TYPE  flag2,flag3,flag4,r,l,space, M[3];
    length(f1, output, M);
    if ( header (f1, output ) != Cube )
    {
        printf("outputFormat: WARNING,Non Cubic output\n");
    }
#if 0
    fprintf(out, "header = %d\n", header(f1, output));
    fprintf(out,"genus = %d %d %f\n", species(f1, output), bodies(f1, output),f1->sinc.d);
    fprintf(out,"%d %d %d %d %d\n", CanonicalRank(f1, output,spin ) , spin, M[0], M[1],M[2]);
#else
    fprintf(out, "header = %lld\n", header(f1, output));
    fprintf(out,"genus = %lld %d %f\n", species(f1, output), bodies(f1, output),f1->sinc.d);
    fprintf(out,"%lld %lld %lld %lld %lld\n", CanonicalRank(f1, output,spin ) , spin, M[0], M[1],M[2]);
#endif
    fprintf (out,"Tensor = {\n");
    flag2 = 0;
    
    for ( r = 0; r < CanonicalRank(f1, output, spin ) ; r++){
        if ( ! flag2 ){
            flag2 =1 ;
        }else
            fprintf (out,",\n");
        fprintf (out,"{\n");
        flag3 = 0;
        
        
        for ( space = 0; space < SPACE ; space++){
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

#if 1
INT_TYPE inputFormat(struct field * f1,char * name,  enum division buffer, enum division input){
    size_t maxRead = 4096;
    char input_line [maxRead];
    double value,lvalue;
    char * inputPt= input_line;;
    //    struct calculation * c2 = malloc( sizeof(struct calculation));
    //    initCalculation(c2);
    
    
    FILE * in = fopen(name, "r");
    if ( in == NULL ){
        printf("failed to load %s\n", name);
        exit(0);
    }
    //    broke = readInput(c2,in);
    //    finalizeInit(c2);
    
    //  getline(&inputPt,&maxRead, in   );
    //  getline(&inputPt,&maxRead, in   );
    getline(&inputPt,&maxRead, in   );
    INT_TYPE ct = 0;
#ifndef MKL
    INT_TYPE head, genus;
    double D;
    INT_TYPE M[3],r1,r,space,flag2,flag3,flag4,l,sp,Nbody;
    
    sscanf(inputPt, "header = %d", &head );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "genus = %d %d %lf", &genus , &Nbody,&D);
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "%d %d %d %d %d", &r1, &sp, &(M[0]), &(M[1]),&(M[2]) );
    getline(&inputPt,&maxRead, in   );
    flag2 = 0;
    printf("header = %d\ngenus = %d %d %f\n%d %d %d %d %d \n",head, genus,Nbody, D,  r1 , sp, M[0], M[1],M[2]);
#else
    INT_TYPE head, genus;
    double D;
    INT_TYPE M[3],r1,r,space,flag2,flag3,flag4,l,sp,Nbody;
    
    sscanf(inputPt, "header = %lld", &head );
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "genus = %lld %lld %lf", &genus , &Nbody,&D);
    getline(&inputPt,&maxRead, in   );
    sscanf(inputPt, "%lld %lld %lld %lld %lld", &r1, &sp, &(M[0]), &(M[1]),&(M[2]) );
    getline(&inputPt,&maxRead, in   );
    flag2 = 0;
    //   printf("header = %lld\ngenus = %lld %lld %f\n%lld %lld %lld %lld %lld \n",head, genus,Nbody, D,  r1 , sp, M[0], M[1],M[2]);
    
#endif
    if ( input == 0 ){
        fclose(in);
        return genus;
    }
    if ( input == 3 ){
        fclose(in);
        return Nbody;
    }
    
    if ( input == 2 ){
        fclose(in);
        return r1;
    }
    
    if ( input == 4 ){
        fclose(in);
        if( Nbody == 2 )
            return sqrt(M[0]);
        else if ( Nbody == 1 )
            return M[0];
        else if ( Nbody == 3 )
            return pow(M[0],1./3.);
        else if ( Nbody == 4 )
            return pow(M[0],1./4.);
    }
    
//    ct = 0;
//    while ( 1 ){
//        if ( EOF == getline(&inputPt,&maxRead, in   ))
//            break;
//        if ( strstr(inputPt, "}") || strstr(inputPt, "{") || strstr(inputPt, ",") ){
//            //skip
//        }else  if ( sscanf (inputPt,"%lf\n",&( value))  ){
//            space = (ct / M[0]) % 3;
//            r = ( ct / ( M[0]*3 ) ) % 3;
//            l = ct % M[0];
//            streams(f1, buffer, sp, space)[ M[space] * r+ l ] = value;
//            ct++;
//        } else if ( sscanf (inputPt,"-%lf\n",&( value)) ){
//            space = (ct / M[0]) % 3;
//            r = ( ct / ( M[0]*3 ) ) % 3;
//            l = ct % M[0];
//            streams(f1, buffer, sp, space)[ M[space] * r+ l ] = value;
//            ct++;
//  //          if ( value > 0 ){
//              printf("wtf\n");
//     //       }
//
//
//        }
//       // printf("%lld %s\n", ct, inputPt);
//    }
//
//    printf("%lld : %lld\n", ct, 3*M[0]*r1);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    for ( r = 0; r < r1 ; r++){
//        if ( ! flag2 ){
//            flag2 =1 ;
//        }else
//            getline(&inputPt,&maxRead, in   );
        getline(&inputPt,&maxRead, in   );
        flag3 = 0;

        for ( space = 0; space < SPACE ; space++){
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
    
    
    
    INT_TYPE length ;
    if ( genus == 1 )//vector!!
        length = M[0];
    else if ( genus == 2 )
        length = sqrt(M[0]);
    else if ( genus == 3 )
        length = pow(M[0],1./3.);
    else if ( genus == 4 )
        length = pow(M[0],1./4.);
    
    else {
        printf("wrong load genus\n");
        exit(0);
    }
    
    
    //    if ( Nbody == 1 ){
    //        exit(0);
    //        //  tOneBand(f1, buffer,length, D,input);
    //    } else if ( Nbody == 2 ) {
    //        tTwoBand(f1, buffer,length, D,input);
    //    } else if ( Nbody == 3 ) {
    //        tThreeBand(f1, buffer,length, D,input);
    //    } else if ( Nbody == 4 ) {
    //        tFourBand(f1, buffer,length, D,input);
    //    }
    fclose(in);
    return 0;
}
#endif



INT_TYPE printOutput ( struct field * f1,INT_TYPE number){
    INT_TYPE i,NC,numC;
    char printCom[MAXSTRING];
    FILE * printOut = stdout;
    double radius,res;
    assignCores(f1,1);
    if ( header == NULL ){
        printf("file\n");
        exit(0);
    }
    size_t ms = MAXSTRING,read;
    char input_line[MAXSTRING];
    char * mask = input_line;
    char name[MAXSTRING];
    if ( f1->mem1->rt->printFlag % 2 ){
    
    }else {
        return 0;
    }
    numC = f1->mem1->rt->monteCarlo;
    NC = f1->mem1->rt->samples;
    //Rs^3 4pi/3 = Volme / Ne
    double Rs = f1->sinc.N1*f1->sinc.d*pow(3./4./pi/f1->body,0.333333333);
    
    {
#ifdef OMP
#pragma omp parallel for private (i,radius,res) schedule(dynamic,1)
#endif
        for ( i=0 ; i <= NC ; i++ )
        {
            radius = f1->sinc.N1*f1->sinc.d*(i*1./NC);
            if ( f1->mem1->rt->runFlag == 7 )
                radius *= 3;
            else if (  ! i )
                continue;
            res =  tComputeRadialPlot( f1 ,number,f1->mem1->rt->runFlag,radius ,numC);
            if ( f1->mem1->rt->runFlag == 7 ){
                res *= pow(f1->sinc.N1*f1->sinc.d,3.)/4./pi;
            }
            else {
            }
            fprintf(printOut,"%d,%d,%d,%d,%12.6f,%12.6f,%12.16f\n",bodies(f1,eigenVectors),(f1->mem1->rt->runFlag)%2,(f1->mem1->rt->runFlag/2)%2, (f1->mem1->rt->runFlag/4)%2,radius,radius/Rs,res);
            fflush(printOut);
        }
    }
    return NC;
}

INT_TYPE printVectorOutput ( struct field * f1,INT_TYPE number){
    INT_TYPE i,j,k,l,NC,numC;
    char printCom[MAXSTRING];
    FILE * printOut = stdout;
    double displacement[3],res;
    assignCores(f1,1);
    //  FILE * header = fopen(f1->mem1->fileList,"r");
    if ( header == NULL ){
        printf("file\n");
        exit(0);
    }
    size_t ms = MAXSTRING,read;
    char input_line[MAXSTRING];
    char * mask = input_line;
    char name[MAXSTRING];
    //  getline(&mask,&ms,header);
    if ( (f1->mem1->rt->printFlag/2)% 2  ){
    }
    else {
        return 0;
    }
    numC = f1->mem1->rt->monteCarlo;
    NC = f1->mem1->rt->samples;
    //   fclose(header);
    {
#ifdef OMP
#pragma omp parallel for private (l,i,j,k,displacement,res) schedule(dynamic,1)
#endif
        
        for ( l = 0; l < NC*NC*NC ; l++)
        {
            i = (l % NC) - (NC-1)/2;
            j = ((l/NC) % NC) - (NC-1)/2;
            k = ((l/(NC*NC)) % NC) - (NC-1)/2;

            displacement[0] = f1->sinc.N1*f1->sinc.d*(i*(1./NC));
            displacement[1] = f1->sinc.N1*f1->sinc.d*(j*(1./NC));
            displacement[2] = f1->sinc.N1*f1->sinc.d*(k*(1./NC));

#ifdef Bill
            if ( i != 0 )
                continue;
#endif
            
            res =  tComputeVectorPlot( f1 ,number,f1->mem1->rt->runFlag,displacement ,numC);
            fprintf(printOut,"%d,%d,%d,%d,%12.6f,%12.6f,%12.6f,%d,%d,%d,%12.16f\n",bodies(f1,eigenVectors),(f1->mem1->rt->runFlag)%2,(f1->mem1->rt->runFlag/2)%2, (f1->mem1->rt->runFlag/4)%2, displacement[0],displacement[1],displacement[2],i,j,k,res);
            fflush(printOut);
        }
    }
    return NC;
}



double tComputeRadialPlot(struct field * f1,INT_TYPE number,  INT_TYPE class,  double radius,INT_TYPE numC )
{
    double res= 0.,err=1.;

    INT_TYPE N1 = f1->sinc.N1;
    INT_TYPE N12 = (N1-1)/2;
    double xl[5],xu[5];
    INT_TYPE space,num=0;

    
    
    
        struct fieldArray fA ;
        fA.f1 = f1;
        fA.number = number;
        
        xl[0] = 0.;
        xu[0] = pi;
        
        xl[1] = 0.;
        xu[1] = 2*pi;
        
        
        
        fA.class = class;
        
        
        
        
        space = 0;
        if ( class % 2 ){
            
            //        xl[space+2] = - pi/f1->sinc.d;
            //        xu[space+2] =   pi/f1->sinc.d;
            xl[space+2] = - N12*f1->sinc.d;
            xu[space+2] =   N12*f1->sinc.d;
            
        }
        else {
            
            xl[space+2] = - N12*f1->sinc.d;
            xu[space+2] =   N12*f1->sinc.d;
            
        }
        space = 1;
        if ( (class/2) % 2 ){
            
            //        xl[space+2] = - pi/f1->sinc.d;
            //        xu[space+2] =   pi/f1->sinc.d;
            xl[space+2] = - N12*f1->sinc.d;
            xu[space+2] =   N12*f1->sinc.d;
            
        }
        else {
            
            xl[space+2] = - N12*f1->sinc.d;
            xu[space+2] =   N12*f1->sinc.d;
            
        }
        space = 2;
        if ( (class/4) % 2 ){
            
            //        xl[space+2] = - pi/f1->sinc.d;
            //        xu[space+2] =   pi/f1->sinc.d;
            xl[space+2] = - N12*f1->sinc.d;
            xu[space+2] =   N12*f1->sinc.d;
            
        }
        else {
            
            xl[space+2] = - N12*f1->sinc.d;
            xu[space+2] =   N12*f1->sinc.d;
            
        }
        
        
        INT_TYPE rank = 0,dim = 5;
    
        if ( bodies (f1, eigenVectors ) == one )
            dim = 2;
        else
            dim  = 5;
#ifdef OMP
        rank = omp_get_thread_num();
        gsl_monte_function G;
        G.f = &evaluateDensityBracket;
        G.dim = dim;
        G.params = &fA;
        
        INT_TYPE iter = 0;
        const gsl_rng_type * T;
        gsl_rng  * r;
        ((struct fieldArray*)(G.params))->radius[rank] = radius;
        gsl_monte_vegas_state * s;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        s = gsl_monte_vegas_alloc(G.dim);
        
        gsl_monte_vegas_integrate ( & G , xl, xu, dim , numC, r,s,&res, &err );
        do
        {
            gsl_monte_vegas_integrate ( & G , xl, xu, dim , numC, r,s,&res, &err );
            //    printf("%f %f %f\n", res,err, gsl_monte_vegas_chisq (s));
        } while (fabs ((gsl_monte_vegas_chisq (s) - 1.0) > 0.1 || max(1.,fabs(res)) > f1->mem1->rt->TARGET)&& num++ < 100);
        
        gsl_monte_vegas_free(s);
        gsl_rng_free(r);
#else
        
        //    {
        //        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
        //
        //        gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
        //                                   &res, &err);
        //        display_results ("vegas warm-up", res, err);
        //
        //        printf ("converging...\n");
        //
        //        do
        //        {
        //            gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
        //                                       &res, &err);
        //            printf ("result = % .6f sigma = % .6f "
        //                    "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        //        }
        //        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
        //
        //        display_results ("vegas final", res, err);
        //
        //        gsl_monte_vegas_free (s);
        //    }
#endif

        return res;
}


double tComputeVectorPlot(struct field * f1,INT_TYPE number,  INT_TYPE class,  double *displacement,INT_TYPE numC )
{
    double res= 0.,err=1.;
    
    INT_TYPE N1 = f1->sinc.N1;
    INT_TYPE N12 = (N1-1)/2;
    double xl[5],xu[5];
    INT_TYPE space,num=0;
    
    struct fieldArray fA ;
    fA.f1 = f1;
    fA.number = number;
    fA.class = class;
    space = 0;
    
    if ( bodies(f1, eigenVectors ) == one ){
        res = evaluateVectorBracket(displacement, 3,(void*)(&fA));
    }
    else {
        
        
        
        if ( class % 2 ){
            xl[space] = - N12*f1->sinc.d;
            xu[space] =   N12*f1->sinc.d;
            //
            //        xl[space] = - pi/f1->sinc.d;
            //        xu[space] =   pi/f1->sinc.d;
        }
        else {
            
            xl[space] = - N12*f1->sinc.d;
            xu[space] =   N12*f1->sinc.d;
            
        }
        space = 1;
        if ( (class/2) % 2 ){
            xl[space] = - N12*f1->sinc.d;
            xu[space] =   N12*f1->sinc.d;
            
            //        xl[space] = - pi/f1->sinc.d;
            //        xu[space] =   pi/f1->sinc.d;
        }
        else {
            
            xl[space] = - N12*f1->sinc.d;
            xu[space] =   N12*f1->sinc.d;
            
        }
        space = 2;
        if ( (class/4) % 2 ){
            xl[space] = - N12*f1->sinc.d;
            xu[space] =   N12*f1->sinc.d;
            
            //        xl[space] = - pi/f1->sinc.d;
            //        xu[space] =   pi/f1->sinc.d;
        }
        else {
            
            xl[space] = - N12*f1->sinc.d;
            xu[space] =   N12*f1->sinc.d;
            
        }
        
        
        
        INT_TYPE rank = 0,dim = 3;
#ifndef APPLE
#ifdef OMP
        rank = omp_get_thread_num();
#else
        rank = 0;
#endif
        gsl_monte_function G;
        G.f = &evaluateVectorBracket;
        G.dim = 3;
        G.params = &fA;
        
        INT_TYPE iter = 0;
        const gsl_rng_type * T;
        gsl_rng  * r;
        ((struct fieldArray*)(G.params))->displacement[rank][0] = displacement[0];
        ((struct fieldArray*)(G.params))->displacement[rank][1] = displacement[1];
        ((struct fieldArray*)(G.params))->displacement[rank][2] = displacement[2];
        gsl_monte_vegas_state * s;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        s = gsl_monte_vegas_alloc(G.dim);
        
#ifdef BILL
        double x0[3];
        x0[0] = 0.;
        x0[1] = 0.;
        x0[2] = 0.;
        res = evaluateVectorBracket(x0,3,&fA );
#else
        
        gsl_monte_vegas_integrate ( & G , xl, xu, dim , numC, r,s,&res, &err );
        do
        {
            gsl_monte_vegas_integrate ( & G , xl, xu, dim , numC, r,s,&res, &err );
            // printf("%f %f %f\n", res,err, gsl_monte_vegas_chisq (s));
        } while ((fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.1 || err/max(1.,fabs(res))  > f1->mem1->rt->TARGET) && num++ < 100);
        
        gsl_monte_vegas_free(s);
        gsl_rng_free(r);
#endif
#else
        
        //    {
        //        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
        //
        //        gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
        //                                   &res, &err);
        //        display_results ("vegas warm-up", res, err);
        //
        //        printf ("converging...\n");
        //
        //        do
        //        {
        //            gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
        //                                       &res, &err);
        //            printf ("result = % .6f sigma = % .6f "
        //                    "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        //        }
        //        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
        //
        //        display_results ("vegas final", res, err);
        //
        //        gsl_monte_vegas_free (s);
        //    }
#endif
    }
    return res;
}

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
        }else
        {
                double periodicPosition;
                periodicPosition = coordinates[space];
                for ( h = 0 ; h < N1 ; h++)
                    pt[space][h] = periodicSinc(lattice,periodicPosition-lattice*(h-N12),N1)/sqrt(lattice);

        }
            
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

double evaluateDensityBracket( double x [], size_t dim , void * params ){
    double y[3],sum=0.;;
    INT_TYPE space,info,cl,cmpl;
    Stream_Type *pt[SPACE];
    INT_TYPE rank = 0;
#ifdef OMP
    rank = omp_get_thread_num();
#endif
    enum division wavefunction,product = squareVector;
    struct field *f1 = ((struct fieldArray *)(params))->f1;
    double radius = ((struct fieldArray*)(params))->radius[rank];
    INT_TYPE class = ((struct fieldArray*)(params))->class;
    INT_TYPE number = ((struct fieldArray*)(params))->number;

    if ( dim == 5 ){
        
        {
            {
                for ( space = 0; space < SPACE ; space++)
                    pt[space] = streams(f1, diagonal1VectorA,rank,space);
                y[0] = x[0+2];
                y[1] = x[1+2];
                y[2] = x[2+2];
                
                f1->sinc.tulip[diagonal1VectorA].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
                
                
                for ( space = 0; space < SPACE ; space++)
                    pt[space] = streams(f1, diagonal1VectorB,rank,space);
                y[0] = cos(x[0])*           radius+x[0+2];
                y[1] = sin(x[0])*cos(x[1])* radius+x[1+2];
                y[2] = sin(x[0])*sin(x[1])* radius+x[2+2];
                
                f1->sinc.tulip[diagonal1VectorB].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
                
                
                
                for ( cmpl = 0; cmpl < 2 ; cmpl++)
                for ( wavefunction = eigenVectors ; wavefunction < eigenVectors+number ; wavefunction++){
//                    if ( f1->sinc.tulip[printOperator].linkNext != nullName || CanonicalRank(f1, printOperator, 0) || CanonicalRank(f1, printOperator, 1)){
//                        tEqua(f1, product,rank, wavefunction,);
//                        tHXpX(rank, f1, printOperator, 0, 1., 0., product, f1->mem1->rt->CONVERGENCE, part(f1,product));
//                    }else {
//                        product = wavefunction;
//                    }

                    
                    if ( bodies ( f1, wavefunction ) == two ){
                        //oneVector :  basisRank  and ONE
                        tMultiplyMP(rank, &info,f1, 1.0, -1, oneVector , rank, 'N', wavefunction, cmpl, 'N', diagonal1VectorB, rank);
                        sum+= sqr(tMultiplyMP(rank,&info, f1, 1.0, -1, nullVector , 0, 'T', oneVector, rank, 'N', diagonal1VectorA, rank));
                    }    if ( bodies ( f1, wavefunction ) == three ){
                        //oneVector : basisRank and BODY - TWO
                        //diagonal2VectorA : 1 and TWO
                        f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
                        tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
                        tMultiplyMP(rank, &info,f1, 1.0, -1, oneVector , rank, 'N', wavefunction, cmpl, 'N', diagonal2VectorA, rank);
                        sum += tMultiplyMP(rank, &info,f1, 1.0, -1, nullVector , 0, 'T', oneVector, rank, 'N', oneVector, rank);
                    }   else   if ( bodies ( f1, wavefunction ) == four ){
                        //oneVector : basisRank and BODY - TWO
                        //diagonal2VectorA : 1 and TWO
                        f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
                        tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
                        tMultiplyMP(rank, &info,f1, 1.0, -1, twoVector , rank, 'N', wavefunction, cmpl, 'N', diagonal2VectorA, rank);
                        sum += tMultiplyMP(rank,&info, f1, 1.0, -1, nullVector , 0, 'T', twoVector, rank, 'N', twoVector, rank);
                    }
                }
            }
        }
    }
    else if (dim == 2 ){
        {
            for ( space = 0; space < SPACE ; space++)
                pt[space] = streams(f1, diagonal1VectorA,rank,space);
            y[0] = cos(x[0])*           radius;
            y[1] = sin(x[0])*cos(x[1])* radius;
            y[2] = sin(x[0])*sin(x[1])* radius;
            
            f1->sinc.tulip[diagonal1VectorA].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
            for ( cmpl = 0; cmpl < 2 ; cmpl++)

            for ( wavefunction = eigenVectors ; wavefunction < eigenVectors+number ; wavefunction++){
//                if ( f1->sinc.tulip[printOperator].linkNext != nullName || CanonicalRank(f1, printOperator, 0) || CanonicalRank(f1, printOperator, 1)){
//                    tHXpX(rank, f1, printOperator, 0, 1., 0., product, f1->mem1->rt->CONVERGENCE, part(f1,product));
//                }else {
//                    product = wavefunction;
//                }
                
                sum += sqr(tMultiplyMP(rank, &info, f1, 1.0, -1, nullVector, rank, 'T', wavefunction, cmpl, 'N', diagonal1VectorA, rank));
            }
        }
    }
    if( ! class)
        return sqr(radius)*sin(x[0])*sum;
    else
        return sin(x[0])*sum;

}


double evaluateVectorBracket( double x [], size_t dim , void * params ){
    double y[3],sum=0.;;
    INT_TYPE space,info,cl,cmpl;
    Stream_Type *pt[SPACE];
    INT_TYPE rank = 0;
#ifdef OMP
    rank = omp_get_thread_num();
#endif
    double *displacement = ((struct fieldArray*)(params))->displacement[rank];
    enum division wavefunction,product = squareVector;
    struct field *f1 = ((struct fieldArray *)(params))->f1;
    INT_TYPE class = ((struct fieldArray*)(params))->class;
    INT_TYPE number = ((struct fieldArray*)(params))->number;

    if ( dim == 3 ){
        
        {
            {
                for ( space = 0; space < SPACE ; space++)
                    pt[space] = streams(f1, diagonal1VectorA,rank,space);
                y[0] = x[0];
                y[1] = x[1];
                y[2] = x[2];
                
                f1->sinc.tulip[diagonal1VectorA].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
                
                
                for ( space = 0; space < SPACE ; space++)
                    pt[space] = streams(f1, diagonal1VectorB,rank,space);
                y[0] = displacement[0]+x[0];
                y[1] = displacement[1]+x[1];
                y[2] = displacement[2]+x[2];
                
                f1->sinc.tulip[diagonal1VectorB].Current[rank] = tFillBasis(pt, y,class,f1->sinc.N1,f1->sinc.d);
                
                for ( cmpl = 0; cmpl < 2 ; cmpl++)

                for ( wavefunction = eigenVectors ; wavefunction < eigenVectors+number ; wavefunction++){
//                    if ( f1->sinc.tulip[printOperator].linkNext != nullName || CanonicalRank(f1, printOperator, 0) || CanonicalRank(f1, printOperator, 1)){
//                        tHXpX(rank, f1, printOperator, 0, 1., 0., product, f1->mem1->rt->CONVERGENCE, part(f1,product));
//                    }else {
//                        product = wavefunction;
//                    }

                    if ( bodies ( f1, wavefunction ) == one ){
                        sum += sqr(tMultiplyMP(rank, &info,f1, 1.0, -1, nullVector , 0, 'T', wavefunction, cmpl, 'N', diagonal1VectorA, rank));
                    } else
                    if ( bodies ( f1, wavefunction ) == two ){
                        //oneVector :  basisRank  and ONE
                        tMultiplyMP(rank, &info,f1, 1.0, -1, oneVector , rank, 'N', wavefunction, cmpl, 'N', diagonal1VectorB, rank);
                        sum+= sqr(tMultiplyMP(rank,&info, f1, 1.0, -1, nullVector , 0, 'N', oneVector, rank, 'N', diagonal1VectorA, rank));
                    }  else   if ( bodies ( f1, wavefunction ) == three ){
                        //oneVector : basisRank and BODY - TWO
                        //diagonal2VectorA : 1 and TWO
                        f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
                        tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
                        tMultiplyMP(rank, &info,f1, 1.0, -1, oneVector , rank, 'N', wavefunction, cmpl, 'N', diagonal2VectorA, rank);
                        sum += tMultiplyMP(rank, &info,f1, 1.0, -1, nullVector , 0, 'T', oneVector, rank, 'N', oneVector, rank);
                    }   else   if ( bodies ( f1, wavefunction ) == four ){
                        //oneVector : basisRank and BODY - TWO
                        //diagonal2VectorA : 1 and TWO
                        f1->sinc.tulip[diagonal2VectorA].Current[rank] = 0;
                        tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
                        tMultiplyMP(rank, &info,f1, 1.0, -1, twoVector , rank, 'N', wavefunction, cmpl, 'N', diagonal2VectorA, rank);
                        sum += tMultiplyMP(rank,&info, f1, 1.0, -1, nullVector , 0, 'T', twoVector, rank, 'N', twoVector, rank);
                    }
                }
            }
        }
    }

    return sum/pow(f1->sinc.d,SPACE);
    
}

INT_TYPE tLoadEigenWeights (struct calculation * c1, char * filename){
    INT_TYPE ct = 0,ct2,number,bod,class,weight,cmpl;
    FILE * in = fopen(filename, "r");
    struct field *f1 = &c1->i.c;
    if ( in == NULL ){
        printf("file of occupations is missing\n");
        exit(0);
    }
    size_t ms = MAXSTRING,read,si;
    char input_line[MAXSTRING];
    char str0[MAXSTRING];
    char * mask = input_line;
    double Occ,Ceg;
    char name[MAXSTRING];
    while (1){
        if (  getline(&mask,&ms,in) > 0 ){
            if ( (!comment(input_line)) && (strlen(input_line) > 1) )  {
                si = sscanf ( input_line, "\"%d\",%d,%d,%lf",&str0,&number,&class, &Occ );
                read = (4== si);
                if ( read && fabs(Occ) > 1e-5){
                    tClear(&c1->i.c, eigenVectors+ct);
                    tClear(&c1->i.c, totalVector);
                    for ( cmpl = 0; cmpl < spins(&c1->i.c, eigenVectors); cmpl++)
                    {
                        sprintf(name,"%s.%d.eigen-%d.%d_mac",c1->cycleName,number,class,cmpl);
                        inputFormat(&c1->i.c, name, totalVector, 1);
                        if ( part(f1, eigenVectors + ct ) >= CanonicalRank(f1, totalVector, cmpl)){
                            tEqua(f1, eigenVectors + ct, cmpl, totalVector, cmpl);
                            printf("%d\t%s\t%f\n", ct, name, fabs(Occ));
                    }
                        else{
                            tScale(&c1->i.c, totalVector, sqrt(fabs(Occ/c1->i.c.Ne))/magnitude(&c1->i.c, totalVector));

                            tCycleDecompostionOneMP(0, f1,totalVector, cmpl, eigenVectors + ct, cmpl, f1->mem1->rt->vCANON, part(f1,eigenVectors + ct), -1);
                            printf("%d\t%s\t%f -> %f\n", ct, name, fabs(Occ),distanceOne(0, f1, totalVector, cmpl, eigenVectors+ct, cmpl));
                        }
                    }
                    tScale(&c1->i.c, eigenVectors+ct, sqrt(fabs(Occ/c1->i.c.Ne))/magnitude(&c1->i.c, eigenVectors+ct));
                    c1->i.c.sinc.tulip[eigenVectors+ct].symmetry = class;
                    ct++;
                }
                }
                if ( ct > c1->i.nStates ){
                    printf("maxed out buffer of states\n");
                    exit(0);
                }
            
        }else {
            break;
        }
    }
    
    return ct;
};



INT_TYPE tLoadEigenWeightsWithConstraints (struct calculation * c1, char * filename, char * constraints ){
    INT_TYPE ct = 0,ct2,number,bod,class,weight,cmpl,space;
    FILE * inC;
    inC = NULL;
    inC = fopen(constraints, "r");

    size_t ms = MAXSTRING,read;
    char input_line[MAXSTRING];
    char * mask = input_line;
    double Occ,Ceg,*pt[SPACE];
    for ( space = 0; space < SPACE ; space++)
        pt[space] = streams(&c1->i.c, diagonal1VectorA,0,space);

    char name[MAXSTRING];
    double xx[4][3],x,y,z;
    INT_TYPE jj=0,jjj,info;
    if ( inC != NULL) {
        while (1){
            if (  getline(&mask,&ms,inC) > 0 ){
                if ( (!comment(input_line)) && (strlen(input_line) > 1) )  {
                    read = (3==sscanf ( input_line, "%lf %lf %lf", &x,&y,&z));
                    if ( read ){
                        printf("\n*add contraint %f %f %f\n", x,y,z);
                        xx[jj][0] = x;
                        xx[jj][1] = y;
                        xx[jj][2] = z;
                        jj++;
                        if( jj == 5 ||  4 < bodies(&c1->i.c, eigenVectors)+jj ){
                            printf("cannot specify that many constraints\n");
                            exit(0);
                        }
                    }
                }
            }else {
                break;
            }
        }
        fclose(inC);
    }
    FILE * in;
    in = NULL;
    in = fopen(filename, "r");
    if ( in == NULL )
        return 0;
    
    enum division highBall = highBallVector + ( 4 - ( c1->i.c.body + jj) ) ;
    
    if ( highBall < highBallVector || highBall > highBallVector4 ){
        printf("oops with highBall!");
        exit(0);
    }
    
    while (1){
        if (  getline(&mask,&ms,in) > 0 ){
            if ( (!comment(input_line)) && (strlen(input_line) > 1) )  {
                read = (4==sscanf ( input_line, "%d %d %d %lf", &ct2,&number,&class, &Occ ));
                if ( read && fabs(Occ) > 1e-5){
                    tClear(&c1->i.c, highBall);
                    for ( cmpl = 0; cmpl < 2; cmpl++)
                    {
                        sprintf(name,"%s.%d.eigen-%d.%d_mac",c1->cycleName,number,class,cmpl);
                        printf("%d\t%s\t%f\n", ct, name, fabs(Occ));
                        inputFormat(&c1->i.c, name, highBall, 1);
                    }
                    tScale(&c1->i.c, highBall, sqrt(fabs(Occ/c1->i.c.Ne))/magnitude(&c1->i.c, highBall));
                    

                    for ( jjj = 0; jjj < jj ; jjj++){
                        c1->i.c.sinc.tulip[diagonal1VectorA].Current[0] = tFillBasis(pt, xx[jjj],c1->rt.runFlag,c1->i.c.sinc.N1,c1->i.c.sinc.d);

                        for ( cmpl = 0; cmpl < 2; cmpl++)
                        {
                            tMultiplyMP(0, &info,&c1->i.c, 1.0, -1, highBall+1+jjj , cmpl, 'N', highBall+jjj, cmpl, 'N', diagonal1VectorA, 0);
                            
                        }
                        
                    }
                    
                    if ( bodies(&c1->i.c, eigenVectors+ct ) == bodies (&c1->i.c, highBall+jj) ){
                        tEquals(&c1->i.c, eigenVectors+ct, highBall+jj);
                        printf(">>>%f\n", magnitude(&c1->i.c, eigenVectors+ct));
                    } else {
                        printf("fix up body counts!");
                        exit(0);
                    }
                    ct++;
                }
            }
            if ( ct > c1->i.nStates ){
                printf("maxed out buffer of states\n");
                exit(0);
            }
            
        }else {
            break;
        }
    }
    fclose(in);
    return ct;
}

