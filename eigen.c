/*
 *  eigen.c
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


#include "eigen.h"

double vale ( struct sortClass * f ){
    INT_TYPE * mm = f->mmm+f->i*SPACE*2 ;
    double value=1.;
    INT_TYPE space;
    for ( space = 0; space < SPACE ; space++)
        if ( f->n1[space] ){
            value *=  f->str[space][mm[2*space] + mm[2*space+1]*f->n1[space]];
            //printf("%f -- %d %d\n", value, space );
        }
    return value;
}

//double yale ( struct sortClass * f ){
//    INT_TYPE * mm = f->mmm+f->i*6 ;
//    return mm[3]+mm[4]*f->nG + mm[5]*f->nG*f->nG;
//}


int sortxComp (const void * elem1, const void * elem2)
{
    double* f = ((double*)elem1);
    double* s = ((double*)elem2);
    double valueF,valueS;
    valueF = f[0];
    valueS = s[0];
    if (valueF > valueS) return  -1;
    if (valueF < valueS) return 1;
    return 0;
}




//int sortx2Comp (const void * elem1, const void * elem2)
//{
//    double ff=1.,ss=1.,fs=1.,sf;
//   // (f+s)**2
//    ff = uw( elem1,elem1);
//    sf = uw(elem2,elem1);
//    fs = uw(elem1,elem2);
//    ss = uw(elem2,elem2);
//
//    if (ff +sf+fs > ss) return  -1;
//    if (ff +sf+fs < ss) return 1;
//    return 0;
//}



int sortComp (const void * elem1, const void * elem2)
{
    struct sortClass* f = ((struct sortClass*)elem1);
    struct sortClass* s = ((struct sortClass*)elem2);
    double valueF,valueS;
    valueF = vale(f);
    valueS = vale(s);
    if (valueF > valueS) return  1;
    if (valueF < valueS) return -1;
    return 0;
}

//int sort2Comp (const void * elem1, const void * elem2)
//{
//    struct sortClass* f = ((struct sortClass*)elem1);
//    struct sortClass* s = ((struct sortClass*)elem2);
//    INT_TYPE valueF,valueS;
//    valueF = yale(f);
//    valueS = yale(s);
//    if (valueF > valueS) return  1;
//    if (valueF < valueS) return -1;
//    return 0;
//}

INT_TYPE tBoot1Construction(struct calculation * c1, struct sinc_label f1, enum division eigen){
    assignCores(f1,1);
    enum bodyType bootBodies = f1.rose[0].body;
    INT_TYPE n1[SPACE];
    length1(f1,n1);
    INT_TYPE N1,rank;
    INT_TYPE space,N2,i,ii,iii,iv,v;
    double * ar,*w;
    DCOMPLEX * arc ;
    INT_TYPE cmplFlag = spins(f1,kinetic)-1;

    
    
    f1.tulip[Iterator].linkNext = kinetic1;
    f1.tulip[kinetic1].linkNext = kinetic2;
    f1.tulip[kinetic2].linkNext = kinetic3;
    f1.tulip[kinetic3].linkNext = kinetic4;
    f1.tulip[kinetic4].linkNext = kinetic5;
    f1.tulip[kinetic5].linkNext = kinetic6;

    if (bootBodies==two){
        f1.tulip[kinetic2].linkNext = nullName;
    }
    else if ( bootBodies == three ){
        f1.tulip[kinetic3].linkNext = nullName;
    }    else if ( bootBodies == four ){
        f1.tulip[kinetic4].linkNext = nullName;
    }    else if ( bootBodies == five ){
        f1.tulip[kinetic5].linkNext = nullName;
    }    else if ( bootBodies == six ){
        f1.tulip[kinetic6].linkNext = nullName;
    }

    
    
    
    
    
    
    
    for ( space = 0 ;space < SPACE ; space++)
        if( f1.rose[space].body != nada){
            
            N1 = n1[space];
            INT_TYPE Nx = N1;//imin(N1,c1->i.bootRestriction);
            struct name_label u = f1.tulip[canonicalBuffersBM];
            enum division em;
            N2 = N1*N1;
            myZero(f1,canonicalBuffersBM,0);
            ar = myStreams(f1, canonicalBuffersBM, 0);
            INT_TYPE part1 = part(f1, canonicalBuffersBM);
            arc = (DCOMPLEX*)myStreams(f1, canonicalBuffersBM, 0);
            w = ar + 4*N2;
//            {//here
//                INT_TYPE j;
//                for ( i = 0; i < N1 ; i++)
//                    for ( j =0 ; j < N1 ; j++)
//                        ar[i*N1+j] = delta(i-j);
//                
//            }
            if ( c1->rt.calcType == electronicStuctureCalculation ){
                INT_TYPE j;

#ifdef BOOTIDENTITY
                    for ( i = 0; i < N1 ; i++)
                        for ( j =0 ; j < N1 ; j++)
                            ar[i*N1+j] = delta(i-j);
#else
                for ( v = 0 ; v < N2 ; v++){
                        if ( f1.tulip[kinetic].spinor == cmpl){
                            cmplFlag = 1;
                            arc[v] = (streams(f1, kinetic,0,space)+space*N2)[v] + I * (streams(f1, kinetic,1,space)+space*N2)[v];
                            printf("%f + I %f\n", creal(arc[v]), cimag(arc[v]));
                        }
                        else
                    {
                            ar[v] = (streams(f1, kinetic,0,space)+space*N2)[v];
                        }
                }
#endif
                }
            
            else if ( c1->rt.calcType == clampProtonElectronCalculation ){
                if ( space < COMPONENT )
                    cblas_dcopy(N2, streams(f1, kinetic,0,space)+space*N2, 1, ar, 1);
                else {
                    cblas_dcopy(N2,streams(f1, kineticMass,0,space), 1, ar, 1);
                    if ( c1->i.oneBody.func.fn == LennardJones){
                    
                   for ( i = 0; i < CanonicalRank(f1, protonRepulsion, 0);i++)
                       cblas_daxpy(N2, 1., streams(f1,protonRepulsion,0,space)+i*N2, 1,ar,1 );
                    
                    printf("LJ %d\n", CanonicalRank(f1, protonRepulsion, 0));
                    }
                    if (0)
                    for ( i = 0 ; i < N1 ; i++) printf("%1.3f,", ar[i*N1+i]);
                    // N12 >= 3
//                    exit(0);
                }
            }
//            if ( c1->i.springFlag ){
//                if ( space == 0 &&  ( c1->rt.runFlag )%2 == 0 ){
//                    cblas_daxpy(N2, 1.0,streams(f1, harmonium,0,0)+0*N2, 1, ar, 1);
//                }
//                if ( space == 1 &&  ( c1->rt.runFlag/2 )%2 == 0 ){
//                    cblas_daxpy(N2, 1.0,streams(f1, harmonium,0,1)+1*N2, 1, ar, 1);
//                }
//                if ( space == 2 &&  ( c1->rt.runFlag/4 )%2 == 0 ){
//                    cblas_daxpy(N2, 1.0,streams(f1, harmonium,0,2)+2*N2, 1, ar, 1);
//                }
//            }
            INT_TYPE complete;
            if ( cmplFlag )
               complete =  tzheev(0, f1, 'V', N1, arc, N1, w);
            else
                complete =   tdsyev(0, f1, 'V', N1, ar, N1, w);
            
            if ( complete ){
                printf("eigensovle failed\n");
                exit(0);
            }
            
            if(0)
            if (  space == 0 || space == COMPONENT ){
                if ( cmplFlag ){
                    for ( i = 0 ; i < N1 ; i++){
                        printf("\n\n%f \n\n", w[i]);

                        for ( ii = 0; ii < N1 ; ii++)
                            printf("\n%d--%f+I%f,",ii+1, creal(arc[i*N1+ii]),cimag(arc[i*N1+ii]));
                    }


                }else {
                    for ( i = 0 ; i < N1 ; i++){
                        printf("\n\n%d %f \n\n",i+1, w[i]);

                        for ( ii = 0; ii < N1 ; ii++)
                            printf("%f,", ar[i*N1+ii]);
                    }

                }
            }

            if ( bootBodies == nada ){

            }   else

            if ( bootBodies == one ){
                
                for ( i = 0  ; i < N1 ; i++)
                    if ( w[i] - w[0] < c1->i.level ){
                        streams(f1,foundationStructure,0,space)[i] = w[i]-w[0];
                        
                        streams(f1,foundationStructure,1,space)[i] = 1;
                    }
                
                if ( cmplFlag ){
                    for ( v = 0 ; v < N2 ; v++){
                        myStreams(f1, bill1+space,0)[v] = creal(arc[v]);
                        myStreams(f1, bill1+space,1)[v] = cimag(arc[v]);
                    }
                } else
                    cblas_dcopy(N2 , ar,1,myStreams(f1, bill1+space,0),1);
            }else if ( bootBodies == two ){
                
                
#ifdef OMP
#pragma omp parallel for private (v,rank,i,ii) schedule(dynamic,1)
#endif
                
                for ( v = 0 ; v < Nx*Nx ; v++){
#ifdef OMP
                    rank = omp_get_thread_num();
#else
                    rank = 0;
#endif
                    
                    i = (v)%Nx;
                    ii = (v/Nx)%Nx;
                    if ( w[i]-w[0] + w[ii] - w[0] < c1->i.level ){

                    
                     //   printf("%d %d %f %f\n", i,ii, c1->i.level, w[i]-w[0] + w[ii] - w[0]);
                        
                        
                        
                        
                        cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,space), 1);
                        f1.tulip[diagonal1VectorA].Current[rank] = 1;
                        cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,space), 1);
                        f1.tulip[diagonal1VectorB].Current[rank] = 1;
                        f1.tulip[diagonalVectorA].Current[rank] = 0;
                        tOuterProductSuOne(f1,space, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonalVectorA, rank);
                        
                        
                        
                        //SA++
                        if ( f1.cat ){
                            f1.tulip[diagonalVectorB].Current[rank] = 0;
                            f1.tulip[diagonalVectorA].Current[rank] = 1;
                            tBuild3IrrOne(rank, f1,space, f1.irrep, diagonalVectorA, rank, diagonalVectorB, rank);
                            if ( cblas_dnrm2(N2, streams(f1, diagonalVectorB,rank,space), 1) > 0.001 ){
                                cblas_dcopy(N2, streams(f1, diagonalVectorB,rank,space), 1, myStreams(f1,bill1+space,0)+v*N2, 1);
                                streams(f1,foundationStructure,0,space)[v] = 0.;
                                for ( em = kinetic1 ; em < kinetic1 + bootBodies ; em++){
                                    tGEMV(rank, f1, space, diagonalVectorA, 0, rank, em, space, 0, diagonalVectorB, 0, rank);
                                    streams(f1,foundationStructure,0,space)[v] += tDOT(rank, f1, space, 1, diagonalVectorA, 0, rank, 1, diagonalVectorB, 0, rank);
                                }
                                if ( streams(f1,foundationStructure,0,space)[v]  < c1->i.level)
                                streams(f1,foundationStructure,1,space)[v] = 1;
                            }
                        } else
                        {
                            cblas_dcopy(N2, streams(f1, diagonalVectorA,rank,space), 1, myStreams(f1,bill1+space,0)+v*N2, 1);
                            streams(f1,foundationStructure,0,space)[v] = w[i]-w[0] + w[ii] - w[0];
                            streams(f1,foundationStructure,1,space)[v] = 1;

                        }
                    
                    }
                    //SA++

                }
                
                
                
                
                
            }else if ( bootBodies == three ){
#ifdef OMP
#pragma omp parallel for private (v,rank,i,ii,iii) schedule(dynamic,1)
#endif
                
                for ( v = 0 ; v < Nx*Nx*Nx ; v++){
#ifdef OMP
                    rank = omp_get_thread_num();
#else
                    rank = 0;
#endif
                    
                    i = (v)%Nx;
                    ii = (v/Nx)%Nx;
                    iii = (v/(Nx*Nx))%Nx;
                    
                    if ( w[i]-w[0] + w[ii] - w[0]+w[iii]-w[0] < c1->i.level )
                    {
                    cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,space), 1);
                    double va = cblas_dnrm2(N1, streams(f1,diagonal1VectorA,rank,space), 1);

                    f1.tulip[diagonal1VectorA].Current[rank] = 1;
                    cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,space), 1);
                    va = cblas_dnrm2(N1, streams(f1,diagonal1VectorB,rank,space), 1);

                    f1.tulip[diagonal1VectorB].Current[rank] = 1;
                    f1.tulip[diagonal2VectorA].Current[rank] = 0;
                    
                    struct name_label d2 = f1.tulip[diagonal2VectorA];
                    
                    
                    tOuterProductSuOne(f1, space,diagonal1VectorA, rank, diagonal1VectorB,rank, diagonal2VectorA, rank);
                    va = cblas_dnrm2(N1*N1, streams(f1,diagonal2VectorA,rank,space), 1);
                    cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,space), 1);
                    
                    f1.tulip[diagonal1VectorC].Current[rank] = 1;
                    f1.tulip[diagonal2VectorA].Current[rank] = 1;

                    f1.tulip[diagonalVectorA].Current[rank] = 0;
                    tOuterProductSuOne(f1,space, diagonal2VectorA, rank, diagonal1VectorC, rank, diagonalVectorA, rank);
                    
                    //SA++
                        if ( f1.cat ){
                            f1.tulip[diagonalVectorB].Current[rank] = 0;
                            f1.tulip[diagonalVectorA].Current[rank] = 1;

                            tBuild3IrrOne(rank, f1,space, f1.irrep, diagonalVectorA, rank, diagonalVectorB, rank);
                        //     printf("\n%d %d %f\n",space, v,cblas_dnrm2(N1*N2, streams(f1, diagonalVectorB,rank,space), 1));
                            if ( cblas_dnrm2(N1*N2, streams(f1, diagonalVectorB,rank,space), 1) > 0.001 ){

                                cblas_dcopy(N1*N2, streams(f1, diagonalVectorB,rank,space), 1, myStreams(f1,bill1+space,0)+v*N1*N2, 1);


                                streams(f1,foundationStructure,0,space)[v] = 0.;
                                for ( em = kinetic1 ; em < kinetic1 + bootBodies ; em++){
                                    tGEMV(rank, f1, space, diagonalVectorA, 0, rank, em, space, 0, diagonalVectorB, 0, rank);
                                    streams(f1,foundationStructure,0,space)[v] += tDOT(rank, f1, space, 1, diagonalVectorA, 0, rank, 1, diagonalVectorB, 0, rank);
                                }
                                if ( streams(f1,foundationStructure,0,space)[v]  < c1->i.level)
                                streams(f1,foundationStructure,1,space)[v] = 1;
                            }
                            //SA++

                        }else
                        {
                            cblas_dcopy(N2*N1, streams(f1, diagonalVectorA,rank,space), 1, myStreams(f1,bill1+space,0)+v*N2*N1, 1);
                            streams(f1,foundationStructure,0,space)[v] = w[i]-w[0] + w[ii] - w[0]+w[iii]-w[0];
                            streams(f1,foundationStructure,1,space)[v] = 1;
                            
                        }
                    }
                }
            }else if ( bootBodies == four ){
                
                
#ifdef OMP
#pragma omp parallel for private (v,rank,i,ii,iii,iv) schedule(dynamic,1)
#endif
                
                for ( v = 0 ; v < Nx*Nx*Nx*Nx ; v++){
#ifdef OMP
                    rank = omp_get_thread_num();
#else
                    rank = 0;
#endif
                    i = (v)%Nx;
                    ii = (v/Nx)%Nx;
                    iii = (v/(Nx*Nx))%Nx;
                    iv = (v/(Nx*Nx*Nx))%Nx;
                    
                    if( w[i]-w[0] + w[ii] - w[0]+w[iii]-w[0] + w[iv]-w[0]< c1->i.level ){
                    cblas_dcopy(N1, ar+i*N1, 1, streams(f1,diagonal1VectorA,rank,space), 1);
                    
                    f1.tulip[diagonal1VectorA].Current[rank] = 1;
                    cblas_dcopy(N1, ar+ii*N1, 1, streams(f1,diagonal1VectorB,rank,space), 1);
                    
                    f1.tulip[diagonal1VectorB].Current[rank] = 1;
                    f1.tulip[diagonal2VectorA].Current[rank] = 0;
                    tOuterProductSuOne(f1, space,diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);

                    cblas_dcopy(N1, ar+iii*N1, 1, streams(f1,diagonal1VectorC,rank,space), 1);
                    //    double u0 = cblas_dnrm2(N1, ar+iii*N1, 1);

                    f1.tulip[diagonal1VectorC].Current[rank] = 1;
                    
                    cblas_dcopy(N1, ar+iv*N1, 1, streams(f1,diagonal1VectorD,rank,space), 1);
                    
                    f1.tulip[diagonal1VectorD].Current[rank] = 1;
                    
                    f1.tulip[diagonal2VectorB].Current[rank] = 0;
                    f1.tulip[diagonalVectorA].Current[rank] = 0;
                    
                    tOuterProductSuOne(f1, space,diagonal1VectorC, rank, diagonal1VectorD, rank, diagonal2VectorB, rank);
                        
                        f1.tulip[diagonal2VectorA].Current[rank] = 1;
                        f1.tulip[diagonal2VectorB].Current[rank] = 1;

                    tOuterProductSuOne(f1,space, diagonal2VectorA, rank, diagonal2VectorB, rank, diagonalVectorA, rank);
                        if (f1.cat ){
                            //SA++
                            f1.tulip[diagonalVectorB].Current[rank] = 0;
                            f1.tulip[diagonalVectorA].Current[rank] = 1;
                            
                            tBuild3IrrOne(rank, f1,space, f1.irrep, diagonalVectorA, rank, diagonalVectorB, rank);
                            // printf("\n%d %d %f\n",space, v,cblas_dnrm2(N1*N2, streams(f1, diagonalVectorB,rank,space), 1));
                            if ( cblas_dnrm2(N2*N2, streams(f1, diagonalVectorB,rank,space), 1) > 0.001 ){
                                
                                cblas_dcopy(N2*N2, streams(f1, diagonalVectorB,rank,space), 1, myStreams(f1,bill1+space,0)+v*N2*N2, 1);
                                
                                
                                streams(f1,foundationStructure,0,space)[v] = 0.;
                                for ( em = kinetic1 ; em < kinetic1 + bootBodies ; em++){
                                    tGEMV(rank, f1, space, diagonalVectorA, 0, rank, em, space, 0, diagonalVectorB, 0, rank);
                                    streams(f1,foundationStructure,0,space)[v] += tDOT(rank, f1, space, 1, diagonalVectorA, 0, rank, 1, diagonalVectorB, 0, rank);
                                }                                if ( streams(f1,foundationStructure,0,space)[v]  < c1->i.level)

                                    streams(f1,foundationStructure,1,space)[v] = 1;
                            }
                            
                            //SA++
                        } else
                        {
                            cblas_dcopy(N2*N2, streams(f1, diagonalVectorA,rank,space), 1, myStreams(f1,bill1+space,0)+v*N2*N2, 1);
                            streams(f1,foundationStructure,0,space)[v] = w[i]-w[0] + w[ii] - w[0]+w[iii]-w[0]+w[iv]-w[0];
                            streams(f1,foundationStructure,1,space)[v] = 1;
                            
                        }
//                       // double u = cblas_dnrm2(N2*N2, streams(f1, diagonalVectorA,rank,space), 1);
//                        cblas_dcopy(N2*N2, streams(f1, diagonalVectorA,rank,space), 1, myStreams(f1,bill1+space,0)+v*N2*N2, 1);
//
//                    streams(f1,foundationStructure,0,space)[v] = w[i]-w[0] + w[ii] - w[0]+w[iii]-w[0] + w[iv]-w[0];
//                    streams(f1,foundationStructure,1,space)[v] = 1;
//                    }
                }
                }
            }
        }
    return 0;
}

INT_TYPE tSortBoot(struct calculation * c1, struct sinc_label f1, enum division eigen){
    assignCores(f1,1);
    enum bodyType bootBodies = f1.rose[0].body;
    INT_TYPE N1,rank;
    INT_TYPE space,N2,i,ii,iii,iv,v;
    double * ar,*w;
    DCOMPLEX * arc ;
    INT_TYPE cmplFlag = 0;
    INT_TYPE tally = 0;
    for ( space = 0 ;space < SPACE ; space++)
        if( f1.rose[space].body != nada){
            
            N1 = vectorLen(f1, space);
            
            tally = 0;
            for ( i = 0  ; i < N1 ; i++)
                if (streams(f1, foundationStructure,1,space)[i] > 0.5 ){
                    if ( i != tally )
                    {
                        cblas_dcopy(N1, myStreams(f1,bill1+space,0)+i*N1, 1, myStreams(f1,bill1+space,0)+(tally)*N1,1);
                        streams(f1, foundationStructure,1,space)[tally] = streams(f1, foundationStructure,1,space)[i];
                        streams(f1, foundationStructure,1,space)[i] = 0.;
                        streams(f1, foundationStructure,0,space)[tally] = streams(f1, foundationStructure,0,space)[i];
                     }
                    tally++;
                }
          //  double tt[N1];
            if ( ! tally )
                return 1;
           // tdgeqr(0, f1, tally, N1, myStreams(f1,bill1+space,0), N1, tt);
            
            for ( i = 0 ; i < tally ; i++)
                cblas_dscal(N1, 1./cblas_dnrm2(N1, myStreams(f1,bill1+space,0)+(i)*N1, 1), myStreams(f1,bill1+space,0)+(i)*N1, 1);
        }
    
    
    return 0;
}
//INT_TYPE tMap (struct calculation * c1 ){
//    struct sinc_label f1 =  c1->i.c.sinc;
//    size_t ms = MAXSTRING;
//    char line0[MAXSTRING];
//    char name[MAXSTRING];
//
//    char *line = line0;
//    INT_TYPE rank = 0;
//    FILE  * list = NULL;
//
//
//    sprintf(name, "%s.body", c1->name);
//    list = fopen(name, "r");
//    if ( list == NULL ){
//        printf("nop");
//        exit(0);
//    }
//    DCOMPLEX one = 1.;
//    int lines= 0,n[6],r[6],space,N1;
//
//    printVector(c1,c1->i.c.sinc,c1->name,c1->name,-1,0, &one);
//
//    getline(&line, &ms, list);
//    while(!feof(list)){
//        enum bodyType outBody = three;
//        if ( outBody == three ){
//            sscanf(line, "(%d,%d)(%d,%d)(%d,%d)", &n[0],&r[0], &n[1],&r[1],&n[2],&r[2]);
//
//            for ( space = 0; space < SPACE ; space++)
//                if ( f1.rose[space].body != nada )
//                {
//                    N1 = f1.rose[space].count1Basis;
//                    cblas_dcopy(N1, streams(f1,f1.user + n[0],0,space)+N1*r[0], 1, streams(f1,diagonal1VectorA,rank,space), 1);
//                    cblas_dcopy(N1, streams(f1,f1.user + n[1],0,space)+N1*r[1], 1, streams(f1,diagonal1VectorB,rank,space), 1);
//                    cblas_dcopy(N1, streams(f1,f1.user + n[2],0,space)+N1*r[2], 1, streams(f1,diagonal1VectorC,rank,space), 1);
//
//                }
//            f1.tulip[diagonal1VectorA].Current[rank] = 1;
//            f1.tulip[diagonal1VectorB].Current[rank] = 1;
//            f1.tulip[diagonal1VectorC].Current[rank] = 1;
//
//
//            f1.tulip[diagonal2VectorA].Current[rank] = 0;
//            f1.tulip[diagonal3VectorA].Current[rank] = 0;
//
//            tOuterProductSu(f1, diagonal1VectorA, rank, diagonal1VectorB, rank, diagonal2VectorA, rank);
//            tOuterProductSu(f1, diagonal2VectorA, rank, diagonal1VectorC, rank, diagonal3VectorA, rank);
//            if ( magnitude(f1, diagonal3VectorA) > 0.01 ){
//                printf("%d %d %d | %d %d %d ", n[0],n[1],n[2],r[0],r[1],r[2]);
//                printf("mag %f\n",magnitude(f1, diagonal3VectorA) );
//                tFilename(c1->name, lines+1, outBody, 0, 0, name);
//                FILE * outVector = fopen(name, "w");
//                outputFormat(f1, outVector, diagonal3VectorA, rank);
//                fclose(outVector);
//                printVector(c1,c1->i.c.sinc,c1->name,c1->name,lines,0, &one);
//            }
//        }
//
//        lines++;
//        getline(&line, &ms, list);
//
//    }
//
//    fclose(list);
//    return 0;
//}



INT_TYPE tSlam (struct sinc_label f1,INT_TYPE allc, enum division vl, double fmax2){
    INT_TYPE tot =0,space,t,n1[SPACE];
    
    
    tot =  tFoundationLevel(f1, nullName, -INFINITY, fmax2, 1, nullName, 0, 0, 0, 0, NULL, 0, 0);
    for ( space = 0 ; space < SPACE ; space++)
        n1[space] = vectorLen(f1, space);
//    if ( allc < tot ){
//        printf("increase foundation %d or decrease levelLevel %f", tot, fmax2);
//        exit(0);
//    }

    INT_TYPE * mmm = malloc(sizeof(INT_TYPE ) * tot * SPACE *2),*mm;
    tot =  tFoundationLevel(f1, nullName, -INFINITY, fmax2, 0, nullName, 0, 0, 0, 0,mmm, 0, 0);

    {
        INT_TYPE n1[SPACE],i;
        {
            INT_TYPE space ;
            for ( space = 0 ; space < SPACE ; space++)
                n1[space]= vectorLen(f1, space);
        }
        struct sortClass * sc = malloc ( sizeof( struct sortClass )*tot); ;
        for ( i = 0; i < tot ; i++){
            for ( space = 0; space < SPACE ; space++){
                sc[i].n1 = n1;
                sc[i].str[space] = streams(f1,foundationStructure,0,space);
            }
            sc[i].i = i;
            sc[i].mmm = mmm;
        }
        qsort(sc, tot, sizeof(struct sortClass), &sortComp);
        
        
        
        for ( t = 0; t < imin(tot,allc) ; t++){
            mm = mmm + 2*SPACE * sc[t].i ;
            printf("%d %d %d -- %d %d %d :: %f\n", mm[0],mm[2],mm[4], mm[1],mm[3],mm[5],vale(sc+t));
            for ( space = 0; space < SPACE ; space++)
                if ( f1.rose[space].body != nada){
                    cblas_dcopy(n1[space], myStreams(f1, bill1+space, 0)+(mm[2*space])*n1[space]+(mm[2*space+1])*n1[space]*n1[space],1,streams(f1,vl+t,0,space),1);
                }
            f1.tulip[vl+t].Current[0] = 1;
            f1.tulip[vl+t].value.value = vale(sc+t);
            f1.tulip[vl+t].value.stage = t;
//            tBuild3Irr(0, f1, f1.irrep, vl+t, 0, copyVector, 0);
//            tEqua(f1, vl+t, 0, copyVector, 0);
            //testSA(f1, vl+t);
        }
        free(sc);
    }
    free(mmm);
    return imin(tot,allc);
}


INT_TYPE tBootManyConstruction (struct calculation * c1, struct sinc_label f1, enum division eigen){
    INT_TYPE sp,cmpl,space;
    //THRESING FLOOR
    INT_TYPE r,i,im;
    INT_TYPE n2[SPACE];
    length(f1, eigen,n2);
    INT_TYPE n1[SPACE];
    for ( space = 0; space< SPACE ; space++)
        n1[space] = sqrt(1.*n2[space]);

    DCOMPLEX * hmat = (DCOMPLEX*)myStreams(f1, matrixHbuild,0), sum ,minus = -1.;
    double * w = (double*)(hmat + n2[0]);
    cmpl = real;
    printf("\nlevel %f\n",c1->i.level);

    for ( i = 0 ;i < cmpl ; i++)
       balance(f1, eigen,i);///NEEDS WORK FOR COMPLEX H's
            for ( space =  0; space < SPACE ; space++ )
                if ( f1.rose[space].body != nada){
                    r = 0;
               //     for ( r = 0; r <part(f1,eigen)/*ASSUME FILLED*/ ; r++)
                    {
                        for ( i = 0; i < n2[space]; i++){
                            hmat[i] =0.;
                            for ( sp = 0 ; sp < cmpl;sp++)
                                if ( sp == 0 )
                                    hmat[i] += streams(f1, eigen,sp,space)[i+r*n2[space]] ;
                                else
                                    hmat[i] += I*streams(f1, eigen,sp,space)[i+r*n2[space]] ;

                        }
                        tzheev (0,f1,'V',n1[space],hmat,n1[space],streams(f1,foundationStructure,0,space)+r*n1[space]);
                        printf("s%d -- %f %f\n", space,streams(f1,foundationStructure,0,space)[0],streams(f1,foundationStructure,0,space)[1] );
                        for ( i = 0; i < n2[space]; i++){
                            myStreams(f1, bill1+space,0)[i+r*n2[space]] = creal(hmat[i]);
                            if ( spins (f1, eigen) > 1 )
                                myStreams(f1, bill1+space,1)[i+r*n2[space]] = cimag(hmat[i]);
                        }
                        for (i = 0; i < n1[space];i++)
                            if (((streams(f1,foundationStructure,0,space)+r*n1[space])[i]) < c1->i.level  ){
                                (streams(f1,foundationStructure,1,space)+r*n1[space])[i] = 1;
                            }

                    }
                }
    return 0;
}



INT_TYPE tFoundationLevel( struct sinc_label  f1, enum division A , double lvlm, double lvlx,INT_TYPE ops,enum division build,INT_TYPE xB, double lvl1, double lvl2, double lvl3,INT_TYPE *mmm, INT_TYPE irrep,double seekPower){
    //GRID
    /// GRID
    ///// GRID
    ////////GRID
    ADDRESS_TYPE mx=1,i,j,k,r1,r2,r3,ii,jj,kk,xx[SPACE];
    INT_TYPE vaMax=0, classicalBasisSize,*mm,space;
    INT_TYPE nG = 1;//tSize(bd);
    INT_TYPE n1[SPACE];
    INT_TYPE r[SPACE];
    INT_TYPE GGG,g[SPACE],iii,flag;
    for ( space = 0 ; space < SPACE ; space++)
        n1[space] = vectorLen(f1, space);
    double value;
    classicalBasisSize = 0;
    //for( lvl = 0; lvl < imin(n1[0],imin(n1[1],n1[2])) ; lvl++)
        
    for ( space = 0; space < SPACE ;space++){
            if ( n1[space] ){
                xx[space] = 1;
                for ( i = 0; i < n1[space] ; i++)
                    if ( streams(f1,foundationStructure,1,space)[i] > 0.5 ){
                        xx[space] = i+1;
                    }
            }else
                xx[space] = 1;
        //printf("xx %d %d\n", xx[space],n1[space]);
    }
    
        for ( space = 0 ; space < SPACE ; space++)
                mx *= nG*xx[space];
    
        for ( iii = 0; iii < mx ; iii++){
            
            GGG = 1;
            for ( space = 0; space < SPACE ; space++)
                {
                    g[space] = (iii/GGG) % nG;
                    GGG *= nG;
                    r[space] = (iii/GGG) % xx[space];
                    GGG *= xx[space];
                }

            
            flag = 1;
            value = 1;
            for ( space = 0; space < SPACE ; space++)
                if ( n1[space] != 0 ) {
                if (  streams(f1,foundationStructure,1,space)[r[space]+g[space]*n1[space]] > 0.5    ){
                    value *= (streams(f1,foundationStructure,0,space)[r[space]+g[space]*n1[space]]);

                } else {
                    flag = 0;
                }
                }
            
            if ( flag )
            if ( vaMax < value )
                vaMax = value;
            
            
            if (flag &&  lvlm <= value && value < lvlx){
                {
                    
                    
                    if ( 1 ){
                        //if ( SPACE <= 3 && ! tIR ( bd,g[0],g[1],g[2],irrep) ){
                        //    continue;
                       // }
                        if ( ops  == 0 ){
#if VERBOSE
                            printf("**%d-->%d : %d %d %d :: %d %d %d :: %f\n",irrep,classicalBasisSize,r[0],r[1],r[2],g[0],g[1],g[2],value);

#endif
                            mm = mmm+classicalBasisSize*SPACE*2;
                            
                            for ( space = 0 ; space < SPACE ; space++){
                                mm[2*space] = r[space];
                                mm[2*space+1] = g[space];
                            }
                        }
                     //   printf("%d-->%d : %d %d %d :: %d %d %d :: %f\n",irrep,classicalBasisSize,r[0],r[1],r[2],g[0],g[1],g[2],value);

                        classicalBasisSize++;
                    }
                }
        }
        }
    if ( ops == 2 )
        return vaMax;
    else
        return classicalBasisSize;
}



INT_TYPE tFilter(struct sinc_label f1, INT_TYPE Ve, INT_TYPE irrep, enum division usr){
    INT_TYPE i,ii,space,cmpl=0,rank;
    assignCores(f1, 1);
    DCOMPLEX norm;
    rank = 0;

    if ( irrep && bodies(f1, usr ) > one ){
        
        printf("Symmetry Adaption -> %d-irrep\n", irrep );
        for( ii= 0 ; ii < Ve ; ii++){
            for ( cmpl = 0 ; cmpl < spins(f1, usr+ii) ; cmpl++){
                
                f1.tulip[copyTwoVector].Current[rank] = 0;
                {
                    if ( f1.cat )
                        tBuild3Irr(rank, f1, irrep, usr+ii, cmpl, copyTwoVector, rank);
                    else
                        tBuildIrr(rank, f1, irrep, usr+ii, cmpl, copyTwoVector, rank);
                }
                
                tCycleDecompostionGridOneMP(-2, f1, copyTwoVector, rank, NULL,usr+ii , cmpl, f1.rt->TARGET, part(f1,usr+ii), f1.rt->powDecompose);
            }
            pMatrixElements( f1, usr+ii, nullName, usr+ii, NULL, &norm);
            tScale(f1, usr+ii, 1./sqrt(cabs(norm)));
        }
//        {
//            rank = 0;
//
//            for ( cmpl = 0 ; cmpl < spins(f1, usr+ii) ; cmpl++){
//                f1.tulip[copyTwoVector].Current[rank] = 0;
//                for ( i = 0 ; i < CanonicalRank(f1, usr+ii, cmpl);i++)
//                {
//                    f1.tulip[copyVector].Current[rank] = 0;
//                    for ( space = 0; space < SPACE ; space++)
//                        xsAdd(1., space, f1, copyVector, rank, f1, usr+ii, i, cmpl);
//                    f1.tulip[copyVector].Current[rank] = 1;
//                    //SA++
//                    if ( f1.cat )
//                        tBuild3Irr(rank, f1, irrep, copyVector, rank, copyTwoVector, rank);
//                    else
//                        tBuildIrr(rank, f1, irrep, copyVector, rank, copyTwoVector, rank);
//                }
//                tCycleDecompostionGridOneMP(-2, f1, copyTwoVector, rank, NULL,usr+ii , cmpl, f1.rt->TARGET, part(f1,usr+ii), f1.rt->powDecompose);
//            }
//            pMatrixElements( f1, usr+ii, nullName, usr+ii, NULL, &norm);
//            tScale(f1, usr+ii, 1./sqrt(cabs(norm)));
//        }
    }
    
#ifdef OMP
#pragma omp parallel for private (ii,rank) schedule(dynamic,1)
#endif
    for ( ii = 0; ii < Ve ; ii++)
    {
#ifdef OMP
        rank = omp_get_thread_num();
#else
        rank = 0;
#endif
        f1.tulip[usr+ii].value.symmetry = tClassify(rank, f1, usr+ii);
    }

    return 0;
}

INT_TYPE tSelect(struct sinc_label  f1, INT_TYPE Ve, INT_TYPE type, enum division usr, enum division usa, INT_TYPE testFlag){
    INT_TYPE sp,info,rank=0,maxEV = f1.maxEV,n,m;
    INT_TYPE stride = maxEV;
    double value;
    
	if ( ! CanonicalRank(f1, usa,0))
        return 0;


    if ( type ) {
        tClear(f1, usr+Ve);
        for ( sp = 0; sp < spins(f1, usr+Ve);sp++){
            f1.tulip[copyVector].Current[rank] = 0;
            tBuildIrr(rank, f1, type, usa, sp, copyVector, rank);
            tCycleDecompostionListOneMP(rank, f1, copyVector, rank, NULL,usr+Ve, sp, f1.rt->vCANON, part(f1,usr+Ve), -1);
        }
        
        value = magnitude(f1, usr+Ve);
        
        if ( isnan(1./value) || isinf(1./value) )
           return 0;
        tScale(f1, usr+Ve, 1./value);

    } else {
	
        value = magnitude(f1, usa);
        
        
       if ( value < 1e-6 || isnan(value) || isinf(value) )
            return 0;
        
        if ( testFlag >= 0 ){
            for ( sp = 0; sp < 1;sp++){

                tEqua(f1, usr+Ve,sp, usa,sp);
            }
            tScale(f1, usr+Ve, 1./value);
            
        }
    }
    if ( testFlag != 1 )
        return 1;

    
    DCOMPLEX HH, * S = (DCOMPLEX*)(myStreams(f1, matrixSbuild, 0));
    myZero(f1, matrixSbuild,0);
    double * ov = myStreams(f1, twoBodyRitz, 0);
    assignCores(f1, 1);
#ifdef OMP
#pragma omp parallel for private (m,n,rank,HH) schedule(dynamic,1)
#endif
    for ( n = 0; n < Ve+1 ; n++)
    {
        
#ifdef OMP
        rank = omp_get_thread_num();
#else
        rank = 0;
#endif
        
        for ( m = 0 ; m <= n   ; m++)    {
            matrixElements(rank, f1, usr+n, nullName, usr+m, &HH, S+n*stride+m);
            S[m*stride+n] = conj(S[n*stride+m]);
        }
        
    }
    assignCores(f1, 0);
    INT_TYPE complete;
    complete = tzheev(0, f1, 'N', Ve+1,S, stride, ov);
    if ( complete ){
        printf("eigensovle failed\n");
        exit(0);
    }

    if ( testFlag ){

        if ( ov[Ve]/ov[0] < f1.rt->TOL && ov[Ve]/ov[0] >0 && ov[0] > 0 ){
                printf("\t**%f**\t%d\n",  ov[Ve]/ov[0],Ve+1 );
                fflush(stdout);

                return 1;
        } else {
        }
    }
    
    return 0;
}


//INT_TYPE tCollect (struct sinc_label  f1, INT_TYPE irrep,enum division usz, INT_TYPE target,double seekPower){
//    INT_TYPE  Ve = 0,Vex;
//    f1.tulip[diagonalVectorA].header = Cube;
//    f1.tulip[diagonalVectorB].header = Cube;
//    struct name_label bd = f1.tulip[build];
//    struct name_label eg = f1.tulip[eigen];
//
//    INT_TYPE space, nG = part(f1,eigen);
//    INT_TYPE flag = 1,ct = 0;
//    double min0= 0,max0= 1+2*tFoundationLevel(f1, build,0,0,2,0,target,1e9,1e9,1e9,NULL,irrep,seekPower) ;
//    double min = min0, max =  max0,vx= 100.,va = 1.;
//    printf("\n\n\t| Seek %d Body vectors in %d components \t|\n", bodies(f1, eigenVectors), SPACE);
//    printf("\t| Greater than Target \t: %6d\t\t|\n",target);
//    
//    ct = tFoundationLevel(f1, build,0,max0,1,0,target,1e9,1e9,1e9,NULL,irrep,seekPower);
//    while(1){
//        if ( flag ){
//            if ( ct < target*va )
//                min = 0.5*(max+min);
//            else
//                max = 0.5*(max+min);
//        }
//        ct = tFoundationLevel(f1, build,0,0.5*(min+max),1,usz,target,1e9,1e9,1e9,NULL,irrep,seekPower);
//        printf("\t| Current \t: %6d \t %f\t\t|\n",ct,0.5*(min+max) );
//        if ( max-min < 1e-9  ){
//            printf("conv");
//            exit(0);
//        }
//        if ( ct >= target && ct <= target*va*vx)
//        {
//            
//            Ve = 0;
//            INT_TYPE *mmm = malloc(sizeof(INT_TYPE ) *2*SPACE*ct),i,*mm;
//            ct = tFoundationLevel(f1, build,0,0.5*(min+max),0,usz,target,1e9,1e9,1e9,mmm,irrep,seekPower);
//            INT_TYPE n1[SPACE];
//            {
//                INT_TYPE space ;
//                for ( space = 0 ; space < SPACE ; space++)
//                    n1[space]= vectorLen(f1, space);
//            }
//            struct sortClass * sc = malloc ( sizeof( struct sortClass )*ct); ;
//            for ( i = 0; i < ct ; i++){
//                for ( space = 0; space < SPACE ; space++){
//                    sc[i].n1[space] = n1[space];
//                    sc[i].str[space] = streams(f1,foundationStructure,0,space);
//                }
//                sc[i].i = i;
//                sc[i].nG = nG;
//                sc[i].mmm = mmm;
//            }
//            
//            qsort(sc, ct, sizeof(struct sortClass), &sortComp);
//            
//            
//            for ( i = 0; i < ct ; i++){
//                mm = mmm+sc[i].i*2*SPACE;
//                tClear(f1,diagonalVectorA);
//                f1.tulip[diagonalVectorA].Current[0] = 1;;
//                
//                for ( space = 0; space < SPACE ; space++)
//                    if ( f1.rose[space].body != nada)
//                    cblas_dcopy(n1[space], myStreams(f1, bill1+space, 0)+(mm[2*space])*n1[space]+(mm[2*space+1])*n1[space]*n1[space],1,streams(f1,diagonalVectorA,0,space),1);
//                
//                //Vex =  tSASplit(f1, irrep, Ve,target, usz, diagonalVectorA);
//                if ( Vex - Ve ){
//                    printf(" |%d| %d %d %d -> %f\n",i, mm[0],mm[2],mm[4],vale(sc+i));
//                }
//                Ve = Vex;
//                if ( Ve > target -2)
//                    break;
//            }
//            free(mmm);
//            free ( sc);
//
////            if ( Ve == target )
//               break;
////            else
////            {
////                printf("instead try floorFlag %d\n", Ve);
////                va *= 1.5;
////                max = max0;
////                printf("\t| increasing buffer region to %1.3f\t|\n", va);
////
////                if ( va > 2 )
////                    exit(2);
////            }
//        }
//        
//    }
////    if ( target != Ve ){
////        printf("ack no !\n");
////        exit(0);
////    }
//    
//    return Ve;
//}


//INT_TYPE tSASplit ( struct sinc_label  f1, INT_TYPE irrep , INT_TYPE Ve , INT_TYPE target,enum division usz, enum division vector){
//    INT_TYPE map[24],nDeg=0,ii;
//    enum bodyType body = f1.rose[0].body;
//
//    if ( body == one || irrep == 0 ){
//        nDeg = 1;
//        map[1] = 0;
//    }else
//
//        if ( body == two ){
//            if ( irrep == 1 ){
//                map[1] = 1;
//                nDeg = 1;
//            }else
//                if ( irrep == 2 ){
//                    nDeg = 1;
//                    map[1] = 2;
//                }
//        }else
//            if ( body== three ){
//                if ( irrep == 1 ){
//                    map[1] = 1;
//                    nDeg = 1;
//                }else
//                    if ( irrep == 2 ){
//                        map[1] = 2;
//                        nDeg = 1;
//                    } else if ( irrep == 3 ){
//                        nDeg = 4;
//                        map[1] = 3;
//                        map[2] = 4;
//                        map[3] = 5;
//                        map[4] = 6;
//                    }
//            }
//            else if ( body == four ){
//                //
//                if ( irrep == 1 ){
//                    nDeg = 1;
//                    map [1] = 1;
//                }else if ( irrep == 2 ){
//                    map[1] = 2;
//                    nDeg = 1;
//                } else if ( irrep == 3 ){
//                    map[1] = 3;
//                    map[2] = 4;
//                    map[3] = 5;
//                    map[4] = 6;
//                    nDeg = 4;
//                } else if ( irrep == 4){
//                    map[1] = 7;
//                    map[2] = 8;
//                    map[3] = 9;
//                    map[4] = 10;
//                    map[5] = 11;
//                    map[6] = 12;
//                    map[7] = 13;
//                    map[8] = 14;
//                    map[9] = 15;
//                    nDeg = 9;
//                }else if ( irrep == 5 ){
//                    map[1] = 16;
//                    map[2] = 17;
//                    map[3] = 18;
//                    map[4] = 19;
//                    map[5] = 20;
//                    map[6] = 21;
//                    map[7] = 22;
//                    map[8] = 23;
//                    map[9] = 24;
//                    nDeg = 9;
//                }
//            }
//
//
//    for ( ii = 1 ; ii <= nDeg ; ii++)
//
//        if ( tSelect(f1, Ve, map[ii]+tPerms(body), usz, vector, 1)){
//            Ve++;
//            if ( Ve > target-2 )
//                return Ve;
//
//        }
//    return Ve;
//}


INT_TYPE tSquareVectors(struct sinc_label f1, INT_TYPE EV2, enum division usz,enum division usr ){
    INT_TYPE i,j;

    for ( i = 0; i < EV2 ; i++)
        for ( j = 0 ; j < EV2 ; j++)
        {
            tClear(f1,usz+i*EV2+j );
            tOuterProductSu(f1, usz+i, 0, usz+j, 0, usr+i*EV2+j, 0);
            tOuterProductSu(f1, usz+i, 1, usz+j, 1, usr+i*EV2+j, 1);
        }    
    return EV2*EV2;
}


INT_TYPE tGreatDivideIteration (INT_TYPE translateFlag ,double sumPart, double realPart, struct sinc_label f1, enum division A , INT_TYPE I1, INT_TYPE I2, enum division usz, INT_TYPE foundation,INT_TYPE nMult, INT_TYPE shift){
    INT_TYPE expon,info;
    INT_TYPE rank ;
    //time_t start_t, lapse_t;
    DCOMPLEX temp,temp2, sum = 0,vhhv,vhv,vhhhv;
    //time(&start_t);
    INT_TYPE iii = 0;
    
    if ( translateFlag ){
        printf (" %1.6f | PSI > + %1.6f H | PSI > \n",sumPart,realPart);
    }else {
        printf (" %1.6f H | PSI >\n", realPart);
    }
    
    for( expon = 1 ; foundation*expon < nMult  ; expon++){
        
//#ifdef OMP
//#pragma omp parallel for private (iii,rank,vhhv,vhv) schedule(dynamic,1)
//#endif
        for ( iii = 0; iii < foundation ; iii++)
        {
            
//#ifdef OMP
//            rank = omp_get_thread_num();
//#else
            rank = 0;
//#endif
            
            
            {
                tEquals(f1,usz+iii+expon*foundation,usz+(expon-1)*foundation+iii );//initialize
//                matrixElements(rank, f1, usz+iii+expon*foundation, nullName, usz+iii+expon*foundation, NULL, &vhhv);
//                tScale(f1, usz+iii+expon*foundation, 1./sqrt(cabs(vhhv)));

                tHXpX(rank, f1, A, translateFlag,sumPart, translateFlag*realPart+ (!translateFlag)*1., 0.0, usz+iii+expon*foundation, f1.rt->TARGET , part(f1,usz+(expon)*foundation+iii),1 == 1);
                
                pMatrixElements( f1, usz+iii+expon*foundation, nullName, usz+iii+expon*foundation, NULL, &vhhv);
                pMatrixElements( f1, usz+iii+expon*foundation, nullName, usz+iii+(expon-1)*foundation, NULL, &vhv);

                f1.tulip[usz+iii+expon*foundation].value.value =creal(vhv)-sumPart ;
                
                
                printf("%d\t uncertainity:  \t %f\n", iii+1,  creal(vhhv) - (creal(vhv))*(creal(vhv)));
                if ( cabs(vhhv) > 0. )
                    tScale(f1, usz+iii+expon*foundation, 1./sqrt(cabs(vhhv)));
                else {
                    printf("oops norms \n");
                    return 1;
                }
                pMatrixElements( f1, usz+iii+expon*foundation, nullName, usz+iii+expon*foundation, NULL, &vhhv);
               if ( fabs(creal(vhhv)-1.0)> 1e-6)
                   printf("check %f\n", creal(vhhv));
            }
            
        }
    }
    
    
    
    
    
    INT_TYPE sum2 = 0;
    expon = 1;
    for ( iii = 0; iii < foundation ; iii++)
        sum2 += CanonicalRank(f1, usz+(expon)*foundation+iii, 0)+CanonicalRank(f1, usz+(expon)*foundation+iii, 1);
    
    printf("-------\n\tAve-Canonical Rank \t| %2.1f \n", sum2*1. / foundation );    
    
    return 0;
}

INT_TYPE tLesserDivideIteration ( struct sinc_label  f1, enum division A , INT_TYPE I1, INT_TYPE I2, enum division usz, INT_TYPE foundation,INT_TYPE nMult, INT_TYPE shift);

INT_TYPE tEdges(struct sinc_label f1, enum division vector){
    INT_TYPE info,spatial;
    DCOMPLEX ov,me;
    enum bodyType bootBodies = f1.rose[0].body;
    double sum[4],totalSum = 0.;

    if ( 1 ){
        //EDGES ALT
        enum block b,bx;
        INT_TYPE iii,jjj=1,dim,irrep;
        sum[0] = 0.;
        sum[1] = 0.;
        sum[2] = 0.;
        sum[3] = 0.;
        //for ( irrep = 0 ;irrep <= 5 ; irrep++)
        {
              //  if ((! c1->i.irrep || f1.tulip[vector].value.symmetry  == irrep) && irrep == c1->i.irrep)
                {
                    bx = tv1;
                    if ( bootBodies == two )
                        bx = tv2;
                    else if ( bootBodies == three )
                        bx = tv3;
                    else if ( bootBodies == four )
                        bx = tv4;
                    for ( b = tv1 ; b <= bx; b++){
                        
                        for ( spatial = 0 ; spatial < 2 ; spatial++){
                            printf("electron %d:%d\t",b,spatial);

                            
#ifdef SPINOR
                            for ( dim = ELEC; dim < SPACE ; dim++)

#else
                            for ( dim = 0; dim < COMPONENT ; dim++)

#endif
                            
                            {
                                tClear(f1,edgeElectronMatrix );
                                if ( spatial == 0 )
                                    tEnd(f1, edgeElectronMatrix, 0, dim);
                                if ( spatial == 1 )
                                    tAlt(f1, edgeElectronMatrix, 0, dim);
                                me = 0.;
                             //   enum division u = edgeElectronMatrix+b;
                              //  struct name_label uu = f1.tulip[u];
                                pMatrixElements( f1, vector, edgeElectronMatrix+b, vector, &me, &ov);
                                printf("%1.8f ", (creal(me/ov)));
                                sum[spatial] += (creal(me/ov));
                                totalSum += sum[spatial];
                            }
                            printf("\n");
                        }
//
                        if ( f1.rt->calcType == clampProtonElectronCalculation)
                            for ( spatial = 0 ; spatial < 2 ; spatial++){
                                printf("proton %d:%d\t",b,spatial);

                                for ( dim = COMPONENT; dim < 2*COMPONENT ; dim++)
                                    if ( f1.rose[dim].body != nada ){

                                    tClear(f1,edgeProtonMatrix );
                                    if ( spatial == 0 )
                                        tEnd(f1, edgeProtonMatrix, 0, dim);
                                    if ( spatial == 1 )
                                        tAlt(f1, edgeProtonMatrix, 0, dim);
                                    me = 0.;
                                    pMatrixElements( f1, vector, edgeProtonMatrix+b, vector, &me, &ov);
                                    printf("%1.8f ", creal(me/ov));
                                    sum[2+spatial] += (creal(me/ov));
                                    totalSum += sum[spatial];

                                }
                               printf("\n");
                            }
                        
                    }
                    printf("\n\n");
                }
        }
    }
    
    if ( totalSum < f1.rt->CAP ){
        return 0;
    }else {
        return cblas_idamax(4, sum, 1)+1;
    }
    return 0;
}

//INT_TYPE tConvergeTest (struct calculation * c1, enum division input){
//
//
//    {
//        struct sinc_label * f1 = &c1->i.c;
//        double sum2,sum;
//        enum division Mat, leftP,ml[1000];
//        INT_TYPE Lane = MaxCore,iii,i,i2,j,j2,rk,s1=0,cat=0,rank,iiii= 0;
//        double Sa [c1->i.bRank  *c1->i.bRank];
//        double Saa [c1->i.bRank  ];
//
//        leftP = input;
//
//
//        do {
//
//            if ( linear == name(f1,leftP) ){
//                Mat = rivers(Lane++, f1, linear, cat );
//                f1.tulip[Mat].block = f1.tulip[leftP].;
//            }else {
//                Mat = leftP;
//            }
//
//            if ( Lane -MaxCore >= MaxCore )
//            {
//                printf("neeed more \n");
//                exit(0);
//            }
//
//            if ( CanonicalRank(f1, Mat, 0)+CanonicalRank(f1, Mat, 1))
//                ml[iiii++] = Mat;
//
//            if (name(f1, leftP) == linear ){
//                cat++;
//
//                if ( cat > f1->Na ){
//                    leftP = f1.tulip[leftP].linkNext;
//                    cat = 0;
//                }
//
//            } else {
//                leftP = f1.tulip[leftP].linkNext;
//            }
//        } while ( leftP != nullName);
//
//
//
//        for ( iii = 0; iii < c1->i.nStates ; iii++){
//            sum = 0.;
//            sum2 = 0.;
//            rk = CanonicalRank(f1, eigenVectors+iii, 0);
//
//            for ( j = 0; j < iiii ; j++)
//            {
//
//#pragma omp parallel for private (i,i2,j2,rank) schedule(dynamic,1)
//                for ( i = 0; i < rk;i++)
//                {
//#if ARES
//
//                    rank = omp_get_thread_num();
//#else
//                    rank = 0;
//#endif
//
//                    //                                if ( Nb){
//                    //                                    tEqua(f1,squareVector, rank, ocean(rank,f1,eigenVectors+iii,i,0),0);
//                    //                                    tpProduct(rank, f1, ml[j], 0, squareVector, rank);
//                    //                                }
//                    //                                else
//                    tMultiplyMP(rank,f1,1., -1,squareVector, rank, 'N', ml[j] , 0,'N', ocean(rank,f1,eigenVectors+iii,i,0), s1);
//
//
//                    Saa[i] = tMultiplyMP(rank,f1,1., -1,nullVector, 0, CDT, squareVector, rank,'N', ocean(rank,f1,eigenVectors+iii,i,0),0);
//
//                    for ( j2 = 0 ; j2 < iiii ; j2++)
//                        for ( i2 = 0 ; i2 <= i   ; i2++){
//                            //                                        if ( Nb){
//                            //                                            tEqua(f1, totalVector, rank, ocean(rank,f1,eigenVectors+iii,i2,0),0);
//                            //                                            tpProduct(rank, f1, ml[j2] , 0, totalVector, rank);
//                            //                                        }
//                            //                                        else
//                            tMultiplyMP(rank,f1,1., -1,totalVector, rank, 'N',ml[j2] , 0,'N', ocean(rank,f1,eigenVectors+iii,i2,0), s1);
//
//                            Sa[i*rk+i2] = tMultiplyMP(rank,f1,1., -1,nullVector, 0, CDT, squareVector,rank,'N', totalVector, rank);
//
//                        }
//                }
//
//                for ( i = 0; i < rk ; i++){
//                    for ( i2 = 0; i2 < rk ; i2++){
//                        sum2 += Sa[i*rk+i2];
//                    }
//                    sum += Saa[i];
//                }
//            }
//
//            printf("Cv%lld \t %f\n", iii+1, sum2-sqr(sum));
//
//        }
//    }
//}



INT_TYPE tEigenCycle (INT_TYPE typer,INT_TYPE minusFlag, struct sinc_label  f1, enum division A ,char permutation,  INT_TYPE Ne, enum division usz, INT_TYPE quantumBasisSize ,INT_TYPE iterations, INT_TYPE foundation, INT_TYPE irrep,INT_TYPE flag,  enum division outputSpace, enum division outputValues){
    INT_TYPE in,gvOut,prevBuild;
    time_t start_t, lapse_t;
    myZero(f1, matrixHbuild, 0);
    myZero(f1, matrixSbuild, 0);
    time(&start_t);
    enum division Mat;
    INT_TYPE cmpl,cmpl2,cmpl3,cat,iii = 0,maxEV = f1.maxEV,rank;
    INT_TYPE stage, maxStage=0,minStage,stride = maxEV;
    double * ritz = myStreams(f1, outputValues, 0);
    enum division el ;
    DCOMPLEX *T  =  (DCOMPLEX *) myStreams(f1, matrixHbuild,0/*CORE RANK*/);
    DCOMPLEX *S  =  (DCOMPLEX *) myStreams(f1, matrixSbuild,0/*CORE RANK*/);
    DCOMPLEX *t  =  T+stride*stride;
    DCOMPLEX *s  =  S+stride*stride;
    
    if ( part(f1, matrixHbuild) * sizeof(double) < 2*stride*stride*sizeof(DCOMPLEX) ){
        printf("ack|n");
        exit(0);
    }
    if ( part(f1, matrixSbuild) * sizeof(double) < 2*stride*stride*sizeof(DCOMPLEX) ){
        printf("ack|n");
        exit(0);
    }

    if ( part(f1, outputValues) < stride ){
        printf("ack|n");
        exit(0);
    }

    INT_TYPE qs,aa[stride];
    INT_TYPE powerMat;

    INT_TYPE info,n,m,s1,g,r,rr,rx,gx,a,a2;
    enum division leftP ;
    prevBuild = 0;
    minStage = f1.tulip[usz].value.stage;
    for ( n = prevBuild; n < quantumBasisSize ; n++)
    {
        if ( magnitude(f1, usz+n) == 0. )
        {
            printf("normaly, I would punch you!\n");
            exit(0);
        }
        maxStage = imax(maxStage ,f1.tulip[usz+n].value.stage);
        minStage = imin(minStage ,f1.tulip[usz+n].value.stage);

        tScale(f1, usz+n, 1./magnitude(f1, usz+n));
    }
    
    assignCores(f1, 1);

        leftP = A;
        
        do {
            Mat = leftP;
        
#if VERBOSE
            if ( Rank(f1,name(f1, Mat))){
                INT_TYPE space;
                for ( space = 0 ;space < SPACE ; space++)
                    printf("%2d", f1.tulip[name(f1,Mat)].space[space].body);
                if (CanonicalRank(f1, name(f1,Mat), 1))
                    printf("\t::\t %d \t| (%d+i%d)\t",name(f1,Mat), CanonicalRank(f1, name(f1,Mat), 0),CanonicalRank(f1, name(f1,Mat), 1));
                else
                    printf("\t::\t %d \t| (%d)\t",name(f1,Mat), CanonicalRank(f1, name(f1,Mat), 0));
                for ( space = 0 ;space < SPACE ; space++)
                    printf("%3d", f1.tulip[Mat].space[space].block);
               // printf("\t%f\n", traceOne(f1,name(f1, Mat),0));
                printf("\n");
            }
            fflush(stdout);
#endif
#ifdef OMP
#pragma omp parallel for private (in,m,rank,n) schedule(dynamic,1)
#endif
            for ( in = 0 ;in < quantumBasisSize * quantumBasisSize; in++)
            {
                
#ifdef OMP
                rank = omp_get_thread_num();
#else
                rank = 0;
#endif
                m = in % quantumBasisSize;
                n = (in/quantumBasisSize) % quantumBasisSize;
                if ( m <= n ){
                    matrixElements(rank, f1, usz+n, Mat, usz+m, T+n*stride+m, S+n*stride+m);
                    T[m*stride+n]= conj(T[n*stride+m] ) ;
                    S[m*stride+n] = conj(S[n*stride+m ]) ;
                }
            }
            
            leftP = f1.tulip[leftP].linkNext;
        } while ( leftP != nullName);
    
    if (minusFlag){
        DCOMPLEX one = -1.;
        cblas_zscal(quantumBasisSize*quantumBasisSize,&one,T,1);
    }
    
    if(0){
        INT_TYPE ii;
        for ( ii = 0 ; ii < quantumBasisSize ; ii++){
            printf("%d %f %f d \n", ii+1, creal(T[ii*stride+ii]), creal(S[ii*stride+ii]));
            printf("%d %f %f x \n", ii+1, creal(T[ii*stride]), creal(S[ii*stride]));
        }
    }
    
    if ( 0){
        printf("\n\nH={");
        for ( m = 0; m < quantumBasisSize ; m++){
            if ( m )
                printf("},{");
            else
                printf("{");
            for ( n =0 ; n < quantumBasisSize ; n++){
                if ( n )
                    printf(",");
                printf("%f+I*%f", creal(T[n*stride+m]),cimag(T[n*stride+m]));
            
            }
        }
        printf("}};");

        
    }
    if ( 0){
        printf("S={");
        for ( m = 0; m < quantumBasisSize ; m++){
            if ( m )
                printf("},{");
            else
                printf("{");
            for ( n =0 ; n < quantumBasisSize ; n++){
                if ( n )
                    printf(",");
                printf("%f+I*%f", creal(S[n*stride+m]),cimag(S[n*stride+m]));
                
            }
        }
        printf("}};\n\n");
        
        
    }


    //    if(0)
//    for ( n = 0; n < quantumBasisSize-foundation; n++)
//        for ( m = 0 ; m <= n   ; m++){
//            if ( part(f1,usz+m) > minRank && part(f1, usz+n) > minRank && n % foundation == m % foundation){
//               // printf("%lld %lld %f %f\n", n/foundation,m/foundation,T[n*stride+m],T[n*stride+m+maxEV*maxEV*2]);
//                T[n*stride+m]= T[n*stride+m+maxEV*maxEV*2];
//                T[m*stride+n] = T[n*stride+m];
//            }
//        }
//    f1->mem1->rt->matrixElementsTotal += countTot;
//    f1->mem1->rt->matrixElementsZero += countLam;
//    time(&lapse_t);
//    f1->mem1->rt->buildTime += difftime(lapse_t, start_t);
//    time(&start_t);]
    qs = quantumBasisSize;
    INT_TYPE x1=0,x2=0,xi=1;
    if ( typer == 1 ){
        x1 = minStage;
        x2 = maxStage;
        xi = 1;
    }else if ( typer == -1 ){
        x1 = maxStage;
        x2 = minStage;
        xi = -1;
    }
    for ( stage = x1 ; (typer != -1 && stage <= x2) || ( typer == -1  && stage >= x2) ; stage += xi){
        qs = 0;
        for ( n = 0; n < quantumBasisSize ; n++)
            switch ( typer ){
                case 1:
                    if (   f1.tulip[usz+n].value.stage <= stage ){
                        aa[qs++] = n;
                    }
                    break;
                case 0:
                    aa[qs++] = n;
                    break;
                case -1:
                    if (   f1.tulip[usz+n].value.stage >= stage ){
                        aa[qs++] = n;
                    }
                    break;
            }
            
        for ( n = 0; n < qs ; n++)
            for ( m = 0; m < qs ; m++){
                t[n*stride + m ] = T[aa[n]*stride + aa[m] ];
                s[n*stride + m ] = S[aa[n]*stride + aa[m] ];
            }

        
        if (1){
            DCOMPLEX one=1. , zero = 0.;
            assignCores(f1, 0);
            INT_TYPE complete;
            char Job = 'V';
//            if ( 0 ){
//                cblas_zcopy(maxEV*stride , S , 1 , S+maxEV*stride , 1);
//                tzheev(0, f1, 'N', quantumBasisSize, S+maxEV*stride, stride, overlap);
//                if (  (overlap[0] > 0.) && overlap[quantumBasisSize-1]/overlap[0] < f1.rt->TOL && flag <= 2)
//                    printf(" Krylov-2\t %f \n",  overlap[quantumBasisSize-1]/overlap[0]);
//                else         if (  (overlap[0] > 0.) && overlap[quantumBasisSize-1]/overlap[0] < f1.rt->TOL && flag == 3)
//                    printf(" Krylov-3 \t %f \n",  overlap[quantumBasisSize-1]/overlap[0]);
//                else         if (  (overlap[0] > 0.) && overlap[quantumBasisSize-1]/overlap[0] < f1.rt->TOL && flag == 4)
//                    printf(" Krylov-4 \t %f \n",  overlap[quantumBasisSize-1]/overlap[0]);
//                else {
//                    printf("Linear dependent! %f\t%1.16f\n",overlap[quantumBasisSize-1],overlap[0]);
//                    return quantumBasisSize;
//                }
//
//                cblas_zcopy(maxEV*stride , S , 1 , S+maxEV*stride , 1);
//            }
            complete = tzhegv (0,f1,Job,qs,t,s,stride,ritz);
            if ( complete ){
                printf("eigensovle failed\n");
                return 0;
            }
            
#if VERBOSE
            printf("eigenSolved\n");
#endif
            fflush(stdout);
        }
        time(&lapse_t);
        //f1->mem1->rt->eigenTime += difftime(lapse_t, start_t);
        time(&start_t);
        
        {       //printf("\nHeader,Number,CEG,Linear,BODY,CLASS,WEIGHT\n");
            for ( iii = 0; iii < imin(qs,Ne) ; iii++)
            {
                f1.tulip[eigenVectors+iii].value.value = ritz[iii];
                printf("%d-Press%d:,%1.15f, %f\n",stage, iii+1, f1.tulip[eigenVectors+iii].value.value,cblas_dznrm2(qs,t+iii*stride, 1));
            }
            //    fflush(stdout);
        }
    }
    if ( quantumBasisSize != qs){
        printf("uhm,  ooops\n");
        exit(0);
    }
    cblas_zcopy(stride*stride, t, 1, T, 1);

//    if (flag == 3 ){
//        if ( Ne > quantumBasisSize ){
//            printf ("warning basis too small\n");
//            exit(0);
//        }
//        el = outputSpace;
//
//
//        INT_TYPE u;
//        printf(" foundation %d < quan %d\n", foundation, quantumBasisSize );
//        for ( u = 0; u < quantumBasisSize ;u++){
//            if ( u < foundation )
//                printf("*");
//            printf ("USZ %d ---%d)\n", u+usz, tPath(f1, usz+u));
//
//
//        }
//
////        for ( u = 0; u < Ne ; u++)
////            f1.tulip[el+u].path = -1;
//
//        INT_TYPE sp,nG= tSize(f1->body),path,fpath,ipath,i,h,xx,six,ii2,jj2,kk2,fx,nm,v;
//        double mag2,*pointers[MaxCore];
//
//        f1.tulip[eigenList].ptRank[0] = 0;
//        f1.tulip[eigenList].purpose = ptObject;
//        {
//            ipath = 0;
//            fpath = 0;
//            for ( v = 0 ; v < nG*nG*nG ; v++ )
//            {
//                ii2 = v % nG;
//                jj2 = (v/nG)%nG;
//                kk2 = (v/(nG*nG))%nG;
//                        if (tIR(f1->body,ii2, jj2 , kk2, irrep))
//                        {
//                            path = v;
//                            fx = 0;
//                            six = 0;//should be clustered..
//                            for ( xx = 0; xx < quantumBasisSize ; xx++)
//                                if ( tPath(f1, usz+xx) == path ){
//                                    if  ( xx < foundation )
//                                        if (! (fx++) ){
//                                            f1.tulip[eigenList].name = usz+xx;
//                                    //        printf("begin ");
//                                        }
//                                    six+= spins(f1, usz+xx)*part(f1, usz+xx);
//                                   // printf ("%d,%d %d\n",  usz+xx, spins(f1, usz+xx),part(f1, usz+xx));
//                                }
//                           // printf(" fx %d\n", fx);
//                            f1.tulip[eigenList].Current[0] = six;
//
//
//                            if ( CanonicalRank(f1, eigenList, 0) )
//                                for ( cmpl = 0; cmpl < spins(f1, usz) ; cmpl++)
//                                {
//                                    printf("%d:+:%d \n", path, CanonicalRank(f1, eigenList, 0));
//
//#ifdef OMP
//#pragma omp parallel for private (iii,rank,rr,h,i,nm,cmpl2,r,sp,mag2) schedule(dynamic,1)
//#endif
//                                    for ( iii = 0; iii < Ne  ; iii++)
//                                    {
//#ifdef OMP
//                                        rank = omp_get_thread_num();
//#else
//                                        rank = 0;
//#endif
//                                        pointers[rank] = myStreams(f1, canonicalBuffersC, rank);
//                                        rr = 0;
//                                        //loop over actual memory order...
//                                        for ( h = 0; h < fx ;h++)
//                                            for ( i = 0; i < quantumBasisSize ; i+=foundation)
//                                            {
//                                                nm = h+fpath + i ;
//                                                if ( nm >= quantumBasisSize ){
//                                                    printf("**%d %d %d %d **", nm,h,fpath, i );
//                                                    exit(0);
//                                                }
//
//                                              //      printf ( "%d:: %d <-",iii+1, nm );
//                                               // fflush(stdout);
//                                                for ( cmpl2 = 0; cmpl2 < spins(f1,usz + nm ) ; cmpl2++)
//                                                    for ( r = 0; r < part(f1,usz + nm ); r++){
//                                                        if ( r < CanonicalRank(f1,usz + nm,cmpl) &&  cmpl2 == cmpl && tPath (f1,usz + nm) == path){
//                                                            pointers[rank][rr++] = (vectors[cmpl]+(iii)*stride)[nm];
//                                                        }else{
//                                                            pointers[rank][rr++] = 0.0;
//                                                        }
//                                                    //    printf ( "->%d\n", rr );
//                                                    }
//                                                //fflush(stdout);
//
//                                            }
//                                     //   printf("%d %d ==%f\n", path, iii+1,cblas_dnrm2 ( rr, pointers[rank],1) );
//                                        if ( rr != CanonicalRank(f1, eigenList, 0) ){
//                                            printf("%d %d exit\n", rr, CanonicalRank(f1, eigenList, 0) );
//                                            exit(0);
//                                        }
//                                        if ( cblas_dnrm2 ( rr, pointers[rank],1) > 1e-6 ){
//                                            if ( rr > part(f1, canonicalBuffersC))
//                                            {
//                                                printf("crap!\n %lld %lld",rr , part(f1, canonicalBuffersC));
//                                                exit(0);
//                                            }
//                                            nm = ipath * Ne + iii;
//                                      //      fflush(stdout);
//
//                                            tCycleDecompostionListOneMP(rank,f1, eigenList, pointers[rank],el+nm, cmpl, f1->mem1->rt->vCANON, part(f1, el+nm), 1.);
////                                            f1.tulip[el+nm].path = path;
//                                            mag2 = 0;
//                                            for ( sp = 0; sp < spins(f1, el+nm) ; sp++)
//                                                mag2 += inner(rank, f1, el+nm, sp);
//                                          //  printf("mag %f\n", sqrt(mag2));
//
//                                            tScale(f1, el+nm, 1./sqrt(mag2));
//                               //             printf("assign %d <= %d\n", nm , path);
//                                        }else {
//
//                                        }
//                                    }
//                                }
//
//                            ipath += 1;
//                            fpath += fx;
//                        }
//        }
//        }
//        if ( fpath != foundation ){
//            printf ("wanring !!! fpath %d foundation %d", fpath, foundation);
//            fflush(stdout);
//            exit(0);
//        }
//
//        time(&lapse_t);
//        printf("\nSpectrum TIME \t %15.0f\n", difftime(lapse_t, start_t));
//        fflush(stdout);
//
//    }else
        if ( flag == 4 )
        {

            double norms,*pointers[MaxCore];
            if (1){
                f1.tulip[eigenList].name = usz;
                f1.tulip[eigenList].species = vector;
                f1.tulip[eigenList].Partition= 0;
                for( g = 0; g < quantumBasisSize ; g++)
                    for ( cmpl = 0 ;cmpl < spins(f1,usz+g) ; cmpl++)
                        for ( r = 0; r < part(f1, usz+g); r++){
                            f1.tulip[eigenList].Partition++;//memory inline. need this.
                        }
                // printf("[[%d %d %lld %lld]]\n", eigenList, usz,f1.tulip[usz].Address, streams(f1, eigenList, 0, 0 )-streams(f1, usz, 0, 0 ));
                for ( cmpl = 0 ;cmpl < spins(f1,eigenVectors) ; cmpl++)
                    
                    for ( powerMat = 1 ; powerMat <= 1 ; powerMat++){
                        
                        el = outputSpace + (powerMat-1)*Ne;
                        //                    printf("begin\n");
                        //                    fflush(stdout);
//#ifdef OMP
//#pragma omp parallel for private (iii,rank,rr,g,r,cmpl2) schedule(dynamic,1)
//#endif
                        for ( iii = 0; iii < imin(quantumBasisSize,Ne) ; iii++)
                        {
                            
//#ifdef OMP
//                            rank = omp_get_thread_num();
//#else
                            rank = 0;
//#endif
                            pointers[rank] = myStreams(f1, canonicalBuffersC, rank);
                            
                            //                        printf( "%d %d\n",iii+1, cblas_idamax(stride, (vectors[cmpl]+(iii)*stride), 1));
                            
                            rr = 0;
                            norms = 0.;
                            for( g = 0; g < quantumBasisSize ; g++){
                                for ( cmpl2 = 0; cmpl2 < spins(f1, usz+g) ; cmpl2++)
                                    for ( r = 0; r < part(f1, usz+g); r++){
                                        if (  r >= CanonicalRank(f1, usz+g, cmpl2)){
                                            pointers[rank][rr++] = 0.0;
                                        }else{
                                            if( cmpl2 == cmpl ){
                                                if ( cmpl == 0 )
                                                    pointers[rank][rr++] = creal((T + iii*stride)[g]);
                                                else if ( cmpl == 1 )
                                                    pointers[rank][rr++] = cimag((T + iii*stride)[g]);
                                            }
                                            else
                                                pointers[rank][rr++] = 0.0;//mask off other *SPIN*( or complex vector).
                                        }
                                          //  printf("%d %d %d %f\n", usz+g,cmpl2,r,myStreams(f1, canonicalBuffersC, rank)[rr-1] );
                                    }
                            }
                            
//                            for ( g = 0 ; g < 35 ; g++)
//                                printf("%1.2f,",pointers[rank][g]);
//                            printf("\n\n");
                            if ( rr > part(f1, canonicalBuffersC))
                            {
                                printf("crap!\n %d %d",rr , part(f1, canonicalBuffersC));
                                fflush(stdout);

                                exit(0);
                            }
                            tCycleDecompostionGridOneMP(-2,f1, eigenList,0, pointers[rank],el+iii   , cmpl, f1.rt->vCANON, part(f1, el+iii), f1.rt->powDecompose);
                        
//                            INT_TYPE jj;
  //                          for ( jj = 0 ;jj <= iii ; jj++)
    //                            printf("_%f_", tMultiplyMP(rank, &info, f1, 1., -1, nullName, 0, 'T', el+jj, 0,'N',el+iii, 0));
                        
                        
                        }
                    }
            }
            time(&lapse_t);
            printf("\nFinal Spectrum TIME \t %15.0f\n", difftime(lapse_t, start_t));
            fflush(stdout);
        }
    
    
    //f1.tulip[matrixHbuild].Current[0] = quantumBasisSize;
    //f1.tulip[matrixSbuild].Current[0] = quantumBasisSize;

    return 0;
}



