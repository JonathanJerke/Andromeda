/*
 *  purity.c
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

#include "purity.h"

INT_TYPE krylov ( struct calculation *c1, struct field f1)
{
    INT_TYPE EV = 0,fi;
    inputFormat(f1.f, "ham.matrix", hamiltonian, 1);
    
    
    f1.f.tulip[Ha].linkNext = hamiltonian;
    f1.f.tulip[Iterator].linkNext = hamiltonian;

    f1.i.qFloor = countLinesFromFile(c1,f1,0,&f1.i.iRank, &f1.i.xRank);
    //count canonical-rank...
    
    
    f1.i.nStates =1  ;
    f1.i.nStates =f1.i.Iterations  ;
    
    iModel(c1,&f1);
    for ( fi =0 ; fi < f1.i.files ; fi++){
        tLoadEigenWeights (c1,f1, f1.i.fileList[fi], &EV,f1.f.user , f1.i.collect);//UNUSUAL!!!
    }
    if (EV == 0 ){
        print(c1,f1,1,0,0,eigenVectors);
        return 1;
    }
    INT_TYPE RdsSize = EV,iterator=0;
    
    if(1){
        printf ("Step \t%d\n", iterator);
        INT_TYPE cmpl,g ;
        //        for ( iii = 0; iii < EV ; iii++){
        //            printf ( "\n Vector \t%d \n", iii+1);
        //            for ( sp = 0 ; sp < spins(f1.f, f1.f.user); sp++)
        //                tCycleDecompostionGridOneMP(-1, f1.f, f1.f.user+iii, sp, NULL, eigenVectors+iii, sp, c1->rt.vCANON, part(f1.f,eigenVectors+iii), -1);
        
        for ( cmpl = 0; cmpl < spins(f1.f, f1.f.user) ; cmpl++){
            tClear(f1.f,totalVector);
            if ( f1.f.cat && f1.i.filter  && f1.i.irrep )
            {
                for( g = 0; g < EV ; g++)
                    tBuild3Irr(0, f1.f, f1.i.irrep, f1.f.user+g, cmpl, totalVector, 0);
            }else {
                for( g = 0; g < EV ; g++)
                    tAddTw(f1.f, totalVector, 0, f1.f.user+g, cmpl);
            }
            tCycleDecompostionGridOneMP(-2, f1.f, totalVector, 0, NULL,eigenVectors , cmpl, c1->rt.vCANON, part(f1.f,eigenVectors), c1->rt.powDecompose);
        }
        double norm = magnitude(f1.f, eigenVectors );
        if ( norm > c1->rt.TARGET ){
            printf("Normed from %f\n", norm );
            tScaleOne(f1.f, eigenVectors, 0, 1/norm);
        }
        else
        {
            print(c1,f1,1,0,0,eigenVectors);
            return 1;
        }
        EV = 1;
        RdsSize = 1;
        tFilter(f1.f, EV, 0, eigenVectors);//classify
        printExpectationValues(f1.f, Ha, eigenVectors);
        fflush(stdout);
        print(c1,f1,1,0,1,eigenVectors);
    }
    
    INT_TYPE flag;
    
    for ( iterator = 1 ; iterator < f1.i.Iterations ; iterator++){
        
        
        
        flag = 1;
        
        
        
        if ( ! tGreatDivideIteration(flag, c1->i.shiftVector[iterator-1][0],c1->i.shiftVector[iterator-1][1],  f1.f,Iterator, 1,0,eigenVectors+RdsSize-EV,EV,2*EV,0)){
            RdsSize += EV;
            
            if(1){
                tFilter(f1.f, EV, !(!f1.i.filter )* f1.i.irrep, eigenVectors+RdsSize-EV);//filter
                printf ("Step \t%d\n", iterator);
                fflush(stdout);
                INT_TYPE iii ;
                for ( iii = 0; iii < EV ; iii++){
                    printf ( "\n Vector \t%d \t %d\n", iii+1, +RdsSize-EV+iii);
                    printExpectationValues(f1.f, Ha, eigenVectors+RdsSize-EV+iii);
                    fflush(stdout);
                    print(c1,f1,0,RdsSize-EV+iii,RdsSize-EV+iii+1,eigenVectors);
                }
            }
        }else {
            break;
        }
        fflush(stdout);
    }
    fModel(&f1.f);
    
    return 0;
}


INT_TYPE purityI ( struct sinc_label f1 ){
    enum division p = f1.purity, o = f1.purityOverlap,t = f1.temp;
    INT_TYPE complete,l,ll,N1,PP,ppp,pp,oo=0,r,rr,k,kk,space,R = f1.purityCanon,rank=0;
    PP = part(f1,f1.purity);
   double * w= myStreams(f1, canonicalBuffersBM, 0);
  // double w[1000];
	 double x=5,m = -5;
    for ( r = 0 ; r < R ; r++)
        for ( space = 0; space < SPACE ; space++){
            N1 = vectorLen(f1, space);
            struct name_label u = f1.tulip[o+R*r+r];
         //	printf("%d %d %d \n", R,part(f1,p+R*r+r),part(f1, o+R*r+r));
	//	fflush(stdout);
	   cblas_dcopy(N1*N1, streams(f1,hamiltonian,0,space)+r*N1*N1, 1, streams(f1,o+R*r+r,0,space)+1*N1*N1, 1);//
         // printf("%d\n", alloc(f1,hamiltonian,space));
//fflush(stdout);
complete =   tdsyev(0, f1, 'V', N1, streams(f1,o+R*r+r,0,space)+1*N1*N1, N1, w);
     

 f1.tulip[p+R*r+r].Current[0] = 1;
            for ( l = 0 ; l < N1 ; l++){
                for ( ll = 0 ; ll < N1 ; ll++)
                    if ( l == ll )
                        streams(f1,p+R*r+r,0,space)[ll*N1+l] = (w[l]-m/R)/pow(x-m,1);
                  else
                       streams(f1,p+R*r+r,0,space)[ll*N1+l] = 0.;
 		printf("%d %d %f\n",l, space,w[l]);
	           }
}
   	if (0) 
	   for ( l = 0; l < N1 ; l++){
                if ( x < w[l])
                    x = w[l];
                if ( m > w[l])
                    m = w[l];
            }
            if ( complete ){
                printf("eigensovle failed\n");
                exit(0);
            }
        
    if(0)
    for ( r = 0 ; r < R ; r++)
        for ( space = 0; space < SPACE ; space++){
            N1 = vectorLen(f1, space);

            f1.tulip[p+R*r+r].Current[0] = 1;
            for ( l = 0 ; l < N1 ; l++){
                for ( ll = 0 ; ll < N1 ; ll++)
                    if ( l == ll )
                        streams(f1,p+R*r+r,0,space)[ll*N1+l] = (w[l]-delta(r)*m)/(x-m);
                    else
                        streams(f1,p+R*r+r,0,space)[ll*N1+l] = 0.;
            }
        }
    //exit(0);
    for ( r = 0 ; r < R ; r++)
        for ( rr = 0 ; rr < R ; rr++){
            f1.tulip[o+R*rr+r].Current[0] = 1;

            for ( space = 0; space < SPACE ; space++){
                N1 = vectorLen(f1, space);
                
                for ( l = 0 ; l < N1 ; l++)
                    for ( ll = 0 ; ll < N1 ; ll++)
                        streams(f1,o+R*rr+r,0,space)[ll*N1+l] = cblas_ddot(N1, streams(f1,o+R*r+r,0,space)+1*N1*N1+l*N1, 1, streams(f1,o+R*rr+rr,0,space)+1*N1*N1+ll*N1, 1);
            }
        }
    return 0;
}

INT_TYPE purityA(struct sinc_label f1){
    enum division p = f1.purity, o = f1.purityOverlap,t = f1.temp;
    INT_TYPE PP,ppp,pp,oo=0,r,rr,k,kk,space,R = f1.purityCanon,rank=0;
    PP = part(f1,f1.purity);
    INT_TYPE xx = part(f1,copyTwo);
    zero(f1, copyTwo, 0);
    for ( r = 0 ; r < R ; r++)
        for ( rr = 0 ; rr < R ; rr++){
            for ( k = 0 ; k < R ; k++)
                for ( kk = 0 ; kk < R ; kk++)
                    for ( space = 0 ; space < SPACE ;space++)
                        for ( pp = 0 ; pp < CanonicalRank(f1, p+R*rr+k, 0);pp++){
                            tGEMM(rank, f1, space, copy,0, 0, o+R*k+kk, oo, 0, p+R*rr+k, pp, 0);
                            for ( ppp = 0 ; ppp < CanonicalRank(f1, p+kk*R+r, 0);ppp++){
                                tGEMM(rank, f1, space, copyTwo, (kk*R+k)*PP*PP + (pp*PP+ppp),0, p+kk*R+r, ppp, 0, copy, 0, 0);
                            }
                        }
            f1.tulip[copyTwo].Current[0] = PP*PP*R*R;
            tCycleDecompostionGridOneMP(-1, f1, copyTwo, 0, NULL, t+R*rr+r, 0, f1.rt->CANON, PP, 1);
        }
    
    for ( r = 0 ; r < R ; r++)
        for ( rr = 0 ; rr < R ; rr++)
            tEqua(f1, p+R*rr+r, 0, t+R*rr+r, 0);
    return 0;
}

INT_TYPE purityB(struct sinc_label f1){
    enum division p = f1.purity, o = f1.purityOverlap,t = f1.temp;
    INT_TYPE PP,ppp,pp,oo=0,r,rr,k,kk,space,R = f1.purityCanon,rank=0;
    PP = part(f1,f1.purity);
    INT_TYPE xx = part(f1,copyTwo);
    zero(f1, copyTwo, 0);
    for ( r = 0 ; r < R ; r++)
        for ( rr = 0 ; rr < R ; rr++){
            tEqua(f1, copyTwo, 0, p+R*rr+r, 0);
            tScaleOne(f1, copyTwo, 0, -2);
            for ( k = 0 ; k < R ; k++)
                for ( kk = 0 ; kk < R ; kk++)
                    for ( space = 0 ; space < SPACE ;space++)
                        for ( pp = 0 ; pp < CanonicalRank(f1, p+R*rr+k, 0);pp++){
                            tGEMM(rank, f1, space, copy,0, 0, o+R*k+kk, oo, 0, p+R*rr+k, pp, 0);
                            for ( ppp = 0 ; ppp < CanonicalRank(f1, p+kk*R+r, 0);ppp++){
                                tGEMM(rank, f1, space, copyTwo, (kk*R+k)*PP*PP + (pp*PP+ppp),0, p+kk*R+r, ppp, 0, copy, 0, 0);
                            }
                        }
            f1.tulip[copyTwo].Current[0] = PP*PP*R*R;
            tCycleDecompostionGridOneMP(-1, f1, copyTwo, 0, NULL, t+R*rr+r, 0, f1.rt->CANON, PP, 1);
            tScaleOne(f1, t+R*rr+r, 0, -1.);
        }
    
    for ( r = 0 ; r < R ; r++)
        for ( rr = 0 ; rr < R ; rr++)
            tEqua(f1, p+R*rr+r, 0, t+R*rr+r, 0);
    return 0;
}


double purityTr ( struct sinc_label f1){
    enum division p = f1.purity, o = f1.purityOverlap;
    double sum = 0.;
    INT_TYPE rank = 0,space,pp,oo=0,r,rr,R = f1.purityCanon,PP = part(f1,f1.purity);
    for ( r = 0 ; r < R ; r++)
        for ( rr = 0 ; rr < R ; rr++){
          //  sum = 0;
            for ( pp = 0 ; pp < CanonicalRank(f1, p+R*r+rr, 0);pp++){
                f1.tulip[copy].Current[0] = 0;
                for ( space = 0 ; space < SPACE ;space++)
                    tGEMM(rank, f1, space, copy,0, 0, o+R*r+rr, oo, 0, p+R*r+rr, pp, 0);
                f1.tulip[copy].Current[0] = 1;

                sum += traceOne(f1,copy,0);
            }
          //  printf("%d %d %f\n",r,rr,sum);
        }
    return sum;
}


int main (INT_TYPE argc , char * argv[]){
    struct calculation c;
    struct field f;
    bootShell(argc, argv,&c,&f);
    // test2();
    
    INT_TYPE space,i,a,plusSize,nStatesTrans=0,nStatesFound=0 ,RdsSize = 0,totalIter = 0;
    FILE * out = stdout;
    struct runTime * rt = & c.rt;
    f.f.rt = rt;
    
    
#ifdef OMP
    if ( c.i.omp > MaxCore ){
        printf("lanes > MaxCore\n");
        exit(0);
    }
    rt->NLanes = c.i.omp;
#pragma omp parallel for private (i)
    for ( i = 0; i < MaxCore ; i++){
        rt->NSlot = omp_get_num_threads();
    }
    if ( rt->NLanes > rt->NSlot ){
        printf("decrease lanes or increase your number of OMP cores\n");
        exit(0);
    }
    
#ifdef MKL
    if ( rt->NSlot < c.i.mkl * c.i.omp )
    {
        printf("not enough slots for mkl*omp\n" );
        exit(0);
    }
    else
    {
        rt->NParallel = c.i.mkl;
    }
#endif
#endif
    
    if ( c.rt.phaseType == buildFoundation ){
        enum division leftP = Ha;
        iModel(&c,&f);
        tClear(f.f, hamiltonian );
        
        tAddTw(f.f, copyTwo, 0, kinetic, 0);
        tAddTw(f.f, copyTwo, 0, linear, 0);
        
        tCycleDecompostionGridOneMP(-1, f.f, copyTwo, 0, NULL, hamiltonian, 0, c.rt.CANON, f.f.purityCanon, 1);
        
        INT_TYPE EV;
        iModel(&c,&f);
        tBoot1Construction(&c,f.f ,build);
        tSortBoot(&c,f.f,build);
        EV =   tSlam(f.f,f.i.qFloor,f.f.user,c.i.level);
        
        
        purityI(f.f);
        int i;       for(i =0; i < 100;i++){
            printf("%f** \n\n", purityTr(f.f));
            if (purityTr(f.f)> 2){
                // printf("%f \n", purityTr(f.f));
                purityA(f.f);
                //printf("%f \n", purityTr(f.f));
            }else{
                purityB(f.f);
                
                // printf("%f \n", purityTr(f.f));
            }
            purityI(f.f);

        }
    }
    return 0;
}
