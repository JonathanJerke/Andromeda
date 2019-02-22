/*
 *  coreUtil.c
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

#include "coreUtil.h"

INT_TYPE spaces( struct field * f1, enum division label){
    if ( bodies(f1,label) == four)
        return SPACE;
    else
        return SPACE;
}


INT_TYPE name ( struct field * f1, enum division label){
    return f1->sinc.tulip[f1->sinc.tulip[f1->sinc.tulip[label].name].name].name;
}

INT_TYPE tPath ( struct field * f1, enum division label){
    return f1->sinc.tulip[label].path;
}


INT_TYPE part ( struct field * f1 , enum division label ){
    
    if ( label > f1->sinc.end ){
        printf("rank past end\n");
    }
    //if ( purpose(f1, label) == Object )
    return f1->sinc.tulip[name(f1,label)].Partition;
    //else if ( purpose(f1, label) == ptObject)
    //    return f1->sinc.tulip[name(f1,label)].Partition-f1->sinc.tulip[label].ptRank[0];
    //else
    //   return 0;
}


INT_TYPE species ( struct field * f1 , enum division label ){
    return f1->sinc.tulip[name(f1,label)].species;
}

enum body bodies ( struct field * f1 , enum division label ){
    return f1->sinc.tulip[name(f1,label)].NBody;
}

enum bodyType bodyType ( struct field * f1 , enum division label ){
    return f1->sinc.tulip[name(f1,label)].TBody;
}


INT_TYPE header ( struct field * f1 , enum division label ){
    return f1->sinc.tulip[name(f1,label)].header;
}

INT_TYPE purpose ( struct field * f1, enum division label){
    return f1->sinc.tulip[label].purpose;
}
INT_TYPE memory ( struct field * f1, enum division label){
    return f1->sinc.tulip[name(f1,label)].memory;
}

//void initPointerTensors(struct field *f1){
//    INT_TYPE sp;
//    f1->sinc.tulip[v2].name = vectors;
//    f1->sinc.tulip[v2].purpose = ptObject;
//    
//    f1->sinc.tulip[u2].name = vectors;
//    f1->sinc.tulip[u2].purpose = ptObject;
//    tClear(f1, v2);
//    tClear(f1,u2);
//    tClear(f1, u);
//    tClear(f1,v);
//    tClear(f1, p);
//    tClear(f1,pp);
//    tClear(f1,pointer);
//    
//    
//    f1->sinc.tulip[u].purpose = ptObject;
//    f1->sinc.tulip[v].purpose = ptObject;
//    f1->sinc.tulip[u].Partition = 1;
//    f1->sinc.tulip[v].Partition = 1;
//    f1->sinc.tulip[u].parallel = 1;
//    f1->sinc.tulip[v].parallel = 1;
//    
//    f1->sinc.tulip[u].species = matrix;
//    f1->sinc.tulip[v].species = matrix;
//    
//    f1->sinc.tulip[p].purpose = ptObject;
//    f1->sinc.tulip[pp].purpose = ptObject;
//    f1->sinc.tulip[q].purpose = ptObject;
//    f1->sinc.tulip[qq].purpose = ptObject;
//    
//    f1->sinc.tulip[pointer].purpose = ptObject;
//    f1->sinc.tulip[pointer].parallel = 1;
//    
//    for ( sp = 0 ; sp < spins(f1, diis) ; sp++){
//        f1->sinc.tulip[p].Current[sp] = f1->sinc.tulip[project].Partition;
//        f1->sinc.tulip[pp].Current[sp] = f1->sinc.tulip[project].Partition;
//        f1->sinc.tulip[q].Current[sp] = f1->sinc.tulip[project].Partition;
//        f1->sinc.tulip[qq].Current[sp] = f1->sinc.tulip[project].Partition;
//        f1->sinc.tulip[u].Current[sp] = 1;
//        f1->sinc.tulip[v].Current[sp] = 1;
//    }
//    
//    f1->sinc.tulip[p].species = matrix;
//    f1->sinc.tulip[pp].species = matrix;
//    f1->sinc.tulip[p].Partition = f1->sinc.tulip[project].Partition;
//    f1->sinc.tulip[pp].Partition = f1->sinc.tulip[project].Partition;
//    f1->sinc.tulip[q].Partition = f1->sinc.tulip[project].Partition;
//    f1->sinc.tulip[qq].Partition = f1->sinc.tulip[project].Partition;
//    f1->sinc.tulip[q].species = matrix;
//    f1->sinc.tulip[qq].species = matrix;
//    
//}

//INT_TYPE defineSpinors (struct field *f1 ){
//    INT_TYPE si,sp,sp2,s;
//
//    for ( s = 0; s < NspinType ; s++)
//        for ( sp = 0 ; sp < NS ; sp++)
//            for (sp2 = 0 ; sp2 < NS ; sp2++)
//                f1->sinc.arraySpin[s][sp][sp2] = -1;;
//
//    si = 0;
//
//
//
//
//#if 0
//    for ( sp = 0 ; sp < NS ; sp++){
//        for (sp2 = 0 ; sp2 < NS ; sp2++){
//            if ( sp <= sp2 ){
//                f1->sinc.arraySpin[full][sp][sp2] = si++;
//                f1->sinc.arraySpin[full][sp2][sp] = f1->sinc.arraySpin[full][sp][sp2];
//            }
//            f1->sinc.arraySpin[coulomb][sp][sp2] = 0;
//        }
//        f1->sinc.arraySpin[none][sp][sp] = 0;
//        f1->sinc.arraySpin[sym][sp][sp] = 0;
//        f1->sinc.arraySpin[diag][sp][sp] = sp;
//    }
//#else
//
////    if ( NS == 2 ){
////        f1->sinc.arraySpin[full][0][0] = 0;
////        f1->sinc.arraySpin[full][1][1] = 1;
////        f1->sinc.arraySpin[full][1][0] = 2;
////        f1->sinc.arraySpin[full][0][1] = 2;
////
////
////        f1->sinc.arraySpin[coulomb][0][0] = 0;
////        f1->sinc.arraySpin[coulomb][1][0] = 0;
////        f1->sinc.arraySpin[coulomb][0][1] = 0;
////        f1->sinc.arraySpin[coulomb][1][1] = 0;
////
////
////        f1->sinc.arraySpin[none][0][0] = 0;
////        f1->sinc.arraySpin[none][1][1] = 0;
////
////        f1->sinc.arraySpin[sym][0][0] = 0;
////        f1->sinc.arraySpin[sym][1][1] = 0;
////
////        f1->sinc.arraySpin[diag][0][0] = 0;
////        f1->sinc.arraySpin[diag][1][1] = 1;
////        return 0;
////
////    }
//
//
////    printf("incompatible spins\n NS %d\n\n",NS);
////    exit(0);
//    return 1;
//
//#endif
//
//
//
//
//    return 0;
//
//}
INT_TYPE * vectorLen ( struct field * f1 , enum division label ){
    
    if ( bodies(f1,label) == one){
        if ( header(f1, label ) == Rds)
            return f1->rds.Basis;
//        else if ( header(f1, label ) == Band)
//            return f1->band.Basis;
        else
            return f1->sinc.Basis;
    }
    
    if ( bodies(f1,label) == two ){
        if ( header(f1, label ) == Rds)
            return f1->rds.Basis2;
//        else
//            if ( header(f1, label ) == Band)
//                return f1->band.Basis2;
            else
                return f1->sinc.Basis2;
    }
    if ( bodies(f1,label) == three){
        if ( header(f1, label ) == Rds)
            return f1->rds.Basis3;
//        else
//            if ( header(f1, label ) == Band)
//                return f1->band.Basis3;
            else
                return f1->sinc.Basis3;
    }
    if ( bodies(f1,label) == four){
        if ( header(f1, label ) == Rds)
            return f1->rds.Basis4;
//        else
//            if ( header(f1, label ) == Band)
//                return f1->band.Basis4;
            else
                return f1->sinc.Basis4;
    }
    
    printf("vectorLen %d %d\n",bodies(f1,label) ,label);
    exit(0);
}

INT_TYPE length ( struct field * f1 , enum division label, INT_TYPE *lens ){
    
    INT_TYPE c=0, c1 = species(f1, label), *vc= vectorLen(f1,label),space;
    
    if ( lens != NULL )
        for ( space = 0 ;space < SPACE ; space++){
            lens[space] = 1;
            for( c = 0 ; c < c1  ; c++)
                lens[space] *= vc[space];
        }
    
    return c1;
}

INT_TYPE alloc ( struct field * f1 , enum division label ){
    INT_TYPE c=0, c1 = species(f1, label), *vc= vectorLen(f1,label),space,lens[SPACE],maxl = 0;;
    
    for ( space = 0 ;space < SPACE ; space++){
        lens[space] = 1;
        for( c = 0 ; c < c1  ; c++)
            lens[space] *= vc[space];
        if ( lens[space] > maxl )
            maxl = lens[space];
    }
    return maxl;
}

INT_TYPE zero ( struct field * f1 , enum division label, INT_TYPE spin ){
    //f1->sinc.tulip[label].Current = 0;
    INT_TYPE i, space,M2[3];
    length(f1, label, M2);
    
    
    
    for ( space = 0; space < SPACE ;space++)
        for ( i = 0; i < M2[space]*part(f1, label) ; i++ ){
            *(streams(f1,label, (spin) , space )+i) = 0.;
        }
    return 0 ;
}

INT_TYPE myZero ( struct field * f1 , enum division label, INT_TYPE spin ){
    //f1->sinc.tulip[label].Current = 0;
    INT_TYPE i, M2[3];
    length(f1, label, M2);
    
    
    
    for ( i = 0; i < M2[0]*part(f1, label) ; i++ ){
        *(myStreams(f1,label, (spin)  )+i) = 0.;
    }
    return 0 ;
}


INT_TYPE tClear ( struct field * f1 , enum division label ){
    INT_TYPE spin ;
    for ( spin = 0; spin < spins(f1,label) ; spin++){
        f1->sinc.tulip[label].Current[spin] = 0;
    }
    
    return 0 ;
}

INT_TYPE CanonicalRank( struct field * f1 , enum division label , INT_TYPE spin ){
    if ( label > f1->sinc.end ){
        printf("Can rank past end\n");
    }
    
    
    if ( spin < spins(f1, name(f1,label)) ){
        
        return f1->sinc.tulip[label].Current[spin];
    }
    else {
        printf( " *****rank : spin out %d (%lld)( %lld >= %lld )\n", label,name(f1, label), spin , spins(f1, name(f1,label)) );
        fflush(stdout);
       // exit(0);
        return 0;
    }
}

INT_TYPE Rank( struct field * f1 , enum division label ){
    INT_TYPE sp,ra=0;;
    for ( sp = 0 ; sp < spins(f1, label);sp++)
        ra += CanonicalRank(f1, label, sp);
    return ra;
}


double getPosition(struct field * f1, INT_TYPE at , INT_TYPE space ){
    return f1->atoms[at].position[space];
}

INT_TYPE spins ( struct field * f1 , enum division label ){    
    enum spinType sp = f1->sinc.tulip[label].spinor;

#ifndef APPLE
#ifdef OMP
    if (f1->sinc.tulip[label].parallel /*== 1*/ )
            return  imax(1,imax(f1->mem1->rt->NCore,(sp == cmpl)*2));
#else
    if (f1->sinc.tulip[label].parallel /*== 1*/ )
        return  imax(1,imax(1,(sp == cmpl)*2));
#endif
//        if (f1->sinc.tulip[label].parallel == 2 )
//            return  imax(f1->mem1->rt->NTraffic,(sp == cmpl)*2);
#endif
    
    
    if ( sp == none )
        return 1;
    else if ( sp == cmpl )
        return 2;
    
    return 0;
//    INT_TYPE ref = 0,mult= 1;
//    }
//    //    if (f1->sinc.tulip[label].permutation == 1 )
//    //        ref =  6;//override parallel..
//    //    if( f1->sinc.tulip[label].scaling == linear )
//    //        mult = f1->Na+1;
//    //    else if ( f1->sinc.tulip[label].scaling == none )
//    //        mult = 1;
//
//    enum spinType sp = f1->sinc.tulip[label].spinor;
//    if ( sp != none ){
//
//    }
//    if ( species(f1,label ) == vector && !(sp == none  )){
//        return mult*imax(NS,ref);
//    }
//    if ( sp == none   ){
//        return mult*imax(1,ref);
//    }
//
//
//    printf("spins down\n");
//    exit(0);
//    return 0;
}


double sumSquare (struct field * f1,  enum division alloy){
    INT_TYPE l,dim ;
    //  double value0 = sqrt(inner(f1,alloy,alloy));
    double norm =0.,product;
    INT_TYPE M2[SPACE];
    length(f1,alloy,M2);
    INT_TYPE iOne = 1,sp;
    for ( sp = 0 ; sp < spins(f1,alloy);sp++)
        for ( l = 0; l < CanonicalRank(f1, alloy,sp) ;l++){
            product = 1.;
            for ( dim = 0; dim < SPACE ; dim++)
                product *= sqr(cblas_dnrm2(M2[dim], streams(f1, alloy,sp,dim)+l*M2[dim],iOne));;
            norm += product ;//* sqr(coeff(f1, alloy , sp,l ));
        }
    return norm;
}


ADDRESS_TYPE  fromBegining( struct field * f1 , enum division label ){
    if ( label > f1->sinc.end ){
        printf("from %d past end\n",label);
    }
    if ( memory (f1, label ) == oneObject ){
        //  printf("-\n");
        printf("ATTEMPTED TO LOAD TRIPLE FROM SINGLE %d\n",label);

        exit(0);
        return f1->sinc.tulip[label].Address ;//3dim untouched
    }else
    if ( label == f1->sinc.end ){
        return f1->sinc.tulip[label].Address;
    }

    else{
        INT_TYPE leng = alloc(f1, label);
        
       // printf("%d %lld -> %lld \n", label,f1->sinc.tulip[label].Address,f1->sinc.tulip[label].Address+ leng*part(f1,label)*spins(f1,label));
        
        return f1->sinc.tulip[label].Address + leng*part(f1,label)*spins(f1,label);
    }
}

ADDRESS_TYPE  fromMyBegining( struct field * f1 , enum division label  ){
    if ( label > f1->sinc.end ){
        printf("from %d past end\n",label);
    }
    if ( memory (f1, label ) == oneObject ){
        INT_TYPE leng = alloc(f1, label);
        //        printf("%lld\n", leng);
        //        fflush(stdout);
        //        printf("%lld\n", spins(f1,label));
        //        fflush(stdout);
        //        printf("%lld\n", part(f1,label));
        //        fflush(stdout);
        
         //printf("%d %lld my %lld \n", label,f1->sinc.tulip[label].myAddress,f1->sinc.tulip[label].myAddress+ leng*part(f1,label)*spins(f1,label));
        
        return f1->sinc.tulip[label].myAddress + leng*part(f1,label)*spins(f1,label);
    }
    if ( label == f1->sinc.end ){
        return f1->sinc.tulip[label].myAddress;
    }
    else{
        printf("ATTEMPTED TO LOAD SINGLE FROM TRIPLE %d\n",label);
        exit(0);
        return f1->sinc.tulip[label].myAddress ;
    }
}



Stream_Type* myStreams ( struct field * f1, enum division label ,INT_TYPE spin ){
    if ( memory(f1, label) != oneObject){
        printf("called my stream %d %d\n",label,purpose(f1, label) );
        exit(0);
    }
    
    
    if ( spin < 0 || spin >= spins(f1, label)){
        printf("\n*my %d %lld\n\n", label, spin );
        exit(0);
    }
    if ( spin != 0  ){
        INT_TYPE leng = alloc(f1, name(f1,label));
        INT_TYPE partit = part(f1, name(f1,label));
        if ( purpose(f1, label ) == ptObject ){
            
            printf("fixme\n");
            exit(0);
        }
        return f1->sinc.rose[SPACE].stream+f1->sinc.tulip[name(f1,label)].myAddress + leng * partit * spin  ;
        
    }
    
    else{
        if ( purpose(f1, label ) == ptObject ){
            INT_TYPE len[SPACE];
            length(f1, name(f1,label), len);
            
            ///        printf("%lld -%lld %lld %lld\n", label,name(f1,label),len[space],f1->sinc.tulip[label].ptRank[spin]);
            // / printf("%lld %lld\n", len[space]*part(f1,label), len[space]*f1->sinc.tulip[label].ptRank[spin]);
            // p/rintf("%lld\n", f1->sinc.tulip[name(f1,label)].Address);
            return f1->sinc.rose[SPACE].stream+f1->sinc.tulip[name(f1,label)].myAddress + len[0]*f1->sinc.tulip[label].ptRank[spin] ;
            
        }else {
            
            return f1->sinc.rose[SPACE].stream+f1->sinc.tulip[name(f1,label)].myAddress ;
        }
    }
    
}


Stream_Type* streams ( struct field * f1, enum division label ,INT_TYPE spin, INT_TYPE space ){
    if ( memory(f1, label) == oneObject){
        return myStreams(f1, label, spin );
    }
    
    if ( spin < 0 || spin >= spins(f1, label)){
        printf("\n* %d %lld %lld\n\n", label, spin , space);
        exit(0);
    }
    if ( space < 0 || space >= SPACE ){
        printf("streams: space out %d %lld\n", label, space);
        exit(0);
    }
    if ( spin != 0  ){
        INT_TYPE leng = alloc(f1, name(f1,label));
        INT_TYPE partit = part(f1, name(f1,label));
        
        if ( purpose(f1, label ) == ptObject ){
            INT_TYPE len[SPACE];
            length(f1, name(f1, label), len);
       //     printf(".");
            return f1->sinc.rose[space].stream+f1->sinc.tulip[name(f1,label)].Address + leng * partit * spin + len[space]*f1->sinc.tulip[label].ptRank[spin] ;
        }
        else
            return f1->sinc.rose[space].stream+f1->sinc.tulip[name(f1,label)].Address + leng * partit * spin  ;
        
    }
    
    else{
        
        if ( purpose(f1, label ) == ptObject ){
            INT_TYPE len[SPACE];
            length(f1, name(f1,label), len);
            
//                    printf("%lld -%lld %lld %lld\n", label,name(f1,label),len[space],f1->sinc.tulip[label].ptRank[spin]);
//             printf("%lld %lld\n", len[space]*part(f1,label), len[space]*f1->sinc.tulip[label].ptRank[spin]);
//printf("%lld %lld\n", f1->sinc.tulip[name(f1,label)].Address,f1->sinc.tulip[interactionExchange].Address);
            return f1->sinc.rose[space].stream+f1->sinc.tulip[name(f1,label)].Address + len[space]*f1->sinc.tulip[label].ptRank[spin] ;
            
        }else {
            //    printf("add %lld\n", f1->sinc.tulip[label].Address);
            
            return f1->sinc.rose[space].stream+f1->sinc.tulip[name(f1,label)].Address ;
        }
        
    }
    

    if ( memory(f1, label) == oneObject){
        return myStreams(f1, label, spin );
    }
    
    if ( spin < 0 || spin >= spins(f1, label)){
        printf("\n* %d %lld %lld\n\n", label, spin , space);
        exit(0);
    }
    if ( space < 0 || space >= SPACE ){
        printf("streams: space out %d %lld\n", label, space);
        exit(0);
    }
    if ( spin != 0  ){
        INT_TYPE leng = alloc(f1, name(f1,label));
        INT_TYPE partit = part(f1, name(f1,label));
        
        if ( purpose(f1, label ) == ptObject ){
            INT_TYPE len[SPACE];
            length(f1, name(f1, label), len);
            
            return f1->sinc.rose[space].stream+f1->sinc.tulip[name(f1,label)].Address + leng * partit * spin + len[space]*f1->sinc.tulip[label].ptRank[spin] ;
        }
        else
            return f1->sinc.rose[space].stream+f1->sinc.tulip[name(f1,label)].Address + leng * partit * spin  ;
        
    }
    
    else{
        
        if ( purpose(f1, label ) == ptObject ){
            INT_TYPE len[SPACE];
            length(f1, name(f1,label), len);
            
            ///        printf("%lld -%lld %lld %lld\n", label,name(f1,label),len[space],f1->sinc.tulip[label].ptRank[spin]);
            // / printf("%lld %lld\n", len[space]*part(f1,label), len[space]*f1->sinc.tulip[label].ptRank[spin]);
            // p/rintf("%lld\n", f1->sinc.tulip[name(f1,label)].Address);
            return f1->sinc.rose[space].stream+f1->sinc.tulip[name(f1,label)].Address + len[space]*f1->sinc.tulip[label].ptRank[spin] ;
            
        }else {
            //    printf("add %lld\n", f1->sinc.tulip[label].Address);
            
            return f1->sinc.rose[space].stream+f1->sinc.tulip[label].Address ;
        }
        
    }
    
}

//category is equal to atomic count , begins at 1

enum division rivers(INT_TYPE rank, struct field * f1, enum division A, INT_TYPE category){
    f1->sinc.tulip[lanes+rank].purpose = ptObject;
    f1->sinc.tulip[lanes+rank].parallel = 0;// = species(f1, A);;
    tClear(f1, lanes+rank);
    f1->sinc.tulip[lanesc+rank].spinor =f1->sinc.tulip[A].spinor;

    f1->sinc.tulip[lanes+rank].name = name(f1,A);
    f1->sinc.tulip[lanes+rank].blockType = f1->sinc.tulip[A].blockType;
    INT_TYPE spin;
    for ( spin = 0 ;spin < spins(f1, A); spin++)
        if ( category == 0 ){
            f1->sinc.tulip[lanes+rank].ptRank[spin] = 0;
            f1->sinc.tulip[lanes+rank].Current[spin] = f1->sinc.tulip[name(f1,A)].stop[spin][category];
#if VERBOSE
            //            if ( A != oneBody)
            printf("%d %d r%d [%d %d] -- %d %d\n",name(f1,A),category,rank, f1->sinc.tulip[lanes+rank].ptRank[spin],f1->sinc.tulip[lanes+rank].Current[spin],f1->sinc.tulip[name(f1,A)].stop[spin][category],0);
#endif
        }
        else{
            f1->sinc.tulip[lanes+rank].ptRank[spin] = f1->sinc.tulip[name(f1,A)].stop[spin][category-1];
            f1->sinc.tulip[lanes+rank].Current[spin] = f1->sinc.tulip[name(f1,A)].stop[spin][category]-f1->sinc.tulip[name(f1,A)].stop[spin][category-1];
#if VERBOSE
            //            if ( A != oneBody)
                printf("%d %d r%d [%d %d] -- %d %d\n",name(f1,A),category,rank, f1->sinc.tulip[lanes+rank].ptRank[spin],f1->sinc.tulip[lanes+rank].Current[spin],f1->sinc.tulip[name(f1,A)].stop[spin][category],f1->sinc.tulip[name(f1,A)].stop[spin][category-1]);
#endif
        }
    
    return lanes+rank;
}

INT_TYPE riversBranch(INT_TYPE rank, struct field * f1, enum division A, INT_TYPE spin, INT_TYPE category){
    INT_TYPE pt , cur;
    
    
    if ( category == 0 ){
        pt = 0;
        cur = f1->sinc.tulip[A].stop[spin][category];
        
    }
    else{
        pt = f1->sinc.tulip[A].stop[spin][category-1];
        cur = f1->sinc.tulip[A].stop[spin][category]-f1->sinc.tulip[A].stop[spin][category-1];
    }
    // printf("%d-> %d\n", A,cur);
    return cur ;
}

enum division ocean(INT_TYPE rank, struct field * f1, enum division A, INT_TYPE l, INT_TYPE spin){
    f1->sinc.tulip[lanesc+rank].purpose = ptObject;
    //    f1->sinc.tulip[lanes+rank].parallel = 0;// = species(f1, A);;
    tClear(f1, lanesc+rank);
    f1->sinc.tulip[lanesc+rank].spinor =f1->sinc.tulip[A].spinor;
    f1->sinc.tulip[lanesc+rank].name = name(f1,A);
    f1->sinc.tulip[lanesc+rank].ptRank[spin] = l+f1->sinc.tulip[A].ptRank[spin];
    f1->sinc.tulip[lanesc+rank].Current[spin] = 1;
    f1->sinc.tulip[lanesc+rank].parallel = f1->sinc.tulip[A].parallel;
    f1->sinc.tulip[lanesc+rank].header = header(f1 ,A);
    f1->sinc.tulip[lanesc+rank].species = species(f1 ,A);
    f1->sinc.tulip[lanesc+rank].NBody = bodies(f1 ,A);
    f1->sinc.tulip[lanesc+rank].memory = f1->sinc.tulip[A].memory;
    f1->sinc.tulip[lanesc+rank].blockType = f1->sinc.tulip[A].blockType;
    return lanesc+rank;
}

double xEqua ( struct field * f1 , enum division targ ,INT_TYPE tspin,struct field * f2 , enum division orig,INT_TYPE ospin ){
    INT_TYPE space;
    INT_TYPE eb = CanonicalRank(f2,orig,ospin);
    INT_TYPE M2[3];
    length(f2, orig, M2);
    
    
    if ( (part(f1,targ) < eb && purpose(f1,targ) != ptObject) || species(f1, targ) != species(f2, orig))
        
    {
        printf("%lld < %lld\n", part(f1,targ),eb);
        printf("%d (%lld)= %d (%lld)\n", targ,tspin, orig,ospin);
        printf("s%lld %lld %lld\n", species(f1, targ) , species(f2, orig),purpose(f1,targ));
        printf("M2 %lld %lld %lld\n", M2[0],M2[1],M2[2]);
        printf("Ranks %lld\n", eb);
        printf("tEqual.. memory\n");
        exit(0);
    }
    
    if ( f1->sinc.N1 != f2->sinc.N1 ){
        if ( bodies(f2, orig) == two )
            xTwoBand(f2, orig,  ospin, f1, targ,tspin,f1->mem1->rt->runFlag );
        else if ( bodies(f2, orig) == three )
            xThreeBand(f2, orig,  ospin, f1, targ,tspin,f1->mem1->rt->runFlag );
        else if ( bodies(f2, orig) == four )
            xFourBand(f2, orig,  ospin, f1, targ,tspin,f1->mem1->rt->runFlag );
    }
    else {
    
    
    
    for ( space = 0; space < SPACE; space++){
        cblas_dcopy(eb*M2[space], streams(f2,orig,ospin,space),1,streams(f1,targ,tspin,space),1);
    }
    }
    f1->sinc.tulip[targ].Current[tspin] = f2->sinc.tulip[orig].Current[ospin];
    f1->sinc.tulip[name(f1,targ)].header = header(f2, name(f2,orig));
    
    return 0;
}


double tEqua ( struct field * f1 , enum division targ ,INT_TYPE tspin, enum division orig,INT_TYPE ospin ){
    return xEqua(f1, targ,tspin, f1, orig,ospin);
}

INT_TYPE tEquals( struct field * f1 , enum division left , enum division right){
    enum spinType spl,spr;
    spl = f1->sinc.tulip[left].spinor;
    spr = f1->sinc.tulip[right].spinor;

    if ( spl > spr ){
        tClear(f1, left);
        tEqua(f1, left, 0, right , 0);
        return 1;
    }else if ( spl == spr ){
        tEqua(f1, left, 0, right , 0);
        if ( spl == cmpl )
            tEqua(f1, left, 1, right , 1);
        return 0;
    }else {
        tEqua(f1, left, 0, right , 0);
        printf("cmpl!\n");
      //  exit(0);
    }
    return 0;
    
}



INT_TYPE tAddTwo( struct field * f1 , enum division left , enum division right){
    enum spinType spl,spr;
    spl = f1->sinc.tulip[left].spinor;
    spr = f1->sinc.tulip[right].spinor;
    
    if ( spl > spr ){
        tClear(f1, left);
        tAddTw(f1, left, 0, right , 0);
        return 1;
    }else if ( spl == spr ){
        tAddTw(f1, left, 0, right , 0);
        if ( spl == cmpl )
            tAddTw(f1, left, 1, right , 1);
        return 0;
    }else {
        tAddTw(f1, left, 0, right , 0);

        printf("cmpl!\n");
       // exit(0);
    }
    return 0;
}

INT_TYPE tScaleOne( struct field * f1, enum division label,INT_TYPE spin, double scalar ){
    
    
    
    
    INT_TYPE L1 = CanonicalRank(f1,label,spin);
    
    double scale = fabs(scalar);
    
    
    INT_TYPE M2[3];
    length(f1,label,M2);
    
    if ( SPACE == 3 ){
        scale = pow ( scale , 1./3.);
        if ( scalar < 0 )
            scale *= -1;
        cblas_dscal(L1*M2[0],scale, streams(f1,label,spin,0),1);
        cblas_dscal(L1*M2[1],scale, streams(f1,label,spin,1),1);
        cblas_dscal(L1*M2[2],scale, streams(f1,label,spin,2),1);
    }else if ( SPACE == 2 ){
        scale = pow ( scale , 1./2.);
        cblas_dscal(L1*M2[1],scale, streams(f1,label,spin,1),1);
        if ( scalar < 0 )
            scale *= -1;
        cblas_dscal(L1*M2[0],scale, streams(f1,label,spin,0),1);
    }else if ( SPACE == 1 ){
        cblas_dscal(L1*M2[0],scalar, streams(f1,label,spin,0),1);
    }

    
    return 0;
}


INT_TYPE tScale( struct field * f1, enum division label, double scalar ){
    INT_TYPE sp;
    for ( sp = 0; sp < spins(f1, label) ; sp++)
        tScaleOne(f1, label, sp,scalar);
    return 0;
}

INT_TYPE tAddTw( struct field* f1 , enum division left, INT_TYPE lspin,  enum division right , INT_TYPE rspin){
    return xAddTw(f1,left, lspin, f1, right, rspin);
}

INT_TYPE xAddTw( struct field* f1 , enum division left, INT_TYPE lspin,struct field* f2 ,  enum division right , INT_TYPE rspin){
    if ( CanonicalRank(f2,right,rspin) ){
        INT_TYPE LL = f1->sinc.tulip[left].Current[lspin];
        INT_TYPE LR = f2->sinc.tulip[right].Current[rspin];
        //printf("++ %lld %lld\n", LL,LR);
        INT_TYPE MM = LL+LR;
        if ( MM > f1->sinc.tulip[left].Partition ){
            //tGetType allocation!
            printf("%lld >  %lld\n",MM, f1->sinc.tulip[left].Partition );
            printf("tAdd more money!\n");
            printf("%d %d\n", left, right);
            return 1;
        }
        INT_TYPE M2[SPACE],space;
        length(f1, left, M2);
        //printf("M2 %d s%d %d\n", M2[0],lspin,rspin);
        for ( space = 0; space < SPACE; space++){
            if ( ( memory(f2,right) == oneObject )  && space ){
                continue;
            }
            
            cblas_dcopy(LR*M2[space], streams(f2,right,rspin,space),1,streams(f1,left,lspin,space)+LL*M2[space],1);
        }
        f1->sinc.tulip[left].Current[lspin] += LR;
    }
    return 0;
}

INT_TYPE tSumMatrices(struct field *f1, enum division sum , enum division mat  ){
    
    INT_TYPE n2[3];
    length(f1, sum,n2);
    if ( bodies ( f1 , sum )== one ){
        
        if (bodies(f1, mat ) == one ){
            tAddTw(f1, sum ,0, mat,0 );
        }
        else {
            printf("eh?");
            //exit(0);
        }
    }else if ( bodies ( f1 , sum )== two && bodyType( f1, sum ) == electron){
        if (bodies(f1, mat ) == two ){
            tAddTwo(f1, sum , mat );
        }else if ( bodies ( f1, mat ) == one ){
            INT_TYPE I1,I2, I3, I4,body,r,space;
            INT_TYPE *n1 = vectorLen(f1, mat);
            double value;
            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
                for ( body = 0 ; body < 2 ; body++){
                    
                    for ( space = 0; space < SPACE ; space++){
                        Stream_Type * stream = streams(f1,sum,0,space)+n2[space]*CanonicalRank(f1, sum, 0);
                        if ( CanonicalRank(f1, sum, 0) > part(f1, sum )){
                            printf("part sum\n");
                            exit(0);
                        }
                        
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)//body 0
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)//body 0
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)//body 1
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)//body 1
                                    {
                                        if ( body == 0 ){
                                            value  = streams(f1,mat,0,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4);
                                        }else {
                                            value  = streams(f1,mat,0,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2);
                                        }
                                        stream[ (I1+I3*n1[space])+ ( I2+I4*n1[space])*n1[space]*n1[space] ] = value;
                                    }
                    }
                    f1->sinc.tulip[sum].Current[0]++;
                    
                    //  tAddTwo(f1, sum , quadCube);
                }
        }else {
            printf("hey!\n");
            exit(0);
        }
    }else if ( bodies ( f1 , sum )== two  && bodyType( f1, sum ) == h2plus){
        if (bodies(f1, mat ) == two ){
            tAddTwo(f1, sum , mat );
        }else if ( bodies ( f1, mat ) == one ){
            INT_TYPE I1,I2, I3, I4,body,r,space;
            INT_TYPE *n1 = vectorLen(f1, mat);
            double value;
            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
            {
                body = 0;
                if ( bodyType ( f1, mat ) == electron)
                    body = 1;

                
                for ( space = 0; space < SPACE ; space++){
                        Stream_Type * stream = streams(f1,sum,0,space)+n2[space]*CanonicalRank(f1, sum, 0);
                        if ( CanonicalRank(f1, sum, 0) > part(f1, sum )){
                            printf("part sum\n");
                            exit(0);
                        }
                        
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)//body 0
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)//body 0
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)//body 1
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)//body 1
                                    {
                                        if ( body == 0 ){
                                            value  = streams(f1,mat,0,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4);
                                        }else {
                                            value  = streams(f1,mat,0,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2);
                                        }
                                        stream[ (I1+I3*n1[space])+ ( I2+I4*n1[space])*n1[space]*n1[space] ] = value;
                                    }
                    }
                    f1->sinc.tulip[sum].Current[0]++;
                    
                    //  tAddTwo(f1, sum , quadCube);
                }
        }else {
            printf("hey!\n");
            exit(0);
        }
    }
    else if ( bodies ( f1, sum ) == three && bodyType(f1, sum ) == electron){
        
        if (bodies(f1, mat ) == three ){
            tAddTwo(f1, sum , mat );
        }else if ( bodies ( f1, mat ) == one ){
            INT_TYPE I1,I2, I3, I4,I5,I6,body,r,space;
            INT_TYPE *n1 = vectorLen(f1, mat);
            double value;
            
            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
                for ( body = 0 ; body < 3 ; body++){
                    for ( space = 0; space < SPACE ; space++){
                        Stream_Type * stream = streams(f1,sum,0,space)+n2[space]*CanonicalRank(f1, sum, 0);
                        if ( CanonicalRank(f1, sum, 0) > part(f1, sum )){
                            printf("part sum\n");
                            exit(0);
                        }
                        
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                                
                                            {
                                                value = 0;
                                                if ( body == 0 ){
                                                    value  = streams(f1,mat,0,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6);
                                                }else if ( body == 1 ) {
                                                    value  = streams(f1,mat,0,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6);
                                                }else if ( body == 2 ) {
                                                    value  = streams(f1,mat,0,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4);
                                                }
                                                stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]] = value;
                                            }
                    }
                    f1->sinc.tulip[sum].Current[0]++;
                    
                }
        }else if ( bodies ( f1, mat ) == two ){
            INT_TYPE I1,I2, I3, I4,I5,I6,pair,r,space,ve;
            INT_TYPE n1[3];
            n1[0] = f1->sinc.N1;
            n1[1] = f1->sinc.N1;
            n1[2] = f1->sinc.N1;

            double value;
            
            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
                for ( pair = 0 ; pair < 3 ; pair++){
                    for ( space = 0; space < SPACE ; space++){
                        Stream_Type * stream = streams(f1,sum,0,space)+n2[space]*CanonicalRank(f1, sum, 0);
                        if ( CanonicalRank(f1, sum, 0) > part(f1, sum )){
                            printf("part sum\n");
                            exit(0);
                        }
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                                
                                            {
                                                value = 0;
                                                if ( pair == 0 ){
                                                    value  = streams(f1,mat,0,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6);
                                                }else if ( pair == 1 ) {
                                                    value  = streams(f1,mat,0,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4);
                                                }else if ( pair == 2 ) {
                                                    value  = streams(f1,mat,0,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2);
                                                }
                                                ve = (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space];
                                                //                                                printf("%f %lld %lld\n", value,ve,n2[space]);
                                                //                                                fflush(stdout);
                                                stream[ ve ] = value;
                                                //                                                printf("x");
                                                //                                                fflush(stdout);
                                                
                                            }
                    }
                    f1->sinc.tulip[sum].Current[0]++;
                }
        }
        else {
            printf("Yo!");
            exit(0);
        }
        
    }    else if ( bodies ( f1, sum ) == three && bodyType(f1, sum ) == h2){
        printf("summatrices\n");
        exit(0);
    }else if (bodies(f1,sum) == four){
        
        if ( bodies ( f1, mat ) == one ){
            INT_TYPE I1,I2, I3, I4,I5,I6,I7,I8,body,r,space;
            INT_TYPE *n1 = vectorLen(f1, mat);
            double value;
            
            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
                for ( body = 0 ; body < 4 ; body++){
                    for ( space = 0; space < SPACE ; space++){
                        Stream_Type * stream = streams(f1,sum,0,space)+n2[space]*CanonicalRank(f1, sum, 0);
                        if ( CanonicalRank(f1, sum, 0) > part(f1, sum )){
                            printf("part sum\n");
                            exit(0);
                        }
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                                for ( I7 = 0 ; I7 < n1[space] ; I7++)
                                                    for ( I8 = 0 ; I8 < n1[space] ; I8++)
                                                        
                                                    {
                                                        value = 0;
                                                        if ( body == 0 ){
                                                            value  = streams(f1,mat,0,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6)*delta(I7-I8);
                                                        }else if ( body == 1 ) {
                                                            value  = streams(f1,mat,0,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6)*delta(I7-I8);
                                                        }else if ( body == 2 ) {
                                                            value  = streams(f1,mat,0,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I7-I8);
                                                        }else if ( body == 3 ) {
                                                            value  = streams(f1,mat,0,space)[ I7*n1[space]+I8 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I5-I6);
                                                        }
                                                        stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
                                                    }
                    }
                    f1->sinc.tulip[sum].Current[0]++;
                    
                }
        }else if ( bodies ( f1, mat ) == two ){
            INT_TYPE I1,I2, I3, I4,I5,I6,I7,I8,pair,r,space;
            INT_TYPE n1[3];
            n1[0] = f1->sinc.N1;
            n1[1] = f1->sinc.N1;
            n1[2] = f1->sinc.N1;
            double value;
            
            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
                for ( pair = 0 ; pair < 6 ; pair++){
                    for ( space = 0; space < SPACE ; space++){
                        Stream_Type * stream = streams(f1,sum,0,space)+n2[space]*CanonicalRank(f1, sum, 0);
                        if ( CanonicalRank(f1, sum, 0) > part(f1, sum )){
                            printf("part sum\n");
                            exit(0);
                        }
                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
                                                for ( I7 = 0 ; I7 < n1[space] ; I7++)
                                                    for ( I8 = 0 ; I8 < n1[space] ; I8++)
                                                        
                                                    {
                                                        //                                            0    e12,     1,3     2,4
                                                        //                                            1    e13,     1,5     2,6
                                                        //                                            2    e23,     3,5     4,6
                                                        //                                            3    e14,     1,7     2,8
                                                        //                                            4    e24,     3,7     4,8
                                                        //                                            5    e34      5,7     6,8
                                                        
                                                        
                                                        if ( pair == 0 ){
                                                            value  = streams(f1,mat,0,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6)*delta(I7-I8);
                                                        }else if ( pair == 1 ) {
                                                            value  = streams(f1,mat,0,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I7-I8);
                                                        }else if ( pair == 2 ) {
                                                            value  = streams(f1,mat,0,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I7-I8);
                                                        }else if ( pair == 3 ) {
                                                            value  = streams(f1,mat,0,space)[ (I1*n1[space]+I7) + (I2*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I5-I6);
                                                        }else if ( pair == 4 ) {
                                                            value  = streams(f1,mat,0,space)[ (I3*n1[space]+I7) + (I4*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I5-I6);
                                                        }else if ( pair == 5 ) {
                                                            value  = streams(f1,mat,0,space)[ (I5*n1[space]+I7) + (I6*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I3-I4);
                                                        }else {
                                                            printf ("rails!\n");
                                                            exit(0);
                                                        }
                                                        stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
                                                    }
                    }
                    f1->sinc.tulip[sum].Current[0]++;
                }
        }
        else {
            printf("Yo!");
            exit(0);
        }
        
    }
    return 0;
}


INT_TYPE tAlt(struct field * f1 , enum division label, INT_TYPE spin , INT_TYPE space1){
    
    INT_TYPE N1 = f1->sinc.N1;
    INT_TYPE N12 = (f1->sinc.N1-1)/2;
    INT_TYPE N2 = f1->sinc.N1*f1->sinc.N1;

    INT_TYPE I1,I2,space;
    INT_TYPE Current    =     f1->sinc.tulip[label].Current[spin]++;
    
    
    if ( f1->sinc.tulip[label].species == vector ){
        printf("reconsider...only one mode represented\n");
    }
    if ( f1->sinc.tulip[label].species == matrix ){
        for ( space = 0; space < SPACE ; space++){
            
            
            Stream_Type * stream = streams(f1,label,spin,space)+Current*N2;
            for ( I1 = 0 ; I1 < N1 ; I1++)
                for ( I2 = 0 ; I2 < N1 ; I2++){
                    stream[I1*N1+I2] =0.;
                }
            for ( I1 = 0 ; I1 < N1 ; I1++)
                for ( I2 = 0 ; I2 < N1 ; I2++){
                    {
                        if ( space == space1)
                            stream[I1*N1+I2] = sign(I2-I1);
                        else{
                            if ( I1 == I2 && I1 == N12 )
                                stream[I1*N1+I2] = 1;
                        }
                    }
                    
                }
        }
    }
    
    
    return 1;
}

INT_TYPE tEnd(struct field * f1 , enum division label, INT_TYPE spin , INT_TYPE space1){
    
    INT_TYPE N1 = f1->sinc.N1;
    INT_TYPE N2 = f1->sinc.N1*f1->sinc.N1;
    INT_TYPE I1,I2,space;
    INT_TYPE Current    =     f1->sinc.tulip[label].Current[spin]++;
    
    
    if ( f1->sinc.tulip[label].species == vector ){}
    if ( f1->sinc.tulip[label].species == matrix ){
        for ( space = 0; space < SPACE ; space++){
            
            
            Stream_Type * stream = streams(f1,label,spin,space)+Current*N2;
            for ( I1 = 0 ; I1 < N1 ; I1++)
                for ( I2 = 0 ; I2 < N1 ; I2++){
                    stream[I1*N1+I2] =0.;
                }
            for ( I1 = 0 ; I1 < N1 ; I1++)
                for ( I2 = 0 ; I2 < N1 ; I2++){
                    {
                        if ( space == space1 )
                        {
                            if ( I1 == 0  )
                                if ( I1 == I2 )
                                    stream[I1*N1+I2] = 1;
                            
                        }
                        else{
                            if ( I1 == I2)
                                stream[I1*N1+I2] = 1;
                        }
                    }
                    
                }
        }
    }
    
    
    return 1;
}

//INT_TYPE tPauli ( struct field * f1  ){
//    tClear(f1, PauliZ);
//    tClear(f1, PauliX);
//    if ( spins(f1,PauliZ) >= 2){
//        tId(f1, PauliZ,f1->sinc.arraySpin[f1->sinc.tulip[PauliZ].spinor][0][0]);
//        tScale(f1, PauliZ,-1);
//        tId(f1, PauliZ,f1->sinc.arraySpin[f1->sinc.tulip[PauliZ].spinor][1][1]);
//    }
//    if (spins(f1,PauliZ)  >= 3 ){
//        tId(f1, PauliX,f1->sinc.arraySpin[f1->sinc.tulip[PauliX].spinor][0][1] );
//    }
//    return 0;
//}

INT_TYPE tId ( struct field *f1 , enum division label,INT_TYPE spin ){
    
    INT_TYPE I1,I2,space;
    INT_TYPE Current ;
    {
        
        if ( f1->sinc.tulip[label].Current[spin] >= f1->sinc.tulip[label].Partition ){
            printf("tryed to add to full array\n");
            exit(0);
            return 0;
        }
        Current =  f1->sinc.tulip[label].Current[spin]++;
    }
    
    {
        if ( f1->sinc.tulip[label].species == vector ){
            INT_TYPE * B1;
            
            B1 = vectorLen(f1,label);
            for ( space = 0; space < SPACE ; space++){
                
                Stream_Type  * stream = streams(f1,label,spin,space)+Current*B1[space];
                for ( I2 = 0 ; I2 < B1[space] ; I2++){
                    stream[I2] = sign(2*I2%2-1);
                }
            }
        }
        
        else if  ( f1->sinc.tulip[label].species == matrix ) {
            INT_TYPE * B1;
            B1 = vectorLen(f1,label);
            
            
            
            for ( space = 0; space < SPACE ; space++){
                
                Stream_Type * stream = streams(f1,label,spin,space)+Current*B1[space]*B1[space];
                for ( I1 = 0 ; I1 < B1[space] ; I1++)
                    for ( I2 = 0 ; I2 < B1[space] ; I2++)
                        stream[I1*B1[space]+I2] =0.;
                
                for ( I1 = 0 ; I1 < B1[space] ; I1++)
                    {
                        
                            stream[I1*B1[space]+I1] = 1;
                        
                    }
            }
        }
        
        
    }
    return 0;
}


INT_TYPE tReplace( struct field *f1 , enum division label,INT_TYPE spin,INT_TYPE space,INT_TYPE l ){
    
    INT_TYPE I1,I2;
    {
        
        if ( f1->sinc.tulip[label].Current[spin] < l )
            return 0;
       }
    {
        if ( f1->sinc.tulip[label].species == vector ){
            INT_TYPE * B1;
            
            B1 = vectorLen(f1,label);
            {
                Stream_Type  * stream = streams(f1,label,spin,space)+l*B1[space];
                for ( I2 = 0 ; I2 < B1[space] ; I2++){
                    stream[I2] = sign(2*I2%2-1);
                }
            }
        }
        
        //        else if ( ( f1->sinc.tulip[label].species == matrix && bodies(f1,label)== two) ||  (f1->sinc.tulip[label].species == quartic && bodies(f1,label) == one)){
        //            for ( space = 0; space < SPACE ; space++){
        //                INT_TYPE * B1;
        //                B1 = vectorLen(f1,label);
        //
        //
        //                Stream_Type * stream = streams(f1,label,spin,space)+Current*B1[space]*B1[space]*B1[space]*B1[space];
        //                for ( I1 = 0 ; I1 < B1[space]*B1[space] ; I1++)
        //                    for ( I2 = 0 ; I2 < B1[space]*B1[space] ; I2++)
        //                        stream[I1*B1[space]*B1[space]+I2] =0.;
        //
        //                for ( I1 = 0 ; I1 < B1[space]*B1[space] ; I1++)
        //                    for ( I2 = 0 ; I2 < B1[space]*B1[space] ; I2++)
        //                    {
        //
        //                        if ( I1==I2  )
        //                            stream[I1*B1[space]*B1[space]+I2] = 1;
        //
        //                    }
        //            }
        //        }
        else if  ( f1->sinc.tulip[label].species == matrix ) {
            INT_TYPE * B1;
            B1 = vectorLen(f1,label);
            
            
            

            {
                Stream_Type * stream = streams(f1,label,spin,space)+l*B1[space]*B1[space];
                for ( I1 = 0 ; I1 < B1[space] ; I1++)
                    for ( I2 = 0 ; I2 < B1[space] ; I2++)
                        stream[I1*B1[space]+I2] =0.;
                
                for ( I1 = 0 ; I1 < B1[space] ; I1++)
                    for ( I2 = 0 ; I2 < B1[space] ; I2++)
                    {
                        
                        if ( I1==I2  )
                            stream[I1*B1[space]+I2] = 1;
                        
                    }
            }
        }
        
        
    }
    return 0;
}

INT_TYPE tBoot ( struct field *f1 , enum division label,INT_TYPE spin ){
    
    INT_TYPE I1,I2,space;
    INT_TYPE Current ;
    {
        
        if ( f1->sinc.tulip[label].Current[spin] >= f1->sinc.tulip[label].Partition )
            return 0;
        Current =  f1->sinc.tulip[label].Current[spin]++;
    }
    {
        if ( f1->sinc.tulip[label].species == vector && bodies(f1,label) == two){
            INT_TYPE B1[3];
            B1[0] = f1->sinc.N1;
            B1[1] = f1->sinc.N1;
            B1[2] = f1->sinc.N1;
            for ( space = 0; space < SPACE ; space++){
                
                Stream_Type  * stream = streams(f1,label,spin,space)+Current*B1[space];
                for ( I1 = 0 ; I1< B1[space] ; I1++)
                    for ( I2 = 0 ; I2 < B1[space] ; I2++){
                        stream[I1*B1[space]+I2] = exp(-abs(I1-(B1[space]-1)/2)*0.01)*exp(-sqr(I2-(B1[space]-1)/2)*0.01);
                }
            }
        }
        else if ( f1->sinc.tulip[label].species == vector && bodies(f1,label) == one){
            INT_TYPE B1[3];
            B1[0] = f1->sinc.N1;
            B1[1] = f1->sinc.N1;
            B1[2] = f1->sinc.N1;
            for ( space = 0; space < SPACE ; space++){
                
                Stream_Type  * stream = streams(f1,label,spin,space)+Current*B1[space];
                    for ( I2 = 0 ; I2 < B1[space] ; I2++){
                        stream[I2] = exp(-sqr(I2-(B1[space]-1)/2)*0.01);
                    }
            }
        }

      
        
    }
    return 0;
}


double vectorElement (struct field * f1, enum division state, INT_TYPE l1,INT_TYPE l2 , INT_TYPE l3 ){
    INT_TYPE spin ;
    double den = 0.;
    double maxDen = 0.;
    INT_TYPE l;
    INT_TYPE * n1 = vectorLen(f1, state);
    for ( spin = 0; spin < spins(f1, state ) ;spin++){
        for ( l = 0; l < CanonicalRank(f1, state, spin ) ; l++)
            den += (streams(f1, state, spin , 0 )[l*n1[0]+l1]* streams(f1, state, spin , 1 )[l*n1[1]+l2]* streams(f1, state, spin , 2 )[l*n1[2]+l3]);
        if ( den > maxDen)
            maxDen = den;
    }
    
    return den;
}

double matrixElement (struct field * f1, enum division label, INT_TYPE i , INT_TYPE i2, INT_TYPE j,INT_TYPE j2, INT_TYPE k , INT_TYPE k2 ){
    INT_TYPE *n1 = f1->sinc.Basis;
    INT_TYPE *N2  = f1->sinc.Basis2;
    
    
    INT_TYPE l;
    INT_TYPE spin = 0;;
    double product,sum = 0;
    for ( l = 0; l < CanonicalRank(f1,label,0); l++){
        product = 1;
        
        product = streams ( f1, label, spin, 0 )[N2[0]*l + n1[0]*(i)+(i2)];
        product *= streams ( f1, label, spin, 1 )[N2[0]*l + n1[1]*(j)+(j2)];
        product *= streams ( f1, label, spin, 2 )[N2[0]*l + n1[2]*(k)+(k2)];
        //   product *= cc ( f1, label, spin, 0 )[N2*l];
        
        sum += product ;//* coeff(f1, label,spin, l);
    }
    return sum;
}

void nuclearArray (struct field * f1,  enum division array,INT_TYPE M1){
    INT_TYPE *n1 = f1->sinc.Basis;
    INT_TYPE N1[3];
    INT_TYPE l = 2;
    length(f1, oneVector, N1);
    INT_TYPE i,k,j,at,mi,mj,mk;
    double * array3d = myStreams(f1, array, 0),x,y,z,d,mini,dis;

    for ( i = 0; i < M1*M1*M1 ; i++)
        array3d[i] = 0.;


    for ( at = 1 ; at <= f1->Na ; at++){
        mini = M1;
        mi = 0;
        mj = 0;
        mk = 0;
        d = 1./ (double)(M1)  * (double)(n1[0]) * f1->sinc.d;
        for ( i = 0; i < M1 ; i++){
            for ( j = 0; j < M1 ; j++){
                for ( k = 0 ; k < M1 ; k++){
                    x = ( i - (M1-1)/2 ) * d ;
                    y = ( j - (M1-1)/2 ) * d;
                    z = ( k - (M1-1)/2 ) * d;

                    dis = sqr( x - f1->atoms[at].position[1] ) + sqr( y - f1->atoms[at].position[2] ) + sqr( z - f1->atoms[at].position[3] );
                    if (  dis < mini ){
                        mini = dis;
                        mi= i;
                        mj = j;
                        mk = k;
                    }
                }
            }
        }
        for ( i = imax(0,mi-l); i < imin(M1,mi+l); i++)
            for ( j = imax(0,mj-l); j < imin(M1,mj+l); j++)
                for ( k = imax(0,mk-l); k < imin(M1,mk+l); k++)
                    if ( sqr(mi-i) + sqr( mj-j)+ sqr(k-mk) < l)
                        array3d[i+j*M1+k*M1*M1] += 1.;
        
    }



}

void vectorArray (struct field * f1, enum division oneVector,  enum division array,INT_TYPE M1){
    INT_TYPE *n1 = f1->sinc.Basis;
    INT_TYPE N1[3],stride[3];
    length(f1, oneVector, N1);
    INT_TYPE i;
    double * array3d = myStreams(f1, array, 0);
    INT_TYPE spin = 0.,space,r;;
    double * basis = myStreams(f1, oneBasis,0);
    double * oneDim = myStreams(f1, oneArray,0);
    double * tempArea = oneDim + 3 * M1;
    if ( part(f1, oneArray ) < 3*M1+M1*M1){
        printf("vectorArray\n");
        exit(0);
    }
    for ( i = 0; i < M1*M1*M1 ; i++)
        array3d[i] = 0.;
    if ( bodies (f1,oneVector) == one && species(f1, oneVector ) == vector){
        stride[0] = n1[0];
        stride[1] = n1[1];
        stride[2] = n1[2];
    }
    else     if ( bodies (f1,oneVector) == one && species(f1, oneVector ) == matrix){
        stride[0] = n1[0]+1;
        stride[1] = n1[1]+1;
        stride[2] = n1[2]+1;
    }else{
        printf("vectorArray2\n");

        exit(0);
    }
    for ( r = 0;  r< CanonicalRank(f1, oneVector, 0); r++){
        
      // oneBasis = (i,ii)  ii in M1
        for ( space = 0; space < SPACE ; space++)
            cblas_dgemv(CblasColMajor, CblasNoTrans,M1,n1[space],1.,
                        basis,M1,streams(f1, oneVector, spin , space ) + r*N1[space] ,stride[space], 0., oneDim+M1*space,1  );
        for ( i = 0; i < M1*M1 ; i++)
            tempArea[i] = 0.;

        cblas_dger(CblasColMajor, M1,M1, 1. , oneDim,1, oneDim+M1,1, tempArea,M1);
        cblas_dger(CblasColMajor, M1*M1,M1, 1. , tempArea,1, oneDim+M1+M1,1, array3d,M1*M1);
    }

}

int double_cmp( const void * a ,const void * b ){
    const double* A = (const double* )a;
    const double* B = (const double* )b;
    if ( *A < * B )
        return 1;
    else if ( *A == *B )
        return 0;
    else return -1;
}


double levelDetermine ( INT_TYPE M1 , double * array , double level){
    double *temp = array+M1*M1*M1,sum,psum;
    INT_TYPE i;
    cblas_dcopy(M1*M1*M1, array, 1, temp, 1);
    qsort(temp, M1*M1*M1, sizeof(double), &double_cmp);
    sum = 0.;
    for ( i = 0 ;i < M1*M1*M1 ; i++)
        sum += temp[i];
    i= 0;
    psum = 0.;
    while (psum < level*sum ){
        psum += temp[i++];
    }
    
    return temp[i-1];
}


INT_TYPE assignCores(struct field * f1, INT_TYPE parallel ){
#ifdef OMP
    INT_TYPE nSlot = f1->mem1->rt->NSlot;
    //    INT_TYPE nTraffic = f1->mem1->rt->NTraffic;
    INT_TYPE nCore = f1->mem1->rt->NCore;

    INT_TYPE omp;
    if ( parallel == 0){
        omp = 1;
    }else if ( parallel /*== 1*/ ){
        omp = nCore;
    }else if ( parallel == 2 ){
    //    omp = nTraffic;
    }else {
        printf (" unknown parallel\n");
        exit(0);
    }
    omp_set_num_threads(omp);
#endif
#ifdef MKL
    mkl_set_num_threads(nSlot/omp);
#endif

    return 0;

}

void printVectorAllocations(struct field *f1){
    INT_TYPE vecG = 3*(f1->sinc.tulip[f1->sinc.end].Address -  f1->sinc.tulip[eigenVectors].Address);
    printf("vector Contribution \t G%f",vecG/1000000000./(sizeof(Stream_Type)));
}
struct basisElement grabBasis (struct field * f1, INT_TYPE component, INT_TYPE particle, INT_TYPE elementIndex){
    struct basisElement boa;
    if ( component == 0 ){
        boa.periodic = f1->mem1->rt->runFlag % 2;
    }  else  if ( component == 1 ){
        boa.periodic = (f1->mem1->rt->runFlag/2) % 2;
    }  else  if ( component == 2 ){
        boa.periodic = (f1->mem1->rt->runFlag/4) % 2;
    }
    boa.basis = SincBasisElement;
    INT_TYPE reducedElementIndex = elementIndex;
    boa.index = reducedElementIndex - (f1->sinc.dims[particle][component]-1)/2;
    boa.length = f1->sinc.d;
    boa.origin = 0.;
    
    boa.auxIndex = f1->sinc.dims[particle][component];
    boa.dim = component;
    // boa.body = particle;
    boa.association = 0;
    boa.auxLength = 1.;
    
    if ( f1->mem1->rt->bodyType == h2plus)
        boa.length /= 1.3;
        
        return boa;
}


#if 0

double xOneBand (INT_TYPE rank,struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic){
    double D = f1->sinc.d;
    INT_TYPE space,i,l,i2,l2,r,aPeriodic;
    INT_TYPE N1 = f2->sinc.N1;
    INT_TYPE N12 = (N1-1)/2;
    INT_TYPE L1 =  f1->sinc.N1;
    INT_TYPE L12 = (L1-1)/2;
    f2->sinc.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);
    
    //    for ( r= 0 ; r < 3 ;r++){
    //    for ( l= 0 ; l < 9 ; l++)
    //        printf("%f:", streams(f1, vector1,s1,0)[l]);
    //        printf("\n");
    //    }
    //        printf("%d %f : %d %f\n", L1, D, N1, f2->sinc.d);
    //
    
    for ( i = 0; i < N1 ; i++)
    {
            
            //build
            for ( space = 0;space < SPACE; space++){
                if ( space == 0 )
                    aPeriodic = periodic%2;
                else if ( space == 1 )
                    aPeriodic = (periodic/2)%2;
                else
                    aPeriodic = (periodic/4)%2;
                
                
                for ( l = 0 ; l < L1 ; l++)
                {
                    if ( aPeriodic ){
                            streams(f1, foundationStructure,rank, space )[l] =
                        periodicSS ( D,D*(l-L12),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i-N12 ) ,f2->sinc.N1);
                            //                            periodicSS ( D,D*(l),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i ) ,f2->sinc.N1)*
                            //                            periodicSS ( D,D*(l2),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i2 ),f2->sinc.N1 );
                            
                        }else {
                            streams(f1, foundationStructure,rank, space )[l] =
                            SS ( D,D*(l-L12), f2->sinc.d, f2->sinc.d*( i-N12 ) );
                            
                        }
                    }
            }
            for ( space = 0;space < SPACE; space++){
                
                for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                    streams(f2, out, s2,space)[r*N1 + (i)] = cblas_ddot(L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1,1);
                    //      printf("%1.3f:", streams(f2, out, s2,space)[r*N1*N1 + (i+i2*N1)]);
                }
                //   printf("\n");
                //
            }
            //  printf("\n");
            
            
            
        }
    f2->sinc.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    
    
    
    return 0.;
}

double xTwoBand (INT_TYPE rank,struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic){
    double D = f1->sinc.d;
    INT_TYPE space,i,l,i2,l2,r,aPeriodic;
    INT_TYPE N1 = f2->sinc.N1;
    INT_TYPE N12 = (N1-1)/2;
    INT_TYPE L1 =  f1->sinc.N1;
    INT_TYPE L12 = (L1-1)/2;
    f2->sinc.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);
    
//    for ( r= 0 ; r < 3 ;r++){
//    for ( l= 0 ; l < 9 ; l++)
//        printf("%f:", streams(f1, vector1,s1,0)[l]);
//        printf("\n");
//    }
//        printf("%d %f : %d %f\n", L1, D, N1, f2->sinc.d);
//
    
    for ( i = 0; i < N1 ; i++)
        for ( i2 = 0; i2 < N1 ; i2++){
            
            
            //build
            for ( space = 0;space < SPACE; space++){
                if ( space == 0 )
                    aPeriodic = periodic%2;
                else if ( space == 1 )
                    aPeriodic = (periodic/2)%2;
                else
                    aPeriodic = (periodic/4)%2;
                
                
                for ( l = 0 ; l < L1 ; l++)
                    for ( l2 = 0 ; l2 < L1 ; l2++){
                        if ( aPeriodic ){
                            streams(f1, foundationStructure,rank, space )[l2*L1+l] =
                            periodicSS ( D,D*(l-L12),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i-N12 ) ,f2->sinc.N1)*
                            periodicSS ( D,D*(l2-L12),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i2-N12 ),f2->sinc.N1 );
//                            periodicSS ( D,D*(l),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i ) ,f2->sinc.N1)*
//                            periodicSS ( D,D*(l2),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i2 ),f2->sinc.N1 );

                        }else {
                            streams(f1, foundationStructure,rank, space )[l2*L1+l] =
                            SS ( D,D*(l-L12), f2->sinc.d, f2->sinc.d*( i-N12 ) )*
                            SS ( D,D*(l2-L12), f2->sinc.d, f2->sinc.d*( i2-N12 ) );
                        }
                    }
            }
            for ( space = 0;space < SPACE; space++){
                
                for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                    streams(f2, out, s2,space)[r*N1*N1 + (i+i2*N1)] = cblas_ddot(L1*L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1*L1,1);
              //      printf("%1.3f:", streams(f2, out, s2,space)[r*N1*N1 + (i+i2*N1)]);
                }
             //   printf("\n");
                //
            }
          //  printf("\n");

            
            
        }
    f2->sinc.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    
    
    
    return 0.;
}


double xThreeBand (INT_TYPE rank,struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic){
    double D = f1->sinc.d;
    INT_TYPE space,i,l,i2,i3,l2,l3,r,aPeriodic;
    INT_TYPE N1 = f2->sinc.N1;
    INT_TYPE N12 = (N1-1)/2;
    INT_TYPE L1 =  f1->sinc.N1;
    INT_TYPE L12 = (L1-1)/2;
    f2->sinc.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);
    for ( i = 0; i < N1 ; i++)
        for ( i2 = 0; i2 < N1 ; i2++)
            for ( i3 = 0; i3 < N1 ; i3++)
            {
                
                
                //build
                for ( space = 0;space < SPACE; space++){
                    if ( space == 0 )
                        aPeriodic = periodic%2;
                    else if ( space == 1 )
                        aPeriodic = (periodic/2)%2;
                    else
                        aPeriodic = (periodic/4)%2;

                    for ( l = 0 ; l < L1 ; l++)
                        for ( l2 = 0 ; l2 < L1 ; l2++)
                            for ( l3 = 0 ; l3 < L1 ; l3++)
                                if ( aPeriodic ){
                                    streams(f1, foundationStructure,rank, space )[(l3*L1*L1+l2*L1+l)] =
                                    periodicSS ( D,D*(l-L12),f1->sinc.N1,f2->sinc.d, f2->sinc.d*( i-N12 ) ,f2->sinc.N1)*
                                    periodicSS ( D,D*(l2-L12),f1->sinc.N1,f2->sinc.d, f2->sinc.d*( i2-N12 ),f2->sinc.N1 )*
                                    periodicSS ( D,D*(l3-L12),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i3-N12 ) ,f2->sinc.N1);
                                }else {
                                    streams(f1, foundationStructure,rank, space )[(l3*L1*L1+l2*L1+l)] =
                                    SS ( D,D*(l-L12), f2->sinc.d, f2->sinc.d*( i-N12 ) )*
                                    SS ( D,D*(l2-L12), f2->sinc.d, f2->sinc.d*( i2-N12 ) )*
                                    SS ( D,D*(l3-L12), f2->sinc.d, f2->sinc.d*( i3-N12 ) );
                                }
                }
                
                for ( space = 0;space < SPACE; space++){
                    
                    for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++)
                        streams(f2, out, s2,space)[r*N1*N1*N1 + (i3*N1*N1+i2*N1+i)] = cblas_ddot(L1*L1*L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1*L1*L1,1);
                    //
                }
                
                
                
            }
    f2->sinc.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    return 0.;
}

double xFourBand (INT_TYPE rank,struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic){
    double D = f1->sinc.d;
    INT_TYPE space,i,l,i2,i3,i4,l2,l3,l4,r,aPeriodic;
    INT_TYPE N1 = f2->sinc.N1;
    INT_TYPE N12 = (N1-1)/2;
    INT_TYPE L1 =  f1->sinc.N1;
    INT_TYPE L12 = (L1-1)/2;
    f2->sinc.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);
    for ( i = 0; i < N1 ; i++)
        for ( i2 = 0; i2 < N1 ; i2++)
            for ( i3 = 0; i3 < N1 ; i3++)
                for ( i4 = 0 ; i4 < N1; i4++)
            {
                
                
                //build
                for ( space = 0;space < SPACE; space++){
                    if ( space == 0 )
                        aPeriodic = periodic%2;
                    else if ( space == 1 )
                        aPeriodic = (periodic/2)%2;
                    else
                        aPeriodic = (periodic/4)%2;

                    for ( l = 0 ; l < L1 ; l++)
                        for ( l2 = 0 ; l2 < L1 ; l2++)
                            for ( l3 = 0 ; l3 < L1 ; l3++)
                                for ( l4 = 0; l4 < L1 ; l4++)
                                if ( aPeriodic ){
                                    streams(f1, foundationStructure,rank, space )[(l4*L1*L1*L1+l3*L1*L1+l2*L1+l)] =
                                    periodicSS ( D,D*(l-L12),f1->sinc.N1,f2->sinc.d, f2->sinc.d*( i-N12 ) ,f2->sinc.N1)*
                                    periodicSS ( D,D*(l2-L12),f1->sinc.N1,f2->sinc.d, f2->sinc.d*( i2-N12 ),f2->sinc.N1 )*
                                    periodicSS ( D,D*(l3-L12),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i3-N12 ) ,f2->sinc.N1)*
                                    periodicSS ( D,D*(l4-L12),f1->sinc.N1, f2->sinc.d, f2->sinc.d*( i4-N12 ) ,f2->sinc.N1);
                                }else {
                                    streams(f1, foundationStructure,rank, space )[(l4*L1*L1*L1+l3*L1*L1+l2*L1+l)] =
                                    SS ( D,D*(l-L12), f2->sinc.d, f2->sinc.d*( i-N12 ) )*
                                    SS ( D,D*(l2-L12), f2->sinc.d, f2->sinc.d*( i2-N12 ) )*
                                    SS ( D,D*(l3-L12), f2->sinc.d, f2->sinc.d*( i3-N12 ) )*
                                    SS ( D,D*(l4-L12), f2->sinc.d, f2->sinc.d*( i4-N12 ) );
                                }
                }
                
                for ( space = 0;space < SPACE; space++){
                    
                    for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++)
                        streams(f2, out, s2,space)[r*N1*N1*N1*N1 + (i4*N1*N1*N1+i3*N1*N1+i2*N1+i)] = cblas_ddot(L1*L1*L1*L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1*L1*L1*L1,1);
                    //
                }
                
                
                
            }
        f2->sinc.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    return 0.;
}

#elif 1




//struct basisElement {
//    enum basisElementType basis;
//    INT_TYPE index;
//    double length;
//    double origin;
//
//    INT_TYPE auxIndex; //for periodic Sincs
//    double auxLength;
//    INT_TYPE association;
//    INT_TYPE dim;
//    INT_TYPE body;
//};
//

double xOneBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2,INT_TYPE oldPeriodic){
    INT_TYPE si,space,i,l,r,rank;
    INT_TYPE N1;
    INT_TYPE L1;
    f2->sinc.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);

    for ( space = 0;space < SPACE; space++){
        N1 = f2->sinc.N1;
        L1 = f1->sinc.N1;
        struct basisElement gl[1][L1];
        struct basisElement gi[1][N1];

        for ( i = 0; i < N1 ; i++)
            gi[0][i] = grabBasis(f2, space, 1, i);

        for ( l = 0 ; l < L1 ; l++)
            gl[0][l] = grabBasis(f1, space, 1, l);
        
//#ifdef OMP
//#pragma omp parallel for private (si,i,rank,l,r) schedule(dynamic,1)
//#endif
        for ( si = 0; si < N1 ; si++){
//#ifdef OMP
//            rank = omp_get_thread_num();
//#else
           rank = 0;
//#endif

       // for ( i = 0; i < N1 ; i++)
            i = si;
            
            //build
            
            for ( l = 0 ; l < L1 ; l++)
            {
                streams(f1, foundationStructure,rank, space )[l] =
                BoB (gl[0][l],gi[0][i] );
            }
            
            
            for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                streams(f2, out, s2,space)[r*N1 + (i)] = cblas_ddot(L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1,1);
            }
        }
    }
    f2->sinc.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    return 0.;
}





double xTwoBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2,INT_TYPE oldPeriodic){
    INT_TYPE space,i,l,i2,l2,r,p,rank;
    INT_TYPE N1,si ;
    INT_TYPE L1 ;
    f2->sinc.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);
    for ( space = 0;space < SPACE; space++){
        N1 = f2->sinc.N1;
        L1 = f1->sinc.N1;
        
        struct basisElement gl[2][L1];
        struct basisElement gi[2][N1];
        
        for ( i = 0; i < N1 ; i++){
            gi[0][i] = grabBasis(f2, space, 1, i);
            gi[1][i] = grabBasis(f2, space, 2, i);
        }
        for ( l = 0 ; l < L1 ; l++){
            gl[0][l] = grabBasis(f1, space, 1, l);
            gl[1][l] = grabBasis(f1, space, 2, l);
        }
        INT_TYPE pt = 2;
        double gx[pt][L1][N1];
        for ( p = 0; p < pt ; p++)
            for ( i = 0; i < N1 ; i++)
                for ( l = 0 ; l < L1 ; l++){
                    gx[p][l][i] = BoB (grabBasis(f1, space, 1+p, l),grabBasis(f2, space, 1+p, i) );
//                  /  printf("%d %d %d %f\n", p,i,l,gx[p][l][i]);
                }
        
                    
                    
        
//#ifdef OMP
//#pragma omp parallel for private (si,i,i2,rank,l,l2,r) schedule(dynamic,1)
//#endif
        for ( si = 0; si < N1*N1 ; si++){
//#ifdef OMP
//            rank = omp_get_thread_num();
//#else
            rank = 0;
//#endif
            i = si% N1;
            i2 = (si/N1)%N1;
        
              //  printf("%d %d %d \n", vector1, i,i2);
                //build
                for ( l = 0 ; l < L1 ; l++)
                    for ( l2 = 0 ; l2 < L1 ; l2++){
                        streams(f1, foundationStructure,rank, space )[l2*L1+l] =
                        gx[0][l][i]*
                        gx[1][l2][i2];
                        
                    }
                
                
                for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                    streams(f2, out, s2,space)[r*N1*N1 + (i+i2*N1)] = cblas_ddot(L1*L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1*L1,1);
                }
            }
    }
    f2->sinc.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    
    
    
    return 0.;
}


double xThreeBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2,INT_TYPE oldPeriodic){
    INT_TYPE si,space,i,l,i2,i3,l2,l3,r,rank;
    INT_TYPE N1 ;
    INT_TYPE L1 ;
    f2->sinc.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);
    
    for ( space = 0;space < SPACE; space++){
        N1 = f2->sinc.N1;
        L1 = f1->sinc.N1;
        struct basisElement gl[3][L1];
        struct basisElement gi[3][N1];
        
        for ( i = 0; i < N1 ; i++){
            gi[0][i] = grabBasis(f2, space, 1, i);
            gi[1][i] = grabBasis(f2, space, 2, i);
            gi[2][i] = grabBasis(f2, space, 3, i);

        }
        for ( l = 0 ; l < L1 ; l++){
            gl[0][l] = grabBasis(f1, space, 1, l);
            gl[1][l] = grabBasis(f1, space, 2, l);
            gl[2][l] = grabBasis(f1, space, 3, l);
        }

        INT_TYPE pt = 3,p;
        double gx[pt][L1][N1];
        for ( p = 0; p < pt ; p++)
            for ( i = 0; i < N1 ; i++)
                for ( l = 0 ; l < L1 ; l++)
                    gx[p][l][i] = BoB (gl[p][l],gi[p][i] );
        
//#ifdef OMP
//#pragma omp parallel for private (si,i,i2,i3,rank,l,l2,l3,r) schedule(dynamic,1)
//#endif
        for ( si = 0; si < N1*N1*N1 ; si++){
//#ifdef OMP
//            rank = omp_get_thread_num();
//#else
            rank = 0;
//#endif

//        for ( i = 0; i < N1 ; i++)
//            for ( i2 = 0; i2 < N1 ; i2++)
//                for ( i3 = 0; i3 < N1 ; i3++)
            i = si% N1;
            i2 = (si/N1)%N1;
            i3 = (si/(N1*N1)) % N1;
                    
                    
                    //build
                    
                    for ( l = 0 ; l < L1 ; l++)
                        for ( l2 = 0 ; l2 < L1 ; l2++)
                            for ( l3 = 0 ; l3 < L1 ; l3++)
                                streams(f1, foundationStructure,rank, space )[(l3*L1*L1+l2*L1+l)] =
                                gx[0][l][i]*
                                gx[1][l2][i2]*
                                gx[2][l3][i3];
                    
                    
                    for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++)
                        streams(f2, out, s2,space)[r*N1*N1*N1 + (i3*N1*N1+i2*N1+i)] = cblas_ddot(L1*L1*L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1*L1*L1,1);
                    //
                    
                    
                    
                    
                }
    }
    f2->sinc.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    return 0.;
}

double xFourBand (struct field *f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2,INT_TYPE oldPeriodic){
    INT_TYPE si,space,i,l,i2,i3,i4,l2,l3,l4,r,rank;
    INT_TYPE N1 ;
    INT_TYPE L1 ;
    f2->sinc.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);
    
    
    for ( space = 0;space < SPACE; space++){
        N1 = f2->sinc.N1;
        L1 = f1->sinc.N1;
        struct basisElement gl[4][L1];
        struct basisElement gi[4][N1];
        
        for ( i = 0; i < N1 ; i++){
            gi[0][i] = grabBasis(f2, space, 1, i);
            gi[1][i] = grabBasis(f2, space, 2, i);
            gi[2][i] = grabBasis(f2, space, 3, i);
            gi[3][i] = grabBasis(f2, space, 4, i);

        }
        for ( l = 0 ; l < L1 ; l++){
            gl[0][l] = grabBasis(f1, space, 1, l);
            gl[1][l] = grabBasis(f1, space, 2, l);
            gl[2][l] = grabBasis(f1, space, 3, l);
            gl[3][l] = grabBasis(f1, space, 4, l);

        }
        INT_TYPE pt = 4,p;
        double gx[pt][L1][N1];
        for ( p = 0; p < pt ; p++)
            for ( i = 0; i < N1 ; i++)
                for ( l = 0 ; l < L1 ; l++)
                    gx[p][l][i] = BoB (gl[p][l],gi[p][i] );


//#ifdef OMP
//#pragma omp parallel for private (si,i,i2,i3,i4,rank,l,l2,l3,l4,r) schedule(dynamic,1)
//#endif
        for ( si = 0; si < N1*N1*N1*N1 ; si++){
//#ifdef OMP
//            rank = omp_get_thread_num();
//#else
            rank = 0;
//#endif
            //        for ( i = 0; i < N1 ; i++)
            //            for ( i2 = 0; i2 < N1 ; i2++)
            //                for ( i3 = 0; i3 < N1 ; i3++)
            i = si% N1;
            i2 = (si/N1)%N1;
            i3 = (si/(N1*N1)) % N1;
            i4 = (si/(N1*N1*N1)) % N1;

                        //build
                        
                        for ( l = 0 ; l < L1 ; l++)
                            for ( l2 = 0 ; l2 < L1 ; l2++)
                                for ( l3 = 0 ; l3 < L1 ; l3++)
                                    for ( l4 = 0; l4 < L1 ; l4++)
                                        streams(f1, foundationStructure,rank, space )[(l4*L1*L1*L1+l3*L1*L1+l2*L1+l)] =
                                        gx[0][l][i]*
                                        gx[1][l2][i2]*
                                        gx[2][l3][i3]*
                                        gx[3][l4][i4];

                        
                        
                        
                        for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++)
                            streams(f2, out, s2,space)[r*N1*N1*N1*N1 + (i4*N1*N1*N1+i3*N1*N1+i2*N1+i)] = cblas_ddot(L1*L1*L1*L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1*L1*L1*L1,1);
                        //
                    }
        
        
        
    }
    f2->sinc.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    return 0.;
}

#endif
