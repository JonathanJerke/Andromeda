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

INT_TYPE spaces( struct sinc_label f1, enum division label){
    if ( bodies(f1,label) == four)
        return SPACE;
    else
        return SPACE;
}


enum division name ( struct sinc_label f1, enum division label){
    return f1.tulip[f1.tulip[f1.tulip[label].name].name].name;
}


INT_TYPE part ( struct sinc_label f1 , enum division label){
    
    if ( label > f1.end || label < 0 ){
        printf("rank past end\n");
    }
    return f1.tulip[name(f1,label)].Partition;
}

INT_TYPE species ( struct sinc_label f1 , enum division label ){
    enum division u = name(f1,label);
    return f1.tulip[u].species;
}



enum bodyType Bodies ( struct sinc_label f1 , enum division label,INT_TYPE space ){
    if ( species ( f1, label ) == vector )
        return f1.rose[space].body ;
    else
        return f1.tulip[name(f1,label)].space[space].body;
        
}

enum bodyType bodies ( struct sinc_label f1 , enum division label ){
    enum bodyType x=nada,u;
    INT_TYPE space;
    for ( space = 0 ; space < SPACE ; space++)
    {
        u =Bodies(f1, label, space);
        if ( u > x )
            x = u ;
    }
    return x;
}

enum particleType particle ( struct sinc_label f1 , enum division label, INT_TYPE space ){
    return f1.rose[space].particle;
}


INT_TYPE header ( struct sinc_label f1 , enum division label ){
    return f1.tulip[name(f1,label)].header;
}

//void initPointerTensors(struct field *f1){
//    INT_TYPE sp;
//    f1.tulip[v2].name = vectors;
//    f1.tulip[v2].purpose = ptObject;
//    
//    f1.tulip[u2].name = vectors;
//    f1.tulip[u2].purpose = ptObject;
//    tClear(f1, v2);
//    tClear(f1,u2);
//    tClear(f1, u);
//    tClear(f1,v);
//    tClear(f1, p);
//    tClear(f1,pp);
//    tClear(f1,pointer);
//    
//    
//    f1.tulip[u].purpose = ptObject;
//    f1.tulip[v].purpose = ptObject;
//    f1.tulip[u].Partition = 1;
//    f1.tulip[v].Partition = 1;
//    f1.tulip[u].parallel = 1;
//    f1.tulip[v].parallel = 1;
//    
//    f1.tulip[u].species = matrix;
//    f1.tulip[v].species = matrix;
//    
//    f1.tulip[p].purpose = ptObject;
//    f1.tulip[pp].purpose = ptObject;
//    f1.tulip[q].purpose = ptObject;
//    f1.tulip[qq].purpose = ptObject;
//    
//    f1.tulip[pointer].purpose = ptObject;
//    f1.tulip[pointer].parallel = 1;
//    
//    for ( sp = 0 ; sp < spins(f1, diis) ; sp++){
//        f1.tulip[p].Current[sp] = f1.tulip[project].Partition;
//        f1.tulip[pp].Current[sp] = f1.tulip[project].Partition;
//        f1.tulip[q].Current[sp] = f1.tulip[project].Partition;
//        f1.tulip[qq].Current[sp] = f1.tulip[project].Partition;
//        f1.tulip[u].Current[sp] = 1;
//        f1.tulip[v].Current[sp] = 1;
//    }
//    
//    f1.tulip[p].species = matrix;
//    f1.tulip[pp].species = matrix;
//    f1.tulip[p].Partition = f1.tulip[project].Partition;
//    f1.tulip[pp].Partition = f1.tulip[project].Partition;
//    f1.tulip[q].Partition = f1.tulip[project].Partition;
//    f1.tulip[qq].Partition = f1.tulip[project].Partition;
//    f1.tulip[q].species = matrix;
//    f1.tulip[qq].species = matrix;
//    
//}

//INT_TYPE defineSpinors (struct field *f1 ){
//    INT_TYPE si,sp,sp2,s;
//
//    for ( s = 0; s < NspinType ; s++)
//        for ( sp = 0 ; sp < NS ; sp++)
//            for (sp2 = 0 ; sp2 < NS ; sp2++)
//                f1.arraySpin[s][sp][sp2] = -1;;
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
//                f1.arraySpin[full][sp][sp2] = si++;
//                f1.arraySpin[full][sp2][sp] = f1.arraySpin[full][sp][sp2];
//            }
//            f1.arraySpin[coulomb][sp][sp2] = 0;
//        }
//        f1.arraySpin[none][sp][sp] = 0;
//        f1.arraySpin[sym][sp][sp] = 0;
//        f1.arraySpin[diag][sp][sp] = sp;
//    }
//#else
//
////    if ( NS == 2 ){
////        f1.arraySpin[full][0][0] = 0;
////        f1.arraySpin[full][1][1] = 1;
////        f1.arraySpin[full][1][0] = 2;
////        f1.arraySpin[full][0][1] = 2;
////
////
////        f1.arraySpin[coulomb][0][0] = 0;
////        f1.arraySpin[coulomb][1][0] = 0;
////        f1.arraySpin[coulomb][0][1] = 0;
////        f1.arraySpin[coulomb][1][1] = 0;
////
////
////        f1.arraySpin[none][0][0] = 0;
////        f1.arraySpin[none][1][1] = 0;
////
////        f1.arraySpin[sym][0][0] = 0;
////        f1.arraySpin[sym][1][1] = 0;
////
////        f1.arraySpin[diag][0][0] = 0;
////        f1.arraySpin[diag][1][1] = 1;
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

INT_TYPE length ( struct sinc_label f1 , enum division label, INT_TYPE *lens ){
    
    INT_TYPE space ;
    for ( space = 0 ; space < SPACE ; space++)
        lens[space] = alloc(f1,label ,space);

    return 0;
}

INT_TYPE vectorLen(struct sinc_label f1, INT_TYPE space){
    if ( f1.rose[space].body == one )
        return f1.rose[space].count1Basis ;
    else if (f1.rose[space].body == two )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis  ;
    else if ( f1.rose[space].body == three )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis  ;
    else if ( f1.rose[space].body == four )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis  ;
    else
        return 0;
}
INT_TYPE outerVectorLen(struct sinc_label f1, enum bodyType bd, INT_TYPE space){
    if ( bd == one )
        return f1.rose[space].count1Basis ;
    else if ( bd == two )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis  ;
    else if ( bd == three )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis  ;
    else if ( bd == four )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis  ;
    else
        return 0;
}


INT_TYPE vector1Len(struct sinc_label f1, INT_TYPE space){
    return f1.rose[space].count1Basis ;
}

void length1(struct sinc_label f1, INT_TYPE *len){
    INT_TYPE space ;
    for ( space = 0 ;space < SPACE;space++)
        len[space] = vector1Len(f1, space);
    
    return;
}

INT_TYPE matrixLen(struct sinc_label f1, enum bodyType body, INT_TYPE space){
    if ( body == one )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis ;
    else if ( body == two )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis  ;
    else if ( body == three )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis  ;
    else if ( body == four )
        return f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis*f1.rose[space].count1Basis  ;
    else
        return 0;
}

INT_TYPE alloc ( struct sinc_label f1 , enum division label ,INT_TYPE space){
    if  ( space <0 || space > SPACE ){
        printf("allocy\n");
        exit(7);
    }
        if ( space == SPACE ){
        return 1;
    }else if ( species(f1, label ) == vector ){
            return vectorLen(f1, space);
        } else if ( species(f1, label ) == matrix ){
            return matrixLen(f1, f1.tulip[name(f1,label)].space[space].body, space);
        }else if ( species(f1, label ) == scalar ){
            return 1;
        }else if ( species(f1, label ) == outerVector ){
            return outerVectorLen(f1, f1.tulip[name(f1,label)].space[space].body, space);
        } else
            return 1;
}

INT_TYPE zero ( struct sinc_label f1 , enum division label, INT_TYPE spin ){
    //f1.tulip[label].Current = 0;
    INT_TYPE i, space,M2[SPACE];
    length(f1, label, M2);
    
    
    for ( space = 0; space < SPACE ;space++)
        if ( f1.rose[space].body != nada){
        Stream_Type * pt = streams(f1,label, (spin) , space );

        for ( i = 0; i < M2[space]*part(f1, label) ; i++ ){
            *(pt+i) = 0.;
        }
    }
    return 0 ;
}

INT_TYPE myZero ( struct sinc_label f1 , enum division label, INT_TYPE spin ){
    //f1.tulip[label].Current = 0;
    INT_TYPE i;

    for ( i = 0; i < part(f1, label) ; i++ ){
        myStreams(f1,label,spin)[i] = 0.;
    }
    return 0 ;
}


INT_TYPE tClear ( struct sinc_label f1 , enum division label ){
    INT_TYPE spin ;
    for ( spin = 0; spin < MaxCore ; spin++){
        f1.tulip[label].Current[spin] = 0;
    }
    return 0 ;
}

//double volume ( struct input * f1 ){
//    return pow( vectorLen(f1->c.sinc,0)*f1->d,3. );
//}

INT_TYPE CanonicalRank( struct sinc_label f1 , enum division label , INT_TYPE spin ){
    if ( label > f1.end ){
        printf("Can rank past end\n");
        fflush(stdout);
        exit(0);
    }
    struct name_label u = f1.tulip[label];
    
    if ( f1.tulip[label].name == label){
        if ( spin < spins(f1, label) ){
            return f1.tulip[label].Current[spin];
        }
            else
                return 0;
        }
    else {
        return f1.tulip[label].Partition;
    }
}

INT_TYPE Rank( struct sinc_label f1 , enum division label ){
    INT_TYPE sp,ra=0;;
    for ( sp = 0 ; sp < spins(f1, label);sp++)
        ra += CanonicalRank(f1, label, sp);
    return ra;
}



INT_TYPE spins ( struct sinc_label f1 , enum division label ){
    enum spinType sp = f1.tulip[label].spinor;

#ifndef APPLE
#ifdef OMP
    if (sp == parallel)
            return  MaxCore;
#endif
#endif
    
    
    if ( sp == real )
        return 1;
    else if ( sp == cmpl )
        return 2;
    else if (sp == parallel )
        return 1;
    return 0;
}


double sumSquare (struct sinc_label  f1,  enum division alloy){
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
                if ( f1.rose[dim].body != nada)
                product *= sqr(cblas_dnrm2(M2[dim], streams(f1, alloy,sp,dim)+l*M2[dim],iOne));;
            norm += product ;//* sqr(coeff(f1, alloy , sp,l ));
        }
    return norm;
}

void assignOneWithPointers( struct sinc_label f1, enum division oneMat , enum particleType particle){
    INT_TYPE space;
    enum block tv ;
    f1.tulip[oneMat].species = matrix;
    {
        for ( space = 0; space < SPACE ; space++)
            if ((f1.rose[space].body != nada ) &&( f1.rose[space].particle == particle || particle == all)){
                f1.tulip[oneMat].space[space].body = one;
            }else {
                f1.tulip[oneMat].space[space].body = nada;
            }
    }

    for ( tv = tv1 ; tv <= tv4 ; tv++){
        f1.tulip[oneMat+tv].Partition = part(f1, oneMat);
        f1.tulip[oneMat+tv].species = matrix;
        f1.tulip[oneMat+tv].name = oneMat;
        if ( tv < tv4 )
            f1.tulip[oneMat+tv].linkNext = oneMat+tv+1;
        for ( space = 0; space < SPACE ; space++)
            if ((f1.rose[space].body != nada ) &&( f1.rose[space].particle == particle || particle == all)){
                f1.tulip[oneMat+tv].space[space].block = tv;
            }else {
                f1.tulip[oneMat+tv].space[space].block = id0;
            }
    }
    
    
    
    return;
}

void assignParticle(struct sinc_label f1, enum division ma, enum particleType pa , enum bodyType ba ){
    INT_TYPE space ;
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.rose[space].body != nada )
        if ( f1.rose[space].particle == pa || pa == all )
            f1.tulip[ma].space[space].body = ba;
}

void assignTwoWithPointers( struct sinc_label f1, enum division twoMat ){
    INT_TYPE space;
    enum block e ;
    enum block maxee = e12;
    if ( bodies(f1, twoMat) == three )
        maxee = e23;
    if ( bodies(f1, twoMat) == four )
        maxee = e34;
    for ( e = e12 ; e <= maxee ; e++){
        f1.tulip[twoMat+e].name = twoMat;
        f1.tulip[twoMat+e].species = matrix;
        for ( space = 0; space < SPACE ; space++){
            f1.tulip[twoMat+e].space[space].block = e;
        }
    }
    return;
}

INT_TYPE sizeofDivision(struct sinc_label f1, enum division head, INT_TYPE space ){
    if ( f1.tulip[head].Partition == 0 )
        return 0;
    if ( f1.tulip[head].memory == bufferAllocation && space < SPACE )
        return 0;
    if ( f1.tulip[head].memory == objectAllocation && space == SPACE )
        return 0;
    if ( f1.tulip[head].memory == noAllocation)
        return 0;
    if ( head == f1.end )
        return 0;
    return alloc( f1, head, space)*part(f1,head)*spins(f1,head);
}


void  fromBeginning( struct sinc_label  f1 ,enum division new, enum division head ){
    if ( new > f1.end ){
        printf("from past end\n");
        fflush(stdout);
        exit(0);
    }
    INT_TYPE space;
    
    if ( head == nullName ){
        
        for ( space = 0; space <= SPACE ; space++)
            f1.tulip[new].space[space].Address = 0;
        
    }else {
        for ( space = 0; space <= SPACE ; space++){
            f1.tulip[new].space[space].Address = f1.tulip[head].space[space].Address;
            if ( f1.tulip[head].space[space].Address == -1 )
            {
                printf("acck!\n");
            }
        }
    }
    if ( f1.tulip[head].memory == objectAllocation){
        for ( space = 0; space < SPACE ; space++)
        if( f1.rose[space].body != nada){
            f1.tulip[new].space[space].Address += sizeofDivision(f1,head,space);
        }
    } else if ( f1.tulip[head].memory == bufferAllocation){
        for ( space = SPACE; space <= SPACE ; space++)
        {
            f1.tulip[new].space[space].Address += sizeofDivision(f1,head,space);
        }
    } else {
        //nothing
    }
//    printf("||%d::", new);
//    for ( space = 0; space <= SPACE ; space++)
//        printf("%lld:", f1.tulip[new].space[space].Address);
//    printf("\n\n");
    return;
}




//Stream_Type* myStreams ( struct field * f1, enum division label ,INT_TYPE spin ){
//    if ( memory(f1, label) != oneObject){
//        printf("called my stream %d %d\n",label,purpose(f1, label) );
//        exit(0);
//    }
//
//
//    if ( spin < 0 || spin >= spins(f1, label)){
//        printf("\n*my %d %lld\n\n", label, spin );
//        exit(0);
//    }
//    if ( spin != 0  ){
//        INT_TYPE leng = alloc(f1, name(f1,label));
//        INT_TYPE partit = part(f1, name(f1,label));
//        if ( purpose(f1, label ) == ptObject ){
//
//            printf("fixme\n");
//            exit(0);
//        }
//        return f1.rose[SPACE].stream+f1.tulip[name(f1,label)].myAddress + leng * partit * spin  ;
//
//    }
//
//    else{
//        if ( purpose(f1, label ) == ptObject ){
//            INT_TYPE len[SPACE];
//            length(f1, name(f1,label), len);
//
//            ///        printf("%lld -%lld %lld %lld\n", label,name(f1,label),len[space],f1.tulip[label].ptRank[spin]);
//            // / printf("%lld %lld\n", len[space]*part(f1,label), len[space]*f1.tulip[label].ptRank[spin]);
//            // p/rintf("%lld\n", f1.tulip[name(f1,label)].Address);
//            return f1.rose[SPACE].stream+f1.tulip[name(f1,label)].myAddress + len[0]*f1.tulip[label].ptRank[spin] ;
//
//        }else {
//
//            return f1.rose[SPACE].stream+f1.tulip[name(f1,label)].myAddress ;
//        }
//    }
//
//}

Stream_Type*  myStreams ( struct sinc_label f1, enum division label ,INT_TYPE spin){
    INT_TYPE space = SPACE;//buffer-space
    if ( spin < 0 || spin >= spins(f1, label)){
        printf("my spins %d %d\n",label ,spin);
        exit(5);
    }

    INT_TYPE leng = alloc(f1, name(f1,label),space);
    INT_TYPE partit = part(f1, name(f1,label));
    if (f1.tulip[name(f1,label)].space[SPACE].Address == -1  ){
        
        exit(0);
    }

    if ( f1.tulip[label].memory == bufferAllocation){
        if ( name(f1,label) != label ){
            return f1.rose[space].stream+f1.tulip[name(f1,label)].space[space].Address + leng * partit * spin + leng*f1.tulip[label].Current[spin] ;
        }
        else{
            return f1.rose[space].stream+f1.tulip[name(f1,label)].space[space].Address + leng * partit * spin  ;
        }
    }else {
        printf("not a buffer!\n %d",label);
        exit(0);
    }
    return NULL;
}


Stream_Type* streams ( struct sinc_label f1, enum division label ,INT_TYPE spin, INT_TYPE space ){
    Stream_Type * uu ;
    if ( spin < 0 || spin >= spins(f1, label)){
        printf("spins %d\n",label);
        exit(5);
    }
    if ( space < 0 || space >= SPACE+1 ){
        printf("space\n");

        exit(4);
        
    }
    if (f1.tulip[name(f1,label)].space[space].Address == -1  ){
        
        exit(0);
    }
    if ( f1.tulip[label].memory == objectAllocation){
        
        INT_TYPE leng = alloc(f1, name(f1,label),space);
        INT_TYPE partit = part(f1, name(f1,label));
        
        if ( name(f1,label) != label ){
            return f1.rose[space].stream+f1.tulip[name(f1,label)].space[space].Address + leng * partit * spin + leng*f1.tulip[label].Current[spin] ;
        }
        else{
             uu =  f1.rose[space].stream+f1.tulip[name(f1,label)].space[space].Address + leng * partit * spin  ;
            return uu;
        }
    } else if ( f1.tulip[label].memory == bufferAllocation){
        myStreams(f1, label, spin);
    }
    
    return NULL;
}

void assignView(INT_TYPE lane, struct sinc_label f1,  enum division A,INT_TYPE part ){
    tClear(f1, lanes+lane);
    f1.tulip[lanes+ lane].spinor =f1.tulip[A].spinor;
    f1.tulip[lanes+ lane].name = name(f1,A);
    f1.tulip[lanes+ lane].header = header(f1 ,A);
    f1.tulip[lanes+ lane].species = species(f1 ,A);
    f1.tulip[lanes+ lane].memory = f1.tulip[A].memory;
//    struct name_label u = f1.tulip[lanes+ lane];
    f1.tulip[lanes+ lane].Partition = part;
    struct name_label u = f1.tulip[lanes+lane];

}

void assignViewBlock(INT_TYPE lane, struct sinc_label f1,  enum division A ){
    INT_TYPE space;
    for ( space = 0 ; space < SPACE ; space++)
        f1.tulip[lanes+lane].space[space].block = f1.tulip[A].space[space].block;

}

enum division ocean(INT_TYPE lane, struct sinc_label f1, INT_TYPE l, INT_TYPE spin){
    enum division A = f1.tulip[lanes+ lane].name;
    f1.tulip[lanes+ lane].Current[spin] = l;
    if ( name(f1,A) != A )
        f1.tulip[lanes+ lane].Current[spin] += f1.tulip[A].Current[spin];

    return lanes+ lane;
}

void xsAdd ( double scalar , INT_TYPE dim ,struct sinc_label f1 , enum division targ ,INT_TYPE tspin,struct sinc_label  f2 , enum division orig,INT_TYPE o,INT_TYPE ospin ){
    INT_TYPE M2 = alloc(f1, orig, dim);
    INT_TYPE N2 = alloc(f1, targ, dim);
    INT_TYPE flag = (N2 == M2),space=dim;

    if ( flag && f2.tulip[orig].memory == objectAllocation && f1.tulip[targ].memory == objectAllocation){
        cblas_dcopy(N2, streams(f2,orig,ospin,space)+N2*o,1,streams(f1,targ,tspin,space)+CanonicalRank(f1, targ, tspin)*N2,1);
        if ( scalar != 1. )
            cblas_dscal(N2, scalar, streams(f1,targ,tspin,space)+CanonicalRank(f1, targ, tspin)*N2,1);
    }
    else {
        printf("add\n");

        exit(0);
    }
}

double xEqua ( struct sinc_label f1 , enum division targ ,INT_TYPE tspin,struct sinc_label  f2 , enum division orig,INT_TYPE ospin ){
    INT_TYPE space,flag=1;
    INT_TYPE eb = CanonicalRank(f2,orig,ospin);
    INT_TYPE M2[SPACE];
    length(f2, orig, M2);
    INT_TYPE N2[SPACE];
    length(f1,targ,N2);
    
    
    
    if ( f2.tulip[orig].memory == objectAllocation && f1.tulip[targ].memory == objectAllocation){
        
        if ( name(f1,targ) == targ  ) {
            
            if ( (part(f1,targ) < eb) )
                
            {
                printf("tEqual.. memory %d %d %d %d\n", targ,part(f1,targ), orig,part(f1,orig) );
                printf("partition %d\n", eb);
                exit(0);
            }
            
            for ( space = 0 ; space < SPACE ; space++)
                if ( f1.rose[space].body != nada)
                flag = flag * (N2[space]==M2[space]);

            
            if ( !flag ){
                if ( bodies(f2, orig) == one )
                    xOneBand(f2, orig,  ospin, f1, targ,tspin,0 );
                if ( bodies(f2, orig) == two )
                    xTwoBand(f2, orig,  ospin, f1, targ,tspin,0 );
                else if ( bodies(f2, orig) == three )
                    xThreeBand(f2, orig,  ospin, f1, targ,tspin,0 );
                else if ( bodies(f2, orig) == four )
                    xFourBand(f2, orig,  ospin, f1, targ,tspin,0 );
            }
            else {
                for ( space = 0; space < SPACE; space++)
                    if ( f1.rose[space].body != nada)
                    {
                    cblas_dcopy(eb*M2[space], streams(f2,orig,ospin,space),1,streams(f1,targ,tspin,space),1);
                }
            }
            struct name_label o = f1.tulip[orig];
            f1.tulip[targ].Current[tspin] = f2.tulip[orig].Current[ospin];
            f1.tulip[name(f1,targ)].header = header(f2, name(f2,orig));
        } else{
            printf("xeq 1");
            exit(9);
        }
    } else {
        printf("xeq 2");

        exit(9);
    }
    return 0;
}


double tEqua ( struct sinc_label f1 , enum division targ ,INT_TYPE tspin, enum division orig,INT_TYPE ospin ){
    return xEqua(f1, targ,tspin, f1, orig,ospin);
}

INT_TYPE tEquals( struct sinc_label f1 , enum division left , enum division right){
    enum spinType spl,spr;
    spl = f1.tulip[left].spinor;
    spr = f1.tulip[right].spinor;

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



INT_TYPE tAddTwo( struct sinc_label f1 , enum division left , enum division right){
    enum spinType spl,spr;
    spl = f1.tulip[left].spinor;
    spr = f1.tulip[right].spinor;
    
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

INT_TYPE tScaleOne( struct sinc_label f1, enum division label,INT_TYPE spin, double scalar ){
    
    if ( scalar == 1. )
        return 0;
    
    
    INT_TYPE L1 = CanonicalRank(f1,label,spin),space;
    
    double scale = fabs(scalar),prod;
    
    
    INT_TYPE M2[SPACE];
    length(f1,label,M2);
    
    prod = 1.;
    
    INT_TYPE dimCount = 0;
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.rose[space].body != nada ){
            dimCount++;
        }
            
    
    for ( space = 0 ; space < SPACE ; space++)
        if ( f1.rose[space].body != nada ){
            cblas_dscal(L1*M2[space],pow(scale,1./dimCount), streams(f1,label,spin,space),1);
        }
    
    if ( scalar < 0 )
        cblas_dscal(L1*M2[0],-1., streams(f1,label,spin,0),1);

    
    
    
    
    
//    
//    if ( SPACE == 3 ){
//        scale = pow ( scale , 1./3.);
//        if ( scalar < 0 )
//            scale *= -1;
//        cblas_dscal(L1*M2[0],scale, streams(f1,label,spin,0),1);
//        cblas_dscal(L1*M2[1],scale, streams(f1,label,spin,1),1);
//        cblas_dscal(L1*M2[2],scale, streams(f1,label,spin,2),1);
//    }else if ( SPACE == 2 ){
//        scale = pow ( scale , 1./2.);
//        cblas_dscal(L1*M2[1],scale, streams(f1,label,spin,1),1);
//        if ( scalar < 0 )
//            scale *= -1;
//        cblas_dscal(L1*M2[0],scale, streams(f1,label,spin,0),1);
//    }else if ( SPACE == 1 ){
//        cblas_dscal(L1*M2[0],scalar, streams(f1,label,spin,0),1);
//    }

    
    return 0;
}


INT_TYPE tScale( struct sinc_label f1, enum division label, double scalar ){
    INT_TYPE sp;
    for ( sp = 0; sp < spins(f1, label) ; sp++)
        tScaleOne(f1, label, sp,scalar);
    return 0;
}

INT_TYPE tAddTw( struct sinc_label  f1 , enum division left, INT_TYPE lspin,  enum division right , INT_TYPE rspin){
    return xAddTw(f1,left, lspin, f1, right, rspin);
}

INT_TYPE xAddTw( struct sinc_label f1 , enum division left, INT_TYPE lspin,struct sinc_label f2 ,  enum division right , INT_TYPE rspin){
    if ( CanonicalRank(f2,right,rspin) ){
        INT_TYPE LL = f1.tulip[left].Current[lspin];
        INT_TYPE LR = f2.tulip[right].Current[rspin];
        //printf("++ %lld %lld\n", LL,LR);
        INT_TYPE MM = LL+LR;
        if ( MM > f1.tulip[left].Partition ){
            //tGetType allocation!
            printf("tAdd more money! %d -> %d\n",left,right);
            exit(8);
        }
        INT_TYPE M2[SPACE],space;
        length(f1, left, M2);
        if ( f1.tulip[right].memory == objectAllocation ){
            for ( space = 0; space < SPACE; space++)
                if (f1.rose[space].body != nada){
                if ( (species(f1, right ) == vector) || (f1.tulip[right].space[space].body != nada && species(f1, right ) == matrix)){
                    cblas_dcopy(LR*M2[space], streams(f2,right,rspin,space),1,streams(f1,left,lspin,space)+LL*M2[space],1);
                }else if ((f1.tulip[right].space[space].body == nada && species(f1, right ) == matrix && f1.tulip[right].space[space].body == id0) ){
                    
                    INT_TYPE i,l,n1 = outerVectorLen(f1,bodies(f1,name(f1,right)),space);//question!
                    for ( i= 0 ; i < LR*M2[space] ; i++)
                        (streams(f1,left,lspin,space)+LL*M2[space])[i] = 0;
                    for ( l = 0; l < LR ; l++)
                        for ( i = 0; i < n1 ; i++)
                            (streams(f1,left,lspin,space)+(LL+l)*M2[space])[i*n1+i] = 1;
                }
            
            
            
                }
        }else if ( f1.tulip[right].memory == bufferAllocation )
        {
            printf("buf\n");
            exit(8);
        }
            else if ( f1.tulip[right].memory == noAllocation ){
            printf("not allocated\n");
            exit(0);
        }
        f1.tulip[left].Current[lspin] += LR;
    }
    return 0;
}

//INT_TYPE tSumMatrices(struct sinc_label f1, enum division sum ,INT_TYPE spin, enum division mat  ){
//
//    INT_TYPE n2[SPACE];
//    length(f1, sum,n2);
//    if ( bodies ( f1 , sum )== one ){
//
//        if (bodies(f1, mat ) == one ){
//            tAddTw(f1, sum ,spin, mat,0 );
//        }
//        else {
//            printf("eh?");
//            //exit(0);
//        }
//    }else if ( bodies ( f1 , sum )== two && f1->mem1->rt->calcType == electronicStuctureCalculation){
//        if (bodies(f1, mat ) == two ){
//            tAddTw(f1, sum ,spin, mat ,0);
//        }else if ( bodies ( f1, mat ) == one ){
//            INT_TYPE I1,I2, I3, I4,body,r,space;
//
//
//            INT_TYPE n1[SPACE];
//            length1(f1,n1);
//
//            double value;
//            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
//                for ( body = 0 ; body < 2 ; body++){
//
//                    for ( space = 0; space < SPACE ; space++){
//                        Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*CanonicalRank(f1, sum, spin);
//                        if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
//                            printf("part sum\n");
//                            exit(0);
//                        }
//
//                        for ( I1 = 0 ; I1 < n1[space] ; I1++)//body 0
//                            for ( I2 = 0 ; I2 < n1[space] ; I2++)//body 0
//                                for ( I3 = 0 ; I3 < n1[space] ; I3++)//body 1
//                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)//body 1
//                                    {
//                                        if ( body == 0 ){
//                                            value  = streams(f1,mat,0,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4);
//                                        }else {
//                                            value  = streams(f1,mat,0,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2);
//                                        }
//                                        stream[ (I1+I3*n1[space])+ ( I2+I4*n1[space])*n1[space]*n1[space] ] = value;
//                                    }
//                    }
//                    f1.tulip[sum].Current[spin]++;
//
//                    //  tAddTwo(f1, sum , quadCube);
//                }
//        }else {
//            printf("hey!\n");
//            exit(0);
//        }
//    }
//    else if ( bodies ( f1, sum ) == three && f1->mem1->rt->calcType == electronicStuctureCalculation){
//
//        if (bodies(f1, mat ) == three ){
//            tAddTw(f1, sum ,spin, mat,0 );
//        }else if ( bodies ( f1, mat ) == one ){
//            INT_TYPE I1,I2, I3, I4,I5,I6,body,r,space;
//            INT_TYPE n1[SPACE];
//            length1(f1,n1);
//            double value;
//
//            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
//                for ( body = 0 ; body < 3 ; body++){
//                    for ( space = 0; space < SPACE ; space++){
//                        Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*CanonicalRank(f1, sum, spin);
//                        if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
//                            printf("part sum\n");
//                            exit(0);
//                        }
//
//                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
//                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
//                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
//                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
//                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
//                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
//
//                                            {
//                                                value = 0;
//                                                if ( body == 0 ){
//                                                    value  = streams(f1,mat,0,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6);
//                                                }else if ( body == 1 ) {
//                                                    value  = streams(f1,mat,0,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6);
//                                                }else if ( body == 2 ) {
//                                                    value  = streams(f1,mat,0,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4);
//                                                }
//                                                stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]] = value;
//                                            }
//                    }
//                    f1.tulip[sum].Current[spin]++;
//
//                }
//        }else if ( bodies ( f1, mat ) == two ){
//            INT_TYPE I1,I2, I3, I4,I5,I6,pair,r,space,ve;
//            INT_TYPE n1[SPACE];
//
//            length1(f1, n1);
//
//            double value;
//
//            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
//                for ( pair = 0 ; pair < 3 ; pair++){
//                    for ( space = 0; space < SPACE ; space++){
//                        Stream_Type * stream = streams(f1,sum,spin,space)+n2[space]*CanonicalRank(f1, sum, spin);
//                        if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
//                            printf("part sum\n");
//                            exit(0);
//                        }
//                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
//                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
//                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
//                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
//                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
//                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
//
//                                            {
//                                                value = 0;
//                                                if ( pair == 0 ){
//                                                    value  = streams(f1,mat,0,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6);
//                                                }else if ( pair == 1 ) {
//                                                    value  = streams(f1,mat,0,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4);
//                                                }else if ( pair == 2 ) {
//                                                    value  = streams(f1,mat,0,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2);
//                                                }
//                                                ve = (I1+I3*n1[space]+I5*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space])*n1[space]*n1[space]*n1[space];
//                                                //                                                printf("%f %lld %lld\n", value,ve,n2[space]);
//                                                //                                                fflush(stdout);
//                                                stream[ ve ] = value;
//                                                //                                                printf("x");
//                                                //                                                fflush(stdout);
//
//                                            }
//                    }
//                    f1.tulip[sum].Current[spin]++;
//                }
//        }
//        else {
//            printf("Yo!");
//            exit(0);
//        }
//
//    }else if (bodies(f1,sum) == four && f1->mem1->rt->calcType == electronicStuctureCalculation){
//
//        if ( bodies ( f1, mat ) == one ){
//            INT_TYPE I1,I2, I3, I4,I5,I6,I7,I8,body,r,space;
//            INT_TYPE n1[SPACE];
//            length1(f1,n1);
//            double value;
//
//            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
//                for ( body = 0 ; body < 4 ; body++){
//                    for ( space = 0; space < SPACE ; space++){
//                        Stream_Type * stream = streams(f1,sum,0,space)+n2[space]*CanonicalRank(f1, sum, 0);
//                        if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
//                            printf("part sum\n");
//                            exit(0);
//                        }
//                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
//                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
//                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
//                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
//                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
//                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
//                                                for ( I7 = 0 ; I7 < n1[space] ; I7++)
//                                                    for ( I8 = 0 ; I8 < n1[space] ; I8++)
//
//                                                    {
//                                                        value = 0;
//                                                        if ( body == 0 ){
//                                                            value  = streams(f1,mat,0,space)[ I1*n1[space]+I2 + r*n1[space]*n1[space] ] * delta(I3-I4)*delta(I5-I6)*delta(I7-I8);
//                                                        }else if ( body == 1 ) {
//                                                            value  = streams(f1,mat,0,space)[ I3*n1[space]+I4 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I5-I6)*delta(I7-I8);
//                                                        }else if ( body == 2 ) {
//                                                            value  = streams(f1,mat,0,space)[ I5*n1[space]+I6 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I7-I8);
//                                                        }else if ( body == 3 ) {
//                                                            value  = streams(f1,mat,0,space)[ I7*n1[space]+I8 + r*n1[space]*n1[space] ] * delta(I1-I2)*delta(I3-I4)*delta(I5-I6);
//                                                        }
//                                                        stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
//                                                    }
//                    }
//                    f1.tulip[sum].Current[0]++;
//
//                }
//        }else if ( bodies ( f1, mat ) == two ){
//            INT_TYPE I1,I2, I3, I4,I5,I6,I7,I8,pair,r,space;
//            INT_TYPE n1[SPACE];
//            length1(f1, n1);
//
//            double value;
//
//            for ( r = 0; r < CanonicalRank(f1, mat , 0 ); r++)
//                for ( pair = 0 ; pair < 6 ; pair++){
//                    for ( space = 0; space < SPACE ; space++){
//                        Stream_Type * stream = streams(f1,sum,0,space)+n2[space]*CanonicalRank(f1, sum, 0);
//                        if ( CanonicalRank(f1, sum, spin) > part(f1, sum )){
//                            printf("part sum\n");
//                            exit(0);
//                        }
//                        for ( I1 = 0 ; I1 < n1[space] ; I1++)
//                            for ( I2 = 0 ; I2 < n1[space] ; I2++)
//                                for ( I3 = 0 ; I3 < n1[space] ; I3++)
//                                    for ( I4 = 0 ; I4 < n1[space] ; I4++)
//                                        for ( I5 = 0 ; I5 < n1[space] ; I5++)
//                                            for ( I6 = 0 ; I6 < n1[space] ; I6++)
//                                                for ( I7 = 0 ; I7 < n1[space] ; I7++)
//                                                    for ( I8 = 0 ; I8 < n1[space] ; I8++)
//
//                                                    {
//                                                        //                                            0    e12,     1,3     2,4
//                                                        //                                            1    e13,     1,5     2,6
//                                                        //                                            2    e23,     3,5     4,6
//                                                        //                                            3    e14,     1,7     2,8
//                                                        //                                            4    e24,     3,7     4,8
//                                                        //                                            5    e34      5,7     6,8
//
//
//                                                        if ( pair == 0 ){
//                                                            value  = streams(f1,mat,0,space)[ (I1*n1[space]+I3) + (I2*n1[space]+I4)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I5-I6)*delta(I7-I8);
//                                                        }else if ( pair == 1 ) {
//                                                            value  = streams(f1,mat,0,space)[ (I1*n1[space]+I5) + (I2*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I7-I8);
//                                                        }else if ( pair == 2 ) {
//                                                            value  = streams(f1,mat,0,space)[ (I3*n1[space]+I5) + (I4*n1[space]+I6)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I7-I8);
//                                                        }else if ( pair == 3 ) {
//                                                            value  = streams(f1,mat,0,space)[ (I1*n1[space]+I7) + (I2*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I3-I4)*delta(I5-I6);
//                                                        }else if ( pair == 4 ) {
//                                                            value  = streams(f1,mat,0,space)[ (I3*n1[space]+I7) + (I4*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I5-I6);
//                                                        }else if ( pair == 5 ) {
//                                                            value  = streams(f1,mat,0,space)[ (I5*n1[space]+I7) + (I6*n1[space]+I8)*n1[space]*n1[space] + r*n1[space]*n1[space]*n1[space]*n1[space] ]*delta(I1-I2)*delta(I3-I4);
//                                                        }else {
//                                                            printf ("rails!\n");
//                                                            exit(0);
//                                                        }
//                                                        stream[ (I1+I3*n1[space]+I5*n1[space]*n1[space]+I7*n1[space]*n1[space]*n1[space])+ (I2+I4*n1[space]+I6*n1[space]*n1[space]+I8*n1[space]*n1[space]*n1[space])*n1[space]*n1[space]*n1[space]*n1[space]] = value;
//                                                    }
//                    }
//                    f1.tulip[sum].Current[0]++;
//                }
//        }
//        else {
//            printf("Yo!");
//            exit(0);
//        }
//
//    }
//    return 0;
//}


INT_TYPE tAlt(struct sinc_label f1 , enum division label, INT_TYPE spin , INT_TYPE space1){
    INT_TYPE N1[SPACE] ;
    length1(f1,N1);
    
    INT_TYPE I1,I2,space;
    INT_TYPE Current    =     f1.tulip[label].Current[spin]++;
    
    
    if ( f1.tulip[label].species == vector ){
        printf("reconsider...only one mode represented\n");
    }
    if ( f1.tulip[label].species == matrix ){
        for ( space = 0; space < SPACE ; space++){
            Stream_Type * stream = streams(f1,label,spin,space)+Current*N1[space]*N1[space];
            for ( I1 = 0 ; I1 < N1[space] ; I1++)
                for ( I2 = 0 ; I2 < N1[space] ; I2++){
                    stream[I1*N1[space]+I2] =0.;
                }
            for ( I1 = 0 ; I1 < N1[space] ; I1++)
                for ( I2 = 0 ; I2 < N1[space] ; I2++){
                    {
                        if ( space == space1)
                            stream[I1*N1[space]+I2] = sign(I2-I1);
                        else{
                            if ( I1 == I2 && I1 == (N1[space]-1)/2 )
                                stream[I1*N1[space]+I2] = 1;
                        }
                    }
                    
                }
        }
    }
    
    
    return 1;
}

INT_TYPE tEnd(struct sinc_label f1 , enum division label, INT_TYPE spin , INT_TYPE space1){
    
    INT_TYPE N1[SPACE] ;
    length1(f1,N1);
    INT_TYPE I1,I2,space;
    INT_TYPE Current    =     f1.tulip[label].Current[spin]++;
    
    
    if ( f1.tulip[label].species == vector ){}
    if ( f1.tulip[label].species == matrix ){
        for ( space = 0; space < SPACE ; space++){
            
            
            Stream_Type * stream = streams(f1,label,spin,space)+Current*N1[space]*N1[space];
            for ( I1 = 0 ; I1 < N1[space] ; I1++)
                for ( I2 = 0 ; I2 < N1[space] ; I2++){
                    stream[I1*N1[space]+I2] =0.;
                }
            for ( I1 = 0 ; I1 < N1[space] ; I1++)
                for ( I2 = 0 ; I2 < N1[space] ; I2++){
                    {
                        if ( space == space1 )
                        {
                            if ( I1 == 0  )
                                if ( I1 == I2 )
                                    stream[I1*N1[space]+I2] = 1;
                            
                        }
                        else{
                            if ( I1 == I2)
                                stream[I1*N1[space]+I2] = 1;
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
//        tId(f1, PauliZ,f1.arraySpin[f1.tulip[PauliZ].spinor][0][0]);
//        tScale(f1, PauliZ,-1);
//        tId(f1, PauliZ,f1.arraySpin[f1.tulip[PauliZ].spinor][1][1]);
//    }
//    if (spins(f1,PauliZ)  >= 3 ){
//        tId(f1, PauliX,f1.arraySpin[f1.tulip[PauliX].spinor][0][1] );
//    }
//    return 0;
//}

INT_TYPE tId ( struct sinc_label f1 , enum division label,INT_TYPE spin ){
    
    INT_TYPE I1,I2,space;
    INT_TYPE Current ;
    {
        
        if ( f1.tulip[label].Current[spin] >= f1.tulip[label].Partition ){
            printf("%d %d\n", label, spin);
            printf("tryed to add to full array\n");
            return 0;
        }
        Current =  f1.tulip[label].Current[spin]++;
    }
    
    {
        if ( f1.tulip[label].species == vector || f1.tulip[label].species == outerVector){
            
            INT_TYPE B1[SPACE];
            length(f1, label, B1);
            for ( space = 0; space < SPACE ; space++)
                if ( f1.rose[space].body != nada)
                {
                
                Stream_Type  * stream = streams(f1,label,spin,space)+Current*B1[space];
                for ( I2 = 0 ; I2 < B1[space] ; I2++){
                    stream[I2] = sign(I2);
                }
            }
        }
        
        else if  ( f1.tulip[label].species == matrix ) {
            INT_TYPE B1;
            

            
            for ( space = 0; space < SPACE ; space++)
                if ( f1.rose[space].body != nada)
                {
                B1 = outerVectorLen(f1,bodies(f1, name(f1,label)), space);
                Stream_Type * stream = streams(f1,label,spin,space)+Current*B1*B1;
                for ( I1 = 0 ; I1 < B1 ; I1++)
                    for ( I2 = 0 ; I2 < B1 ; I2++)
                        stream[I1*B1+I2] =0.;
                for ( I1 = 0 ; I1 < B1 ; I1++)
                    {
                        
                        stream[I1*B1+I1] = 1;
                        
                    }

                }
        }
        
        
    }
    return 0;
}

INT_TYPE tOv ( struct sinc_label f1 , enum division label,INT_TYPE spin ){
    
    INT_TYPE I1,I2,space;
    INT_TYPE Current ;
    {
        
        if ( f1.tulip[label].Current[spin] >= f1.tulip[label].Partition ){
            printf("%d %d\n", label, spin);
            printf("tryed to add to full array\n");
            return 0;
        }
        Current =  f1.tulip[label].Current[spin]++;
    }
    
    {
        if ( f1.tulip[label].species == vector || f1.tulip[label].species == outerVector){
            
            INT_TYPE B1[SPACE];
            length(f1, label, B1);
            for ( space = 0; space < SPACE ; space++)
                if ( f1.rose[space].body != nada)
                {
                    
                    Stream_Type  * stream = streams(f1,label,spin,space)+Current*B1[space];
                    for ( I2 = 0 ; I2 < B1[space] ; I2++){
                        stream[I2] = sign(I2);
                    }
                }
        }
        
        else if  ( f1.tulip[label].species == matrix && bodies(f1, label) == one) {
            INT_TYPE B1;
            
            
            
            for ( space = 0; space < SPACE ; space++)
                if ( f1.rose[space].body != nada)
                {
                    B1 = outerVectorLen(f1,bodies(f1, name(f1,label)), space);
                    Stream_Type * stream = streams(f1,label,spin,space)+Current*B1*B1;
                    for ( I1 = 0 ; I1 < B1 ; I1++)
                        for ( I2 = 0 ; I2 < B1 ; I2++)
                            stream[I1*B1+I2] =BoB(f1.rose[space].basisList[I1], f1.rose[space].basisList[I2]);
                    
                }
        }
        
        
    }
    return 0;
}


INT_TYPE tReplace( struct sinc_label f1 , enum division label,INT_TYPE spin,INT_TYPE space,INT_TYPE l ){
    
    INT_TYPE I1,I2;
    {
        
        if ( f1.tulip[label].Current[spin] < l )
            return 0;
       }
    {
        if ( f1.tulip[label].species == vector ){
            INT_TYPE B1[SPACE];
            length(f1, label, B1);
            {
                Stream_Type  * stream = streams(f1,label,spin,space)+l*B1[space];
                for ( I2 = 0 ; I2 < B1[space] ; I2++){
                    stream[I2] = sign(I2%2);
                }
            }
        }
        
        //        else if ( ( f1.tulip[label].species == matrix && bodies(f1,label)== two) ||  (f1.tulip[label].species == quartic && bodies(f1,label) == one)){
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
        else if  ( f1.tulip[label].species == matrix ) {
            INT_TYPE B1[SPACE];
            length(f1, label,B1);
            
            

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

INT_TYPE tBoot ( struct sinc_label f1 , enum division label,INT_TYPE spin ){
    
    INT_TYPE I1,I2,space;
    INT_TYPE Current ;
    {
        
        if ( f1.tulip[label].Current[spin] >= f1.tulip[label].Partition )
            return 0;
        Current =  f1.tulip[label].Current[spin]++;
    }
    {
        if ( f1.tulip[label].species == vector && bodies(f1,label) == two){
            INT_TYPE B1[SPACE];
            length(f1, label, B1);
            for ( space = 0; space < SPACE ; space++){
                
                Stream_Type  * stream = streams(f1,label,spin,space)+Current*B1[space];
                for ( I1 = 0 ; I1< B1[space] ; I1++)
                    for ( I2 = 0 ; I2 < B1[space] ; I2++){
                        stream[I1*B1[space]+I2] = exp(-abs(I1-(B1[space]-1)/2)*0.01)*exp(-sqr(I2-(B1[space]-1)/2)*0.01);
                }
            }
        }
        else if ( f1.tulip[label].species == vector && bodies(f1,label) == one){
            INT_TYPE B1[SPACE];
            length(f1, label, B1);
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


//double vectorElement (struct field * f1, enum division state, INT_TYPE l1,INT_TYPE l2 , INT_TYPE l3 ){
//    INT_TYPE spin ,space;
//    double den = 0.;
//    double maxDen = 0.;
//    INT_TYPE l;
//    INT_TYPE n1[SPACE];
//    for ( space = 0; space < SPACE ; space++)
//        n1[space] = vectorLen(f1, space);
//    for ( spin = 0; spin < spins(f1, state ) ;spin++){
//        for ( l = 0; l < CanonicalRank(f1, state, spin ) ; l++)
//            den += (streams(f1, state, spin , 0 )[l*n1[0]+l1]* streams(f1, state, spin , 1 )[l*n1[1]+l2]* streams(f1, state, spin , 2 )[l*n1[2]+l3]);
//        if ( den > maxDen)
//            maxDen = den;
//    }
//
//    return den;
//}

double matrixElement (struct sinc_label f1, enum division label, INT_TYPE i , INT_TYPE i2, INT_TYPE j,INT_TYPE j2, INT_TYPE k , INT_TYPE k2 ){
    INT_TYPE N2[SPACE],n1[SPACE];
    INT_TYPE space;
    length1(f1, n1);
    for ( space = 0; space < 3 ; space++){
        N2[space] = n1[space]*n1[space];
    }

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

//void nuclearArray (struct input * f1,  enum division array,INT_TYPE M1){
//    INT_TYPE n1[SPACE] ;
//    length1(f1->c.sinc,n1);
//    INT_TYPE N1[SPACE];
//    INT_TYPE l = 2;
//    length(f1->c.sinc, oneVector, N1);
//    INT_TYPE i,k,j,at,mi,mj,mk;
//    double * array3d = myStreams(f1->c.sinc, array, 0),x,y,z,d,mini,dis;
//
//    for ( i = 0; i < M1*M1*M1 ; i++)
//        array3d[i] = 0.;
//
//
//    for ( at = 1 ; at <= f1->Na ; at++){
//        mini = M1;
//        mi = 0;
//        mj = 0;
//        mk = 0;
//        d = 1./ (double)(M1)  * (double)(n1[0]) * f1->d;
//        for ( i = 0; i < M1 ; i++){
//            for ( j = 0; j < M1 ; j++){
//                for ( k = 0 ; k < M1 ; k++){
//                    x = ( i - (M1-1)/2 ) * d ;
//                    y = ( j - (M1-1)/2 ) * d;
//                    z = ( k - (M1-1)/2 ) * d;
//
//                    dis = sqr( x - f1->atoms[at].position[1] ) + sqr( y - f1->atoms[at].position[2] ) + sqr( z - f1->atoms[at].position[3] );
//                    if (  dis < mini ){
//                        mini = dis;
//                        mi= i;
//                        mj = j;
//                        mk = k;
//                    }
//                }
//            }
//        }
//        for ( i = imax(0,mi-l); i < imin(M1,mi+l); i++)
//            for ( j = imax(0,mj-l); j < imin(M1,mj+l); j++)
//                for ( k = imax(0,mk-l); k < imin(M1,mk+l); k++)
//                    if ( sqr(mi-i) + sqr( mj-j)+ sqr(k-mk) < l)
//                        array3d[i+j*M1+k*M1*M1] += 1.;
//
//    }
//
//
//
//}

void vectorArray (struct sinc_label f1, enum division oneVector,  enum division array,INT_TYPE M1){
    INT_TYPE n1[SPACE];
    length1(f1,n1);
    INT_TYPE space,N1[SPACE],stride[SPACE];
    length(f1, oneVector, N1);
    INT_TYPE i;
    double * array3d = myStreams(f1, array, 0);
    INT_TYPE spin = 0,r;;
    double * basis = myStreams(f1, oneBasis,0);
    double * oneDim = myStreams(f1, oneArray,0);
    double * tempArea = oneDim + 3 * M1;
    if ( part(f1, oneArray ) < 3*M1+M1*M1){
        printf("vectorArray\n");
        exit(0);
    }
    for ( i = 0; i < n1[0]*n1[1]*n1[2] ; i++)
        array3d[i] = 0.;
    if ( bodies (f1,oneVector) == one && species(f1, oneVector ) == vector){
        stride[0] = 1;
        stride[1] = 1;
        stride[2] = 1;
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

INT_TYPE  countLinesFromFile( struct field f1,INT_TYPE location, INT_TYPE *ir ){
    INT_TYPE fi,lines = 0;
    size_t ms = MAXSTRING;
    char line0[MAXSTRING];
    char name[MAXSTRING];
    char name2[MAXSTRING];

    char *line = line0;
    INT_TYPE FIT ;
    if ( location == 0 )
        FIT = f1.i.files ;
    else if ( location == 1 )
        FIT = f1.i.filesVectorOperator ;
    else
        exit(1);
    
        for ( fi =0 ; fi < FIT; fi++){
            if ( location == 0 )
                strcpy(name,f1.i.fileList[fi] );
            else if ( location == 1)
                strcpy(name,f1.i.fileVectorOperator[fi] );
            else
                exit(0);
            FILE * fp = fopen(name,"r");
            if ( fp == NULL ) {
                printf("cannot find file: %s\n", name);
                exit(0);
            }
            getline(&line, &ms, fp);
            while(!feof(fp))
            {
                if (! comment(line))
                {
                   // printf("%s\n", line);
                  //  fflush(stdout);
                    tFromReadToFilename(NULL, line,  name2, 0,0);
                   // printf("%s\n", name2);
                   // fflush(stdout);

                    *ir = imax(*ir,inputFormat(f1.f, name2, nullName, 2));
                   // printf("%d \n", c1->i.iRank );
                   // fflush(stdout);

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
    return lines;
}


INT_TYPE assignCores(struct sinc_label f1, INT_TYPE parallel ){
#ifdef OMP
    INT_TYPE nSlot = f1.rt->NSlot;
#ifdef MKL
    INT_TYPE nParallel = f1.rt->NParallel;
#endif
    INT_TYPE nLanes = f1.rt->NLanes;

    INT_TYPE omp;
    if ( parallel == 0){
        omp = 1;
    }else if ( parallel ){
        omp = nLanes;
    }
    omp_set_num_threads(omp);
#endif
#ifdef MKL
    if ( omp == 0 )
        mkl_set_num_threads(nParallel*nLanes);
    else
        mkl_set_num_threads(nParallel);
#endif

    return 0;

}

//double lattice ( struct input * f1, INT_TYPE space ){
//    return f1->d;
//}

void printVectorAllocations(struct sinc_label f1){
    INT_TYPE space;
    for ( space = 0; space <= SPACE ; space++){
        ADDRESS_TYPE vecG = 3*(f1.tulip[f1.end].space[space].Address -  f1.tulip[eigenVectors].space[space].Address);
        printf("%d --vector Contribution \t G%f",space,vecG/1000000000./(sizeof(Stream_Type)));
    }
}

struct basisElement defineSincBasis (enum noteType note, enum componentType component, enum basisElementType basis, double lattice , double origin, INT_TYPE count1, INT_TYPE elementIndex ){
    struct basisElement boa;
    
    
    boa.note = note;
    
    boa.type = component;
    
    boa.basis = basis;
    
    boa.length = lattice;
    
    boa.origin = origin;
    
    boa.grid = count1;
    
    boa.index = 0;
    boa.index2 = 0;
    
    if ( boa.type > 3 ){
        boa.grid /= 2;//allocation should be 2*actual grid.!
    }
    if ( boa.type > 3 ) {
        if ( elementIndex >= boa.grid   ) {
            //periodic domain...
            INT_TYPE n,ct = 0 ;
            for ( n = - boa.grid ; n < boa.grid ; n+= 2){//see Ewald 7.5 MATHEMATICA
                if ( ct == elementIndex - boa.grid ){
                    boa.index2 = n;
                   // boa.type += 3;
                    break;
                }
                ct++;
            }
        }else {
            boa.index = elementIndex;
            //not boosted..
        }
        
    } else {
        boa.index = elementIndex - (boa.grid-1)/2 ;
    }
    return boa;

}

struct basisElement defineGaussBasis (enum noteType note, enum componentType component, enum basisElementType basis, double lattice , double origin, INT_TYPE count1, INT_TYPE elementIndex ){
    struct basisElement boa;
    
    
    boa.note = note;
    
    boa.type = component;
    
    boa.basis = basis;
    
    boa.length = lattice;
    
    boa.origin = origin;
    
    boa.grid = count1;
    
    boa.index = elementIndex;
    boa.index2 = 0;
    
    return boa;
    
}

struct basisElement grabBasis (struct sinc_label f1, INT_TYPE space, INT_TYPE particle, INT_TYPE elementIndex){
    if (  0 <= elementIndex && 0 <= space && space < SPACE && electron<= particle && particle < PARTICLE+1){
        if ( elementIndex < f1.rose[space].count1Basis  )
            return f1.rose[space].basisList[elementIndex];
    }
    
    
    exit(0);
}
        
//double interpolateOne( struct field * f1 , double * position ,enum division vec){
//    double *pt[SPACE];
//    INT_TYPE n1[SPACE],info,space;
//    if ( vec == copyVector){
//        printf("cpV");
//        exit(1);
//    }
//    for ( space = 0; space < SPACE ; space++)
//        pt[space] = streams(f1, copyVector, 0, space);
//    f1.tulip[copyVector].Current[0] = 1;
//    length1(f1, n1);
//    
//    
//    //tFillBasis( pt,position, 0, n1[3], f1.rose[3].lattice);
//    return tMultiplyMP(0, &info, f1, 1., -1, nullName, 0, 'T', copyVector, 0, 'N', vec, 0);
//}

struct basisElement transformBasis( INT_TYPE flip, double scale, struct basisElement ba ){
    ba.length *= scale;

    if ( flip ) {
        ba.origin *= -1.;
        ba.index *= -1;
        return ba;
    }
    
    return ba;
    
}


#if 0

double xOneBand (INT_TYPE rank,struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic){
    double D = f1.d;
    INT_TYPE space,i,l,i2,l2,r,aPeriodic;
    INT_TYPE N1 = f2.N1;
    INT_TYPE N12 = (N1-1)/2;
    INT_TYPE L1 =  f1.N1;
    INT_TYPE L12 = (L1-1)/2;
    f2.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);
    
    //    for ( r= 0 ; r < 3 ;r++){
    //    for ( l= 0 ; l < 9 ; l++)
    //        printf("%f:", streams(f1, vector1,s1,0)[l]);
    //        printf("\n");
    //    }
    //        printf("%d %f : %d %f\n", L1, D, N1, f2.d);
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
                        printf("ask me\n");
                        exit(0);
                            
                        }else {
                            streams(f1, foundationStructure,rank, space )[l] =
                            SS ( D,D*(l-L12), f2.d, f2.d*( i-N12 ) );
                            
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
    f2.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    
    
    
    return 0.;
}

double xTwoBand (INT_TYPE rank,struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic){
    double D = f1.d;
    INT_TYPE space,i,l,i2,l2,r,aPeriodic;
    INT_TYPE N1 = f2.N1;
    INT_TYPE N12 = (N1-1)/2;
    INT_TYPE L1 =  f1.N1;
    INT_TYPE L12 = (L1-1)/2;
    f2.tulip[out].Current[s2] = 0;
    zero(f1,out,s2);
    
//    for ( r= 0 ; r < 3 ;r++){
//    for ( l= 0 ; l < 9 ; l++)
//        printf("%f:", streams(f1, vector1,s1,0)[l]);
//        printf("\n");
//    }
//        printf("%d %f : %d %f\n", L1, D, N1, f2.d);
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
                            printf("ask me\n");
                            exit(0);

//                            streams(f1, foundationStructure,rank, space )[l2*L1+l] =
//                            periodicSS ( D,D*(l-L12),f1.N1, f2.d, f2.d*( i-N12 ) ,f2.N1)*
//                            periodicSS ( D,D*(l2-L12),f1.N1, f2.d, f2.d*( i2-N12 ),f2.N1 );
//                            periodicSS ( D,D*(l),f1.N1, f2.d, f2.d*( i ) ,f2.N1)*
//                            periodicSS ( D,D*(l2),f1.N1, f2.d, f2.d*( i2 ),f2.N1 );

                        }else {
                            streams(f1, foundationStructure,rank, space )[l2*L1+l] =
                            SS ( D,D*(l-L12), f2.d, f2.d*( i-N12 ) )*
                            SS ( D,D*(l2-L12), f2.d, f2.d*( i2-N12 ) );
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
    f2.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    
    
    
    return 0.;
}


double xThreeBand (INT_TYPE rank,struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic){
    double D = f1.d;
    INT_TYPE space,i,l,i2,i3,l2,l3,r,aPeriodic;
    INT_TYPE N1 = f2.N1;
    INT_TYPE N12 = (N1-1)/2;
    INT_TYPE L1 =  f1.N1;
    INT_TYPE L12 = (L1-1)/2;
    f2.tulip[out].Current[s2] = 0;
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
                                    printf("ask me\n");
                                    exit(0);

//                                    streams(f1, foundationStructure,rank, space )[(l3*L1*L1+l2*L1+l)] =
//                                    periodicSS ( D,D*(l-L12),f1.N1,f2.d, f2.d*( i-N12 ) ,f2.N1)*
//                                    periodicSS ( D,D*(l2-L12),f1.N1,f2.d, f2.d*( i2-N12 ),f2.N1 )*
//                                    periodicSS ( D,D*(l3-L12),f1.N1, f2.d, f2.d*( i3-N12 ) ,f2.N1);
                                }else {
                                    streams(f1, foundationStructure,rank, space )[(l3*L1*L1+l2*L1+l)] =
                                    SS ( D,D*(l-L12), f2.d, f2.d*( i-N12 ) )*
                                    SS ( D,D*(l2-L12), f2.d, f2.d*( i2-N12 ) )*
                                    SS ( D,D*(l3-L12), f2.d, f2.d*( i3-N12 ) );
                                }
                }
                
                for ( space = 0;space < SPACE; space++){
                    
                    for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++)
                        streams(f2, out, s2,space)[r*N1*N1*N1 + (i3*N1*N1+i2*N1+i)] = cblas_ddot(L1*L1*L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1*L1*L1,1);
                    //
                }
                
                
                
            }
    f2.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    return 0.;
}

double xFourBand (INT_TYPE rank,struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct field * f2, enum division out,INT_TYPE s2, INT_TYPE periodic){
    double D = f1.d;
    INT_TYPE space,i,l,i2,i3,i4,l2,l3,l4,r,aPeriodic;
    INT_TYPE N1 = f2.N1;
    INT_TYPE N12 = (N1-1)/2;
    INT_TYPE L1 =  f1.N1;
    INT_TYPE L12 = (L1-1)/2;
    f2.tulip[out].Current[s2] = 0;
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
                                    printf("ask me\n");
                                    exit(0);

//                                    streams(f1, foundationStructure,rank, space )[(l4*L1*L1*L1+l3*L1*L1+l2*L1+l)] =
//                                    periodicSS ( D,D*(l-L12),f1.N1,f2.d, f2.d*( i-N12 ) ,f2.N1)*
//                                    periodicSS ( D,D*(l2-L12),f1.N1,f2.d, f2.d*( i2-N12 ),f2.N1 )*
//                                    periodicSS ( D,D*(l3-L12),f1.N1, f2.d, f2.d*( i3-N12 ) ,f2.N1)*
//                                    periodicSS ( D,D*(l4-L12),f1.N1, f2.d, f2.d*( i4-N12 ) ,f2.N1);
                                }else {
                                    streams(f1, foundationStructure,rank, space )[(l4*L1*L1*L1+l3*L1*L1+l2*L1+l)] =
                                    SS ( D,D*(l-L12), f2.d, f2.d*( i-N12 ) )*
                                    SS ( D,D*(l2-L12), f2.d, f2.d*( i2-N12 ) )*
                                    SS ( D,D*(l3-L12), f2.d, f2.d*( i3-N12 ) )*
                                    SS ( D,D*(l4-L12), f2.d, f2.d*( i4-N12 ) );
                                }
                }
                
                for ( space = 0;space < SPACE; space++){
                    
                    for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++)
                        streams(f2, out, s2,space)[r*N1*N1*N1*N1 + (i4*N1*N1*N1+i3*N1*N1+i2*N1+i)] = cblas_ddot(L1*L1*L1*L1, streams(f1, foundationStructure,rank, space ),1,streams(f1, vector1,s1,space)+r*L1*L1*L1*L1,1);
                    //
                }
                
                
                
            }
        f2.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
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

double xOneBand (struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct sinc_label f2, enum division out,INT_TYPE s2,INT_TYPE oldPeriodic){
    INT_TYPE sl,space,i,l,r,rank;
    INT_TYPE n1[SPACE],N1;
    length1(f1,n1);
    INT_TYPE n2[SPACE];
    length1(f2,n2);
    
    INT_TYPE L1;
    f2.tulip[out].Current[s2] = 0;
    zero(f2,out,s2);
    
    for ( space = 0;space < SPACE; space++)
        if ( f1.rose[space].body != nada){
            N1 = n2[space];
            L1 = n1[space];
            
            for ( sl = 0; sl < L1 ; sl++){
                rank = 0;
                l = sl;
                
                //build
#ifdef OMP
#pragma omp parallel for private (i)
#endif
                for ( i = 0 ; i < N1 ; i++)
                {
                    myStreams(f2, bandBasis,rank )[i] =
                    BoB (grabBasis(f1, space, f1.rose[space].particle,l),grabBasis(f2, space, f2.rose[space].particle, i) );
                }
#ifdef OMP
#pragma omp parallel for private (r)
#endif
                for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                    cblas_daxpy(N1, (streams(f1, vector1,s1,space)+r*L1)[l], myStreams(f2, bandBasis, rank), 1, streams(f2, out, s2,space)+r*N1, 1);
                }
            }
        }
    f2.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
//    double ov ;
//    matrixElements(rank, f2, out, nullName, out, NULL, &ov);
//
//    double x, xv[SPACE];
//    for ( x = -10 ; x <= 10 ; x+=0.1){
//        for (space = 0; space < SPACE; space++)
//            xv[space] = 0.;
//        xv[3] = x;
//        printf("{%f,%f,%f},\n", x, interpolateOne(f2, xv, out),interpolateOne(f1,xv, vector1) );
//    }
//
    return 0.;
}

double xTwoBand (struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct sinc_label  f2, enum division out,INT_TYPE s2,INT_TYPE oldPeriodic){
    INT_TYPE space,i,l,i2,l2,r,rank,sl;
    INT_TYPE N1,si ;
    INT_TYPE L1 ;
    f2.tulip[out].Current[s2] = 0;
    zero(f2,out,s2);
    INT_TYPE n1[SPACE];
    length1(f1,n1);
    INT_TYPE n2[SPACE];
    length1(f2,n2);

    for ( space = 0;space < SPACE; space++)
    if ( f1.rose[space].body != nada)
    {
        N1 = n2[space];
        L1 = n1[space];
        
        for ( sl = 0; sl < L1*L1 ; sl++){
            rank = 0;
            l = sl% L1;
            l2 = (sl/L1)%L1;

            //build
#ifdef OMP
#pragma omp parallel for private (si,i,i2)
#endif
            for ( si = 0 ; si < N1*N1 ; si++)
            {
                i = si% N1;
                i2 = (si/N1)%N1;

                myStreams(f2, bandBasis,rank )[si] =
                BoB (grabBasis(f1, space, f1.rose[space].particle,l),grabBasis(f2, space, f2.rose[space].particle, i) )*
                BoB (grabBasis(f1, space, f1.rose[space].particle,l2),grabBasis(f2, space, f2.rose[space].particle, i2) );
            }
#ifdef OMP
#pragma omp parallel for private (r)
#endif
            for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                cblas_daxpy(N1*N1, (streams(f1, vector1,s1,space)+r*L1*L1)[sl], myStreams(f2, bandBasis, rank), 1, streams(f2, out, s2,space)+r*N1*N1, 1);
            }
        }
    }
    f2.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    return 0.;
}


double xThreeBand (struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct sinc_label  f2, enum division out,INT_TYPE s2,INT_TYPE oldPeriodic){
    INT_TYPE si,space,i,l,i2,i3,l2,l3,r,rank,sl;
    INT_TYPE N1 ;
    INT_TYPE L1 ;
    f2.tulip[out].Current[s2] = 0;
    zero(f2,out,s2);
    INT_TYPE n1[SPACE];
    length1(f1,n1);
    INT_TYPE n2[SPACE];
    length1(f2,n2);

    for ( space = 0;space < SPACE; space++)
        if ( f1.rose[space].body != nada){
            N1 = n2[space];
            L1 = n1[space];
            
            for ( sl = 0; sl < L1*L1*L1 ; sl++){
                rank = 0;
                l = sl% L1;
                l2 = (sl/L1)%L1;
                l3 = (sl/(L1*L1))%L1;

                //build
#ifdef OMP
#pragma omp parallel for private (si,i,i2,i3)
#endif
                for ( si = 0 ; si < N1*N1*N1 ; si++)
                {
                    i = si% N1;
                    i2 = (si/N1)%N1;
                    i3 = (si/(N1*N1))%N1;

                    myStreams(f2, bandBasis,rank )[si] =
                    BoB (grabBasis(f1, space, f1.rose[space].particle,l),grabBasis(f2, space, f2.rose[space].particle, i) )*
                    BoB (grabBasis(f1, space, f1.rose[space].particle,l2),grabBasis(f2, space, f2.rose[space].particle, i2) )*
                    BoB (grabBasis(f1, space, f1.rose[space].particle,l3),grabBasis(f2, space, f2.rose[space].particle, i3) );
                    
                }
#ifdef OMP
#pragma omp parallel for private (r)
#endif
                for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                    cblas_daxpy(N1*N1*N1, (streams(f1, vector1,s1,space)+r*L1*L1*L1)[sl], myStreams(f2, bandBasis, rank), 1, streams(f2, out, s2,space)+r*N1*N1*N1, 1);
                }
            }
        }
    f2.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    return 0.;
}

double xFourBand (struct sinc_label f1, enum division vector1 ,INT_TYPE s1, struct sinc_label f2, enum division out,INT_TYPE s2,INT_TYPE oldPeriodic){
    INT_TYPE si,space,i,l,i2,i3,i4,l2,l3,l4,r,rank,sl;
    INT_TYPE N1 ;
    INT_TYPE L1 ;
    f2.tulip[out].Current[s2] = 0;
    zero(f2,out,s2);
    INT_TYPE n1[SPACE];
    length1(f1,n1);
    INT_TYPE n2[SPACE];
    length1(f2,n2);

    
    for ( space = 0;space < SPACE; space++)
        if ( f1.rose[space].body != nada){
            N1 = n2[space];
            L1 = n1[space];
            
            for ( sl = 0; sl < L1*L1*L1*L1 ; sl++){
                rank = 0;
                l = sl% L1;
                l2 = (sl/L1)%L1;
                l3 = (sl/(L1*L1))%L1;
                l4 = (sl/(L1*L1*L1))%L1;

                //build
#ifdef OMP
#pragma omp parallel for private (si,i,i2,i3,i4)
#endif
                for ( si = 0 ; si < N1*N1*N1*N1 ; si++)
                {
                    i = si% N1;
                    i2 = (si/N1)%N1;
                    i3 = (si/(N1*N1))%N1;
                    i4 = (si/(N1*N1*N1))%N1;

                    myStreams(f2, bandBasis,rank )[si] =
                    BoB (grabBasis(f1, space, f1.rose[space].particle,l),grabBasis(f2, space, f2.rose[space].particle, i) )*
                    BoB (grabBasis(f1, space, f1.rose[space].particle,l2),grabBasis(f2, space, f2.rose[space].particle, i2) )*
                    BoB (grabBasis(f1, space, f1.rose[space].particle,l3),grabBasis(f2, space, f2.rose[space].particle, i3) )*
                    BoB (grabBasis(f1, space, f1.rose[space].particle,l4),grabBasis(f2, space, f2.rose[space].particle, i4) );

                    
                }
#ifdef OMP
#pragma omp parallel for private (r)
#endif
                for ( r = 0 ; r < CanonicalRank(f1, vector1, s1); r++){
                    cblas_daxpy(N1*N1*N1*N1, (streams(f1, vector1,s1,space)+r*L1*L1*L1*L1)[sl], myStreams(f2, bandBasis, rank), 1, streams(f2, out, s2,space)+r*N1*N1*N1*N1, 1);
                }
            }
        }
    f2.tulip[out].Current[s2] = CanonicalRank(f1, vector1, s1);
    
    return 0.;
}

#endif


