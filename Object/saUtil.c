/**
 *  saUtil.c
 *
 *
 *  Copyright 2021 Jonathan Jerke and Bill Poirier.
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

#include "saUtil.h"

/**
 *translate commands given into strides of a component array
 */
bodyType commandSA(  bodyType bd, inta act,   blockType cl,   blockType bl,inta perm[], inta op[]){
    
    perm[0] = 0;
    perm[1] = 1;
    perm[2] = 2;
    perm[3] = 3;
    perm[4] = 4;
    perm[5] = 5;

    
    
    switch ( bd ){
        case one:
            perm[0]     = 0;
            op[0]    = 0;
            return one;
        case two:
            if ( act == 1 ){
                perm[0] = 0;
                perm[1] = 1;
            }else if ( act ==2 ){
                perm[0] = 1;
                perm[1] = 0;
            }
            if ( cl == tv1){
                switch(bl){
                    case id0:
                //1-body
                    case tv1:
                        op[0] = 0;
                        op[1] = 1;
                        return one;
                    case tv2:
                        op[0] = 1;
                        op[1] = 0;
                        return one;
                //1-body
                    case e12://if 2-body --on first like tv1
                        op[0] = 0;
                        op[1] = 1;
                        return two;
                    default:
                        break;

                }
            }else if ( cl == tv2 ){//2-body--on second like tv2
                switch(bl){
                    case id0:
                    case tv1:
                        return one;
                    case tv2:
                        return one;
                    case e12:
                        op[0] = 1;
                        op[1] = 0;
                        return two;
                    default:
                        break;

                }
            }else if ( cl == e12){
                op[0] = 0;
                op[1] = 1;
                return two;
            }
            case three:
//            1, 2, 3,//1///tv1//e12
//            1, 3, 2,//2///e13
//            2, 1, 3,//3///tv2
//            3, 1, 2,//4// e23-1
//            2, 3, 1,//5// e23
//            3, 2, 1//6////tv3
            switch (act) {
                case 1:
                    perm[0] = 0;
                    perm[1] = 1;
                    perm[2] = 2;
                    break;
                case 2:
                    perm[0] = 0;
                    perm[1] = 2;
                    perm[2] = 1;
                    break;
                case 3:
                    perm[0] = 1;
                    perm[1] = 0;
                    perm[2] = 2;
                    break;
                case 4:
                    perm[0] = 2;
                    perm[1] = 0;
                    perm[2] = 1;
                    break;
                case 5:
                    perm[0] = 1;
                    perm[1] = 2;
                    perm[2] = 0;
                    break;
                case 6:
                    perm[0] = 2;
                    perm[1] = 1;
                    perm[2] = 0;
                    break;
                default:
                    break;

                
            }
            if ( cl == tv1){
                    switch(bl){
                        case id0:
                    //1-body
                        case tv1:
                            op[0] = 0;
                            op[1] = 1;
                            op[2] = 2;
                            return one;
                        case tv2:
                            op[0] = 1;
                            op[1] = 0;
                            op[2] = 2;
                            return one;
                        case tv3:
                            op[0] = 2;
                            op[1] = 1;
                            op[2] = 0;
                            return one;

                    //1-body
                        case e12://if 2-body --on first like tv1
                            op[0] = 0;
                            op[1] = 1;
                            op[2] = 2;
                            return two;
                        case e13:
                            op[0] = 0;
                            op[1] = 2;
                            op[2] = 1;
                            return two;
                        case e23:
                            op[0] = 1;
                            op[1] = 2;
                            op[2] = 0;
                            return two;
                        default:
                            break;

                    }
                }else if ( cl == tv2 ){//2-body--on second like tv2
                    switch(bl){
                        case id0:
                        case tv1:
                            return one;
                        case tv2:
                            return one;
                        case e12:
                            op[0] = 1;
                            op[1] = 0;
                            op[2] = 2;
                            return two;
                        case e13:
                            op[0] = 2;
                            op[1] = 0;
                            op[2] = 1;
                            return two;
                        case e23:
                            op[0] = 2;
                            op[1] = 1;
                            op[2] = 0;
                            return two;
                        default:
                            break;

                    }
                }else if ( cl == e12){
                    switch(bl){
                        case e12:
                            op[0] = 0;
                            op[1] = 1;
                            op[2] = 2;
                            return two;
                        case e13:
                            op[0] = 0;
                            op[1] = 2;
                            op[2] = 1;
                            return two;
                        case e23:
                            op[0] = 1;
                            op[1] = 2;
                            op[2] = 0;
                            return two;
                        default:
                            break;


                    }
                }
        case four:
            perm[0] = 0;
            perm[1] = 1;
            perm[2] = 2;
            perm[3] = 3;
            if ( cl == e12){
                switch(bl){
                    case e12:
                        op[0] = 0;
                        op[1] = 1;
                        op[2] = 2;
                        op[3] = 3;
                        return two;
                    case e13:
                        op[0] = 0;
                        op[1] = 2;
                        op[2] = 1;
                        op[3] = 3;
                        return two;
                    case e23:
                        op[0] = 1;
                        op[1] = 2;
                        op[2] = 0;
                        op[3] = 3;
                        return two;
                    case e14:
                        op[0] = 0;
                        op[1] = 2;
                        op[2] = 3;
                        op[3] = 1;
                        return two;
                    case e24:
                        op[0] = 2;
                        op[1] = 0;
                        op[2] = 3;
                        op[3] = 1;
                        return two;
                    case e34:
                        op[0] = 2;
                        op[1] = 3;
                        op[2] = 0;
                        op[3] = 1;
                        return two;
                    default:
                        printf("code!");
                        exit(0);
                        break;
                }
            }
            break;
            default:
            exit(0);
                break;

    }
    return 0;
}


/**
 *Character table elements
 */
double tGetProjection (  bodyType bd , inta irrep , inta op ){
    op--;
    inta nsyp=0 ,msyp=0;
        const static double syp2 [] = {
            1.,1.,
            1.,-1.
    
        };
        const static double syp3 [] = {
            1.,1.,1., 1., 1., 1.,
            1.,-1.,-1.,1.,1.,-1.,
            2.,0.,0.,-1.,-1.,0.
        };
    const static double syp4 [] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., -1., -1., 1., 1., -1., -1., 1., 1.,
        -1., -1., 1., 1., -1., -1., 1., 1., -1., -1., 1., 1., -1., -1., 1.,
        2., 0, 0, -1., -1., 0, 0, 2., -1., 0, 0, -1., -1., 0, 0, -1., 2., 0,
        0, -1., -1., 0, 0, 2., 3., -1., -1., 0, 0, -1., -1., -1., 0, 1., 1.,
        0, 0, 1., -1., 0, -1., 1., 1., 0, 0, -1., 1., -1., 3., 1., 1., 0, 0,
        1., 1., -1., 0, -1., -1., 0, 0, -1., 1., 0, -1., -1., -1., 0, 0, 1.,
        -1., -1.};

        if ( bd == one ){
            return 1;
        }
        else    if ( bd == two ){
            nsyp = 2;
            msyp = 2;
        }
        else if ( bd == three ){
            nsyp = 6;
            msyp = 3;
        }else if ( bd == four ){
            nsyp = 24;
            msyp = 5;
        }
        else {
            printf("unsupported body\n");
            exit(0);
        }
    
    
        if ( irrep <= 0 || irrep > msyp )
        {
            printf("irrep outside bounds\n %d", irrep);
            exit(0);
        }
        if ( op < 0 || op >= nsyp ){
            printf("group action outside bounds\n");
            exit(0);
        }
    
    
        if ( bd == two ){
            return syp2[(irrep-1)*nsyp+op]/nsyp;
        }
        else if ( bd == three ){
            return syp3[(irrep-1)*nsyp+op]/nsyp;
        }
        else if ( bd == four ){
            return syp4[(irrep-1)*nsyp+op]/nsyp;
        }
        return 0.;
}

/**
 *Defining particular vectors in degenerate irreps
 */
double tGetVector(  bodyType bd , inta type , inta perm ){
    perm--;
    type--;
    inta nsyp=0 ,msyp=0;
    const static double syp2 [] = {
        sr2,sr2,
        sr2,-sr2
        
    };
    const static double syp3 [] = {
        /********/
        sr6,sr6,sr6, sr6, sr6, sr6,//A1
        sr6,sr6,sr6,-sr6,-sr6,-sr6,//A2
        hf    , -hf     ,0.      ,hf        , -hf       , 0.,
        hf*sr3,hf*sr3   , -sr3   ,hf*sr3    ,hf*sr3     ,-sr3,
        hf    , -hf     ,0.      ,-hf       , hf        , 0.,
        hf*sr3,hf*sr3   , -sr3   ,-hf*sr3   ,-hf*sr3    ,sr3,
    };
    
    const static double syp4[] = {
        /*1: A1*/
        sr6*hf,sr6*hf,sr6*hf,sr6*hf,    sr6*hf,sr6*hf,sr6*hf,sr6*hf,    sr6*hf,sr6*hf,sr6*hf,sr6*hf,
        sr6*hf,sr6*hf,sr6*hf,sr6*hf,    sr6*hf,sr6*hf,sr6*hf,sr6*hf,    sr6*hf,sr6*hf,sr6*hf,sr6*hf,
        
        /*2: A2*/
        sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,    sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,    sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,
        sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,    sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,    sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,
        
        /*3: E-1*/
        sr6,    0.,0.,-sr6*hf,      -sr6*hf,0.,0., sr6,       -sr6*hf,0.,0.,-sr6*hf,
        -sr6*hf,0.,0.,-sr6*hf,       sr6,   0.,0.,-sr6*hf,    -sr6*hf,0.,0.,sr6,
        
        /*4: E-2*/
        0.,sr6,-sr6*hf,0.,          0.,-sr6*hf, sr6,0.,    0.,-sr6*hf,-sr6*hf,0.,
        0.,-sr6*hf,-sr6*hf,0.,      0.,sr6, -sr6*hf,0.,    0.,-sr6*hf,sr6,0.,
        
        /*5: E-3*/
        0.,0.,sr2*hf,0.,             0.,-sr2*hf,0.,0.,       0.,-sr2*hf ,sr2*hf,0.,
        0.,sr2*hf ,-sr2*hf,0.,        0.,0.,-sr2*hf,0,        0.,sr2*hf,0.,0.,
        
        /*6: E-4*/
        0.,0.,0.,sr2*hf,             -sr2*hf,0.,0.,0.,         -sr2*hf,0.,0.,sr2*hf,
        sr2*hf,0.,0.,-sr2*hf,       0.,0.,0.,-sr2*hf,           sr2*hf,0.,0.,0.,
        
        /*7: T1-1*/
        hf*sr2/sr3,-hf*sr6,-hf*sr6,0.,      0.,-hf*sr6,-hf*sr6,-hf*sr6,     0.,hf*sr6,hf*sr6,0.,
        0., hf*sr6,-hf*sr6,0.,              -hf*sr6, hf*sr6, hf*sr6,0.,     0.,-hf*sr6,hf*sr6,-hf*sr6,
        
        /*8: T1-2*/
        0. , sr3,-hf*hf*hf*sr3,-hf*hf*hf/sr3,                               -hf*hf*hf/sr3,-hf*hf*hf*sr3, -hf*sr3, -hf*sr3,
        hf*hf*hf/sr3,hf*hf*hf*sr3,hf*hf*hf*sr3,hf*hf*hf/sr3,                  hf*hf*hf/sr3,hf*hf*hf*sr3,-hf*hf*hf*sr3,-hf*hf*hf/sr3,
        hf*hf*sr3,-hf*hf*sr3,hf*hf*hf*sr3,hf*hf*hf/sr3,                       -hf*hf*hf/sr3,-hf*hf*hf*sr3,-hf*hf*sr3,  hf*hf*sr3,
        
        /*9: T1-3*/
        0.,0.,hf*hf*hf/sr3/sr7,-hf*hf*hf*sr7/sr3/sr3/sr3,
        -hf*hf*hf*sr7/sr3/sr3/sr3,-hf*hf*hf*sr7/sr3,-hf*sr3*sr7, hf*sr3*sr7,
        -sr3*hf*hf*hf/sr7,1./sr3*hf*hf*hf*sr7,-hf*hf*hf*sr3*sr7/sr5/sr5, hf*hf*hf*sr7/sr3/sr3/sr3,
        -sr3*hf*hf*hf/sr7,-hf*hf*hf*sr3*sr7/sr5/sr5,-hf*hf*hf/sr3*sr7,sr3*hf*hf*hf/sr7,
        hf*hf*sr7/sr3, hf*hf*sr3*sr7, hf*hf*hf*sr7/sr3,  hf*hf*hf*sr7/sr3/sr3/sr3,
        hf*hf*hf*sr3/sr7,-11.*hf*hf*hf*sr3*sr7, hf*hf*sr3*sr7,  -hf*hf*sr3*sr7/sr5/sr5,
        
        
        /*10: T1-4*/
        0.,0.,0.,hf*sr2*sr7/sr5/sr3,
        -sr5*sr2*sr7/sr3, -hf*sr5*sr2*sr7/sr3/sr3/sr3, sr5*sr2*sr3*sr7,-sr5*sr2*sr3*sr7,
        0.,-hf*sr6*sr7/sr5, sr5*sr2*sr7/sr3, -hf * sr7*sr3*sr5*sr2,
        -hf*sr5*sr6/sr7,-sr3*sr7*sr5/hf/sr2,hf*sr6*sr7/sr5,0.,
        sr3*sr7*sr5/sr2, sr7*sr5*sr2/sr3,   hf*sr7*sr5*sr2/sr3/sr3/sr3, sr7*sr5*sr2/sr3,
        -hf*sr5*sr6/sr7,sr5*sr2*sr3*sr7,-sr3*sr7*sr5/hf/sr2,-sr5*sr2*sr3*sr7,
        
        /*11: T1-5*/
        0.,0.,0.,0.,
        hf/sr3/sr3*sr2*sr5,-hf/sr3/sr3*sr2*sr5,sr3*sr3*sr2*sr5,-sr3*sr3*sr2*sr5,
        -sr6*sr6*sr2/sr5,sr6*sr6*sr2/sr5,-sr5*sr3*sr3/sr2,sr5*sr3*sr3/sr2,
        -sr2*sr5*sr3*sr3,sr2*sr5*sr3*sr3,sr6*sr6*sr2/sr5,-sr6*sr6*sr2/sr5,
        sr5*sr3*sr3/sr2,-sr5*sr3*sr3/sr2,-sr2*sr5*sr6*sr6,sr6*sr6*sr2*sr5,
        -sr2*sr5*sr3*sr3,sr2*sr5*sr3*sr3,sr2*sr5*sr3*sr3,-sr2*sr5*sr3*sr3,
        
        /*12: T1-6*/
        0.,0.,0.,0.,
        0.,0.,sr3*sr3/sr2,-sr3*sr3/sr2,
        -sr6*sr6*sr2, sr6*sr6*sr2,sr6*sr6*sr2,-sr6*sr6*sr2,
        -sr6*sr6*sr2,sr6*sr6*sr2,-sr3*sr3*sr2,sr3*sr3*sr2,
        sr3*sr3*sr2,-sr3*sr3*sr2,sr6*sr6*sr2,-sr6*sr6*sr2,
        sr3*sr3*sr2,-sr3*sr3*sr2,-sr3*sr3*sr2,sr3*sr3*sr2,
        
        /*13: T1-7*/
        0.,0.,0.,0.,
        0.,0.,0.,0.,
        hf*sr6/sr5,-hf*sr5*sr6,-hf*sr6/sr5,hf*sr6*sr5,
        -hf*sr2*sr5/sr3,hf*sr2*sr5/sr3,-sr6*sr5,-sr6*sr5,
        sr6*sr5,sr6*sr5,hf*sr2*sr5/sr3,-hf*sr2*sr5/sr3,
        sr5*sr6,sr5*sr6,-sr5*sr6,-sr5*sr6,
        
        /*14: T1-8*/
        0.,0.,0.,0.,
        0.,0.,0.,0.,
        0.,sr5,0,-sr5,
        hf*sr5,-hf*sr5,-hf*sr5,-hf*sr5,
        hf*sr5,hf*sr5,-hf*sr5,hf*sr5,
        hf*sr5,hf*sr5,-hf*sr5,-hf*sr5,
        
        /*15: T1-9*/
        0.,0.,0.,0.,
        0.,0.,0.,0.,
        0.,0.,0.,0.,
        hf*sr3,-hf*sr3,-hf*sr3,hf*sr3,
        hf*sr3,-hf*sr3,hf*sr3,-hf*sr3,
        -hf*sr3,hf*sr3,hf*sr3,-hf*sr3,
        
        /*16: T2-1*/
        hf*sr2/sr3,hf*sr2*sr3,hf*sr2*sr3,0.,
        0.,hf*sr2*sr3,hf*sr2*sr3,-hf*sr2*sr3,
        0.,-hf*sr2*sr3,-hf*sr2*sr3,0.,
        0.,-hf*sr2*sr3,hf*sr2*sr3,0.,
        -hf*sr2*sr3,-hf*sr2*sr3,-hf*sr2*sr3,0.,
        0.,hf*sr2*sr3,-hf*sr2*sr3,-hf*sr2*sr3,
        
        /*17: T2-2*/
        0.,sr3,-hf*hf*hf*sr3,hf*hf*hf/sr3,
        hf*hf*hf/sr3,-hf*hf*hf*sr3,-hf*sr3,hf*sr3,
        -hf*hf*hf/sr3,hf*hf*hf*sr3,hf*hf*hf*sr3,-hf*hf*hf/sr3,
        -hf*hf*hf/sr3,hf*hf*hf*sr3,-hf*hf*hf*sr3,hf*hf*hf/sr3,
        -hf*hf*sr3,-hf*hf*sr3,hf*hf*hf*sr3,-hf*hf*hf/sr3,
        hf*hf*hf/sr3,-hf*hf*hf*sr3,-hf*hf*sr3,-hf*hf*sr3,
        
        /*18: T2-3*/
        0.,0.,hf*hf*hf/sr3/sr7,hf*hf*hf*sr7/sr3/sr3/sr3,
        hf*hf*hf*sr7/sr3/sr3/sr3,-hf*hf*hf*sr7/sr3,-hf*sr3*sr7,-hf*sr3*sr7,
        sr3/sr7*hf*hf*hf, sr7/sr3*hf*hf*hf, -hf*hf*hf*sr3*sr7/sr5/sr5, -hf*hf*hf*sr7/sr3/sr3/sr3,
        hf*hf*hf*sr3/sr7,-hf*hf*hf*sr3*sr7/sr5/sr5,-hf*hf*hf*sr7/sr3,-hf*hf*hf*sr3/sr7,
        -hf*hf*sr7/sr3, hf*hf*sr3*sr7,  hf*hf*hf*sr7/sr3, - hf*hf*hf*sr7/sr3/sr3/sr3,
        -sr3/sr7*hf*hf*hf, -11.*hf*hf*hf*sr3*sr7, hf*hf*sr3*sr7, hf*hf*sr3*sr7/sr5/sr5,
        
        
        /*19: T2-4*/
        0.,0.,0.,hf*sr2*sr7/sr3/sr5,
        -sr7*sr5*sr2/sr3,hf*sr7*sr5*sr2/sr3/sr3/sr3,    - sr7*sr3*sr5*sr2,- sr7*sr3*sr5*sr2,
        0.,hf/sr5*sr6*sr7,  -sr7*sr5*sr2/sr3, -hf*sr7*sr3*sr5*sr2,
        -hf*sr5*sr6/sr7,sr7*sr3*sr5/sr2/hf,-hf*sr6*sr7/sr5,0.,
        sr7*sr3*sr5/sr2,-sr7*sr5*sr2/sr3,-hf*sr7*sr5*sr2/sr3/sr3/sr3,sr7*sr5*sr2/sr3,
        -hf*sr5*sr6/sr7,-sr7*sr3*sr5*sr2,sr7*sr3*sr5/sr2/hf,-sr7*sr3*sr5*sr2,
        
        /*20: T2-5*/
        0.,0.,0.,0.,
        sr5*sr2*hf/sr3/sr3,sr5*sr2*hf/sr3/sr3,-sr2*sr5*sr3*sr3,-sr2*sr5*sr3*sr3,
        -sr6*sr6*sr2/sr5,-sr6*sr6*sr2/sr5,sr3*sr3*sr5/sr2,sr3*sr3*sr5/sr2,
        -sr3*sr3*sr2*sr5, -sr3*sr3*sr2*sr5,-sr6*sr6*sr2/sr5,-sr6*sr6*sr2/sr5,
        sr3*sr3*sr5/sr2,sr3*sr3*sr5/sr2,sr6*sr6*sr2*sr5,sr6*sr6*sr2*sr5,
        -sr3*sr3*sr2*sr5, -sr3*sr3*sr2*sr5,-sr3*sr3*sr2*sr5, -sr3*sr3*sr2*sr5,
        
        
        /*21: T2-6*/
        0.,0.,0.,0.,
        0.,0.,sr3*sr3/sr2,sr3*sr3/sr2,
        sr6*sr6*sr2,sr6*sr6*sr2,sr6*sr6*sr2,sr6*sr6*sr2,
        sr6*sr6*sr2,sr6*sr6*sr2,-sr3*sr3*sr2,-sr3*sr3*sr2,
        -sr3*sr3*sr2,-sr3*sr3*sr2,sr6*sr6*sr2,sr6*sr6*sr2,
        -sr3*sr3*sr2,-sr3*sr3*sr2,-sr3*sr3*sr2,-sr3*sr3*sr2,
        
        /*22: T2-7*/
        0.,0.,0.,0.,
        0.,0.,0.,0.,
        hf*sr6/sr5,hf*sr6*sr5,hf*sr6/sr5,hf*sr5*sr6,
        -sr2*sr5/sr3*hf, -sr2*sr5/sr3*hf,sr5*sr6,-sr5*sr6,
        sr5*sr6,-sr5*sr6, -sr2*sr5/sr3*hf, -sr2*sr5/sr3*hf,
        sr5*sr6,-sr5*sr6,sr5*sr6,-sr5*sr6,
        
        /*23: T2-8*/
        0.,0.,0.,0.,
        0.,0.,0.,0.,
        0.,sr5,0.,sr5,
        -hf*sr5,-hf*sr5,-hf*sr5,hf*sr5,
        -hf*sr5,hf*sr5,-hf*sr5,-hf*sr5,
        -hf*sr5,hf*sr5,-hf*sr5,hf*sr5,
        
        /*24: T2-9*/
        0.,0.,0.,0.,
        0.,0.,0.,0.,
        0.,0.,0.,0.,
        hf*sr3,hf*sr3,hf*sr3,hf*sr3,
        hf*sr3,hf*sr3,-hf*sr3,-hf*sr3,
        -hf*sr3,-hf*sr3,-hf*sr3,-hf*sr3
        
        //////////////////////////////////
        /////////////////////////////////
        ////////////////////////////////
        ///////////////////////////////
        //////////////////////////////
        
        
    };
    
    
    if ( bd == one )
    {
        nsyp = 1;
        msyp = 1;
    }
    else  if ( bd == two ){
        nsyp = 2;
        msyp = 2;
    }
    else if ( bd == three ){
        nsyp = 6;
        msyp = 6;
    }else if ( bd == four ){
        nsyp = 24;
        msyp = 24;
    }
    else {
        printf("bod\n");
        exit(0);
    }
    
    
    if ( type < 0 || type >= msyp )
    {
        printf("he\n %d", type);
        exit(0);
    }
    if ( perm < 0 || perm >= nsyp ){
        printf("hm\n");
        exit(0);
    }
    
    
    if ( bd == two ){
        return syp2[type*nsyp+perm]*sr2;
    }
    else if ( bd == three ){
        return syp3[type*nsyp+perm]*sr6;
    }
    else if ( bd == four ){
        return syp4[type*nsyp+perm]*sr6*hf;
    }
    return 0.;
};

/**
 *Classification sub-routine
 *
 */
inta tClassifyComponents(   sinc_label  f1 , double * up, double * entropy){
    
    if ( bodies(f1,eigenVectors ) == one ){
        return 1;
    }
    double entr,sum;
    inta nGroup=0,xt,irrep;
    
    if ( bodies(f1, eigenVectors ) == two ){
        nGroup = 2;
    }
    else if ( bodies ( f1, eigenVectors ) == three ){
        nGroup = 3;
    }
    else if ( bodies (f1, eigenVectors ) == four  ){
        nGroup = 5;
    }
    xt=0;
    entr = 0.;
    sum = 0.;
    for ( irrep = 0 ; irrep <= nGroup ; irrep++ ){
        sum += fabs(up[irrep]);
        if ( fabs(up[irrep])> fabs(up[xt]))
            xt = irrep;
    }
    for ( irrep = 0 ; irrep <= nGroup ; irrep++ ){
        if ( fabs(up[irrep]) > 1e-6 ){
            entr += -(fabs(up[irrep])/sum)*log(fabs(up[irrep])/sum);
        }
    }
    *entropy = entr;
//    printf("\n");

    return xt;
}

/**
 *Classification routine
 *
*/
inta tClassify( sinc_label  f1 ,   division label){
    double up[48],entropy;
    if ( bodies(f1, label ) == one ){
        f1.name[label].value.symmetry = 0;
        f1.name[label].value.value2 = 0;
        
        return 0;
        
    }
    inta i,irrep;
    for ( i = 0; i < 48 ; i++)
        up[i] = 0.;
    tTabulateInnerProjection( f1, label, up);
    irrep =  tClassifyComponents(f1, up,&entropy);
    f1.name[label].value.value2 =entropy;

    return irrep;
}



/**
 *Irrep projection routine
 *
 *Center of 'filter'
 *
*/
inta tBuildIrr ( inta rank,   sinc_label  f1, inta meta , division origin, inta ospin,division targ , inta tspin){
    inta irrep;
    if ( ! allowQ(f1.rt,blockPermutationsblock)){
        printf("blockPermutationsblock Allow!\n");
        fflush(stdout);
        exit(0);
    }

    if ( meta == 0 || bodies(f1,origin ) == one ){
        tEqua(f1, targ, tspin, origin, ospin);
        return 0;
    }
    if (! CanonicalRank(f1, origin, ospin ))
        return 0;
    
    inta o,space,i,nPerm=tPerms(bodies(f1,origin));
    f1.name[permutationVector].Current[rank] = 1;
    f1.name[permutation2Vector].Current[rank] = 1;

    if ( CanonicalRank(f1, origin, ospin)*nPerm > part(f1,targ)){
        printf("Irr\n");
        printf("%d %d %d\n", CanonicalRank(f1, origin, ospin), part(f1, origin), part(f1, targ));
        exit(0);
    }
    
        for ( i = 1; i <= nPerm ; i++){
            irrep = meta;
            f1.name[permutation2Vector].Current[rank] = 0;
            for ( o = 0; o < CanonicalRank(f1, origin, ospin);o++){
                for ( space = 0; space < SPACE ; space++)
                    if ( f1.canon[space].body != nada )
                        xsEqu(1., space, f1, permutationVector, 0, rank, space, f1, origin, o, ospin);
                f1.name[permutation2Vector].Current[rank] = 0;
                tPermute(rank,f1, i, permutationVector, rank, permutation2Vector, rank);
                tScaleOne(f1, permutation2Vector, rank, tGetProjection(bodies(f1, origin), irrep, i));
                tAddTw(f1, targ, tspin, permutation2Vector, rank);
            }
        }
    
    return 0;
}

/**
 * group action of leftChar by dimension
 */
inta tPermuteOne(inta rank,   sinc_label  f1, inta dim, inta leftChar ,   division left, inta l, inta lspin,   division equals, inta e, inta espin){
    inta at1[] = {1};
    inta at2[] = {
        1, 2,
        2, 1
    };
    inta at3[] = {
        1, 2, 3,//1///tv1//e12
        1, 3, 2,//2///e13
        2, 1, 3,//3///tv2
        3, 1, 2,//4// e23-1
        2, 3, 1,//5// e23
        3, 2, 1//6////tv3
    };
    inta at4[] = {
        1, 2, 3, 4,
        1, 2, 4, 3,
        1, 3, 2, 4,
        1, 4, 2, 3,
        1, 3, 4, 2,
        1, 4, 3, 2,
        2, 1, 3, 4,
        2, 1, 4, 3,
        3, 1, 2, 4,
        4, 1, 2, 3,
        3, 1, 4, 2,
        4, 1, 3, 2,
        2, 3, 1, 4,
        2, 4, 1, 3,
        3, 2, 1, 4,
        4, 2, 1, 3,
        3, 4, 1, 2,
        4, 3, 1, 2,
        2, 3, 4, 1,
        2, 4, 3, 1,
        3, 2, 4, 1,
        4, 2, 3, 1,
        3, 4, 2, 1,
        4, 3, 2, 1};
    inta at5[] = {
        1, 2, 3, 4, 5,//tv1//e12
        1, 2, 3, 5, 4,
        1, 2, 4, 3, 5,
        1, 2, 5, 3, 4,
        1, 2, 4, 5, 3,
        1, 2, 5, 4, 3,
        1, 3, 2, 4, 5,//e13
        1, 3, 2, 5, 4,
        1, 4, 2, 3, 5,
        1, 5, 2, 3, 4,
        1, 4, 2, 5, 3,
        1, 5, 2, 4, 3,
        1, 3, 4, 2, 5,
        1, 3, 5, 2, 4,
        1, 4, 3, 2, 5,//e14
        1, 5, 3, 2, 4,
        1, 4, 5, 2, 3,
        1, 5, 4, 2, 3,
        1, 3, 4, 5, 2,
        1, 3, 5, 4, 2,
        1, 4, 3, 5, 2,
        1, 5, 3, 4, 2,//e15
        1, 4, 5, 3, 2,
        1, 5, 4, 3, 2,
        2, 1, 3, 4, 5,//tv2
        2, 1, 3, 5, 4,
        2, 1, 4, 3, 5,
        2, 1, 5, 3, 4,
        2, 1, 4, 5, 3,
        2, 1, 5, 4, 3,
        3, 1, 2, 4, 5,
        3, 1, 2, 5, 4,
        4, 1, 2, 3, 5,
        5, 1, 2, 3, 4,//e23 : -1
        4, 1, 2, 5, 3,
        5, 1, 2, 4, 3,
        3, 1, 4, 2, 5,
        3, 1, 5, 2, 4,//e24 : -1
        4, 1, 3, 2, 5,
        5, 1, 3, 2, 4,
        4, 1, 5, 2, 3,
        5, 1, 4, 2, 3,
        3, 1, 4, 5, 2,
        3, 1, 5, 4, 2,
        4, 1, 3, 5, 2,
        5, 1, 3, 4, 2,//e25 : -1
        4, 1, 5, 3, 2,
        5, 1, 4, 3, 2,
        2, 3, 1, 4, 5,
        2, 3, 1, 5, 4,
        2, 4, 1, 3, 5,
        2, 5, 1, 3, 4,
        2, 4, 1, 5, 3,//e24 : 1 //3, 1, 5, 2, 4
        2, 5, 1, 4, 3,
        3, 2, 1, 4, 5,//tv3
        3, 2, 1, 5, 4,
        4, 2, 1, 3, 5,
        5, 2, 1, 3, 4,
        4, 2, 1, 5, 3,
        5, 2, 1, 4, 3,
        3, 4, 1, 2, 5,
        3, 5, 1, 2, 4,//e34
        4, 3, 1, 2, 5,
        5, 3, 1, 2, 4,
        4, 5, 1, 2, 3,
        5, 4, 1, 2, 3,
        3, 4, 1, 5, 2, //e34 : 1 // 3, 5, 1, 2, 4
        3, 5, 1, 4, 2, //e35 : 1 // 3, 5, 1, 4, 2
        4, 3, 1, 5, 2,
        5, 3, 1, 4, 2,
        4, 5, 1, 3, 2, //e45 : 1 // 3, 5, 4, 1, 2
        5, 4, 1, 3, 2,
        2, 3, 4, 1, 5,
        2, 3, 5, 1, 4,
        2, 4, 3, 1, 5,
        2, 5, 3, 1, 4,
        2, 4, 5, 1, 3,
        2, 5, 4, 1, 3,
        3, 2, 4, 1, 5,
        3, 2, 5, 1, 4,
        4, 2, 3, 1, 5,//tv4
        5, 2, 3, 1, 4,
        4, 2, 5, 1, 3,
        5, 2, 4, 1, 3,
        3, 4, 2, 1, 5,
        3, 5, 2, 1, 4,
        4, 3, 2, 1, 5,
        5, 3, 2, 1, 4,
        4, 5, 2, 1, 3,
        5, 4, 2, 1, 3,
        3, 4, 5, 1, 2,
        3, 5, 4, 1, 2,//e45 : -1
        4, 3, 5, 1, 2,
        5, 3, 4, 1, 2,
        4, 5, 3, 1, 2,
        5, 4, 3, 1, 2,
        2, 3, 4, 5, 1,//e23:1
        2, 3, 5, 4, 1,
        2, 4, 3, 5, 1,
        2, 5, 3, 4, 1,//e25 : 5, 1, 3, 4, 2
        2, 4, 5, 3, 1,
        2, 5, 4, 3, 1,
        3, 2, 4, 5, 1,
        3, 2, 5, 4, 1,
        4, 2, 3, 5, 1,
        5, 2, 3, 4, 1,//tv5
        4, 2, 5, 3, 1,
        5, 2, 4, 3, 1,
        3, 4, 2, 5, 1,
        3, 5, 2, 4, 1,
        4, 3, 2, 5, 1,
        5, 3, 2, 4, 1,
        4, 5, 2, 3, 1,
        5, 4, 2, 3, 1,
        3, 4, 5, 2, 1,
        3, 5, 4, 2, 1,
        4, 3, 5, 2, 1,
        5, 3, 4, 2, 1,
        4, 5, 3, 2, 1,
        5, 4, 3, 2, 1
    };
    inta at6[] = {1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 6, 5, 1, 2, 3, 5, 4, 6, 1, 2, 3, 6, 4,
        5, 1, 2, 3, 5, 6, 4, 1, 2, 3, 6, 5, 4, 1, 2, 4, 3, 5, 6, 1, 2, 4, 3,
        6, 5, 1, 2, 5, 3, 4, 6, 1, 2, 6, 3, 4, 5, 1, 2, 5, 3, 6, 4, 1, 2, 6,
        3, 5, 4, 1, 2, 4, 5, 3, 6, 1, 2, 4, 6, 3, 5, 1, 2, 5, 4, 3, 6, 1, 2,
        6, 4, 3, 5, 1, 2, 5, 6, 3, 4, 1, 2, 6, 5, 3, 4, 1, 2, 4, 5, 6, 3, 1,
        2, 4, 6, 5, 3, 1, 2, 5, 4, 6, 3, 1, 2, 6, 4, 5, 3, 1, 2, 5, 6, 4, 3,
        1, 2, 6, 5, 4, 3, 1, 3, 2, 4, 5, 6, 1, 3, 2, 4, 6, 5, 1, 3, 2, 5, 4,
        6, 1, 3, 2, 6, 4, 5, 1, 3, 2, 5, 6, 4, 1, 3, 2, 6, 5, 4, 1, 4, 2, 3,
        5, 6, 1, 4, 2, 3, 6, 5, 1, 5, 2, 3, 4, 6, 1, 6, 2, 3, 4, 5, 1, 5, 2,
        3, 6, 4, 1, 6, 2, 3, 5, 4, 1, 4, 2, 5, 3, 6, 1, 4, 2, 6, 3, 5, 1, 5,
        2, 4, 3, 6, 1, 6, 2, 4, 3, 5, 1, 5, 2, 6, 3, 4, 1, 6, 2, 5, 3, 4, 1,
        4, 2, 5, 6, 3, 1, 4, 2, 6, 5, 3, 1, 5, 2, 4, 6, 3, 1, 6, 2, 4, 5, 3,
        1, 5, 2, 6, 4, 3, 1, 6, 2, 5, 4, 3, 1, 3, 4, 2, 5, 6, 1, 3, 4, 2, 6,
        5, 1, 3, 5, 2, 4, 6, 1, 3, 6, 2, 4, 5, 1, 3, 5, 2, 6, 4, 1, 3, 6, 2,
        5, 4, 1, 4, 3, 2, 5, 6, 1, 4, 3, 2, 6, 5, 1, 5, 3, 2, 4, 6, 1, 6, 3,
        2, 4, 5, 1, 5, 3, 2, 6, 4, 1, 6, 3, 2, 5, 4, 1, 4, 5, 2, 3, 6, 1, 4,
        6, 2, 3, 5, 1, 5, 4, 2, 3, 6, 1, 6, 4, 2, 3, 5, 1, 5, 6, 2, 3, 4, 1,
        6, 5, 2, 3, 4, 1, 4, 5, 2, 6, 3, 1, 4, 6, 2, 5, 3, 1, 5, 4, 2, 6, 3,
        1, 6, 4, 2, 5, 3, 1, 5, 6, 2, 4, 3, 1, 6, 5, 2, 4, 3, 1, 3, 4, 5, 2,
        6, 1, 3, 4, 6, 2, 5, 1, 3, 5, 4, 2, 6, 1, 3, 6, 4, 2, 5, 1, 3, 5, 6,
        2, 4, 1, 3, 6, 5, 2, 4, 1, 4, 3, 5, 2, 6, 1, 4, 3, 6, 2, 5, 1, 5, 3,
        4, 2, 6, 1, 6, 3, 4, 2, 5, 1, 5, 3, 6, 2, 4, 1, 6, 3, 5, 2, 4, 1, 4,
        5, 3, 2, 6, 1, 4, 6, 3, 2, 5, 1, 5, 4, 3, 2, 6, 1, 6, 4, 3, 2, 5, 1,
        5, 6, 3, 2, 4, 1, 6, 5, 3, 2, 4, 1, 4, 5, 6, 2, 3, 1, 4, 6, 5, 2, 3,
        1, 5, 4, 6, 2, 3, 1, 6, 4, 5, 2, 3, 1, 5, 6, 4, 2, 3, 1, 6, 5, 4, 2,
        3, 1, 3, 4, 5, 6, 2, 1, 3, 4, 6, 5, 2, 1, 3, 5, 4, 6, 2, 1, 3, 6, 4,
        5, 2, 1, 3, 5, 6, 4, 2, 1, 3, 6, 5, 4, 2, 1, 4, 3, 5, 6, 2, 1, 4, 3,
        6, 5, 2, 1, 5, 3, 4, 6, 2, 1, 6, 3, 4, 5, 2, 1, 5, 3, 6, 4, 2, 1, 6,
        3, 5, 4, 2, 1, 4, 5, 3, 6, 2, 1, 4, 6, 3, 5, 2, 1, 5, 4, 3, 6, 2, 1,
        6, 4, 3, 5, 2, 1, 5, 6, 3, 4, 2, 1, 6, 5, 3, 4, 2, 1, 4, 5, 6, 3, 2,
        1, 4, 6, 5, 3, 2, 1, 5, 4, 6, 3, 2, 1, 6, 4, 5, 3, 2, 1, 5, 6, 4, 3,
        2, 1, 6, 5, 4, 3, 2, 2, 1, 3, 4, 5, 6, 2, 1, 3, 4, 6, 5, 2, 1, 3, 5,
        4, 6, 2, 1, 3, 6, 4, 5, 2, 1, 3, 5, 6, 4, 2, 1, 3, 6, 5, 4, 2, 1, 4,
        3, 5, 6, 2, 1, 4, 3, 6, 5, 2, 1, 5, 3, 4, 6, 2, 1, 6, 3, 4, 5, 2, 1,
        5, 3, 6, 4, 2, 1, 6, 3, 5, 4, 2, 1, 4, 5, 3, 6, 2, 1, 4, 6, 3, 5, 2,
        1, 5, 4, 3, 6, 2, 1, 6, 4, 3, 5, 2, 1, 5, 6, 3, 4, 2, 1, 6, 5, 3, 4,
        2, 1, 4, 5, 6, 3, 2, 1, 4, 6, 5, 3, 2, 1, 5, 4, 6, 3, 2, 1, 6, 4, 5,
        3, 2, 1, 5, 6, 4, 3, 2, 1, 6, 5, 4, 3, 3, 1, 2, 4, 5, 6, 3, 1, 2, 4,
        6, 5, 3, 1, 2, 5, 4, 6, 3, 1, 2, 6, 4, 5, 3, 1, 2, 5, 6, 4, 3, 1, 2,
        6, 5, 4, 4, 1, 2, 3, 5, 6, 4, 1, 2, 3, 6, 5, 5, 1, 2, 3, 4, 6, 6, 1,
        2, 3, 4, 5, 5, 1, 2, 3, 6, 4, 6, 1, 2, 3, 5, 4, 4, 1, 2, 5, 3, 6, 4,
        1, 2, 6, 3, 5, 5, 1, 2, 4, 3, 6, 6, 1, 2, 4, 3, 5, 5, 1, 2, 6, 3, 4,
        6, 1, 2, 5, 3, 4, 4, 1, 2, 5, 6, 3, 4, 1, 2, 6, 5, 3, 5, 1, 2, 4, 6,
        3, 6, 1, 2, 4, 5, 3, 5, 1, 2, 6, 4, 3, 6, 1, 2, 5, 4, 3, 3, 1, 4, 2,
        5, 6, 3, 1, 4, 2, 6, 5, 3, 1, 5, 2, 4, 6, 3, 1, 6, 2, 4, 5, 3, 1, 5,
        2, 6, 4, 3, 1, 6, 2, 5, 4, 4, 1, 3, 2, 5, 6, 4, 1, 3, 2, 6, 5, 5, 1,
        3, 2, 4, 6, 6, 1, 3, 2, 4, 5, 5, 1, 3, 2, 6, 4, 6, 1, 3, 2, 5, 4, 4,
        1, 5, 2, 3, 6, 4, 1, 6, 2, 3, 5, 5, 1, 4, 2, 3, 6, 6, 1, 4, 2, 3, 5,
        5, 1, 6, 2, 3, 4, 6, 1, 5, 2, 3, 4, 4, 1, 5, 2, 6, 3, 4, 1, 6, 2, 5,
        3, 5, 1, 4, 2, 6, 3, 6, 1, 4, 2, 5, 3, 5, 1, 6, 2, 4, 3, 6, 1, 5, 2,
        4, 3, 3, 1, 4, 5, 2, 6, 3, 1, 4, 6, 2, 5, 3, 1, 5, 4, 2, 6, 3, 1, 6,
        4, 2, 5, 3, 1, 5, 6, 2, 4, 3, 1, 6, 5, 2, 4, 4, 1, 3, 5, 2, 6, 4, 1,
        3, 6, 2, 5, 5, 1, 3, 4, 2, 6, 6, 1, 3, 4, 2, 5, 5, 1, 3, 6, 2, 4, 6,
        1, 3, 5, 2, 4, 4, 1, 5, 3, 2, 6, 4, 1, 6, 3, 2, 5, 5, 1, 4, 3, 2, 6,
        6, 1, 4, 3, 2, 5, 5, 1, 6, 3, 2, 4, 6, 1, 5, 3, 2, 4, 4, 1, 5, 6, 2,
        3, 4, 1, 6, 5, 2, 3, 5, 1, 4, 6, 2, 3, 6, 1, 4, 5, 2, 3, 5, 1, 6, 4,
        2, 3, 6, 1, 5, 4, 2, 3, 3, 1, 4, 5, 6, 2, 3, 1, 4, 6, 5, 2, 3, 1, 5,
        4, 6, 2, 3, 1, 6, 4, 5, 2, 3, 1, 5, 6, 4, 2, 3, 1, 6, 5, 4, 2, 4, 1,
        3, 5, 6, 2, 4, 1, 3, 6, 5, 2, 5, 1, 3, 4, 6, 2, 6, 1, 3, 4, 5, 2, 5,
        1, 3, 6, 4, 2, 6, 1, 3, 5, 4, 2, 4, 1, 5, 3, 6, 2, 4, 1, 6, 3, 5, 2,
        5, 1, 4, 3, 6, 2, 6, 1, 4, 3, 5, 2, 5, 1, 6, 3, 4, 2, 6, 1, 5, 3, 4,
        2, 4, 1, 5, 6, 3, 2, 4, 1, 6, 5, 3, 2, 5, 1, 4, 6, 3, 2, 6, 1, 4, 5,
        3, 2, 5, 1, 6, 4, 3, 2, 6, 1, 5, 4, 3, 2, 2, 3, 1, 4, 5, 6, 2, 3, 1,
        4, 6, 5, 2, 3, 1, 5, 4, 6, 2, 3, 1, 6, 4, 5, 2, 3, 1, 5, 6, 4, 2, 3,
        1, 6, 5, 4, 2, 4, 1, 3, 5, 6, 2, 4, 1, 3, 6, 5, 2, 5, 1, 3, 4, 6, 2,
        6, 1, 3, 4, 5, 2, 5, 1, 3, 6, 4, 2, 6, 1, 3, 5, 4, 2, 4, 1, 5, 3, 6,
        2, 4, 1, 6, 3, 5, 2, 5, 1, 4, 3, 6, 2, 6, 1, 4, 3, 5, 2, 5, 1, 6, 3,
        4, 2, 6, 1, 5, 3, 4, 2, 4, 1, 5, 6, 3, 2, 4, 1, 6, 5, 3, 2, 5, 1, 4,
        6, 3, 2, 6, 1, 4, 5, 3, 2, 5, 1, 6, 4, 3, 2, 6, 1, 5, 4, 3, 3, 2, 1,
        4, 5, 6, 3, 2, 1, 4, 6, 5, 3, 2, 1, 5, 4, 6, 3, 2, 1, 6, 4, 5, 3, 2,
        1, 5, 6, 4, 3, 2, 1, 6, 5, 4, 4, 2, 1, 3, 5, 6, 4, 2, 1, 3, 6, 5, 5,
        2, 1, 3, 4, 6, 6, 2, 1, 3, 4, 5, 5, 2, 1, 3, 6, 4, 6, 2, 1, 3, 5, 4,
        4, 2, 1, 5, 3, 6, 4, 2, 1, 6, 3, 5, 5, 2, 1, 4, 3, 6, 6, 2, 1, 4, 3,
        5, 5, 2, 1, 6, 3, 4, 6, 2, 1, 5, 3, 4, 4, 2, 1, 5, 6, 3, 4, 2, 1, 6,
        5, 3, 5, 2, 1, 4, 6, 3, 6, 2, 1, 4, 5, 3, 5, 2, 1, 6, 4, 3, 6, 2, 1,
        5, 4, 3, 3, 4, 1, 2, 5, 6, 3, 4, 1, 2, 6, 5, 3, 5, 1, 2, 4, 6, 3, 6,
        1, 2, 4, 5, 3, 5, 1, 2, 6, 4, 3, 6, 1, 2, 5, 4, 4, 3, 1, 2, 5, 6, 4,
        3, 1, 2, 6, 5, 5, 3, 1, 2, 4, 6, 6, 3, 1, 2, 4, 5, 5, 3, 1, 2, 6, 4,
        6, 3, 1, 2, 5, 4, 4, 5, 1, 2, 3, 6, 4, 6, 1, 2, 3, 5, 5, 4, 1, 2, 3,
        6, 6, 4, 1, 2, 3, 5, 5, 6, 1, 2, 3, 4, 6, 5, 1, 2, 3, 4, 4, 5, 1, 2,
        6, 3, 4, 6, 1, 2, 5, 3, 5, 4, 1, 2, 6, 3, 6, 4, 1, 2, 5, 3, 5, 6, 1,
        2, 4, 3, 6, 5, 1, 2, 4, 3, 3, 4, 1, 5, 2, 6, 3, 4, 1, 6, 2, 5, 3, 5,
        1, 4, 2, 6, 3, 6, 1, 4, 2, 5, 3, 5, 1, 6, 2, 4, 3, 6, 1, 5, 2, 4, 4,
        3, 1, 5, 2, 6, 4, 3, 1, 6, 2, 5, 5, 3, 1, 4, 2, 6, 6, 3, 1, 4, 2, 5,
        5, 3, 1, 6, 2, 4, 6, 3, 1, 5, 2, 4, 4, 5, 1, 3, 2, 6, 4, 6, 1, 3, 2,
        5, 5, 4, 1, 3, 2, 6, 6, 4, 1, 3, 2, 5, 5, 6, 1, 3, 2, 4, 6, 5, 1, 3,
        2, 4, 4, 5, 1, 6, 2, 3, 4, 6, 1, 5, 2, 3, 5, 4, 1, 6, 2, 3, 6, 4, 1,
        5, 2, 3, 5, 6, 1, 4, 2, 3, 6, 5, 1, 4, 2, 3, 3, 4, 1, 5, 6, 2, 3, 4,
        1, 6, 5, 2, 3, 5, 1, 4, 6, 2, 3, 6, 1, 4, 5, 2, 3, 5, 1, 6, 4, 2, 3,
        6, 1, 5, 4, 2, 4, 3, 1, 5, 6, 2, 4, 3, 1, 6, 5, 2, 5, 3, 1, 4, 6, 2,
        6, 3, 1, 4, 5, 2, 5, 3, 1, 6, 4, 2, 6, 3, 1, 5, 4, 2, 4, 5, 1, 3, 6,
        2, 4, 6, 1, 3, 5, 2, 5, 4, 1, 3, 6, 2, 6, 4, 1, 3, 5, 2, 5, 6, 1, 3,
        4, 2, 6, 5, 1, 3, 4, 2, 4, 5, 1, 6, 3, 2, 4, 6, 1, 5, 3, 2, 5, 4, 1,
        6, 3, 2, 6, 4, 1, 5, 3, 2, 5, 6, 1, 4, 3, 2, 6, 5, 1, 4, 3, 2, 2, 3,
        4, 1, 5, 6, 2, 3, 4, 1, 6, 5, 2, 3, 5, 1, 4, 6, 2, 3, 6, 1, 4, 5, 2,
        3, 5, 1, 6, 4, 2, 3, 6, 1, 5, 4, 2, 4, 3, 1, 5, 6, 2, 4, 3, 1, 6, 5,
        2, 5, 3, 1, 4, 6, 2, 6, 3, 1, 4, 5, 2, 5, 3, 1, 6, 4, 2, 6, 3, 1, 5,
        4, 2, 4, 5, 1, 3, 6, 2, 4, 6, 1, 3, 5, 2, 5, 4, 1, 3, 6, 2, 6, 4, 1,
        3, 5, 2, 5, 6, 1, 3, 4, 2, 6, 5, 1, 3, 4, 2, 4, 5, 1, 6, 3, 2, 4, 6,
        1, 5, 3, 2, 5, 4, 1, 6, 3, 2, 6, 4, 1, 5, 3, 2, 5, 6, 1, 4, 3, 2, 6,
        5, 1, 4, 3, 3, 2, 4, 1, 5, 6, 3, 2, 4, 1, 6, 5, 3, 2, 5, 1, 4, 6, 3,
        2, 6, 1, 4, 5, 3, 2, 5, 1, 6, 4, 3, 2, 6, 1, 5, 4, 4, 2, 3, 1, 5, 6,
        4, 2, 3, 1, 6, 5, 5, 2, 3, 1, 4, 6, 6, 2, 3, 1, 4, 5, 5, 2, 3, 1, 6,
        4, 6, 2, 3, 1, 5, 4, 4, 2, 5, 1, 3, 6, 4, 2, 6, 1, 3, 5, 5, 2, 4, 1,
        3, 6, 6, 2, 4, 1, 3, 5, 5, 2, 6, 1, 3, 4, 6, 2, 5, 1, 3, 4, 4, 2, 5,
        1, 6, 3, 4, 2, 6, 1, 5, 3, 5, 2, 4, 1, 6, 3, 6, 2, 4, 1, 5, 3, 5, 2,
        6, 1, 4, 3, 6, 2, 5, 1, 4, 3, 3, 4, 2, 1, 5, 6, 3, 4, 2, 1, 6, 5, 3,
        5, 2, 1, 4, 6, 3, 6, 2, 1, 4, 5, 3, 5, 2, 1, 6, 4, 3, 6, 2, 1, 5, 4,
        4, 3, 2, 1, 5, 6, 4, 3, 2, 1, 6, 5, 5, 3, 2, 1, 4, 6, 6, 3, 2, 1, 4,
        5, 5, 3, 2, 1, 6, 4, 6, 3, 2, 1, 5, 4, 4, 5, 2, 1, 3, 6, 4, 6, 2, 1,
        3, 5, 5, 4, 2, 1, 3, 6, 6, 4, 2, 1, 3, 5, 5, 6, 2, 1, 3, 4, 6, 5, 2,
        1, 3, 4, 4, 5, 2, 1, 6, 3, 4, 6, 2, 1, 5, 3, 5, 4, 2, 1, 6, 3, 6, 4,
        2, 1, 5, 3, 5, 6, 2, 1, 4, 3, 6, 5, 2, 1, 4, 3, 3, 4, 5, 1, 2, 6, 3,
        4, 6, 1, 2, 5, 3, 5, 4, 1, 2, 6, 3, 6, 4, 1, 2, 5, 3, 5, 6, 1, 2, 4,
        3, 6, 5, 1, 2, 4, 4, 3, 5, 1, 2, 6, 4, 3, 6, 1, 2, 5, 5, 3, 4, 1, 2,
        6, 6, 3, 4, 1, 2, 5, 5, 3, 6, 1, 2, 4, 6, 3, 5, 1, 2, 4, 4, 5, 3, 1,
        2, 6, 4, 6, 3, 1, 2, 5, 5, 4, 3, 1, 2, 6, 6, 4, 3, 1, 2, 5, 5, 6, 3,
        1, 2, 4, 6, 5, 3, 1, 2, 4, 4, 5, 6, 1, 2, 3, 4, 6, 5, 1, 2, 3, 5, 4,
        6, 1, 2, 3, 6, 4, 5, 1, 2, 3, 5, 6, 4, 1, 2, 3, 6, 5, 4, 1, 2, 3, 3,
        4, 5, 1, 6, 2, 3, 4, 6, 1, 5, 2, 3, 5, 4, 1, 6, 2, 3, 6, 4, 1, 5, 2,
        3, 5, 6, 1, 4, 2, 3, 6, 5, 1, 4, 2, 4, 3, 5, 1, 6, 2, 4, 3, 6, 1, 5,
        2, 5, 3, 4, 1, 6, 2, 6, 3, 4, 1, 5, 2, 5, 3, 6, 1, 4, 2, 6, 3, 5, 1,
        4, 2, 4, 5, 3, 1, 6, 2, 4, 6, 3, 1, 5, 2, 5, 4, 3, 1, 6, 2, 6, 4, 3,
        1, 5, 2, 5, 6, 3, 1, 4, 2, 6, 5, 3, 1, 4, 2, 4, 5, 6, 1, 3, 2, 4, 6,
        5, 1, 3, 2, 5, 4, 6, 1, 3, 2, 6, 4, 5, 1, 3, 2, 5, 6, 4, 1, 3, 2, 6,
        5, 4, 1, 3, 2, 2, 3, 4, 5, 1, 6, 2, 3, 4, 6, 1, 5, 2, 3, 5, 4, 1, 6,
        2, 3, 6, 4, 1, 5, 2, 3, 5, 6, 1, 4, 2, 3, 6, 5, 1, 4, 2, 4, 3, 5, 1,
        6, 2, 4, 3, 6, 1, 5, 2, 5, 3, 4, 1, 6, 2, 6, 3, 4, 1, 5, 2, 5, 3, 6,
        1, 4, 2, 6, 3, 5, 1, 4, 2, 4, 5, 3, 1, 6, 2, 4, 6, 3, 1, 5, 2, 5, 4,
        3, 1, 6, 2, 6, 4, 3, 1, 5, 2, 5, 6, 3, 1, 4, 2, 6, 5, 3, 1, 4, 2, 4,
        5, 6, 1, 3, 2, 4, 6, 5, 1, 3, 2, 5, 4, 6, 1, 3, 2, 6, 4, 5, 1, 3, 2,
        5, 6, 4, 1, 3, 2, 6, 5, 4, 1, 3, 3, 2, 4, 5, 1, 6, 3, 2, 4, 6, 1, 5,
        3, 2, 5, 4, 1, 6, 3, 2, 6, 4, 1, 5, 3, 2, 5, 6, 1, 4, 3, 2, 6, 5, 1,
        4, 4, 2, 3, 5, 1, 6, 4, 2, 3, 6, 1, 5, 5, 2, 3, 4, 1, 6, 6, 2, 3, 4,
        1, 5, 5, 2, 3, 6, 1, 4, 6, 2, 3, 5, 1, 4, 4, 2, 5, 3, 1, 6, 4, 2, 6,
        3, 1, 5, 5, 2, 4, 3, 1, 6, 6, 2, 4, 3, 1, 5, 5, 2, 6, 3, 1, 4, 6, 2,
        5, 3, 1, 4, 4, 2, 5, 6, 1, 3, 4, 2, 6, 5, 1, 3, 5, 2, 4, 6, 1, 3, 6,
        2, 4, 5, 1, 3, 5, 2, 6, 4, 1, 3, 6, 2, 5, 4, 1, 3, 3, 4, 2, 5, 1, 6,
        3, 4, 2, 6, 1, 5, 3, 5, 2, 4, 1, 6, 3, 6, 2, 4, 1, 5, 3, 5, 2, 6, 1,
        4, 3, 6, 2, 5, 1, 4, 4, 3, 2, 5, 1, 6, 4, 3, 2, 6, 1, 5, 5, 3, 2, 4,
        1, 6, 6, 3, 2, 4, 1, 5, 5, 3, 2, 6, 1, 4, 6, 3, 2, 5, 1, 4, 4, 5, 2,
        3, 1, 6, 4, 6, 2, 3, 1, 5, 5, 4, 2, 3, 1, 6, 6, 4, 2, 3, 1, 5, 5, 6,
        2, 3, 1, 4, 6, 5, 2, 3, 1, 4, 4, 5, 2, 6, 1, 3, 4, 6, 2, 5, 1, 3, 5,
        4, 2, 6, 1, 3, 6, 4, 2, 5, 1, 3, 5, 6, 2, 4, 1, 3, 6, 5, 2, 4, 1, 3,
        3, 4, 5, 2, 1, 6, 3, 4, 6, 2, 1, 5, 3, 5, 4, 2, 1, 6, 3, 6, 4, 2, 1,
        5, 3, 5, 6, 2, 1, 4, 3, 6, 5, 2, 1, 4, 4, 3, 5, 2, 1, 6, 4, 3, 6, 2,
        1, 5, 5, 3, 4, 2, 1, 6, 6, 3, 4, 2, 1, 5, 5, 3, 6, 2, 1, 4, 6, 3, 5,
        2, 1, 4, 4, 5, 3, 2, 1, 6, 4, 6, 3, 2, 1, 5, 5, 4, 3, 2, 1, 6, 6, 4,
        3, 2, 1, 5, 5, 6, 3, 2, 1, 4, 6, 5, 3, 2, 1, 4, 4, 5, 6, 2, 1, 3, 4,
        6, 5, 2, 1, 3, 5, 4, 6, 2, 1, 3, 6, 4, 5, 2, 1, 3, 5, 6, 4, 2, 1, 3,
        6, 5, 4, 2, 1, 3, 3, 4, 5, 6, 1, 2, 3, 4, 6, 5, 1, 2, 3, 5, 4, 6, 1,
        2, 3, 6, 4, 5, 1, 2, 3, 5, 6, 4, 1, 2, 3, 6, 5, 4, 1, 2, 4, 3, 5, 6,
        1, 2, 4, 3, 6, 5, 1, 2, 5, 3, 4, 6, 1, 2, 6, 3, 4, 5, 1, 2, 5, 3, 6,
        4, 1, 2, 6, 3, 5, 4, 1, 2, 4, 5, 3, 6, 1, 2, 4, 6, 3, 5, 1, 2, 5, 4,
        3, 6, 1, 2, 6, 4, 3, 5, 1, 2, 5, 6, 3, 4, 1, 2, 6, 5, 3, 4, 1, 2, 4,
        5, 6, 3, 1, 2, 4, 6, 5, 3, 1, 2, 5, 4, 6, 3, 1, 2, 6, 4, 5, 3, 1, 2,
        5, 6, 4, 3, 1, 2, 6, 5, 4, 3, 1, 2, 2, 3, 4, 5, 6, 1, 2, 3, 4, 6, 5,
        1, 2, 3, 5, 4, 6, 1, 2, 3, 6, 4, 5, 1, 2, 3, 5, 6, 4, 1, 2, 3, 6, 5,
        4, 1, 2, 4, 3, 5, 6, 1, 2, 4, 3, 6, 5, 1, 2, 5, 3, 4, 6, 1, 2, 6, 3,
        4, 5, 1, 2, 5, 3, 6, 4, 1, 2, 6, 3, 5, 4, 1, 2, 4, 5, 3, 6, 1, 2, 4,
        6, 3, 5, 1, 2, 5, 4, 3, 6, 1, 2, 6, 4, 3, 5, 1, 2, 5, 6, 3, 4, 1, 2,
        6, 5, 3, 4, 1, 2, 4, 5, 6, 3, 1, 2, 4, 6, 5, 3, 1, 2, 5, 4, 6, 3, 1,
        2, 6, 4, 5, 3, 1, 2, 5, 6, 4, 3, 1, 2, 6, 5, 4, 3, 1, 3, 2, 4, 5, 6,
        1, 3, 2, 4, 6, 5, 1, 3, 2, 5, 4, 6, 1, 3, 2, 6, 4, 5, 1, 3, 2, 5, 6,
        4, 1, 3, 2, 6, 5, 4, 1, 4, 2, 3, 5, 6, 1, 4, 2, 3, 6, 5, 1, 5, 2, 3,
        4, 6, 1, 6, 2, 3, 4, 5, 1, 5, 2, 3, 6, 4, 1, 6, 2, 3, 5, 4, 1, 4, 2,
        5, 3, 6, 1, 4, 2, 6, 3, 5, 1, 5, 2, 4, 3, 6, 1, 6, 2, 4, 3, 5, 1, 5,
        2, 6, 3, 4, 1, 6, 2, 5, 3, 4, 1, 4, 2, 5, 6, 3, 1, 4, 2, 6, 5, 3, 1,
        5, 2, 4, 6, 3, 1, 6, 2, 4, 5, 3, 1, 5, 2, 6, 4, 3, 1, 6, 2, 5, 4, 3,
        1, 3, 4, 2, 5, 6, 1, 3, 4, 2, 6, 5, 1, 3, 5, 2, 4, 6, 1, 3, 6, 2, 4,
        5, 1, 3, 5, 2, 6, 4, 1, 3, 6, 2, 5, 4, 1, 4, 3, 2, 5, 6, 1, 4, 3, 2,
        6, 5, 1, 5, 3, 2, 4, 6, 1, 6, 3, 2, 4, 5, 1, 5, 3, 2, 6, 4, 1, 6, 3,
        2, 5, 4, 1, 4, 5, 2, 3, 6, 1, 4, 6, 2, 3, 5, 1, 5, 4, 2, 3, 6, 1, 6,
        4, 2, 3, 5, 1, 5, 6, 2, 3, 4, 1, 6, 5, 2, 3, 4, 1, 4, 5, 2, 6, 3, 1,
        4, 6, 2, 5, 3, 1, 5, 4, 2, 6, 3, 1, 6, 4, 2, 5, 3, 1, 5, 6, 2, 4, 3,
        1, 6, 5, 2, 4, 3, 1, 3, 4, 5, 2, 6, 1, 3, 4, 6, 2, 5, 1, 3, 5, 4, 2,
        6, 1, 3, 6, 4, 2, 5, 1, 3, 5, 6, 2, 4, 1, 3, 6, 5, 2, 4, 1, 4, 3, 5,
        2, 6, 1, 4, 3, 6, 2, 5, 1, 5, 3, 4, 2, 6, 1, 6, 3, 4, 2, 5, 1, 5, 3,
        6, 2, 4, 1, 6, 3, 5, 2, 4, 1, 4, 5, 3, 2, 6, 1, 4, 6, 3, 2, 5, 1, 5,
        4, 3, 2, 6, 1, 6, 4, 3, 2, 5, 1, 5, 6, 3, 2, 4, 1, 6, 5, 3, 2, 4, 1,
        4, 5, 6, 2, 3, 1, 4, 6, 5, 2, 3, 1, 5, 4, 6, 2, 3, 1, 6, 4, 5, 2, 3,
        1, 5, 6, 4, 2, 3, 1, 6, 5, 4, 2, 3, 1, 3, 4, 5, 6, 2, 1, 3, 4, 6, 5,
        2, 1, 3, 5, 4, 6, 2, 1, 3, 6, 4, 5, 2, 1, 3, 5, 6, 4, 2, 1, 3, 6, 5,
        4, 2, 1, 4, 3, 5, 6, 2, 1, 4, 3, 6, 5, 2, 1, 5, 3, 4, 6, 2, 1, 6, 3,
        4, 5, 2, 1, 5, 3, 6, 4, 2, 1, 6, 3, 5, 4, 2, 1, 4, 5, 3, 6, 2, 1, 4,
        6, 3, 5, 2, 1, 5, 4, 3, 6, 2, 1, 6, 4, 3, 5, 2, 1, 5, 6, 3, 4, 2, 1,
        6, 5, 3, 4, 2, 1, 4, 5, 6, 3, 2, 1, 4, 6, 5, 3, 2, 1, 5, 4, 6, 3, 2,
        1, 6, 4, 5, 3, 2, 1, 5, 6, 4, 3, 2, 1, 6, 5, 4, 3, 2, 1};

    inta *at;
    inta nd = (bodies(f1,left));
    inta np = tPerms(bodies(f1,left));
    leftChar--;

    if ( leftChar < 0 || leftChar >= np ){
        printf("oops...%dretry\n %d %d", leftChar,left, equals);
        exit(1);
    }
    switch (bodies(f1,left)){
        case nada:
            return 0;
        case one :
            at = at1+nd*(leftChar);
            break;
        case two:
            at = at2+nd*(leftChar);
            break;
        case three:
            at = at3+nd*(leftChar);
            break;
        case four:
            at = at4+nd*(leftChar);
            break;
        case five:
            at = at5+nd*(leftChar);
            break;
        case six:
            at = at6+nd*(leftChar);
            break;
        }
    inta v[12],u[12],cb,cb2,cb1,ll = vector1Len(f1, dim),prev,d;
    floata * str = streams(f1, left, lspin, dim)+ l * alloc(f1, left, dim);
    floata * out = streams(f1, equals, espin, dim)+e*alloc(f1, equals, dim);
//    printf("%d:%d %d %d\n",leftChar, at[0],at[1],at[2]);
//    fflush(stdout);

    cb1 = Power(ll,nd);
    for ( cb = 0 ; cb < cb1; cb++)
    {
        
        
        prev = 1;
        for ( d = 0 ; d < nd ; d++){
            v[d] = (cb/prev)%ll;
            prev *= ll;
        }
        
        for ( d = 0 ; d < nd ; d++){
            u[d] = v[at[d]-1];//OTHER
            //u[at[d]-1] = v[d];//ORIGINAL

           //WHAT DOES THIS MEAN?
            //3 1 2
            
            //CHOOSE:
            //OTHER
            //u[1] = v[3]
            //u[2] = v[1]
            //u[3] = v[2]
//
//OR
//          ORIGINAL
            // u[3] = v[1]
            // u[1] = v[2]
            // u[2] = v[3]
            
            //THE OTHER IS CORRECT!!
            
        }
        
        cb2 = 0;
        prev = 1;
        for ( d = 0 ; d < nd ; d++){
            cb2 += prev * u[d];
            prev *= ll;
        }
        out[cb2] = str[cb];

    }
    
    return 0;
}



/**
 * group action of leftChar
*/
inta tPermute(inta rank,   sinc_label f1, inta leftChar ,   division left, inta lspin,   division equals, inta espin){
    inta l,space;
    
    for ( l = 0; l < CanonicalRank(f1,left,lspin); l++){
        for ( space = 0; space < SPACE ; space++)
            if ( f1.canon[space].body != nada)
            tPermuteOne(rank, f1, space, leftChar, left, l, lspin, equals,f1.name[equals].Current[espin], espin);
        f1.name[equals].Current[espin]++;
    }

    return 0;
}



/**
 *Conduct inner products with group action
*/
inta tAllCompPermMultiplyMP( sinc_label  f1 ,   division left ,inta lspin,   division right ,inta rspin, double * sequ){
    inta rank = 0;
    inta dim,l,r;
    double prod;
    if ( CanonicalRank(f1, left, lspin ) * CanonicalRank(f1, right, rspin ) == 0)
        return 0;
    
    if ( bodies(f1, left ) != bodies(f1,right)){
        printf("tGetType real!\n");
        exit(0);
    }
    inta i,nPerm=tPerms(bodies(f1,left));
    
    for ( i = 1; i <= nPerm ; i++){
        sequ[i] = 0.;
    }
    
    if ( CanonicalRank(f1, left, lspin ) * CanonicalRank(f1, right, rspin ) == 0){
        return 0;
    }

    for ( i = 1; i <= nPerm ; i++){
        sequ[i] =0 ;
        for ( l = 0; l < CanonicalRank(f1,left,lspin) ; l++)
            for ( r = 0; r < CanonicalRank(f1,left,rspin) ; r++){
                prod = 1;
                for ( dim = 0; dim < SPACE ; dim++)
                    if ( f1.canon[dim].body != nada)
                        prod *= tDOT(rank, f1, dim, i, left, l, lspin, CDT, right, r, rspin);
                sequ[i] += prod;
            }
    }
    return nPerm;
}

/**
 *Manage inner products with group action
 */
inta tTabulateInnerProjection( sinc_label  f1 ,   division vec, double *up){

    if ( bodies(f1,vec ) == one ){
        up[0] = 1;
        return 1;
    }
    inta g,p,nPerm=0,nGroup=0;
    
    
    if ( bodies(f1, vec ) == two ){
        nPerm = 2;
        nGroup = 2;
    }
    else if ( bodies ( f1, vec ) == three ){
        nPerm = 6;
        nGroup = 3;
        
    }
    else if ( bodies (f1, vec ) == four  ){
        nPerm = 24;
        nGroup = 5;
    }else {
        printf("opps\n");
        fflush(stdout);
        exit(0);
    }
    double buff[25];

    tAllCompPermMultiplyMP( f1, vec, 0, vec,0, buff);
    for ( g = 1; g <= nGroup ; g++)
        for ( p = 1; p <= nPerm ; p++){
            up[g]  += tGetProjection(bodies(f1,vec), g, p)*buff[p];
        }
    if ( CanonicalRank(f1, vec ,1)){
        tAllCompPermMultiplyMP( f1, vec, 1, vec,1, buff);
        for ( g = 1; g <= nGroup ; g++)
            for ( p = 1; p <= nPerm ; p++)
                up[g]  += tGetProjection(bodies(f1,vec), g, p)*buff[p];
    }
    return nGroup;
}

/**
 *Number of irreps
 */
inta tSize(  bodyType bd){
    switch (bd) {
        case nada:
            return 0;
        case one:
            return 1;
        case two:
            return 2;
        case three:
            return 3;
        case four:
            return 5;
        case five:
            return 1;
        case six:
            return 1;
    }
}


/**
 *Number of group actions
 */
inta tPerms(  bodyType bd){
    switch (bd) {
        case nada:
            return 0;
        case one:
            return 1;
        case two:
            return 2;
        case three:
            return 6;
        case four:
            return 24;
        case five:
            return 120;
        case six:
            return 720;
    }
}



