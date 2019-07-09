/*
 *  saUtil.c
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

#include "saUtil.h"

INT_TYPE tFil ( struct sinc_label  f1, enum division A, enum division v , INT_TYPE * i ){
    INT_TYPE n1[3],space;
    length1(f1,n1);
    if ( bodies (f1, v ) == two ){
        for ( space = 0; space < SPACE ; space++){
            cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space]*n1[space],1, streams(f1,diagonal1VectorA,0,space),1);
        }
        for ( space = 0; space < SPACE ; space++){
            cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space+3]*n1[space],1, streams(f1,diagonal1VectorB,0,space),1);
        }
        f1.tulip[diagonal1VectorA].Current[0] = 1;
        f1.tulip[diagonal1VectorB].Current[0] = 1;

        f1.tulip[v].Current[0] = 0;

        tOuterProductSu(f1, diagonal1VectorA, 0, diagonal1VectorB, 0, v, 0);

    }else
        if ( bodies (f1, v ) == three ){
            for ( space = 0; space < SPACE ; space++){
                cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space]*n1[space],1, streams(f1,diagonal1VectorA,0,space),1);
            }
            for ( space = 0; space < SPACE ; space++){
                cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space+3]*n1[space],1, streams(f1,diagonal1VectorB,0,space),1);
            }
            for ( space = 0; space < SPACE ; space++){
                cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space+6]*n1[space],1, streams(f1,diagonal1VectorC,0,space),1);
            }
            f1.tulip[diagonal1VectorA].Current[0] = 1;
            f1.tulip[diagonal1VectorB].Current[0] = 1;
            f1.tulip[diagonal1VectorC].Current[0] = 1;
            f1.tulip[diagonal2VectorA].Current[0] = 0;
            f1.tulip[v].Current[0] = 0;

            tOuterProductSu(f1, diagonal1VectorA, 0, diagonal1VectorB, 0, diagonal2VectorA, 0);
            tOuterProductSu(f1, diagonal2VectorA, 0, diagonal1VectorC, 0, v, 0);

        }else     if ( bodies (f1, v ) == four ){
            for ( space = 0; space < SPACE ; space++){
                cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space]*n1[space],1, streams(f1,diagonal1VectorA,0,space),1);
            }
            for ( space = 0; space < SPACE ; space++){
                cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space+3]*n1[space],1, streams(f1,diagonal1VectorB,0,space),1);
            }
            for ( space = 0; space < SPACE ; space++){
                cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space+6]*n1[space],1, streams(f1,diagonal1VectorC,0,space),1);
            }
            for ( space = 0; space < SPACE ; space++){
                cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space+9]*n1[space],1, streams(f1,diagonal1VectorD,0,space),1);
            }
            f1.tulip[diagonal1VectorA].Current[0] = 1;
            f1.tulip[diagonal1VectorB].Current[0] = 1;
            f1.tulip[diagonal1VectorC].Current[0] = 1;
            f1.tulip[diagonal1VectorD].Current[0] = 1;
            f1.tulip[diagonal2VectorA].Current[0] = 0;
            f1.tulip[diagonal2VectorB].Current[0] = 0;
            f1.tulip[v].Current[0] = 0;


            tOuterProductSu(f1, diagonal1VectorA, 0, diagonal1VectorB, 0, diagonal2VectorA, 0);
            tOuterProductSu(f1, diagonal1VectorC, 0, diagonal1VectorD, 0, diagonal2VectorB, 0);
            tOuterProductSu(f1, diagonal2VectorA, 0, diagonal2VectorB, 0, v, 0);

        }

    return 0;
}




//INT_TYPE tInnerTest( struct field * f1, enum division A ,enum division B){
//    char c;
//    INT_TYPE n1[3],space,i[100],j[100],nPerm,ii;
//    double seq[100];
//    length1(f1,n1);
//    if ( bodies ( f1, A ) != one || bodies (f1, B ) != one ){
//        printf ("ment for onebody\n");
//        exit(0);
//    }
//    
//    for ( space = 0; space < SPACE ; space++)
//        tdsyev (0,f1,'V',n1[space],streams(f1,A,0,space),n1[space],streams(f1,B,0,space));
//    for ( ii = 0; ii < 100 ;ii++)
//        i[ii] = ii/3 % n1[0];
//    
//    tFil(f1, A, copyThreeVector, i);
//    
//    if ( bodies ( f1, eigenVectors )== three ){
//        for ( c = 0; c < 5 ; c++){
//            tClear(f1, copyFourVector);
//            printf("\n\nCHAR %c\n", 'A'+c);
//            tPermute(0, f1, 'A'+c, copyThreeVector, 0, copyFourVector, 0);
//            nPerm = tAllCompPermMultiplyMP(0, f1, copyThreeVector, 0, copyFourVector, 0, seq);
//            for ( ii = 0 ; ii < nPerm; ii++)
//                printf("%lld\t%f\n",ii+1,seq[ii]);
//        }
//    }else
//        if ( bodies ( f1, eigenVectors )== four ){
//            for ( c = 0; c < 23 ; c++){
//                tClear(f1,copyFourVector);
//                printf("\n\nCHAR %c\n", 'a'+c);
//                tPermute(0, f1, 'a'+c, copyThreeVector, 0, copyFourVector, 0);
//                nPerm = tAllCompPermMultiplyMP(0, f1, copyThreeVector, 0, copyFourVector, 0, seq);
//                for ( ii = 0 ; ii < nPerm; ii++)
//                    printf("%lld\t%f\n",ii+1,seq[ii]);
//            }
//        }
//    
//    
//    if ( bodies ( f1, eigenVectors )== three ){
//        for ( c = 'A'; c <= 'E' ; c++){
//            if ( c == 'A' ){
//                //            train[0] = 'T';//(1)            123
//                //            train[1] = 'A';//(123)          231
//                //            train[2] = 'B';//(123).(123)    312
//                //            train[3] = 'C';//(12)           213
//                //            train[4] = 'D';//(13)           321
//                //            train[5] = 'E';//(23)           132
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+3];
//                    j[space+3] = i[space+6];
//                    j[space+6] = i[space];
//                    
//                }
//            }else
//                if ( c == 'B' ){
//                    //            train[0] = 'T';//(1)            123
//                    //            train[1] = 'A';//(123)          231
//                    //            train[2] = 'B';//(123).(123)    312
//                    //            train[3] = 'C';//(12)           213
//                    //            train[4] = 'D';//(13)           321
//                    //            train[5] = 'E';//(23)           132
//                    for ( space = 0 ; space < SPACE ; space++){
//                        j[space] = i[space+6];
//                        j[space+3] = i[space];
//                        j[space+6] = i[space+3];
//                        
//                    }
//                }
//                else            if ( c == 'C' ){
//                    //            train[0] = 'T';//(1)            123
//                    //            train[1] = 'A';//(123)          231
//                    //            train[2] = 'B';//(123).(123)    312
//                    //            train[3] = 'C';//(12)           213
//                    //            train[4] = 'D';//(13)           321
//                    //            train[5] = 'E';//(23)           132
//                    for ( space = 0 ; space < SPACE ; space++){
//                        j[space] = i[space+3];
//                        j[space+3] = i[space];
//                        j[space+6] = i[space+6];
//                        
//                    }
//                }else             if ( c == 'D' ){
//                    //            train[0] = 'T';//(1)            123
//                    //            train[1] = 'A';//(123)          231
//                    //            train[2] = 'B';//(123).(123)    312
//                    //            train[3] = 'C';//(12)           213
//                    //            train[4] = 'D';//(13)           321
//                    //            train[5] = 'E';//(23)           132
//                    for ( space = 0 ; space < SPACE ; space++){
//                        j[space] = i[space+6];
//                        j[space+3] = i[space+3];
//                        j[space+6] = i[space];
//                        
//                    }
//                }else             if ( c == 'E' ){
//                    //            train[0] = 'T';//(1)            123
//                    //            train[1] = 'A';//(123)          231
//                    //            train[2] = 'B';//(123).(123)    312
//                    //            train[3] = 'C';//(12)           213
//                    //            train[4] = 'D';//(13)           321
//                    //            train[5] = 'E';//(23)           132
//                    for ( space = 0 ; space < SPACE ; space++){
//                        j[space] = i[space];
//                        j[space+3] = i[space+6];
//                        j[space+6] = i[space+3];
//                        
//                    }
//                }
//            
//            
//            
//            
//            
//            
//            tFil(f1, A, copyFourVector, j);
//            nPerm = tAllCompPermMultiplyMP(0, f1, copyThreeVector, 0, copyFourVector, 0, seq);
//            for ( ii = 0 ; ii < nPerm; ii++)
//                printf("%lld\t%f\n",ii+1,seq[ii]);
//            
//        }
//    }
//    if ( bodies ( f1, eigenVectors )== four ){
//        for ( c = 'a'; c <= 'w' ; c++){
//            
//            // 0,0,0,0 : i,j,k,l   :: T 0
//            // 0,2,1,2 : i,j,l,k   :: 'a'24
//            // 1,3,1,2 : i,k,j,l   :: 'b'20
//            // 1,1,0,0 : i,k,l,j   :: 'c'4
//            // 0,3,1,0 : i,l,j,k   :: 'd'12
//            // 1,2,1,3 : i,l,k,j   :: 'e'22
//            // 1,0,0,0 : j,i,k,l   :: 'f'2
//            // 1,2,1,2 : j,i,l,k   :: 'g' 23
//            // 0,3,1,2 : j,k,i,l   :: 'h'19
//            // 0,1,0,0 : j,k,l,i   :: 'i'3
//            // 1,3,1,0 : j,l,i,k   :: 'j'13
//            // 0,2,1,3 : j,l,k,i   :: 'k'21
//            // 0,2,1,1 : k,i,j,l   :: 'l'15
//            // 1,1,1,0 : k,i,l,j   :: 'm'10
//            // 1,2,1,1 : k,j,i,l   :: 'n'16
//            // 0,1,1,0 : k,j,l,i   :: 'o'9
//            // 0,2,0,0 : k,l,i,j   :: 'p'5
//            // 1,1,0,1 : k,l,j,i   :: 'q'14
//            // 0,3,0,0 : l,i,j,k   :: 'r'7
//            // 1,3,1,1 : l,i,k,j   :: 's'18
//            // 1,3,0,0 : l,j,i,k   :: 't'8
//            // 0,3,1,1 : l,j,k,i   :: 'u'17
//            // 1,2,0,0 : l,k,i,j   :: 'v'6
//            // 1,2,1,0 : l,k,j,i   :: 'w'11
//            if ( c == 'a' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space];
//                    j[space+3] = i[space+3];
//                    j[space+6] = i[space+9];
//                    j[space+9] = i[space+6];
//                }
//            }else            if ( c == 'b' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space];
//                    j[space+3] = i[space+6];
//                    j[space+6] = i[space+3];
//                    j[space+9] = i[space+9];
//                }
//            }else            if ( c == 'c' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space];
//                    j[space+3] = i[space+6];
//                    j[space+6] = i[space+9];
//                    j[space+9] = i[space+3];
//                }
//            }else            if ( c == 'd' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space];
//                    j[space+3] = i[space+9];
//                    j[space+6] = i[space+3];
//                    j[space+9] = i[space+6];
//                }
//            }else            if ( c == 'e' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space];
//                    j[space+3] = i[space+9];
//                    j[space+6] = i[space+6];
//                    j[space+9] = i[space+3];
//                }
//            }else            if ( c == 'f' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+3];
//                    j[space+3] = i[space];
//                    j[space+6] = i[space+6];
//                    j[space+9] = i[space+9];
//                }
//            }else            if ( c == 'g' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+3];
//                    j[space+3] = i[space];
//                    j[space+6] = i[space+9];
//                    j[space+9] = i[space+6];
//                }
//            }else            if ( c == 'h' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+3];
//                    j[space+3] = i[space+6];
//                    j[space+6] = i[space];
//                    j[space+9] = i[space+9];
//                }
//            }else            if ( c == 'i' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+3];
//                    j[space+3] = i[space+6];
//                    j[space+6] = i[space+9];
//                    j[space+9] = i[space];
//                }
//            }else            if ( c == 'j' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+3];
//                    j[space+3] = i[space+9];
//                    j[space+6] = i[space];
//                    j[space+9] = i[space+6];
//                }
//            }else            if ( c == 'k' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+3];
//                    j[space+3] = i[space+9];
//                    j[space+6] = i[space+6];
//                    j[space+9] = i[space];
//                }
//            }else            if ( c == 'l' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+6];
//                    j[space+3] = i[space];
//                    j[space+6] = i[space+3];
//                    j[space+9] = i[space+9];
//                }
//            }else            if ( c == 'm' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+6];
//                    j[space+3] = i[space];
//                    j[space+6] = i[space+9];
//                    j[space+9] = i[space+3];
//                }
//            }else            if ( c == 'n' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+6];
//                    j[space+3] = i[space+3];
//                    j[space+6] = i[space];
//                    j[space+9] = i[space+9];
//                }
//            }else            if ( c == 'o' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+6];
//                    j[space+3] = i[space+3];
//                    j[space+6] = i[space+9];
//                    j[space+9] = i[space];
//                }
//            }else            if ( c == 'p' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+6];
//                    j[space+3] = i[space+9];
//                    j[space+6] = i[space];
//                    j[space+9] = i[space+3];
//                }
//            }else            if ( c == 'q' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+6];
//                    j[space+3] = i[space+9];
//                    j[space+6] = i[space+3];
//                    j[space+9] = i[space];
//                }
//            }else            if ( c == 'r' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+9];
//                    j[space+3] = i[space];
//                    j[space+6] = i[space+3];
//                    j[space+9] = i[space+6];
//                }
//            }else            if ( c == 's' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+9];
//                    j[space+3] = i[space];
//                    j[space+6] = i[space+6];
//                    j[space+9] = i[space+3];
//                }
//            }else            if ( c == 't' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+9];
//                    j[space+3] = i[space+3];
//                    j[space+6] = i[space];
//                    j[space+9] = i[space+6];
//                }
//            }else            if ( c == 'u' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+9];
//                    j[space+3] = i[space+3];
//                    j[space+6] = i[space+6];
//                    j[space+9] = i[space];
//                }
//            }else            if ( c == 'v' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+9];
//                    j[space+3] = i[space+6];
//                    j[space+6] = i[space];
//                    j[space+9] = i[space+3];
//                }
//            }else            if ( c == 'w' ){
//                for ( space = 0 ; space < SPACE ; space++){
//                    j[space] = i[space+9];
//                    j[space+3] = i[space+6];
//                    j[space+6] = i[space+3];
//                    j[space+9] = i[space];
//                }
//            }
//            
//            
//            
//            
//            
//            tFil(f1, A, copyFourVector, j);
//            nPerm = tAllCompPermMultiplyMP(0, f1, copyThreeVector, 0, copyFourVector, 0, seq);
//            for ( ii = 0 ; ii < nPerm; ii++)
//                printf("%lld\t%f\n",ii+1,seq[ii]);
//            
//            
//            
//            
//        }
//        
//        
//        
//        
//        
//        
//        
//        
//        
//        
//        
//        
//    }
//    
//    return 0;
//}


double deg(struct sinc_label f1, INT_TYPE cl ){
    double deg = 0.;
    INT_TYPE Ns = bodies(f1, eigenVectors);
    if ( Ns == 4 ){
        if ( cl == 2 )
            deg = 5.;//5
        else if ( cl ==3 )
            deg = 0.5;//2
        else if ( cl == 4 )
            deg = 1.;//9
        //2^4 = 16 = 5+2+9
    }else
        if ( Ns == 3 ){
            if( cl == 2 )
                deg = 4.;
            else if ( cl == 3 )
                deg = 1.;//8
            //2^3 =  8 = 4+4
        }else
            if ( Ns == 2 ){//
                if ( cl == 2 )
                    deg = 3.;
                else if ( cl == 1 )
                    deg = 1.;
            }
            else
                if ( Ns == 1 ){
                    deg = 2.;
                }
    return deg;
}

double tGetProjection (enum bodyType bd , INT_TYPE irrep , INT_TYPE op ){
    //THESE ARE PROJECTIONS based on enumerated operations.
    //i.e. for 3-electrons
    //Sum[(Op[x],{x,1,6}].v = v_A1
    //
    
    INT_TYPE nsyp=0 ,msyp=0;
        const static double syp2 [] = {
            0.5,0.5,
            0.5,-0.5
    
        };
        const static double syp3 [] = {
            1./6,1./6,1./6, 1./6, 1./6, 1./6,//A1-36 of them
            1./6,1./6,1./6,-1./6,-1./6,-1./6,//A2-36 of them
            2./6,1./6,1./6,0.,0.,0.//EE-144 of them
        };
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
            printf("bod\n");
            exit(0);
        }
    
    
        if ( irrep <= 0 || irrep > msyp )
        {
            printf("he\n %lld", irrep);
            exit(0);
        }
        if ( op < 0 || op >= nsyp ){
            printf("hm\n");
            exit(0);
        }
    
    
        if ( bd == two ){
            return syp2[(irrep-1)*nsyp+op];
        }
        else if ( bd == three ){
            return syp3[(irrep-1)*nsyp+op];
        }
//        else if ( bd == four ){
//            return syp4[(irrep-1)*nsyp+op];
//        }
        return 0.;
}





//Bill's work,  3component breakdown
//double tGetType(enum bodyType bd , INT_TYPE type , INT_TYPE perm ){
//
//    INT_TYPE nsyp=0 ,msyp=0;
//    const static double syp2 [] = {
//        sr2,sr2,
//        sr2,-sr2
//
//    };
//    const static double syp3 [] = {
//        /********/
//        sr6,sr6,sr6, sr6, sr6, sr6,//A1
//        sr6,sr6,sr6,-sr6,-sr6,-sr6,//A2
//
//
//
//
////        hf*sr3,hf*sr3   , -sr3   ,hf*sr3    ,hf*sr3     ,-sr3,
////        hf    , -hf     ,0.      ,-hf       , hf        , 0.,
////        hf*sr3,hf*sr3   , -sr3   ,-hf*sr3   ,-hf*sr3    ,sr3,
////        hf    , -hf     ,0.      ,hf        , -hf       , 0.,
//        //////////////////////////////////
//        /////////////////////////////////
//        ////////////////////////////////
//        ///////////////////////////////
//        //////////////////////////////
//    };
//
//    const static double syp4[] = {
//        /*1: A1*/
//        sr6*hf,sr6*hf,sr6*hf,sr6*hf,    sr6*hf,sr6*hf,sr6*hf,sr6*hf,    sr6*hf,sr6*hf,sr6*hf,sr6*hf,
//        sr6*hf,sr6*hf,sr6*hf,sr6*hf,    sr6*hf,sr6*hf,sr6*hf,sr6*hf,    sr6*hf,sr6*hf,sr6*hf,sr6*hf,
//
//        /*2: A2*/
//        sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,    sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,    sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,
//        sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,    sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,    sr6*hf,-sr6*hf,-sr6*hf,sr6*hf,
//
//        /*3: E-1*/
//        sr6,    0.,0.,-sr6*hf,      -sr6*hf,0.,0., sr6,       -sr6*hf,0.,0.,-sr6*hf,
//        -sr6*hf,0.,0.,-sr6*hf,       sr6,   0.,0.,-sr6*hf,    -sr6*hf,0.,0.,sr6,
//
//        /*4: E-2*/
//        0.,sr6,-sr6*hf,0.,          0.,-sr6*hf, sr6,0.,    0.,-sr6*hf,-sr6*hf,0.,
//        0.,-sr6*hf,-sr6*hf,0.,      0.,sr6, -sr6*hf,0.,    0.,-sr6*hf,sr6,0.,
//
//        /*5: E-3*/
//        0.,0.,sr2*hf,0.,             0.,-sr2*hf,0.,0.,       0.,-sr2*hf ,sr2*hf,0.,
//        0.,sr2*hf ,-sr2*hf,0.,        0.,0.,-sr2*hf,0,        0.,sr2*hf,0.,0.,
//
//        /*6: E-4*/
//        0.,0.,0.,sr2*hf,             -sr2*hf,0.,0.,0.,         -sr2*hf,0.,0.,sr2*hf,
//        sr2*hf,0.,0.,-sr2*hf,       0.,0.,0.,-sr2*hf,           sr2*hf,0.,0.,0.,
//
//        /*7: T1-1*/
//        hf*sr2/sr3,-hf*sr6,-hf*sr6,0.,      0.,-hf*sr6,-hf*sr6,-hf*sr6,     0.,hf*sr6,hf*sr6,0.,
//        0., hf*sr6,-hf*sr6,0.,              -hf*sr6, hf*sr6, hf*sr6,0.,     0.,-hf*sr6,hf*sr6,-hf*sr6,
//
//        /*8: T1-2*/
//        0. , sr3,-hf*hf*hf*sr3,-hf*hf*hf/sr3,                               -hf*hf*hf/sr3,-hf*hf*hf*sr3, -hf*sr3, -hf*sr3,
//        hf*hf*hf/sr3,hf*hf*hf*sr3,hf*hf*hf*sr3,hf*hf*hf/sr3,                  hf*hf*hf/sr3,hf*hf*hf*sr3,-hf*hf*hf*sr3,-hf*hf*hf/sr3,
//        hf*hf*sr3,-hf*hf*sr3,hf*hf*hf*sr3,hf*hf*hf/sr3,                       -hf*hf*hf/sr3,-hf*hf*hf*sr3,-hf*hf*sr3,  hf*hf*sr3,
//
//        /*9: T1-3*/
//        0.,0.,hf*hf*hf/sr3/sr7,-hf*hf*hf*sr7/sr3/sr3/sr3,
//        -hf*hf*hf*sr7/sr3/sr3/sr3,-hf*hf*hf*sr7/sr3,-hf*sr3*sr7, hf*sr3*sr7,
//        -sr3*hf*hf*hf/sr7,1./sr3*hf*hf*hf*sr7,-hf*hf*hf*sr3*sr7/sr5/sr5, hf*hf*hf*sr7/sr3/sr3/sr3,
//        -sr3*hf*hf*hf/sr7,-hf*hf*hf*sr3*sr7/sr5/sr5,-hf*hf*hf/sr3*sr7,sr3*hf*hf*hf/sr7,
//        hf*hf*sr7/sr3, hf*hf*sr3*sr7, hf*hf*hf*sr7/sr3,  hf*hf*hf*sr7/sr3/sr3/sr3,
//        hf*hf*hf*sr3/sr7,-11.*hf*hf*hf*sr3*sr7, hf*hf*sr3*sr7,  -hf*hf*sr3*sr7/sr5/sr5,
//
//
//        /*10: T1-4*/
//        0.,0.,0.,hf*sr2*sr7/sr5/sr3,
//        -sr5*sr2*sr7/sr3, -hf*sr5*sr2*sr7/sr3/sr3/sr3, sr5*sr2*sr3*sr7,-sr5*sr2*sr3*sr7,
//        0.,-hf*sr6*sr7/sr5, sr5*sr2*sr7/sr3, -hf * sr7*sr3*sr5*sr2,
//        -hf*sr5*sr6/sr7,-sr3*sr7*sr5/hf/sr2,hf*sr6*sr7/sr5,0.,
//        sr3*sr7*sr5/sr2, sr7*sr5*sr2/sr3,   hf*sr7*sr5*sr2/sr3/sr3/sr3, sr7*sr5*sr2/sr3,
//        -hf*sr5*sr6/sr7,sr5*sr2*sr3*sr7,-sr3*sr7*sr5/hf/sr2,-sr5*sr2*sr3*sr7,
//
//        /*11: T1-5*/
//        0.,0.,0.,0.,
//        hf/sr3/sr3*sr2*sr5,-hf/sr3/sr3*sr2*sr5,sr3*sr3*sr2*sr5,-sr3*sr3*sr2*sr5,
//        -sr6*sr6*sr2/sr5,sr6*sr6*sr2/sr5,-sr5*sr3*sr3/sr2,sr5*sr3*sr3/sr2,
//        -sr2*sr5*sr3*sr3,sr2*sr5*sr3*sr3,sr6*sr6*sr2/sr5,-sr6*sr6*sr2/sr5,
//        sr5*sr3*sr3/sr2,-sr5*sr3*sr3/sr2,-sr2*sr5*sr6*sr6,sr6*sr6*sr2*sr5,
//        -sr2*sr5*sr3*sr3,sr2*sr5*sr3*sr3,sr2*sr5*sr3*sr3,-sr2*sr5*sr3*sr3,
//
//        /*12: T1-6*/
//        0.,0.,0.,0.,
//        0.,0.,sr3*sr3/sr2,-sr3*sr3/sr2,
//        -sr6*sr6*sr2, sr6*sr6*sr2,sr6*sr6*sr2,-sr6*sr6*sr2,
//        -sr6*sr6*sr2,sr6*sr6*sr2,-sr3*sr3*sr2,sr3*sr3*sr2,
//        sr3*sr3*sr2,-sr3*sr3*sr2,sr6*sr6*sr2,-sr6*sr6*sr2,
//        sr3*sr3*sr2,-sr3*sr3*sr2,-sr3*sr3*sr2,sr3*sr3*sr2,
//
//        /*13: T1-7*/
//        0.,0.,0.,0.,
//        0.,0.,0.,0.,
//        hf*sr6/sr5,-hf*sr5*sr6,-hf*sr6/sr5,hf*sr6*sr5,
//        -hf*sr2*sr5/sr3,hf*sr2*sr5/sr3,-sr6*sr5,-sr6*sr5,
//        sr6*sr5,sr6*sr5,hf*sr2*sr5/sr3,-hf*sr2*sr5/sr3,
//        sr5*sr6,sr5*sr6,-sr5*sr6,-sr5*sr6,
//
//        /*14: T1-8*/
//        0.,0.,0.,0.,
//        0.,0.,0.,0.,
//        0.,sr5,0,-sr5,
//        hf*sr5,-hf*sr5,-hf*sr5,-hf*sr5,
//        hf*sr5,hf*sr5,-hf*sr5,hf*sr5,
//        hf*sr5,hf*sr5,-hf*sr5,-hf*sr5,
//
//        /*15: T1-9*/
//        0.,0.,0.,0.,
//        0.,0.,0.,0.,
//        0.,0.,0.,0.,
//        hf*sr3,-hf*sr3,-hf*sr3,hf*sr3,
//        hf*sr3,-hf*sr3,hf*sr3,-hf*sr3,
//        -hf*sr3,hf*sr3,hf*sr3,-hf*sr3,
//
//        /*16: T2-1*/
//        hf*sr2/sr3,hf*sr2*sr3,hf*sr2*sr3,0.,
//        0.,hf*sr2*sr3,hf*sr2*sr3,-hf*sr2*sr3,
//        0.,-hf*sr2*sr3,-hf*sr2*sr3,0.,
//        0.,-hf*sr2*sr3,hf*sr2*sr3,0.,
//        -hf*sr2*sr3,-hf*sr2*sr3,-hf*sr2*sr3,0.,
//        0.,hf*sr2*sr3,-hf*sr2*sr3,-hf*sr2*sr3,
//
//        /*17: T2-2*/
//        0.,sr3,-hf*hf*hf*sr3,hf*hf*hf/sr3,
//        hf*hf*hf/sr3,-hf*hf*hf*sr3,-hf*sr3,hf*sr3,
//        -hf*hf*hf/sr3,hf*hf*hf*sr3,hf*hf*hf*sr3,-hf*hf*hf/sr3,
//        -hf*hf*hf/sr3,hf*hf*hf*sr3,-hf*hf*hf*sr3,hf*hf*hf/sr3,
//        -hf*hf*sr3,-hf*hf*sr3,hf*hf*hf*sr3,-hf*hf*hf/sr3,
//        hf*hf*hf/sr3,-hf*hf*hf*sr3,-hf*hf*sr3,-hf*hf*sr3,
//
//        /*18: T2-3*/
//        0.,0.,hf*hf*hf/sr3/sr7,hf*hf*hf*sr7/sr3/sr3/sr3,
//        hf*hf*hf*sr7/sr3/sr3/sr3,-hf*hf*hf*sr7/sr3,-hf*sr3*sr7,-hf*sr3*sr7,
//        sr3/sr7*hf*hf*hf, sr7/sr3*hf*hf*hf, -hf*hf*hf*sr3*sr7/sr5/sr5, -hf*hf*hf*sr7/sr3/sr3/sr3,
//        hf*hf*hf*sr3/sr7,-hf*hf*hf*sr3*sr7/sr5/sr5,-hf*hf*hf*sr7/sr3,-hf*hf*hf*sr3/sr7,
//        -hf*hf*sr7/sr3, hf*hf*sr3*sr7,  hf*hf*hf*sr7/sr3, - hf*hf*hf*sr7/sr3/sr3/sr3,
//        -sr3/sr7*hf*hf*hf, -11.*hf*hf*hf*sr3*sr7, hf*hf*sr3*sr7, hf*hf*sr3*sr7/sr5/sr5,
//
//
//        /*19: T2-4*/
//        0.,0.,0.,hf*sr2*sr7/sr3/sr5,
//        -sr7*sr5*sr2/sr3,hf*sr7*sr5*sr2/sr3/sr3/sr3,    - sr7*sr3*sr5*sr2,- sr7*sr3*sr5*sr2,
//        0.,hf/sr5*sr6*sr7,  -sr7*sr5*sr2/sr3, -hf*sr7*sr3*sr5*sr2,
//        -hf*sr5*sr6/sr7,sr7*sr3*sr5/sr2/hf,-hf*sr6*sr7/sr5,0.,
//        sr7*sr3*sr5/sr2,-sr7*sr5*sr2/sr3,-hf*sr7*sr5*sr2/sr3/sr3/sr3,sr7*sr5*sr2/sr3,
//        -hf*sr5*sr6/sr7,-sr7*sr3*sr5*sr2,sr7*sr3*sr5/sr2/hf,-sr7*sr3*sr5*sr2,
//
//        /*20: T2-5*/
//        0.,0.,0.,0.,
//        sr5*sr2*hf/sr3/sr3,sr5*sr2*hf/sr3/sr3,-sr2*sr5*sr3*sr3,-sr2*sr5*sr3*sr3,
//        -sr6*sr6*sr2/sr5,-sr6*sr6*sr2/sr5,sr3*sr3*sr5/sr2,sr3*sr3*sr5/sr2,
//        -sr3*sr3*sr2*sr5, -sr3*sr3*sr2*sr5,-sr6*sr6*sr2/sr5,-sr6*sr6*sr2/sr5,
//        sr3*sr3*sr5/sr2,sr3*sr3*sr5/sr2,sr6*sr6*sr2*sr5,sr6*sr6*sr2*sr5,
//        -sr3*sr3*sr2*sr5, -sr3*sr3*sr2*sr5,-sr3*sr3*sr2*sr5, -sr3*sr3*sr2*sr5,
//
//
//        /*21: T2-6*/
//        0.,0.,0.,0.,
//        0.,0.,sr3*sr3/sr2,sr3*sr3/sr2,
//        sr6*sr6*sr2,sr6*sr6*sr2,sr6*sr6*sr2,sr6*sr6*sr2,
//        sr6*sr6*sr2,sr6*sr6*sr2,-sr3*sr3*sr2,-sr3*sr3*sr2,
//        -sr3*sr3*sr2,-sr3*sr3*sr2,sr6*sr6*sr2,sr6*sr6*sr2,
//        -sr3*sr3*sr2,-sr3*sr3*sr2,-sr3*sr3*sr2,-sr3*sr3*sr2,
//
//        /*22: T2-7*/
//        0.,0.,0.,0.,
//        0.,0.,0.,0.,
//        hf*sr6/sr5,hf*sr6*sr5,hf*sr6/sr5,hf*sr5*sr6,
//        -sr2*sr5/sr3*hf, -sr2*sr5/sr3*hf,sr5*sr6,-sr5*sr6,
//        sr5*sr6,-sr5*sr6, -sr2*sr5/sr3*hf, -sr2*sr5/sr3*hf,
//        sr5*sr6,-sr5*sr6,sr5*sr6,-sr5*sr6,
//
//        /*23: T2-8*/
//        0.,0.,0.,0.,
//        0.,0.,0.,0.,
//        0.,sr5,0.,sr5,
//        -hf*sr5,-hf*sr5,-hf*sr5,hf*sr5,
//        -hf*sr5,hf*sr5,-hf*sr5,-hf*sr5,
//        -hf*sr5,hf*sr5,-hf*sr5,hf*sr5,
//
//        /*24: T2-9*/
//        0.,0.,0.,0.,
//        0.,0.,0.,0.,
//        0.,0.,0.,0.,
//        hf*sr3,hf*sr3,hf*sr3,hf*sr3,
//        hf*sr3,hf*sr3,-hf*sr3,-hf*sr3,
//        -hf*sr3,-hf*sr3,-hf*sr3,-hf*sr3
//
//        //////////////////////////////////
//        /////////////////////////////////
//        ////////////////////////////////
//        ///////////////////////////////
//        //////////////////////////////
//
//
//    };
//
//    if ( bd == one ){
//        nsyp = 1;
//        msyp = 1;
//    }
//    else    if ( bd == two ){
//        nsyp = 2;
//        msyp = 2;
//    }
//    else if ( bd == three ){
//        nsyp = 6;
//        msyp = 6;
//    }else if ( bd == four ){
//        nsyp = 24;
//        msyp = 24;
//    }
//    else {
//        printf("bod\n");
//        exit(0);
//    }
//
//
//    if ( type <= 0 || type > msyp )
//    {
//        printf("he\n %lld", type);
//        exit(0);
//    }
//    if ( perm < 0 || perm >= nsyp ){
//        printf("hm\n");
//        exit(0);
//    }
//
//
//    if ( bd == two ){
//        return syp2[(type-1)*nsyp+perm];
//    }
//    else if ( bd == three ){
//        return syp3[(type-1)*nsyp+perm];
//    }
//    else if ( bd == four ){
//        return syp4[(type-1)*nsyp+perm];
//    }
//    return 0.;
//};

INT_TYPE tClassifyComponents( struct sinc_label  f1 , double * up, double * entropy){
    
    if ( bodies(f1,eigenVectors ) == one ){
        return 1;
    }
    double entr,sum;
    INT_TYPE nGroup=0,xt,irrep;
    
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
    if ( entr < f1.rt->maxEntropy ){
//        printf("**\n");

        return xt;

    }
//    printf("\n");

    return 0;
}

INT_TYPE tClassify(INT_TYPE rank, struct sinc_label  f1 , enum division label){
    double up[48],entropy;
    if ( bodies(f1, label ) == one ){
        f1.tulip[label].value.symmetry = 0;
        f1.tulip[label].value.value2 = 0;
        
        return 0;
        
    }
    INT_TYPE i,irrep;
    for ( i = 0; i < 48 ; i++)
        up[i] = 0.;
    tTabulateInnerProjection(rank, f1, label, up);
    irrep =  tClassifyComponents(f1, up,&entropy);
    f1.tulip[label].value.value2 =entropy;

    return irrep;
}





INT_TYPE tBuildIrr ( INT_TYPE rank, struct sinc_label  f1, char meta , enum division origin, INT_TYPE ospin, enum division targ , INT_TYPE tspin){
    INT_TYPE irrep;

    if ( meta == 0 || bodies(f1,origin ) == one ){
        tEqua(f1, targ, tspin, origin, ospin);
        return 0;
    }
    if (! CanonicalRank(f1, origin, ospin ))
        return 0;
    
    INT_TYPE i,nPerm=0;
    char train[24];
    if ( bodies(f1, origin ) == two ){
        nPerm = 2;
        train[0] = 'N';
        train[1] = 'T';
    }
    if ( bodies(f1, origin ) == three ){
        train[0] = 'N';
        train[1] = 'A';
        train[2] = 'B';
        train[3] = 'C';
        train[4] = 'D';
        train[5] = 'E';
        nPerm = 6;
    }
    else if ( bodies(f1, origin ) == four ){
        //

        train[0] = 'N';
        train[1] = 'a';
        train[2] = 'b';
        train[3] = 'c';
        train[4] = 'd';
        train[5] = 'e';
        train[6] = 'f';
        train[7] = 'g';
        train[8] = 'h';
        train[9] = 'i';
        train[10] = 'j';
        train[11] = 'k';
        train[12] = 'l';
        train[13] = 'm';
        train[14] = 'n';
        train[15] = 'o';
        train[16] = 'p';
        train[17] = 'q';
        train[18] = 'r';
        train[19] = 's';
        train[20] = 't';
        train[21] = 'u';
        train[22] = 'v';
        train[23] = 'w';
        nPerm = 24;
    }
    if ( CanonicalRank(f1, origin, ospin)*nPerm > part(f1,targ)){
        printf("Irr\n");
        printf("%d %d %d\n", CanonicalRank(f1, origin, ospin), part(f1, origin), part(f1, targ));
        exit(0);
    }
    
        for ( i = 0; i < nPerm ; i++){
            irrep = meta;
            f1.tulip[permutation2Vector].Current[rank] = 0;
            tPermute(rank,f1, train[i], origin, ospin, permutation2Vector, rank);
            tScaleOne(f1, permutation2Vector, rank, tGetProjection(bodies(f1, origin), irrep, i));
            tAddTw(f1, targ, tspin, permutation2Vector, rank);
        }
    
    return 0;
}

char matrixAction ( enum bodyType bd, enum block bk, char direction){
    
    //action on right vector...direction = 1
    //action on product vector ... direction = -1
    
    if ( bd == two ){
        switch (bk){
            case tv1 :
                if ( direction == 1 )
                    return 'T';
                else
                    return 'T';
            case iv1 :
                if ( direction == 1 )
                    return 'T';
                else
                    return 'T';

            case tv2 :
                if ( direction == 1 )
                    return 'N';
                else
                    return 'N';
                
            case iv2 :
                if ( direction == 1 )
                    return 'N';
                else
                    return 'N';

            case e12:
                return 'N';
                
        }
        
        
    } else if ( bd == three ){
        //            train[0] = 'T';//(1)            123
        //            train[1] = 'A';//(123)          231
        //            train[2] = 'B';//(123).(123)    312
        //            train[3] = 'C';//(12)           213
        //            train[4] = 'D';//(13)           321
        //            train[5] = 'E';//(23)           132

        switch ( bk){
            case tv1:
                return 'N';
            case tv2:
                return 'C';
            case tv3 :
                return 'D';
            case iv1:
                return 'N';
            case iv2:
                return 'C';
            case iv3 :
                return 'D';
            case e12 :
                return 'N';
            case e13 :
                return 'E';
            case e23 :
                if ( direction == 1)
                    return 'A';
                else
                    return 'B';
        }
    }
    else if ( bd == four ){
        // 0,0,0,0 : i,j,k,l   :: T 0
        // 0,2,1,2 : i,j,l,k   :: 'a'24
        // 1,3,1,2 : i,k,j,l   :: 'b'20
        // 1,1,0,0 : i,k,l,j   :: 'c'4
        // 0,3,1,0 : i,l,j,k   :: 'd'12
        // 1,2,1,3 : i,l,k,j   :: 'e'22
        // 1,0,0,0 : j,i,k,l   :: 'f'2
        // 1,2,1,2 : j,i,l,k   :: 'g' 23
        // 0,3,1,2 : j,k,i,l   :: 'h'19
        // 0,1,0,0 : j,k,l,i   :: 'i'3
        // 1,3,1,0 : j,l,i,k   :: 'j'13
        // 0,2,1,3 : j,l,k,i   :: 'k'21
        // 0,2,1,1 : k,i,j,l   :: 'l'15
        // 1,1,1,0 : k,i,l,j   :: 'm'10
        // 1,2,1,1 : k,j,i,l   :: 'n'16
        // 0,1,1,0 : k,j,l,i   :: 'o'9
        // 0,2,0,0 : k,l,i,j   :: 'p'5
        // 1,1,0,1 : k,l,j,i   :: 'q'14
        // 0,3,0,0 : l,i,j,k   :: 'r'7
        // 1,3,1,1 : l,i,k,j   :: 's'18
        // 1,3,0,0 : l,j,i,k   :: 't'8
        // 0,3,1,1 : l,j,k,i   :: 'u'17
        // 0,2,1,0 : l,k,i,j   :: 'v'6
        // 1,2,1,0 : l,k,j,i   :: 'w'11
        switch ( bk ){
            case tv1:
                return 'N';
            case tv2:
                return 'f';
            case tv3:
                return 'n';
            case tv4:
                return 'u';
            case iv1:
                return 'N';
            case iv2:
                return 'f';
            case iv3:
                return 'n';
            case iv4:
                return 'u';
            case e12:
                return 'N';
            case e13:
                return 'b';
            case e14:
                return 'e';
            case e23:
                if ( direction == 1)
                    return 'i';
                else
                    return 'r';
            case e24:
                if ( direction == 1 )
                    return 'j';
                else
                    return 'm';
            case e34:
                return 'p';
                
        }
    }

    return 'N';
}

INT_TYPE tPermuteOne(INT_TYPE rank, struct sinc_label  f1, INT_TYPE dim, char leftChar , enum division left, INT_TYPE l, INT_TYPE lspin, enum division equals, INT_TYPE espin){
    INT_TYPE LN2 = alloc(f1,left,dim);
    INT_TYPE cur = CanonicalRank(f1, equals, espin);
    INT_TYPE N1 = vector1Len(f1, dim),flagTranspose,flagTranspose2,flagTranspose3,flagTranspose4,space,A,B,AA,BB;
    double *array[6];
    {
        INT_TYPE bs;
        for ( bs = 0; bs < 6 ; bs++)
            array[bs] = myStreams( f1, tensorBuffers+bs, rank);
    }
    double * pleft;

    if ( part(f1, equals ) < CanonicalRank(f1, equals, espin ) ){
        printf("too small %d\n",equals);
        exit(0);
    }
    if ( bodies(f1, eigenVectors) == one ){
            cblas_dcopy(LN2, streams(f1, left, lspin, dim)+l*LN2, 1, streams(f1, equals, espin, dim)+cur*LN2, 1);
    }
    
    if ( bodies(f1, eigenVectors) == two)
    {
        if ( leftChar == 'N' ){//the Transpose is IMPLICIT!!!!
            cblas_dcopy(LN2, streams(f1, left, lspin, dim)+l*LN2, 1, streams(f1, equals, espin, dim)+cur*LN2, 1);
        }else {
            transpose(N1, N1,streams(f1, left, lspin, dim)+l*LN2,streams(f1, equals, espin, dim)+cur*LN2);
        }
    }
    
    
    else if ( bodies(f1, eigenVectors) == three)
    {

        
        //
        //            train[0] = 'T';//(1)            123
        //            train[1] = 'A';//(123)          231
        //            train[2] = 'B';//(123).(123)    312
        //            train[3] = 'C';//(12)           213
        //            train[4] = 'D';//(13)           321
        //            train[5] = 'E';//(23)           132
        
        flagTranspose3 = 0;
      //  for ( space = 0; space < SPACE ; space++){
            A = N1*N1;
            B = N1;
       // }
        
        if ( leftChar == 'N' ){
            flagTranspose = 0;
            flagTranspose2 = 0;
            cblas_dcopy(LN2, streams(f1, left, lspin, dim)+l*LN2, 1, streams(f1, equals, espin, dim)+cur*LN2, 1);
            return 0;
        }
        else  if ( leftChar == 'A' ){
            flagTranspose = 1;// a | b c
            
            flagTranspose2 = 0;
                A = N1;
                B = N1*N1;
        }else  if ( leftChar == 'B' ){
            flagTranspose = 1; // a b | c
            flagTranspose2 = 0;
                A = N1*N1;
                B = N1;
        } else if ( leftChar == 'C' ){
            flagTranspose = 0;
            flagTranspose2 = 1;
        } else  if ( leftChar == 'D' ){
            flagTranspose = 1; // a b | c
            flagTranspose2 = 1;
                A = N1*N1;
                B = N1;
            
        }else  if ( leftChar == 'E' ){
            flagTranspose = 1; // a b | c
            flagTranspose2 = 0;
            flagTranspose3 = 1;
            A = N1*N1;
                B = N1;
        }
        
        else {
            printf("unknown flag %c\n",leftChar);
            exit(0);
        }
        
        
    
            INT_TYPE bs,o;
    
                {
                    space = dim;
            //
                    pleft = streams( f1, left, lspin,space )+l*LN2;
                    
                    bs = 0;
                    
                    if (flagTranspose2){
                        for ( o = 0 ; o < N1 ;o++){
                            transpose(N1, N1,pleft+N1*N1*o,array[bs]+N1*N1*o);
                            
                        }
                        pleft = array[bs++];
                    }
                    if ( flagTranspose  ){
                        transpose(A, B,pleft ,array[bs]);
                        
                        pleft = array[bs++];
                        
                    }
                    
                    if ( flagTranspose3){
                        for ( o = 0 ; o < N1 ;o++){
                            transpose(N1, N1,pleft+N1*N1*o,array[bs]+N1*N1*o);
                        }
                        pleft = array[bs++];
                        
                    }
                    cblas_dcopy(LN2, pleft, 1, streams(f1,equals, espin, space)+cur*LN2 ,1);
                    
                }
            }
            
            
            
        
        
        

    else if (bodies ( f1, eigenVectors ) == four ){
        {
            
                flagTranspose = 0;
                flagTranspose2 = 0;
                flagTranspose3 = 0;
                flagTranspose4 = 0;
                
                    A = N1*N1;
                    B = N1*N1;
                
                    AA = N1;
                BB = N1*N1*N1;
            
                
                
                // 0,0,0,0 : i,j,k,l   :: 'N' 0
                // 0,2,1,2 : i,j,l,k   :: 'a'24
                // 1,3,1,2 : i,k,j,l   :: 'b'20
                // 1,1,0,0 : i,k,l,j   :: 'c'4
                // 0,3,1,0 : i,l,j,k   :: 'd'12
                // 1,2,1,3 : i,l,k,j   :: 'e'22
                // 1,0,0,0 : j,i,k,l   :: 'f'2
                // 1,2,1,2 : j,i,l,k   :: 'g' 23
                // 0,3,1,2 : j,k,i,l   :: 'h'19
                // 0,1,0,0 : j,k,l,i   :: 'i'3
                // 1,3,1,0 : j,l,i,k   :: 'j'13
                // 0,2,1,3 : j,l,k,i   :: 'k'21
                // 0,2,1,1 : k,i,j,l   :: 'l'15
                // 1,1,1,0 : k,i,l,j   :: 'm'10
                // 1,2,1,1 : k,j,i,l   :: 'n'16
                // 0,1,1,0 : k,j,l,i   :: 'o'9
                // 0,2,0,0 : k,l,i,j   :: 'p'5
                // 1,1,0,1 : k,l,j,i   :: 'q'14
                // 0,3,0,0 : l,i,j,k   :: 'r'7
                // 1,3,1,1 : l,i,k,j   :: 's'18
                // 1,3,0,0 : l,j,i,k   :: 't'8
                // 0,3,1,1 : l,j,k,i   :: 'u'17
                // 0,2,1,0 : l,k,i,j   :: 'v'6
                // 1,2,1,0 : l,k,j,i   :: 'w'11
                
                
                
                if ( leftChar == 'N' ){//(i,j,k,l)
                    cblas_dcopy(LN2, streams(f1, left, lspin, dim)+l*LN2, 1, streams(f1, equals, espin, dim)+cur*LN2, 1);
                    return 0;

                } else if ( leftChar == 'a' ){//(i,j,l,k)
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    A = N1*N1;
                        B = N1*N1;
                        AA = N1*N1;
                        BB = N1*N1;
                    
                }else if ( leftChar == 'b' ){//(i,j,k,l)->(j,i,k,l)->(l,j,i,k)->(j,l,i,k)->(i,k,j,l)->(i,k,j,l)
                    flagTranspose = 1;
                    flagTranspose2 = 1;
                        A = N1*N1*N1;
                        B = N1;
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                        AA = N1*N1;
                        BB = N1*N1;
                    
                    
                }else  if ( leftChar == 'c' ){//(i,k,l,j)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                        A = N1;
                        B = N1*N1*N1;
                }else  if ( leftChar == 'd' ){//(i,l,j,k)
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    
                        A = N1*N1*N1;
                        B = N1;
                }else  if ( leftChar == 'e' ){//(i,l,k,j) 6
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                        A = N1*N1;
                        B = N1*N1;
                    
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    
                        AA = N1*N1*N1;
                        BB = N1;
                    
                }else  if ( leftChar == 'f' ){//(j,i,k,l)
                    flagTranspose = 1; // a b | c
                }else  if ( leftChar == 'g' ){//(j,i,l,k)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                        A = N1*N1;
                        B = N1*N1;
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                        AA = N1*N1;
                        BB = N1*N1;
                    
                }else  if ( leftChar == 'h' ){//(j,k,i,l)
                    flagTranspose = 0; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A = N1*N1*N1;
                        B = N1;
                    }
                    
                    flagTranspose3 = 1; // a b | c
                    flagTranspose4 = 1;
                        AA= N1*N1;
                        BB = N1*N1;
                    
                }else  if ( leftChar == 'i' ){//(j,k,l,i)
                    flagTranspose = 0; // a b | c
                    flagTranspose2 = 1;
                    
                        A = N1;
                        B = N1*N1*N1;
                }
                else  if ( leftChar == 'j' ){//(j,l,i,k)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                        A = N1*N1*N1;
                    B = N1;
                    
                }else  if ( leftChar == 'k' ){//(j,l,k,i)12
                    flagTranspose2 = 2;
                        A = N1*N1;
                        B = N1*N1;
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                        AA = N1*N1*N1;
                        BB = N1;
                    
                }else  if ( leftChar == 'l' ){//(k,i,j,l)
                    flagTranspose = 0; // a b | c
                    flagTranspose2 = 1;
                    
                        A = N1*N1;
                        B = N1*N1;
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    
                        AA = N1;
                        BB = N1*N1*N1;
                    
                }else  if ( leftChar == 'm' ){//(k,i,l,j)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    
                        A = N1;
                        B = N1*N1*N1;
                }else  if ( leftChar == 'n' ){//(k,j,i,l)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                        A = N1*N1;
                        B = N1*N1;
                    
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                        AA = N1;
                        BB = N1*N1*N1;
                    
                    
                }else  if ( leftChar == 'o' ){//(k,j,l,i)
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                        A = N1;
                        B = N1*N1*N1;
                }else  if ( leftChar == 'p' ){//(k,l,i,j)
                    flagTranspose2 = 1;
                    
                        A = N1*N1;
                        B = N1*N1;
                }else  if ( leftChar == 'q' ){//i,j,k,l->j,i,k,l->i,k,l,j->k,l,j,i->(k,l,j,i)//18
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                        A = N1;
                        B = N1*N1*N1;
                    
                    flagTranspose4 = 1;
                    
                        AA = N1;
                        BB = N1*N1*N1;
            
                }
                else  if ( leftChar == 'r' ){//(l,i,j,k)
                    flagTranspose2 = 1;
                    
                        A = N1*N1*N1;
                        B = N1;
                    
                }
                else  if ( leftChar == 's' ){//(l,i,k,j)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                        A = N1*N1*N1;
                        B = N1;
                    flagTranspose3 = 1; // a b | c
                    flagTranspose4 = 1;
                        AA = N1;
                        BB = N1*N1*N1;
                    
                }
                else  if ( leftChar == 't' ){//(l,j,i,k)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                        A = N1*N1*N1;
                        B = N1;
                    
                }
                else  if ( leftChar == 'u' ){//(l,j,k,i)
                    flagTranspose2 = 1;
                    
                        A = N1*N1*N1;
                        B = N1;
                    
                    
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                        AA = N1;
                        BB = N1*N1*N1;
                    
                }
                // 0,2,1,0 : l,k,i,j   :: 'v'6
                else  if ( leftChar == 'v' ){//i,j,k,l->j,i,k,l->k,l,j,i->(l,k,i,j)
                    flagTranspose2 = 1;
                    
                        A = N1*N1;
                        B = N1*N1;
                    
                    flagTranspose3 = 1;
                    
                }
                else  if ( leftChar == 'w' ){//(l,k,j,i)//END
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    
                        A = N1*N1;
                        B = N1*N1;
                
                }
                
                
                
                
                else {
                    printf("unknown flag %c\n",leftChar);
                    exit(0);
                }
            
                {
                    INT_TYPE bs,o;
                    {
                        {
                            space = dim;
                            bs = 0;
                            //
                            pleft = streams( f1, left, lspin,space )+l*LN2;
                            
                            
                            if (flagTranspose){
                                for ( o = 0 ; o < N1*N1 ;o++)
                                    transpose(N1, N1,pleft+N1*N1*o,array[bs]+N1*N1*o);
                                pleft = array[bs++];
                            }
                            
                            if ( flagTranspose2  ){
                                transpose(A, B,pleft,array[bs]);
                                pleft = array[bs++];
                            }
                            
                            if ( flagTranspose3){
                                for ( o = 0 ; o < N1*N1 ;o++)
                                    transpose(N1, N1,pleft+N1*N1*o ,array[bs]+N1*N1*o);
                                pleft = array[bs++];
                                
                                
                            }
                            
                            if ( flagTranspose4  ){
                                transpose(AA, BB,pleft,array[bs]);
                                pleft= array[bs++];
                                
                            }
                            cblas_dcopy(LN2, pleft, 1, streams(f1,equals, espin, space)+cur*LN2 ,1);
                            
                        }
                        
                        
                        
                    }
                    
                }
                

                
                
        }

    }
    
    return 0;
}


INT_TYPE tPermute(INT_TYPE rank, struct sinc_label f1, char leftChar , enum division left, INT_TYPE lspin, enum division equals, INT_TYPE espin){
    INT_TYPE l,space;
    
    for ( l = 0; l < CanonicalRank(f1,left,lspin); l++){
        for ( space = 0; space < SPACE ; space++)
            tPermuteOne(rank, f1, space, leftChar, left, l, lspin, equals, espin);
        f1.tulip[equals].Current[espin]++;
    }

    return 0;
}




INT_TYPE tAllCompPermMultiplyMP( INT_TYPE rank, struct sinc_label  f1 , enum division left ,INT_TYPE lspin, enum division right ,INT_TYPE rspin, double * sequ){
    
    if ( CanonicalRank(f1, left, lspin ) * CanonicalRank(f1, right, rspin ) == 0)
        return 0;
    
    if ( bodies(f1, left ) != bodies(f1,right)){
        printf("tGetType real!\n");
        exit(0);
    }
    char train[24];
    INT_TYPE nPerm,i,info ;
    nPerm = 1;

    if ( bodies(f1,left) == one )
        return 0;
    else
    if ( bodies(f1, left ) == two ){
        //
        train[0] = 'N';
        train[1] = 'T';//(12)
        nPerm = 2;
    }
    else if ( bodies(f1, left ) == three ){
        //
        train[0] = 'N';//(1)            123
        train[1] = 'A';//(123)          231
        train[2] = 'B';//(123).(123)    312
        train[3] = 'C';//(12)           213
        train[4] = 'D';//(13)           321
        train[5] = 'E';//(23)           132
        nPerm = 6;
    }
    else if ( bodies(f1, left ) == four ){
        //
        train[0] = 'N';
        train[1] = 'a';
        train[2] = 'b';
        train[3] = 'c';
        train[4] = 'd';
        train[5] = 'e';
        train[6] = 'f';
        train[7] = 'g';
        train[8] = 'h';
        train[9] = 'i';
        train[10] = 'j';
        train[11] = 'k';
        train[12] = 'l';
        train[13] = 'm';
        train[14] = 'n';
        train[15] = 'o';
        train[16] = 'p';
        train[17] = 'q';
        train[18] = 'r';
        train[19] = 's';
        train[20] = 't';
        train[21] = 'u';
        train[22] = 'v';
        train[23] = 'w';
        nPerm = 24;
    }

    for ( i = 0; i < nPerm ; i++){
        sequ[i] = tMultiplyMP(rank, &info,f1,1. , -1, nullVector, 0, train[i], left, lspin, 'N', right, rspin);
    }
    
    return nPerm;
}


INT_TYPE tTabulateInnerProjection( INT_TYPE rank, struct sinc_label  f1 , enum division vec, double *up){

    if ( bodies(f1,vec ) == one ){
        up[0] = 1;
        return 1;
    }
    INT_TYPE g,p,nPerm=0,nGroup=0;
    
    
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
    double buff[720];

    tAllCompPermMultiplyMP(rank, f1, vec, 0, vec,0, buff);
    for ( g = 1; g <= nGroup ; g++)
        for ( p = 0; p < nPerm ; p++)
            up[g]  += tGetProjection(bodies(f1,vec), g, p)*buff[p];

    
    tAllCompPermMultiplyMP(rank, f1, vec, 1, vec,1, buff);
    for ( g = 1; g <= nGroup ; g++)
        for ( p = 0; p < nPerm ; p++)
            up[g]  += tGetProjection(bodies(f1,vec), g, p)*buff[p];

    return nGroup;
}



//INT_TYPE tTest ( enum bodyType bd ){
//    INT_TYPE a,b,nX=0;
//
//    if ( bd == two )
//        nX = 4;
//    if ( bd == three )
//        nX = 9;
//    if ( bd == four )
//        nX = 29;
//
//
//    for ( a = 1  ; a <= nX  ; a++)
//        for ( b = 1  ; b <= nX  ; b++)
//            printf("%d\t%d\t%f\n", a,b, tGetInner(bd, a, b));
////    if ( bd == three )
////        for ( a = 1 ; a <= 3 ; a++)
////            for ( b = 0 ; b < 6 ; b++){
////                if ( b % 3 == 0 )
////                    printf("\n");
////
////                printf("%f,", tGetIrrep(bd, a, b));
////
////            }
////
////
////    if ( bd == four )
////    for ( a = 1 ; a <= 5 ; a++)
////        for ( b = 0 ; b < 24 ; b++){
////            if ( b % 4 == 0 )
////                printf("\n");
////
////            printf("%f,", tGetIrrep(bd, a, b));
////
////        }
//
//    return 0;
//}


INT_TYPE tSize(enum bodyType bd){
    
    INT_TYPE nG;
    if ( bd == two )
        nG = 2;
    else if ( bd == three )
        nG = 3;
    else  if ( bd == four )
        nG = 5;
    else
        nG = 1;
    return nG;
}

INT_TYPE tPerms(enum bodyType bd){
    
    INT_TYPE nP;
    if ( bd == two )
        nP = 2;
    else if ( bd == three )
        nP = 6;
    else  if ( bd == four )
        nP = 24;
    else
        nP = 1;
    return nP;
}
