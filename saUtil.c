/*
 *  saUtil.c
 *
 *
 *  Copyright 2018 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University
 *  and the Robert A. Welch Foundation.
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

INT_TYPE tFil ( struct field * f1, enum division A, enum division v , INT_TYPE * i ){
    INT_TYPE n1[3],space;
    n1[0] = vectorLen(f1,A)[0];
    n1[1] = vectorLen(f1,A)[1];
    n1[2] = vectorLen(f1,A)[2];
    if ( bodies (f1, v ) == two ){
        for ( space = 0; space < SPACE ; space++){
            cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space]*n1[space],1, streams(f1,diagonal1VectorA,0,space),1);
        }
        for ( space = 0; space < SPACE ; space++){
            cblas_dcopy(n1[space], streams(f1,A,0,space)+i[space+3]*n1[space],1, streams(f1,diagonal1VectorB,0,space),1);
        }
        f1->sinc.tulip[diagonal1VectorA].Current[0] = 1;
        f1->sinc.tulip[v].Current[0] = 0;
        
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
            f1->sinc.tulip[diagonal1VectorA].Current[0] = 1;
            f1->sinc.tulip[diagonal1VectorB].Current[0] = 1;
            f1->sinc.tulip[diagonal1VectorC].Current[0] = 1;
            f1->sinc.tulip[diagonal2VectorA].Current[0] = 0;
            f1->sinc.tulip[v].Current[0] = 0;
            
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
            f1->sinc.tulip[diagonal1VectorA].Current[0] = 1;
            f1->sinc.tulip[diagonal1VectorB].Current[0] = 1;
            f1->sinc.tulip[diagonal1VectorC].Current[0] = 1;
            f1->sinc.tulip[diagonal1VectorD].Current[0] = 1;
            f1->sinc.tulip[diagonal2VectorA].Current[0] = 0;
            f1->sinc.tulip[diagonal2VectorB].Current[0] = 0;
            f1->sinc.tulip[v].Current[0] = 0;
            
            
            tOuterProductSu(f1, diagonal1VectorA, 0, diagonal1VectorB, 0, diagonal2VectorA, 0);
            tOuterProductSu(f1, diagonal1VectorC, 0, diagonal1VectorD, 0, diagonal2VectorB, 0);
            tOuterProductSu(f1, diagonal2VectorA, 0, diagonal2VectorB, 0, v, 0);
            
        }
    
    return 0;
}




INT_TYPE tInnerTest( struct field * f1, enum division A ,enum division B){
    char c;
    INT_TYPE n1[3],space,i[100],j[100],nPerm,ii;
    double seq[100];
    n1[0] = vectorLen(f1,A)[0];
    n1[1] = vectorLen(f1,A)[1];
    n1[2] = vectorLen(f1,A)[2];
    if ( bodies ( f1, A ) != one || bodies (f1, B ) != one ){
        printf ("ment for onebody\n");
        exit(0);
    }
    
    for ( space = 0; space < SPACE ; space++)
        tdsyev (0,f1,'V',n1[space],streams(f1,A,0,space),n1[space],streams(f1,B,0,space));
    for ( ii = 0; ii < 100 ;ii++)
        i[ii] = ii/3 % n1[0];
    
    tFil(f1, A, copyThreeVector, i);
    
    
    INT_TYPE jj;
    double ma [ 24 * 24];
    
    
    
    if ( bodies ( f1, eigenVectors )== three ){
        for ( c = 0; c < 5 ; c++){
            printf("\n\nCHAR %c\n", 'A'+c);
            tPermute(0, f1, 'A'+c, copyThreeVector, 0, copyFourVector, 0);
            nPerm = tAllCompPermMultiplyMP(0, f1, copyThreeVector, 0, copyFourVector, 0, seq);
            for ( ii = 0 ; ii < nPerm; ii++)
                printf("%lld\t%f\n",ii+1,seq[ii]);
        }
    }else
        if ( bodies ( f1, eigenVectors )== four ){
            for ( c = 0; c < 23 ; c++){
                printf("\n\nCHAR %c\n", 'a'+c);
                tPermute(0, f1, 'a'+c, copyThreeVector, 0, copyFourVector, 0);
                nPerm = tAllCompPermMultiplyMP(0, f1, copyThreeVector, 0, copyFourVector, 0, seq);
                for ( ii = 0 ; ii < nPerm; ii++)
                    printf("%lld\t%f\n",ii+1,seq[ii]);
            }
        }
    
    
    if ( bodies ( f1, eigenVectors )== three ){
        for ( c = 'A'; c <= 'E' ; c++){
            if ( c == 'A' ){
                //            train[0] = 'T';//(1)            123
                //            train[1] = 'A';//(123)          231
                //            train[2] = 'B';//(123).(123)    312
                //            train[3] = 'C';//(12)           213
                //            train[4] = 'D';//(13)           321
                //            train[5] = 'E';//(23)           132
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+3];
                    j[space+3] = i[space+6];
                    j[space+6] = i[space];
                    
                }
            }else
                if ( c == 'B' ){
                    //            train[0] = 'T';//(1)            123
                    //            train[1] = 'A';//(123)          231
                    //            train[2] = 'B';//(123).(123)    312
                    //            train[3] = 'C';//(12)           213
                    //            train[4] = 'D';//(13)           321
                    //            train[5] = 'E';//(23)           132
                    for ( space = 0 ; space < SPACE ; space++){
                        j[space] = i[space+6];
                        j[space+3] = i[space];
                        j[space+6] = i[space+3];
                        
                    }
                }
                else            if ( c == 'C' ){
                    //            train[0] = 'T';//(1)            123
                    //            train[1] = 'A';//(123)          231
                    //            train[2] = 'B';//(123).(123)    312
                    //            train[3] = 'C';//(12)           213
                    //            train[4] = 'D';//(13)           321
                    //            train[5] = 'E';//(23)           132
                    for ( space = 0 ; space < SPACE ; space++){
                        j[space] = i[space+3];
                        j[space+3] = i[space];
                        j[space+6] = i[space+6];
                        
                    }
                }else             if ( c == 'D' ){
                    //            train[0] = 'T';//(1)            123
                    //            train[1] = 'A';//(123)          231
                    //            train[2] = 'B';//(123).(123)    312
                    //            train[3] = 'C';//(12)           213
                    //            train[4] = 'D';//(13)           321
                    //            train[5] = 'E';//(23)           132
                    for ( space = 0 ; space < SPACE ; space++){
                        j[space] = i[space+6];
                        j[space+3] = i[space+3];
                        j[space+6] = i[space];
                        
                    }
                }else             if ( c == 'E' ){
                    //            train[0] = 'T';//(1)            123
                    //            train[1] = 'A';//(123)          231
                    //            train[2] = 'B';//(123).(123)    312
                    //            train[3] = 'C';//(12)           213
                    //            train[4] = 'D';//(13)           321
                    //            train[5] = 'E';//(23)           132
                    for ( space = 0 ; space < SPACE ; space++){
                        j[space] = i[space];
                        j[space+3] = i[space+6];
                        j[space+6] = i[space+3];
                        
                    }
                }
            
            
            
            
            
            
            tFil(f1, A, copyFourVector, j);
            nPerm = tAllCompPermMultiplyMP(0, f1, copyThreeVector, 0, copyFourVector, 0, seq);
            for ( ii = 0 ; ii < nPerm; ii++)
                printf("%lld\t%f\n",ii+1,seq[ii]);
            
        }
    }
    if ( bodies ( f1, eigenVectors )== four ){
        for ( c = 'a'; c <= 'w' ; c++){
            
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
            // 1,2,0,0 : l,k,i,j   :: 'v'6
            // 1,2,1,0 : l,k,j,i   :: 'w'11
            if ( c == 'a' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space];
                    j[space+3] = i[space+3];
                    j[space+6] = i[space+9];
                    j[space+9] = i[space+6];
                }
            }else            if ( c == 'b' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space];
                    j[space+3] = i[space+6];
                    j[space+6] = i[space+3];
                    j[space+9] = i[space+9];
                }
            }else            if ( c == 'c' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space];
                    j[space+3] = i[space+6];
                    j[space+6] = i[space+9];
                    j[space+9] = i[space+3];
                }
            }else            if ( c == 'd' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space];
                    j[space+3] = i[space+9];
                    j[space+6] = i[space+3];
                    j[space+9] = i[space+6];
                }
            }else            if ( c == 'e' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space];
                    j[space+3] = i[space+9];
                    j[space+6] = i[space+6];
                    j[space+9] = i[space+3];
                }
            }else            if ( c == 'f' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+3];
                    j[space+3] = i[space];
                    j[space+6] = i[space+6];
                    j[space+9] = i[space+9];
                }
            }else            if ( c == 'g' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+3];
                    j[space+3] = i[space];
                    j[space+6] = i[space+9];
                    j[space+9] = i[space+6];
                }
            }else            if ( c == 'h' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+3];
                    j[space+3] = i[space+6];
                    j[space+6] = i[space];
                    j[space+9] = i[space+9];
                }
            }else            if ( c == 'i' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+3];
                    j[space+3] = i[space+6];
                    j[space+6] = i[space+9];
                    j[space+9] = i[space];
                }
            }else            if ( c == 'j' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+3];
                    j[space+3] = i[space+9];
                    j[space+6] = i[space];
                    j[space+9] = i[space+6];
                }
            }else            if ( c == 'k' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+3];
                    j[space+3] = i[space+9];
                    j[space+6] = i[space+6];
                    j[space+9] = i[space];
                }
            }else            if ( c == 'l' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+6];
                    j[space+3] = i[space];
                    j[space+6] = i[space+3];
                    j[space+9] = i[space+9];
                }
            }else            if ( c == 'm' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+6];
                    j[space+3] = i[space];
                    j[space+6] = i[space+9];
                    j[space+9] = i[space+3];
                }
            }else            if ( c == 'n' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+6];
                    j[space+3] = i[space+3];
                    j[space+6] = i[space];
                    j[space+9] = i[space+9];
                }
            }else            if ( c == 'o' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+6];
                    j[space+3] = i[space+3];
                    j[space+6] = i[space+9];
                    j[space+9] = i[space];
                }
            }else            if ( c == 'p' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+6];
                    j[space+3] = i[space+9];
                    j[space+6] = i[space];
                    j[space+9] = i[space+3];
                }
            }else            if ( c == 'q' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+6];
                    j[space+3] = i[space+9];
                    j[space+6] = i[space+3];
                    j[space+9] = i[space];
                }
            }else            if ( c == 'r' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+9];
                    j[space+3] = i[space];
                    j[space+6] = i[space+3];
                    j[space+9] = i[space+6];
                }
            }else            if ( c == 's' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+9];
                    j[space+3] = i[space];
                    j[space+6] = i[space+6];
                    j[space+9] = i[space+3];
                }
            }else            if ( c == 't' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+9];
                    j[space+3] = i[space+3];
                    j[space+6] = i[space];
                    j[space+9] = i[space+6];
                }
            }else            if ( c == 'u' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+9];
                    j[space+3] = i[space+3];
                    j[space+6] = i[space+6];
                    j[space+9] = i[space];
                }
            }else            if ( c == 'v' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+9];
                    j[space+3] = i[space+6];
                    j[space+6] = i[space];
                    j[space+9] = i[space+3];
                }
            }else            if ( c == 'w' ){
                for ( space = 0 ; space < SPACE ; space++){
                    j[space] = i[space+9];
                    j[space+3] = i[space+6];
                    j[space+6] = i[space+3];
                    j[space+9] = i[space];
                }
            }
            
            
            
            
            
            tFil(f1, A, copyFourVector, j);
            nPerm = tAllCompPermMultiplyMP(0, f1, copyThreeVector, 0, copyFourVector, 0, seq);
            for ( ii = 0 ; ii < nPerm; ii++)
                printf("%lld\t%f\n",ii+1,seq[ii]);
            
            
            
            
        }
        
        
        
        
        
        
        
        
        
        
        
        
    }
    
    return 0;
}


double deg(struct field *f1, INT_TYPE cl ){
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


double get(enum body bd , INT_TYPE type , INT_TYPE i ){
    
    INT_TYPE nsyp=0 ,msyp=0;
    const static double syp2 [] = {
        sr2,sr2,
        sr2,-sr2,
        1.,1.,
        1.,-1.
    };
    const static double syp3 [] = {
        /********/
        sr6,sr6,sr6, sr6, sr6, sr6,
        sr6,sr6,sr6,-sr6,-sr6,-sr6,
        hf*sr3,hf*sr3, -sr3,hf*sr3,hf*sr3,-sr3,
        hf, -hf,0.,-hf, hf , 0.,
        hf*sr3,hf*sr3, -sr3,-hf*sr3,-hf*sr3,sr3,
        hf, -hf,0.,hf, -hf , 0.,
        //////////////////////////////////
        /////////////////////////////////
        ////////////////////////////////
        ///////////////////////////////
        //////////////////////////////

        /* 7:: A1 */
        1.,1.,1., 1., 1., 1.,
        
        /* 8:: A2*/
        1.,1.,1.,-1.,-1.,-1.,

        /* 9 :: E*/
        2., -1., -1.,0.,0.,0.
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
        -hf*sr3,-hf*sr3,-hf*sr3,-hf*sr3,
        
        //////////////////////////////////
        /////////////////////////////////
        ////////////////////////////////
        ///////////////////////////////
        //////////////////////////////
        /* 25:: A1*/
        1.,1.,1.,1.,
        1.,1.,1.,1.,
        1.,1.,1.,1.,
        1.,1.,1.,1.,
        1.,1.,1.,1.,
        1.,1.,1.,1.,
        
        
        /* 26:: A2*/
        1.,-1.,-1.,1.,
        1.,-1.,-1.,1.,
        -1.,-1.,1.,1.,
        1.,1.,-1.,-1.,
        1.,-1.,-1.,-1.,
        1.,-1.,1.,1.,

        /* 27:: E */
        2.,0.,0.,-1.,
        -1.,0.,0.,2.,
        0.,0.,-1.,-1.,
        -1.,-1.,0.,0.,
        2.,0.,0.,0.,
        -1.,0.,-1.,2.,

        /* 28:: T1 */
        3.,-1.,-1.,0.,
        0.,-1.,-1.,-1.,
        1.,1.,0.,0.,
        0.,0.,-1.,1.,
        -1.,1.,1.,1.,
        0.,-1.,0.,-1.,

        /* 29:: T2 */
        3.,1.,1.,0.,
        0.,1.,1.,-1.,
        -1.,-1.,0.,0.,
        0.,0.,1.,-1.,
        -1.,-1.,-1.,-1.,
        0.,1.,0.,-1.

      /*
       
       * 6 sigd
       / 6 S4
       : 3 C2
         8 C3
       ()1 E
       
            (){i, j, k, l},   *{i, j, l, k},  *{i, k, j, l},   {i, k, l, j},
            {i, l, j, k},   *{i, l, k, j},  *{j, i, k, l},   :{j, i, l, k},
            /{j, k, i, l},   /{j, k, l, i},  {j, l, i, k},   {j, l, k, i},
            {k, i, j, l},   {k, i, l, j},   *{k, j, i, l},     /{k, j, l, i},
            :{k, l, i, j},  /{k, l, j, i},   /{l, i, j,k},    /{l, i, k, j},
            {l, j, i, k},   *{l, j, k, i},   {l, k, i, j},   :{l, k,j, i}}
       */
    };

    if ( bd == one ){
        nsyp = 1;
        msyp = 1;
    }
    else    if ( bd == two ){
        nsyp = 2;
        msyp = 4;
    }
    else if ( bd == three ){
        nsyp = 6;
        msyp = 9;
    }else if ( bd == four ){
        nsyp = 24;
        msyp = 29;
    }
    else {
        printf("bod\n");
        exit(0);
    }
    
    
    if ( type <= 0 || type > msyp )
    {
        printf("he\n %lld", type);
        exit(0);
    }
    if ( i < 0 || i >= nsyp ){
        printf("hm\n");
        exit(0);
    }
    
    
    if ( bd == two ){
        return syp2[(type-1)*nsyp+i];
    }
    else if ( bd == three ){
        return syp3[(type-1)*nsyp+i];
    }
    else if ( bd == four ){
        return syp4[(type-1)*nsyp+i];
    }
    return 0.;
};
INT_TYPE tClassifyComponents( struct field * f1 , double * up, double * entropy){
    
    if ( bodies(f1,eigenVectors ) == one ){
        return 1;
    }
    double entr,sum;
    INT_TYPE nPerm=0,nGroup=0,xt,type;
    
    if ( bodies(f1, eigenVectors ) == two ){
        nGroup = 2;
    }
    else if ( bodies ( f1, eigenVectors ) == three ){
        nPerm = 6;
        nGroup = 3;//JLJ:HERE
    }
    else if ( bodies (f1, eigenVectors ) == four  ){
        nPerm = 24;
        nGroup = 5;//JLJ:HERE
    }
    xt=0;
    entr = 0.;
    sum = 0.;
    for ( type = 0 ; type < nGroup ; type++ ){
        sum += fabs(up[type]);
        if ( fabs(up[type])> fabs(up[xt]))
            xt = type;
    }
    for ( type = 0 ; type < nGroup ; type++ ){
        if ( fabs(up[type]) > 1e-6 ){
            entr += -(fabs(up[type])/sum)*log(fabs(up[type])/sum);
        }
//        printf("%1.3f,", up[type]);
    }
    *entropy = entr;
    if ( entr < f1->mem1->rt->maxEntropy ){
//        printf("**\n");

        return xt+1;

    }
//    printf("\n");

    return 0;
}

INT_TYPE tClassify(INT_TYPE rank, struct field * f1 , enum division label){
    double up[48],entropy;
    INT_TYPE i,type;
    for ( i = 0; i < 48 ; i++)
        up[i] = 0.;
    
    tTabulateProjection(rank, f1, label, label, up);
    type =  tClassifyComponents(f1, up,&entropy);
    f1->sinc.tulip[label].value2 =entropy;
    return type;
}


INT_TYPE tBuildIrr ( INT_TYPE rank, struct field * f1, char type/*1..6 (3b)*/ , enum division origin, INT_TYPE ospin, enum division targ , INT_TYPE tspin){
    INT_TYPE map[24];

    if ( type == 0 ){
        tEqua(f1, targ, tspin, origin, ospin);
        return 0;
    }
    if (! CanonicalRank(f1, origin, ospin ))
        return 0;
    
    if ( origin == diagonalVector || targ == diagonalVector ){
        printf("nop!");
        exit(0);
    }
    INT_TYPE i,nPerm=0,nDeg=0,d;
    char train[24];
    double sum;
    if ( bodies(f1, origin ) == two ){
        nPerm = 2;
        nDeg = 1;
        train[0] = 'T';
        train[1] = 'N';
        map[1] = 1;
        if ( type == 2 ){
            nDeg = 1;
            map[1] = 2;
        }
    }
    if ( bodies(f1, origin ) == three ){
        train[0] = 'T';
        train[1] = 'A';
        train[2] = 'B';
        train[3] = 'C';
        train[4] = 'D';
        train[5] = 'E';
        nPerm = 6;
        map[1] = 1;
        if ( type == 2 ){
            map[2] = 2;
            nDeg = 1;
        } else {
            nDeg = 4;
            map[1] = 3;
            map[2] = 4;
            map[3] = 5;
            map[4] = 6;
        }
    }
    else if ( bodies(f1, origin ) == four ){
        //
        nDeg = 1;
        map [1] = 1;
        if ( type == 2 ){
            map[1] = 2;
            nDeg = 1;
        } else if ( type == 3 ){
            map[1] = 3;
            map[2] = 4;
            map[3] = 5;
            map[4] = 6;
            nDeg = 4;
        } else if ( type == 4){
            map[1] = 7;
            map[2] = 8;
            map[3] = 9;
            map[4] = 10;
            map[5] = 11;
            map[6] = 12;
            map[7] = 13;
            map[8] = 14;
            map[9] = 15;
            nDeg = 9;
        }else if ( type == 5 ){
            map[1] = 16;
            map[2] = 17;
            map[3] = 18;
            map[4] = 19;
            map[5] = 20;
            map[6] = 21;
            map[7] = 22;
            map[8] = 23;
            map[9] = 24;
            nDeg = 9;
        }
        train[0] = 'T';
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
        exit(0);
    }
    for ( i = 0; i < nPerm ; i++){
        sum = 0.;
        for ( d = 1 ; d <= nDeg ; d++)
            sum += get(bodies(f1, origin ),map[d],i);
        
        
        if ( fabs(sum) > 1e-6  ){
            f1->sinc.tulip[diagonalVector].Current[rank] = 0;
            tPermute(rank,f1, train[i], origin, ospin, diagonalVector, rank);
//		printf("+%d\n", CanonicalRank(f1,diagonalVector,rank));
            tScaleOne(f1, diagonalVector, rank, sum);
            tAddTw(f1, targ, tspin, diagonalVector, rank);
//	printf( "%d %f\n", i, sum);
 
       }
    }
    
    return 0;
}


INT_TYPE tPermute(INT_TYPE rank, struct field * f1, char leftChar , enum division left, INT_TYPE lspin, enum division equals, INT_TYPE espin){
    INT_TYPE LN2[SPACE];
    length(f1, left, LN2);
    
    INT_TYPE l,flagTranspose,flagTranspose2,flagTranspose3,flagTranspose4,space,*N1, A[SPACE],B[SPACE],AA[SPACE],BB[SPACE];
    double *array[4];
    if ( part(f1, equals ) < CanonicalRank(f1, left, lspin ) ){
        printf("too small\n");
        exit(0);
    }
    
    
    if ( bodies(f1, v) == two)
    {
        N1 = f1->sinc.Basis;
        
        
        
        
        if ( leftChar == 'T' ){
            tEqua(f1, equals, espin, left, lspin );
            return 0;
        }else {
            for ( l = 0 ; l < CanonicalRank(f1, left, lspin ) ; l++){
                for ( space = 0; space < SPACE ;space++){
                    
                    transpose(N1[space], N1[space],streams(f1, left, lspin, space)+l*LN2[space],streams(f1, equals, espin, space)+l*LN2[space]);
                    
                    
                    
                }
            }
        }
    }
    
    else if ( bodies(f1, v) == three)
    {
        N1 = f1->sinc.Basis;
        
        
        //
        //            train[0] = 'T';//(1)            123
        //            train[1] = 'A';//(123)          231
        //            train[2] = 'B';//(123).(123)    312
        //            train[3] = 'C';//(12)           213
        //            train[4] = 'D';//(13)           321
        //            train[5] = 'E';//(23)           132
        
        flagTranspose3 = 0;
        for ( space = 0; space < SPACE ; space++){
            A[space] = N1[space]*N1[space];
            B[space] = N1[space];
        }
        
        if ( leftChar == 'T' ){
            flagTranspose = 0;
            flagTranspose2 = 0;
            tEqua(f1, equals, espin, left, lspin );
            return 0;
        }
        else  if ( leftChar == 'A' ){
            flagTranspose = 1;// a | b c
            
            flagTranspose2 = 0;
            for ( space = 0; space < SPACE ; space++){
                A[space] = N1[space];
                B[space] = N1[space]*N1[space];
            }
        }else  if ( leftChar == 'B' ){
            flagTranspose = 1; // a b | c
            flagTranspose2 = 0;
            for ( space = 0; space < SPACE ; space++){
                A[space] = N1[space]*N1[space];
                B[space] = N1[space];
            }
        } else if ( leftChar == 'C' ){
            flagTranspose = 0;
            flagTranspose2 = 1;
        } else  if ( leftChar == 'D' ){
            flagTranspose = 1; // a b | c
            flagTranspose2 = 1;
            for ( space = 0; space < SPACE ; space++){
                A[space] = N1[space]*N1[space];
                B[space] = N1[space];
            }
        }else  if ( leftChar == 'E' ){
            flagTranspose = 1; // a b | c
            flagTranspose2 = 0;
            flagTranspose3 = 1;
            for ( space = 0; space < SPACE ; space++){
                A[space] = N1[space]*N1[space];
                B[space] = N1[space];
            }
        }
        
        else {
            printf("unknown flag %c\n",leftChar);
            exit(0);
        }
        if ( ( species(f1, tensorBuffers ) >= vector && bodies(f1,tensorBuffers)  == three && part(f1, tensorBuffers)>0) ){
            for ( space =0 ; space < spaces(f1,tensorBuffers) ; space++)
                array[space] = streams( f1, tensorBuffers, rank,space );
        }else {
            printf("oopss\n");
            exit(0);
        }
        
        
        double * pleft[SPACE];
        {
            INT_TYPE bs,o,l;
            for ( l = 0 ; l < CanonicalRank(f1, left, lspin ) ; l++){
                for ( space = 0; space < SPACE ;space++){
                    //
                    pleft[space] = streams( f1, left, lspin,space )+l*LN2[space];
                    
                    bs = 0;
                    
                    if (flagTranspose2){
                        for ( o = 0 ; o < N1[space] ;o++){
                            transpose(N1[space], N1[space],pleft[space]+N1[space]*N1[space]*o,array[bs]+N1[space]*N1[space]*o);
                            
                        }
                        pleft[space] = array[bs++];
                    }
                    if ( flagTranspose  ){
                        transpose(A[space], B[space],pleft[space] ,array[bs]);
                        
                        pleft[space] = array[bs++];
                        
                    }
                    
                    if ( flagTranspose3){
                        for ( o = 0 ; o < N1[space] ;o++){
                            transpose(N1[space], N1[space],pleft[space]+N1[space]*N1[space]*o,array[bs]+N1[space]*N1[space]*o);
                        }
                        pleft[space] = array[bs++];
                        
                    }
                    cblas_dcopy(LN2[space], pleft[space], 1, streams(f1,equals, espin, space)+l*LN2[space] ,1);
                    
                }
            }
            
            
            
        }
        
        
    }
    else if (bodies ( f1, v ) == four ){
        {
            {
                double * pleft[SPACE];
                
                INT_TYPE r;
                N1 = f1->sinc.Basis;
                flagTranspose = 0;
                flagTranspose2 = 0;
                flagTranspose3 = 0;
                flagTranspose4 = 0;
                
                for ( space = 0; space < SPACE ; space++){
                    A[space] = N1[space]*N1[space];
                    B[space] = N1[space]*N1[space];
                }
                for ( space = 0; space < SPACE ; space++){
                    AA[space] = N1[space];
                    BB[space] = N1[space]*N1[space]*N1[space];
                }
                
                
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
                
                
                
                if ( leftChar == 'T' ){//(i,j,k,l)
                    
                } else if ( leftChar == 'a' ){//(i,j,l,k)
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space];
                        B[space] = N1[space]*N1[space];
                    }
                    
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space]*N1[space];
                        BB[space] = N1[space]*N1[space];
                    }
                    
                }else if ( leftChar == 'b' ){//(i,j,k,l)->(j,i,k,l)->(l,j,i,k)->(j,l,i,k)->(i,k,j,l)->(i,k,j,l)
                    flagTranspose = 1;
                    flagTranspose2 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space]*N1[space];
                        B[space] = N1[space];
                    }
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space]*N1[space];
                        BB[space] = N1[space]*N1[space];
                    }
                    
                    
                }else  if ( leftChar == 'c' ){//(i,k,l,j)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space];
                        B[space] = N1[space]*N1[space]*N1[space];
                    }
                }else  if ( leftChar == 'd' ){//(i,l,j,k)
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space]*N1[space];
                        B[space] = N1[space];
                    }
                }else  if ( leftChar == 'e' ){//(i,l,k,j) 6
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space];
                        B[space] = N1[space]*N1[space];
                    }
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space]*N1[space]*N1[space];
                        BB[space] = N1[space];
                    }
                    
                }else  if ( leftChar == 'f' ){//(j,i,k,l)
                    flagTranspose = 1; // a b | c
                }else  if ( leftChar == 'g' ){//(j,i,l,k)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space];
                        B[space] = N1[space]*N1[space];
                    }
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space]*N1[space];
                        BB[space] = N1[space]*N1[space];
                    }
                    
                }else  if ( leftChar == 'h' ){//(j,k,i,l)
                    flagTranspose = 0; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space]*N1[space];
                        B[space] = N1[space];
                    }
                    
                    flagTranspose3 = 1; // a b | c
                    flagTranspose4 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space]*N1[space];
                        BB[space] = N1[space]*N1[space];
                    }
                    
                }else  if ( leftChar == 'i' ){//(j,k,l,i)
                    flagTranspose = 0; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space];
                        B[space] = N1[space]*N1[space]*N1[space];
                    }
                }
                else  if ( leftChar == 'j' ){//(j,l,i,k)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space]*N1[space];
                        B[space] = N1[space];
                    }
                }else  if ( leftChar == 'k' ){//(j,l,k,i)12
                    flagTranspose2 = 2;
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space];
                        B[space] = N1[space]*N1[space];
                    }
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space]*N1[space]*N1[space];
                        BB[space] = N1[space];
                    }
                    
                }else  if ( leftChar == 'l' ){//(k,i,j,l)
                    flagTranspose = 0; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space];
                        B[space] = N1[space]*N1[space];
                    }
                    
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space];
                        BB[space] = N1[space]*N1[space]*N1[space];
                    }
                    
                }else  if ( leftChar == 'm' ){//(k,i,l,j)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space];
                        B[space] = N1[space]*N1[space]*N1[space];
                    }
                }else  if ( leftChar == 'n' ){//(k,j,i,l)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space];
                        B[space] = N1[space]*N1[space];
                    }
                    
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space];
                        BB[space] = N1[space]*N1[space]*N1[space];
                    }
                    
                    
                }else  if ( leftChar == 'o' ){//(k,j,l,i)
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space];
                        B[space] = N1[space]*N1[space]*N1[space];
                    }
                }else  if ( leftChar == 'p' ){//(k,l,i,j)
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space];
                        B[space] = N1[space]*N1[space];
                    }
                }else  if ( leftChar == 'q' ){//i,j,k,l->j,i,k,l->i,k,l,j->k,l,j,i->(k,l,j,i)//18
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space];
                        B[space] = N1[space]*N1[space]*N1[space];
                    }
                    flagTranspose4 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space];
                        BB[space] = N1[space]*N1[space]*N1[space];
                    }
                    
                }
                else  if ( leftChar == 'r' ){//(l,i,j,k)
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space]*N1[space];
                        B[space] = N1[space];
                    }
                }
                else  if ( leftChar == 's' ){//(l,i,k,j)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space]*N1[space];
                        B[space] = N1[space];
                    }
                    flagTranspose3 = 1; // a b | c
                    flagTranspose4 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space];
                        BB[space] = N1[space]*N1[space]*N1[space];
                    }
                    
                }
                else  if ( leftChar == 't' ){//(l,j,i,k)
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space]*N1[space];
                        B[space] = N1[space];
                    }
                }
                else  if ( leftChar == 'u' ){//(l,j,k,i)
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space]*N1[space];
                        B[space] = N1[space];
                    }
                    
                    flagTranspose3 = 1;
                    flagTranspose4 = 1;
                    for ( space = 0; space < SPACE ; space++){
                        AA[space] = N1[space];
                        BB[space] = N1[space]*N1[space]*N1[space];
                    }
                }
                // 0,2,1,0 : l,k,i,j   :: 'v'6
                else  if ( leftChar == 'v' ){//i,j,k,l->j,i,k,l->k,l,j,i->(l,k,i,j)
                    flagTranspose2 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space];
                        B[space] = N1[space]*N1[space];
                    }
                    flagTranspose3 = 1;
                    
                }
                else  if ( leftChar == 'w' ){//(l,k,j,i)//END
                    flagTranspose = 1; // a b | c
                    flagTranspose2 = 1;
                    flagTranspose3 = 1;
                    
                    for ( space = 0; space < SPACE ; space++){
                        A[space] = N1[space]*N1[space];
                        B[space] = N1[space]*N1[space];
                    }
                }
                
                
                
                
                else {
                    printf("unknown flag %c\n",leftChar);
                    exit(0);
                }
                
                if ( ( species(f1, tensorBuffers ) >= vector && bodies(f1,tensorBuffers)  >= four && part(f1, tensorBuffers)>=1) ){
                    for ( space =0 ; space < SPACE ; space++)
                        array[space] = streams( f1, tensorBuffers, rank,space);
                    array[3] = streams( f1, tensorBuffers, rank,0)+LN2[0];
                    
                }else {
                    printf("oopss\n");
                    exit(0);
                }
                
                
                {
                    INT_TYPE bs,o;
                    for ( l = 0 ; l < CanonicalRank(f1, left, lspin ) ; l++){
                        
                        for ( space = 0; space < spaces(f1,left) ;space++){
                            bs = 0;
                            //
                            pleft[space] = streams( f1, left, lspin,space )+l*LN2[space];
                            
                            
                            if (flagTranspose){
                                for ( o = 0 ; o < N1[space]*N1[space] ;o++)
                                    transpose(N1[space], N1[space],pleft[space]+N1[space]*N1[space]*o,array[bs]+N1[space]*N1[space]*o);
                                pleft[space] = array[bs++];
                            }
                            
                            if ( flagTranspose2  ){
                                transpose(A[space], B[space],pleft[space],array[bs]);
                                pleft[space] = array[bs++];
                            }
                            
                            if ( flagTranspose3){
                                for ( o = 0 ; o < N1[space]*N1[space] ;o++)
                                    transpose(N1[space], N1[space],pleft[space]+N1[space]*N1[space]*o ,array[bs]+N1[space]*N1[space]*o);
                                pleft[space] = array[bs++];
                                
                                
                            }
                            
                            if ( flagTranspose4  ){
                                transpose(AA[space], BB[space],pleft[space],array[bs]);
                                pleft[space] = array[bs++];
                                
                            }
                            cblas_dcopy(LN2[space], pleft[space], 1, streams(f1,equals, espin, space)+l*LN2[space] ,1);
                            
                        }
                        
                        
                        
                    }
                    
                }
                
                
                
                
            }
        }
    }
    
    f1->sinc.tulip[equals].Current[espin] = CanonicalRank(f1, left, lspin);
    return 0;
    
}




INT_TYPE tAllCompPermMultiplyMP( INT_TYPE rank, struct field * f1 , enum division left ,INT_TYPE lspin, enum division right ,INT_TYPE rspin, double * sequ){
    if ( bodies(f1, left ) != bodies(f1,right)){
        printf("get real!\n");
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
        train[0] = 'T';
        train[1] = 'N';//(12)
        nPerm = 2;
    }
    else if ( bodies(f1, left ) == three ){
        //
        train[0] = 'T';//(1)            123
        train[1] = 'A';//(123)          231
        train[2] = 'B';//(123).(123)    312
        train[3] = 'C';//(12)           213
        train[4] = 'D';//(13)           321
        train[5] = 'E';//(23)           132
        nPerm = 6;
    }
    else if ( bodies(f1, left ) == four ){
        //
        train[0] = 'T';
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

    for ( i = 0; i < nPerm ; i++)
        sequ[i] = tMultiplyMP(rank, &info,f1,1. , -1, nullVector, 0, train[i], left, lspin, 'N', right, rspin);
        
    
    return nPerm;
}


//INT_TYPE tAddUpComponents( INT_TYPE rank, struct field * f1 , enum division left , enum division right ,  double *up){
//
////    if ( bodies(f1, left ) > two )
////        printf ( "auktung!\n");
//
//    if ( bodies(f1,left ) == one ){
//        return 1;
//    }
//    double rMA[24*24], iMA[24*24], buff[24];
//    INT_TYPE sumj,i,j,i2,i1,nPerm=0,nGroup=0,Pe[24][24];
//    char c,c0;
//
//    for ( i = 0; i < 24*24 ; i++)
//    {
//        rMA[i] = 0.;
//        iMA[i] = 0.;
//    }
//
//    for ( i = 0 ; i < 24 ; i++)
//        for ( j = 0 ; j < 24 ; j++)
//            Pe[i][j] = 0;
//
//
//    if ( bodies(f1, left ) == two ){
//        nPerm = 2;
//        nGroup = 2;
//        c0 = 'N';
//
//
//
//        Pe[0][0] = 1;
//
//        Pe[1][0] = 2;
//
//
//
//    }
//    else if ( bodies ( f1, left ) == three ){
//        nPerm = 6;
//        nGroup = 3;//JLJ:HERE
//        c0 = 'A';
//        Pe[0][0] = 1;
//
//        Pe[1][0] = 2;
//
//        Pe[2][0] = 3;
//        Pe[2][1] = 4;
//        Pe[2][2] = 5;
//        Pe[2][3] = 6;
//
//
////        for ( i = 0 ; i < 6 ; i++)
////            Pe[i][0] = i+1;
//
//
//    }
//    else if ( bodies (f1, left ) == four  ){
//        nPerm = 24;
//        nGroup = 24;//JLJ:HERE
//        c0 = 'a';
////        for ( i = 0 ; i < 24 ; i++)
////            //for ( j = 0; j < 24 ; j++)
////                Pe[i][0] = i+1;
//
//
//        Pe[0][0] = 1;
//
//        Pe[1][0] = 2;
//
//        Pe[2][0] = 3;
//        Pe[2][1] = 4;
//        Pe[2][2] = 5;
//        Pe[2][3] = 6;
//
//        Pe[3][0] = 7;
//        Pe[3][1] = 8;
//        Pe[3][2] = 9;
//        Pe[3][3] = 10;
//        Pe[3][4] = 11;
//        Pe[3][5] = 12;
//        Pe[3][6] = 13;
//        Pe[3][7] = 14;
//        Pe[3][8] = 15;
//
//        Pe[4][0] = 16;
//        Pe[4][1] = 17;
//        Pe[4][2] = 18;
//        Pe[4][3] = 19;
//        Pe[4][4] = 20;
//        Pe[4][5] = 21;
//        Pe[4][6] = 22;
//        Pe[4][7] = 23;
//        Pe[4][8] = 24;
//
//    }else {
//        printf("opps\n");
//        exit(0);
//    }
//    double cmpl,cmpl2;
//    enum division alloy = right;
//    enum division alloyBak;
//    if ( species (f1, alloy ) == vector && bodies(f1,alloy) == one){
//        alloyBak = trainVector;
//    }
//    else     if ( species (f1, alloy ) == vector && bodies(f1,alloy) == two){
//        alloyBak = trainVector2;
//    }
//    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == one){
//
//        f1->sinc.tulip[trainMatrix].header = header(f1,alloy);
//
//        alloyBak = trainMatrix;
//    }else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == two){
//
//        f1->sinc.tulip[trainMatrix2].header = header(f1,alloy);
//
//        alloyBak = trainMatrix2;
//    }
//    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == three){
//
//        f1->sinc.tulip[trainVector3].header = header(f1,alloy);
//
//        alloyBak = trainVector3;
//    }
//    else if ( species(f1, alloy ) == vector && bodies(f1,alloy) == four){
//
//        f1->sinc.tulip[trainVector4].header = header(f1,alloy);
//
//        alloyBak = trainVector4;
//    }
//    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == three){
//
//        f1->sinc.tulip[trainMatrix3].header = header(f1,alloy);
//
//        alloyBak = trainMatrix3;
//    }
//    else if ( species(f1, alloy ) == matrix && bodies(f1,alloy) == four){
//
//        f1->sinc.tulip[trainMatrix4].header = header(f1,alloy);
//
//        alloyBak = trainMatrix4;
//    }
//    else if ( species(f1, alloy ) == quartic ){
//        f1->sinc.tulip[trainQuartic].header = header(f1,alloy);
//
//        alloyBak = trainQuartic;
//    } else{
//        printf("cd: tObject species %d \n", alloy);
//
//        exit(0);
//
//    }
//    for ( cmpl = 0; cmpl < 2 ; cmpl++){
//        cmpl2 = cmpl;
//       // for ( cmpl2 = 0; cmpl2 < 2; cmpl2++){
//
//        for (c = 0; c < nPerm ; c++){
//            if ( c ){
//                tPermute(rank, f1,(c-1)+c0, right, cmpl, alloyBak, rank);
//                tAllCompPermMultiplyMP(rank, f1, left, cmpl2, alloyBak,rank, buff);
//            }
//
//            else {
//                tAllCompPermMultiplyMP(rank, f1, left, cmpl2, right,cmpl, buff);
//            }
//            if ( cmpl == 0 && cmpl2 == 0 )
//                cblas_daxpy(nPerm, 1.0, buff, 1, rMA+c*nPerm, 1);
////            if ( cmpl == 1 && cmpl2 == 0 )
////                cblas_daxpy(nPerm, 1.0, buff, 1, iMA+c*nPerm, 1);
////            if ( cmpl == 0 && cmpl2 == 1 )
////                cblas_daxpy(nPerm, -1.0, buff, 1, iMA+c*nPerm, 1);
//            if ( cmpl == 1 && cmpl2 == 1 )
//                cblas_daxpy(nPerm, 1.0, buff, 1, rMA+c*nPerm, 1);
//        }
//    }
////    printf("mag %f\n", magnitude(f1, left));
////    for ( i = 0; i < 2; i++){
////        for ( j = 0 ; j < 2 ; j++)
////            printf("(%1.2f)", rMA[i*nPerm+j]);
////        printf("\n");
////    }
//    for (  i = 0 ; i < nGroup ; i++){
//        sumj = 0;
//        for ( j = 0; j < nPerm ; j++)
//            if ( Pe[i][j] ){
//                sumj++;
//                for ( i1 = 0 ; i1 < nPerm ; i1++)
//                    for ( i2 = 0; i2 < nPerm ; i2++)
//                        up[i] += get(bodies(f1, left ), Pe[i][j], i1)* rMA[i1*nPerm+i2]* get(bodies(f1, left ), Pe[i][j], i2);
//        }
//        up[i] /= nPerm;
//    }
//
//    return nGroup;
//}


INT_TYPE tTabulateProjection( INT_TYPE rank, struct field * f1 , enum division left , enum division right ,  double *up){
    
    if ( bodies(f1,left ) == one ){
        return 1;
    }
    INT_TYPE i,g,p,nPerm=0,nGroup=0;
    
    
    if ( bodies(f1, left ) == two ){
        nPerm = 2;
        nGroup = 2;
    }
    else if ( bodies ( f1, left ) == three ){
        nPerm = 6;
        nGroup = 3;
        
    }
    else if ( bodies (f1, left ) == four  ){
        nPerm = 24;
        nGroup = 5;
    }else {
        printf("opps\n");
        exit(0);
    }

    double buff[24];
    DCOMPLEX gup[24];
    for ( i = 0; i< 24 ;i++)
        gup[i] = 0.;

    tAllCompPermMultiplyMP(rank, f1, left, 0, right,0, buff);
    for ( g = 0; g < nGroup ; g++)
        for ( p = 0; p < nPerm ; p++)
            gup[g] += get(bodies(f1,right), g+nPerm+1, p)*buff[p];

    tAllCompPermMultiplyMP(rank, f1, left, 1, right,1, buff);
    for ( g = 0; g < nGroup ; g++)
        for ( p = 0; p < nPerm ; p++)
            gup[g] += get(bodies(f1,right), g+nPerm+1, p)*buff[p];

    tAllCompPermMultiplyMP(rank, f1, left, 0, right,1, buff);
    for ( g = 0; g < nGroup ; g++)
        for ( p = 0; p < nPerm ; p++)
            gup[g] += I*get(bodies(f1,right), g+nPerm+1, p)*buff[p];

    tAllCompPermMultiplyMP(rank, f1, left, 1, right,0, buff);
    for ( g = 0; g < nGroup ; g++)
        for ( p = 0; p < nPerm ; p++)
            gup[g] += -I*get(bodies(f1,right), g+nPerm+1, p)*buff[p];
    
    for ( g = 0; g < 24 ; g++)
        up[g] = cabs(gup[g])/nPerm;
        
    return nGroup;
}


INT_TYPE tVeto ( enum body bd , INT_TYPE type , INT_TYPE *a){
    if ( bd == two ){
        if ( type == 2 )
            if ( a[0] != a[1] )
                return 1;
    } else if ( bd == three ){
        if ( type == 1 )
            return 1;
        if (type == 2 ){
            if ( a[0] != a[1] && a[1] != a[2] && a[2] != a[0])
                return 1;
        }
        if ( type ==3  ){
            if ( a[0] != a[1] || a[1] != a[2] || a[2] != a[0])
                return 1;
        }
    }else if ( bd == four ){
        if ( type == 1 )
            return 1;
        if ( type == 2)
            if ( (a[0] != a[1] + a[1] != a[2] + a[2] != a[0] + a[0] != a[3] + a[1] != a[3] + a[2] != a[3] ) == 0)
                return 1;
        if ( type == 3 )
            if ( ( (a[0] != a[1]) + (a[1] != a[2]) + (a[2] != a[0]) + (a[0] != a[3]) + (a[1] != a[3]) + (a[2] != a[3])) < 3)
                return 1;
        if ( type == 4 )
            if ( ( (a[0] != a[1]) + (a[1] != a[2]) + (a[2] != a[0]) + (a[0] != a[3]) + (a[1] != a[3]) + (a[2] != a[3])) < 2)
                return 1;
        if ( type == 5 )
            if ( ( (a[0] != a[1]) + (a[1] != a[2]) + (a[2] != a[0]) + (a[0] != a[3]) + (a[1] != a[3]) + (a[2] != a[3])) < 6)
                return 1;
    }
    
    return 0;
}
