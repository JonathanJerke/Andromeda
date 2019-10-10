/*
 *  system.h
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


//HERE LIES THE external CODE level changes!

#ifndef system_h
#define system_h

#define COMPONENT 3
#define PARTICLE 1
#define SPACE (COMPONENT * PARTICLE)
#define MAXATOM 4


//Currently under a dozen actual blocks possible, there is extra.
#define BlockCount 16

//#define APPLE

#ifdef APPLE
#define BIT_INT
#define MAXSTRING 1
#define SUPERMAXSTRING 10
#define MAXFILE 1
//APPLE
#define MaxCore 2
//APPLE
#endif





#ifndef APPLE
//NOT APPLE!!
#define MAXSTRING 1024
#define MAXFILE 100
#define OMP
//SELECT BIT-length
//int:
#define BIT_INT
//long int:         #define BIT_LONG
//long long int:    #define MKL

//USE libgslcblas else use libblas
#define GSL_CBLAS



//NEEDED
#define GSL_LIB

//UNIX
#define MaxCore 4
//UNIX
#ifndef GSL_LIB
    THIS PROGRAM NEEDS INTEGRATION ROUTINES
#endif


//NOT APPLE!!
#endif



#endif /* system_h */
