/*
 *  system.h
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


//HERE LIES THE external CODE level changes!


#ifndef system_h
#define system_h

//few fixed memory arrays
#define MAXATOM 4

//#define OMP
//#define MKL
#define GSL_LIB


#define MaxCore 2





//need NOT define.
//#define PARAMETER_PATH ""
#ifndef GSL_LIB
    THIS PROGRAM NEEDS INTEGRATION ROUTINES
#endif


#endif /* system_h */
