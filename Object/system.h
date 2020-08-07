/**
 *  system.h
 *
 *
 *  Copyright 2020 Jonathan Jerke and Bill Poirier.
 *  We acknowledge the generous support of Texas Tech University,
 *  the Robert A. Welch Foundation, and the Army Research Office.
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

/**
 *Andromeda: a few-body plane wave calculator
 *
 *v9.0
 *quantumGalaxies.org
 *
 *Jonathan Jerke
 *Bill Poirier
 *
 *Texas Tech University
*/

#ifndef system_h
#define system_h
#define ASTER_FLAT
///to compile with acceleration in APPLE, not for distribution
//#define APPLE
///Number of components, no limit
#define SPACE 12
///Number of 'atoms' under geometry,  really no reason to to keep this big
#define MAXATOM 16
///All internal numbers are in au, including Hartrees.
///input in Angstroms unless 'Angstroms 0' set
#define a0  0.52917721
///natural
#define pi  3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
///Standard Position of identity group action
#define CDT 1
///switch
#define VERBOSE 0
///Maximum number of input terms at prompt, at minor cost to increase
#define MAXTERM 1000
///Maximum bodies per component currently supported
#define MAXBODY 3
///Switch to make stuff complex at compiler level,not working yet
//#define COMPLEXME
///really just vistidual, but some places still use it
#define COMPONENT 3
///basically overwritten, will remove
#define MAX_PARAM_FUNC 4
///The block-memory commands act to negate memory allocations,  this is the maximum number of blocks.  Leave this alone, unless you add blocks.
#define BLOCK_COUNT 6

#ifdef APPLE
    #define SPHERE
    #define BIT_INT
    #define MAXSTRING 1
    #define SUPERMAXSTRING 10
    #define MAX_FILE 1
    #define MAX_CORE 2
#else
///for including omp.h
    #define OMP
///intel MKL
    #define MKL

///Probably too much, but dont care,
    #define MAXSTRING 1024
///Probably too much, but dont care,
    #define SUPERMAXSTRING 2048
///Maximum number of .mac files that are loadable,  (not vectors, ---> files of .vector)
    #define MAX_FILE 10000
///Maybe necessary for integrating momentums in GaussianSinc-basis--otherwise, not used.
//#define GSL_LIB
///Set true if you want to use GSL cblas
//#define GSL_CBLAS
///Normal super computer size, could be more or less,
    #define MAX_CORE 72
#endif

#endif /* system_h */
