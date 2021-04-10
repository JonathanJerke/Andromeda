/**
 *  system.h
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

/**
 *Andromeda: a few-body plane wave calculator
 *
 *v9.8
 *quantumGalaxies.org
 *
 *Jonathan Jerke
 *Bill Poirier
 *
 *Texas Tech University
*/

#ifndef system_h
#define system_h
#define MOMENTUM_INTERVAL 51
#define MOMENTUM_INTERVAL_RE 51
#define MOMENTUM_INTERVAL_IM 51
//#define SEPARATE_ONE_BODY
//#define VERBOSE_ALS
#define writeHDF5
#define readHDF5

///These additions do work now,  be careful about legacy
#define READ_FAST
#define WRITE_FAST
///to preserve info, may need to incrementally bootup.

#define MODULARIZE_INPUT
#define MODULARIZE_OUTPUT


#define RAND_FOUNDATION
//#define CHERRY_PICKER
///to compile with acceleration in APPLE, not for distribution
///Number of components, no limit
#define SPACE 3
///Number of 'atoms' under geometry,  really no reason to to keep this big
#define MAXATOM 3
///All internal numbers are in au, including Hartrees.
///input in Angstroms unless 'Angstroms 0' set
#define a0  0.52917721
///natural
#define pi  3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
///Standard Position of identity group action
#define CDT 1
///switch
#define VERBOSE 0
///print info on ALS decomposition
//#define VERBOSE_ALS
///Maximum number of input terms at prompt, at minor cost to increase
#define MAXTERM 100

///Maximum bodies per component currently supported
///only coded for 4bodies with PBC
#define MAXBODY 4

///basically a fundamental on ## bodies in interaction
#define MAX_SPLIT 2

///Switch to make stuff complex at compiler level,not working yet
//#define COMPLEXME
///really just vistidual, but some places still use it
#define COMPONENT 3
///basically overwritten, will remove
#define MAX_PARAM_FUNC 4
///The block-memory commands act to negate memory allocations,  this is the maximum number of blocks.  Leave this alone, unless you add blocks.
#define BLOCK_COUNT 24

#ifdef APPLE
    #define SPHERE
    #define BIT_INT
    #define MAXSTRING 76
    #define SUPERMAXSTRING 200
    #define MAX_FILE 1000
    #define MAX_CORE 12
#else
///for including omp.h
    #define OMP
///intel MKL

#ifndef MKL
/// gnu compatible
    #define ATLAS
    #define LAPACKE
    #define BIT_INT
#endif
///Probably too much, but dont care,
    #define MAXSTRING 76
///Probably too much, but dont care,
    #define SUPERMAXSTRING 200
///Maximum number of .mac files that are loadable,  (not vectors, rather files of .vector)
    #define MAX_FILE 1000
///Maybe necessary for integrating momentums in GaussianSinc-basis--otherwise, not used.
//#define GSL_LIB
///Set true if you want to use GSL cblas
//#define GSL_CBLAS
///Normal super computer size, could be more or less,
    #define MAX_CORE 12
#endif

#endif /* system_h */
