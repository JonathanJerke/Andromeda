/**
 *  constants.h
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
 *v9.6
 *quantumGalaxies.org
 *
 *Jonathan Jerke
 *Bill Poirier
 *
 *Texas Tech University
 *Welch Foundation
 *Army Research Office
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "system.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef APPLE
#include "Accelerate/Accelerate.h"
#include "complex.h"
typedef double __complex__ DCOMPLEX;
typedef __CLPK_integer  inta;
typedef __CLPK_doublereal floata;
typedef unsigned long ADDRESS_TYPE;
#else

typedef double floata;

#ifdef  OMP
#include "omp.h"
#endif

#ifdef GSL_LIB
///May not need all of these,
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#endif

#ifdef MKL
#include "mkl.h"
#include "mkl_lapacke.h"
#include "mkl_types.h"
typedef MKL_INT inta;
typedef MKL_INT ADDRESS_TYPE;
typedef double __complex__ DCOMPLEX;
typedef MKL_Complex16 DCOMPLEX_PRIME;
#endif


#ifdef ATLAS
#include "cblas-atlas.h"
#define LAPACK_COL_MAJOR 102
#endif


#ifdef LAPACKE
#include "lapacke.h"
//#include "lapacke_utils.h"
#endif

#ifdef GSL_CBLAS
#include "gsl/gsl_blas.h"
#endif

#ifdef BIT_INT
typedef long long int ADDRESS_TYPE;
typedef int inta;
typedef double __complex__ DCOMPLEX;
typedef double __complex__ DCOMPLEX_PRIME;
#endif

#endif

#include "complex.h"

#ifdef BIT_LONG
typedef long long int ADDRESS_TYPE;
#define lapack_int long int
typedef lapack_int inta;
typedef double __complex__ DCOMPLEX;
typedef double __complex__ DCOMPLEX_PRIME;
#endif


#ifdef COMPLEXME
typedef DCOMPLEX mea;
#else
typedef double mea;
#endif

/**
 *Component by component body type
 */
enum bodyType {
    nada,
    one,
    two,
    three,
    four,
    five,
    six
};

/**
 *not used.
 */
enum particleType{
    all
};






/**
 *Command run by job
 */
enum phaseType{
    buildFoundation,
    productKrylov,
    solveRitz,
    buildTraces,
    formOCSB,
    iterateOCSB
};

/**
 *For future
 */
enum calculationType{
    nullCalculation,
    electronicStuctureCalculation
};

/**
 *TO turn off printing for bootup.
 */
enum bootType{
    noMatrices,
    fullMatrices
};

/**
 *Interaction type labed, i.e. tv1 = act on particle 1 of component in question.
 */
enum blockType{
    id0,
    tv1,
    tv2,
    tv3,
    tv4,
    tv5,
    tv6,
    e12,
    e13,
    e23,
    e14,
    e24,
    e34,
    e15,
    e25,
    e35,
    e45,
    e16,
    e26,
    e36,
    e46,
    e56
};

/**
 *Controling aspects of Hamiltonian Terms and more generally
 *the concept of a division.
 */
enum genusType{
    scalar,
    vector,
    matrix,
    outerVector,
    diagonalMatrix,
    eikon,
    eikonDiagonal,
    eikonOffDiagonal,
    eikonSemiDiagonal,
    eikonDeriv,
    eikonKinetic,
    eikonConstant,
    eikonLinear,
    eikonSpring,
    eikonElement,
    eikonOuter
};

/**
 *The function type, defined in InverseLaplaceTransform
 */
enum functionType{
    nullFunction,
    Pseudo,
    Yukawa,
    Coulomb,
    Morse,
    LennardJones,
    LDA,
    BLYP,
    Gaussian
};


/**
 *For counting number of spins associated with a division
 */
enum spinType {
    none,
    real,
    cmpl,
    parallel
};


/**
 *Exclusively canon[space<SPACE] for objectAllocation
 *and canon[SPACE] for bufferAllocation
 *just one of those redundancies
 */
enum memoryType {
    noAllocation,
    objectAllocation,
    bufferAllocation
};

/**
 *Printed at runtime, 'y'
 */

enum basisElementType {
    nullBasisElement,
    SincBasisElement,
    GaussianBasisElement,
    DiracDeltaElement,
    StateBasisElement,
    overlapBasisElement
};

/**
 *canon parameter, for labeling component
 */
enum componentType {
    nullComponent,
    spatialComponent1,
    spatialComponent2,
    spatialComponent3,
    periodicComponent1,
    periodicComponent2,
    periodicComponent3,
};


/**
 *For periodic boundary conditions
 */
enum noteType {
    nullNote,
    interactionCell
};

/**
 *The Lebequse metric of integration of beta integral
 *
 *These controls can create most any kind of quantum operator.
 */
enum metricType {
    dirac,
    separateDirac,
    interval,
    semiIndefinite,
    pureInterval,
    pureSemiIndefinite
};


















/**
 *Enumerate memory controls
 *
 *This is control by negation, setting a positive switch 'blockMemory 1' will turn off allocaiton
 *of totalVector, which is the primary buffer.
 *
 *TrainVectors will allow ALS to run
 *CopyVectors are only used in exploratory works.
 *Transfer Basis will allow one to boot from a previous (smaller) basis.
 */
enum blockMemoryType{
    passBlock,
    blockTotalVectorBlock,
    blockTrainVectorsblock,
    blockCopyBlock,
    blockTransferBasisblock,
    blockMatrixElementsblock,
    blockPermutationsblock,
    blockParallelMultiplyblock,
    blockParallelMatrixElementblock,
    blockParallelPermuteblock,
    blockPrintStuffblock,
    blockTotalVectorParallelBlock,
    blockComponentblock,
    blockDiagonalMatrixblock
};


/**
 *Direct
*/

///B
///block: 3,4, 5, 6,7, 9

///b
///block: 3, 4, 5, 6, 7,  9

///C
///block 1,2,3,4,6,7,9

///c
///block 1,2,3,4,6,7,8,9

///D
///block 3,4,5,6,7,9


/**
 *Collect
*/

///B
///block: 3, 4,5, 7,  9

///b
///block: 3, 4, 5, 7, 9

///C
///collect 1
///block 3,4, 6,7,9

///c
///collect 1
///block 3,4,6,7,8,9

///D
///block 3,4,5,6,7,9

/**
 *CollectFilter at output
*/

///B
///block: 3, 4,5, 7, 9

///b
///filter 4
///block: 3, 4, 5, 7, 9

///C
///collect 1
///block 3,4, 6,7,9

///c
///collect 1
///block 3,4,6,7,8,9

///D
///filter 4
///block 3,4,5,6,7,9


/**
 *Filter post creation of vector
*/

///B
///block: 3,4, 5, 6,7,  9

///b
///block: 3, 4, 5, 7,  9
///filter 4

///C
///block 1,2,3,4,6,7,9

///c
///block 1,2,3,4,6,7,9

///D
///block 3,4,5,7,9
///filter 4


/**
 *InputFilter
 *Collect and filter at input
*/

///B
///block: 3,4,5,7,9

///b
///block: 3,4,5,7,9

///C
///filter 3
///collect 1
///block 3,4,7,9

///c
///filter 3
///collect 1
///block 3,4,7,9

///D
///filter 3
///block 3,4,5,7,9

















/**
 *Divisions define and indicate data.
 *Lined up with 500 so a division number can be easily read.
 */
enum division{
    nullName,
    nullOverlap,
    nullVector,
    nullMatrix,
    buffer,
    bra,
    bra2,
    bra3,
    bra4,
    bra5,
    bra6,
    bra7,
    lowdinVec,
    lowdinMatrix,
    overlap,
    Ha, 
    Iterator,
    printOperator,
    copy,
    copyTwo , 
    copyThree,
    copyFour,
    vectorDiagonalMatrix,
    eigen,
    PauliX,
    PauliY,
    PauliZ,
    trainQuartic,
    cycleVector,
    cycleMatrix,
    cycleQuartic,
    directBak, 
    oneBody,
    entropyVector,
    entropyUnit,
    vectors3,
	basis,
    subBasis,
    quad,
    quad2,
    bandBasis,
    quadCube,
    quad2Cube,
    diagonalCube,
    copyVector,
    copyTwoVector,
    copyThreeVector,
    copyFourVector,
    scalarTemp,
    productVector,
    permutationVector,
    permutation2Vector,
    permutation3Vector,
    permutation4Vector,
    eikonBuffer,
    basisBuffers,
    canonicalBuffers,
    guideBuffer,
    trackBuffer,
    canonicalVector,
    vectorCubeBuffers,
    diagonalQuad,
    diDiagonalQuad,
    diagonal3VectorA,
    diagonal2VectorA,
    diagonal2VectorB,
    diagonal1VectorA,
    diagonalVectorA,
    diagonalVectorB,
    diagonal1VectorB,
    diagonal1VectorC,
    diagonal1VectorD,
    canonicalBuffersB,
    canonicalBuffersC,
    canonicalmvVector,
    canonicalmv2Vector,
    canonicalmv3Vector,
    canonicaldotVector,
    canonicaldot2Vector,
    canonicaldot3Vector,
    canonicalvvVector,
    canonicalvv2Vector,
    canonicalvv3Vector,
    canonicalmeVector,
    canonicalme2Vector,
    canonicalme3Vector,
    multiplyVector,
    component,
    componentTotal,
    CanonicalBuffers,
    oneVector,
    twoVector,
    twoBodyRitz,
    totalVector,
    squareTwoVector,
    dsyBuffers,
    matrixSbuild,
    matrixHbuild,
    vectorHbuild,
    twoBodyProjector,
    eigenVectors
};

typedef enum bodyType bodyType;
typedef enum particleType particleType;
typedef enum phaseType phaseType;
typedef enum calculationType calculationType;
typedef enum bootType bootType;
typedef enum blockType blockType;
typedef enum genusType genusType;
typedef enum functionType functionType;
typedef enum spinType spinType;
typedef enum memoryType memoryType;
typedef enum basisElementType basisElementType;
typedef enum componentType componentType;
typedef enum noteType noteType;
typedef enum metricType metricType;
typedef enum blockMemoryType blockMemoryType;
typedef enum division division;

typedef struct function_label function_label;
typedef struct atom_label atom_label;
typedef struct basisElement_label basisElement_label;
typedef struct canon_label canon_label;
typedef struct metric_label metric_label;
typedef struct term_label term_label;
typedef struct runTime runTime;
typedef struct space_label space_label;
typedef struct value_label value_label;
typedef struct name_label name_label;
typedef struct nameDispenser_label nameDispenser_label;
typedef struct sinc_label sinc_label;
typedef struct input_label input_label;
typedef struct field field;
typedef struct general_index general_index;
typedef struct general_2index general_2index;
typedef struct input input;
typedef struct calculation calculation;

struct function_label{
    ///number of canonical ranks of function per chain element
    inta interval;
    ///degree of off-diagonal
    inta contr;
    ///function type
    functionType fn;
    ///allowances for function parameters, not used
    double param[MAX_PARAM_FUNC];
};

struct atom_label {
    ///count begins with 1, 3 components allocated
    double position[3+1];
    ///ion charge
    inta Z;
};

struct basisElement_label {
    ///the type of of basis in component
    basisElementType basis;
    ///the type of a particular component
    componentType component;
    ///for periodic systems
    noteType note;
    ///lattice index
    inta index;
    ///boost for periodicBoostComponent
    inta index2;
    double length;
    ///the position on the axis of the left edge of the grid
    double origin;
    ///L , number of grid sites in 1 dimension
    inta grid;
};

struct dimensions_label {
    ///1D lattice spacing
    double lattice;
    ///in [0,1], larger values will stage basis into momentum more
    double attack;
    ///spatial posiiton of the anchor
    double origin;
    ///the fraction of grid for which origin attached; internally anchor is transformed to 0.
    double anchor;
};

struct canon_label {
    ///actual pointer , developers ignore this level
    floata *stream ;
    ///type of basis element
    basisElementType basis;
    ///the size of the component group
    componentType component;
    ///particle count on component, nada for null components throughout
    bodyType body;
    ///an ID allowing for grouping of components
    inta label;
    ///L , the length of the 1D basis, 'count1Basis'
    inta count1Basis;
    ///L+, the increment of L under 'count1Stage'
    inta count1Inc;
    ///the ID of the component in the group
    inta space;
    ///all the double informaton
    struct dimensions_label particle[MAXBODY+1];
};

struct metric_label {
    ///future
    inta pow[SPACE];
    ///future
    inta powB[SPACE];
    ///future
    inta deriv[SPACE];
    ///a potential and inner parameters
    function_label fn;
    ///metric type, i.e. point, interval, or indefinite interval
    metricType metric;
    ///lower and upper bounds,  where -1 => infinity
    double beta[2];//lower and upper bound
};

struct term_label {
    ///optional, loaded for diagonalMatrix
    char filename[MAXSTRING];
    ///internal potential parameters
    function_label func;
    ///description of term
    char desc[16];
    ///external multiplication of term
    double scalar ;
    ///internal indicator of interaction type, i.e. 'kinetic','spring',...
    inta type;
    ///multiply type and particle address, i.e. e2= 2 and 7 = e12
    blockType bl ;
    ///Symmetry Adpated action on vector multiply
    inta act;
    ///rescaling particle 1 of twoBody interaction
    double adjustOne;
    ///ID of component group to address
    inta label;
    ///invert flag on particle 1
    inta invert;
    ///new term flag, may be 0 or 1
    inta headFlag ;
    ///tieing oneBody field to atomic geometry and Z
    inta atom;
    ///for state elements
    inta bra;
    ///for state elements
    inta ket;
    ///multiply each Gaussian Term by more G(0)'s
    inta embed;
    ///the metric of the beta
    struct metric_label mu;
};

struct runTime {
    ///maxiumum swing in canonical rank
    inta dynamic;
    ///scripting level changes to memory allocations, mostly inert
    blockMemoryType memBlock[BLOCK_COUNT];
    ///number of memory allocated processing lanes for OMP
    inta NLanes;
#ifdef OMP
    ///number of OMP cores total
    inta NSlot;
#endif
#ifdef MKL
    ///number of cores per lane
    inta NParallel;
#endif
    ///standard tolerance for canonical decomposition, has units
    double TOLERANCE ;
    ///relative Tolerance for canonical decomposition, does not have units
    double relativeTOLERANCE;
    ///how small vectors can get
    double THRESHOLD;
    ///COLLECT MAX CONDITION number
    double XCONDITION;
    ///Beylkin's condition parameter
    double ALPHA ;
    ///controlled by 'maxCycle'
    inta MAX_CYCLE;
    ///inert
    calculationType calcType;
    ///ritz,krylov,...
    phaseType phaseType;
};

struct space_label{
    ///actual memory address of allocated named element, or buffer,
    ADDRESS_TYPE Address;
    ///when not vector, it will allocate based on this parameter, unless its an eikon
    bodyType body;
    ///the type/address of the interaction iff a matrix
    blockType block;
    ///for SA action
    inta act;
    ///will invert particle 1
    inta invert;
    ///for explicit matrix elements, called state elements
    inta bra;
    ///for explicit matrix elements, called state elements
    inta ket;
};

struct value_label{
    ///vector title
    char title[MAXSTRING];
    ///a way to label the iteration number
    inta stage;
    ///the energy
    double value;
    ///the Symmetry Adaption entropy
    double value2;
    ///the best irrep description
    inta symmetry;
};

struct name_label {
    ///pointer iff name != name
    division name;
    ///vector, matrix,...
    genusType species;
    ///parallel, complex, real
    spinType spinor;
    ///another term
    division linkNext;
    ///lump into a term, defined like link, in so far that it adds
    division chainNext;
    ///equal lineups in chain will multiply together
    inta multId;
    ///within an SOP collective...direct sum
    division loopNext;
    ///Allocated number of elements of type defined wherein
    inta Partition;
    ///extra level of protection of memory types
    memoryType memory;
    ///various numerical values
    value_label value;
    ///allocations, commands, and internal memory pointers
    space_label space[SPACE+1];
    ///dynamic canonical length of SOP-object
    inta Current[MAX_CORE];
    ///if name != name , then Begin is a pointer.
    inta Begin[MAX_CORE];
};

struct nameDispenser_label{
    ///Current label index, will return this one
    inta currLabel;
    ///defined by 'names', will control max number of name_labels allocated
    inta maxLabel;
    ///division of block of names
    division head;
};

struct sinc_label {
    ///many 1-term twoBody SOP vectors
    struct nameDispenser_label eikonLabels;
    ///many 0-term nullBody SOP vectors
    struct nameDispenser_label nullLabels;
    ///control the process of defining terms
    bootType boot;
    ///count of various matrix sizes
    inta maxEV;
    ///Irreducible representation called for, may be zero.
    inta irrep;
    ///flags on booted memory
    inta bootedMemory;
    ///SOP component definitions
    canon_label canon[SPACE+1];
    ///names that pass and control SOP meaning
    name_label *name;
    ///used for storing the foundation
    division user;
    ///for loading in vectors as operators
    division vectorOperator;
    ///last name
    division end;
    ///complexity / parallel , exactly the number of allocated versions in terms of MAX_CORE
    spinType cmpl;
    ///pointer to runTime environment
    struct runTime * rt;
};

struct input_label {
    ///range of canonRank
    inta flex;
    ///0 - no terms, -1 - all terms, n -- nth term
    inta OpIndex;
    ///used in slightly antiquited SA
    inta body;
    ///used in slightly antiquited SA
    inta irrep;
    ///number of files for 'sum' load
    inta matrices ;
    ///list of filenames for 'sum' load
    char matrixList [MAX_FILE][MAXSTRING];
    ///number of files for 'sum' load
    inta files ;
    ///list of filenames for 'vector' load
    char fileList [MAX_FILE][MAXSTRING];
    ///number of files for 'operator' load
    inta filesVectorOperator ;
    ///list of filenaems for 'operator' load
    char fileVectorOperator [MAX_FILE][MAXSTRING];
    ///may be 1 or 2
    inta Iterations;
    ///number of vectors io
    inta nStates;
    ///number of ranks at input
    inta iRank ;
    ///number of ranks for processing
    inta canonRank;
    ///max ranks for (all) krylovs
    inta xRank;
    ///max number of internal vectors, for foundation
    inta qFloor;
    ///flag for operator loaded by 'operator'
    inta nOperator;
    ///+1-> filter on input to program; +2-> filter at input. ; +4 -> filter after processing vectors 
    inta filter;
    ///switch 0/1 control max condition of input of krylov vectors
    inta collect;
    ///secondary level controlling complexity of numbers
    inta cmpl;
};

struct field {
    ///for boot up
    struct input_label i;
    ///for most processing
    struct sinc_label f;
};

struct input {
    ///number of empty names
    inta numNames;
    ///number of twoBody names
    inta numVectors;
    ///switch for loading krylov matrix elements
    inta build;
    ///number of chain elements, not terms exactly
    inta termNumber ;
    ///*Term input container
    struct term_label terms[MAXTERM];
    ///avoid printing to disk lower iterations
    inta minIterationPrint;
    ///dir for 'control'
    char controlPath[MAXSTRING];
    ///shiftVector[0] + H shiftVector[1]
    double shiftVector[2];
    ///switch shift of H
    inta shiftFlag ;
#ifdef OMP
    inta omp;
#endif

#ifdef MKL
    inta mkl;
#endif
    ///wave-let momentum domain ,  counts from 0,1,2,3,...
    inta SymmetrizedGaussianLevel;
    ///length scale of cloud
    floata SymmetrizedGaussianWidth;
    ///for krylov, max rank per term
    inta Lambda;
    ///Gb of of total allocated ram cap
    inta RAMmax ;
    ///0,1 flag for switching geometry to Angstroms
    inta Angstroms;
    ///atom geometry
    struct atom_label atoms[MAXATOM+1];
    ///number of atoms
    inta Na;
    ///ocsb number of partitions
    inta nocsb;
    ///ocsb partition number
    inta iocsb;
};

struct calculation {
    ///the name of the computation leading script '*Body ----'
	char name[MAXSTRING];
    ///general input parameters
	struct input i;
    ///runtime parameters
    struct runTime rt;
};

#endif




