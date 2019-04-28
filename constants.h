/*
 *  constants.h
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

//VERSION 5.7

#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "system.h"
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#define VERBOSE 0 

#ifdef APPLE
#include "Accelerate/Accelerate.h"
#include "complex.h"
typedef double __complex__ DCOMPLEX;
typedef __CLPK_integer  INT_TYPE;
typedef __CLPK_doublereal Stream_Type;
typedef unsigned long ADDRESS_TYPE;
#else


typedef double Stream_Type;

#include "omp.h"
#include "complex.h"

#ifdef GSL_LIB
#define MAXINT 1000
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
typedef MKL_INT INT_TYPE;
typedef MKL_INT ADDRESS_TYPE;
typedef double __complex__ DCOMPLEX;
typedef MKL_Complex16 DCOMPLEX_PRIME;
#else
typedef long ADDRESS_TYPE;
#include "omp.h"
#include "lapacke.h"
#include "lapacke_utils.h"
#include "gsl/gsl_blas.h"
typedef lapack_int INT_TYPE;
typedef double __complex__ DCOMPLEX;
typedef double __complex__ DCOMPLEX_PRIME;


#endif


#endif

#define MAXSTRING 1024
#define MaxParamFunc 4

#define Mag  0.0000021271911715764754 // Hartree/Tesla
//input in Angstroms...
#define a0  0.52917721  //A
//all internal numbers are in au
#define Ry 13.60569253  //eV
//used in spin-orbit coupling
#define Alpha 0.0072973525698  //unitless

#define pi  3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
        
#define iNS 16


#define MAX_PROTOTYPE_ATOMS 100
#define MAX_ANGULAR_SYMMETRY 2
#define MAX_LABEL 3
#define NO_SINC_LABEL_NAME 2

////////////////////////
/////////
////

#define NO_ACTION 'N'
#define CDT 'T'
#define MaxSpin 2
#define nSAG 5

//pseudo-potentials


enum bodyType {
    nada,
    one,
    two,
    three,
    four
};

enum particleType{
    nullParticle,
    electron,
    proton,
    all,
};

enum phaseType{
    buildFoundation,//0
    productKrylov,//1
    solveRitz,//2
    decomposeTensor,//3
    bandStage,//4
    frameDensity, //5
    svdOperation //6
};

enum calculationType{
    nullCalculation,
    electronicStuctureCalculation,
    clampProtonElectronCalculation,
    protonsElectronsCalculation,
    svdCalculation,
};

enum block{
    id0,
    tv1,//kinetic + potential block 1
    tv2,//1-b block 2
    tv3,
    tv4,
    iv1,
    iv2,
    iv3,
    iv4,
    e12,//ee interaction between label 1 and 2
    e13,
    e23,//changed..
    e14,
    e24,
    e34
};

enum genus{
    scalar,
    vector,
    matrix,
    outerVector,
};

enum shape{
    Cube,
    Rds,
    RdsBasis,
    Band,
    BandBasis,
};

struct atom_label {
    INT_TYPE Z;
    INT_TYPE label;
};

struct atom {
	double position[4];
	struct atom_label label;
};

enum spinType {
    none,
    real,
    cmpl,
    parallel
};

enum memoryType {
    noAllocation,
    objectAllocation,
    bufferAllocation
};

enum functionType{
    nullFunction,//0
    Pseudo,//1
    Yukawa,//2
    Coulomb,//3
    Morse,//4
    LJ,//5
    LDA,//6
    BLYP//7
};

enum basisElementType {
    nullBasisElement,
    SincBasisElement,
    GaussianBasisElement,
    DiracDeltaElement
};

enum componentType {
    nullComponent,
    spatialComponent1,
    spatialComponent2,
    spatialComponent3,
    periodicComponent1,
    periodicComponent2,
    periodicComponent3,
};

struct basisElement {
    enum basisElementType basis;
    enum componentType type;
    INT_TYPE index;
    double length;
    double origin;
    INT_TYPE grid;
};


struct canon {
    Stream_Type *stream ;
    //vector types
    enum basisElementType basis;
    enum componentType component;
    enum bodyType body;
    enum particleType particle;
    INT_TYPE count1Basis;
    double lattice;
    double origin;
};































enum division{
    nullName,
    kineticMass,
    kineticMass1,
    kineticMass2,
    kineticMass3,
    kineticMass4,
    kinetic,//0
    kinetic1,
    kinetic2,
    kinetic3,
    kinetic4,
    linear,
    external1,
    external2,
    external3,
    external4,//121,122,123,124
    interactionExchangePlus,//oneBody
    interaction1Plus,
    interaction2Plus,
    interaction3Plus,
    interaction4Plus,
    interactionExchangeMinus,//oneBody
    interaction1Minus,
    interaction2Minus,
    interaction3Minus,
    interaction4Minus,
    interactionEwald,
    interaction1Ewald,
    interaction2Ewald,
    interaction3Ewald,
    interaction4Ewald,
    shortenEwald,
    shorten1Ewald,
    shorten2Ewald,
    shorten3Ewald,
    shorten4Ewald,
    shortenPlus,
    shorten1Plus,
    shorten2Plus,
    shorten3Plus,
    shorten4Plus,
    shortenMinus,
    shorten1Minus,
    shorten2Minus,
    shorten3Minus,
    shorten4Minus,
    interactionTwoAcrossDimensions,
    interactionTAD11,
    interactionTAD12,
    interactionTAD13,
    interactionTAD14,
    interactionTAD21,
    interactionTAD22,
    interactionTAD23,
    interactionTAD24,
    interactionTAD31,
    interactionTAD32,
    interactionTAD33,
    interactionTAD34,
    interactionTAD41,
    interactionTAD42,
    interactionTAD43,
    interactionTAD44,
    shortTwoAcrossDimensions,
    shortTAD11,
    shortTAD12,
    shortTAD13,
    shortTAD14,
    shortTAD21,
    shortTAD22,
    shortTAD23,
    shortTAD24,
    shortTAD31,
    shortTAD32,
    shortTAD33,
    shortTAD34,
    shortTAD41,
    shortTAD42,
    shortTAD43,
    shortTAD44,
    protonRepulsion,
    proton1,
    proton2,
    proton3,
    proton4,
    X,//
    X1,
    X2,
    X3,
    X4,
    harmonium,
    harmonium1,
    harmonium2,
    harmonium3,
    harmonium4,//98,99,100,101
    vectorMomentum,
    vectorMomentum1,
    vectorMomentum2,
    vectorMomentum3,
    vectorMomentum4,
    interactionExchange,//22
    interaction12,
    interaction13,
    interaction23,
    interaction14,
    interaction24,
    interaction34,
    interactionExchangeB,//22
    interaction12B,
    interaction13B,
    interaction23B,
    interaction14B,
    interaction24B,
    interaction34B,
    edgeElectronMatrix,//63
    edgeEMatrix1,//63
    edgeEMatrix2,//63
    edgeEMatrix3,//63
    edgeEMatrix4,//63
    edgeProtonMatrix,//63
    edgePMatrix1,//63
    edgePMatrix2,//63
    edgePMatrix3,//63
    edgePMatrix4,//63
    hartree,
    nullScalar,//7
    nullVector,//8
    nullMatrix,//9
    forces,//10
    inversion,
    overlap,
    overlap1,
    overlap2,
    overlap3,
    overlap4,
    resolveBufferMatrix,//11
    distanceBufferVector,//12
    distanceBufferMatrix,//13
    distanceBufferScalar,//14
    maxBufferVector,//15
    squareVector,//16
    diagonal,//18
    square,//19
    project,//20
    interactionDirect,//21
    buffer,//24,
    oneByOneBuffer,
    Ha, //25
    Iterator,
    printOperator,
    copy,//17
    copyTwo , //23
    copyThree,//51
    copyFour,//52
    copyFive,
    eigen,//64
    bill1,
    bill2,
    bill3,
    bill4,
    bill5,
    bill6,
    PauliX,
    PauliY,
    PauliZ,//73,74,75
    trainQuartic,//76,77,78
    cycleVector,
    cycleMatrix,
    cycleQuartic,//79,80,81
    copySix,//82
    directBak, //83
    oneBody,
    entropyVector,
    entropyUnit,//84,85
    vectors3,//90
    basis,
    subBasis,//92,93
    edges,//94
    edges1,
    edges2,
    edges3,
    edgesA1,
    edgesA2,
    edgesA3,
    quad,
    quad2,//102,103
    bandBasis,//105
    quadCube,
    quad2Cube,
    diagonalCube,//106,107,108
    canonicalBuffers,
    copyVector,
    copyTwoVector,
    copyThreeVector,//113
    copyFourVector,
    productVector,
    permutationVector,
    permutation2Vector,
    permutation3Vector,
    permutation4Vector,
    build,
    squareTwo,//
    foundationBasis,
    basisBuffers,
    tensorBuffers,
    tensorBuffers2,
    tensorBuffers3,
    tensorBuffers4,
    tensorBuffers5,
    tensorBuffers6,
    guideBuffer,
    trackBuffer,
    vectorCubeBuffers,//125,126,127
    diagonalQuad,
    diDiagonalQuad,
    diagonal2VectorA,
    diagonal2VectorB,
    diagonal1VectorA,
    diagonalVectorA,
    diagonalVectorB,
    diagonal1VectorB,
    diagonal1VectorC,
    diagonal1VectorD,
    canonicalBuffersB,
    canonicalBuffersBM,
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
    complement,//null
    complement1,
    complement2,
    complement3,
    complementTwo,//null
    complementTwo1,
    complementTwo2,
    complementTwo3,
    eigenList,
    canonicalBuffersB2,
    oneVector,
    twoVector,
    twoBodyRitz,
    conditionOverlapNumbers,
    manyElectronsRitz,//
    foundationStructure,
    foundationEquals,
    totalVector,
    totalFuzzyVector,
    squareTwoVector,//
    dsyBuffers,
    oneArray,
    oneBasis,
    threeArray,
    seconds,
    MomentumDot,
    matrixSbuild,
    matrixHbuild,
    vectorHbuild,
    twoBodyProjector,
    lanes,
    lanes2,
    lanes3,
    lanes4,
    lanes5,
    lanes6,
    lanes7,
    lanes8,
    lanes9,
    lanes10,
    eigenVectors
};

struct space_label{
    ADDRESS_TYPE Address;
    
    //matrix types
    enum bodyType body;
    enum block block;
    //not vector types
};

struct value_label{
    double value;
    double value2;
    INT_TYPE symmetry;
};

struct name_label {
    enum division name;//pointer iff name != name
    enum genus species;
    enum shape header;
    enum spinType spinor;
    enum division linkNext;
    INT_TYPE Partition;
    enum memoryType memory;
    struct value_label value;
    struct space_label space[SPACE+1];
    INT_TYPE Current[MaxCore];//if name != name , then Current is a pointer.
};

struct sinc_label {
    INT_TYPE maxEV;
    double d;//input lattice spacing
    struct canon rose[SPACE+1];
    struct name_label *tulip;//vectors
    enum division user;
    enum division vectorOperator;
    enum division end;

};

struct function_label{
    INT_TYPE interval;
    enum   functionType fn;
    double param[MaxParamFunc];
};

struct interaction_label{
    INT_TYPE num;
    struct function_label func;
};

struct field {
    struct atom atoms[MAXATOM+1];
    INT_TYPE Na;
    INT_TYPE Ne;
    enum bodyType body;
	struct sinc_label sinc;
    struct interaction_label twoBody;
    struct interaction_label oneBody;
    struct MEM * mem1;
};

struct general_index{//one dimension
    struct basisElement bra;
    struct basisElement ket;
    INT_TYPE action ;//0 = no derivative.  1 derivative.
    double x0;//gaussain1 position in 1-d
    INT_TYPE pointer;

    //line in sand...
    double b0;
    double b1;
    INT_TYPE l0;
    INT_TYPE l1;
    double x1;
    
    
	double d;
	double n;//function index.
	double m;
    
};

struct general_2index{//one dimension
	struct general_index i[2];
    double beta;
    double momentumShift;
    INT_TYPE point;
    INT_TYPE powSpace;//r^2*pow in Gaussian term!!
    INT_TYPE gaussianAccelerationFlag;//pre-integrated selection
    INT_TYPE realFlag;//1 == real , else == complex
    
    //for element calculations
    struct function_label * fl;
};

struct runTime {
    INT_TYPE runFlag;
#ifdef OMP
    INT_TYPE NLanes;
    INT_TYPE NSlot;

    double position[MaxCore][6];
#endif
    
#ifdef MKL
    INT_TYPE NParallel;
#endif
    
    INT_TYPE monteCarlo;
    INT_TYPE samples;
    INT_TYPE printFlag;

    double maxEntropy;
    double TARGET ;
    double ALPHA ; //condition
    double CANON ; //threshold
    double TOL;
    double vCANON;
    double targetCondition;
    
    enum bodyType body;
    enum calculationType calcType;
    enum phaseType phaseType;
};

struct MEM {
    INT_TYPE bootedMemory;
    INT_TYPE files ;
    char fileList [MAXSTRING][MAXSTRING];
    
    INT_TYPE filesVectorOperator ;
    char fileVectorOperator [MAXSTRING][MAXSTRING];
    struct runTime * rt;
};


struct input {
    double minClamp;
    double maxClamp;
    double orgClamp;
    INT_TYPE complexType ;
    double level;
    INT_TYPE magFlag;
    double  mag ;
    INT_TYPE cycleStep;
    INT_TYPE lookBack;
    INT_TYPE filter;
    int irrep;
#ifdef OMP
    INT_TYPE omp;
#endif

#ifdef MKL
    INT_TYPE mkl;
#endif

    INT_TYPE l2;
    INT_TYPE nTargets;
    double scalar;
    double turn;
    double param1;
    double param2;
    INT_TYPE interval;
    INT_TYPE sectors;
    double attack;
//    INT_TYPE cSA[nSAG*nSAG*nSAG];
    INT_TYPE side;
    
//    INT_TYPE hartreeFockFlag;
    INT_TYPE vectorOperatorFlag;
    INT_TYPE M1;
    double vectorMomentum;
    double springConstant;
    INT_TYPE springFlag;
    INT_TYPE OCSBflag;
    INT_TYPE potentialFlag;
    INT_TYPE RAMmax ;
    INT_TYPE decomposeRankMatrix;
    INT_TYPE iCharge;
	INT_TYPE epi;
    INT_TYPE around;
	double d;
    double D;
    INT_TYPE outputFlag;
    INT_TYPE cycles;
    INT_TYPE Iterations;
    INT_TYPE charge;
    INT_TYPE canonRank;
    INT_TYPE heliumFlag;
 //   INT_TYPE paths;
    INT_TYPE nStates;
    INT_TYPE iRank ;
    INT_TYPE bRank;
    INT_TYPE dRank;
    INT_TYPE qFloor;
    INT_TYPE Angstroms;
    struct field c;

};


struct calculation {
	char name[MAXSTRING];
	struct input i;
    struct runTime rt;
    struct MEM mem;
};

#endif




