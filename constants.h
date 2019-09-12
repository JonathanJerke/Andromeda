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
 
//VERSION 7.3

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

#define CDT 1
#define MaxSpin 2

enum bodyType {
    nada,
    one,
    two,
    three,
    four,
    five,
    six
};

enum particleType{
    nullParticle,
    electron,
    proton,
    pair,
    all,
};

enum phaseType{
    buildFoundation,//0
    productKrylov,//1
    solveRitz,//2
    svdOperation, //3
    distillMatrix,//4
    reportMatrix//5
};

enum calculationType{
    nullCalculation,
    electronicStuctureCalculation,
    clampProtonElectronCalculation,
    protonsElectronsCalculation,
    svdCalculation,
};

enum bootType{
    noMatrices,
    fullMatrices
};

enum block{
    id0,
    //nada
    tv1,
    tv2,
    tv3,
    tv4,
    tv5,
    tv6,
    //ee interaction between label 1 and 2
    e12,
    e13,
    e23,
    //four
    e14,
    e24,
    e34,
    //five
    e15,
    e25,
    e35,
    e45,
    //six
    e16,
    e26,
    e36,
    e46,
    e56
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

enum functionType{
    nullFunction,//0
    Pseudo,//1
    Yukawa,//2
    Coulomb,//3
    Morse,//4
    LennardJones,//5
    LDA,//6
    BLYP,//7
    Gaussian//8
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


enum basisElementType {
    nullBasisElement,
    SincBasisElement,
    GaussianBasisElement,
    DiracDeltaElement
};

enum componentType {
    nullComponent,
    spatialComponent1,//1
    spatialComponent2,//2
    spatialComponent3,//3
    periodicComponent1,//4
    periodicComponent2,//5
    periodicComponent3,//6
    periodicSumComponent1,//7
    periodicSumComponent2,//8
    periodicSumComponent3,//9
    periodicBoostComponent1,//10
    periodicBoostComponent2,//11
    periodicBoostComponent3,//12
};

enum noteType {
    nullNote,//0
    interactionCell//1
};

struct basisElement {
    enum basisElementType basis;
    enum componentType type;
    enum noteType note;
    INT_TYPE index;
    INT_TYPE index2;//boost for periodicBoostComponent##
    double length;
    double origin;
    INT_TYPE grid;
};

struct canon {
    Stream_Type *stream ;
    struct basisElement *basisList;
    //vector types
    enum basisElementType basis;
    enum componentType component;
    enum bodyType body;
    enum particleType particle;
    INT_TYPE count1Basis;
};

enum metricType {
    dirac,
    separateDirac,
    interval,
    semiIndefinite,
};

struct metric_label {
    double pow[SPACE];
    double powB[SPACE];
    double deriv[SPACE];
    struct function_label fn;
    enum metricType metric;
    double beta[2];//lower and upper bound
    //beta here...  -beta^2 is the exponent
};


//blockMemory:
enum blockMemoryType{
    passBlock,//
    blockHamiltonianBlock,//1
    blockTrainHamiltonianBlock,//2
    blockTrainingHamiltonianBlock,//3
    blockFoundationBlock,//4
    blockBuildHamiltonianBlock,//5
    blockEigenDecomposeBlock//6
};































































enum division{
    nullName,
    nullScalar,//7
    nullVector,//8
    nullMatrix,//9
    hamiltonian,
    trainHamiltonian,
    h12,
    h13,
    h23,
    h14,
    h24,
    h34,//five
    h15,
    h25,
    h35,
    h45,//six
    h16,
    h26,
    h36,
    h46,
    h56,
	kineticMass,
    kineticMass1,
    kineticMass2,
    kineticMass3,
    kineticMass4,
    kineticMass5,
    kineticMass6,
    kinetic,//0
    kinetic1,
    kinetic2,
    kinetic3,
    kinetic4,
    kinetic5,
    kinetic6,
    linear,
    external1,
    external2,
    external3,
    external4,//121,122,123,124
    external5,
    external6,
    interactionExchangePlus,//oneBody
    interaction1Plus,
    interaction2Plus,
    interaction3Plus,
    interaction4Plus,
    interaction5Plus,
    interaction6Plus,
    interactionExchangeMinus,//oneBody
    interaction1Minus,
    interaction2Minus,
    interaction3Minus,
    interaction4Minus,
    interaction5Minus,
    interaction6Minus,
    interactionEwald,
    interaction12Ewald,
    interaction13Ewald,
    interaction23Ewald,
    interaction14Ewald,
    interaction24Ewald,
    interaction34Ewald,
    interaction15Ewald,
    interaction25Ewald,
    interaction35Ewald,
    interaction45Ewald,
    interaction16Ewald,
    interaction26Ewald,
    interaction36Ewald,
    interaction46Ewald,
    interaction56Ewald,
    intercellularSelfEwald,
    intercellularSelf1Ewald,
    intercellularSelf2Ewald,
    intercellularSelf3Ewald,
    intercellularSelf4Ewald,
    intercellularSelf5Ewald,
    intercellularSelf6Ewald,
    intracellularSelfEwald,
    intracellularSelf1Ewald,
    intracellularSelf2Ewald,
    intracellularSelf3Ewald,
    intracellularSelf4Ewald,
    intracellularSelf5Ewald,
    intracellularSelf6Ewald,
    jelliumElectron,
    jellium1Electron,
    jellium2Electron,
    jellium3Electron,
    jellium4Electron,
    jellium5Electron,
    jellium6Electron,
    shortenPlus,
    shorten1Plus,
    shorten2Plus,
    shorten3Plus,
    shorten4Plus,
    shorten5Plus,
    shorten6Plus,
    shortenMinus,
    shorten1Minus,
    shorten2Minus,
    shorten3Minus,
    shorten4Minus,
    shorten5Minus,
    shorten6Minus,
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
    interactionTAD15,
    interactionTAD25,
    interactionTAD35,
    interactionTAD45,
    interactionTAD16,
    interactionTAD26,
    interactionTAD36,
    interactionTAD46,
    interactionTAD56,
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
    shortTAD15,
    shortTAD25,
    shortTAD35,
    shortTAD45,
    shortTAD16,
    shortTAD26,
    shortTAD36,
    shortTAD46,
    shortTAD56,
    protonRepulsion,
    proton1,
    proton2,
    proton3,
    proton4,
    proton5,
    proton6,
    X,//
    X1,
    X2,
    X3,
    X4,
    X5,
    X6,
    harmonium,
    harmonium1,
    harmonium2,
    harmonium3,
    harmonium4,//98,99,100,101
    harmonium5,
    harmonium6,
    vectorMomentum,
    vectorMomentum1,
    vectorMomentum2,
    vectorMomentum3,
    vectorMomentum4,
    vectorMomentum5,
    vectorMomentum6,
    interactionExchange,//22
    interaction12,
    interaction13,
    interaction23,
    interaction14,
    interaction24,
    interaction34,
    interaction15,
    interaction25,
    interaction35,
    interaction45,
    interaction16,
    interaction26,
    interaction36,
    interaction46,
    interaction56,
    interactionExchangeB,//22
    interaction12B,
    interaction13B,
    interaction23B,
    interaction14B,
    interaction24B,
    interaction34B,
    interaction15B,
    interaction25B,
    interaction35B,
    interaction45B,
    interaction16B,
    interaction26B,
    interaction36B,
    interaction46B,
    interaction56B,
    edgeElectronMatrix,//63
    edgeEMatrix1,//63
    edgeEMatrix2,//63
    edgeEMatrix3,//63
    edgeEMatrix4,//63
    edgeEMatrix5,//63
    edgeEMatrix6,//63
    edgeProtonMatrix,//63
    edgePMatrix1,//63
    edgePMatrix2,//63
    edgePMatrix3,//63
    edgePMatrix4,//63
    edgePMatrix5,//63
    edgePMatrix6,//63
    hartree,
    forces,//10
    inversion,
    overlap,
    overlap1,
    overlap2,
    overlap3,
    overlap4,
    overlap5,
    overlap6,
    tempOneMatrix,
    tempTwoMatrix,
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
    canonicalBuffers0,
    copyVector,
    copyTwoVector,
    copyThreeVector,//113
    copyFourVector,
    scalarTemp,
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
    guideBuffer0,
    trackBuffer0,
    vectorCubeBuffers,//125,126,127
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
    eigenVectors
};

//struct divisionRange {
//    enum division begin;
//    enum division end;
//};


struct space_label{
    ADDRESS_TYPE Address;
    
    //matrix types
    enum bodyType body;
    enum block block;
    //not vector types
};

struct value_label{
    char title [MAXSTRING];
    INT_TYPE stage;
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
    enum bootType boot;
    INT_TYPE maxEV;
    INT_TYPE irrep;
    INT_TYPE cat;
    INT_TYPE bootedMemory;
    struct canon rose[SPACE+1];
    struct name_label *tulip;//vectors
    enum division user;
    enum division vectorOperator;
    enum division end;
    enum division purity;

#ifdef PURITY
    enum division purityOverlap;
    enum division temp;
    INT_TYPE purityCanon;
#endif
    enum spinType cmpl;
    struct runTime * rt;
};


struct input_label {
    INT_TYPE body;
    INT_TYPE irrep;
    INT_TYPE cat;
    INT_TYPE files ;
    INT_TYPE filesVectorOperator ;
    char fileVectorOperator [MAXFILE][MAXSTRING];
    char fileList [MAXFILE][MAXSTRING];
    INT_TYPE epi;
    INT_TYPE around;
    double d;
    double D;
    double attack;
    INT_TYPE Iterations;
    INT_TYPE nStates;
    INT_TYPE iRank ;
    INT_TYPE bRank;
    INT_TYPE xRank;
    INT_TYPE qFloor;
    INT_TYPE nOperator;
    INT_TYPE filter;
    INT_TYPE collect;
    INT_TYPE cmpl;
};


struct field {
    struct input_label i;
    struct sinc_label f;
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
    enum spinType realFlag;//1 == real , else == complex

    
    double d;
    double n;//function index.
    double m;
    
};

struct general_2index{//one dimension
	struct general_index i[2];
    double beta;
    double momentumShift;
    INT_TYPE point;
    INT_TYPE pow2[2];
    INT_TYPE powSpace;//r^2*pow in Gaussian term!!
    INT_TYPE gaussianAccelerationFlag;//pre-integrated selection
    enum spinType realFlag;//1 == real , else == complex
    
    //for element calculations
    struct function_label * fl;
};

struct runTime {
    INT_TYPE powDecompose;
    INT_TYPE runFlag;
    enum blockMemoryType memBlock[BlockCount];

#ifdef OMP
    INT_TYPE NLanes;
    INT_TYPE NSlot;
    double position[MaxCore][6];
#else
    INT_TYPE NLanes;
#endif
    
#ifdef MKL
    INT_TYPE NParallel;
#endif
//    INT_TYPE samples;

    double maxEntropy;
    double TARGET ;
    double ALPHA ; //condition
    double CANON ; //threshold
    double TOL;
    double vCANON;
    double EWALD;
    double targetCondition;
    double CAP;
    enum calculationType calcType;
    enum phaseType phaseType;
};

struct input {
    double shiftVector[100][2];
    INT_TYPE barrier;
    double massElectron;//electron
    double massProton;//protons
    double massClampPair;//clamped second particle
    INT_TYPE shiftFlag ;
    double realPart ;
    double minClamp;
    double maxClamp;
    double orgClamp;
    double level;
    INT_TYPE magFlag;
    double  mag ;
#ifdef OMP
    INT_TYPE omp;
#endif

#ifdef MKL
    INT_TYPE mkl;
#endif
    INT_TYPE canonRank;
    double scalar;
    double turn;
    double param1;
    double param2;
    INT_TYPE interval;
    INT_TYPE M1;
    double vectorMomentum;
    double springConstant;
    INT_TYPE springFlag;
  //  INT_TYPE OCSBflag;
    INT_TYPE RAMmax ;
  //  INT_TYPE bootRestriction;
    INT_TYPE decomposeRankMatrix;
    INT_TYPE Angstroms;
    struct atom atoms[MAXATOM+1];
    INT_TYPE Na;
    struct interaction_label twoBody;
    struct interaction_label oneBody;
};


struct calculation {
	char name[MAXSTRING];
	struct input i;
    struct runTime rt;
};

#endif




