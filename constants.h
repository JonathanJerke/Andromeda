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


#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "system.h"
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define VERBOSE 0 

#ifdef APPLE
#include "Accelerate/Accelerate.h"
#include "complex.h"
typedef double __complex__ DCOMPLEX;
typedef __CLPK_integer INT_TYPE;
typedef __CLPK_doublereal Stream_Type;
typedef __CLPK_integer ADDRESS_TYPE;

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

#define hardMaxCore 72
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

//////////////  Quick and dirty defintion which is appropriately FAST to keep memory and in buffer...
/////////
////
//2147483647

#define CDT 'T'
#define MaxSpin 2
#define nSAG 5

//pseudo-potentials

struct canon {
    Stream_Type *stream ;// C_F*F
};

//enum symmetry {
//    nullSymmetry,
//    symmetricalMatrices,
//};
//
enum body {
    nada,
    one,
    two,
    three,
    four
};

enum bodyType{
    electron,
    proton,
    h2plus,//for labeling purposes
    h2//for labeling purposes
};

enum block{
    tv1,//kinetic + potential block 1
    tv2,//1-b block 2
    tv3,
    tv4,
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
    cubic,
    quartic,
};

enum shape{
    Cube,
    Rds,
    RdsBasis,
    Band,
    BandBasis,
};


enum purposeType{
    nullObject,
    Object,//trainable
    tObject,//sets to train to Object
    vObject,//in stack, for building new vectors...unset to nullify list element
    ptObject//views / pointers to Obj
};

enum memoryType{
    threeObject,
    oneObject
};


struct atom_label {
    INT_TYPE iZ;

    INT_TYPE Z;
    INT_TYPE Cc;
    INT_TYPE label;
};

struct atom {
	double position[4];
	struct atom_label label;
};

#define NspinType 2
enum spinType {
    none,
    cmpl
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
    DiracDelta
};


struct basisElement {
    INT_TYPE periodic;
    enum basisElementType basis;
    INT_TYPE index;
    double length;
    double origin;
    
    INT_TYPE auxIndex; //for periodic Sincs
};




















































enum division{
    kineticMass,
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
    interactionExchange,//22
    interaction12,
    interaction13,
    interaction23,
    interaction14,
    interaction24,
    interaction34, //-157
    interactionExchangePlus,
    interaction12Plus,
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
    hartree,
    nullScalar,//7
    nullVector,//8
    nullMatrix,//9
    forces,//10
    overlap,
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
    Ha, //25
    Iterator,
    printOperator,
    highBallVector,
    highBallVector2,
    highBallVector3,
    highBallVector4,
    secondVector,//48
    maxBuffer4Vector,//49
    maxBuffer5Vector,//50
    copy,//17
    copyTwo , //23
    copyThree,//51
    copyFour,//52
    copyFive,
    edgeMatrix,//63
    eigen,//64
    trainVector,
    trainVector2,
    trainVector3,
    trainVector4,
    trainMatrix,
    trainMatrix2,
    trainMatrix3,
    trainMatrix4,
    spBuffer,//71,72
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
    init,
    initValence,
    initLumos ,//95,96,97
    quad,
    quad2,//102,103
    bandBasis,//104
    nullName,//105
    quadCube,
    quad2Cube,
    diagonalCube,//106,107,108
    canonicalBuffers,
    eigenBuffers,//109,110
    copyVector,
    copyCube,
    diagonalVector,
    copyTwoVector,
    copyThreeVector,//113
    copyFourVector,
    productVector,
    permutationVector,
    build,
    squareTwo,//
    foundationBasis,
    diDiagonalVector,
    diDiagonal,//118,119,120
    basisBuffers,
    tensorBuffers,
    tensorBuffers2,
    tensorBuffers3,
    tensorBuffers4,
    tensorBuffers5,
    tensorBuffers6,
    vectorCubeBuffers,//125,126,127
    diagonalQuad,
    diDiagonalQuad,
    rtDensity,
    build2,//128,129,130,131
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
    canonicalBuffersBX,
    canonicalBuffersC,
    canonicalBuffersD,//
    complement,
    complementTwo,
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
    lanes11,lanes12,lanes13,lanes14,lanes15,lane16,lanes17,lanes18,lanes19,lanes20,lanes21,lanes22,lanes23,lanes24,lanes25,lanes26,lanes27,lanes28,lanes29,lanes30,lanes31,lane32,lanes33,lanes34,lanes35,lanes36,
    lanesb,lanesb2,lanesb3,lanesb4,lanesb5,lanesb6,lanesb7,lanesb8,lanesb9,lanesb10,lanesb11,lanesb12,lanesb13,lanesb14,lanesb15,lanesb16,lanesb17,lanesb18,lanesb19,lanesb20,lanesb21,lanesb22,lanesb23,lanesb24,lanesb25,lanesb26,lanesb27,lanesb28,lanesb29,lanesb30,lanesb31,lanesb32,lanesb33,lanesb34,lanesb35,lanesb36,lanesc,lanesc2,lanesc3,lanesc4,lanesc5,lanesc6,lanesc7,lanesc8,lanesc9,lanesc10,lanesc11,lanesc12,lanesc13,lanesc14,lanesc15,lanesc16,lanesc17,lanesc18,lanesc19,lanesc20,lanesc21,lanesc22,lanesc23,lanesc24,lanesc25,lanesc26,lanesc27,lanesc28,lanesc29,lanesc30,lanesc31,lanesc32,lanesc33,lanesc34,lanesc35,lanesc36,lanesd,lanesd2,lanesd3,lanesd4,lanesd5,lanesd6,lanesd7,lanesd8,lanesd9,lanesd10,lanesd11,lanesd12,lanesd13,lanesd14,lanesd15,lanesd16,lanesd17,lanesd18,lanesd19,lanesd20,lanesd21,lanesd22,lanesd23,lanesd24,lanesd25,lanesd26,lanesd27,lanesd28,lanesd29,lanesd30,lanesd31,lanesd32,lanesd33,lanesd34,lanesd35,lanesd36,

    
    
    
    
    eigenVectors
};


struct name_label {
    INT_TYPE sigN[3];
    INT_TYPE Partition;
    INT_TYPE Current[MaxCore];
    INT_TYPE ptRank[MaxCore];
    ADDRESS_TYPE Address;
    ADDRESS_TYPE myAddress;
    INT_TYPE path ;
    enum body NBody;
    enum bodyType TBody;
    double value;
    double value2;
    double value3;

    INT_TYPE parallel;
    INT_TYPE stop[MaxSpin][MAXATOM+1];//want to get rid of this scaling monster
    int symmetry;
    int symmetry2;
    enum division name;
    enum genus species;
    enum shape header;
    enum purposeType purpose;
    enum memoryType memory;
    enum block blockType;
    enum spinType spinor;
    enum division linkNext;
};

struct sinc_label {
    double volume;
    INT_TYPE maxEV;

    INT_TYPE N1;
    INT_TYPE dims[7][3];
    double d;//input lattice spacing
    double dd[7];//multiplying global lattice length, d, at fixed ratios.
    struct canon rose[SPACE+1];
    struct name_label *tulip;//many names,  Myrid
    enum division density;
    enum division end;
    INT_TYPE Basis[3];
    INT_TYPE Basis2[3];
    INT_TYPE Basis3[3];
    INT_TYPE Basis4[3];
};

struct rds_label{
    INT_TYPE flag;
    INT_TYPE Basis[3];
    INT_TYPE Basis2[3];
    INT_TYPE Basis3[3];
    INT_TYPE Basis4[3];
};

//struct band_label{
//    INT_TYPE Basis[3];
//    INT_TYPE Basis2[3];
//    INT_TYPE Basis3[3];//unused
//    INT_TYPE Basis4[3];//unused
//    double D;
//};

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
    enum body body;
    INT_TYPE ir;
	struct sinc_label sinc;
    struct rds_label rds;
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
    double alpha;
    double momentumShift;
    INT_TYPE point;
    INT_TYPE powSpace;//r^2*pow in Gaussian term!!
    INT_TYPE N1;
    INT_TYPE body;
    INT_TYPE periodic;//1==periodic basis//2==periodic external field
    INT_TYPE gaussianAccelerationFlag;
    INT_TYPE realFlag;
    struct function_label * fl;
};

struct runTime {
    INT_TYPE runFlag;

#ifdef OMP
    INT_TYPE NCore;
    INT_TYPE NSlot;

    double position[MaxCore][6];
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
    
    enum body body;
    enum bodyType bodyType;

};

struct MEM {
    INT_TYPE bootedMemory;
    char fileList [MAXSTRING];
    char constraintFile [MAXSTRING];
    char densityName[MAXSTRING];
    struct runTime * rt;
};


struct input {
    INT_TYPE theory;//1 normal schrodinger.  2 N-body theory.
    INT_TYPE cycleStep;
    INT_TYPE lookBack;
    INT_TYPE filter;
    int irrep;
#ifdef OMP
    INT_TYPE omp;
#endif
    double seekPower;
    INT_TYPE l2;
    INT_TYPE group;
    INT_TYPE nTargets;
    double scalar;
    double turn;
    double param1;
    double param2;
    INT_TYPE interval;
    INT_TYPE sectors;
    double attack;
    INT_TYPE cSA[nSAG*nSAG*nSAG];
    INT_TYPE side;
    
    INT_TYPE hartreeFockFlag;
    INT_TYPE densityFlag;
    enum body bodyDensity;
    INT_TYPE M1;
    double vectorMomentum;
    double springConstant;
    INT_TYPE springFlag;
    INT_TYPE potentialFlag;
    INT_TYPE RAMmax ;
    INT_TYPE decomposeRankMatrix;
    INT_TYPE iCharge;
	struct field c;
	INT_TYPE epi;
	double d;
    INT_TYPE outputFlag;
    INT_TYPE cycles;
    INT_TYPE Iterations;
    INT_TYPE charge;
    INT_TYPE canonRank;
    INT_TYPE heliumFlag;
    INT_TYPE paths;
    INT_TYPE nStates;
    INT_TYPE iRank ;
    INT_TYPE bRank;
    INT_TYPE dRank;
    INT_TYPE qFloor;
    INT_TYPE Angstroms;
};


struct calculation {
	char name[MAXSTRING];
    char cycleName[MAXSTRING];
	struct input i;
    struct runTime rt;
    struct MEM mem;
};

#endif




