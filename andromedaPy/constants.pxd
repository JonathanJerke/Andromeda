###-------------###
###constants.pxd###
###-------------###


#*  Copyright 2021 Jonathan Jerke and Bill Poirier.
#*  We acknowledge the generous support of Texas Tech University,
#*  the Robert A. Welch Foundation, and the Army Research Office.


#*   *   This file is part of Andromeda.

#*   *   Andromeda is free software: you can redistribute it and/or modify
#*   *   it under the terms of the GNU General Public License as published by
#*   *   the Free Software Foundation, either version 3 of the License.

#*   *   Andromeda is distributed in the hope that it will be useful,
#*   *   but WITHOUT ANY WARRANTY; without even the implied warranty of
#*   *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#*   *   GNU General Public License for more details.

#*   *   You should have received a copy of the GNU General Public License
#*   *   along with Andromeda.  If not, see <https://www.gnu.org/licenses/>.

include "system.pxi"

ctypedef int inta
ctypedef double floata
ctypedef floata mea
ctypedef int ADDRESS_TYPE

cdef extern from "../Object/constants.h":


    cdef enum bodyType:
        nada,
        one,
        two,
        three,
        four,
        five,
        six

    cdef enum particleType:
        all

    cdef enum phaseType:
        buildFoundation,
        productKrylov,
        solveRitz,
        buildTraces,
        formOCSB,
        iterateOCSB

    cdef enum calculationType:
        nullCalculation,
        electronicStuctureCalculation

    cdef enum bootType:
        noMatrices,
        fullMatrices

    cdef enum blockType:
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

    cdef enum genusType:
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

    cdef enum functionType:
        nullFunction,
        Pseudo,
        Yukawa,
        Coulomb,
        Morse,
        LennardJones,
        LDA,
        BLYP,
        Gaussian

    cdef enum spinType:
        none,
        real,
        cmpl,
        parallel

    cdef enum memoryType:
        noAllocation,
        objectAllocation,
        bufferAllocation

    cdef enum basisElementType:
        nullBasisElement,
        SincBasisElement,
        GaussianBasisElement,
        DiracDeltaElement,
        StateBasisElement,
        overlapBasisElement

    cdef enum componentType:
        nullComponent,
        spatialComponent1,
        spatialComponent2,
        spatialComponent3,
        periodicComponent1,
        periodicComponent2,
        periodicComponent3,

    cdef enum noteType:
        nullNote,
        interactionCell

    cdef enum metricType:
        dirac,
        separateDirac,
        interval,
        semiIndefinite,
        pureInterval,
        pureSemiIndefinite

    cdef enum blockMemoryType:
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
        blockDiagonalMatrixblock,
        blockCompressionBlock

    cdef enum division:
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
        canonicalBuffersD,
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

    cdef struct function_label:
        inta interval
        inta contr
        functionType fn
        double param[4]

    cdef struct atom_label:
        double position[3+1]
        inta Z

    cdef struct basisElement_label:
        basisElementType basis
        componentType component
        noteType note
        inta index
        inta index2
        double length
        double origin
        inta grid

    cdef struct dimensions_label:
        double lattice
        double attack
        double origin
        double anchor
    
    cdef struct canon_label:
        floata *stream
        basisElementType basis
        componentType component
        bodyType body
        inta label
        inta count1Basis
        inta count1Inc
        inta space
        dimensions_label particle[3+1]

    cdef struct metric_label:
        inta pow[SPACE]
        inta powB[SPACE]
        inta deriv[SPACE]
        function_label fn
        metricType metric
        double beta[2]

    cdef struct term_label:
        char filename[MAXSTRING]
        function_label func
        char desc[16]
        double scalar
        inta type
        blockType bl
        inta act
        double adjustOne
        inta label
        inta invert
        inta headFlag
        inta atom
        inta bra
        inta ket
        inta embed
        metric_label mu

    cdef struct runTime:
    	inta boost
        inta dynamic
        blockMemoryType memBlock[16]
        inta NLanes
        inta NSlot
        inta NParallel
        double TOLERANCE
        double relativeTOLERANCE
        double THRESHOLD
        double XCONDITION
        double ALPHA
        inta MAX_CYCLE
        calculationType calcType
        phaseType phaseType

    cdef struct space_label:
        ADDRESS_TYPE Address
        bodyType body
        blockType block
        inta act
        inta invert
        inta bra
        inta ket

    cdef struct value_label:
        char title[MAXSTRING]
        inta stage
        double value
        double value2
        inta symmetry

    cdef struct name_label:
        division name
        genusType species
        spinType spinor
        division linkNext
        division chainNext
        inta multId
        division loopNext
        inta Partition
        memoryType memory
        value_label value
        space_label space[SPACE+1]
        inta Current[4]
        inta Begin[4]

    cdef struct nameDispenser_label:
        inta currLabel
        inta maxLabel
        division head

    cdef struct sinc_label:
        nameDispenser_label eikonLabels
        nameDispenser_label nullLabels
        bootType boot
        inta maxEV
        inta irrep
        inta bootedMemory
        canon_label canon[SPACE+1]
        name_label *name
        division user
        division vectorOperator
        division end
        spinType cmpl
        runTime * rt
    
    cdef struct input_label:
        inta flex
        inta OpIndex
        inta body
        inta irrep
        inta matrices
        char matrixList [1000][MAXSTRING]
        inta files
        char fileList [1000][MAXSTRING]
        inta filesVectorOperator
        char fileVectorOperator [1000][MAXSTRING]
        inta Iterations
        inta nStates
        inta iRank
        inta canonRank
        inta xRank
        inta qFloor
        inta nOperator
        inta filter
        inta collect
        inta cmpl

    cdef struct field:
        input_label i
        sinc_label f

    cdef struct input:
        inta numNames
        inta numVectors
        inta build
        inta termNumber
        term_label terms[100]
        inta minIterationPrint
        char controlPath[MAXSTRING]
        double shiftVector[2]
        inta shiftFlag
#       inta omp
#       inta mkl
        inta SymmetrizedGaussianLevel
        floata SymmetrizedGaussianWidth
        inta Lambda
        inta RAMmax
        inta Angstroms
        atom_label atoms[12+1]
        inta Na
        inta nocsb
        inta iocsb

    cdef struct calculation:
        char name[MAXSTRING]
        input i
        runTime rt

