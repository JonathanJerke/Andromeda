/*
 *  ioPrint.h
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

#ifndef ioPrint_h
#define ioPrint_h
#include "mAls.h"
#include"input.h"
#include"Model.h"
struct fieldArray {
    struct field * f1;
    double radius[MaxCore];
    double displacement[MaxCore][SPACE];
    INT_TYPE class;
    INT_TYPE number;
};
INT_TYPE print(struct calculation *c,INT_TYPE lv,enum division eigenVectors );
INT_TYPE printVector (struct calculation *c,char * name,char * vectorName,INT_TYPE iv, INT_TYPE irrep, DCOMPLEX * vector);
double evaluateDensityBracket( double x [], size_t dim , void * params );
double evaluateVectorBracket( double x [], size_t dim , void * params );
INT_TYPE tSymmetryClass (char * cycleName, char * read , char * filename,INT_TYPE cmplFlag, INT_TYPE cmpl);
void tFromReadToToken (char * read , char * token);
void outputFormat(struct field  * f1, FILE * out, enum division output ,INT_TYPE spin);
double inputFormat(struct field * f1,char * name,  enum division buffer, enum division input);
double tComputeRadialPlot(struct field * f1, INT_TYPE number, INT_TYPE class,  double radius,INT_TYPE numC );
double tComputeVectorPlot(struct field * f1,INT_TYPE number,  INT_TYPE class,  double *displacement,INT_TYPE numC );
INT_TYPE ioStoreMatrix(struct field * f1, enum division op, INT_TYPE spin, char * filename, INT_TYPE ioIn );
INT_TYPE printOutput ( struct field * f1,INT_TYPE number);
INT_TYPE printVectorOutput ( struct field * f1,INT_TYPE number);
INT_TYPE printFaceOutput ( struct field * f1,INT_TYPE number);

INT_TYPE tFillBasis(Stream_Type ** pt/*3 vectors*/, double * coordinates/*3 numbers*/, INT_TYPE class,INT_TYPE N1,double lattice);
INT_TYPE tLoadEigenWeights ( struct calculation * c1, char * filename , enum division input);
INT_TYPE tLoadEigenWeightsWithConstraints (struct calculation * c1, char * filename, char * constraints );

void tFilename (char * cycleName, INT_TYPE count, INT_TYPE body ,INT_TYPE IRREP, INT_TYPE cmpl, char * filename);
DCOMPLEX tFromReadToFilename (char * cycleName, char * read , char * filename,INT_TYPE cmplFlag, INT_TYPE cmpl);

#endif /* ioPrint_h */
