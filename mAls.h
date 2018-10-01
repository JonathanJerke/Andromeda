/*
 *  mAls.h
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

#ifndef mAls_h
#define mAls_h
#include "constants.h"
#include "coreUtil.h"
#include "coreMath.h"
#include "interfaceMath.h"
double tCycleMultiplyMP ( INT_TYPE rank,struct field * f1 , char mc,enum division mat,INT_TYPE ms, enum division vec, INT_TYPE vs,   enum division alloy ,  INT_TYPE spin ,double tolerance, INT_TYPE maxRun, double power);
double tInnerVectorListMP( INT_TYPE rank, struct field * f1 , enum division origin, double * coeff, enum division vector,INT_TYPE spin );
double tInnerListMP( INT_TYPE rank, struct field * f1 , enum division origin, double * coeff );
double matrixElements ( INT_TYPE rank,struct field * f1 , enum division uec,INT_TYPE us, char mc, enum division mat,INT_TYPE ms, enum division vec,INT_TYPE vs);
double canonicalMultiplyMP( INT_TYPE rank,struct field * f1 , char mc,enum division mat,INT_TYPE ms, enum division vec, INT_TYPE vs,   enum division alloy ,  INT_TYPE spin ,double tolerance);
INT_TYPE normalize (struct field * f1,  enum division alloy, INT_TYPE spin, INT_TYPE space);
INT_TYPE spread (struct field * f1, enum division origin, INT_TYPE os, enum division alloy, INT_TYPE spin, INT_TYPE space, Stream_Type * output,Stream_Type * output2);
INT_TYPE balance (struct field * f1,  enum division alloy, INT_TYPE spin);
double canonicalDecompositionMP( INT_TYPE rank,struct field * f1 , enum division origin,INT_TYPE os,   enum division alloy ,  INT_TYPE spin ,double tolerance);
double canonicalListDecompositionMP( INT_TYPE rank,struct field * f1 , Stream_Type * cofact, enum division origin,INT_TYPE os,   enum division alloy ,  INT_TYPE spin ,double tolerance);
double tCycleDecompostionOneMP ( INT_TYPE rank, struct field * f1 , enum division origin, INT_TYPE os, enum division alloy,INT_TYPE spin,  double tolerance , INT_TYPE maxRun , double power  );
double tCycleDecompostionListOneMP ( INT_TYPE rank, struct field * f1 , enum division origin, double * coeff, enum division alloy,INT_TYPE spin,  double tolerance , INT_TYPE maxRun , double power  );
INT_TYPE tOuterProductSu( struct field * f1,enum division vector , INT_TYPE a, enum division vector2,INT_TYPE b, enum division proj, INT_TYPE c);
//double tHermDev(INT_TYPE rank,  struct field * f1,  enum division H, double dt, double X0,enum division right);
double tMultiplyMP(INT_TYPE rank,  INT_TYPE * info,struct field * f1,double scale, INT_TYPE beta,  enum division equals, INT_TYPE espin ,char leftChar, enum division left, INT_TYPE lspin, char rightChar,enum division right, INT_TYPE rspin);
double distanceOne(INT_TYPE rank,struct field * f1 , enum division alloy , INT_TYPE spin , enum division alloyBak, INT_TYPE spin2);
double inner(INT_TYPE rank,struct field * f1 , enum division alloy , INT_TYPE spin );
double magnitude ( struct field * f1 , enum division alloy );
void tHXpX (  INT_TYPE rank, struct field * f1 , enum division left,INT_TYPE shiftFlag, double product, double productCmpl, enum division equals ,  double tolerance , INT_TYPE maxRun  );
double positioningElectrons2 (INT_TYPE rank, struct field * f1 , enum division oneVector, enum division wavefunction,double x1, double y1, double z1, double x2, double y2, double z2);
INT_TYPE ready ( struct calculation * c );
INT_TYPE tConstructDensity(struct calculation * calc , INT_TYPE ct );
#endif /* mAls_h */
