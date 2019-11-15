/*
 *  mAls.h
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

#ifndef mAls_h
#define mAls_h
#include "constants.h"
#include "coreUtil.h"
#include "coreMath.h"
#include "interfaceMath.h"
#include "saUtil.h"
#include "eigen.h"
INT_TYPE completeOverlap (INT_TYPE rank, struct sinc_label  f1, INT_TYPE dim,enum division vector,INT_TYPE v,INT_TYPE spin, enum division ov , INT_TYPE v2,INT_TYPE sp2);
INT_TYPE analyze (struct sinc_label  f1 , enum division term,INT_TYPE sp );
void pContract ( INT_TYPE rank,struct sinc_label  *f1, enum division mat,INT_TYPE ms, enum division vector ,INT_TYPE vs1, enum division vector2,INT_TYPE vs2);
void tContract ( INT_TYPE rank,struct sinc_label  f1, enum division mat,INT_TYPE ms, enum division vector ,INT_TYPE vs1, enum division vector2,INT_TYPE vs2);
double tCycleMultiplyMP ( INT_TYPE rank,struct sinc_label  f1 ,INT_TYPE begin, INT_TYPE end, char mc,enum division mat,INT_TYPE ms, enum division vec, INT_TYPE vs,   enum division alloy ,  INT_TYPE spin ,double tolerance, INT_TYPE maxRun, double power);
double tInnerVectorListMP( INT_TYPE rank, struct sinc_label  f1 , enum division origin, double * coeff, enum division vector,INT_TYPE spin );
double tInnerListMP( INT_TYPE rank, struct sinc_label  f1 , enum division origin, double * coeff );
void matrixElements ( INT_TYPE rank,struct sinc_label f1 , enum division bra, enum division mat, enum division ket,DCOMPLEX *ME,DCOMPLEX *OV );
void pMatrixElements ( struct sinc_label  f1 , enum division bra, enum division mat, enum division ket,DCOMPLEX *ME,DCOMPLEX *OV );
void pOverlap ( INT_TYPE rank,struct sinc_label  f1 , enum division bra,INT_TYPE b1, INT_TYPE b2,INT_TYPE sp, enum division ket,INT_TYPE k1, INT_TYPE k2,INT_TYPE sp2,DCOMPLEX *OV );
double canonicalMultiplyMP( INT_TYPE rank,struct sinc_label  f1 , INT_TYPE begin, INT_TYPE end,char mc,enum division mat,INT_TYPE ms, enum division vec, INT_TYPE vs,   enum division alloy ,  INT_TYPE spin ,double tolerance);
INT_TYPE normalize (struct sinc_label  f1,  enum division alloy,INT_TYPE l3,INT_TYPE l4, INT_TYPE spin, INT_TYPE space);
double distanceFrac1 (struct sinc_label  f1 ,enum division alloy, INT_TYPE a1,INT_TYPE a2,INT_TYPE os, enum division alloyBak,INT_TYPE b1, INT_TYPE b2,INT_TYPE os2);
INT_TYPE spread (struct sinc_label  f1, enum division origin, INT_TYPE l1,INT_TYPE l2,INT_TYPE os, enum division alloy,INT_TYPE l3 , INT_TYPE l4, INT_TYPE spin, INT_TYPE space, Stream_Type * output,Stream_Type * output2);
INT_TYPE pSpread (struct sinc_label f1, enum division origin, INT_TYPE l1, INT_TYPE l2, INT_TYPE os, enum division alloy, INT_TYPE l3, INT_TYPE l4, INT_TYPE spin, INT_TYPE space, Stream_Type * output,Stream_Type * output2);
INT_TYPE balance (struct sinc_label f1,  enum division alloy, INT_TYPE spin);
double canonicalListDecompositionMP( INT_TYPE rank,struct sinc_label f1 , Stream_Type * cofact, enum division origin,INT_TYPE os,   enum division alloy ,  INT_TYPE spin ,double tolerance,double magn,INT_TYPE preferred);
double tCycleDecompostionChromaticOneMP ( struct sinc_label  f1 , enum division origin,INT_TYPE os, double * coeff, enum division alloy,INT_TYPE spin,  double tolerance , INT_TYPE maxRun , double power  );
double tCycleDecompostionGridOneMP ( INT_TYPE rank, struct sinc_label  f1 , enum division origin,INT_TYPE os, double * coeff, enum division alloy,INT_TYPE spin,  double tolerance , INT_TYPE maxRun , double power  );
double tCycleDecompostionListOneMP ( INT_TYPE rank, struct sinc_label  f1 , enum division origin,INT_TYPE os, double * coeff, enum division alloy,INT_TYPE spin,  double tolerance , INT_TYPE maxRun , double power  );
INT_TYPE tOuterProductSu( struct sinc_label  f1,enum division vector , INT_TYPE a, enum division vector2,INT_TYPE b, enum division proj, INT_TYPE c);
double tMultiplyMP(INT_TYPE rank,  INT_TYPE * info,struct sinc_label  f1,double scale, INT_TYPE beta,  enum division equals, INT_TYPE espin ,char leftChar, enum division left, INT_TYPE lspin, char rightChar,enum division right, INT_TYPE rspin);
double tMultiplyOne (INT_TYPE rank, struct sinc_label  f1,INT_TYPE dim,  enum division equals,INT_TYPE e, INT_TYPE espin , enum division left,INT_TYPE l,INT_TYPE lspin, enum division right,INT_TYPE r, INT_TYPE rspin);
double distance1(struct sinc_label  f1 ,enum division alloy ,INT_TYPE sp,  enum division alloyBak,INT_TYPE sp2);
double distance(struct sinc_label  f1 , enum division alloy , enum division alloyBak);
double inner(struct sinc_label  f1 , enum division alloy, INT_TYPE os );
double magnitude ( struct sinc_label  f1 , enum division alloy );
double pMagnitude ( struct sinc_label * f1 , enum division alloy );
INT_TYPE sortTerms (struct sinc_label  f1 , enum division term,INT_TYPE sp,enum division sorted,INT_TYPE sps );

void tHXpX (  INT_TYPE rank, struct sinc_label  f1 , enum division left,INT_TYPE shiftFlag, double sum,double product, double productCmpl, enum division equals ,  double tolerance , INT_TYPE maxRun,INT_TYPE solo  );
void pHXpX (  INT_TYPE rank, struct sinc_label  *f1 , enum division left,INT_TYPE shiftFlag, double sum,double product, double productCmpl, enum division equals ,  double tolerance , INT_TYPE maxRun,INT_TYPE solo  );

double positioningElectrons2 (INT_TYPE rank, struct sinc_label  f1 , enum division oneVector, enum division wavefunction,double x1, double y1, double z1, double x2, double y2, double z2);
INT_TYPE pReady ( struct sinc_label *f1 );
INT_TYPE ready ( struct sinc_label f1 );
INT_TYPE bootedQ ( struct sinc_label f1);

INT_TYPE xConstructFoundation (struct sinc_label calc , enum division usr, INT_TYPE UR, struct sinc_label calc2, enum division usz, INT_TYPE UZ ,INT_TYPE mx);
INT_TYPE printExpectationValues (struct sinc_label  f1 , enum division ha  , enum division vector);
INT_TYPE pPrintExpectationValues (struct sinc_label * f1 , enum division ha  , enum division vector);

INT_TYPE tOuterProductSuOne( struct sinc_label  f1,INT_TYPE space,enum division vector , INT_TYPE a, enum division vector2,INT_TYPE b, enum division proj, INT_TYPE c);
double canonicalGridDecompositionMP( INT_TYPE rank,struct sinc_label  f1 , Stream_Type * cofact, enum division origin,INT_TYPE l1,INT_TYPE l2,INT_TYPE os,   enum division alloy ,INT_TYPE l3,INT_TYPE l4,  INT_TYPE spin ,double tolerance,double magn, INT_TYPE preferred);
INT_TYPE tGEMV (INT_TYPE rank,  struct sinc_label  f1,INT_TYPE dim,  enum division equals,INT_TYPE e, INT_TYPE espin,enum division left,INT_TYPE l,INT_TYPE lspin, enum division right, INT_TYPE r,INT_TYPE rspin );
INT_TYPE tGEMM (INT_TYPE rank,  struct sinc_label  f1,INT_TYPE dim,  enum division equals, INT_TYPE e,INT_TYPE espin,enum division left,INT_TYPE l,INT_TYPE lspin, enum division right, INT_TYPE r,INT_TYPE rspin );
INT_TYPE tGEVV (INT_TYPE rank,  struct sinc_label  f1,INT_TYPE dim,  enum division equals, INT_TYPE espin,INT_TYPE leftChar, enum division left,INT_TYPE l,INT_TYPE lspin, INT_TYPE rightChar, enum division right, INT_TYPE r,INT_TYPE rspin );
double tDOT (INT_TYPE rank,  struct sinc_label  f1,INT_TYPE dim,char leftChar, enum division left,INT_TYPE l,INT_TYPE lspin, char rightChar, enum division right, INT_TYPE r,INT_TYPE rspin );
double traceOne( struct sinc_label  f1 , enum division label , INT_TYPE spin );

#endif /* mAls_h */
