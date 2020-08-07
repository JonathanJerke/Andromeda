/**
 *  mAls.h
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

#ifndef mAls_h
#define mAls_h
#include "constants.h"
#include "coreUtil.h"
#include "coreMath.h"
#include "interfaceMath.h"
#include "saUtil.h"
#include "eigen.h"

inta canonicalRankDecomposition( sinc_label  f1 , floata * cofact,floata *GG, division origin,inta l1,inta l2,inta os, division alloy ,inta l3 , inta l4,  inta spin ,double tolerance,double relativeTolerance, double condition,double threshold, inta maxCycle);
double AsterCanonicalRankDecomposition ( inta rank,  sinc_label  f1 , double * cofact, division origin,inta os, division alloy,inta spin,  double tolerance ,  double relativeTolerance, double condition,double threshold, inta maxCycle , inta canon );
double printExpectationValues ( calculation *c,   sinc_label  f1 , division ha  , division vector);
double tMatrixElements ( inta rank,  sinc_label  f1 ,   division bra, inta bspin,  division mat, inta mspin,  division ket, inta kspin );
inta tOuterProductSu( sinc_label  f1,  division vector , inta a, division vector2,inta b, division proj, inta c);
double magnitude ( sinc_label  f1 , division alloy, inta spin);
inta tOuterProductSuOne( sinc_label  f1,inta space, division vector , inta a,   division vector2,inta b,   division proj, inta c);
inta tGEMV (inta rank, sinc_label  f1,inta dim, division equals,inta e, inta espin, division left,inta l,inta lspin, division right, inta r,inta rspin );
inta tGEVV (inta rank,  sinc_label  f1,inta dim, division equals,inta e, inta espin, division left,inta l,inta lspin, division right, inta r,inta rspin );
double tDOT (inta rank, sinc_label  f1,inta dim,char leftChar, division left,inta l,inta lspin, char rightChar, division right, inta r,inta rspin );
inta tHX(  inta rank, sinc_label f1 ,division left, inta l, inta im, double prod, division ket , inta k, inta sp2,   division oket, inta o,inta ospin );
void tHXpY (  inta rank,   sinc_label f1 ,  division bra,   division left,inta shiftFlag,  division right ,  double tolerance ,double relativeTolerance,double condition,double threshold, inta maxCycle, inta canon,inta X1);
#endif /* mAls_h */
