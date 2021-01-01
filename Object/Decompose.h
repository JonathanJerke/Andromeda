/**
*  BoostFibonacci.h
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


#ifndef Decompose_h
#define Decompose_h
#include "constants.h"
#include "coreUtil.h"
#include "coreMath.h"
#include "interfaceMath.h"

floata canonicalRankDecomposition( sinc_label  f1 , floata * cofact,inta G,floata *GG, division origin,inta l1,inta l2,inta os, inta neo,division alloy ,inta l3 , inta l4,  inta spin ,double tolerance,double relativeTolerance, double condition,double maxCondition, inta maxCycle);

double CanonicalRankDecomposition ( sinc_label  f1 , double * cofact, division origin,inta os, division alloy,inta spin,  double tolerance ,  double relativeTolerance, double condition,double threshold, inta maxCycle ,double maxCondition, inta canon,inta X1 );

#endif /* Decompose_h */