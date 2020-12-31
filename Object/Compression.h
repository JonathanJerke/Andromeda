/**
*  Compression.h
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

#ifndef Compression_h
#define Compression_h
#include "constants.h"
#include "coreUtil.h"
#include "coreMath.h"
#include "interfaceMath.h"



floata canonicalRankCompression( inta ** spatial, floata * cofact,sinc_label f1, inta G,floata *GG, division origin,inta l1,inta l2,inta os, sinc_label  f2 ,inta neo,division alloy ,inta l3 , inta l4, inta spin ,double tolerance,double relativeTolerance, double condition,double maxCondition, inta maxCycle);

double CanonicalRankCompression ( inta ** spatial, double * cofact, sinc_label f1, division origin,inta os, sinc_label  f2 , division alloy,inta spin, double tolerance ,  double relativeTolerance, double condition,double threshold, inta maxCycle ,double maxCondition, inta canon,inta X1 );

#endif /* Compression_h */
