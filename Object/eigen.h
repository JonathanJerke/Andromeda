/**
 *  eigen.h
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


#ifndef eigen_h
#define eigen_h
#include "constants.h"
#include "coreUtil.h"
#include "saUtil.h"
#include "interfaceMath.h"

inta tSelect(  sinc_label  f1, inta Ve, inta type,   division usr,  inta testFlag);
inta tFilter(  sinc_label  f1, inta Ve, inta type,   division usr);
inta tBuildMatrix (inta minusFlag,   sinc_label  f1,   division A ,   division usz, inta quantumBasisSize);
inta tSolveMatrix (inta typer,   sinc_label  f1,inta Ne,  division usz, inta quantumBasisSize,   division outputValues);
floata tComponent( sinc_label f1, division hamiltonian, inta space, inta index);
floata tComponentPoint(calculation * c, sinc_label f1, division hamiltonian, inta space, inta index);
#endif /* eigen_h */
