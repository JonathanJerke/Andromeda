/**
 *  coreForce.h
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

#ifndef coreForce_h
#define coreForce_h

#include "constants.h"
#include "coreMath.h"
#include "coreUtil.h"
#include "ioPrint.h"

void getDescription (   function_label *fn ,double scalar,FILE * outString);
void getMetric (   metric_label mu ,FILE * outString);

mea BoB (  basisElement_label b1,   basisElement_label b2);
inta buildMetric( double ,inta Z,   function_label func,inta am,   metric_label * metric);
inta buildPairWisePotential(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,  blockType bl, division pair,   inta particle1 ,inta embed, inta overline,  spinType cmpl,   metric_label mu);
inta buildExternalPotential(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,  blockType bl, division single,   inta particle1,inta embed, inta overline,  spinType cmpl,   metric_label mu,inta a);
inta buildConstant(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl);
inta buildDeriv(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl);
inta buildKinetic(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl);
inta buildClampKinetic(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,  inta particle1, inta overline,   spinType cmpl);
inta buildSpring(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl);
inta buildLinear(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl);
inta buildElement(  calculation *c1,   sinc_label *f1,double scalar,inta invert,inta act,   blockType bl,  division single,inta label, inta overline,   spinType cmpl,inta bra, inta ket);
inta assignDiagonalMatrix(calculation *c1,field *f, char * filename, division single );
#endif /* coreForce_h */
