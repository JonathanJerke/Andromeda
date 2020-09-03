/**
 *  saUtil.h
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


#ifndef saUtil_h
#define saUtil_h

#include "constants.h"
#include "coreUtil.h"
#include "mAls.h"
#define sr7 0.3779644730092272272
#define sr6 0.4082482904638630164
#define sr5 0.4472135954999579393
#define hf  0.5000000000000000000
#define sr3 0.5773502691896257645
#define sr2 0.7071067811865475244

void testSAAgain (   sinc_label f1 ,   division vector );
inta tPerms(  bodyType bd);
void tTestSA (  bodyType bd, inta n);
inta tSize(  bodyType bd);
inta tPaths(  bodyType bd , inta irrep );
inta nEqua(  bodyType bd, inta *a );
void tIntr (   bodyType bd , inta eq , double * a);
inta tInnerTest(   sinc_label f1,   division A ,  division B);
inta tSA (  bodyType bd, inta X, inta Y, inta Z, inta T );
double tGetProjection(   bodyType bd, inta type , inta perm );
double tGetVector(  bodyType bd , inta type , inta perm );
inta tClassifyComponents(   sinc_label  , double * up, double * entropy );
inta tClassify(inta rank,   sinc_label  ,   division label);
inta tBuildIrr ( inta rank,   sinc_label , inta meta ,   division origin, inta ospin,   division targ , inta tspin);
inta tPermuteOne(inta rank,   sinc_label , inta dim, inta leftChar ,   division left, inta l, inta lspin,   division equals,inta e, inta espin);
inta tPermute(inta rank,   sinc_label , inta leftChar ,   division left, inta lspin,   division equals, inta espin);
inta tAllCompPermMultiplyMP( inta rank,   sinc_label  f1 ,   division left ,inta lspin,   division right ,inta rspin, double * sequ);
  bodyType commandSA(  bodyType bd, inta act,   blockType cl,   blockType bl,inta perm[], inta op[]);
inta tAddUpComponents( inta rank,   sinc_label  f1 ,   division left ,   division right ,  double *up);
inta tTabulateInnerProjection( inta rank,   sinc_label  f1 ,   division vec, double *up);
double tTestInner (   bodyType bd, inta i , inta j );
inta tTest (   bodyType bd );
inta tSizeUp(inta rank,   sinc_label  f1 , inta type,   division label);
inta testSA (   sinc_label f1 ,   division vector );
inta irreps (   bodyType bd, inta type );
inta mapir(  bodyType bd , inta class );
#endif /* saUtil_h */
