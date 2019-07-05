/*
 *  Model.h
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

#ifndef Model_h
#define Model_h

#include "constants.h"
#include "coreUtil.h"
#include "coreMath.h"
#include "interfaceMath.h"
#include "mAls.h"
#include "coreForce.h"
#include "main.h"
struct calculation initCal (void);
struct calculation gas (void );
INT_TYPE iModel( struct calculation * c1, struct field *f1);
INT_TYPE fModel( struct sinc_label *f1);
void resetExternal(struct calculation * i, INT_TYPE number, double scale );
struct field initField (void );
INT_TYPE singleModel( struct calculation * c1, struct field *f);
INT_TYPE multModel( struct calculation * c1,INT_TYPE nv, struct field * v,struct field *f);
#endif /* Model_h */
