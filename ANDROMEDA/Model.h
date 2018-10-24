/*
 *  Model.h
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
INT_TYPE iModel( struct calculation * c1);
INT_TYPE fModel( struct calculation * c1);
void resetExternal(struct calculation * i, INT_TYPE number, double scale );

#endif /* Model_h */
