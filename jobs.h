/*
 *  jobs.h
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

#ifndef jobs_h
#define jobs_h

#include "Model.h"
#include "input.h"
#include "eigen.h"
#ifndef APPLE
#include "ioPrint.h"
#endif

INT_TYPE loadDensity(struct calculation * c1 );

INT_TYPE foundation (struct calculation *c1 );
INT_TYPE krylov (struct calculation *c1 );
INT_TYPE ritz (struct calculation *c1 );
INT_TYPE decompose (struct calculation *c1 );

#endif /* jobs_h */
