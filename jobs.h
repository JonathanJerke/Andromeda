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

INT_TYPE sumTo2(struct sinc_label f1, INT_TYPE nn,enum division mat,INT_TYPE ms, enum division sum,INT_TYPE spin);
INT_TYPE sumTo3(struct sinc_label f1, enum division mat,INT_TYPE ms, enum division sum,INT_TYPE spin);
INT_TYPE sumTo4(struct sinc_label f1, enum division mat,INT_TYPE ms, enum division sum,INT_TYPE spin);

INT_TYPE countHam ( struct calculation *c1 , struct field f1 );
INT_TYPE foundation1 (struct calculation *c1 , struct field f1);
INT_TYPE foundationM (struct calculation *c1 , struct field f1);
INT_TYPE decompose ( struct calculation *c1, struct field f1);
INT_TYPE krylov (struct calculation *c1  , struct field f1);
INT_TYPE ritz (struct calculation *c1  , struct field f1);
INT_TYPE frameEwald( struct calculation * c , struct field f);

#endif /* jobs_h */
