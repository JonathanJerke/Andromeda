/**
 *  jobs.h
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

#ifndef jobs_h
#define jobs_h

#include "Model.h"
#include "input.h"
#include "eigen.h"
#ifndef APPLE
#include "ioPrint.h"
#endif

inta foundationB ( calculation *c1  ,  field f1);
inta foundationS ( calculation *c1  ,  field f1);
double singlekrylov ( calculation *c1  ,  field f1);
inta ritz ( calculation *c1  ,  field f1);
int run (inta argc , char * argv[]);

#endif /* jobs_h */
