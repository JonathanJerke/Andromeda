/*
 *  input.h
 *
 *
 *  Copyright 2019 Jonathan Jerke and Bill Poirier.
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


#ifndef input_h
#define input_h

#include"constants.h"
#include"coreMath.h"
#include <dirent.h>
INT_TYPE resetA( struct runTime *rt);
INT_TYPE allowQ( struct runTime *  f1, enum blockMemoryType a );
INT_TYPE blockA( struct runTime *  f1, enum blockMemoryType a );
INT_TYPE rotateGeometry (struct calculation * c, double * u  );
INT_TYPE getInitialGeneral(struct calculation * c, const char * input_line );
INT_TYPE initCalculation(struct calculation * c );
INT_TYPE destructCalculation ( struct calculation * c );
INT_TYPE control ( const char * line );
INT_TYPE comment ( const char * line );
INT_TYPE getCore(struct calculation * c, const char * input_line );
INT_TYPE readInput(struct calculation *c ,struct input_label *f1, FILE * in);
INT_TYPE setOutputSchedule( struct calculation * c, const char * segement );
INT_TYPE  getParam ( struct calculation * c,struct input_label* f1, const char * input_line );
INT_TYPE  getGeometry(struct calculation * c,const char * input_line );
INT_TYPE modGeometry(struct calculation * c,const char * input_line );
INT_TYPE getBasis(struct calculation * c,  char * input_line );
INT_TYPE getInputOutput(struct calculation * c,struct input_label *f1,const char * input_line );
INT_TYPE finalizeInit(struct calculation * c );
INT_TYPE modifyGeometry(struct calculation * c,const char * input_line );
INT_TYPE getBuildParameters(struct calculation * c,const char * input_line );
INT_TYPE defineNaturalBasis ( struct calculation * c , char * input_line );
INT_TYPE assignNaturalBasis ( struct calculation * c, char * input_line );
INT_TYPE bootShell (INT_TYPE argc , char * argv[],struct calculation * c1, struct field *f);
INT_TYPE estSize ( INT_TYPE interval );
#endif /* input_h */
