/**
 *  input.h
 *
 *
 *  Copyright 2021 Jonathan Jerke and Bill Poirier.
 *  Ongoing support for this program is coordinated through quantumgalaxies.org.
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


#ifndef input_h
#define input_h

#include"constants.h"
#include"coreMath.h"

inta resetA(   runTime *rt);
inta allowQ(   runTime *  f1,   blockMemoryType a );
inta blockA(   runTime *  f1,   blockMemoryType a );
inta rotateGeometry (  calculation * c, double * u  );
inta getInitialGeneral(  calculation * c, const char * input_line );
inta initCalculation(  calculation * c );
inta destructCalculation (   calculation * c );
inta control ( const char * line );
inta comment ( const char * line );
inta getCore(  calculation * c, const char * input_line );
inta readInput(  calculation *c ,  field *f1, FILE * in);
inta getParam (   calculation * c,  input_label* f1, const char * input_line );
inta getGeometry(  calculation * c,const char * input_line );
inta getInputOutput(  calculation * c,  field *f1,const char * input_line );
inta finalizeInit(  calculation * c );
inta modifyGeometry(  calculation * c,const char * input_line );
inta getBuildParameters(  calculation * c,const char * input_line );
inta bootShell (inta argc , char * argv[], calculation * c1, field *f);
inta readShell (inta argc , char * argv[], calculation * c1, field *f);

#endif /* input_h */
