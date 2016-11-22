/*
    This file is part of atus-pro package.

    atus-pro is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atus-pro is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with atus-pro.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
    Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
    (ZARM - Center of Applied Space Technology and Microgravity, Germany, http://www.zarm.uni-bremen.de/)

    Public use and modification of this code are allowed provided that the following papers are cited:

    Marojević, Želimir, Ertan Göklü, und Claus Lämmerzahl. "ATUS-PRO: A FEM-based solver for the time-dependent and stationary Gross–Pitaevskii equation",
    Computer Physics Communications, Vol 202, 2016, p. 216--232. doi:10.1016/j.cpc.2015.12.004.

    W. Bangerth and D. Davydov and T. Heister and L. Heltai and G. Kanschat and M. Kronbichler and M. Maier and B. Turcksin and D. Wells.
    "The \texttt{deal.II} Library, Version 8.4", Journal of Numerical Mathematics, vol 24, 2016.

    The authors would be grateful for all information and/or comments regarding the use of the code.
*/


/**
* @file gen_params.cpp
* @author Želimir Marojević
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <fstream>

#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

const long N1[] = {4, 4, 4};
const long NA = 10000;
const long Ndmu = 40;
const double dmu = .1;
double omega[] = {0.5, 0.5, 0.5};

int main( int argc, char *argv[] )
{
  //const char * HomePath = getenv ("HOME");

  string tmpstr;

  char base_folder[255];
  char folder[255];
  char filename[255];
  char shellcmd[255];
  char datetime[255];

  time_t rawtime;
  struct tm *timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(datetime, 255, "%F_%R", timeinfo);

#if POTENTIAL==1
  sprintf( base_folder, "HTRAP_%s", datetime );
#endif
#if POTENTIAL==2
  sprintf( base_folder, "GOST_%s", datetime );
  omega[0] = 0.5039287608;
#endif

#if DIMENSION==2
  for ( long i1 = 0; i1 < N1[0]; i1++ )
  {
    for ( long j1 = 0; j1 < N1[1]; j1++ )
    {
      sprintf( folder, "%ld_%ld", i1, j1 );
      sprintf( shellcmd, "mkdir -p %s/%s", base_folder, folder );
      system( shellcmd );

      sprintf( filename, "%s/%s/params.xml", base_folder, folder );
      FILE *fh = fopen( filename, "w" );

      XMLPrinter printer(fh);
      printer.OpenElement( "PARAMETER" );
      printer.OpenElement( "filename" );
      printer.PushText("final.bin");
      printer.CloseElement();
      printer.OpenElement( "guess_fct" );
      printer.PushText( "exp(-x*x-y*y)" );
      printer.CloseElement();

      // problem related parameter
      printer.OpenElement( "PHYSICS" );
      printer.OpenElement( "QN1" );
      tmpstr = to_string(i1) + "," + to_string(j1);
      printer.PushText(tmpstr.c_str());
      printer.CloseElement();
      printer.OpenElement( "omega" );
      tmpstr = to_string(omega[0]) + "," + to_string(omega[1]) + "," + to_string(omega[2]);
      printer.PushText(tmpstr.c_str());
      printer.CloseElement();
      printer.OpenElement( "gs_1" );
      printer.PushText("1,1");
      printer.CloseElement();
      printer.OpenElement( "mu" );
      printer.PushText("5,5");
      printer.CloseElement();
      printer.CloseElement(); // close PHYSICS

      // mesh related parameter
      printer.OpenElement( "MESH" );
      printer.OpenElement( "xrange" );
      printer.PushText("0,30");
      printer.CloseElement();
      printer.OpenElement( "yrange" );
      printer.PushText("-15,15");
      printer.CloseElement();
      printer.OpenElement( "zrange" );
      printer.PushText("-15,15");
      printer.CloseElement();
      printer.OpenElement( "global_refinements" );
      printer.PushText( "10" );
      printer.CloseElement();
      printer.CloseElement(); // close MESH

      // algorithm related parameter
      printer.OpenElement( "ALGORITHM" );
      printer.OpenElement( "ti" );
      printer.PushText( "1" );
      printer.CloseElement();
      printer.OpenElement( "NA" );
      tmpstr = to_string(NA);
      printer.PushText( tmpstr.c_str() );
      printer.CloseElement();
      printer.OpenElement( "NK" );
      printer.PushText( "10" );
      printer.CloseElement();
      printer.OpenElement( "dmu" );
      tmpstr = to_string(dmu);
      printer.PushText( tmpstr.c_str() );
      printer.CloseElement();
      printer.OpenElement( "Ndmu" );
      tmpstr = to_string(Ndmu);
      printer.PushText( tmpstr.c_str() );
      printer.CloseElement();
      printer.OpenElement( "dt" );
      printer.PushText( "0.001" );
      printer.CloseElement();
      printer.OpenElement( "epsilon" );
      printer.PushText( "1e-5,1e-10" );
      printer.CloseElement();
      printer.OpenElement( "df" );
      printer.PushText( "1" );
      printer.CloseElement();
      printer.CloseElement(); // close ALGORITHM
      printer.CloseElement(); // close PARAMETER
    }
  }
#endif
  return 0;
}
