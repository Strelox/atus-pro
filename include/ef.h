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
 * @file ef.h
 * @brief The file ef.h contains the Airy, Hermite polynome, and the eigenfunctions for the harmonic oscillator in 2d polar coordinates.
 * @author Želimir Marojević
 */

#ifndef _EF_H_
#define _EF_H_

#include "gsl/gsl_sf.h"
#include "math.h"

typedef double (*LINEF)(const double);

double AiryEF_1( const double x )
{
  return 1.4261 * gsl_sf_airy_Ai(x - 2.33811, 0);
}

double AiryEF_2( const double x )
{
  return 1.24516 * gsl_sf_airy_Ai(x - 4.08795, 0);
}

double AiryEF_3( const double x )
{
  return 1.1558 * gsl_sf_airy_Ai(x - 5.52056, 0);
}

double AiryEF_4( const double x )
{
  return 1.09787 * gsl_sf_airy_Ai(x - 6.78671, 0);
}

double AiryEF_5( const double x )
{
  return 1.05559 * gsl_sf_airy_Ai(x - 7.94413, 0);
}

double AiryEF_6( const double x )
{
  return 1.02258 * gsl_sf_airy_Ai(x - 9.02265, 0);
}

double AiryEF_7( const double x )
{
  return 0.995649 * gsl_sf_airy_Ai(x - 10.0402, 0);
}

double AiryEF_8( const double x )
{
  return 0.97301 * gsl_sf_airy_Ai(x - 11.0085, 0);
}

double AiryEF_9( const double x )
{
  return 0.953543 * gsl_sf_airy_Ai(x - 11.936, 0);
}

double AiryEF_10( const double x )
{
  return 0.93651 * gsl_sf_airy_Ai(x - 12.8288, 0);
}

LINEF AIRYEF[] = {&AiryEF_1, &AiryEF_2, &AiryEF_3, &AiryEF_4, &AiryEF_5, &AiryEF_6, &AiryEF_7, &AiryEF_8, &AiryEF_9, &AiryEF_10};

double R_polar_HO_0_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x);
  return retval;
}

double R_polar_HO_0_1( const double x )
{
  double retval = exp(-0.5000000000e0 * x) * sqrt(x);
  return retval;
}

double R_polar_HO_0_2( const double x )
{
  double retval = exp(-0.5000000000e0 * x) * x;
  return retval;
}

double R_polar_HO_0_3( const double x )
{
  double retval = exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1);
  return retval;
}

double R_polar_HO_0_4( const double x )
{
  double retval = exp(-0.5000000000e0 * x) * x * x;
  return retval;
}

double R_polar_HO_0_5( const double x )
{
  double retval = exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1);
  return retval;
}

double R_polar_HO_0_6( const double x )
{
  double retval = exp(-0.5000000000e0 * x) * pow(x, 0.3e1);
  return retval;
}

double R_polar_HO_0_7( const double x )
{
  double retval = exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1);
  return retval;
}

double R_polar_HO_0_8( const double x )
{
  double retval = exp(-0.5000000000e0 * x) * pow(x, 0.4e1);
  return retval;
}

double R_polar_HO_0_9( const double x )
{
  double retval = exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1);
  return retval;
}

double R_polar_HO_1_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x) - 0.1e1 * exp(-0.5000000000e0 * x) * x;
  return retval;
}

double R_polar_HO_1_1( const double x )
{
  double retval = 0.2e1 * exp(-0.5000000000e0 * x) * sqrt(x) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1);
  return retval;
}

double R_polar_HO_1_2( const double x )
{
  double retval = 0.3e1 * exp(-0.5000000000e0 * x) * x - 0.1e1 * exp(-0.5000000000e0 * x) * x * x;
  return retval;
}

double R_polar_HO_1_3( const double x )
{
  double retval = 0.4e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1);
  return retval;
}

double R_polar_HO_1_4( const double x )
{
  double retval = 0.5e1 * exp(-0.5000000000e0 * x) * x * x - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1);
  return retval;
}

double R_polar_HO_1_5( const double x )
{
  double retval = 0.6e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1);
  return retval;
}

double R_polar_HO_1_6( const double x )
{
  double retval = 0.7e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1);
  return retval;
}

double R_polar_HO_1_7( const double x )
{
  double retval = 0.8e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1);
  return retval;
}

double R_polar_HO_1_8( const double x )
{
  double retval = 0.9e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1);
  return retval;
}

double R_polar_HO_1_9( const double x )
{
  double retval = 0.10e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1);
  return retval;
}

double R_polar_HO_2_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x) - 0.2e1 * exp(-0.5000000000e0 * x) * x + 0.5000000000e0 * exp(-0.5000000000e0 * x) * x * x;
  return retval;
}

double R_polar_HO_2_1( const double x )
{
  double retval = 0.3e1 * exp(-0.5000000000e0 * x) * sqrt(x) - 0.3e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1);
  return retval;
}

double R_polar_HO_2_2( const double x )
{
  double retval = 0.6e1 * exp(-0.5000000000e0 * x) * x - 0.4e1 * exp(-0.5000000000e0 * x) * x * x + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1);
  return retval;
}

double R_polar_HO_2_3( const double x )
{
  double retval = 0.10e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) - 0.5e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1);
  return retval;
}

double R_polar_HO_2_4( const double x )
{
  double retval = 0.15e2 * exp(-0.5000000000e0 * x) * x * x - 0.6e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1);
  return retval;
}

double R_polar_HO_2_5( const double x )
{
  double retval = 0.21e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.7e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1);
  return retval;
}

double R_polar_HO_2_6( const double x )
{
  double retval = 0.28e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.8e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1);
  return retval;
}

double R_polar_HO_2_7( const double x )
{
  double retval = 0.36e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.9e1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1);
  return retval;
}

double R_polar_HO_2_8( const double x )
{
  double retval = 0.45e2 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.10e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1);
  return retval;
}

double R_polar_HO_2_9( const double x )
{
  double retval = 0.55e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.11e2 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1);
  return retval;
}

double R_polar_HO_3_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x) - 0.3e1 * exp(-0.5000000000e0 * x) * x + 0.1500000000e1 * exp(-0.5000000000e0 * x) * x * x - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1);
  return retval;
}

double R_polar_HO_3_1( const double x )
{
  double retval = 0.4e1 * exp(-0.5000000000e0 * x) * sqrt(x) - 0.6e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) + 0.2e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1);
  return retval;
}

double R_polar_HO_3_2( const double x )
{
  double retval = 0.10e2 * exp(-0.5000000000e0 * x) * x - 0.10e2 * exp(-0.5000000000e0 * x) * x * x + 0.2500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1);
  return retval;
}

double R_polar_HO_3_3( const double x )
{
  double retval = 0.20e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) - 0.15e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) + 0.3e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1);
  return retval;
}

double R_polar_HO_3_4( const double x )
{
  double retval = 0.35e2 * exp(-0.5000000000e0 * x) * x * x - 0.21e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.3500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1);
  return retval;
}

double R_polar_HO_3_5( const double x )
{
  double retval = 0.56e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.28e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.4e1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1);
  return retval;
}

double R_polar_HO_3_6( const double x )
{
  double retval = 0.84e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.36e2 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.4500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1);
  return retval;
}

double R_polar_HO_3_7( const double x )
{
  double retval = 0.120e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.45e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.5e1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1);
  return retval;
}

double R_polar_HO_3_8( const double x )
{
  double retval = 0.165e3 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.55e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.5500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1);
  return retval;
}

double R_polar_HO_3_9( const double x )
{
  double retval = 0.220e3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.66e2 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.6e1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1);
  return retval;
}

double R_polar_HO_4_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x) - 0.4e1 * exp(-0.5000000000e0 * x) * x + 0.3e1 * exp(-0.5000000000e0 * x) * x * x - 0.6666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1);
  return retval;
}

double R_polar_HO_4_1( const double x )
{
  double retval = 0.5e1 * exp(-0.5000000000e0 * x) * sqrt(x) - 0.10e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) + 0.5e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.8333333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1);
  return retval;
}

double R_polar_HO_4_2( const double x )
{
  double retval = 0.15e2 * exp(-0.5000000000e0 * x) * x - 0.20e2 * exp(-0.5000000000e0 * x) * x * x + 0.7500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1);
  return retval;
}

double R_polar_HO_4_3( const double x )
{
  double retval = 0.35e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) - 0.35e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) + 0.1050000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.1166666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1);
  return retval;
}

double R_polar_HO_4_4( const double x )
{
  double retval = 0.70e2 * exp(-0.5000000000e0 * x) * x * x - 0.56e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.14e2 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.1333333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1);
  return retval;
}

double R_polar_HO_4_5( const double x )
{
  double retval = 0.126e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.84e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.18e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.1500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1);
  return retval;
}

double R_polar_HO_4_6( const double x )
{
  double retval = 0.210e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.120e3 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.2250000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.1666666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1);
  return retval;
}

double R_polar_HO_4_7( const double x )
{
  double retval = 0.330e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.165e3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.2750000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.1833333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1);
  return retval;
}

double R_polar_HO_4_8( const double x )
{
  double retval = 0.495e3 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.220e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.33e2 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.2e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1);
  return retval;
}

double R_polar_HO_4_9( const double x )
{
  double retval = 0.715e3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.286e3 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.39e2 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.2166666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.4166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1);
  return retval;
}

double R_polar_HO_5_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x) - 0.5e1 * exp(-0.5000000000e0 * x) * x + 0.5e1 * exp(-0.5000000000e0 * x) * x * x - 0.1666666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.2083333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1);
  return retval;
}

double R_polar_HO_5_1( const double x )
{
  double retval = 0.6e1 * exp(-0.5000000000e0 * x) * sqrt(x) - 0.15e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) + 0.10e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.2500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.2500000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1);
  return retval;
}

double R_polar_HO_5_2( const double x )
{
  double retval = 0.21e2 * exp(-0.5000000000e0 * x) * x - 0.35e2 * exp(-0.5000000000e0 * x) * x * x + 0.1750000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.3500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.2916666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1);
  return retval;
}

double R_polar_HO_5_3( const double x )
{
  double retval = 0.56e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) - 0.70e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) + 0.28e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.4666666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.3333333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1);
  return retval;
}

double R_polar_HO_5_4( const double x )
{
  double retval = 0.126e3 * exp(-0.5000000000e0 * x) * x * x - 0.126e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.42e2 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.6e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.3750000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1);
  return retval;
}

double R_polar_HO_5_5( const double x )
{
  double retval = 0.252e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.210e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.60e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.7500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.4166666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1);
  return retval;
}

double R_polar_HO_5_6( const double x )
{
  double retval = 0.462e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.330e3 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.8250000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.9166666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.4583333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1);
  return retval;
}

double R_polar_HO_5_7( const double x )
{
  double retval = 0.792e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.495e3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.110e3 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.11e2 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.5000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1);
  return retval;
}

double R_polar_HO_5_8( const double x )
{
  double retval = 0.1287e4 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.715e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.143e3 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.13e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.5416666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1);
  return retval;
}

double R_polar_HO_5_9( const double x )
{
  double retval = 0.2002e4 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.1001e4 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.182e3 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.1516666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.5833333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) - 0.8333333333e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1);
  return retval;
}

double R_polar_HO_6_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x) - 0.6e1 * exp(-0.5000000000e0 * x) * x + 0.7500000000e1 * exp(-0.5000000000e0 * x) * x * x - 0.3333333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.6250000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.5000000000e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1);
  return retval;
}

double R_polar_HO_6_1( const double x )
{
  double retval = 0.7e1 * exp(-0.5000000000e0 * x) * sqrt(x) - 0.21e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) + 0.1750000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.5833333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.8750000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.5833333333e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1);
  return retval;
}

double R_polar_HO_6_2( const double x )
{
  double retval = 0.28e2 * exp(-0.5000000000e0 * x) * x - 0.56e2 * exp(-0.5000000000e0 * x) * x * x + 0.35e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.9333333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.1166666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.6666666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1);
  return retval;
}

double R_polar_HO_6_3( const double x )
{
  double retval = 0.84e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) - 0.126e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) + 0.63e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.14e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.1500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.7500000000e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1);
  return retval;
}

double R_polar_HO_6_4( const double x )
{
  double retval = 0.210e3 * exp(-0.5000000000e0 * x) * x * x - 0.252e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.105e3 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.20e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.1875000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.8333333333e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1);
  return retval;
}

double R_polar_HO_6_5( const double x )
{
  double retval = 0.462e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.462e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.165e3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.2750000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.2291666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.9166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1);
  return retval;
}

double R_polar_HO_6_6( const double x )
{
  double retval = 0.924e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.792e3 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.2475000000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.3666666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.2750000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) - 0.1000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1);
  return retval;
}

double R_polar_HO_6_7( const double x )
{
  double retval = 0.1716e4 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.1287e4 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.3575000000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.4766666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.3250000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) - 0.1083333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1);
  return retval;
}

double R_polar_HO_6_8( const double x )
{
  double retval = 0.3003e4 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.2002e4 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.5005000000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.6066666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.3791666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) - 0.1166666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2);
  return retval;
}

double R_polar_HO_6_9( const double x )
{
  double retval = 0.5005e4 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.3003e4 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.6825000000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.7583333333e2 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.4375000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) - 0.1250000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) + 0.1388888889e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1);
  return retval;
}

double R_polar_HO_7_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x) - 0.7e1 * exp(-0.5000000000e0 * x) * x + 0.1050000000e2 * exp(-0.5000000000e0 * x) * x * x - 0.5833333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.1458333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.1750000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.9722222222e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1);
  return retval;
}

double R_polar_HO_7_1( const double x )
{
  double retval = 0.8e1 * exp(-0.5000000000e0 * x) * sqrt(x) - 0.28e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) + 0.28e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.1166666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.2333333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.2333333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.1111111111e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1);
  return retval;
}

double R_polar_HO_7_2( const double x )
{
  double retval = 0.36e2 * exp(-0.5000000000e0 * x) * x - 0.84e2 * exp(-0.5000000000e0 * x) * x * x + 0.63e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.21e2 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.3500000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.3000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.1250000000e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1);
  return retval;
}

double R_polar_HO_7_3( const double x )
{
  double retval = 0.120e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) - 0.210e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) + 0.126e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.35e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.5e1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.3750000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.1388888889e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1);
  return retval;
}

double R_polar_HO_7_4( const double x )
{
  double retval = 0.330e3 * exp(-0.5000000000e0 * x) * x * x - 0.462e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.231e3 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.55e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.6875000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.4583333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.1527777778e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1);
  return retval;
}

double R_polar_HO_7_5( const double x )
{
  double retval = 0.792e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.924e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.396e3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.8250000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.9166666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.5500000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.1666666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1);
  return retval;
}

double R_polar_HO_7_6( const double x )
{
  double retval = 0.1716e4 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.1716e4 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.6435000000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.1191666667e3 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.1191666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) - 0.6500000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) + 0.1805555556e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2);
  return retval;
}

double R_polar_HO_7_7( const double x )
{
  double retval = 0.3432e4 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.3003e4 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.1001e4 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.1668333333e3 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.1516666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) - 0.7583333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) + 0.1944444444e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1);
  return retval;
}

double R_polar_HO_7_8( const double x )
{
  double retval = 0.6435e4 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.5005e4 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.1501500000e4 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.2275000000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.1895833333e2 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) - 0.8750000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) + 0.2083333333e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2);
  return retval;
}

double R_polar_HO_7_9( const double x )
{
  double retval = 0.11440e5 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.8008e4 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.2184e4 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.3033333333e3 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.2333333333e2 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) + 0.2222222222e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1) - 0.1984126984e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.23e2 / 0.2e1);
  return retval;
}

double R_polar_HO_8_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x) - 0.8e1 * exp(-0.5000000000e0 * x) * x + 0.14e2 * exp(-0.5000000000e0 * x) * x * x - 0.9333333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.2916666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.4666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.3888888889e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.1587301587e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1);
  return retval;
}

double R_polar_HO_8_1( const double x )
{
  double retval = 0.9e1 * exp(-0.5000000000e0 * x) * sqrt(x) - 0.36e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) + 0.42e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.21e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.5250000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.7000000000e0 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.5000000000e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.1785714286e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1);
  return retval;
}

double R_polar_HO_8_2( const double x )
{
  double retval = 0.45e2 * exp(-0.5000000000e0 * x) * x - 0.120e3 * exp(-0.5000000000e0 * x) * x * x + 0.105e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.42e2 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.8750000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.1e1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.6250000000e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) - 0.1984126984e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1);
  return retval;
}

double R_polar_HO_8_3( const double x )
{
  double retval = 0.165e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) - 0.330e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) + 0.231e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.77e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.1375000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.1375000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.7638888889e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) - 0.2182539683e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1);
  return retval;
}

double R_polar_HO_8_4( const double x )
{
  double retval = 0.495e3 * exp(-0.5000000000e0 * x) * x * x - 0.792e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.462e3 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.132e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.2062500000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.1833333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.9166666667e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) - 0.2380952381e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2);
  return retval;
}

double R_polar_HO_8_5( const double x )
{
  double retval = 0.1287e4 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.1716e4 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.858e3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.2145000000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.2979166667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.2383333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.1083333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) - 0.2579365079e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1);
  return retval;
}

double R_polar_HO_8_6( const double x )
{
  double retval = 0.3003e4 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.3432e4 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.1501500000e4 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.3336666667e3 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.4170833333e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) - 0.3033333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) + 0.1263888889e0 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) - 0.2777777778e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2);
  return retval;
}

double R_polar_HO_8_7( const double x )
{
  double retval = 0.6435e4 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.6435e4 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.2502500000e4 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.5005000000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.5687500000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) - 0.3791666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) + 0.1458333333e0 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) - 0.2976190476e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.23e2 / 0.2e1);
  return retval;
}

double R_polar_HO_8_8( const double x )
{
  double retval = 0.12870e5 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.11440e5 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.4004e4 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.728e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.7583333333e2 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) - 0.4666666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) + 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2) - 0.3174603175e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.12e2);
  return retval;
}

double R_polar_HO_8_9( const double x )
{
  double retval = 0.24310e5 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.19448e5 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.6188e4 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.1031333333e4 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.9916666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) - 0.5666666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) + 0.1888888889e0 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1) - 0.3373015873e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.23e2 / 0.2e1) + 0.2480158730e-4 * exp(-0.5000000000e0 * x) * pow(x, 0.25e2 / 0.2e1);
  return retval;
}

double R_polar_HO_9_0( const double x )
{
  double retval = exp(-0.5000000000e0 * x) - 0.9e1 * exp(-0.5000000000e0 * x) * x + 0.18e2 * exp(-0.5000000000e0 * x) * x * x - 0.14e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.5250000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.1050000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.1166666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.7142857143e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.2232142857e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1);
  return retval;
}

double R_polar_HO_9_1( const double x )
{
  double retval = 0.10e2 * exp(-0.5000000000e0 * x) * sqrt(x) - 0.45e2 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) + 0.60e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.35e2 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.1050000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.1750000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.1666666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.8928571429e-2 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.2480158730e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1);
  return retval;
}

double R_polar_HO_9_2( const double x )
{
  double retval = 0.55e2 * exp(-0.5000000000e0 * x) * x - 0.165e3 * exp(-0.5000000000e0 * x) * x * x + 0.165e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.77e2 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.1925000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.2750000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.2291666667e0 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) - 0.1091269841e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) + 0.2728174603e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2);
  return retval;
}

double R_polar_HO_9_3( const double x )
{
  double retval = 0.220e3 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1 / 0.2e1) - 0.495e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) + 0.396e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.154e3 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.33e2 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.4125000000e1 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.3055555556e0 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) - 0.1309523810e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) + 0.2976190476e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1);
  return retval;
}

double R_polar_HO_9_4( const double x )
{
  double retval = 0.715e3 * exp(-0.5000000000e0 * x) * x * x - 0.1287e4 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) + 0.858e3 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.286e3 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.5362500000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.5958333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.3972222222e0 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) - 0.1547619048e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) + 0.3224206349e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2);
  return retval;
}

double R_polar_HO_9_5( const double x )
{
  double retval = 0.2002e4 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1 / 0.2e1) - 0.3003e4 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) + 0.1716e4 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.5005000000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.8341666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.8341666667e1 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.5055555556e0 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) - 0.1805555556e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) + 0.3472222222e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.23e2 / 0.2e1);
  return retval;
}

double R_polar_HO_9_6( const double x )
{
  double retval = 0.5005e4 * exp(-0.5000000000e0 * x) * pow(x, 0.3e1) - 0.6435e4 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) + 0.3217500000e4 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) - 0.8341666667e3 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) + 0.1251250000e3 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) - 0.1137500000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) + 0.6319444444e0 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) - 0.2083333333e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2) + 0.3720238095e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.12e2);
  return retval;
}

double R_polar_HO_9_7( const double x )
{
  double retval = 0.11440e5 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1 / 0.2e1) - 0.12870e5 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) + 0.5720e4 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) - 0.1334666667e4 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) + 0.182e3 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) - 0.1516666667e2 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) + 0.7777777778e0 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) - 0.2380952381e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1) + 0.3968253968e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.23e2 / 0.2e1) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.25e2 / 0.2e1);
  return retval;
}

double R_polar_HO_9_8( const double x )
{
  double retval = 0.24310e5 * exp(-0.5000000000e0 * x) * pow(x, 0.4e1) - 0.24310e5 * exp(-0.5000000000e0 * x) * pow(x, 0.5e1) + 0.9724e4 * exp(-0.5000000000e0 * x) * pow(x, 0.6e1) - 0.2062666667e4 * exp(-0.5000000000e0 * x) * pow(x, 0.7e1) + 0.2578333333e3 * exp(-0.5000000000e0 * x) * pow(x, 0.8e1) - 0.1983333333e2 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1) + 0.9444444444e0 * exp(-0.5000000000e0 * x) * pow(x, 0.10e2) - 0.2698412698e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2) + 0.4216269841e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.12e2) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2);
  return retval;
}

double R_polar_HO_9_9( const double x )
{
  double retval = 0.48620e5 * exp(-0.5000000000e0 * x) * pow(x, 0.9e1 / 0.2e1) - 0.43758e5 * exp(-0.5000000000e0 * x) * pow(x, 0.11e2 / 0.2e1) + 0.15912e5 * exp(-0.5000000000e0 * x) * pow(x, 0.13e2 / 0.2e1) - 0.3094e4 * exp(-0.5000000000e0 * x) * pow(x, 0.15e2 / 0.2e1) + 0.357e3 * exp(-0.5000000000e0 * x) * pow(x, 0.17e2 / 0.2e1) - 0.2550000000e2 * exp(-0.5000000000e0 * x) * pow(x, 0.19e2 / 0.2e1) + 0.1133333333e1 * exp(-0.5000000000e0 * x) * pow(x, 0.21e2 / 0.2e1) - 0.3035714286e-1 * exp(-0.5000000000e0 * x) * pow(x, 0.23e2 / 0.2e1) + 0.4464285714e-3 * exp(-0.5000000000e0 * x) * pow(x, 0.25e2 / 0.2e1) - 0.2755731922e-5 * exp(-0.5000000000e0 * x) * pow(x, 0.27e2 / 0.2e1);
  return retval;
}

LINEF R_POLAR_HO[] = {&R_polar_HO_0_0, &R_polar_HO_0_1, &R_polar_HO_0_2, &R_polar_HO_0_3, &R_polar_HO_0_4, &R_polar_HO_0_5, &R_polar_HO_0_6, &R_polar_HO_0_7, &R_polar_HO_0_8, &R_polar_HO_0_9, &R_polar_HO_1_0, &R_polar_HO_1_1, &R_polar_HO_1_2, &R_polar_HO_1_3, &R_polar_HO_1_4, &R_polar_HO_1_5, &R_polar_HO_1_6, &R_polar_HO_1_7, &R_polar_HO_1_8, &R_polar_HO_1_9, &R_polar_HO_2_0, &R_polar_HO_2_1, &R_polar_HO_2_2, &R_polar_HO_2_3, &R_polar_HO_2_4, &R_polar_HO_2_5, &R_polar_HO_2_6, &R_polar_HO_2_7, &R_polar_HO_2_8, &R_polar_HO_2_9, &R_polar_HO_3_0, &R_polar_HO_3_1, &R_polar_HO_3_2, &R_polar_HO_3_3, &R_polar_HO_3_4, &R_polar_HO_3_5, &R_polar_HO_3_6, &R_polar_HO_3_7, &R_polar_HO_3_8, &R_polar_HO_3_9, &R_polar_HO_4_0, &R_polar_HO_4_1, &R_polar_HO_4_2, &R_polar_HO_4_3, &R_polar_HO_4_4, &R_polar_HO_4_5, &R_polar_HO_4_6, &R_polar_HO_4_7, &R_polar_HO_4_8, &R_polar_HO_4_9, &R_polar_HO_5_0, &R_polar_HO_5_1, &R_polar_HO_5_2, &R_polar_HO_5_3, &R_polar_HO_5_4, &R_polar_HO_5_5, &R_polar_HO_5_6, &R_polar_HO_5_7, &R_polar_HO_5_8, &R_polar_HO_5_9, &R_polar_HO_6_0, &R_polar_HO_6_1, &R_polar_HO_6_2, &R_polar_HO_6_3, &R_polar_HO_6_4, &R_polar_HO_6_5, &R_polar_HO_6_6, &R_polar_HO_6_7, &R_polar_HO_6_8, &R_polar_HO_6_9, &R_polar_HO_7_0, &R_polar_HO_7_1, &R_polar_HO_7_2, &R_polar_HO_7_3, &R_polar_HO_7_4, &R_polar_HO_7_5, &R_polar_HO_7_6, &R_polar_HO_7_7, &R_polar_HO_7_8, &R_polar_HO_7_9, &R_polar_HO_8_0, &R_polar_HO_8_1, &R_polar_HO_8_2, &R_polar_HO_8_3, &R_polar_HO_8_4, &R_polar_HO_8_5, &R_polar_HO_8_6, &R_polar_HO_8_7, &R_polar_HO_8_8, &R_polar_HO_8_9, &R_polar_HO_9_0, &R_polar_HO_9_1, &R_polar_HO_9_2, &R_polar_HO_9_3, &R_polar_HO_9_4, &R_polar_HO_9_5, &R_polar_HO_9_6, &R_polar_HO_9_7, &R_polar_HO_9_8, &R_polar_HO_9_9};


double HermiteH_0( const double x )
{
  return 1.0;
}

double HermiteH_1( const double x )
{
  double retval = 2 * x;
  return retval;
}

double HermiteH_2( const double x )
{
  double retval = -2 + 4 * x * x;
  return retval;
}

double HermiteH_3( const double x )
{
  double retval = 0.8e1 * pow(x, 0.3e1) - 0.12e2 * x;
  return retval;
}

double HermiteH_4( const double x )
{
  double retval = 0.12e2 + 0.16e2 * pow(x, 0.4e1) - 0.48e2 * x * x;
  return retval;
}

double HermiteH_5( const double x )
{
  double retval = 0.32e2 * pow(x, 0.5e1) - 0.160e3 * pow(x, 0.3e1) + 0.120e3 * x;
  return retval;
}

double HermiteH_6( const double x )
{
  double retval = -0.120e3 + 0.64e2 * pow(x, 0.6e1) - 0.480e3 * pow(x, 0.4e1) + 0.720e3 * x * x;
  return retval;
}

double HermiteH_7( const double x )
{
  double retval = 0.128e3 * pow(x, 0.7e1) - 0.1344e4 * pow(x, 0.5e1) + 0.3360e4 * pow(x, 0.3e1) - 0.1680e4 * x;
  return retval;
}

double HermiteH_8( const double x )
{
  double retval = 0.1680e4 + 0.256e3 * pow(x, 0.8e1) - 0.3584e4 * pow(x, 0.6e1) + 0.13440e5 * pow(x, 0.4e1) - 0.13440e5 * x * x;
  return retval;
}

double HermiteH_9( const double x )
{
  double retval = 0.512e3 * pow(x, 0.9e1) - 0.9216e4 * pow(x, 0.7e1) + 0.48384e5 * pow(x, 0.5e1) - 0.80640e5 * pow(x, 0.3e1) + 0.30240e5 * x;
  return retval;
}

LINEF HERMITEH[] = {&HermiteH_0, &HermiteH_1, &HermiteH_2, &HermiteH_3, &HermiteH_4, &HermiteH_5, &HermiteH_6, &HermiteH_7, &HermiteH_8, &HermiteH_9 };

#endif