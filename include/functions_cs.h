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
 * @file functions_cs.h
 * @brief The file functions_cs.h contains classes for eigenfunctions and real potentials.
 * @author Želimir Marojević
 */

#ifndef _MY_FUNCTIONS
#define _MY_FUNCTIONS

#include <deal.II/base/function.h>
#include "ef.h"

using namespace std;
using namespace dealii;

namespace BreedSolver_cs
{
  class myexception: public exception
  {
    virtual const char *what() const throw()
    {
      return "Quantum number out of range.";
    }
  } myex;

  template <int dim>
  class CEigenfunctions : public Function<dim>
  {
  public:
    CEigenfunctions ( unsigned QN[3], const std::vector<double> &a ) : Function<dim>()
    {
      if ( QN[0] < 0 || QN[0] > 9 ) throw myex;
      if ( QN[1] < 0 || QN[1] > 9 ) throw myex;
      if ( QN[2] < 0 || QN[2] > 9 ) throw myex;
      m_QNx = QN[0];
      m_QNy = QN[1];
      m_QNz = QN[2];
      m_fakx = a[0];
      m_faky = a[1];
      m_fakz = a[2];
    }
    virtual double value ( const Point<dim> &p, const unsigned int  component = 0) const;

  private:
    unsigned int m_QNx; /**< quantum number of the gravitational eigenfunctions */
    unsigned int m_QNy; /**< quantum number of the eigenfunctions for the harmonic trap in polar coordinates */
    unsigned int m_QNz; /**< angular momentum quantum number */
    double m_fakx; /**< gravitational acceleration \f$ \beta \f$ */
    double m_faky; /**< harmonic trapping frequency \f$ \omega \f$ */
    double m_fakz; /**< not used by _cs programs */
  };

  template <int dim>
  double CEigenfunctions<dim>::value( const Point<dim> &p, const unsigned int component ) const
  {
    double retval = (*AIRYEF[m_QNx])(pow(m_fakx, 1.0 / 3.0) * p(0)) * (*R_POLAR_HO[10 * m_QNy + m_QNz])(m_faky * m_faky * p(1) * p(1));
    //double retval = (*AIRYEF[m_QNx])(pow(m_fakx,1.0/3.0)*p(0)) * (*HERMITEH[m_QNy])(sqrt(m_faky)*p(1))*exp(-0.5*m_faky*p(1)*p(1));
    return retval;
  }

  /*************************************************************************************************/

  template <int dim>
  class CPotential : public Function<dim>
  {
  public:
    CPotential ( const std::vector<double> &a, const long l  ) : Function<dim>()
    {
      m_fakx = a[0];
      m_faky = a[1];
      m_fakz = a[2];
      m_m = double(l * l);
    }
    virtual double value ( const Point<dim> &p, const unsigned int  component = 0) const;

    double m_fakx; /**< gravitational acceleration \f$ \beta \f$ */
    double m_faky; /**< harmonic trapping frequency \f$ \omega \f$ */
    double m_fakz; /**< not used by _cs programs */
    double m_m; /**< angular momentum quantum number */
  };

  /** Computes the value of the potential.
   * It is not necessary to invoke this function directly.
   * @param[in] p a spatial point (p(0) ≡ z, p(1) ≡ r)
   * @return \f$ V(\rho,z) \f$
   */
  template <int dim>
  double CPotential<dim>::value( const Point<dim> &p, const unsigned int component ) const
  {
    double retval = 0;
    double rq = p(1) * p(1);
    if ( m_m == 0 )
    {
      retval = m_fakx * p(0) + m_faky * m_faky * rq;
    }
    else
    {
      if ( rq > 0 )
        retval = m_fakx * p(0) + m_faky * m_faky * rq + m_m / rq;
      else
        retval = 0;
    }
    return retval;
  }
}
#endif