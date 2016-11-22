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
 * @file functions.h
 * @brief The file functions.h contains classes for real eigenfunctions and real potentials.
 * @author Želimir Marojević
 */

#ifndef _MY_FUNCTIONS
#define _MY_FUNCTIONS

#include <deal.II/base/function.h>
#include <exception>
#include "ef.h"

using namespace std;
using namespace dealii;

namespace BreedSolver
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
    virtual double value ( const Point<dim> &p, const unsigned component = 0) const;

  private:
    unsigned m_QNx; /**< quantum number of the eigenfunctions of the harmonic trap or the gravitational trap */
    unsigned m_QNy; /**< quantum number of the eigenfunctions of the harmonic trap for the \f$ y \f$ coordinate */
    unsigned m_QNz; /**< quantum number of the eigenfunctions of the harmonic trap for the \f$ z \f$ coordinate */
    double m_fakx; /**< \f$ \omega_x \f$ or \f$ \beta \f$ */
    double m_faky; /**< \f$ \omega_y \f$ */
    double m_fakz; /**< \f$ \omega_z \f$ */
  };

  template <int dim>
  double CEigenfunctions<dim>::value( const Point<dim> &p, const unsigned component ) const
  {
    double retval;

#if DIMENSION==2
#if POTENTIAL==1
    retval = (*HERMITEH[m_QNx])(sqrt(m_fakx) * p(0)) * exp(-0.5 * m_fakx * p(0) * p(0)) * (*HERMITEH[m_QNy])(sqrt(m_faky) * p(1)) * exp(-0.5 * m_faky * p(1) * p(1));
#endif
#if POTENTIAL==2
    retval = (*AIRYEF[m_QNx])(pow(m_fakx, 1.0 / 3.0) * p(0)) * (*HERMITEH[m_QNy])(sqrt(m_faky) * p(1)) * exp(-0.5 * m_faky * p(1) * p(1));
#endif
#endif

#if DIMENSION==3
#if POTENTIAL==1
    retval = (*HERMITEH[m_QNx])(sqrt(m_fakx) * p(0)) * exp(-0.5 * m_fakx * p(0) * p(0)) * (*HERMITEH[m_QNy])(sqrt(m_faky) * p(1)) * exp(-0.5 * m_faky * p(1) * p(1)) * (*HERMITEH[m_QNz])(sqrt(m_fakz) * p(2)) * exp(-0.5 * m_fakz * p(2) * p(2));
#endif
#if POTENTIAL==2
    retval = (*AIRYEF[m_QNx])(pow(m_fakx, 1.0 / 3.0) * p(0)) * (*HERMITEH[m_QNy])(sqrt(m_faky) * p(1)) * exp(-0.5 * m_faky * p(1) * p(1)) * (*HERMITEH[m_QNz])(sqrt(m_fakz) * p(2)) * exp(-0.5 * m_fakz * p(2) * p(2));
#endif
#endif
    return retval;
  }

  /*************************************************************************************************/

  template <int dim>
  class CPotential : public Function<dim>
  {
  public:
    CPotential ( const std::vector<double> &a ) : Function<dim>()
    {
      m_fakx = a[0];
      m_faky = a[1];
      m_fakz = a[2];
    }
    virtual double value ( const Point<dim> &p, const unsigned component = 0) const;

    double m_fakx; /**< \f$ \omega_x \f$ or \f$ \beta \f$ */
    double m_faky; /**< \f$ \omega_y \f$ */
    double m_fakz; /**< \f$ \omega_z \f$ */
  };

  /** Computes the value of the potential.
   * It is not necessary to invoke this function directly.
   * @param[in] p a spatial point (p(0) ≡ x, p(1) ≡ y, p(2) ≡ z)
   * @return \f$ V(x,y) \f$ or \f$ V(x,y,z) \f$
   */
  template <int dim>
  double CPotential<dim>::value( const Point<dim> &p, const unsigned component ) const
  {
    double retval;

#if DIMENSION==2
#if POTENTIAL==1
    retval = m_fakx * m_fakx * p(0) * p(0) + m_faky * m_faky * p(1) * p(1);
#endif
#if POTENTIAL==2
    retval = m_fakx * p(0) + m_faky * m_faky * p(1) * p(1);
#endif
#endif

#if DIMENSION==3
#if POTENTIAL==1
    retval = m_fakx * m_fakx * p(0) * p(0) + m_faky * m_faky * p(1) * p(1) + m_fakz * m_fakz * p(2) * p(2);
#endif
#if POTENTIAL==2
    retval = m_fakx * p(0) + m_faky * m_faky * p(1) * p(1) + m_fakz * m_fakz * p(2) * p(2);
#endif
#endif
    return retval;
  }
}
#endif