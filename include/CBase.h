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


#ifndef __CBASE_H__
#define __CBASE_H__

#include <mpi.h>
#include <cstdint>
#include "ref_pt_list.h"

enum Status { SUCCESS, FAILED, ZERO_SOL, CONTINUE };

template <int dim>
class CBase
{
public:
  CBase( const string & );
  virtual ~CBase() {};

  int find_ortho_min( bool = true );

  const double l2norm_t() const;

  virtual void compute_contributions() = 0;

  MPI_Comm mpi_communicator;
protected:
  void screening();

  ///
  double m_t[dim];
  /// initial guess for the initial position if function space
  double m_t_guess[dim];

  double m_xmin, m_xmax;
  double m_ymin, m_ymax;
  double m_zmin, m_zmax;

  /// current residual of the rhs
  double m_res;
  /// residual of the previous step
  double m_res_old;
  double m_resp;
  double m_res_over_resp;
  /// initial value for the initial position in function space
  double m_ti;

  /// damping factor of the newton step
  double m_df;

  std::vector<double> m_epsilon;

  /// bifurcation parameter \f$ \mu \f$
  double m_mu;
  /// \f$ \Delta\mu \f$
  double m_dmu;
  /// non linearity parameter \f$ \gamma \f$
  std::vector<double> m_gs;
  /// parameters for the potential
  std::vector<double> m_omega;
  /// particle number
  double m_N;
  /// final L_2 norm of the L_2 gradient
  double m_final_error;

  /// counter for the inner iteration loop
  unsigned m_counter;
  /// initial global refinement of the grid
  unsigned m_global_refinement;
  /// total number of cells including ghost cells and artificial cells
  unsigned m_total_no_cells;
  /// total number of active cells
  unsigned m_total_no_active_cells;
  /// determines the frequency of intermediate outputs
  unsigned m_NA;
  /// number of \f$\Delta\mu\f$ increases
  unsigned m_Ndmu;
  /// specifies the quantum numbers for the initial guess
  unsigned m_QN1[3];

  /// true if m_myrank is zero
  bool m_root;
  /// MPI rank
  int m_rank;

  /// output stream for benchmark statistics
  ofstream m_computing_timer_log;
  /// object which tracks run times of different functions
  TimerOutput m_computing_timer;
  /// reference to the parameter handler
  MyParameterHandler m_ph;
  /// output stream which redirects only root outputs to std::out
  ConditionalOStream pcout;

  MyUtils::ref_pt_list<dim> m_ref_pt_list;
  MyUtils::ref_pt_list<dim> m_ref_pt_list_tmp;
  string m_guess_str; /// a string containing the function
  string m_filename;
};


/** Default constructor
 */
template <int dim>
CBase<dim>::CBase( const string &xml_filename )
  :
  mpi_communicator(MPI_COMM_WORLD),
  m_computing_timer_log("benchmark.txt"),
  m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times ),
  m_ph(xml_filename),
  pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
{
  try
  {
    m_filename = m_ph.Get_Parameter( "filename" );
    m_guess_str = m_ph.Get_Parameter( "guess_fct" );

    m_omega = m_ph.Get_Physics("omega");
    m_gs = m_ph.Get_Physics("gs_1");
    m_QN1[0] = int(m_ph.Get_Physics("QN1", 0));
    m_QN1[1] = int(m_ph.Get_Physics("QN1", 1));
    m_QN1[2] = int(m_ph.Get_Physics("QN1", 2));

    m_xmin = m_ph.Get_Mesh("xrange", 0);
    m_xmax = m_ph.Get_Mesh("xrange", 1);
    m_ymin = m_ph.Get_Mesh("yrange", 0);
    m_ymax = m_ph.Get_Mesh("yrange", 1);
    m_zmin = m_ph.Get_Mesh("zrange", 0);
    m_zmax = m_ph.Get_Mesh("zrange", 1);
    m_global_refinement = m_ph.Get_Mesh("global_refinements", 0);

    m_ti = m_ph.Get_Algorithm("ti", 0);
    m_epsilon = m_ph.Get_Algorithm("epsilon");
    m_t[0] = m_ti;
    m_t[1] = m_ti;
    m_t_guess[0] = m_ti;
    m_t_guess[1] = m_ti;

    m_NA = int(m_ph.Get_Algorithm("NA", 0));
    m_Ndmu = m_ph.Get_Algorithm("Ndmu", 0);
    m_dmu = m_ph.Get_Algorithm("dmu", 0);
    m_df = m_ph.Get_Algorithm("df", 0);
  }
  catch ( const std::string info )
  {
    std::cerr << info << endl;
    MPI_Abort( mpi_communicator, 0 );
  }

  m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
  MPI_Comm_rank(mpi_communicator, &m_rank);

  m_counter = 0;
  m_final_error = 0;
}

template <int dim>
const double CBase<dim>::l2norm_t() const
{
  double retval = 0;
  for ( int i = 0; i < dim; i++ )
    retval += m_t[i] * m_t[i];
  return sqrt(retval);
}

template <int dim>
void CBase<dim>::screening()
{
  m_ref_pt_list_tmp.reset( 5, 20 );

  typename vector<MyUtils::ref_pt_list_item<dim>>::iterator it = m_ref_pt_list_tmp.m_list.begin();
  for ( ; it != m_ref_pt_list_tmp.m_list.end(); it++ )
  {
    size_t iter = 0;
    int status;
    double size;

    gsl_vector *x = gsl_vector_alloc(dim);
    gsl_vector *ss = gsl_vector_alloc(dim);
    gsl_vector_set_all (ss, 1.0);

    gsl_multimin_function my_func;

    my_func.n = dim;
    my_func.f = fun<DIMENSION>;
    my_func.params = this;

    for ( int i = 0; i < dim; i++ )
      gsl_vector_set (x, i, (*it).ti[i]);

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, dim);

    gsl_multimin_fminimizer_set (s, &my_func, x, ss );
    do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);
      if (status) break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-10);
    }
    while ( status == GSL_CONTINUE && iter < 10000 );
    //printf ("status = %s\n", gsl_strerror (status));

    (*it).f = s->fval;
    for ( int i = 0; i < dim; i++ )
    {
      (*it).t[i] = gsl_vector_get (s->x, i);
    }
    //it.DumpItem( std::cout );

    gsl_vector_free (x);
    gsl_vector_free (ss);
    gsl_multimin_fminimizer_free (s);
  }
  m_ref_pt_list_tmp.remove_zero_ref_points();
}

/**
 * @param[in] bool
 */
template <int dim>
int CBase<dim>::find_ortho_min( bool bdel_f_non_zero )
{
  compute_contributions();

  m_computing_timer.enter_section(__func__);

  if ( m_root )
  {
    double l2_norm_t_old = 0;
    double min = std::numeric_limits<double>::max();
    for ( int i = 0; i < dim; i++ )
      l2_norm_t_old += m_t[i] * m_t[i];
    l2_norm_t_old = sqrt(l2_norm_t_old);

    CBase<dim>::screening();

    m_ref_pt_list_tmp.remove_zero_ref_points();
    m_ref_pt_list_tmp.remove_duplicates();
    if (bdel_f_non_zero) m_ref_pt_list_tmp.remove_f_non_zero();
    //m_ref_pt_list_tmp.Dump( std::cout );

    typename vector<MyUtils::ref_pt_list_item<dim>>::iterator it = m_ref_pt_list_tmp.m_list.begin();
    for ( ; it != m_ref_pt_list_tmp.m_list.end(); it++ )
    {
      double tmp1 = fabs((*it).l2norm_t() - l2_norm_t_old);
      if ( tmp1 < min )
      {
        min = tmp1;
        for ( int i = 0; i < dim; i++ )
          m_t[i] = (*it).t[i];
      }
    }

    //WARNING: noch nicht getestet
    if ( dim == 2 )
      if ( m_t[0] / fabs(m_t[0]) < 0 )
        for ( int i = 0; i < dim; i++ ) m_t[i] *= -1.0;
  }

  int retval = m_ref_pt_list_tmp.m_list.empty();
  MPI_Bcast( m_t, dim, MPI_DOUBLE, 0, mpi_communicator );
  MPI_Bcast( &retval, dim, MPI_INT, 0, mpi_communicator );
  m_computing_timer.exit_section();
  return retval;
}

#endif