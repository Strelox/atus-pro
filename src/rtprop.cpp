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
 * @file rtprop.cpp Real time propagation for the Gross--Pitaevskii equation (Cartesian coordinates)
 * @brief Real time propagation for the Gross--Pitaevskii equation (Cartesian coordinates)
 * @author Želimir Marojević
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include "mpi.h"
#include "config.h"
#include "my_table.h"
#include "functions.h"
#include "MyParameterHandler.h"

///namespace realtime_propagation
namespace realtime_propagation
{
  using namespace std;
  using namespace dealii;

  /** BLA1
   */
  template <int dim>
  class MySolver
  {
  public:
    MySolver( const std::string & );
    ~MySolver();

    void run ();
    double Particle_Number( LA::MPI::Vector & );
    void Expectation_value_momentum( LA::MPI::Vector &, double * );
    void Expectation_value_position( LA::MPI::Vector &, double * );
  protected:
    void make_grid();
    void setup_system();
    void assemble_system();
    void assemble_rhs();

    void DoIter();
    void solve();
    void output_results ( string );
    void load( string );
    void save( string );

    MyParameterHandler m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector system_rhs;
    LA::MPI::Vector newton_update;
    LA::MPI::Vector m_Psi; // Psi(t)
    LA::MPI::Vector m_Psi_t; // Psi(t)
    LA::MPI::Vector m_workspace;
    LA::MPI::Vector m_workspace_ng;
    Vector<double> m_error_per_cell;

    ConditionalOStream pcout;
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;

    bool m_root;

    double m_gs;
    double m_t;
    double m_dt;
    std::vector<double> m_omega;
    double m_res;

    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;
    double m_zmin;
    double m_zmax;

    int m_rank;

    /// total number of outputs
    unsigned m_NA;
    /// number of intermediate time steps between two outputs
    unsigned m_NK;
    /// number of initial global grid refinements
    unsigned m_global_refinement;

    MyTable m_table;
    string m_filename;
  };

  /** Default constructor
   */
  template <int dim>
  MySolver<dim>::MySolver ( const std::string &xml_filename )
    :
    m_ph(xml_filename),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(MY_FE_DEGREE), 2),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    m_computing_timer_log("benchmark.txt"),
    m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput::cpu_and_wall_times )
  {
    try
    {
      m_filename = m_ph.Get_Parameter( "filename" );

      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1", 0);

      m_xmin = m_ph.Get_Mesh("xrange", 0);
      m_xmax = m_ph.Get_Mesh("xrange", 1);
      m_ymin = m_ph.Get_Mesh("yrange", 0);
      m_ymax = m_ph.Get_Mesh("yrange", 1);
      m_zmin = m_ph.Get_Mesh("zrange", 0);
      m_zmax = m_ph.Get_Mesh("zrange", 1);

      m_NA = int(m_ph.Get_Algorithm("NA", 0));
      m_NK = int(m_ph.Get_Algorithm("NK", 0));
      m_dt = m_ph.Get_Algorithm("dt", 0);
    }
    catch ( const std::string info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }

    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_rank);
  }

  /** Default Destructor
   */
  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }

  /** Computes the particle number of the dof function vec
   */
  template <int dim>
  double MySolver<dim>::Particle_Number( LA::MPI::Vector &vec )
  {
    m_computing_timer.enter_section(__func__);
    double tmp1 = 0.0;

    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; cell++)
    {
      if ( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for ( unsigned qp = 0; qp < n_q_points; qp++ )
          tmp1 += fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
      }
    }

    double retval;
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
    return retval;
  }

  template <int dim>
  void MySolver<dim>::make_grid ()
  {
    m_computing_timer.enter_section(__func__);
#if DIMENSION==2
    Point<dim, double> pt1( m_xmin, m_ymin );
    Point<dim, double> pt2( m_xmax, m_ymax );
#endif
#if DIMENSION==3
    Point<dim, double> pt1( m_xmin, m_ymin, m_zmin );
    Point<dim, double> pt2( m_xmax, m_ymax, m_zmax );
#endif

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);
    m_computing_timer.exit_section();
  }

  /**
   */
  template <int dim>
  void MySolver<dim>::setup_system()
  {
    m_computing_timer.enter_section(__func__);
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_t.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    newton_update.reinit(locally_owned_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(2), constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

    m_computing_timer.exit_section();
  }

  /** Compute the expectation value of the position operator \f$ \hat{x} \f$.
   * \f[
   * \int_{\Omega} vec^* \hat{x} vec dV
   * \f]
   * @param[in] vec parallel distributed dof vector
   * @returns retval expects a double array of size 3
   */
  template <int dim>
  void MySolver<dim>::Expectation_value_position( LA::MPI::Vector &vec, double *retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0, 0, 0}, JxWxn;

    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));
    Point<dim> spacept;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      if ( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for ( unsigned qp = 0; qp < n_q_points; qp++ )
        {
          JxWxn = fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
          spacept = fe_values.quadrature_point(qp);
          tmp[0] += spacept[0] * JxWxn;
          tmp[1] += spacept[1] * JxWxn;
#if dim == 3
          tmp[2] += spacept[2] * JxWxn;
#endif
        }
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }

  /** Computes the expectation value of the momentum operator \f$ \hat{p} \f$ in real space.
   * \f[
   * \int_{\Omega} vec^* \hat{p} vec dx
   * \f]
   * @param[in] vec parallel distributed dof vector
   * @param[out] retval expects a double array of size 3
   */
  template <int dim>
  void MySolver<dim>::Expectation_value_momentum( LA::MPI::Vector &vec, double *retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0, 0, 0}, JxW;

    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> vec_grads(n_q_points, vector<Tensor<1, dim>>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; cell++)
    {
      if ( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        fe_values.get_function_gradients( vec, vec_grads );

        for ( unsigned qp = 0; qp < n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          tmp[0] += JxW * (vec_vals[qp][0] * vec_grads[qp][1][0] - vec_vals[qp][1] * vec_grads[qp][0][0]);
          tmp[1] += JxW * (vec_vals[qp][0] * vec_grads[qp][1][1] - vec_vals[qp][1] * vec_grads[qp][0][1]);
#if dim == 3
          tmp[2] += JxW * (vec_vals[qp][0] * vec_grads[qp][1][2] - vec_vals[qp][1] * vec_grads[qp][0][2]);
#endif
        }
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }

  /** Output of data to vtu files
   * @param[in] path specifies the output folder
   */
  template <int dim>
  void MySolver<dim>::output_results ( string path )
  {
    m_computing_timer.enter_section(__func__);
    string filename;

    vector<std::string> solution_names;

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.add_data_vector (subdomain, "subdomain");
    data_out.build_patches ();

    filename = path + "solution-" + to_string(m_t) + ".vtu";
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );

    static vector<pair<double, string>> times_and_names;
    times_and_names.push_back (pair<double, string> (m_t, filename));
    ofstream pvd_output ("solution.pvd");
    data_out.write_pvd_record (pvd_output, times_and_names);

    m_computing_timer.exit_section();
  }

  /** Stores the triangulation and the vector m_Psi in the deal.ii native format.
   * @param[in] filename specifies the total path and file name of the output file
   */
  template <int dim>
  void MySolver<dim>::save( string filename )
  {
    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_Psi);
    triangulation.save( filename.c_str() );
  }

  /** loads the stored triangulation and the vector m_Psi from the deal.ii native format
   * @param[in] filename specifies the total path and file name of the input file
   */
  template <int dim>
  void MySolver<dim>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    setup_system();
    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(m_Psi);
  }

  template <int dim>
  void MySolver<dim>::assemble_system ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    BreedSolver::CPotential<dim> Potential ( m_omega );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_matrix = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> Psi_grad(n_q_points, vector<Tensor<1, dim>>(2));
    vector<Vector<double>> Psi_t(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> Psi_t_grad(n_q_points, vector<Tensor<1, dim>>(2));

    double JxW, pot = 0, tmp1a, tmp1b, tmp2, sum_re, sum_req, sum_im, sum_imq;

    const double fak2 = 0.5 * m_dt;
    const double fak4 = 0.25 * m_gs * m_dt;
    const double fak8 = 0.125 * m_gs * m_dt;


    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; cell++)
    {
      if ( cell->is_locally_owned() )
      {
        cell_matrix = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi, Psi);
        fe_values.get_function_gradients(m_Psi, Psi_grad);
        fe_values.get_function_values(m_Psi_t, Psi_t);
        fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

        for ( unsigned qp = 0; qp < n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          pot = Potential.value(fe_values.quadrature_point(qp));

          sum_re = Psi[qp][0] + Psi_t[qp][0];
          sum_im = Psi[qp][1] + Psi_t[qp][1];
          sum_req = sum_re * sum_re;
          sum_imq = sum_im * sum_im;
          tmp1a = fak8 * (sum_req + 3 * sum_imq);
          tmp1b = fak8 * (sum_imq + 3 * sum_req);
          tmp2 = fak4 * sum_re * sum_im;

          for ( unsigned i = 0; i < dofs_per_cell; i++ )
          {
            for ( unsigned j = 0; j < dofs_per_cell; j++ )
            {
              cell_matrix(i, j) += JxW * ((1.0 - tmp2) * fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) + (1.0 + tmp2) * fe_values[it].value(i, qp) * fe_values[it].value(j, qp)
                                          - tmp1a * fe_values[rt].value(i, qp) * fe_values[it].value(j, qp) - fak2 * (fe_values[rt].gradient(i, qp) * fe_values[it].gradient(j, qp) + pot * fe_values[rt].value(i, qp) * fe_values[it].value(j, qp))
                                          + tmp1b * fe_values[it].value(i, qp) * fe_values[rt].value(j, qp) + fak2 * (fe_values[it].gradient(i, qp) * fe_values[rt].gradient(j, qp) + pot * fe_values[it].value(i, qp) * fe_values[rt].value(j, qp)));
            }
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, local_dof_indices, system_matrix);
      }
    }
    system_matrix.compress(VectorOperation::add);
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::assemble_rhs ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    BreedSolver::CPotential<dim> Potential ( m_omega );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_rhs = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> Psi_grad(n_q_points, vector<Tensor<1, dim>>(2));
    vector<Vector<double>> Psi_t(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> Psi_t_grad(n_q_points, vector<Tensor<1, dim>>(2));

    double JxW, pot = 0, tmp1, sum_re, sum_im;

    const double fak2 = 0.5 * m_dt;
    const double fak8 = 0.125 * m_gs * m_dt;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      if ( cell->is_locally_owned() )
      {
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi, Psi);
        fe_values.get_function_gradients(m_Psi, Psi_grad);
        fe_values.get_function_values(m_Psi_t, Psi_t);
        fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

        for ( unsigned qp = 0; qp < n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          pot = Potential.value(fe_values.quadrature_point(qp));

          sum_re = Psi[qp][0] + Psi_t[qp][0];
          sum_im = Psi[qp][1] + Psi_t[qp][1];
          tmp1 = fak8 * (sum_re * sum_re + sum_im * sum_im);

          for ( unsigned i = 0; i < dofs_per_cell; i++ )
          {
            cell_rhs(i) += JxW * (-fak2 * ((Psi_grad[qp][1] + Psi_t_grad[qp][1]) * fe_values[rt].gradient(i, qp) + pot * sum_im * fe_values[rt].value(i, qp))
                                  + fak2 * ((Psi_grad[qp][0] + Psi_t_grad[qp][0]) * fe_values[it].gradient(i, qp) + pot * sum_re * fe_values[it].value(i, qp))
                                  + (Psi_t[qp][0] - Psi[qp][0]) * fe_values[rt].value(i, qp) - tmp1 * sum_im * fe_values[rt].value(i, qp)
                                  + (Psi_t[qp][1] - Psi[qp][1]) * fe_values[it].value(i, qp) + tmp1 * sum_re * fe_values[it].value(i, qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_rhs, local_dof_indices, system_rhs);
      }
    }
    system_rhs.compress(VectorOperation::add);

    m_workspace = system_rhs;
    VectorTools::integrate_difference ( dof_handler,  m_workspace, ZeroFunction<dim>(), m_error_per_cell,  QGauss<dim>(fe.degree + 2), VectorTools::L2_norm);
    tmp1 = m_error_per_cell.l2_norm();
    tmp1 = tmp1 * tmp1;
    MPI_Allreduce( &tmp1, &m_res, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_res = sqrt(m_res);
    m_computing_timer.exit_section();
  }

  /**
   * This procedure solves
   * @returns m_workspace contains the correction for the newton step
   * @param none
   */
  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);

    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, newton_update, system_rhs);

    constraints.distribute (newton_update);
    m_computing_timer.exit_section();
  }

  /**
   * This procedure implements the Newton method. We correct out initial guess m_Psi_t as long as the error of right hand side has not reached the threshold \f$ 10^{-16} \f$.
   */
  template<int dim>
  void MySolver<dim>::DoIter()
  {
    m_computing_timer.enter_section(__func__);

    m_Psi_t = 0;
    m_res = 0;
    assemble_rhs();
    pcout << "m_res = " << m_res << endl;
    do
    {
      assemble_system();
      solve();

      m_workspace_ng = m_Psi_t;
      m_workspace_ng.add( -1, newton_update );
      constraints.distribute(m_workspace_ng);
      m_Psi_t = m_workspace_ng;

      assemble_rhs();
      pcout << "m_res = " << m_res << endl;
    }
    while ( m_res > 1e-16 );
    m_t += m_dt;

    m_Psi = m_Psi_t;
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    double N;

    load( m_filename );

    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    pcout << "min_cell_diameter = " << min_cell_diameter << "\n";
    pcout << "max_cell_diameter = " << max_cell_diameter << "\n";
    pcout << "dt/dx^2 == " << m_dt / (min_cell_diameter * min_cell_diameter) << endl;

    double p[] = {0, 0, 0};
    double pos[] = {0, 0, 0};

    system_rhs = m_Psi;
    constraints.distribute(system_rhs);
    m_Psi = system_rhs;

    output_results("");

    N = Particle_Number(m_Psi);
    //pcout << "N == " << N << endl;
    Expectation_value_position( m_Psi, pos );
    Expectation_value_momentum( m_Psi, p );

    pcout << "t == " << m_t << endl;
    pcout << "N == " << N << endl;
    pcout << "p == " << p[0] / N << ", " << p[1] / N << ", " << p[2] / N << endl;
    pcout << "pos == " << pos[0] / N << ", " << pos[1] / N << ", " << pos[2] / N << endl;

    for ( unsigned i = 1; i <= m_NA; i++ )
    {
      for ( unsigned j = 1; j <= m_NK; j++ )
      {
        pcout << "t == " << m_t << endl;
        DoIter();
      }

      N = Particle_Number(m_Psi);
      Expectation_value_position( m_Psi, pos );
      Expectation_value_momentum( m_Psi, p );

      pcout << "N == " << N << endl;
      pcout << "p == " << p[0] / N << ", " << p[1] / N << ", " << p[2] / N << endl;
      pcout << "pos == " << pos[0] / N << ", " << pos[1] / N << ", " << pos[2] / N << endl;

      output_results("");
    }
  }
} // end of namespace

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    realtime_propagation::MySolver<DIMENSION> solver("params.xml");
    solver.run();
  }
  return EXIT_SUCCESS;
}
