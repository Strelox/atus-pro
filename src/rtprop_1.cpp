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

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/function_parser.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include "config.h"
#include "my_table.h"
#include "functions.h"
#include "MyParameterHandler.h"

///namespace realtime_propagation_1
namespace realtime_propagation_1
{
  using namespace std;
  using namespace dealii;

  class MySolver
  {
  public:
    MySolver( const std::string & );
    ~MySolver();

    void run ();
    double Particle_Number( Vector<double> & );
    void Expectation_value_momentum( Vector<double> &, double & );
    void Expectation_value_position( Vector<double> &, double & );

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
    Triangulation<1> triangulation;
    FESystem<1> fe;
    DoFHandler<1> dof_handler;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> system_rhs;
    Vector<double> newton_update;
    Vector<double> m_Psi; // Psi(t)
    Vector<double> m_Psi_t; // Psi(t)
    Vector<double> m_workspace;
    Vector<double> m_error_per_cell;

    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;

    double m_gs;
    double m_t;
    double m_dt;
    std::vector<double> m_omega;
    double m_res;

    double m_xmin;
    double m_xmax;

    /// total number of outputs
    unsigned m_NA;
    /// number of intermediate time steps between two outputs
    unsigned m_NK;

    unsigned m_global_refinement;

    MyTable m_table;
    string m_filename;
  };

  /** Default constructor
   */
  MySolver::MySolver ( const std::string &xml_filename )
    :
    m_ph(xml_filename),
    triangulation (),
    fe (FE_Q<1>(MY_FE_DEGREE), 2),
    dof_handler (triangulation),
    m_computing_timer_log("benchmark.txt"),
    m_computing_timer( m_computing_timer_log, TimerOutput::summary, TimerOutput::cpu_and_wall_times )
  {
    try
    {
      m_filename = m_ph.Get_Parameter( "filename" );

      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1", 0);

      m_xmin = m_ph.Get_Mesh("xrange", 0);
      m_xmax = m_ph.Get_Mesh("xrange", 1);

      m_global_refinement = m_ph.Get_Algorithm("global_refinements", 0);
      m_NA = int(m_ph.Get_Algorithm("NA", 0));
      m_NK = int(m_ph.Get_Algorithm("NK", 0));
      m_dt = m_ph.Get_Algorithm("dt", 0);
    }
    catch ( const std::string info )
    {
      std::cerr << info << endl;
      exit(0);
    }
  }

  /** Default constructor
   */
  MySolver::~MySolver ()
  {
    dof_handler.clear ();
  }

  /** Computes the particle number of the dof function vec
   */
  double MySolver::Particle_Number( Vector<double> &vec )
  {
    m_computing_timer.enter_section(__func__);
    double retval = 0;

    const QGauss<1>  quadrature_formula(fe.degree + 1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec, vec_vals );
      for ( unsigned qp = 0; qp < n_q_points; qp++ )
        retval += fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
    }
    m_computing_timer.exit_section();
    return retval;
  }

  void MySolver::make_grid ()
  {
    m_computing_timer.enter_section(__func__);

    Point<1, double> pt1( m_xmin );
    Point<1, double> pt2( m_xmax );

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);

    m_computing_timer.exit_section();
  }

  void MySolver::setup_system()
  {
    m_computing_timer.enter_section(__func__);
    dof_handler.distribute_dofs (fe);
    DoFRenumbering::component_wise (dof_handler);

    m_Psi.reinit (dof_handler.n_dofs());
    m_workspace.reinit (dof_handler.n_dofs());
    m_Psi_t.reinit (dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    newton_update.reinit(dof_handler.n_dofs());
    m_error_per_cell.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);
    m_computing_timer.exit_section();
  }

  void MySolver::Expectation_value_position( Vector<double> &vec, double &retval )
  {
    m_computing_timer.enter_section(__func__);
    double JxWxn;
    retval = 0;

    const QGauss<1>  quadrature_formula(fe.degree + 1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));
    Point<1> spacept;

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec, vec_vals );
      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        JxWxn = fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
        spacept = fe_values.quadrature_point(qp);
        retval += spacept[0] * JxWxn;
      }
    }
    m_computing_timer.exit_section();
  }

  void MySolver::Expectation_value_momentum( Vector<double> &vec, double &retval )
  {
    m_computing_timer.enter_section(__func__);

    retval = 0;

    const QGauss<1>  quadrature_formula(fe.degree + 1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, 1>>> vec_grads(n_q_points, vector<Tensor<1, 1>>(2));

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec, vec_vals );
      fe_values.get_function_gradients( vec, vec_grads );

      for (unsigned q_point = 0; q_point < n_q_points; ++q_point)
      {
        retval += fe_values.JxW(q_point) * (vec_vals[q_point][0] * vec_grads[q_point][1][0] - vec_vals[q_point][1] * vec_grads[q_point][0][0]);
      }
    }
    m_computing_timer.exit_section();
  }

  /** Output of data to vtu files
   * @param[in] path specifies the output folder
   */
  void MySolver::output_results ( string path )
  {
    m_computing_timer.enter_section(__func__);

    vector<std::string> solution_names;

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.build_patches ();

    string filename = path + "solution-" + to_string(m_t) + ".gnuplot";
    ofstream output (filename);
    data_out.write_gnuplot ( output );

    m_computing_timer.exit_section();
  }

  void MySolver::save( string filename )
  {
    Vector<double> tmp( 2 * dof_handler.n_dofs());
    tmp = 0;

    for ( unsigned i = 0; i < dof_handler.n_dofs(); i++ )
    {
      tmp[i] = m_Psi[i];
    }

    ofstream out(filename);
    tmp.block_write(out);
  }

  void MySolver::load( string filename )
  {
    ifstream in(filename);
    m_Psi.block_read(in);
  }

  void MySolver::assemble_system ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<1> quadrature_formula(fe.degree + 1);

    CPotential<1> Potential ( m_omega );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_matrix = 0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, 1>>> Psi_grad(n_q_points, vector<Tensor<1, 1>>(2));
    vector<Vector<double>> Psi_t(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, 1>>> Psi_t_grad(n_q_points, vector<Tensor<1, 1>>(2));

    double JxW, pot = 0, tmp1a, tmp1b, tmp2, sum_re, sum_req, sum_im, sum_imq;

    const double fak2 = 0.5 * m_dt;
    const double fak4 = 0.25 * m_gs * m_dt;
    const double fak8 = 0.125 * m_gs * m_dt;

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
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

        for (unsigned i = 0; i < dofs_per_cell; i++ )
        {
          for (unsigned j = 0; j < dofs_per_cell; j++ )
          {
            cell_matrix(i, j) += JxW * ((1.0 - tmp2) * fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) + (1.0 + tmp2) * fe_values[it].value(i, qp) * fe_values[it].value(j, qp)
                                        - tmp1a * fe_values[rt].value(i, qp) * fe_values[it].value(j, qp) - fak2 * (fe_values[rt].gradient(i, qp) * fe_values[it].gradient(j, qp) + pot * fe_values[rt].value(i, qp) * fe_values[it].value(j, qp))
                                        + tmp1b * fe_values[it].value(i, qp) * fe_values[rt].value(j, qp) + fak2 * (fe_values[it].gradient(i, qp) * fe_values[rt].gradient(j, qp) + pot * fe_values[it].value(i, qp) * fe_values[rt].value(j, qp)));
          }
        }
      }

      cell->get_dof_indices (local_dof_indices);
      for (unsigned i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned j = 0; j < dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
      }
    }

    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(2), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, newton_update, system_rhs);

    m_computing_timer.exit_section();
  }

  void MySolver::assemble_rhs ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<1> quadrature_formula(fe.degree + 1);

    CPotential<1> Potential (m_omega);

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_rhs = 0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, 1>>> Psi_grad(n_q_points, vector<Tensor<1, 1>>(2));
    vector<Vector<double>> Psi_t(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, 1>>> Psi_t_grad(n_q_points, vector<Tensor<1, 1>>(2));

    double JxW, pot = 0, tmp1, sum_re, sum_im;

    const double fak2 = 0.5 * m_dt;
    const double fak8 = 0.125 * m_gs * m_dt;

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi, Psi);
      fe_values.get_function_gradients(m_Psi, Psi_grad);
      fe_values.get_function_values(m_Psi_t, Psi_t);
      fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

      cell->get_dof_indices (local_dof_indices);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        JxW = fe_values.JxW(qp);
        pot = Potential.value(fe_values.quadrature_point(qp));

        sum_re = Psi[qp][0] + Psi_t[qp][0];
        sum_im = Psi[qp][1] + Psi_t[qp][1];
        tmp1 = fak8 * (sum_re * sum_re + sum_im * sum_im);

        for (unsigned i = 0; i < dofs_per_cell; i++ )
        {
          system_rhs(local_dof_indices[i]) += JxW * (-fak2 * ((Psi_grad[qp][1] + Psi_t_grad[qp][1]) * fe_values[rt].gradient(i, qp) + pot * sum_im * fe_values[rt].value(i, qp))
                                                     + fak2 * ((Psi_grad[qp][0] + Psi_t_grad[qp][0]) * fe_values[it].gradient(i, qp) + pot * sum_re * fe_values[it].value(i, qp))
                                                     + (Psi_t[qp][0] - Psi[qp][0]) * fe_values[rt].value(i, qp) - tmp1 * sum_im * fe_values[rt].value(i, qp)
                                                     + (Psi_t[qp][1] - Psi[qp][1]) * fe_values[it].value(i, qp) + tmp1 * sum_re * fe_values[it].value(i, qp));
        }
      }
    }
    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(2), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, newton_update, system_rhs);

    m_workspace = system_rhs;
    VectorTools::integrate_difference ( dof_handler,  m_workspace, ZeroFunction<1>(2), m_error_per_cell,  QGauss<1>(fe.degree + 1), VectorTools::L2_norm);
    m_res = m_error_per_cell.l2_norm();
    m_computing_timer.exit_section();
  }

  void MySolver::solve ()
  {
    m_computing_timer.enter_section(__func__);

    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (newton_update, system_rhs);

    m_computing_timer.exit_section();
  }

  void MySolver::DoIter()
  {
    m_computing_timer.enter_section(__func__);

    m_Psi_t = 0;
    m_res = 0;
    assemble_rhs();
    //cout << "m_res = " << m_res << endl;
    do
    {
      assemble_system();
      solve();

      m_Psi_t.add( -1, newton_update );

      assemble_rhs();
      //cout << "m_res = " << m_res << endl;
    }
    while ( m_res > 1e-16 );
    m_t += m_dt;

    m_Psi = m_Psi_t;
    m_computing_timer.exit_section();
  }

  void MySolver::run()
  {
    double N;

    make_grid();
    setup_system();

    load( m_filename );

    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    cout << "min_cell_diameter = " << min_cell_diameter << "\n";
    cout << "max_cell_diameter = " << max_cell_diameter << "\n";
    cout << "dt/dx^2 == " << m_dt / (min_cell_diameter * min_cell_diameter) << endl;

    double p;
    double pos;

    output_results("");

    N = Particle_Number(m_Psi);
    //cout << "N == " << N << endl;
    Expectation_value_position( m_Psi, pos );
    Expectation_value_momentum( m_Psi, p );

    cout << "t == " << m_t << endl;
    cout << "N == " << N << endl;
    cout << "p == " << p / N << endl;
    cout << "pos == " << pos / N << endl;

    for ( unsigned i = 1; i <= m_NA; i++ )
    {
      for ( unsigned j = 1; j <= m_NK; j++ )
      {
        cout << "t == " << m_t << endl;
        DoIter();
      }

      N = Particle_Number(m_Psi);
      Expectation_value_position( m_Psi, pos );
      Expectation_value_momentum( m_Psi, p );

      cout << "N == " << N << endl;
      cout << "p == " << p / N << endl;
      cout << "pos == " << pos / N << endl;

      output_results("");
    }
  }
} // end of namespace

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  realtime_propagation_1::MySolver solver("params.xml");
  solver.run();
  return EXIT_SUCCESS;
}
