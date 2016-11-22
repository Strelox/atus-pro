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
* @file breed_sob_1.cpp
* @author Želimir Marojević
* @brief This program solves the stationary Gross--Pitaevskii equation for a real wave function in 1D using the sobolev gradient method.
*
* A wider description of the method can be found in the book of John Neuberger, "Sobolev Gradients and Differential Equations", Springer.
* A shorter description among others by Kazemi and Eckart can be found on arxiv: http://arxiv.org/abs/0906.3206
*/

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <fstream>
#include <iostream>
#include <cmath>
#include <array>

#include "config.h"
#include "functions.h"
#include "my_table.h"
#include "MyParameterHandler.h"

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

  class MySolver
  {
  public:
    MySolver( const std::string );
    virtual ~MySolver();

    void run ();

  protected:

    double m_xmin, m_xmax;
    double m_res, m_res_old, m_resp;
    double m_final_error;
    double m_N;
    double m_df;
    double m_mu;
    double m_gs;
    vector<double> m_omega;
    vector<double> m_epsilon;

    unsigned m_counter;
    unsigned m_global_refinement;
    unsigned m_total_no_cells;
    unsigned m_total_no_active_cells;
    unsigned m_NA;

    int DoIter( string = "" );

    void make_grid();
    void setup_system();
    void assemble_rhs();
    void assemble_system();
    void compute_Psi_sob();
    void compute_mu();
    void save( string );
    //void load( string );
    void Project_gradient();
    void Interpolate_R_to_C();

    double Particle_Number( const Vector<double> & );

    void solve();
    void compute_E_lin( const Vector<double> &, double &, double &, double & );

    void output_results ( string, string = "step" );

    MyParameterHandler m_ph;
    Triangulation<1> triangulation;
    FE_Q<1> fe;
    DoFHandler<1> dof_handler;
    FESystem<1> fe_2;
    DoFHandler<1> dof_handler_2;

    SparsityPattern m_sparsity_pattern, m_sparsity_pattern_2;
    SparseMatrix<double> m_system_matrix, m_system_matrix_2;

    Vector<double> m_system_rhs, m_system_rhs_2;
    Vector<double> m_Psi;
    Vector<double> m_Psi_C;
    Vector<double> m_sob_grad;
    Vector<double> m_Psi_sob;
    Vector<double> m_error_per_cell;

    string m_guess_str;
    MyTable m_table;
    MyTable m_results;
  };

  /**
   * Constructor
   */
  MySolver::MySolver ( const std::string xmlfilename )
    :
    m_ph(xmlfilename),
    triangulation (typename Triangulation<1>::MeshSmoothing(Triangulation<1>::limit_level_difference_at_vertices | Triangulation<1>::smoothing_on_refinement | Triangulation<1>::smoothing_on_coarsening), true),
    fe (MY_FE_DEGREE),
    dof_handler (triangulation),
    fe_2 (FE_Q<1>(MY_FE_DEGREE), 2),
    dof_handler_2 (triangulation)
  {
    try
    {
      m_guess_str = m_ph.Get_Parameter( "guess_fct" );

      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1", 0);

      m_xmin = m_ph.Get_Mesh("xrange", 0);
      m_xmax = m_ph.Get_Mesh("xrange", 1);
      m_global_refinement = m_ph.Get_Mesh("global_refinements", 0);

      m_epsilon = m_ph.Get_Algorithm("epsilon");
      m_NA = int(m_ph.Get_Algorithm("NA", 0));
      m_df = m_ph.Get_Algorithm("df", 0);
    }
    catch ( const std::string info )
    {
      std::cerr << info << endl;
      exit(0);
    }

    m_counter = 0;
    m_final_error = 0;
  }

  MySolver::~MySolver ()
  {
    dof_handler.clear();
    dof_handler_2.clear();
  }

  void MySolver::compute_E_lin( const Vector<double> &vec, double &T, double &N, double &W )
  {
    CPotential<1> Potential( m_omega );
    const QGauss<1>  quadrature_formula(fe.degree + 1);
    FEValues<1> fe_values (fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals(n_q_points);
    vector<Tensor<1, 1>> grad(n_q_points);

    T = 0;
    N = 0;
    W = 0;

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell != endc; cell++ )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(vec, vals);
      fe_values.get_function_gradients(vec, grad);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double vec_val_q = vals[qp] * vals[qp];
        double JxW = fe_values.JxW(qp);
        T += JxW * ( grad[qp] * grad[qp] + Potential.value(fe_values.quadrature_point(qp)) * vec_val_q );
        N += JxW * vec_val_q;
        W += JxW * vec_val_q * vec_val_q;
      }
    }
  }

  double MySolver::Particle_Number( const Vector<double> &vec )
  {
    double retval = 0;

    const QGauss<1>  quadrature_formula(fe.degree + 1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals(n_q_points);

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell != endc; cell++ )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec, vals );

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
        retval += fe_values.JxW(qp) * (vals[qp] * vals[qp]);
    }
    return retval;
  }

  /**
  * Project_gradient() computes the projection of the sobolev gradient onto the subspace of constant norm.
  * This procedure orthogonalizes the sobolev gradient with respect to the wave function, so that
  * the norm of the wave function m_Psi is not changed for sufficiently small step lenght.
  */
  void MySolver::Project_gradient()
  {
    double sum[] = {0, 0};

    const QGauss<1> quadrature_formula(fe.degree + 1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals_1(n_q_points);
    vector<double> vals_2(n_q_points);
    vector<double> vals_3(n_q_points);

    typename DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell != endc; cell++ )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( m_Psi_sob, vals_1 );
      fe_values.get_function_values( m_sob_grad, vals_2 );
      fe_values.get_function_values( m_Psi, vals_3 );

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        sum[0] += JxW * (vals_3[qp] * vals_2[qp]);
        sum[1] += JxW * (vals_3[qp] * vals_1[qp]);
      }
    }
    m_sob_grad.add( -sum[0] / sum[1], m_Psi_sob );
  }

  /**
  * assemble_rhs () assembles the right hand side of the system of linear equations which is necessary to compute the sobolev gradient.
  * Mainly we compute the integral of the L_2 gradient multplied by the test functions.
  */
  void MySolver::assemble_rhs ()
  {
    CPotential<1> Potential( m_omega );
    const QGauss<1> quadrature_formula(fe.degree + 1);

    m_system_rhs = 0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    vector<double> vals(n_q_points);
    vector<Tensor<1, 1>> grads(n_q_points);

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell != endc; cell++ )
    {
      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);

      fe_values.get_function_values(m_Psi, vals);
      fe_values.get_function_gradients(m_Psi, grads);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        double Q1 = Potential.value(fe_values.quadrature_point(qp)) + m_gs * (vals[qp] * vals[qp]);

        for ( unsigned i = 0; i < dofs_per_cell; ++i)
          m_system_rhs(local_dof_indices[i]) += JxW * (grads[qp] * fe_values.shape_grad(i, qp) + Q1 * vals[qp] * fe_values.shape_value(i, qp));
      }
    }

    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, m_system_matrix, m_sob_grad, m_system_rhs);
  }

  /**
  * assemble_system () assembles the matrix of the operator (1-Laplace) in the weak formulation.
  */
  void MySolver::assemble_system ()
  {
    CPotential<1> Potential( m_omega );
    const QGauss<1> quadrature_formula(fe.degree + 1);

    m_system_matrix = 0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell != endc; cell++ )
    {
      cell_matrix = 0;
      fe_values.reinit (cell);

      cell->get_dof_indices (local_dof_indices);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        //double Q1 = Potential.value(fe_values.quadrature_point(qp));

        for ( unsigned i = 0; i < dofs_per_cell; i++ )
          for ( unsigned j = 0; j < dofs_per_cell; j++ )
            cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
      }

      for ( unsigned i = 0; i < dofs_per_cell; i++ )
      {
        for ( unsigned j = 0; j < dofs_per_cell; j++ )
          m_system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
      }
    }

    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, m_system_matrix, m_sob_grad, m_system_rhs);
  }

  /**
  * computes the projection of m_Psi from the L_2 space onto the space W^1,2.
  */
  void MySolver::compute_Psi_sob ()
  {
    const QGauss<1> quadrature_formula(fe.degree + 1);

    m_system_matrix = 0;
    m_system_rhs = 0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    vector<double> vals(n_q_points);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell != endc; ++cell )
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi, vals);

      cell->get_dof_indices (local_dof_indices);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        double val = vals[qp];

        for ( unsigned i = 0; i < dofs_per_cell; i++ )
        {
          cell_rhs(i) = JxW * val * fe_values.shape_value(i, qp);
          for ( unsigned j = 0; j < dofs_per_cell; j++ )
            cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));

        }
      }

      for ( unsigned i = 0; i < dofs_per_cell; i++ )
      {
        m_system_rhs(local_dof_indices[i]) += cell_rhs(i);
        for ( unsigned j = 0; j < dofs_per_cell; j++ )
          m_system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
      }
    }

    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, m_system_matrix, m_Psi_sob, m_system_rhs);

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(m_system_matrix);
    A_direct.vmult (m_Psi_sob, m_system_rhs);
  }

  void MySolver::compute_mu()
  {
    CPotential<1> Potential( m_omega );
    const QGauss<1>  quadrature_formula(fe.degree + 1);
    FEValues<1> fe_values (fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals(n_q_points);
    vector<Tensor<1, 1>> grads(n_q_points);

    m_mu = 0;
    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell != endc; cell++ )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi, vals);
      fe_values.get_function_gradients(m_Psi, grads);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double uq = vals[qp] * vals[qp];
        m_mu += fe_values.JxW(qp) * (grads[qp] * grads[qp] + (Potential.value(fe_values.quadrature_point(qp)) + m_gs * uq) * uq);
      }
    }
  }

  void MySolver::solve ()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(m_system_matrix);
    A_direct.vmult (m_sob_grad, m_system_rhs);
  }

  void MySolver::make_grid ()
  {
    Point<1, double> pt1(m_xmin);
    Point<1, double> pt2(m_xmax);

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);
  }

  void MySolver::setup_system ()
  {
    dof_handler.distribute_dofs (fe);

    m_Psi.reinit (dof_handler.n_dofs());
    m_Psi_sob.reinit (dof_handler.n_dofs());
    m_system_rhs.reinit(dof_handler.n_dofs());
    m_sob_grad.reinit (dof_handler.n_dofs());
    m_error_per_cell.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);

    m_sparsity_pattern.copy_from(dsp);
    m_system_matrix.reinit (m_sparsity_pattern);

    // stuff for the second dof handler
    dof_handler_2.distribute_dofs (fe_2);

    m_Psi_C.reinit (dof_handler_2.n_dofs());
    m_system_rhs_2.reinit(dof_handler_2.n_dofs());

    DynamicSparsityPattern dsp_2(dof_handler_2.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler_2, dsp_2);

    m_sparsity_pattern_2.copy_from(dsp_2);
    m_system_matrix_2.reinit (m_sparsity_pattern_2);
  }

  void MySolver::output_results ( string path, string prefix )
  {
    string filename = path + prefix + "-" + Utilities::int_to_string (m_counter, 5) + ".gnuplot";

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi, "Psi");
    data_out.add_data_vector (m_sob_grad, "L_2 gradient");
    data_out.add_data_vector (m_error_per_cell, "error_per_cell");
    data_out.build_patches ();

    ofstream output (filename.c_str());
    data_out.write_gnuplot (output);
  }

  int MySolver::DoIter ( string path )
  {
    int retval = Status::SUCCESS;

    m_table.clear();

    assemble_rhs();

    m_res = 0;
    m_counter = 0;
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "- " << m_counter << endl;

      assemble_system();
      solve();

      compute_Psi_sob();

      Project_gradient();

      m_res = m_sob_grad.l2_norm();

      m_Psi.add( -m_df*0.001, m_sob_grad);

      //compute_mu(m_mu);
      m_N = Particle_Number(m_Psi);

      if ( fabs(m_N - 1) > 1e-5 ) m_Psi *= 1 / sqrt(m_N);

      assemble_rhs();

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      if ( m_counter % m_NA == 0 ) output_results(path);

      columns &cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      //m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      //m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      cout << m_table;
      cout << "m_final_error =" << m_final_error << endl;
      if ( m_res < m_epsilon[1] )
      {
        retval = Status::SUCCESS;
        break;
      }

      m_counter++;
    }
    while ( true );

    m_N = Particle_Number(m_Psi);

    if ( m_N < 1e-5 ) retval = Status::ZERO_SOL;

    string filename = path + "log.csv";
    m_table.dump_2_file(filename);

    return retval;
  }

  void MySolver::run()
  {
    int status;
    double T, N, W;

    make_grid();
    setup_system();

    unsigned QN[3] = {};
    CEigenfunctions<1> Ef1( QN, m_omega );
    VectorTools::interpolate (dof_handler, Ef1, m_Psi );

    output_results("", "guess");

    compute_E_lin( m_Psi, T, N, W );
    m_Psi *= sqrt(1 / N);

    cout << setprecision(9);
    cout << "T = " << T << endl;
    cout << "N = " << N << endl;
    cout << "W = " << W << endl;
    //cout << "m_mu = " << m_mu << endl;

    status = DoIter("");

    cout << "m_final_error of m_Psi: " << m_final_error << endl;

    if ( status == Status::SUCCESS )
    {
      Interpolate_R_to_C();
      output_results("","final");
    }

    ofstream ofs("log.txt");
    ofs << m_table;
  }

  /**
  * computes the interpolation of the real wave function m_Psi onto a FE space representing a complex wave function with imaginary part zero.
  */
  void MySolver::Interpolate_R_to_C()
  {
    const QGauss<1> quadrature_formula(fe.degree + 1);
    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    m_system_rhs_2 = 0;
    m_system_matrix_2 = 0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_JxW_values);
    FEValues<1> fe_values_2 (fe_2, quadrature_formula, update_values | update_JxW_values);

    const unsigned dofs_per_cell = fe_2.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    vector<double> vals(n_q_points);
    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    typename DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    typename DoFHandler<1>::active_cell_iterator cell_2 = dof_handler_2.begin_active();
    for ( ; cell != endc; ++cell, ++cell_2 )
    {
      cell_rhs = 0;
      cell_matrix = 0;

      fe_values.reinit (cell);
      fe_values_2.reinit (cell_2);
      fe_values.get_function_values(m_Psi, vals);

      cell_2->get_dof_indices (local_dof_indices);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values_2.JxW(qp);
        double tmp1 = vals[qp];

        for ( unsigned i = 0; i < dofs_per_cell; i++ )
        {
          cell_rhs(i) += JxW * tmp1 * fe_values_2[rt].value(i, qp);
          for ( unsigned j = 0; j < dofs_per_cell; j++ )
          {
            cell_matrix(i, j) += JxW * (fe_values_2[rt].value(i, qp) * fe_values_2[rt].value(j, qp) + fe_values_2[it].value(i, qp) * fe_values_2[it].value(j, qp));
          }
        }

        for ( unsigned i = 0; i < dofs_per_cell; i++ )
        {
          m_system_rhs_2(local_dof_indices[i]) += cell_rhs(i);
          for ( unsigned j = 0; j < dofs_per_cell; j++ )
          {
            m_system_matrix_2.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
          }
        }
      }
    }

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(m_system_matrix_2);
    A_direct.vmult (m_Psi_C, m_system_rhs_2);

    ofstream out("Cfinal.bin");
    m_Psi_C.block_write(out);
  }

  void MySolver::save( string filename )
  {
  }
} // end of namespace

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    BreedSolver::MySolver solver("params.xml");
    solver.run ();
  }
  return EXIT_SUCCESS;
}