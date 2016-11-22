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
 * @file my_table.h
 * @brief A simple class which manages a table with double entries.
 * @author Želimir Marojević
 */

#ifndef _MYTABLE_
#define _MYTABLE_

#include <vector>
#include <string>
#include <map>

using namespace std;

typedef map<string, double> columns;
typedef vector<columns> table;

class MyTable
{
public:
  MyTable();
  ~MyTable();

  friend ostream &operator<<( ostream &, MyTable & );

  void clear();
  columns &new_line();
  void insert( columns &, const string, const double );
  void dump_2_file( const string );

  static const string COUNTER;
  static const string RES;
  static const string RESP;
  static const string RES_OVER_RESP;
  static const string PARTICLE_NUMBER;
  static const string MU;
  static const string GS;
  static const string l2norm_t;
  static const string t1;
  static const string t2;
  static const string total_no_cells;
  static const string total_no_active_cells;
  static const string STEPSIZE;
  static const string STATUS;
  static const string time;
  static const string ev_position_x;
  static const string ev_position_y;
  static const string ev_position_z;
  static const string ev_momentum_x;
  static const string ev_momentum_y;
  static const string ev_momentum_z;
  static const string ev_difference;

  void save_txt( string );

  table m_table;
protected:
  vector<string> m_order;
};

ostream &operator<<( ostream &stream, MyTable &obj );
#endif