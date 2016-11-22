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
 * @file my_table.cpp
 * @brief A simple class which manages a table with double entries.
 * @author Želimir Marojević
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "my_table.h"

const string MyTable::COUNTER = "Counter";
const string MyTable::RES = "res";
const string MyTable::RESP = "resp";
const string MyTable::RES_OVER_RESP = "|res / resp|";
const string MyTable::PARTICLE_NUMBER = "N";
const string MyTable::MU = "mu";
const string MyTable::GS = "gs";
const string MyTable::l2norm_t = "l2 norm t";
const string MyTable::t1 = "t1";
const string MyTable::t2 = "t2";
const string MyTable::total_no_cells = "total no of cells";
const string MyTable::total_no_active_cells = "total no of active cells";
const string MyTable::STEPSIZE = "stepsize";
const string MyTable::STATUS = "status";
const string MyTable::time = "time";
const string MyTable::ev_position_x = "expectation value position x";
const string MyTable::ev_position_y = "expectation value position y";
const string MyTable::ev_position_z = "expectation value position z";
const string MyTable::ev_momentum_x = "expectation value momentum x";;
const string MyTable::ev_momentum_y = "expectation value momentum y";;
const string MyTable::ev_momentum_z = "expectation value momentum z";;

MyTable::MyTable()
{
}

MyTable::~MyTable()
{
  clear();
}

void MyTable::clear()
{
  table::iterator row = m_table.begin();
  for ( ; row != m_table.end(); ++row )
  {
    (*row).clear();
  }
  m_table.clear();
  m_order.clear();
}

columns &MyTable::new_line()
{
  m_order.clear();
  columns cols;
  m_table.push_back(cols);
  columns &ref = m_table.back();

  return ref;
}

void MyTable::insert(columns &cols, const string str, const double val )
{
  cols[str] = val;
  m_order.push_back(str);
}

void MyTable::save_txt( string path )
{
  ofstream out( path );
  out << *this;
}

ostream &operator<<( ostream &stream, MyTable &obj )
{
  columns &cols = obj.m_table.back();

  for ( unsigned int i = 0; i < obj.m_order.size(); i++ )
  {
    stream << setw(14) << std::left << obj.m_order[i] << " == " << cols[obj.m_order[i]] << endl;
  }
  stream << "\n";
  return stream;
}

void MyTable::dump_2_file( const string path )
{
  if ( m_table.size() == 0 ) return;

  unsigned int c, r;
  ofstream ofs( path );

  ofs << "# ";
  for ( c = 0; c < m_order.size() - 1; c++ )
  {
    ofs << m_order[c] << ";";
  }
  ofs << m_order[c] << "\n";

  for ( r = 0; r < m_table.size(); r++ )
  {
    columns &cols = m_table[r];
    for ( c = 0; c < m_order.size() - 1; c++ )
    {
      ofs << cols[m_order[c]] << ";";
    }
    ofs << cols[m_order[c]] << "\n";
  }
  ofs.close();
}