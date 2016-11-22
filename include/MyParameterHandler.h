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


#include "tinyxml2.h"
#include <string>
#include <vector>
#include <map>
#include <list>

#ifndef __class_MyParameterHandler__
#define __class_MyParameterHandler__

class MyParameterHandler
{
public:
  MyParameterHandler( const std::string );
  virtual ~MyParameterHandler() {};

  std::string Get_Parameter( const std::string );

  void Set_Physics( const std::string, const std::vector<double> & );
  double Get_Physics( const std::string, const int );
  std::vector<double> Get_Physics( const std::string );

  double Get_Mesh( const std::string, const int );
  std::vector<double> Get_Mesh( const std::string );

  double Get_Algorithm( const std::string, const int );
  std::vector<double> Get_Algorithm( const std::string );

  void SaveXMLFile( const std::string & );

  int Get_NA();
  int Get_NK();

protected:
  void populate_vconstants( const std::string, std::map<std::string, std::vector<double>> & );
  void populate_parameter();

  tinyxml2::XMLDocument m_xml_doc;
  tinyxml2::XMLNode *m_pRoot;

  std::map<std::string, std::vector<double>> m_map_physics;
  std::map<std::string, std::vector<double>> m_map_mesh;
  std::map<std::string, std::vector<double>> m_map_algorithm;
  std::map<std::string, std::string> m_map_parameter;
};

#endif
