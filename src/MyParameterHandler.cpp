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


#include "MyParameterHandler.h"
#include "strtk.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <stdexcept>

MyParameterHandler::MyParameterHandler( const std::string filename )
{
  m_xml_doc.LoadFile( filename.c_str() );
  m_pRoot = m_xml_doc.FirstChild();

  populate_vconstants( "PHYSICS", m_map_physics );
  populate_vconstants( "MESH", m_map_mesh );
  populate_vconstants( "ALGORITHM", m_map_algorithm );

  populate_parameter();
}

void MyParameterHandler::populate_vconstants( const std::string section, std::map<std::string, std::vector<double>> &mymap )
{
  mymap.clear();

  std::string tmp, str;
  std::vector<std::string> vec;
  tinyxml2::XMLNode *pconstants = m_pRoot->FirstChildElement(section.c_str());

  double val;
  while ( pconstants != nullptr )
  {
    tinyxml2::XMLNode *pelem = pconstants->FirstChild();
    while ( pelem != nullptr )
    {
      vec.clear();
      str = pelem->ToElement()->Name();
      tmp = pelem->ToElement()->GetText();
      strtk::parse(tmp, ",", vec);

      for ( auto i : vec )
      {
        try
        {
          val = stod(i);
        }
        catch ( const std::invalid_argument &ia )
        {
          std::cerr << "Error Parsing xml file: Unable to convert " << i << " to double for element <" << str << "> in section " << section << '\n';
        }
        mymap[str].push_back(val);
      }
      pelem = pelem->NextSibling();
    }
    pconstants = pconstants->NextSiblingElement(section.c_str());
  }
}

void MyParameterHandler::populate_parameter()
{
  m_map_parameter.clear();

  std::string tmp, str;
  tinyxml2::XMLNode *pRoot = m_pRoot->FirstChild();

  while ( pRoot != nullptr )
  {
    str = pRoot->ToElement()->Name();
    //std::transform(str.begin(), str.end(),str.begin(), ::toupper);

    if ( str == "PHYSICS" || str == "MESH" || str == "ALGORITHM" )
    {
      pRoot = pRoot->NextSibling();
      continue;
    }

    tmp = pRoot->ToElement()->GetText();
    m_map_parameter[str] = tmp;
    pRoot = pRoot->NextSibling();
  }

//  for( auto i : m_map_parameter )
//    std::cout << i.first << ", " << i.second << std::endl;
}

/*
void MyParameterHandler::Setup_muParser( mu::Parser& mup )
{
  mup.ClearConst();
  //mup.ClearFun();

  for( auto i : m_map_constants )
  {
    //std::cout << i.first << ", " << i.second << std::endl;
    mup.DefineConst( i.first.c_str(), i.second );
  }

  mup.DefineFun("Heaviside", Heaviside, false);
  mup.DefineFun("rect", rect, false);
  mup.DefineFun("sign", sign, false);
}
*/

std::string MyParameterHandler::Get_Parameter( const std::string k )
{
  auto it = m_map_parameter.find(k);
  if ( it == m_map_parameter.end() ) throw std::string( "Error: Could not find the key: " + k + " in section PARAMETER." );
  return (*it).second;
}

void MyParameterHandler::Set_Physics( const std::string k, const std::vector<double> &newdata )
{
  std::string newstr;
  auto it = m_map_physics.find(k);
  if ( it == m_map_physics.end() ) throw std::string( "Error: Could not find the key " + k + " in section PHYSICS." );

  m_map_physics.at(k) = newdata;

  for ( unsigned i = 0; i < newdata.size() - 1; i++ )
    newstr += std::to_string(newdata[i]) + ",";
  newstr += std::to_string(newdata[newdata.size() - 1]);

  tinyxml2::XMLNode *pnode = m_pRoot->FirstChildElement("PHYSICS");

  while ( pnode != nullptr )
  {
    tinyxml2::XMLNode *pelem = pnode->FirstChild();
    while ( pelem != nullptr )
    {
      if ( k == pelem->ToElement()->Name() )
      {
        pelem->ToElement()->SetText(newstr.c_str());
        break;
      }
      pelem = pelem->NextSibling();
    }
    pnode = pnode->NextSiblingElement("PHYSICS");
  }
}

std::vector<double> MyParameterHandler::Get_Physics( const std::string k )
{
  std::vector<double> retval;
  auto it = m_map_physics.find(k);
  if ( it == m_map_physics.end() ) throw std::string( "Error: Could not find the key " + k + " in section PHYSICS." );
  retval = (*it).second;
  return retval;
}

double MyParameterHandler::Get_Physics( const std::string k, const int p )
{
  auto it = m_map_physics.find(k);
  if ( it == m_map_physics.end() ) throw std::string( "Error: Could not find the key " + k + " in section PHYSICS." );
  return (*it).second[p];
}

std::vector<double> MyParameterHandler::Get_Mesh( const std::string k )
{
  std::vector<double> retval;
  auto it = m_map_mesh.find(k);
  if ( it == m_map_mesh.end() ) throw std::string( "Error: Could not find the key " + k + " in section MESH." );
  retval = (*it).second;
  return retval;
}

double MyParameterHandler::Get_Mesh( const std::string k, const int p )
{
  auto it = m_map_mesh.find(k);
  if ( it == m_map_mesh.end() ) throw std::string( "Error: Could not find the key " + k + " in section MESH." );
  return (*it).second[p];
}

std::vector<double> MyParameterHandler::Get_Algorithm( const std::string k )
{
  std::vector<double> retval;
  auto it = m_map_algorithm.find(k);
  if ( it == m_map_algorithm.end() ) throw std::string( "Error: Could not find the key " + k + " in section ALGORITHM." );
  retval = (*it).second;
  return retval;
}

double MyParameterHandler::Get_Algorithm( const std::string k, const int p )
{
  auto it = m_map_algorithm.find(k);
  if ( it == m_map_algorithm.end() ) throw std::string( "Error: Could not find the key " + k + " in section ALGORITHM." );
  return (*it).second[p];
}

int MyParameterHandler::Get_NA()
{
  int retval = 10;
  auto it = m_map_algorithm.find("NA");
  if ( it != m_map_algorithm.end() ) retval = int((*it).second[0]);
  return retval;
}

int MyParameterHandler::Get_NK()
{
  int retval = 10;
  auto it = m_map_algorithm.find("NK");
  if ( it != m_map_algorithm.end() ) retval = int((*it).second[0]);
  return retval;
}

void MyParameterHandler::SaveXMLFile( const std::string &filename )
{
  m_xml_doc.SaveFile( filename.c_str() );
}
