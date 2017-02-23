
#ifndef _MY_FUNCTIONS
#define _MY_FUNCTIONS

#include <deal.II/base/function.h>
#include "eigenfunctions_AiryAi.h"
#include "eigenfunctions_HO.h"
#include "eigenfunctions_Polar_HO.h"

using namespace std;
using namespace dealii;

template <int dim>
class CEigenfunctions : public Function<dim>
{
  public:
    CEigenfunctions ( unsigned QN[3], std::vector<double> a ) : Function<dim>() 
    { 
      m_QNx=QN[0]; m_QNy=QN[1]; m_QNz=QN[2]; 
      m_fakx=a[0]; m_faky=a[1]; m_fakz=a[2]; 
    }
    virtual double value ( const Point<dim> &p, const unsigned component = 0) const;
  
    void change_params( double a[3] )
    {
      m_fakx=a[0]; m_faky=a[1]; m_fakz=a[2]; 
    }
  private:
    unsigned int m_QNx;
    unsigned int m_QNy;
    unsigned int m_QNz;
    double m_fakx;
    double m_faky;
    double m_fakz;
};

template <int dim>
double CEigenfunctions<dim>::value( const Point<dim> &p, const unsigned component ) const
{
  double retval;
  
  switch(dim)
  {
    case 1:
            #if POTENTIAL==1
              retval = (*EF_HO[m_QNx])(m_fakx,p(0));
            #endif
            #if POTENTIAL==2
              retval = (*AIRYEF[m_QNx])(pow(m_fakx,1.0/3.0),p(0));
            #endif
    break;
    case 2:
            #if POTENTIAL==1
              retval = (*EF_HO[m_QNx])(m_fakx,p(0)) * (*EF_HO[m_QNy])(m_faky,p(1));
            #endif
            #if POTENTIAL==2
              retval = (*AIRYEF[m_QNx])(pow(m_fakx,1.0/3.0),p(0)) * (*EF_HO[m_QNy])(m_faky,p(1));
            #endif
    break;
    case 3:
            #if POTENTIAL==1
              retval = (*EF_HO[m_QNx])(m_fakx,p(0)) * (*EF_HO[m_QNy])(m_faky,p(1)) * (*EF_HO[m_QNz])(m_fakz,p(2));
            #endif
            #if POTENTIAL==2
              retval = (*AIRYEF[m_QNx])(pow(m_fakx,1.0/3.0),p(0)) * (*EF_HO[m_QNy])(m_faky,p(1)) * (*EF_HO[m_QNz])(m_fakz,p(2));
            #endif
    break;
  }
return retval;
}

/*************************************************************************************************/

template <int dim>
class CPotential : public Function<dim>
{
  public:
    CPotential ( std::vector<double> a ) : Function<dim>() { m_fakx=a[0]; m_faky=a[1]; m_fakz=a[2];  }
    virtual double value ( const Point<dim> &p, const unsigned component = 0) const;

    double m_fakx;
    double m_faky;
    double m_fakz;
};
  
/*************************************************************************************************/
template <int dim>
double CPotential<dim>::value( const Point<dim> &p, const unsigned component ) const
{
  double retval;

  switch( dim )
  {
    case 1: 
            #if POTENTIAL==1
              retval = m_fakx*m_fakx*p(0)*p(0);
            #endif
            #if POTENTIAL==2
              retval = m_fakx*p(0);
            #endif
    break;
    case 2: 
            #if POTENTIAL==1
              retval = m_fakx*m_fakx*p(0)*p(0) + m_faky*m_faky*p(1)*p(1);
            #endif
            #if POTENTIAL==2
              retval = m_fakx*p(0) + m_faky*m_faky*p(1)*p(1);
            #endif
    break;
    case 3: 
            #if POTENTIAL==1
              retval = m_fakx*m_fakx*p(0)*p(0) + m_faky*m_faky*p(1)*p(1) + m_fakz*m_fakz*p(2)*p(2);
            #endif
            #if POTENTIAL==2
              retval = m_fakx*p(0) + m_faky*m_faky*p(1)*p(1) + m_fakz*m_fakz*p(2)*p(2);
            #endif
    break;
  }
return retval;
}
#endif