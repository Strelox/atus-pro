
#ifndef _MY_FUNCTIONS
#define _MY_FUNCTIONS

#include <deal.II/base/function.h>
#include "eigenfunctions_AiryAi.h"
#include "eigenfunctions_Polar_HO.h"

using namespace std;
using namespace dealii;
   
template <int dim>
class CEigenfunctions : public Function<dim>
{
  public:
    CEigenfunctions ( unsigned int QN[3], std::vector<double> a ) : Function<dim>() 
    { 
      m_QNx=QN[0]; m_QNy=QN[1]; m_QNz=QN[2]; 
      m_fakx=a[0]; m_faky=a[1]; m_fakz=a[2]; 
    }
    virtual double value ( const Point<dim> &p, const unsigned int  component = 0) const;
  
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
double CEigenfunctions<dim>::value( const Point<dim> &p, const unsigned int component ) const
{
  double retval = (*AIRYEF[m_QNx])(pow(m_fakx,1.0/3.0),p(0)) * (*EF_PHO[m_QNy+51*m_QNz])(m_faky,p(1));
return retval;
}

/*************************************************************************************************/

template <int dim>
class CPotential : public Function<dim>
{
  public:
    CPotential ( const std::vector<double> a, const long l  ) : Function<dim>() { m_fakx=a[0]; m_faky=a[1]; m_fakz=a[2]; m_m = double(l*l); }
    virtual double value ( const Point<dim> &p, const unsigned int  component = 0) const;

    double m_fakx;
    double m_faky;
    double m_fakz;
    double m_m;
};
  
/*************************************************************************************************/
template <int dim>
double CPotential<dim>::value( const Point<dim> &p, const unsigned int component ) const
{
  double retval=0;
  double rq = p(1)*p(1);
  if( m_m == 0 )
  {
    retval = m_fakx*p(0) + m_faky*m_faky*rq;
  }
  else
  {
    if( rq > 0 ) 
      retval = m_fakx*p(0) + m_faky*m_faky*rq + m_m/rq;
    else
      retval = 0;
  }
return retval;
}
#endif