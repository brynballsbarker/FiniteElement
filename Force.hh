//---------------------------------------------------------------------------//
/*!
 * \file Force.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_FORCE_HH
#define FEA_FORCE_HH

#include "Element.hh"

#include <array>
#include <vector>

namespace FEA
{

//---------------------------------------------------------------------------//
/*!
 * \class Force
 *
 * \brief Force function. 
 */
class Force
{
  public:

    // Constructor.
    Force( const int f_order,
           const double c );

    // Get the value of f at a point x.
    double evaluate( const double& coord );
         
    bool checkLocal();

  private:

    // Order of f.
    int d_f_order;

    // Value of c.
    int d_c;

    // defined on local
    bool d_local;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_FORCE_HH

//---------------------------------------------------------------------------//
// end Force.hh
//---------------------------------------------------------------------------//
