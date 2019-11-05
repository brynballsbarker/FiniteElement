//---------------------------------------------------------------------------//
/*!
 * \file TrueSol.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_TRUESOL_HH
#define FEA_TRUESOL_HH

#include "Element.hh"

#include <array>
#include <vector>

namespace FEA
{

//---------------------------------------------------------------------------//
/*!
 * \class TrueSol
 *
 * \brief True solution to system.
 */
class TrueSol
{
  public:

    // Constructor.
    TrueSol( const int f_order,
             const double c );

    // Get the value of u at a point x.
    double evaluate( const double& coord );

    bool checkLocal();

  private:

    // Order of f.
    int d_f_order;

    // Value of c.
    int d_c;

    // Defined on local.
    bool d_local;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_TRUESOL_HH

//---------------------------------------------------------------------------//
// end TrueSol.hh
//---------------------------------------------------------------------------//
