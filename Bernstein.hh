//---------------------------------------------------------------------------//
/*!
 * \file Bernstein.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_BERNSTEIN_HH
#define FEA_BERNSTEIN_HH

#include "Helpers.hh"

#include <array>
#include <vector>

namespace FEA
{

//---------------------------------------------------------------------------//
/*!
 * \class Bernstein
 *
 * \brief True solution to system.
 */
class Bernstein
{
  public:

    // Constructor.
    Bernstein( const int p_val,
             const int a_val );

    // Get the value of shpae func at a point x.
    double evaluate( const double& coord );

    // Get the derivate at point x.
    double evaluateDeriv( const double& coord );

  private:

    // Order of f.
    int d_p;

    // Value of c.
    int d_a;

    // p choose c.
    int d_choose;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_BERNSTEIN_HH

//---------------------------------------------------------------------------//
// end Bernstein.hh
//---------------------------------------------------------------------------//
