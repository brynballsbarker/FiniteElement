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
 * \brief Bernstein basis function.
 */
class Bernstein
{
  public:

    // Constructor.
    Bernstein( const int p_val,
             const int a_val );

    // Get the value of shape func at a point x.
    double evaluate( const double& coord );

    // Get the derivate at point x.
    double evaluateDeriv( const double& coord );

    // Get the second derivate at point x.
    double evaluateSecondDeriv( const double& coord );

  private:

    // Shape function order..
    int d_p;

    // Which shape function.
    int d_a;

    // p choose a.
    int d_choose;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_BERNSTEIN_HH

//---------------------------------------------------------------------------//
// end Bernstein.hh
//---------------------------------------------------------------------------//
