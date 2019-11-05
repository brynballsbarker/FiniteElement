//---------------------------------------------------------------------------//
/*!
 * \file Bernoulli.hh
 */
//---------------------------------------------------------------------------//

#ifndef FEA_BERNOULLI_HH
#define FEA_BERNOULLI_HH

#include "Element.hh"

#include <array>
#include <vector>

namespace FEA
{

//---------------------------------------------------------------------------//
/*!
 * \class Bernoulli
 *
 * \brief True solution to system.
 */
class Bernoulli
{
  public:

    // Constructor.
    Bernoulli( const int p_val,
             const double a_val );

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

    // Defined on local
    bool d_local;

};

} 

//---------------------------------------------------------------------------//

#endif // end FEA_BERNOULLI_HH

//---------------------------------------------------------------------------//
// end Bernoulli.hh
//---------------------------------------------------------------------------//
